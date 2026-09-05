#!/usr/bin/env python3
"""Local leap-frog startup regression for the 500 SOL + 500 MEOH fixture.

Uses GROMACS long-format TRR dumps (nine significant digits for float32).
Requires the project's TIP4P/2005 / rigid TraPPE-UA topology and initial GRO.
All generated input, logs and results live in a new output directory.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess


def frames(path):
    result = []
    pattern = re.compile(r'\s*(x|v|f|box)\[\s*\d+\]=\{([^}]+)\}')
    for line in path.read_text().splitlines():
        if re.search(r'frame \d+:', line):
            result.append({key: [] for key in ('x', 'v', 'f', 'box')})
        match = pattern.match(line)
        if match:
            result[-1][match[1]].append([float(v) for v in match[2].split(',')])
    return result


def norm(v):
    return math.sqrt(sum(c*c for c in v))


def sub(a, b):
    return [x-y for x, y in zip(a, b)]


def image(v, box):
    return [c-box[d][d]*math.floor(c/box[d][d]+0.5) for d, c in enumerate(v)]


def com(frame, start, mass, key):
    vectors = frame[key][start:start+len(mass)]
    if key == 'x':
        vectors = [[a+b for a, b in zip(vectors[0], image(sub(v, vectors[0]), frame['box']))]
                   for v in vectors]
    return [sum(m*v[d] for m, v in zip(mass, vectors))/sum(mass) for d in range(3)]


def kinetic(frame, molecules):
    return sum(0.5*m*sum(c*c for c in frame['v'][start+i])
               for start, mass in molecules for i, m in enumerate(mass))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gmx', required=True)
    parser.add_argument('--gro', type=Path, required=True)
    parser.add_argument('--top', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.gro = args.gro.resolve()
    args.top = args.top.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    lines = args.gro.read_text().splitlines()
    assert int(lines[1]) == 3500
    molecules = [(4*i, [15.9994, 1.008, 1.008, 0]) for i in range(500)]
    molecules += [(2000+3*i, [15.035, 15.9994, 1.008]) for i in range(500)]
    for start, mass in molecules:
        names = ['OW', 'HW1', 'HW2', 'MW'] if start < 2000 else ['C', 'O', 'HO']
        assert [line[10:15].strip() for line in lines[2+start:2+start+len(mass)]] == names
    template = Path(__file__).with_name('fixed-molecular-com-force-validation.mdp').read_text()
    for option, value in {'nsteps': 20, 'nstxout': 1, 'nstvout': 1,
                          'nstenergy': 1, 'nstlog': 1}.items():
        template = re.sub(rf'(?m)^{option}\s*=.*$', f'{option} = {value}', template)
    template += '\nld-seed = 123456\n'
    env = dict(os.environ)
    env.pop('GMX_FIXED_MOLECULAR_COM', None)
    env.pop('GMX_FIXED_MOLECULAR_COM_FORCE_FILE', None)
    env['GMX_PRINT_LONGFORMAT'] = '1'

    def run(command, log, runenv):
        with log.open('w') as out:
            subprocess.run([args.gmx]+command, cwd=args.output, env=runenv,
                           stdout=out, stderr=subprocess.STDOUT, check=True)

    all_frames = {}
    for continuation in ('yes', 'no'):
        mdp = args.output / f'{continuation}.mdp'
        mdp.write_text(re.sub(r'(?m)^continuation\s*=.*$', f'continuation = {continuation}', template))
        tpr = args.output / f'{continuation}.tpr'
        run(['grompp', '-f', str(mdp), '-c', str(args.gro), '-p', str(args.top),
             '-o', str(tpr), '-po', str(args.output/f'{continuation}-processed.mdp'),
             '-maxwarn', '1'],  # Expected: comm-mode=None; the feature removes molecular translation.
            args.output/f'{continuation}-grompp.log', env)
        for enabled in (False, True):
            label = continuation + ('-on' if enabled else '-off')
            outdir = args.output/label
            outdir.mkdir()
            runenv = dict(env)
            if enabled:
                runenv['GMX_FIXED_MOLECULAR_COM'] = '1'
            run(['mdrun', '-s', str(tpr), '-deffnm', str(outdir/'run'), '-ntomp', '8',
                 '-update', 'cpu', '-nb', 'gpu', '-pme', 'gpu', '-notunepme'],
                outdir/'console.log', runenv)
            run(['dump', '-f', str(outdir/'run.trr')], outdir/'trajectory.txt', env)
            all_frames[label] = frames(outdir/'trajectory.txt')
            expected = '3000' if enabled else '6000'
            assert re.search(r'nrdf:\s*'+expected+r'\b', (outdir/'run.log').read_text())

    targets = [com(all_frames['yes-off'][0], start, mass, 'x') for start, mass in molecules]
    report = {}
    for continuation in ('yes', 'no'):
        trajectory = all_frames[continuation+'-on']
        assert len(trajectory) == 21
        errors, velocities, finite_difference, geometry = [], [], [], []
        for frame_index, frame in enumerate(trajectory):
            for index, (start, mass) in enumerate(molecules):
                errors.append(norm(image(sub(com(frame, start, mass, 'x'), targets[index]), frame['box'])))
                velocities.append(norm(com(frame, start, mass, 'v')))
                if frame_index:
                    # Output v at t_n is v_(n-1/2), the velocity connecting x_(n-1) to x_n.
                    previous = trajectory[frame_index-1]
                    for atom in range(start, start+3):
                        dx = image(sub(frame['x'][atom], previous['x'][atom]), frame['box'])
                        finite_difference.append(norm(sub([c/0.002 for c in dx], frame['v'][atom])))
                # continuation=yes deliberately trusts rounded input geometry at frame zero.
                if frame_index or continuation == 'no':
                    lengths = [0.09572, 0.09572, 0.15139] if start < 2000 else [0.143, 0.1948205404, 0.0945]
                    for (a, b), length in zip(((0, 1), (0, 2), (1, 2)), lengths):
                        distance = norm(image(sub(frame['x'][start+a], frame['x'][start+b]), frame['box']))
                        geometry.append(abs(distance-length))
        initial_velocity = max(norm(com(trajectory[0], start, mass, 'v')) for start, mass in molecules)
        report[continuation] = dict(frames=len(trajectory), max_com_nm=max(errors),
            rms_com_nm=norm(errors)/math.sqrt(len(errors)), max_com_velocity_nm_ps=max(velocities),
            initial_max_com_velocity_nm_ps=initial_velocity,
            max_leapfrog_finite_difference_nm_ps=max(finite_difference),
            max_rigid_distance_error_nm=max(geometry), initial_kinetic_kj_mol=kinetic(trajectory[0], molecules))
        assert max(errors) < 3e-6, report
        assert initial_velocity < 1e-6, report
        assert max(velocities) < 1e-3, report
        assert max(finite_difference) < 2e-3, report
        assert max(geometry) < 3e-6, report
    off = all_frames['yes-off'][0]
    translation = sum(0.5*sum(mass)*norm(com(off, start, mass, 'v'))**2 for start, mass in molecules)
    expected = kinetic(off, molecules)-translation
    report['initial_energy'] = dict(unprojected_kj_mol=kinetic(off, molecules),
        removed_translation_kj_mol=translation, expected_projected_kj_mol=expected)
    assert abs(report['yes']['initial_kinetic_kj_mol']-expected) < 0.01
    (args.output/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
