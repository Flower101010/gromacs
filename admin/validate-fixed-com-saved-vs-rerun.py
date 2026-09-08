#!/usr/bin/env python3
"""Compare live fixed-COM mapped forces with independent mdrun -rerun forces.

This local test creates a short constrained trajectory whose every TRR frame
contains coordinates, box and forces. It sums rerun atom forces in exactly the
same whole-molecule topology order as GMX_FIXED_MOLECULAR_COM_FORCE_FILE.
It also repeats rerun to establish the local GPU roundoff scale.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess


def dump_frames(gmx, trr, output):
    with output.open('w') as log:
        subprocess.run([str(gmx), 'dump', '-f', str(trr)], stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    frames = []
    frame_header = re.compile(r'^.* frame \d+:$', re.MULTILINE)
    vector = re.compile(r'\s*([xf])\[\s*\d+\]=\{(.*?)\}', re.DOTALL)
    # dump writes its status/exit quote on the same stream and it can split a
    # long vector line. Its single-line messages have no force data.
    text = re.sub(r'\n(?:GROMACS reminds you:.*|Thanx for Using GROMACS.*)\n', '\n',
                  output.read_text())
    positions = list(frame_header.finditer(text))
    for index, header in enumerate(positions):
        frame_text = text[header.end():positions[index + 1].start() if index + 1 < len(positions) else None]
        frame = {'x': [], 'f': []}
        for match in vector.finditer(frame_text):
            frame[match[1]].append([float(v) for v in ''.join(match[2].split()).split(',')])
        # gmx dump prints the box matrix first using x[0..2] labels.
        if len(frame['x']) == 3503:
            frame['x'] = frame['x'][3:]
        frames.append(frame)
    return frames


def molecular_forces(frame):
    forces = frame['f']
    assert len(forces) == 3500
    mapped = []
    for molecule in range(500):
        mapped.extend([sum(forces[4*molecule + atom][d] for atom in range(4)) for d in range(3)])
    for molecule in range(500):
        start = 2000 + 3*molecule
        mapped.extend([sum(forces[start + atom][d] for atom in range(3)) for d in range(3)])
    return mapped


def maximum_and_rms(first, second):
    difference = [a-b for a, b in zip(first, second)]
    return max(abs(x) for x in difference), math.sqrt(sum(x*x for x in difference) / len(difference))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gmx', type=Path, required=True)
    parser.add_argument('--gro', type=Path, required=True)
    parser.add_argument('--top', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=20)
    args = parser.parse_args()
    if args.steps < 1:
        raise SystemExit('--steps must be positive')
    args.gmx, args.gro, args.top, args.output = (p.resolve() for p in
                                                   (args.gmx, args.gro, args.top, args.output))
    args.output.mkdir(parents=True, exist_ok=False)
    template = Path(__file__).with_name('fixed-molecular-com-force-validation.mdp').read_text()
    for key, value in {'nsteps': args.steps, 'nstxout': 1, 'nstvout': 0,
                       'nstfout': 1, 'nstenergy': 0, 'nstlog': 0}.items():
        template = re.sub(rf'(?m)^{key}\s*=.*$', f'{key} = {value}', template)
    (args.output / 'run.mdp').write_text(template + '\nld-seed = 123456\n')
    clean_env = {k: v for k, v in os.environ.items() if not k.startswith('GMX_FIXED_MOLECULAR_COM')}

    def run(command, cwd, env):
        with (cwd / 'console.log').open('w') as log:
            subprocess.run([str(args.gmx)] + command, cwd=cwd, env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True)

    production = args.output / 'production'
    production.mkdir()
    run(['grompp', '-f', str(args.output / 'run.mdp'), '-c', str(args.gro), '-p', str(args.top),
         '-o', 'run.tpr', '-maxwarn', '1'], production, clean_env)
    production_env = dict(clean_env, GMX_FIXED_MOLECULAR_COM='1',
                          GMX_FIXED_MOLECULAR_COM_FORCE_FILE=str(production / 'saved-force.bin'))
    run(['mdrun', '-s', 'run.tpr', '-deffnm', 'run', '-ntomp', '8', '-update', 'cpu',
         '-nb', 'gpu', '-pme', 'gpu', '-notunepme'], production, production_env)
    force_bytes = (production / 'saved-force.bin').read_bytes()
    magic, version, bead_count = struct.unpack_from('=8sIQ', force_bytes)
    assert (magic, version, bead_count) == (b'GMXCOMF1', 1, 1000)
    saved = list(struct.iter_unpack('=d' + 'f' * (3 * bead_count), force_bytes[20:]))
    expected_frames = args.steps + 1
    assert len(saved) == expected_frames

    reruns = []
    for name in ('rerun-a', 'rerun-b'):
        directory = args.output / name
        directory.mkdir()
        run(['mdrun', '-s', str(production / 'run.tpr'), '-rerun', str(production / 'run.trr'),
             '-deffnm', 'run', '-ntomp', '8', '-nb', 'gpu', '-pme', 'gpu', '-notunepme'],
            directory, clean_env)
        frames = dump_frames(args.gmx, directory / 'run.trr', directory / 'dump.txt')
        assert len(frames) == expected_frames
        assert all(len(frame['x']) == 3500 and len(frame['f']) == 3500 for frame in frames)
        reruns.append([molecular_forces(frame) for frame in frames])

    saved_values = [list(record[1:]) for record in saved]
    saved_vs_rerun = [maximum_and_rms(a, b) for a, b in zip(saved_values, reruns[0])]
    rerun_vs_rerun = [maximum_and_rms(a, b) for a, b in zip(reruns[0], reruns[1])]
    flat_saved = [value for frame in saved_values for value in frame]
    flat_rerun = [value for frame in reruns[0] for value in frame]
    mean_saved = [sum(frame[i] for frame in saved_values) / expected_frames for i in range(3 * bead_count)]
    mean_rerun = [sum(frame[i] for frame in reruns[0]) / expected_frames for i in range(3 * bead_count)]
    mean_difference = [a-b for a, b in zip(mean_saved, mean_rerun)]
    report = {
        'frames': expected_frames,
        'beads': 1000,
        'saved_vs_rerun_max_abs_kj_mol_nm': max(value[0] for value in saved_vs_rerun),
        'saved_vs_rerun_rms_kj_mol_nm': math.sqrt(sum(value[1]**2 for value in saved_vs_rerun) / expected_frames),
        'rerun_repeat_max_abs_kj_mol_nm': max(value[0] for value in rerun_vs_rerun),
        'rerun_repeat_rms_kj_mol_nm': math.sqrt(sum(value[1]**2 for value in rerun_vs_rerun) / expected_frames),
        'largest_force_component_kj_mol_nm': max(max(abs(value) for value in flat_saved),
                                                  max(abs(value) for value in flat_rerun)),
        'saved_vs_rerun_max_relative_to_largest_component':
            max(value[0] for value in saved_vs_rerun) /
            max(max(abs(value) for value in flat_saved), max(abs(value) for value in flat_rerun)),
        'time_mean_max_abs_difference_kj_mol_nm': max(abs(value) for value in mean_difference),
        'time_mean_rms_difference_kj_mol_nm': math.sqrt(sum(value*value for value in mean_difference)
                                                         / len(mean_difference)),
        'per_frame_saved_vs_rerun_max_abs_kj_mol_nm': [value[0] for value in saved_vs_rerun],
        'times_ps': [record[0] for record in saved],
    }
    # This is a diagnostic bound for this exact single-precision GPU fixture.
    assert report['saved_vs_rerun_max_abs_kj_mol_nm'] < 0.1, report
    (args.output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
