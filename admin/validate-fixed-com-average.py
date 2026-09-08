#!/usr/bin/env python3
"""Local real-timestep online/offline averaging regression (GPU NB/PME, CPU update).

Requires the existing fixed-COM GRO/topology fixture. All outputs go to a new directory.
Uses only the Python standard library. The optional diagnostic records do_force's
actual nine virial components, not the constraint-inclusive EDR pressure virial.
"""
import argparse
import json
import os
from pathlib import Path
import re
import struct
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('gmx', 'gro', 'top', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    args.gmx, args.gro, args.top, args.output = (
        p.resolve() for p in (args.gmx, args.gro, args.top, args.output))
    args.output.mkdir(parents=True, exist_ok=False)
    env = {k: v for k, v in os.environ.items() if not k.startswith('GMX_FIXED_MOLECULAR_COM')}
    env['GMX_FIXED_MOLECULAR_COM'] = '1'
    template = Path(__file__).with_name('fixed-molecular-com-force-validation.mdp').read_text()
    # Nonzero init-step verifies that skip is run-relative, not an absolute MD step.
    template += '\ninit-step = 10\nld-seed = 123456\n'

    def command(cmd, directory, runenv, success=True):
        with (directory / (cmd[0] + '.console')).open('w') as log:
            result = subprocess.run([str(args.gmx)] + cmd, cwd=directory, env=runenv,
                                    stdout=log, stderr=subprocess.STDOUT)
        assert (result.returncode == 0) == success, (directory, cmd)

    report = {}
    streams = {}
    for label, skip, stride, force_output, energy_stride in (
            ('all', 0, 1, 1, 1), ('skip', 3, 1, 1, 1), ('sparse', 2, 3, 1, 100),
            ('empty', 9, 1, 1, 1), ('no_trajectory', 0, 1, 0, 1),
            ('average_off', 0, 1, 1, 1), ('average_off_repeat', 0, 1, 1, 1)):
        directory = args.output / label
        directory.mkdir()
        mdp = template
        for key, value in {'nsteps': 8, 'nstfout': force_output,
                           'nstcalcenergy': energy_stride}.items():
            mdp = re.sub(rf'(?m)^{key}\s*=.*$', f'{key} = {value}', mdp)
        (directory / 'run.mdp').write_text(mdp)
        command(['grompp', '-f', 'run.mdp', '-c', str(args.gro), '-p', str(args.top),
                 '-o', 'run.tpr', '-maxwarn', '1'], directory, env)
        runenv = dict(env)
        if force_output:
            runenv['GMX_FIXED_MOLECULAR_COM_FORCE_FILE'] = str(directory / 'force.bin')
        if not label.startswith('average_off'):
            runenv.update(GMX_FIXED_MOLECULAR_COM_AVERAGE_FILE=str(directory / 'average.txt'),
                          GMX_FIXED_MOLECULAR_COM_AVERAGE_SKIP_STEPS=str(skip),
                          GMX_FIXED_MOLECULAR_COM_AVERAGE_STRIDE=str(stride),
                          GMX_FIXED_MOLECULAR_COM_AVERAGE_VIRIAL_FILE=str(directory / 'virial.txt'))
        command(['mdrun', '-s', 'run.tpr', '-deffnm', 'run', '-ntomp', '8',
                 '-update', 'cpu', '-nb', 'gpu', '-pme', 'gpu', '-notunepme'],
                directory, runenv)
        if force_output:
            binary = (directory / 'force.bin').read_bytes()
            magic, version, beads = struct.unpack_from('=8sIQ', binary)
            assert (magic, version) == (b'GMXCOMF1', 1)
            records = list(struct.iter_unpack('=d' + 'f' * (3 * beads), binary[20:]))
            assert len(records) == 9
            streams[label] = records
        if label.startswith('average_off'):
            assert not (directory / 'average.txt').exists()
            continue
        lines = [line.split() for line in (directory / 'average.txt').read_text().splitlines()
                 if line and not line.startswith('#')]
        count = int(next(line[1] for line in lines if line[0] == 'samples'))
        selected = list(range(skip, 9, stride))
        assert count == len(selected)
        instantaneous = [list(map(float, line.split())) for line in
                         (directory / 'virial.txt').read_text().splitlines()]
        assert [int(row[0]) for row in instantaneous] == selected
        force = [float(v) for line in lines if line[0] == 'force' for v in line[2:]]
        virial = [float(v) for line in lines if line[0] == 'virial' for v in line[1:]]
        if not count:
            assert not force and not virial
            report[label] = {'samples': 0}
            continue
        # Independent offline reduction of the existing binary estimator, all components.
        reference = streams['all'] if not force_output else records
        offline_force = [sum(reference[i][j+1] for i in selected) / count
                         for j in range(len(force))]
        offline_virial = [sum(row[j+2] for row in instantaneous) / count for j in range(9)]
        force_error = max(abs(a-b) for a, b in zip(force, offline_force))
        virial_error = max(abs(a-b) for a, b in zip(virial, offline_virial))
        # Separate GPU trajectories need a roundoff tolerance; same-run means do not.
        assert force_error < (0.1 if not force_output else 1e-10), force_error
        assert virial_error < 1e-10, virial_error
        if not force_output:
            assert not (directory / 'run.trr').exists()
        report[label] = dict(samples=count, max_force_error=force_error,
                             max_virial_error=virial_error, sampled_relative_steps=selected)

    difference = max(abs(a-b) for row_a, row_b in zip(streams['all'], streams['average_off'])
                     for a, b in zip(row_a, row_b))
    repeat_difference = max(abs(a-b) for row_a, row_b in
                            zip(streams['average_off'], streams['average_off_repeat'])
                            for a, b in zip(row_a, row_b))
    assert difference < 0.1 and repeat_difference < 0.1, (difference, repeat_difference)
    report['average_off'] = dict(max_instantaneous_force_difference=difference,
                                off_repeat_max_difference=repeat_difference)
    # Fail loudly for invalid configuration, before entering production MD.
    for label, updates in (
            ('bad_skip', {'GMX_FIXED_MOLECULAR_COM_AVERAGE_SKIP_STEPS': '-1'}),
            ('bad_stride', {'GMX_FIXED_MOLECULAR_COM_AVERAGE_STRIDE': '0'}),
            ('no_fixed_com', {'GMX_FIXED_MOLECULAR_COM': None})):
        directory = args.output / label
        directory.mkdir()
        runenv = dict(env, GMX_FIXED_MOLECULAR_COM_AVERAGE_FILE=str(directory / 'average.txt'))
        for key, value in updates.items():
            if value is None:
                runenv.pop(key)
            else:
                runenv[key] = value
        command(['mdrun', '-s', str(args.output / 'all' / 'run.tpr'), '-ntomp', '8',
                 '-update', 'cpu', '-nb', 'gpu', '-pme', 'gpu'], directory, runenv, success=False)
        expected = {'bad_skip': 'AVERAGE_SKIP_STEPS must be an integer',
                    'bad_stride': 'AVERAGE_STRIDE must be an integer',
                    'no_fixed_com': 'averaging requires GMX_FIXED_MOLECULAR_COM'}[label]
        assert expected in (directory / 'mdrun.console').read_text()
    (args.output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
