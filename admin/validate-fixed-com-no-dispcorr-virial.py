#!/usr/bin/env python3
"""Validate online no-DispCorr configurational virial in normal fixed-COM MD.

Runs normal fixed-COM NVT with DispCorr=EnerPres and DispCorr=no. At step zero
both evaluate the same coordinates, so after online subtraction the tensors
must agree within the independent GPU repeat floor. At later steps, GPU
roundoff can make independently run trajectories diverge; those runs still
verify each online average against its own instantaneous diagnostic samples.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def data_rows(path):
    return [[float(value) for value in line.split()] for line in path.read_text().splitlines()
            if line and not line.startswith('#')]


def read_average_virial(path):
    return [float(value) for line in path.read_text().splitlines()
            if line.startswith('virial ') for value in line.split()[1:]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gmx', type=Path, required=True)
    parser.add_argument('--gro', type=Path, required=True)
    parser.add_argument('--top', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=20)
    args = parser.parse_args()
    if args.steps < 0:
        raise SystemExit('--steps must be non-negative')
    args.gmx, args.gro, args.top, args.output = (path.resolve() for path in
                                                   (args.gmx, args.gro, args.top, args.output))
    args.output.mkdir(parents=True, exist_ok=False)
    base_mdp = Path(__file__).with_name('fixed-molecular-com-force-validation.mdp').read_text()
    for key, value in {'nsteps': args.steps, 'nstxout': 0, 'nstvout': 0,
                       'nstfout': 0, 'nstenergy': 0, 'nstlog': 0}.items():
        base_mdp = re.sub(rf'(?m)^{key}\s*=.*$', f'{key} = {value}', base_mdp)
    base_mdp += '\nld-seed = 123456\n'
    base_env = {key: value for key, value in os.environ.items()
                if not key.startswith('GMX_FIXED_MOLECULAR_COM')}

    def run(command, directory, env):
        result = subprocess.run([str(args.gmx)] + command, cwd=directory, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        (directory / (command[0] + '.console')).write_text(result.stdout)
        result.check_returncode()

    results = {}
    for label, setting in (('enerpres', 'EnerPres'), ('no', 'no'), ('no-repeat', 'no')):
        directory = args.output / label
        directory.mkdir()
        mdp = re.sub(r'(?mi)^DispCorr\s*=.*$', f'DispCorr = {setting}', base_mdp)
        (directory / 'run.mdp').write_text(mdp)
        run(['grompp', '-f', 'run.mdp', '-c', str(args.gro), '-p', str(args.top),
             '-o', 'run.tpr', '-maxwarn', '1'], directory, base_env)
        average = directory / 'average.txt'
        diagnostic = directory / 'diagnostic.txt'
        env = dict(base_env, GMX_FIXED_MOLECULAR_COM='1',
                   GMX_FIXED_MOLECULAR_COM_AVERAGE_FILE=str(average),
                   GMX_FIXED_MOLECULAR_COM_AVERAGE_STRIDE='1',
                   GMX_FIXED_MOLECULAR_COM_AVERAGE_VIRIAL_FILE=str(diagnostic))
        run(['mdrun', '-s', 'run.tpr', '-deffnm', 'run', '-ntomp', '8', '-update', 'cpu',
             '-nb', 'gpu', '-pme', 'gpu', '-notunepme'], directory, env)
        diagnostic_rows = data_rows(diagnostic)
        assert len(diagnostic_rows) == args.steps + 1 and all(len(row) == 11 for row in diagnostic_rows)
        results[label] = dict(diagnostic=diagnostic_rows, average=read_average_virial(average))
        offline = [sum(row[2 + 3*i + j] for row in diagnostic_rows) / len(diagnostic_rows)
                   for i in range(3) for j in range(3)]
        assert max(abs(left-right) for left, right in zip(offline, results[label]['average'])) < 1e-10

    difference = [[left-right for left, right in zip(a[2:], b[2:])]
                  for a, b in zip(results['enerpres']['diagnostic'], results['no']['diagnostic'])]
    repeat_difference = [[left-right for left, right in zip(a[2:], b[2:])]
                         for a, b in zip(results['no']['diagnostic'], results['no-repeat']['diagnostic'])]
    flat = [value for row in difference for value in row]
    report = {
        'frames': args.steps + 1,
        'online_no_dispcorr_enerpres_vs_no_max_abs_kj_mol': max(abs(value) for value in flat),
        'online_no_dispcorr_enerpres_vs_no_rms_kj_mol':
            (sum(value*value for value in flat) / len(flat))**0.5,
        'no_dispcorr_repeat_max_abs_kj_mol':
            max(abs(value) for row in repeat_difference for value in row),
        'no_dispcorr_repeat_rms_kj_mol':
            (sum(value*value for row in repeat_difference for value in row)
             / sum(len(row) for row in repeat_difference))**0.5,
        'average_tensor_max_abs_difference_kj_mol':
            max(abs(left-right) for left, right in zip(results['enerpres']['average'], results['no']['average'])),
    }
    assert report['online_no_dispcorr_enerpres_vs_no_max_abs_kj_mol'] < 0.2, report
    (args.output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
