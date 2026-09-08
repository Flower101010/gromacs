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


def dump_frames(gmx, trr, output, natoms):
    with output.open('w') as log:
        dump_env = dict(os.environ, GMX_PRINT_LONGFORMAT='1')
        subprocess.run([str(gmx), 'dump', '-f', str(trr)], stdout=log,
                       stderr=subprocess.STDOUT, env=dump_env, check=True)
    frames = []
    frame_header = re.compile(r'^.* frame \d+:$', re.MULTILINE)
    vector = re.compile(r'\s*([xf])\s*\[\s*[\d\s]+\]\s*=\{(.*?)\}', re.DOTALL)
    box_vector = re.compile(r'\s*box\[\s*\d+\]=\{(.*?)\}', re.DOTALL)
    # dump writes its status/exit quote on the same stream and it can split a
    # long vector line. Its single-line messages have no force data.
    text = re.sub(r'\n(?:GROMACS reminds you:.*|Thanx for Using GROMACS.*)\n', '\n',
                  output.read_text())
    positions = list(frame_header.finditer(text))
    for index, header in enumerate(positions):
        frame_text = text[header.end():positions[index + 1].start() if index + 1 < len(positions) else None]
        frame = {'x': [], 'f': [], 'box': []}
        for match in vector.finditer(frame_text):
            frame[match[1]].append([float(v) for v in ''.join(match[2].split()).split(',')])
        # gmx dump prints the box matrix first using x[0..2] labels.
        if len(frame['x']) == natoms + 3:
            frame['x'] = frame['x'][3:]
        frame['box'] = [[float(v) for v in ''.join(match.group(1).split()).split(',')]
                        for match in box_vector.finditer(frame_text)]
        frames.append(frame)
    return frames


def molecular_ranges(gro):
    lines = gro.read_text().splitlines()
    natoms = int(lines[1])
    assert len(lines) >= natoms + 3
    ranges, start, residue = [], 0, lines[2][:5]
    for atom in range(1, natoms):
        next_residue = lines[2 + atom][:5]
        if next_residue != residue:
            ranges.append((start, atom))
            start, residue = atom, next_residue
    ranges.append((start, natoms))
    return natoms, ranges


def molecular_forces(frame, ranges):
    forces = frame['f']
    assert len(forces) == ranges[-1][1]
    mapped = []
    for start, end in ranges:
        mapped.extend([sum(forces[atom][d] for atom in range(start, end)) for d in range(3)])
    return mapped


def maximum_and_rms(first, second):
    difference = [a-b for a, b in zip(first, second)]
    return max(abs(x) for x in difference), math.sqrt(sum(x*x for x in difference) / len(difference))


def signed_xyz_mean(first, second, bead_count):
    return [sum(first[frame][3*bead + axis] - second[frame][3*bead + axis]
                for frame in range(len(first)) for bead in range(bead_count)) / (len(first) * bead_count)
            for axis in range(3)]


def block_sem(values, blocks=10):
    """Conservative contiguous-block SEM for each force component."""
    block_size = len(values) // blocks
    assert block_size > 0
    means = [[sum(values[block*block_size + sample][component] for sample in range(block_size))
              / block_size for block in range(blocks)] for component in range(len(values[0]))]
    return [math.sqrt(sum((value - sum(component_means) / blocks)**2 for value in component_means)
                      / (blocks - 1)) / math.sqrt(blocks)
            for component_means in means]


def percentile(values, p):
    values = sorted(values)
    return values[round(p * (len(values) - 1))]


def periodic_coordinate_maximum(first, second):
    maximum = 0.0
    for left, right in zip(first, second):
        assert len(left['box']) == len(right['box']) == 3
        for left_atom, right_atom in zip(left['x'], right['x']):
            for axis in range(3):
                box_length = left['box'][axis][axis]
                delta = left_atom[axis] - right_atom[axis]
                delta -= box_length * round(delta / box_length)
                maximum = max(maximum, abs(delta))
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gmx', type=Path, required=True)
    parser.add_argument('--gro', type=Path, required=True)
    parser.add_argument('--top', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=20)
    parser.add_argument('--stride', type=int, default=1,
                        help='write and compare every Nth MD step (default: 1)')
    args = parser.parse_args()
    if args.steps < 1 or args.stride < 1 or args.steps % args.stride:
        raise SystemExit('--steps and --stride must be positive, and steps must be divisible by stride')
    args.gmx, args.gro, args.top, args.output = (p.resolve() for p in
                                                   (args.gmx, args.gro, args.top, args.output))
    natoms, ranges = molecular_ranges(args.gro)
    expected_beads = len(ranges)
    args.output.mkdir(parents=True, exist_ok=False)
    template = Path(__file__).with_name('fixed-molecular-com-force-validation.mdp').read_text()
    for key, value in {'nsteps': args.steps, 'nstxout': args.stride, 'nstvout': 0,
                       'nstfout': args.stride, 'nstenergy': 0, 'nstlog': 0}.items():
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
    assert (magic, version, bead_count) == (b'GMXCOMF1', 1, expected_beads)
    saved = list(struct.iter_unpack('=d' + 'f' * (3 * bead_count), force_bytes[20:]))
    expected_frames = args.steps // args.stride + 1
    assert len(saved) == expected_frames

    production_frames = dump_frames(args.gmx, production / 'run.trr', production / 'dump.txt', natoms)
    assert len(production_frames) == expected_frames
    assert all(len(frame['x']) == natoms and len(frame['f']) == natoms for frame in production_frames)
    production_mapped = [molecular_forces(frame, ranges) for frame in production_frames]
    # A genuine MD force single point: the ON case initializes the fixed-COM
    # constraint machinery, but -nsteps 0 performs no integration/update.
    # This tests the actual feature, rather than only its rerun environment.
    single_points = []
    for name, enabled in (('single-point-off', False), ('single-point-on', True)):
        directory = args.output / name
        directory.mkdir()
        point_env = dict(clean_env)
        if enabled:
            point_env['GMX_FIXED_MOLECULAR_COM'] = '1'
        run(['mdrun', '-s', str(production / 'run.tpr'), '-deffnm', 'run', '-nsteps', '0',
             '-ntomp', '8', '-update', 'cpu', '-nb', 'gpu', '-pme', 'gpu', '-notunepme'],
            directory, point_env)
        frames = dump_frames(args.gmx, directory / 'run.trr', directory / 'dump.txt', natoms)
        assert len(frames) == 1
        assert len(frames[0]['x']) == natoms and len(frames[0]['f']) == natoms
        single_points.append(molecular_forces(frames[0], ranges))
    reruns, rerun_frames = [], []
    # ON here enters the fixed-COM setup but rerun does no update/projection.
    # It is a direct same-coordinate test for force-array contamination.
    for name, enabled in (('rerun-off', False), ('rerun-on', True)):
        directory = args.output / name
        directory.mkdir()
        rerun_env = dict(clean_env)
        if enabled:
            rerun_env['GMX_FIXED_MOLECULAR_COM'] = '1'
        run(['mdrun', '-s', str(production / 'run.tpr'), '-rerun', str(production / 'run.trr'),
             '-deffnm', 'run', '-ntomp', '8', '-nb', 'gpu', '-pme', 'gpu', '-notunepme'],
            directory, rerun_env)
        frames = dump_frames(args.gmx, directory / 'run.trr', directory / 'dump.txt', natoms)
        assert len(frames) == expected_frames
        assert all(len(frame['x']) == natoms and len(frame['f']) == natoms for frame in frames)
        reruns.append([molecular_forces(frame, ranges) for frame in frames])
        rerun_frames.append(frames)

    saved_values = [list(record[1:]) for record in saved]
    saved_vs_production = [maximum_and_rms(a, b) for a, b in zip(saved_values, production_mapped)]
    saved_vs_rerun = [maximum_and_rms(a, b) for a, b in zip(saved_values, reruns[0])]
    single_point_on_off = maximum_and_rms(single_points[0], single_points[1])
    rerun_on_off = [maximum_and_rms(a, b) for a, b in zip(reruns[0], reruns[1])]
    flat_saved = [value for frame in saved_values for value in frame]
    flat_rerun = [value for frame in reruns[0] for value in frame]
    mean_saved = [sum(frame[i] for frame in saved_values) / expected_frames for i in range(3 * bead_count)]
    mean_rerun = [sum(frame[i] for frame in reruns[0]) / expected_frames for i in range(3 * bead_count)]
    mean_difference = [a-b for a, b in zip(mean_saved, mean_rerun)]
    sem = block_sem(saved_values)
    sem_ratios = [abs(delta) / error for delta, error in zip(mean_difference, sem) if error != 0]
    report = {
        'frames': expected_frames,
        'stride_steps': args.stride,
        'beads': bead_count,
        'saved_vs_same_live_force_buffer_max_abs_kj_mol_nm': max(value[0] for value in saved_vs_production),
        'saved_vs_same_live_force_buffer_rms_kj_mol_nm':
            math.sqrt(sum(value[1]**2 for value in saved_vs_production) / expected_frames),
        'production_vs_rerun_off_coordinate_minimum_image_max_abs_nm':
            periodic_coordinate_maximum(production_frames, rerun_frames[0]),
        'saved_vs_rerun_max_abs_kj_mol_nm': max(value[0] for value in saved_vs_rerun),
        'saved_vs_rerun_rms_kj_mol_nm': math.sqrt(sum(value[1]**2 for value in saved_vs_rerun) / expected_frames),
        'fixed_com_on_vs_off_rerun_max_abs_kj_mol_nm': max(value[0] for value in rerun_on_off),
        'fixed_com_on_vs_off_rerun_rms_kj_mol_nm': math.sqrt(sum(value[1]**2 for value in rerun_on_off) / expected_frames),
        'fixed_com_on_vs_off_live_single_point_max_abs_kj_mol_nm': single_point_on_off[0],
        'fixed_com_on_vs_off_live_single_point_rms_kj_mol_nm': single_point_on_off[1],
        'largest_force_component_kj_mol_nm': max(max(abs(value) for value in flat_saved),
                                                  max(abs(value) for value in flat_rerun)),
        'saved_vs_rerun_max_relative_to_largest_component':
            max(value[0] for value in saved_vs_rerun) /
            max(max(abs(value) for value in flat_saved), max(abs(value) for value in flat_rerun)),
        'time_mean_max_abs_difference_kj_mol_nm': max(abs(value) for value in mean_difference),
        'time_mean_rms_difference_kj_mol_nm': math.sqrt(sum(value*value for value in mean_difference)
                                                         / len(mean_difference)),
        'signed_mean_delta_xyz_kj_mol_nm': signed_xyz_mean(saved_values, reruns[0], bead_count),
        'signed_mean_delta_vector_norm_kj_mol_nm': math.sqrt(sum(value*value for value in signed_xyz_mean(saved_values, reruns[0], bead_count))),
        'block_sem_blocks': 10,
        'saved_vs_rerun_mean_difference_over_block_sem_median': percentile(sem_ratios, 0.5),
        'saved_vs_rerun_mean_difference_over_block_sem_p95': percentile(sem_ratios, 0.95),
        'saved_vs_rerun_mean_difference_over_block_sem_max': max(sem_ratios),
        'per_frame_saved_vs_rerun_max_abs_kj_mol_nm': [value[0] for value in saved_vs_rerun],
        'times_ps': [record[0] for record in saved],
    }
    # This is a diagnostic bound for this exact single-precision GPU fixture.
    assert report['saved_vs_same_live_force_buffer_max_abs_kj_mol_nm'] < 0.02, report
    assert report['production_vs_rerun_off_coordinate_minimum_image_max_abs_nm'] < 2e-6, report
    assert report['saved_vs_rerun_max_abs_kj_mol_nm'] < 0.1, report
    assert report['fixed_com_on_vs_off_rerun_max_abs_kj_mol_nm'] < 0.1, report
    assert report['fixed_com_on_vs_off_live_single_point_max_abs_kj_mol_nm'] < 0.1, report
    (args.output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
