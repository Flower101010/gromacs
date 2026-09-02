#!/usr/bin/env python3
"""Generate a 3-D generic-pull reference using whole molecular COM groups.

This is deliberately separate from CGWorkflow.  Zero-mass virtual sites may
remain in an index group: grompp uses the topology masses, so they contribute
zero to the pull COM while retaining the complete molecular membership.
"""

from __future__ import annotations

import argparse
from pathlib import Path


AXES = (("1", "0", "0"), ("0", "1", "0"), ("0", "0", "1"))


def read_molecules(gro: Path) -> list[list[int]]:
    lines = gro.read_text(encoding="ascii").splitlines()
    natoms = int(lines[1])
    molecules: list[list[int]] = []
    previous: tuple[int, str] | None = None
    for line in lines[2 : 2 + natoms]:
        key = (int(line[:5]), line[5:10].strip())
        atom = int(line[15:20])
        if key != previous:
            molecules.append([])
            previous = key
        molecules[-1].append(atom)
    if len(molecules) != 1000:
        raise ValueError(f"Expected 1000 molecules, found {len(molecules)}")
    return molecules


def write_index(path: Path, molecules: list[list[int]]) -> None:
    with path.open("w", encoding="ascii") as out:
        out.write("[ System ]\n")
        for atom in range(1, max(molecules[-1]) + 1):
            out.write(f"{atom}")
            out.write("\n" if atom % 15 == 0 else " ")
        out.write("\n")
        for i, atoms in enumerate(molecules, 1):
            out.write(f"[ MOL_{i:04d} ]\n")
            out.write(" ".join(map(str, atoms)) + "\n")


def write_mdp(path: Path, molecules: list[list[int]], nsteps: int) -> None:
    lines = [
        "; Whole-molecule mass-weighted COM generic-pull reference",
        "integrator = md", "dt = 0.002", f"nsteps = {nsteps}", "continuation = yes",
        "constraint-algorithm = lincs", "constraints = all-bonds", "lincs-iter = 2", "lincs-order = 6",
        "lincs-warnangle = 30", "cutoff-scheme = Verlet", "nstlist = 20", "rlist = 1.4",
        "coulombtype = PME", "rcoulomb = 1.4", "vdwtype = Cut-off", "rvdw = 1.4", "DispCorr = EnerPres",
        "pme-order = 4", "fourierspacing = 0.12", "ewald-rtol = 1e-5", "tcoupl = v-rescale",
        "tc-grps = System", "tau-t = 0.5", "ref-t = 298.15", "pcoupl = no", "gen-vel = no",
        "comm-mode = None", "nstcalcenergy = 100", "nstenergy = 0", "nstlog = 0", "nstxout = 0",
        "nstvout = 0", "nstfout = 0", "nstxout-compressed = 0", "pull = yes",
        f"pull-ngroups = {len(molecules)}", f"pull-ncoords = {3 * len(molecules)}", "pull-constr-tol = 1e-8",
        "pull-nstxout = 0", "pull-nstfout = 0", "pull-fout-average = no", "",
    ]
    for group, atoms in enumerate(molecules, 1):
        lines += [f"pull-group{group}-name = MOL_{group:04d}", f"pull-group{group}-pbcatom = {atoms[0]}"]
        for axis, vec in enumerate(AXES, 1):
            coordinate = 3 * (group - 1) + axis
            lines += [
                f"pull-coord{coordinate}-type = constraint", f"pull-coord{coordinate}-geometry = direction-periodic",
                f"pull-coord{coordinate}-groups = 0 {group}", f"pull-coord{coordinate}-vec = {' '.join(vec)}",
                f"pull-coord{coordinate}-origin = 0 0 0", f"pull-coord{coordinate}-start = yes",
            ]
        lines.append("")
    path.write_text("\n".join(lines), encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nsteps", type=int, default=50000)
    args = parser.parse_args()
    molecules = read_molecules(args.gro)
    args.output.mkdir(parents=True, exist_ok=True)
    write_index(args.output / "whole-molecule.ndx", molecules)
    write_mdp(args.output / "whole-molecule-pull.mdp", molecules, args.nsteps)


if __name__ == "__main__":
    main()
