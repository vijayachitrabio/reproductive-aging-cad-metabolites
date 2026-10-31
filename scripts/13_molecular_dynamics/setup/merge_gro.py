#!/usr/bin/env python3
from pathlib import Path
import sys

def read_gro(path):
    lines = Path(path).read_text().splitlines()
    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    if len(atoms) != natoms:
        raise ValueError(f"{path}: expected {natoms} atoms, found {len(atoms)}")
    return atoms, lines[2 + natoms]

def renumber(line, number):
    return f"{line[:15]}{number:5d}{line[20:]}"

protein, box = read_gro(sys.argv[1])
ligand, _ = read_gro(sys.argv[2])
merged = [renumber(line, i) for i, line in enumerate(protein + ligand, 1)]
Path(sys.argv[3]).write_text(
    "C2 FADS1 complex\n" + f"{len(merged):5d}\n" +
    "\n".join(merged) + "\n" + box + "\n"
)
