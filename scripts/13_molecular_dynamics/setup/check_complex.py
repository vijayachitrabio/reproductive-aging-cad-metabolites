#!/usr/bin/env python3
"""Checkpoint: verify a merged protein+ligand complex.gro is sane before spending
any compute on equilibration/production.

Checks:
  1. Ligand internal geometry unchanged from the original docked pose (centroid-
     subtracted RMSD between the ligand atoms in complex.gro and the standalone
     ligand .gro) -- catches accidental corruption/reordering during the merge.
  2. Protein-ligand centroid separation in complex.gro -- catches a merge in
     which ligand and protein were translated independently and ended up
     far apart instead of docked together.
"""
import sys

def read_gro(path):
    with open(path) as f:
        lines = f.read().splitlines()
    natoms = int(lines[1].strip())
    atoms = []
    for line in lines[2:2 + natoms]:
        resname = line[5:10].strip()
        x, y, z = float(line[20:28]), float(line[28:36]), float(line[36:44])
        atoms.append((resname, x, y, z))
    return atoms

def centroid(coords):
    n = len(coords)
    return (sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n)

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def rmsd_centroid_subtracted(a_coords, b_coords):
    ca, cb = centroid(a_coords), centroid(b_coords)
    sq = 0.0
    for (x1, y1, z1), (x2, y2, z2) in zip(a_coords, b_coords):
        dx, dy, dz = (x1 - ca[0]) - (x2 - cb[0]), (y1 - ca[1]) - (y2 - cb[1]), (z1 - ca[2]) - (z2 - cb[2])
        sq += dx*dx + dy*dy + dz*dz
    return (sq / len(a_coords)) ** 0.5

complex_path, standalone_ligand_path = sys.argv[1], sys.argv[2]

complex_atoms = read_gro(complex_path)
lig_complex = [(x, y, z) for (resname, x, y, z) in complex_atoms if resname == "MOL"]
prot_complex = [(x, y, z) for (resname, x, y, z) in complex_atoms if resname != "MOL"]

lig_standalone_raw = read_gro(standalone_ligand_path)
lig_standalone = [(x, y, z) for (_, x, y, z) in lig_standalone_raw]

assert len(lig_complex) == len(lig_standalone), \
    f"ligand atom count mismatch: complex has {len(lig_complex)}, standalone has {len(lig_standalone)}"

rmsd = rmsd_centroid_subtracted(lig_complex, lig_standalone)
prot_c = centroid(prot_complex)
lig_c = centroid(lig_complex)
sep = sum((a - b) ** 2 for a, b in zip(prot_c, lig_c)) ** 0.5

print(f"protein atoms: {len(prot_complex)}")
print(f"ligand atoms:  {len(lig_complex)}")
print(f"protein centroid (nm): {prot_c}")
print(f"ligand centroid  (nm): {lig_c}")
print(f"protein-ligand centroid separation (nm): {sep:.4f}")
print(f"ligand internal-geometry RMSD vs docked pose (nm): {rmsd:.6f}")
print()

ok = True
if rmsd > 0.01:
    print(f"FAIL: ligand geometry changed by {rmsd:.4f} nm during merge -- expected ~0")
    ok = False
if sep > 3.0:
    print(f"FAIL: protein-ligand centroid separation is {sep:.2f} nm -- ligand is NOT in/near the pocket")
    ok = False

if ok:
    print("PASS: ligand geometry preserved and ligand sits near the protein.")
    sys.exit(0)
else:
    sys.exit(1)
