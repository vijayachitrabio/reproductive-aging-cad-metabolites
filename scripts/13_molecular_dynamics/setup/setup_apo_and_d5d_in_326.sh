#!/bin/bash
MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
set -e
source ~/miniforge3/bin/activate md
FF=charmm36-feb2026_cgenff-5.0.ff
ROOT=${MD_ROOT}/fads1_md

# ============================================================
# APO (no ligand)
# ============================================================
mkdir -p "$ROOT/apo_v2"
cd "$ROOT/apo_v2"
cp "$ROOT/apo/FADS1.pdb" .
rm -rf "./$FF"; cp -r "$ROOT/$FF" "./$FF"

gmx pdb2gmx -f FADS1.pdb -o FADS1_processed.gro -water tip3p -ff "${FF%.ff}" -ignh
gmx editconf -f FADS1_processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f "$ROOT/mdp_files/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

gmx grompp -f "$ROOT/mdp_files/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -nb gpu

gmx grompp -f "$ROOT/mdp_files/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu

gmx grompp -f "$ROOT/mdp_files/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu

gmx grompp -f "$ROOT/mdp_files/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o apo_md_0_10.tpr -maxwarn 1
gmx mdrun -deffnm apo_md_0_10 -nb gpu -pme gpu -bonded gpu -update gpu

echo "=== APO_V2 DONE ==="

# ============================================================
# T-COMPLEX
# ============================================================
mkdir -p "$ROOT/T_complex_v2"
cd "$ROOT/T_complex_v2"
cp "$ROOT/T_complex/FADS1.pdb" .
cp -r "$ROOT/T_complex/T_LIG.acpype" .
rm -rf "./$FF"; cp -r "$ROOT/$FF" "./$FF"

gmx pdb2gmx -f FADS1.pdb -o FADS1_processed.gro -water tip3p -ff "${FF%.ff}" -ignh

python3 "$SCRIPT_DIR/merge_gro.py" FADS1_processed.gro T_LIG.acpype/T_LIG_GMX.gro complex.gro

# patch topology: insert ligand include after forcefield include, add to [molecules]
python3 - << 'PYEOF'
with open("topol.top") as f:
    content = f.read()
content = content.replace(
    '#include "./charmm36-feb2026_cgenff-5.0.ff/forcefield.itp"\n',
    '#include "./charmm36-feb2026_cgenff-5.0.ff/forcefield.itp"\n\n; Include ligand parameters\n#include "T_LIG.acpype/T_LIG_GMX.itp"\n'
)
content = content.rstrip() + "\nT_LIG               1\n"
with open("topol.top", "w") as f:
    f.write(content)
PYEOF

gmx editconf -f complex.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f "$ROOT/mdp_files/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

gmx grompp -f "$ROOT/mdp_files/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -nb gpu

gmx grompp -f "$ROOT/mdp_files/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu

gmx grompp -f "$ROOT/mdp_files/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu

gmx grompp -f "$ROOT/mdp_files/md_100ns.mdp" -c npt.gro -t npt.cpt -p topol.top -o T_md_0_100.tpr -maxwarn 1
gmx mdrun -deffnm T_md_0_100 -nb gpu -pme gpu -bonded gpu -update gpu

echo "=== T_COMPLEX_V2 DONE ==="
