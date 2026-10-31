#!/bin/bash
#SBATCH --job-name=c2_rebuild_step1
#SBATCH --account=uthpc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
MD_SCRIPTS="${MD_SCRIPTS:-$SLURM_SUBMIT_DIR}"   # this repository folder (scripts/13_molecular_dynamics/...)
source ~/miniconda3/bin/activate md_2026
set -e

cd ${MD_ROOT}/fads1_md_c2/C2_complex

echo "=== Step 1: merge protein.gro (unboxed) + C2_LIG.gro (docked pose) — same coordinate frame ==="
python3 "$MD_SCRIPTS/merge_gro.py" protein.gro C2_LIG.gro complex_unboxed.gro

echo "=== Step 2: box the MERGED complex (cubic, matching C1 exactly: 13.95 nm) ==="
gmx editconf -f complex_unboxed.gro -o complex.gro -bt cubic -box 13.95 13.95 13.95 -c

echo "=== Step 3: checkpoint — ligand geometry + protein-ligand proximity ==="
python3 "$MD_SCRIPTS/check_complex.py" complex.gro C2_LIG.gro
