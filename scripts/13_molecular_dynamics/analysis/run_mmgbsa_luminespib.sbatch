#!/bin/bash
#SBATCH --job-name=c1_mmpbsa
#SBATCH --partition=main
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=c1_mmpbsa_%j.log

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
MD_SCRIPTS="${MD_SCRIPTS:-$SLURM_SUBMIT_DIR}"   # this repository folder (scripts/13_molecular_dynamics/...)
source /etc/profile.d/modules.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mmpbsa

cd ${MD_ROOT}/fads1_md_c1/C1_complex
MPI_NP=8 bash "$MD_SCRIPTS/run_mmgbsa.sh" . C1_md_0_100 mmpbsa.in ~/miniconda3/envs/md_2026/bin/gmx

echo "C1 MMPBSA JOB FINISHED"
