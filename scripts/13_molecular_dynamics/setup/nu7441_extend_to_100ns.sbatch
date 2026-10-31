#!/bin/bash
#SBATCH --job-name=c2_extend_100ns
#SBATCH --account=uthpc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
source ~/miniconda3/bin/activate md_2026
set -e
cd ${MD_ROOT}/fads1_md_c2/C2_complex

echo "=== Continuing C2 production from 50ns checkpoint to 100ns ==="
gmx mdrun -s C2_md_0_100.tpr -cpi C2_md_0_50.cpt -deffnm C2_md_0_50 \
    -nb gpu -pme gpu -ntmpi 1 -ntomp 8

echo "=== C2 EXTENSION TO 100NS FINISHED ==="
