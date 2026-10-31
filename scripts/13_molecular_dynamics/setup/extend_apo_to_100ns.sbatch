#!/bin/bash
#SBATCH --job-name=apo_extend
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=apo_extend_%j.log

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
source /etc/profile.d/modules.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate md_2026

cd ${MD_ROOT}/fads1_md_apo

# Extend from 10ns to 100ns (add 90ns = 90000 ps) via checkpoint continuation --
# no rerun of the existing 10ns, picks up exactly where it left off.
gmx convert-tpr -s apo_md_0_10.tpr -extend 90000 -o apo_md_0_100.tpr

gmx mdrun -s apo_md_0_100.tpr -cpi apo_md_0_10.cpt -deffnm apo_md_0_10 \
    -nb gpu -pme gpu -bonded gpu -update gpu -ntomp 8

echo "APO EXTENSION TO 100NS FINISHED"
