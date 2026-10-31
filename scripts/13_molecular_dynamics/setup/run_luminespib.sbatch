#!/bin/bash
#SBATCH --job-name=c1_full
#SBATCH --account=uthpc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
source ~/miniconda3/bin/activate md_2026
set -e

cd ${MD_ROOT}/fads1_md_c1/C1_complex
rm -f nvt.* npt.* C1_md_0_100.* \#*\#

echo "=== NVT ==="
gmx grompp -f ../mdp_files/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -nb gpu -pme gpu -ntomp 8

echo "=== NPT ==="
gmx grompp -f ../mdp_files/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt -nb gpu -pme gpu -ntomp 8

echo "=== PRODUCTION 100ns ==="
gmx grompp -f ../mdp_files/md_100ns.mdp -c npt.gro -t npt.cpt -p topol.top -o C1_md_0_100.tpr -maxwarn 1
gmx mdrun -deffnm C1_md_0_100 -nb gpu -pme gpu -ntomp 8
