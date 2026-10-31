#!/bin/bash
#SBATCH --job-name=c2_rebuild_full
#SBATCH --account=uthpc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=1-12:00:00
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
source ~/miniconda3/bin/activate md_2026
set -e

cd ${MD_ROOT}/fads1_md_c2/C2_complex

echo "=== Solvate (TIP3P) ==="
gmx solvate -cp complex.gro -cs spc216.gro -p topol.top -o solvated.gro

echo "=== Neutralize (genion) ==="
gmx grompp -f ${MD_ROOT}/fads1_md_c1/mdp_files/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ions.tpr -o ions.gro -p topol.top -neutral -pname NA -nname CL

echo "=== Energy minimization ==="
gmx grompp -f ${MD_ROOT}/fads1_md_c1/mdp_files/minim.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -nb gpu -ntomp 8

echo "=== NVT equilibration ==="
gmx grompp -f ${MD_ROOT}/fads1_md_c1/mdp_files/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -nb gpu -pme gpu -ntomp 8

echo "=== NPT equilibration ==="
gmx grompp -f ${MD_ROOT}/fads1_md_c1/mdp_files/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt -nb gpu -pme gpu -ntomp 8

echo "=== PRODUCTION 50 ns ==="
gmx grompp -f ${MD_ROOT}/fads1_md_c1/mdp_files/md_50ns.mdp -c npt.gro -t npt.cpt -p topol.top -o C2_md_0_50.tpr -maxwarn 1
gmx mdrun -deffnm C2_md_0_50 -nb gpu -pme gpu -ntomp 8

echo "=== DONE ==="
