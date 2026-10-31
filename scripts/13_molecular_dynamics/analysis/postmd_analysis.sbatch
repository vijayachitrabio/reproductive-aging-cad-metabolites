#!/bin/bash
# Post-MD analysis of the apo FADS1 and FADS1-NU7441 (C2) 100-ns trajectories,
# matching the pipeline used for D5D-IN-326 (T). MM-GBSA is run for C2 only
# because apo FADS1 has no ligand.
#SBATCH --job-name=fads1_postmd
#SBATCH --account=uthpc
#SBATCH --partition=main
#SBATCH --time=5-00:00:00
#SBATCH --output=%x_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
set -eo pipefail
source "$HOME/miniconda3/bin/activate" md_2026
set -u
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
GMX="$HOME/miniconda3/envs/md_2026/bin/gmx"
MMPBSA="$HOME/miniconda3/envs/mmpbsa/bin/gmx_MMPBSA"

run_apo() {
    cd "${MD_ROOT}/fads1_md_apo"
    test -s apo_md_0_10.tpr; test -s apo_md_0_10.xtc
    printf 'q\n' | "$GMX" make_ndx -f apo_md_0_10.gro -o apo_index.ndx > apo_make_ndx.log 2>&1
    printf 'Protein\nProtein\n' | "$GMX" trjconv -s apo_md_0_10.tpr -f apo_md_0_10.xtc -o apo_md_center.xtc -pbc mol -center -n apo_index.ndx
    printf 'Protein\nProtein\n' | "$GMX" trjconv -s apo_md_0_10.tpr -f apo_md_center.xtc -o apo_md_fit.xtc -fit rot+trans -n apo_index.ndx
    printf 'Protein\n' | "$GMX" trjconv -s apo_md_0_10.tpr -f apo_md_fit.xtc -o apo_start.pdb -dump 0 -n apo_index.ndx
    printf 'Protein\n' | "$GMX" trjconv -s apo_md_0_10.tpr -f apo_md_fit.xtc -o apo_end.pdb -dump 100000 -n apo_index.ndx
    printf '4\n4\n' | "$GMX" rms -s apo_md_0_10.tpr -f apo_md_fit.xtc -tu ns -n apo_index.ndx -o apo_protein_rmsd.xvg
    printf '1\n' | "$GMX" rmsf -s apo_md_0_10.tpr -f apo_md_fit.xtc -n apo_index.ndx -res -o apo_rmsf.xvg
    printf '1\n' | "$GMX" gyrate -s apo_md_0_10.tpr -f apo_md_fit.xtc -n apo_index.ndx -o apo_gyrate.xvg
    printf '1\n1\n' | "$GMX" hbond-legacy -s apo_md_0_10.tpr -f apo_md_fit.xtc -tu ns -n apo_index.ndx -num apo_hbnum.xvg
    printf '1\n1\n' | "$GMX" hbond-legacy -s apo_md_0_10.tpr -f apo_md_fit.xtc -tu ns -n apo_index.ndx -hbm apo_hbond_matrix.xpm
    "$GMX" sasa -f apo_md_fit.xtc -s apo_md_0_10.tpr -n apo_index.ndx -o apo_area.xvg -odg apo_dgsolv.xvg -tu ns -surface Protein
    printf '4\n4\n' | "$GMX" covar -f apo_md_fit.xtc -s apo_md_0_10.tpr -n apo_index.ndx -o apo_eigenval.xvg -tu ns -v apo_eigenvec.trr
    printf '4\n4\n' | "$GMX" anaeig -v apo_eigenvec.trr -f apo_md_fit.xtc -eig apo_eigenval.xvg -s apo_md_0_10.tpr -n apo_index.ndx -first 1 -last 2 -2d apo_2dproj.xvg -comp apo_eigcomp.xvg -rmsf apo_eigrmsf.xvg -tu ns
}

run_c2() {
    cd "${MD_ROOT}/fads1_md_c2/C2_complex"
    test -s C2_md_0_100.tpr; test -s C2_md_0_100.xtc
    # C2 defaults match the completed MD workflow: 1=Protein, 4=Backbone, 13=MOL.
    printf '1 | 13\nname 20 Protein_LIG\n13 & ! a H*\nname 21 LIG_Heavy\nq\n' | "$GMX" make_ndx -f npt.gro -o c2_index_full.ndx > c2_make_ndx_full.log 2>&1
    printf 'Protein_LIG\nProtein_LIG\n' | "$GMX" trjconv -s C2_md_0_100.tpr -f C2_md_0_100.xtc -o c2_md_center.xtc -pbc mol -center -n c2_index_full.ndx
    printf 'Protein_LIG\nProtein_LIG\n' | "$GMX" trjconv -s C2_md_0_100.tpr -f c2_md_center.xtc -o c2_md_fit.xtc -fit rot+trans -n c2_index_full.ndx
    printf 'Protein_LIG\n' | "$GMX" trjconv -s C2_md_0_100.tpr -f c2_md_fit.xtc -o c2_start.pdb -dump 0 -n c2_index_full.ndx
    printf 'Protein_LIG\n' | "$GMX" trjconv -s C2_md_0_100.tpr -f c2_md_fit.xtc -o c2_end.pdb -dump 100000 -n c2_index_full.ndx
    printf '4\n4\n' | "$GMX" rms -s C2_md_0_100.tpr -f c2_md_fit.xtc -tu ns -n c2_index_full.ndx -o c2_protein_rmsd.xvg
    printf '4\n21\n' | "$GMX" rms -s C2_md_0_100.tpr -f c2_md_fit.xtc -tu ns -n c2_index_full.ndx -o c2_lig_rmsd_vs_avg.xvg
    printf '4\n21\n' | "$GMX" rms -s c2_start.pdb -f c2_md_fit.xtc -tu ns -n c2_index_full.ndx -o c2_lig_rmsd_vs_start.xvg
    printf '1\n' | "$GMX" rmsf -s C2_md_0_100.tpr -f c2_md_fit.xtc -n c2_index_full.ndx -res -o c2_rmsf.xvg
    printf '1\n' | "$GMX" gyrate -s C2_md_0_100.tpr -f c2_md_fit.xtc -n c2_index_full.ndx -o c2_gyrate.xvg
    printf '1\n13\n' | "$GMX" hbond-legacy -s C2_md_0_100.tpr -f c2_md_fit.xtc -tu ns -n c2_index_full.ndx -num c2_hbnum.xvg
    printf '1\n13\n' | "$GMX" hbond-legacy -s C2_md_0_100.tpr -f c2_md_fit.xtc -tu ns -n c2_index_full.ndx -hbm c2_hbond_matrix.xpm
    "$GMX" sasa -f c2_md_fit.xtc -s C2_md_0_100.tpr -n c2_index_full.ndx -o c2_area.xvg -odg c2_dgsolv.xvg -tu ns -surface Protein_LIG
    cat > c2_ie.mdp <<'EOF'
integrator      = md
dt              = 0.002
nsteps          = 0
nstlist         = 10
cutoff-scheme   = Verlet
pbc             = xyz
coulombtype     = PME
rcoulomb        = 1.2
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
rlist           = 1.2
energygrps      = Protein MOL
EOF
    "$GMX" grompp -f c2_ie.mdp -c npt.gro -t npt.cpt -p topol.top -n c2_index_full.ndx -o c2_ie.tpr -maxwarn 5
    "$GMX" mdrun -deffnm c2_ie -rerun C2_md_0_100.xtc -nb cpu -ntmpi 1 -ntomp "$SLURM_CPUS_PER_TASK" -v
    printf 'Coul-SR:Protein-MOL\nLJ-SR:Protein-MOL\n\n' | "$GMX" energy -f c2_ie.edr -o c2_interaction_energy.xvg
    printf '4\n4\n' | "$GMX" covar -f c2_md_fit.xtc -s C2_md_0_100.tpr -n c2_index_full.ndx -o c2_eigenval.xvg -tu ns -v c2_eigenvec.trr
    printf '4\n4\n' | "$GMX" anaeig -v c2_eigenvec.trr -f c2_md_fit.xtc -eig c2_eigenval.xvg -s C2_md_0_100.tpr -n c2_index_full.ndx -first 1 -last 2 -2d c2_2dproj.xvg -comp c2_eigcomp.xvg -rmsf c2_eigrmsf.xvg -tu ns
    cat > c2_mmpbsa.in <<'EOF'
&general
  startframe=8001, endframe=10000, interval=10,
  verbose=2, keep_files=0,
  PBRadii=2,
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=0,
  print_res="within 6",
/
EOF
    set +u
    source "$HOME/miniconda3/bin/activate" mmpbsa
    set -u
    export PATH="$HOME/miniconda3/envs/md_2026/bin:$PATH"
    printf 'q\n' | "$GMX" make_ndx -f npt.gro -o c2_index_mmpbsa.ndx > c2_make_ndx_mmpbsa.log 2>&1
    "$MMPBSA" -O -i c2_mmpbsa.in -cs C2_md_0_100.tpr -ci c2_index_mmpbsa.ndx -cg 1 13 -ct c2_md_fit.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA_C2_md_0_100.dat -do FINAL_DECOMP_MMPBSA_C2_md_0_100.dat -eo FINAL_RESULTS_MMPBSA_C2_md_0_100.csv -nogui
}

echo "FADS1 post-MD analysis started: $(date -Is)"
run_apo
run_c2
echo "FADS1 post-MD analysis completed: $(date -Is)"
