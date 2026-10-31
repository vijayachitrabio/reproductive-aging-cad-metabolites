#!/bin/bash
# FADS1-Luminespib (C1) post-MD analysis: RMSD, RMSF, radius of gyration, H-bonds,
# SASA, interaction energy, PCA. Group numbers for this system
# (from gmx make_ndx on index_mmpbsa.ndx):
#   1=Protein  4=Backbone  13=MOL(ligand)
# Two custom groups added below: 20=Protein_LIG, 21=LIG_Heavy.
# Do not stop on the first failing step -- these are independent analyses.
MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
set +e
cd ${MD_ROOT}/fads1_md_c1/C1_complex

TPR=C1_md_0_100.tpr
XTC=C1_md_0_100.xtc
NDX=index_full.ndx

echo "=== Build extended index (Protein_LIG, LIG_Heavy) ==="
printf "1 | 13\nname 20 Protein_LIG\n13 & ! a H*\nname 21 LIG_Heavy\nq\n" | \
    gmx make_ndx -n index_mmpbsa.ndx -o "$NDX" > make_ndx_full.log 2>&1
grep -E "^\s*(20|21)\s" make_ndx_full.log

echo "=== PBC-center then rot+trans fit ==="
printf "Protein_LIG\nProtein_LIG\n" | gmx trjconv -s "$TPR" -f "$XTC" -o md_center.xtc -pbc mol -center -n "$NDX"
printf "Protein_LIG\nProtein_LIG\n" | gmx trjconv -s "$TPR" -f md_center.xtc -o md_fit.xtc -fit rot+trans -n "$NDX"

echo "=== Start/end reference frames (0 ns / 100 ns) ==="
printf "Protein_LIG\n" | gmx trjconv -s "$TPR" -f md_fit.xtc -o c1_start.pdb -dump 0 -n "$NDX"
printf "Protein_LIG\n" | gmx trjconv -s "$TPR" -f md_fit.xtc -o c1_end.pdb -dump 100000 -n "$NDX"

echo "=== RMSD: protein backbone (self-consistency, fit=Backbone, calc=Backbone) ==="
printf "4\n4\n" | gmx rms -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -o c1_protein_rmsd.xvg

echo "=== RMSD: ligand heavy atoms, fit=Backbone / calc=LIG_Heavy, vs whole-traj avg ==="
printf "4\n21\n" | gmx rms -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -o c1_lig_rmsd_vs_avg.xvg

echo "=== RMSD: ligand heavy atoms vs the t=0 start frame specifically ==="
printf "4\n21\n" | gmx rms -s c1_start.pdb -f md_fit.xtc -tu ns -n "$NDX" -o c1_lig_rmsd_vs_start.xvg

echo "=== RMSF: per-residue, protein ==="
printf "1\n" | gmx rmsf -s "$TPR" -f md_fit.xtc -n "$NDX" -res -o c1_rmsf.xvg

echo "=== Radius of gyration: protein ==="
printf "1\n" | gmx gyrate -s "$TPR" -f md_fit.xtc -n "$NDX" -o c1_gyrate.xvg

echo "=== H-bonds: protein-ligand (count + matrix) ==="
printf "1\n13\n" | gmx hbond -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -num c1_hbnum.xvg
printf "1\n13\n" | gmx hbond -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -hbm c1_hbond_matrix.xpm

echo "=== SASA: protein+ligand surface, plus solvation free energy estimate ==="
gmx sasa -f md_fit.xtc -s "$TPR" -n "$NDX" -o c1_area.xvg -odg c1_dgsolv.xvg -tu ns -surface Protein_LIG

echo "=== Interaction energy: rerun with energygrps=Protein MOL ==="
gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n "$NDX" -o ie.tpr -maxwarn 5
gmx mdrun -deffnm ie -rerun "$XTC" -nb cpu -v
echo "-- available energy terms (diagnostic) --"
printf "\n" | gmx energy -f ie.edr -o /tmp/energy_terms_dump.xvg 2>&1 | grep -E "^\s*[0-9]+\s"
printf "Coul-SR:Protein-MOL\nLJ-SR:Protein-MOL\n\n" | gmx energy -f ie.edr -o c1_interaction_energy.xvg

echo "=== PCA: covariance + eigenvector analysis on protein backbone ==="
printf "4\n4\n" | gmx covar -f md_fit.xtc -s "$TPR" -n "$NDX" -o c1_eigenval.xvg -tu ns -v c1_eigenvec.trr
printf "4\n4\n" | gmx anaeig -v c1_eigenvec.trr -f md_fit.xtc -eig c1_eigenval.xvg -s "$TPR" -n "$NDX" \
    -first 1 -last 2 -2d c1_2dproj.xvg -comp c1_eigcomp.xvg -rmsf c1_eigrmsf.xvg -tu ns

echo "C1 ANALYSIS PIPELINE FINISHED"
