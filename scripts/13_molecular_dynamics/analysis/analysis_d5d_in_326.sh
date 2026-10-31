#!/bin/bash
# FADS1-D5D-IN-326 (T) post-MD analysis: same pipeline as analysis_luminespib.sh.
# Group numbers for this system (gmx make_ndx on npt.gro): 1=Protein 4=Backbone
# 13=MOL (ligand, 43 atoms). Uses gmx hbond-legacy, because the newer gmx hbond
# module does not recognise the ligand donors/acceptors.
MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
set +e
cd ${MD_ROOT}/fads1_md/T_complex

TPR=T_md_0_100.tpr
XTC=T_md_0_100.xtc
NDX=index_full.ndx

echo "=== Build base + extended index (Protein_LIG, LIG_Heavy) ==="
printf "1 | 13\nname 20 Protein_LIG\n13 & ! a H*\nname 21 LIG_Heavy\nq\n" | \
    gmx make_ndx -f npt.gro -o "$NDX" > make_ndx_full.log 2>&1
grep -E "^\s*(20|21)\s" make_ndx_full.log

echo "=== PBC-center then rot+trans fit ==="
printf "Protein_LIG\nProtein_LIG\n" | gmx trjconv -s "$TPR" -f "$XTC" -o md_center.xtc -pbc mol -center -n "$NDX"
printf "Protein_LIG\nProtein_LIG\n" | gmx trjconv -s "$TPR" -f md_center.xtc -o md_fit.xtc -fit rot+trans -n "$NDX"

echo "=== Start/end reference frames (0 ns / 100 ns) ==="
printf "Protein_LIG\n" | gmx trjconv -s "$TPR" -f md_fit.xtc -o t_start.pdb -dump 0 -n "$NDX"
printf "Protein_LIG\n" | gmx trjconv -s "$TPR" -f md_fit.xtc -o t_end.pdb -dump 100000 -n "$NDX"

echo "=== RMSD: protein backbone (self-consistency, fit=Backbone, calc=Backbone) ==="
printf "4\n4\n" | gmx rms -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -o t_protein_rmsd.xvg

echo "=== RMSD: ligand heavy atoms, fit=Backbone / calc=LIG_Heavy, vs whole-traj avg ==="
printf "4\n21\n" | gmx rms -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -o t_lig_rmsd_vs_avg.xvg

echo "=== RMSD: ligand heavy atoms vs the t=0 start frame specifically ==="
printf "4\n21\n" | gmx rms -s t_start.pdb -f md_fit.xtc -tu ns -n "$NDX" -o t_lig_rmsd_vs_start.xvg

echo "=== RMSF: per-residue, protein ==="
printf "1\n" | gmx rmsf -s "$TPR" -f md_fit.xtc -n "$NDX" -res -o t_rmsf.xvg

echo "=== Radius of gyration: protein ==="
printf "1\n" | gmx gyrate -s "$TPR" -f md_fit.xtc -n "$NDX" -o t_gyrate.xvg

echo "=== H-bonds: protein-ligand (count + matrix) ==="
printf "1\n13\n" | gmx hbond-legacy -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -num t_hbnum.xvg
printf "1\n13\n" | gmx hbond-legacy -s "$TPR" -f md_fit.xtc -tu ns -n "$NDX" -hbm t_hbond_matrix.xpm

echo "=== SASA: protein+ligand surface, plus solvation free energy estimate ==="
gmx sasa -f md_fit.xtc -s "$TPR" -n "$NDX" -o t_area.xvg -odg t_dgsolv.xvg -tu ns -surface Protein_LIG

echo "=== Interaction energy: rerun with energygrps=Protein MOL ==="
gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n "$NDX" -o ie.tpr -maxwarn 5
gmx mdrun -deffnm ie -rerun "$XTC" -nb cpu -v
echo "-- available energy terms (diagnostic) --"
printf "\n" | gmx energy -f ie.edr -o /tmp/energy_terms_dump_t.xvg 2>&1 | grep -E "^\s*[0-9]+\s"
printf "Coul-SR:Protein-MOL\nLJ-SR:Protein-MOL\n\n" | gmx energy -f ie.edr -o t_interaction_energy.xvg

echo "=== PCA: covariance + eigenvector analysis on protein backbone ==="
printf "4\n4\n" | gmx covar -f md_fit.xtc -s "$TPR" -n "$NDX" -o t_eigenval.xvg -tu ns -v t_eigenvec.trr
printf "4\n4\n" | gmx anaeig -v t_eigenvec.trr -f md_fit.xtc -eig t_eigenval.xvg -s "$TPR" -n "$NDX" \
    -first 1 -last 2 -2d t_2dproj.xvg -comp t_eigcomp.xvg -rmsf t_eigrmsf.xvg -tu ns

echo "T ANALYSIS PIPELINE FINISHED"
