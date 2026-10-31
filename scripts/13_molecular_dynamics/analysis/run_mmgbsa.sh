#!/bin/bash
# Generic MM-GBSA + per-residue decomposition runner.
# Usage: run_mmgbsa.sh <system_dir> <deffnm> <mmpbsa_in> <gmx_bin>
# Assumes: <system_dir>/<deffnm>.tpr and .xtc exist (finished production),
# topol.top is in <system_dir>, ligand residue name is "MOL" (ACPYPE default,
# confirmed for both T_LIG and C1_LIG), receptor is "Protein".
# <gmx_bin> is the full path to the gmx/gmx_mpi binary in the PRODUCTION env
# (not the isolated mmpbsa env, which has no GROMACS installed) — e.g.
# ~/miniforge3/envs/md/bin/gmx on GCP, ~/miniconda3/envs/md_2026/bin/gmx_mpi on Rocket.
set -e
SYSDIR="$1"
DEFFNM="$2"
MMPBSA_IN="$3"
GMX_BIN="$4"

cd "$SYSDIR"

# gmx_MMPBSA internally shells out to a bare "gmx" for its own trajectory
# processing, which fails in the isolated mmpbsa env (no GROMACS installed
# there on purpose, to avoid touching the production env). Append the
# production env's bin dir so "gmx" resolves, while python/gmx_MMPBSA
# itself still resolve to the (PATH-first) mmpbsa env.
export PATH="$PATH:$(dirname "$GMX_BIN")"

# gmx_MMPBSA's sander re-evaluation is NOT periodic-boundary-aware -- feeding
# it the raw production trajectory causes intermittent BOND-energy overflow
# ("Some energy terms are undefined") on frames where a molecule wraps across
# the box edge. Fix: PBC-correct (whole molecules, centered) before use, same
# as the separate RMSD/RMSF/etc analysis pipeline already does. Reuse its
# md_fit.xtc if present; otherwise generate a PBC-corrected trajectory here.
TRAJ="${DEFFNM}.xtc"
if [ -f md_fit.xtc ]; then
    TRAJ="md_fit.xtc"
elif [ -f md_pbc.xtc ]; then
    TRAJ="md_pbc.xtc"
else
    echo "No PBC-corrected trajectory found -- generating md_pbc.xtc from ${DEFFNM}.xtc"
    printf "Protein\nSystem\n" | "$GMX_BIN" trjconv -s "${DEFFNM}.tpr" -f "${DEFFNM}.xtc" \
        -o md_pbc.xtc -pbc mol -center
    TRAJ="md_pbc.xtc"
fi
echo "Using PBC-corrected trajectory: $TRAJ"

# Build index file with only the default groups (System, Protein, MOL, etc.)
"$GMX_BIN" make_ndx -f npt.gro -o index_mmpbsa.ndx <<< "q" > make_ndx.log 2>&1

PROT_GRP=$(grep -oP '^\s*\K[0-9]+(?=\s+Protein\s+:)' make_ndx.log | head -1)
LIG_GRP=$(grep -oP '^\s*\K[0-9]+(?=\s+MOL\s+:)' make_ndx.log | head -1)

if [ -z "$PROT_GRP" ] || [ -z "$LIG_GRP" ]; then
    echo "FATAL: could not auto-detect Protein/MOL group numbers from make_ndx.log"
    grep -E "Protein|MOL" make_ndx.log
    exit 1
fi
echo "Using receptor group $PROT_GRP (Protein), ligand group $LIG_GRP (MOL)"

# MPI_NP > 1 fails for this workload -- multiple ranks race on
# shared sander scratch files during the per-frame ligand-only calculation,
# producing spurious "undefined energy terms" failures. Always run serial.
MPI_NP="${MPI_NP:-1}"
if [ "$MPI_NP" -gt 1 ]; then
    echo "WARNING: MPI_NP=$MPI_NP requested, but MPI parallelism is confirmed to cause" \
         "race-condition failures for this workload. Forcing serial (MPI_NP=1)."
fi
RUNNER="gmx_MMPBSA"

$RUNNER -O -i "$MMPBSA_IN" \
    -cs "${DEFFNM}.tpr" \
    -ci index_mmpbsa.ndx \
    -cg "$PROT_GRP" "$LIG_GRP" \
    -ct "$TRAJ" \
    -cp topol.top \
    -o "FINAL_RESULTS_MMPBSA_${DEFFNM}.dat" \
    -do "FINAL_DECOMP_MMPBSA_${DEFFNM}.dat" \
    -eo "FINAL_RESULTS_MMPBSA_${DEFFNM}.csv" \
    -nogui

echo "MMPBSA DONE FOR ${DEFFNM}"
