#!/bin/bash
# FADS1-D5D-IN-326 (T): post-MD analysis, then MM-GBSA, run sequentially on a
# 4-vCPU cloud VM without SLURM.
MD_ROOT="${MD_ROOT:-$HOME}"   # parent directory of the per-system MD run folders
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
set +e
cd ${MD_ROOT}/fads1_md/T_complex

echo "########## STAGE 1: post-MD analysis ##########"
source ~/miniforge3/bin/activate md
bash "$SCRIPT_DIR/analysis_d5d_in_326.sh"

echo "########## STAGE 2: MM-PBSA ##########"
source ~/miniforge3/bin/activate mmpbsa
export PATH="$PATH:$HOME/miniforge3/envs/md/bin"
bash "$SCRIPT_DIR/run_mmgbsa.sh" . T_md_0_100 mmpbsa.in ~/miniforge3/envs/md/bin/gmx

echo "T_ALL_PIPELINE_FINISHED"
