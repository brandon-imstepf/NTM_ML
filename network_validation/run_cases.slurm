#!/bin/bash -l
# Run numeric NTM cases on ROMULUS as a SLURM array, CHUNK rows per array task.
# Uses MATLAB R2022b at /opt/matlab/r2022b on the compute nodes (node01 and node07 checked).
#
#   sbatch --array=1-40 run_cases.slurm manifests/numeric.csv 1
#
# Smoke test first (4 time steps, separate output folder; sbatch passes env vars):
#   SNTM_SMOKE=1 SNTM_OUT_DIR=$PWD/outputs_smoke sbatch --array=1 run_cases.slurm manifests/smoke_numeric.csv 1
#
# init_rescale comparison (interior grid at 0.002, plus the Fig. 8 setting at 0.002 and 0.01):
#   sbatch --array=1-26 run_cases.slurm manifests/rescale_numeric.csv 1
#   sbatch --array=1-130%40 run_cases.slurm manifests/rescale_surrogate.csv 1
#
# Tighter ODE tolerances for rows whose steady-state solve fails at the defaults:
#   SNTM_RELTOL=1e-6 SNTM_ABSTOL=1e-8 sbatch --array=28-30,34,36-38 run_cases.slurm manifests/numeric.csv 1
#
# Array size = ceil(rows / CHUNK). Rerunning skips cases whose output exists.
# Surrogate cases run the same way (GP and NN-int. need the Statistics Toolbox):
#   sbatch --array=1-200 run_cases.slurm manifests/surrogate.csv 1
#SBATCH --job-name=sntm
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out

set -eo pipefail
MANIFEST=${1:?manifest csv}
CHUNK=${2:-1}

VAULT=/home/bimstepf/endor/vault
export NTM_ML_ROOT=$VAULT/NTM_ML
export NTM_SIM_ROOT=$VAULT/Tau_Transport-jtorok_dev
export SNTM_OUT_DIR=${SNTM_OUT_DIR:-$VAULT/NTM_ML/network_validation/outputs}

# MATLAB R2022b on the compute nodes (found by probe_cluster.slurm). Its
# modulefile lives in /act/modulefiles, which SLURM jobs do not see, so call
# the executable by path. Override with MATLAB_BIN=... if it moves.
MATLAB_BIN=${MATLAB_BIN:-/opt/matlab/r2022b/bin/matlab}
[[ -x "$MATLAB_BIN" ]] || { echo "MATLAB not found at $MATLAB_BIN on $(hostname)" >&2; exit 1; }

FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
LAST=$(( FIRST + CHUNK - 1 ))
cd "$NTM_ML_ROOT/network_validation"
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID rows=$FIRST-$LAST cpus=$SLURM_CPUS_PER_TASK"

"$MATLAB_BIN" -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); run_sntm_cases('$MANIFEST', $FIRST, $LAST)"
