#!/bin/bash -l
# Run numeric NTM cases on ROMULUS as a SLURM array, CHUNK rows per array task.
# MATLAB (R2022b) exists only on the oldnodes partition (node01-05), via `module load matlab`.
#
#   MATLAB_BIN=/path/to/bin/matlab sbatch --array=1-40 run_cases.slurm manifests/numeric.csv 1
#
# Smoke test first (4 time steps, separate output folder; sbatch passes env vars):
#   MATLAB_BIN=/path/to/bin/matlab SNTM_SMOKE=1 SNTM_OUT_DIR=$PWD/outputs_smoke sbatch --array=1 run_cases.slurm manifests/smoke_numeric.csv 1
#
# Array size = ceil(rows / CHUNK). Rerunning skips cases whose output exists.
# Surrogate cases need Statistics Toolbox and run locally (run_local.m).
#SBATCH --job-name=sntm
#SBATCH --partition=oldnodes   # MATLAB module exists only on node01-05
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

# MATLAB: full path if MATLAB_BIN is set (inside SLURM jobs the module system
# inherits the head node's MODULEPATH and cannot see matlab), else the module.
if [[ -z "${MATLAB_BIN:-}" ]]; then
    module load matlab
    MATLAB_BIN=$(command -v matlab) || { echo "matlab not found; set MATLAB_BIN" >&2; exit 1; }
fi

FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
LAST=$(( FIRST + CHUNK - 1 ))
cd "$NTM_ML_ROOT/network_validation"
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID rows=$FIRST-$LAST cpus=$SLURM_CPUS_PER_TASK"

"$MATLAB_BIN" -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); run_sntm_cases('$MANIFEST', $FIRST, $LAST)"
