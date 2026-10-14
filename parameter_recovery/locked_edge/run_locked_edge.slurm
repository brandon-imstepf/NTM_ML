#!/bin/bash -l
# Single-edge locked-parameter datasets on ROMULUS, one locked set per array task.
# Array size = number of data rows in the sets file.
#
#   sbatch --array=1-9 run_locked_edge.slurm locked_sets.csv
#
# Smoke test first (3 x 3 grid, outputs_smoke/; sbatch passes env vars):
#   LOCKED_SMOKE=1 sbatch --array=1 run_locked_edge.slurm locked_sets.csv
#
# Rerunning skips sets whose output CSV exists.
#SBATCH --job-name=locked_edge
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%A_%a.out

set -eo pipefail
SETS=${1:?sets csv}

VAULT=/home/bimstepf/endor/vault
export NTM_SIM_ROOT=$VAULT/Tau_Transport-jtorok_dev

# Its modulefile is not visible to SLURM jobs, so call MATLAB by path
MATLAB_BIN=${MATLAB_BIN:-/opt/matlab/r2022b/bin/matlab}
[[ -x "$MATLAB_BIN" ]] || { echo "MATLAB not found at $MATLAB_BIN on $(hostname)" >&2; exit 1; }

cd "$VAULT/NTM_ML/parameter_recovery/locked_edge"
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID"

"$MATLAB_BIN" -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); make_locked_edge_data('$SETS', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_TASK_ID)"
