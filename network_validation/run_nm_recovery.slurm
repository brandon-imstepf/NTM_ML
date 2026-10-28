#!/bin/bash -l
# Nelder-Mead single-parameter recovery on the network (nelder_mead_network_recovery.m),
# one parameter per array task: 1 = gamma1, 2 = lambda1, 3 = delta, 4 = epsilon.
# Each task runs its simulations one after the other (about 5 min each for MLP,
# 40 at most), so a task takes up to about 4 h.
#
#   sbatch --array=1-4 run_nm_recovery.slurm g0.0005_l0.05_d20_e80
#   sbatch --array=1-4 run_nm_recovery.slurm g0.0005_l0.05_d20_e80 mlp l2f
#
# Rerunning after a timeout resumes from the simulations saved in outputs_nm/.
#SBATCH --job-name=nm
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out

set -eo pipefail
TRUTH=${1:?truth_id from manifests/recovery_truths.csv}
METHOD=${2:-mlp}
CONV=${3:-l2f}
PARAMS=(gamma1 lambda1 delta epsilon)
PARAM=${PARAMS[$((SLURM_ARRAY_TASK_ID - 1))]}

VAULT=/home/bimstepf/endor/vault
export NTM_ML_ROOT=$VAULT/NTM_ML
export NTM_SIM_ROOT=$VAULT/Tau_Transport-jtorok_dev

# MATLAB R2022b by path (see run_cases.slurm)
MATLAB_BIN=${MATLAB_BIN:-/opt/matlab/r2022b/bin/matlab}
[[ -x "$MATLAB_BIN" ]] || { echo "MATLAB not found at $MATLAB_BIN on $(hostname)" >&2; exit 1; }

cd "$NTM_ML_ROOT/network_validation"
echo "host=$(hostname) truth=$TRUTH param=$PARAM method=$METHOD conv=$CONV"

"$MATLAB_BIN" -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK}); nelder_mead_network_recovery('$TRUTH', '$PARAM', '$METHOD', '$CONV')"
