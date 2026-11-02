#!/bin/bash -l
# Nelder-Mead single-parameter recovery at 60 days (nelder_mead_network_recovery.m),
# MLP surrogate, l2f, against the 60-day truths of manifests/recovery_truth_mlp_l2f.csv.
# One array task per truth and parameter: task = 4 * (truth - 1) + parameter, with
# truth = row of the l2f truths in manifests/recovery_truths.csv (1-8) and
# parameter 1 = gamma1, 2 = lambda1, 3 = delta, 4 = epsilon.
#
#   GRID   base (default trange), dt0.45 or dt0.225 (SNTM_DT_MAX in days); the
#          truth is outputs_t60_<GRID>/ and the searches simulate on the same grid
#   TRUTH  numeric (the real test) or mlp (the MLP's own run: loss 0 at the truth,
#          tests the search alone)
#   START  mid (default): the middle of the training range, first step 5% of the
#          start coordinate; near: 0.98 x the true value, first step 0.005 decades
#          (gamma1, lambda1) or 1 (delta, epsilon), as in nm_failfast_tests.m
#
#   sbatch --array=1-32 run_nm_t60.slurm dt0.225 numeric
#   sbatch --array=1-32 run_nm_t60.slurm dt0.225 mlp near
#
# Writes outputs_nm/t60_<GRID>[_near]/<TRUTH>_truth/nm_l2f__mlp__<truth_id>__<param>.mat.
# Rerunning after a timeout resumes from the saved simulations.
#SBATCH --job-name=nm60
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out

set -eo pipefail
GRID=${1:?grid: base, dt0.45 or dt0.225}
TRUTH=${2:?truth: numeric or mlp}
START=${3:-mid}
PARAMS=(gamma1 lambda1 delta epsilon)
FIRST_STEPS=(0.005 0.005 1 1)
P=$(( (SLURM_ARRAY_TASK_ID - 1) % 4 ))
PARAM=${PARAMS[$P]}
TRUTH_ROW=$(( (SLURM_ARRAY_TASK_ID - 1) / 4 + 1 ))

VAULT=/home/bimstepf/endor/vault
export NTM_ML_ROOT=${NTM_ML_ROOT:-$VAULT/NTM_ML}
export NTM_SIM_ROOT=${NTM_SIM_ROOT:-$VAULT/Tau_Transport-jtorok_dev}
cd "$NTM_ML_ROOT/network_validation"

TRUTH_ID=$(awk -F, '$2 == "l2f" {print $1}' manifests/recovery_truths.csv | sed -n "${TRUTH_ROW}p")
[[ -n "$TRUTH_ID" ]] || { echo "no l2f truth in row $TRUTH_ROW" >&2; exit 1; }
export SNTM_T_END=60
case $GRID in
  base) unset SNTM_DT_MAX ;;
  dt0.45) export SNTM_DT_MAX=0.45 ;;
  dt0.225) export SNTM_DT_MAX=0.225 ;;
  *) echo "unknown grid $GRID" >&2; exit 1 ;;
esac
export NM_TRUTH_FILE=$PWD/outputs_t60_$GRID/rec_l2f__${TRUTH}__ECl__$TRUTH_ID.mat
SUFFIX=
case $START in
  mid) ;;
  near)
    # true value of the searched parameter: columns gamma1, lambda1, delta, epsilon = 3-6
    TRUE_VALUE=$(awk -F, -v id="$TRUTH_ID" -v c=$((P + 3)) '$1 == id && $2 == "l2f" {print $c}' manifests/recovery_truths.csv)
    export NM_START=$(awk -v x="$TRUE_VALUE" 'BEGIN {printf "%.10g", 0.98 * x}')
    export NM_STEP=${FIRST_STEPS[$P]}
    SUFFIX=_near ;;
  *) echo "unknown start $START" >&2; exit 1 ;;
esac
export NM_OUT_DIR=$PWD/outputs_nm/t60_$GRID$SUFFIX/${TRUTH}_truth
[[ -f "$NM_TRUTH_FILE" ]] || { echo "truth not found: $NM_TRUTH_FILE" >&2; exit 1; }

# MATLAB R2022b by path (see run_cases.slurm)
MATLAB_BIN=${MATLAB_BIN:-/opt/matlab/r2022b/bin/matlab}
[[ -x "$MATLAB_BIN" ]] || { echo "MATLAB not found at $MATLAB_BIN on $(hostname)" >&2; exit 1; }

echo "host=$(hostname) grid=$GRID truth=$TRUTH truth_id=$TRUTH_ID param=$PARAM dt_max=${SNTM_DT_MAX:-none} start=$START ${NM_START:+($NM_START, step $NM_STEP)}"
echo "truth file $NM_TRUTH_FILE"
"$MATLAB_BIN" -batch "maxNumCompThreads(${SLURM_CPUS_PER_TASK:-1}); nelder_mead_network_recovery('$TRUTH_ID', '$PARAM', 'mlp', 'l2f')"
