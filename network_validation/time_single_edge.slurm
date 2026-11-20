#!/bin/bash -l
# Time the numeric single-edge data generator on a cluster node (time_single_edge.m):
# the first N samples of NetworkTransportModel_SingleEdge_AllParams.m, in a parfor
# on as many workers as CPUs requested (the generator itself uses 4), extrapolated
# to 10^6 samples. The result is printed at the end of the log and written to
# results/time_single_edge_<host>_<workers>w.csv.
#   sbatch time_single_edge.slurm            # 400 samples on 4 workers
#   sbatch --cpus-per-task=8 time_single_edge.slurm 800
#SBATCH --job-name=time_se
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out

set -eo pipefail
N=${1:-400}

VAULT=/home/bimstepf/endor/vault
export NTM_ML_ROOT=$VAULT/NTM_ML
export NTM_SIM_ROOT=${NTM_SIM_ROOT:-$VAULT/Tau_Transport-jtorok_dev}

MATLAB_BIN=${MATLAB_BIN:-/opt/matlab/r2022b/bin/matlab}
[[ -x "$MATLAB_BIN" ]] || { echo "MATLAB not found at $MATLAB_BIN on $(hostname)" >&2; exit 1; }

cd "$NTM_ML_ROOT/network_validation"
echo "host=$(hostname) cpus=$SLURM_CPUS_PER_TASK samples=$N"
grep -m1 "model name" /proc/cpuinfo || true

"$MATLAB_BIN" -batch "time_single_edge($N, $SLURM_CPUS_PER_TASK)"
