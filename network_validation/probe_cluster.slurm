#!/bin/bash -l
# Diagnostic job: reports what a compute node offers for running the NTM
# simulations (hardware, account, modules, MATLAB, toolboxes, licence).
# Submit once per node and send back the probe_*.out files:
#
#   sbatch probe_cluster.slurm                                   # default partition, wherever SLURM puts it
#   for n in node01 node02 node03 node04 node05; do sbatch -p oldnodes -w $n probe_cluster.slurm; done
#   sbatch -p himem probe_cluster.slurm
#
#SBATCH --job-name=probe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=probe_%j_%N.out

section() { printf '\n==================== %s ====================\n' "$1"; }
VAULT=/home/bimstepf/endor/vault

section "NODE"
echo "host: $(hostname)   partition: $SLURM_JOB_PARTITION   job: $SLURM_JOB_ID"
grep PRETTY_NAME /etc/os-release; uname -r
lscpu | grep -E '^(Model name|CPU\(s\)|Thread|Core|Socket)'
free -g | head -2
df -h /home "$VAULT" /tmp 2>&1 | sort -u

section "ACCOUNT AND FILES"
id 2>&1
echo "HOME=$HOME"
ls -d "$VAULT"/* 2>&1
( cd "$VAULT/NTM_ML" 2>&1 && git log -1 --oneline 2>&1 && ls network_validation 2>&1 )
touch "$VAULT/NTM_ML/network_validation/.write_test" 2>&1 && echo "write access: yes" && rm -f "$VAULT/NTM_ML/network_validation/.write_test"

section "SLURM ENVIRONMENT"
env | grep -E '^SLURM_(CPUS|MEM|JOB_PARTITION|NODELIST|SUBMIT_HOST)' | sort

section "MODULE SYSTEM"
type module 2>&1 | head -1
echo "MODULEPATH=$MODULEPATH"
module avail 2>&1 | head -40
echo "--- module load matlab (inherited environment):"
module load matlab 2>&1; command -v matlab || echo "matlab not on PATH"
echo "--- module load matlab (fresh login shell, clean environment):"
env -i HOME="$HOME" bash -lc 'echo "MODULEPATH=$MODULEPATH"; module load matlab 2>&1; command -v matlab || echo "matlab not on PATH"'

section "SEARCH FOR MATLAB"
CANDS=$( { command -v matlab
           env -i HOME="$HOME" bash -lc 'module load matlab >/dev/null 2>&1; command -v matlab'
           ls -d /opt/*/bin/matlab /opt/*/*/bin/matlab /opt/*/*/*/bin/matlab \
                 /usr/local/*/bin/matlab /usr/local/*/*/bin/matlab /usr/local/*/*/*/bin/matlab \
                 /share/*/*/bin/matlab /apps/*/*/bin/matlab /software/*/*/bin/matlab \
                 /cm/shared/apps/*/*/bin/matlab 2>/dev/null
           find /opt /usr/local /share /apps /software /cm -maxdepth 6 -type f -path '*/bin/matlab' 2>/dev/null
         } | xargs -r -n1 readlink -f 2>/dev/null | sort -u )
echo "${CANDS:-no matlab executable found}"
echo "MLM_LICENSE_FILE=$MLM_LICENSE_FILE  LM_LICENSE_FILE=$LM_LICENSE_FILE"

for M in $CANDS; do
    ROOT=$(dirname "$(dirname "$M")")
    section "MATLAB AT $ROOT"
    cat "$ROOT/VersionInfo.xml" 2>/dev/null | grep -o '<version>.*</version>\|<release>.*</release>'
    echo "--- toolbox folders:"; ls "$ROOT/toolbox" 2>&1 | tr '\n' ' '; echo
    echo "--- licence files:"; ls -l "$ROOT/licenses" 2>&1
    if grep -qE '<release>R20(19|2[0-9])' "$ROOT/VersionInfo.xml" 2>/dev/null; then
        echo "--- launching MATLAB (licence checkout and the functions the NTM needs):"
        ( cd "$VAULT/NTM_ML/network_validation" && timeout 900 "$M" -batch "probe_matlab" 2>&1 )
        echo "MATLAB exit code: $?"
    else
        echo "--- skipping launch (release older than R2019a has no -batch)"
    fi
done

section "OTHER"
command -v python3 && python3 --version
echo "done: $(date)"
