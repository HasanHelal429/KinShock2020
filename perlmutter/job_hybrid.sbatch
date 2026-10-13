#!/bin/bash -l
# One HYBRID run (H-PICShock deck: kinetic ions, fluid electrons) on one GPU.
#
#   sbatch -A m5032_g -q debug -t 00:30:00 --export=ALL,RUN_DIR=<abs run dir> \
#          perlmutter/job_hybrid.sbatch
#
# H-PICShock (github.com/HasanHelal429/H-PICShock) generates these decks from its own
# config schema, which kinshock cannot load, so submit_chain/run_warpx do not apply. This
# keeps their launch invariants: cd into the run dir (diags/ lands there), the profile with
# the CUDA 12 shim (_common.sh: pm_profile), the commit-A binary, output in run.log, and a
# refusal to overwrite existing diags.
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -G 1
#SBATCH -J hybrid
#SBATCH -o slurm-hybrid-%j.out

set -euo pipefail
KINSHOCK_PM="${KINSHOCK_PM:-/pscratch/sd/h/hhelal/KinShock2020/perlmutter}"
export KINSHOCK_PM
# shellcheck disable=SC1090
source "$KINSHOCK_PM/_common.sh"
: "${RUN_DIR:?set RUN_DIR to the run directory}"
cd "$RUN_DIR"
if compgen -G "diags/*" > /dev/null; then
    echo "job_hybrid: $RUN_DIR/diags already has output -- refusing to overwrite" >&2; exit 1
fi
decks=(inputs_*); deck="${decks[0]}"
bin="$(pm_binary "${BINARY:-A}")"
pm_profile
export OMP_NUM_THREADS=16 MPICH_OFI_NIC_POLICY=GPU
echo "job_hybrid: $(basename "$RUN_DIR") deck=$deck bin=$bin"
set +e
srun --cpu-bind=cores "$bin" "$deck" > run.log 2>&1
rc=$?
set -e
echo "job_hybrid: exit $rc"
tail -3 run.log
exit $rc
