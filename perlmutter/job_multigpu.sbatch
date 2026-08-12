#!/bin/bash -l
# Single multi-GPU WarpX run. ONE JOB = ONE RUN (no array).
#
# WHY THIS EXISTS ALONGSIDE job.sbatch. job.sbatch is hard-wired `-n 1 -G 1` because every
# S_phase config sets max_grid_size = n_cell, i.e. ONE box, and one box cannot be split
# across ranks. The production R1_paper_* decks decompose into several boxes and need one
# rank per box: boxes must be >= ranks or a rank sits idle, the same rule
# scripts/launch.sh enforces on chablis.
#
# RANKS AND GPUS COME FROM THE SUBMIT LINE. The #SBATCH values below are only a default
# (4, matching a full Perlmutter GPU node); sbatch command-line flags override in-file
# directives, so `sbatch -n 2 -G 2 ...` gives a 2-GPU run with no edit here. **Match the
# rank count to the deck's box count** -- amr.max_grid_size = n_cell/<boxes>.
#
#   4 boxes (mgs = n_cell/4) -> -n 4 -G 4, and QOS `regular`: gpu_shared is capped at
#                               gres/gpu=2, so anything above 2 GPUs needs a whole node.
#   2 boxes (mgs = n_cell/2) -> -n 2 -G 2, QOS `shared` (billed per GPU).
#
# Everything else is deliberately NOT duplicated: this sources _common.sh and calls
# run_warpx, so the cd-before-launch invariant, the progress logger, the refuse-to-
# overwrite-diags guard and the post-run --verify all behave exactly as they do for the
# single-GPU path.
#
# Submit with:
#   sbatch -A <acct> -q regular -t 48:00:00 -n 4 -G 4 \
#          --export=ALL,KINSHOCK_PM=<repo>/perlmutter,RUN=runs/R1_phase/<ID> \
#          perlmutter/job_multigpu.sbatch
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 32
#SBATCH -G 4
#SBATCH -J kinshockmg
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

set -euo pipefail

# shellcheck disable=SC1090
source "${KINSHOCK_PM}/_common.sh"

: "${RUN:?RUN not set -- pass RUN=runs/<phase>_phase/<ID> via sbatch --export}"

echo "job_multigpu: ranks=${SLURM_NTASKS:-?} gpus=${SLURM_GPUS_ON_NODE:-?} run=$RUN"
run_warpx "$RUN"
