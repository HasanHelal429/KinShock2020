#!/bin/bash -l
# Generic single-GPU WarpX job body. ONE ARRAY TASK = ONE RUN.
#
# Do not sbatch this directly -- perlmutter/submit.sh supplies -A/-q/-t/--array and the
# RUNS/BINARY/WORKROOT environment, all of which come from perlmutter/site.conf.
#
# Resources are per-task because these runs are single-GPU by construction: every
# S_phase config sets numerics.max_grid_size = n_cell, i.e. ONE box, and one box cannot
# be split across ranks. Asking for more GPUs would leave them idle.
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -G 1
#SBATCH -J kinshock
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err

set -euo pipefail

# shellcheck disable=SC1090
source "${KINSHOCK_PM}/_common.sh"

# One run per line; the array index selects the line. A file rather than an exported
# variable because sbatch --export is comma-separated and mangles anything containing
# spaces -- see the note in submit.sh.
IDX="${SLURM_ARRAY_TASK_ID:-0}"
SPEC="$(sed -n "$((IDX + 1))p" "${RUNLIST:?RUNLIST not set - submit via perlmutter/submit.sh}")"
[[ -n "$SPEC" ]] || { echo "job: no run at index $IDX in $RUNLIST" >&2; exit 2; }

# Each entry is "<run_dir_relative_to_repo>[:<label>]". A label means "run this config in
# a scratch working directory under that name" -- how the A/B test gets its replicate runs
# without a second copy of an identical config, which would be a lie about what varies.
RUN_REL="${SPEC%%:*}"
LABEL="${SPEC#*:}"
[[ "$LABEL" == "$SPEC" ]] && LABEL=""

run_warpx "$RUN_REL" "$LABEL"
