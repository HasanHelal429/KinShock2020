#!/bin/bash
# scripts/launch.sh -- the ONE correct way to start a WarpX run in this repo.
#
# The generated deck sets no `diag*.file_prefix`, so WarpX writes plotfiles to `diags/`
# RELATIVE TO THE LAUNCH CWD. Launching two runs from the repo root makes them share
# ./diags/ and clobber each other (WarpX leaves .old.NNNN rename files as the tell) --
# that cost a rerun on the R2/R3 controls. This script always cd's into the run dir
# first, so diag1/diag_fields/reducedfiles land under runs/<ID>/diags/.
#
# Also applies the benchmarked thread settings (RESULTS 2026-07-23: near-linear to 8
# cores, ~1.8x vs 4; max_grid_size/tiling/sort_intervals were neutral-to-negative).
#
# Usage:  scripts/launch.sh [options] <run_dir> [-- <warpx args>]
#
#   -j, --threads N    OMP_NUM_THREADS (default 8)
#   -g, --gpu [LIST]   run on GPU(s) LIST (default 0; e.g. "0" or "0,1"): selects the CUDA
#                      build, 1 thread, pins CUDA_VISIBLE_DEVICES so unlisted cards stay
#                      free for other users, and runs ONE MPI RANK PER DEVICE via mpirun.
#                      Two guards, both measured (RESULTS 2026-08-04): refuses without
#                      numerics.max_grid_size, since AMReX's default decomposition costs
#                      6.5x on GPU (1.21x vs 7.89x over 8 CPU threads); and refuses if the
#                      decomposition gives fewer boxes than ranks, since the extra ranks
#                      would idle and the run would just look half as fast as expected.
#                      Two cards measured 1.77x (91% parallel efficiency) on R1_paper_470eV.
#   -w, --warpx PATH   WarpX binary (default $KINSHOCK_WARPX, else the repo's usual build)
#   -b, --background   detach and return immediately (prints the PID)
#   -L, --logger       start scripts/run_progress_logger.py (DEFAULT ON)
#       --no-logger    do NOT start the progress logger (rarely what you want:
#                      a run with no progress.log leaves no wall-clock record)
#       --every-pct P  logger checkpoint every P percent of max_step (default 10)
#       --poll S       logger poll interval in seconds (default 30)
#   -f, --force        launch even though diags/ already holds output (see below)
#   -n, --dry-run      print what would run, change nothing
#
# Anything after `--` is appended as ParmParse overrides, e.g.
#   scripts/launch.sh runs/R0_phase/R0 -- max_step=20
# NOTE overrides are not reflected in config.yaml, so `make_inputs.py --verify` will
# flag them afterwards. Use them for smoke tests, not for physics.
#
# Stdout/stderr go to <run_dir>/run.log (gitignored; findings belong in RESULTS.md).
#
# Refuses to start when diags/ is already populated, because relaunching overwrites a
# finished run's plotfiles in place -- pass --force once you mean it, or move the old
# diags/ aside. Run `make_inputs.py <run_dir> --check` first if the deck may be stale.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WARPX="${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
WARPX_CUDA="${KINSHOCK_WARPX_CUDA:-/home/hhelal/warpx-cda/build_cuda1d/bin/warpx.1d}"
THREADS=8
GPU=""          # empty = CPU; otherwise the device index
BACKGROUND=0
LOGGER=1        # progress logger is ON by default -- every run gets a progress.log
EVERY_PCT=""
POLL=""
FORCE=0
DRYRUN=0
RUN_DIR=""
EXTRA=()

die() { echo "launch: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --)              shift; EXTRA=("$@"); break ;;
        -j|--threads)    THREADS="${2:-}"; shift 2 ;;
        -g|--gpu)        # optional device LIST; bare -g means device 0. "0,1" runs one MPI
                         # rank per device via mpirun.
                         if [[ "${2:-}" =~ ^[0-9]+(,[0-9]+)*$ ]]; then GPU="$2"; shift 2
                         else GPU=0; shift; fi ;;
        -w|--warpx)      WARPX="${2:-}";   shift 2 ;;
        -b|--background) BACKGROUND=1; shift ;;
        -L|--logger)     LOGGER=1;     shift ;;
        --no-logger)     LOGGER=0;     shift ;;
        --every-pct)     EVERY_PCT="$2"; shift 2 ;;
        --poll)          POLL="$2";      shift 2 ;;
        -f|--force)      FORCE=1;      shift ;;
        -n|--dry-run)    DRYRUN=1;     shift ;;
        -h|--help)       sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
        -*)              die "unknown option '$1' (try --help)" ;;
        *)               [[ -n "$RUN_DIR" ]] && die "one run_dir at a time (got '$RUN_DIR' and '$1')"
                         RUN_DIR="$1"; shift ;;
    esac
done

[[ -n "$RUN_DIR" ]] || die "usage: scripts/launch.sh [options] <run_dir>"
[[ -d "$RUN_DIR" ]] || die "no such run dir: $RUN_DIR"
RUN_DIR="$(cd "$RUN_DIR" && pwd)"                       # absolute: we are about to cd
[[ -f "$RUN_DIR/config.yaml" ]] || die "$RUN_DIR has no config.yaml -- is it a run dir?"

# Exactly one deck, so we never guess which input file was meant.
shopt -s nullglob
DECKS=("$RUN_DIR"/inputs_*)
shopt -u nullglob
case ${#DECKS[@]} in
    0) die "no deck in $RUN_DIR -- run: python scripts/make_inputs.py $RUN_DIR" ;;
    1) DECK="$(basename "${DECKS[0]}")" ;;
    *) die "$RUN_DIR has ${#DECKS[@]} decks (${DECKS[*]##*/}) -- keep one" ;;
esac

# --- GPU mode -----------------------------------------------------------------------
# Selects the CUDA build, drops to one thread (OMP is irrelevant on device), pins
# CUDA_VISIBLE_DEVICES so unlisted cards stay free for other users on this shared box, and
# runs one MPI rank per listed device. Placed after deck discovery because the box-vs-rank
# check reads amr.n_cell / amr.max_grid_size out of the DECK.
MPI_PREFIX=()
if [[ -n "$GPU" ]]; then
    [[ "$WARPX" == "${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}" ]] \
        && WARPX="$WARPX_CUDA"          # -w wins if given explicitly
    THREADS=1
    export CUDA_VISIBLE_DEVICES="$GPU"
    NRANKS=$(awk -F, '{print NF}' <<< "$GPU")     # one rank per device

    # The single biggest GPU footgun: without amr.max_grid_size AMReX picks a CPU-friendly
    # ~235-box decomposition and the GPU runs 6.5x slower -- 1.21x over 8 CPU threads
    # instead of 7.89x (RESULTS 2026-08-04). A silent multi-day mistake, so refuse.
    grep -qE '^\s*max_grid_size\s*:' "$RUN_DIR/config.yaml" \
        || die "GPU mode needs numerics.max_grid_size in $RUN_DIR/config.yaml
     (n_cell for one GPU, n_cell/2 for two). Without it AMReX's default decomposition costs
     6.5x on GPU. Add it to config.yaml and regenerate the deck, not as an override."

    if [[ $NRANKS -gt 1 ]]; then
        command -v mpirun >/dev/null || die "-g $GPU needs mpirun for $NRANKS ranks"
        MPI_PREFIX=(mpirun -np "$NRANKS")
        # Boxes must be >= ranks or a rank owns nothing and idles -- which looks like a
        # working run at half the expected speed. Read the DECK, so this reflects what
        # WarpX will actually be given rather than what the config meant.
        NC=$(awk -F= '/^amr\.n_cell/{gsub(/ /,"",$2); print $2}'        "${DECKS[0]}")
        MGS=$(awk -F= '/^amr\.max_grid_size/{gsub(/ /,"",$2); print $2}' "${DECKS[0]}")
        NBOX=$(( (NC + MGS - 1) / MGS ))
        [[ $NBOX -ge $NRANKS ]] || die "decomposition gives $NBOX box(es) for $NRANKS ranks
     (amr.n_cell=$NC, amr.max_grid_size=$MGS). Ranks beyond the first would idle. Set
     numerics.max_grid_size to about n_cell/$NRANKS and regenerate."
        echo "launch: $NRANKS MPI ranks, $NBOX boxes (n_cell=$NC, max_grid_size=$MGS)"
    fi
    if command -v nvidia-smi >/dev/null; then
        nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv,noheader \
                   --id="$GPU" | sed 's/^/launch: gpu /'
    fi
fi
# After the GPU block, so it validates the binary actually about to run.
[[ -x "$WARPX" ]] || die "WarpX binary not executable: $WARPX (set --warpx or \$KINSHOCK_WARPX)"

if [[ -d "$RUN_DIR/diags" ]] && compgen -G "$RUN_DIR/diags/*" >/dev/null; then
    if [[ $FORCE -eq 1 ]]; then
        echo "launch: --force, writing into the existing $RUN_DIR/diags"
    elif [[ $DRYRUN -eq 1 ]]; then
        echo "launch: WARNING $RUN_DIR/diags already has output -- a real launch would"
        echo "launch:         refuse this without --force."
    else
        die "$RUN_DIR/diags already has output -- relaunching overwrites it in place.
     Move it aside, or pass --force if that is what you want."
    fi
fi

echo "launch: $(basename "$RUN_DIR")  deck=$DECK  threads=$THREADS"
echo "launch: cwd=$RUN_DIR  (so diags/ lands here, not in the repo root)"
echo "launch: ${MPI_PREFIX[*]:-}${MPI_PREFIX:+ }$WARPX $DECK ${EXTRA[*]:-} > run.log 2>&1"
if [[ $DRYRUN -eq 1 ]]; then
    if [[ $LOGGER -eq 1 ]]; then
        echo "launch: progress logger WOULD start -> $(basename "$RUN_DIR")/progress.log (every=${EVERY_PCT:-10}%% poll=${POLL:-30}s)"
    else
        echo "launch: progress logger DISABLED (--no-logger) — run leaves no wall-clock record"
    fi
    echo "launch: --dry-run, nothing started."
    exit 0
fi

cd "$RUN_DIR"                                           # THE POINT OF THIS SCRIPT
export OMP_NUM_THREADS="$THREADS" OMP_PROC_BIND=spread OMP_PLACES=cores

start_logger() {   # after WarpX, so run.log exists (the logger waits for it anyway)
    [[ $LOGGER -eq 1 ]] || return 0
    local largs=()
    [[ -n "$EVERY_PCT" ]] && largs+=(--every-pct "$EVERY_PCT")
    [[ -n "$POLL"      ]] && largs+=(--poll "$POLL")
    nohup python "$ROOT/scripts/run_progress_logger.py" "$RUN_DIR" \
        ${largs[@]+"${largs[@]}"} > "$RUN_DIR/logger.out" 2>&1 &
    echo "launch: progress logger pid $! -> $(basename "$RUN_DIR")/progress.log"
}

if [[ $BACKGROUND -eq 1 ]]; then
    nohup ${MPI_PREFIX[@]+"${MPI_PREFIX[@]}"} "$WARPX" "$DECK" ${EXTRA[@]+"${EXTRA[@]}"} > run.log 2>&1 &
    echo "launch: warpx pid $! -> $(basename "$RUN_DIR")/run.log"
    start_logger
    echo "launch: tail -f $RUN_DIR/run.log"
else
    start_logger
    echo "launch: running in the foreground; tail -f $RUN_DIR/run.log"
    exec ${MPI_PREFIX[@]+"${MPI_PREFIX[@]}"} "$WARPX" "$DECK" ${EXTRA[@]+"${EXTRA[@]}"} > run.log 2>&1
fi
