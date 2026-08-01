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
#   scripts/launch.sh runs/R0 -- max_step=20
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
THREADS=8
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
[[ -x "$WARPX" ]] || die "WarpX binary not executable: $WARPX (set --warpx or \$KINSHOCK_WARPX)"

# Exactly one deck, so we never guess which input file was meant.
shopt -s nullglob
DECKS=("$RUN_DIR"/inputs_*)
shopt -u nullglob
case ${#DECKS[@]} in
    0) die "no deck in $RUN_DIR -- run: python scripts/make_inputs.py $RUN_DIR" ;;
    1) DECK="$(basename "${DECKS[0]}")" ;;
    *) die "$RUN_DIR has ${#DECKS[@]} decks (${DECKS[*]##*/}) -- keep one" ;;
esac

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
echo "launch: $WARPX $DECK ${EXTRA[*]:-} > run.log 2>&1"
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
    nohup "$WARPX" "$DECK" ${EXTRA[@]+"${EXTRA[@]}"} > run.log 2>&1 &
    echo "launch: warpx pid $! -> $(basename "$RUN_DIR")/run.log"
    start_logger
    echo "launch: tail -f $RUN_DIR/run.log"
else
    start_logger
    echo "launch: running in the foreground; tail -f $RUN_DIR/run.log"
    exec "$WARPX" "$DECK" ${EXTRA[@]+"${EXTRA[@]}"} > run.log 2>&1
fi
