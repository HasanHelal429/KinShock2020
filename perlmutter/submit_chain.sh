#!/bin/bash -l
# Submit a run as a chain of restartable single-GPU jobs (for runs longer than one QOS
# walltime). The run's config must set diagnostics.checkpoint_intervals, which makes the
# deck write periodic checkpoints and exit cleanly on SIGUSR1.
#
#   perlmutter/submit_chain.sh runs/IM_phase/im_470eV_t03 --segments 3
#   perlmutter/submit_chain.sh <run_dir> --dry                 # print, submit nothing
#   perlmutter/submit_chain.sh <run_dir> --resume --segments 1 # extend an existing chain
#
# Options: --segments N (3)  --time HH:MM:SS (48:00:00)  --qos Q (shared)
#          --signal-lead S (600: seconds before walltime that SIGUSR1 is sent)
#          --work DIR ($PSCRATCH/kinshock_runs/<phase>/<id>)  --binary A|B ($SWEEP_BUILD)
#
# The run executes in WORK, not in the repo: a multi-day chain must not live inside a
# checkout (or worktree) that can move under it. WORK gets the deck, config.yaml, the
# run README and a snapshot of src/ scripts/ perlmutter/ in WORK/.chain/code, which is
# what every segment runs -- and the provenance of the code that rendered the deck.
set -euo pipefail
KINSHOCK_PM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export KINSHOCK_PM
REPO="$(cd "$KINSHOCK_PM/.." && pwd)"
# shellcheck disable=SC1090
source "$KINSHOCK_PM/_common.sh"

RUN="${1:?usage: $0 <run_dir> [--segments N] [--time T] [--qos Q] [--work DIR] [--resume] [--dry]}"; shift
SEGS=3; WALL=48:00:00; QOS=shared; LEAD=600; WORK=""; RESUME=0; DRY=0
BIN="$SWEEP_BUILD"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --segments)    SEGS="$2"; shift 2 ;;
        --time)        WALL="$2"; shift 2 ;;
        --qos)         QOS="$2"; shift 2 ;;
        --signal-lead) LEAD="$2"; shift 2 ;;
        --work)        WORK="$2"; shift 2 ;;
        --binary)      BIN="$2"; shift 2 ;;
        --resume)      RESUME=1; shift ;;
        --dry)         DRY=1; shift ;;
        *) echo "unknown option '$1'" >&2; exit 2 ;;
    esac
done

SRC="$(cd "$REPO/$RUN" 2>/dev/null && pwd)" || { echo "submit_chain: no run dir '$RUN' under $REPO" >&2; exit 2; }
[[ -f "$SRC/config.yaml" ]] || { echo "submit_chain: no config.yaml in $SRC" >&2; exit 2; }
DECK="$(ls "$SRC"/inputs_* | head -1)"
grep -q "^chk.format *= *checkpoint" "$DECK" && grep -q "^warpx.break_signals" "$DECK" || {
    echo "submit_chain: $DECK has no checkpoint diag / break signal -- set" >&2
    echo "              diagnostics.checkpoint_intervals in config.yaml and regenerate" >&2; exit 2; }
pm_profile > /dev/null 2>&1 || true    # the system python3 is 3.6; the profile's venv is not
python3 "$REPO/scripts/make_inputs.py" "$SRC" --check > /dev/null || {
    echo "submit_chain: $DECK is stale against config.yaml -- regenerate it first" >&2; exit 2; }

PHASE="$(basename "$(dirname "$SRC")")"; ID="$(basename "$SRC")"
WORK="${WORK:-${PSCRATCH:?}/kinshock_runs/$PHASE/$ID}"
if [[ $RESUME -eq 0 ]] && compgen -G "$WORK/diags/*" > /dev/null; then
    echo "submit_chain: $WORK/diags already has output. --resume to extend that chain," >&2
    echo "              or move it aside for a fresh start." >&2; exit 2
fi
if [[ $RESUME -eq 1 ]]; then
    cmp -s "$DECK" "$WORK/$(basename "$DECK")" || {
        echo "submit_chain: --resume but the repo deck differs from $WORK's -- refusing to" >&2
        echo "              continue a run under a different deck" >&2; exit 2; }
fi

echo "run   $SRC"
echo "work  $WORK"
echo "chain $SEGS x $WALL on -q $QOS, SIGUSR1 ${LEAD}s before walltime, binary $BIN ($(pm_binary "$BIN"))"
if [[ $DRY -eq 0 ]]; then
    mkdir -p "$WORK/.chain"
    if [[ $RESUME -eq 0 ]]; then
        cp "$DECK" "$SRC/config.yaml" "$WORK/"
        [[ -f "$SRC/README.md" ]] && cp "$SRC/README.md" "$WORK/"
        rm -rf "$WORK/.chain/code"; mkdir -p "$WORK/.chain/code/perlmutter"
        rsync -a --exclude __pycache__ "$REPO/src" "$REPO/scripts" "$WORK/.chain/code/"
        cp "$KINSHOCK_PM"/{_common.sh,site.conf,job_chain.sbatch} "$WORK/.chain/code/perlmutter/"
        git -C "$REPO" rev-parse HEAD > "$WORK/.chain/code/COMMIT" 2>/dev/null || true
        git -C "$REPO" status --porcelain >> "$WORK/.chain/code/COMMIT" 2>/dev/null || true
        : > "$WORK/.chain/jobids"
    fi
fi

prev="$(tail -1 "$WORK/.chain/jobids" 2>/dev/null || true)"
for i in $(seq 1 "$SEGS"); do
    dep=(); [[ -n "$prev" ]] && dep=(--dependency="afterany:$prev")
    cmd=(sbatch --parsable -A "$NERSC_ACCOUNT" -q "$QOS" -t "$WALL" -J "c_$ID"
         --signal="USR1@$LEAD" "${dep[@]}"
         -o "$WORK/.chain/slurm-%j.out" -e "$WORK/.chain/slurm-%j.out"
         --export="ALL,CHAIN_WORK=$WORK,BINARY=$BIN"
         "$WORK/.chain/code/perlmutter/job_chain.sbatch")
    echo "+ ${cmd[*]}"
    if [[ $DRY -eq 1 ]]; then prev="DRY$i"; continue; fi
    prev="$("${cmd[@]}")"
    echo "$prev" >> "$WORK/.chain/jobids"
    echo "  -> job $prev"
done
[[ $DRY -eq 1 ]] && echo "(--dry: nothing submitted)"
exit 0
