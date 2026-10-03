#!/bin/bash -l
# One segment of a chained, restartable single-GPU run. Do not sbatch this directly --
# perlmutter/submit_chain.sh supplies -A/-q/-t/--signal/--dependency and CHAIN_WORK.
#
# Each segment restarts from the newest checkpoint in $CHAIN_WORK/diags (none -> fresh
# start), runs until WarpX reaches max_step or until Slurm's SIGUSR1 (sent --signal
# seconds before walltime) makes it write a checkpoint and exit (warpx.break_signals,
# rendered by deck.render when diagnostics.checkpoint_intervals is set). The next
# segment, queued with afterany, picks up from there.
#
# Everything runs from a code snapshot in $CHAIN_WORK/.chain/code taken at submit time,
# so a multi-day chain does not depend on the checkout it was submitted from.
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -G 1

set -euo pipefail
shopt -s nullglob
: "${CHAIN_WORK:?CHAIN_WORK not set - submit via perlmutter/submit_chain.sh}"
CODE="$CHAIN_WORK/.chain/code"
export KINSHOCK_PM="$CODE/perlmutter"
# shellcheck disable=SC1090
source "$KINSHOCK_PM/_common.sh"

cd "$CHAIN_WORK"                       # diags/ lands here (see run_warpx in _common.sh)
LOG="$CHAIN_WORK/.chain/chain.log"
JOB="${SLURM_JOB_ID:-nojob}"
note() { echo "$(date +%Y-%m-%dT%H:%M:%S) job=$JOB $*" | tee -a "$LOG"; }

decks=(inputs_*); deck="${decks[0]}"
max_step="$(awk -F= '/^max_step/{gsub(/ /,"",$2); v=$2} END{print v}' "$deck")"
bin="$(pm_binary "${BINARY:-$SWEEP_BUILD}")"

last_step() {  # last completed step in a WarpX log (0 if none)
    grep -aoE "^STEP [0-9]+ ends" "$1" 2>/dev/null | tail -1 | awk '{print $2}' | grep . || echo 0
}
cancel_rest() {  # cancel the chain's later, still-queued segments
    local ids; ids="$(grep -v "^$JOB\$" .chain/jobids 2>/dev/null | tr '\n' ' ')"
    [[ -n "${ids// /}" ]] && scancel $ids 2>/dev/null || true
}
chks_by_step() {  # checkpoint dirs (not *.bad), oldest first
    local d
    for d in diags/chk[0-9]*; do
        [[ "$d" =~ /chk([0-9]+)$ ]] && echo "$((10#${BASH_REMATCH[1]})) $d"
    done | sort -n | cut -d' ' -f2
}
newest_chk() { chks_by_step | tail -1; }

# Already finished (an earlier segment reached max_step)?
for f in run_seg*.log; do
    [[ -f "$f" ]] || continue
    if [[ "$(last_step "$f")" -ge "$max_step" ]]; then
        note "run already complete (step $max_step in $f) - nothing to do"
        cancel_rest; exit 0
    fi
done

pm_profile
export MPICH_OFI_NIC_POLICY=GPU
export OMP_NUM_THREADS=16

# Up to two attempts: if the newest checkpoint is unreadable (a segment killed mid-write),
# mark it bad and fall back to the one before it.
rc=1
for attempt in 1 2; do
    chk="$(newest_chk)"
    segs=(run_seg*.log); seg=$(( ${#segs[@]} + 1 ))
    seglog="$(printf 'run_seg%02d.log' "$seg")"
    ln -sfn "$seglog" run.log          # the progress logger and --verify read run.log
    restart=()
    start=0
    if [[ -n "$chk" ]]; then
        restart=("amr.restart=$chk")
        start="${chk##*chk}"; start=$((10#$start))
    fi
    note "segment $seg start step=$start chk=${chk:-none} bin=$bin"

    logger_pid=""
    python3 "$CODE/scripts/run_progress_logger.py" "$CHAIN_WORK" >> .chain/logger.out 2>&1 &
    logger_pid=$!

    set +e
    srun --cpu-bind=cores "$bin" "$deck" "${restart[@]}" > "$seglog" 2>&1
    rc=$?
    set -e
    sleep 35; kill "$logger_pid" 2>/dev/null || true

    end="$(last_step "$seglog")"
    note "segment $seg exit rc=$rc end step=$end (max_step $max_step)"
    if [[ "$end" -le "$start" && -n "$chk" && "$attempt" -eq 1 ]]; then
        note "no step advanced from $chk - marking it bad and retrying from the previous one"
        mv "$chk" "$chk.bad"
        continue
    fi
    break
done

# Keep the newest three checkpoints; each holds every particle (~GBs).
chks_by_step | head -n -3 | while read -r old; do rm -rf "$old"; done

if [[ "$end" -ge "$max_step" ]]; then
    note "reached max_step - verifying deck against config and cancelling the rest of the chain"
    python3 "$CODE/scripts/make_inputs.py" "$CHAIN_WORK" --verify 2>&1 | tee -a "$LOG" || true
    cancel_rest
    exit 0
fi
if [[ "$end" -le "$start" ]]; then
    note "NO PROGRESS this segment - cancelling the rest of the chain; see $seglog"
    tail -15 "$seglog" >> "$LOG" || true
    cancel_rest
    exit 1
fi
exit 0
