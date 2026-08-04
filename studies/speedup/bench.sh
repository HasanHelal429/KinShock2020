#!/bin/bash
# studies/speedup/bench.sh -- time a SHORT WarpX run and report tail-averaged s/step.
#
# WHY NOT scripts/launch.sh. launch.sh cd's into the run dir so diags land there, which
# is exactly right for a physics run and exactly wrong for a benchmark: seven benchmark
# points would all write into runs/<ID>/diags and clobber the staged production run.
# This script cd's into its own scratch dir under studies/speedup/out/<label> instead,
# and suppresses plotfiles (they are 277 MiB each here, and disk is at 94%).
#
# It mirrors launch.sh's thread environment exactly (OMP_PROC_BIND=spread,
# OMP_PLACES=cores) so numbers are comparable to the R1_paper_470eV_pilot baseline.
#
# THE REPORTED NUMBER is the mean of WarpX's own per-step "This step = X s", with the
# first 20% of steps dropped. The warm-up matters: the pilot's step 1 took 0.241 s
# against an asymptote of 0.111, and its rate was still moving at step 2500.
#
#   bench.sh <label> <run_dir> <threads> <warpx_bin> <max_step> <cap_s> [parmparse...]
#
# Exits 0 even if WarpX aborts or hits the cap -- a failed point is a RESULT (that is
# the whole question for the implicit lever), so the caller decides what it means.
# `status` in the summary line is one of ok / capped / aborted.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [[ $# -lt 6 ]]; then
    sed -n '2,22p' "${BASH_SOURCE[0]}" >&2
    exit 2
fi

LABEL="$1"; RUN_DIR="$2"; THREADS="$3"; BIN="$4"; MAX_STEP="$5"; CAP="$6"; shift 6
EXTRA=("$@")

RUN_DIR="$(cd "$RUN_DIR" && pwd)"
[[ -x "$BIN" ]] || { echo "bench: not executable: $BIN" >&2; exit 2; }

shopt -s nullglob
DECKS=("$RUN_DIR"/inputs_*)
shopt -u nullglob
[[ ${#DECKS[@]} -eq 1 ]] || { echo "bench: want exactly one deck in $RUN_DIR, got ${#DECKS[@]}" >&2; exit 2; }
DECK="${DECKS[0]}"

OUT="$ROOT/studies/speedup/out/$LABEL"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"

export OMP_NUM_THREADS="$THREADS" OMP_PROC_BIND=spread OMP_PLACES=cores

{
    echo "label      $LABEL"
    echo "deck       $DECK"
    echo "binary     $BIN"
    echo "threads    $THREADS"
    echo "max_step   $MAX_STEP"
    echo "cap_s      $CAP"
    echo "extra      ${EXTRA[*]:-(none)}"
    echo "started    $(date -Is)"
    echo "load_before $(cut -d' ' -f1 /proc/loadavg)"
} > meta.txt

# diag*.intervals=0 kills the 277 MiB plotfiles; the EP/PN reduced diags stay (they
# cost 11 s per 50k steps and are what the GPU-vs-CPU agreement check compares).
T0=$(date +%s.%N)
timeout --signal=INT "$CAP" \
    "$BIN" "$DECK" \
    max_step="$MAX_STEP" \
    diag1.intervals=0 \
    diag_fields.intervals=0 \
    ${EXTRA[@]+"${EXTRA[@]}"} > run.log 2>&1
RC=$?
T1=$(date +%s.%N)

STATUS=ok
[[ $RC -eq 124 || $RC -eq 130 ]] && STATUS=capped
[[ $RC -ne 0 && $STATUS == ok ]] && STATUS=aborted

python3 - "$OUT" "$LABEL" "$THREADS" "$STATUS" "$RC" "$T0" "$T1" <<'PY'
import re, sys, statistics as st
out, label, threads, status, rc, t0, t1 = sys.argv[1:]
steps = [float(m) for m in re.findall(r"This step = ([0-9.eE+-]+) s", open(f"{out}/run.log", errors="replace").read())]
wall = float(t1) - float(t0)
if not steps:
    line = f"{label:34s} thr={threads:>3s}  n=0  NO STEPS  status={status} rc={rc} wall={wall:.0f}s"
else:
    tail = steps[len(steps) // 5:] or steps          # drop the first 20% (warm-up)
    line = (f"{label:34s} thr={threads:>3s}  n={len(steps):6d}  "
            f"s/step={st.mean(tail):8.5f}  median={st.median(tail):8.5f}  "
            f"status={status} wall={wall:.0f}s")
open(f"{out}/result.txt", "w").write(line + "\n")
print(line)
PY
exit 0
