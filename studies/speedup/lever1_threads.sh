#!/bin/bash
# Lever 1: OMP thread count.
#
# The pilot ran at 8 threads (run.log:12) on a 32-core box. RESULTS 2026-07-23 measured
# "near-linear to 8 cores, ~1.8x vs 4" but never went ABOVE 8, and the standing note is
# that WarpX collapses ~20x somewhere above 12. This scan finds the actual knee.
#
# Tile supply is not the limit: the deck runs 235 grids over 30000 cells (112-128 cells
# each, run.log:28), so even 24 threads gets ~10 boxes apiece. The two max_grid_size
# variants test whether box GRANULARITY rather than count is what bites at high counts.
#
# Every point is capped at 400 s, so a collapsed high-thread point costs 400 s and not
# an hour. A capped point still reports a valid s/step from the steps it did finish.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="${1:-$HERE/../../runs/R1_paper_470eV}"
BIN="${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
STEPS=1500
CAP=400

echo "=== Lever 1: OMP thread scan ==="
echo "deck $RUN   binary $BIN   steps $STEPS   cap ${CAP}s/point"
echo "baseline for reference: 0.1113 s/step at 8 threads (R1_paper_470eV_pilot)"
echo

for T in 4 8 12 16 20 24; do
    "$HERE/bench.sh" "l1_thr${T}" "$RUN" "$T" "$BIN" "$STEPS" "$CAP"
done

# Granularity variants at the two thread counts most likely to be limited by it.
"$HERE/bench.sh" "l1_thr16_mgs64" "$RUN" 16 "$BIN" "$STEPS" "$CAP" amr.max_grid_size=64
"$HERE/bench.sh" "l1_thr24_mgs64" "$RUN" 24 "$BIN" "$STEPS" "$CAP" amr.max_grid_size=64

echo
echo "=== Lever 1 summary ==="
cat "$HERE"/out/l1_*/result.txt
