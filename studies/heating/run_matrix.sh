#!/bin/bash
# Grid-heating mitigation pilot matrix -- runs/heat_phase/, five variants, SERIALLY.
#
# WHY SERIAL. Each variant uses BOTH GPUs (one MPI rank per card, amr.max_grid_size = 15000).
# Two variants at once would halve each other's cards and make the timings meaningless --
# and timing is half the point, since the cost of shape 3 is an estimate, not a measurement.
#
# WHAT IT ANSWERS. R1_paper_470eV heats the far upstream 10 -> 115 eV over 141 t_ab, taking
# beta_0 from 0.196 to 2.25 against Table I's 0.2 (RESULTS 2026-08-05). Holding beta_0 under
# 1.0 needs ~2.6x less heating; under 0.40 needs ~10x. Three of the levers have never been
# measured for HEATING on this deck:
#   h1  filter_npass 8  -- 31% NOISE reduction was measured; heating effect unknown
#   h2  particle_shape 3 -- expected large, cost estimated at +40%, both unverified
#   h4  ppc 400          -- 1/N_ppc is textbook, unverified here
# h3 tests whether the two cheap levers compose. h0 is the control.
#
# METRIC: dT_0 in the far upstream over the 30 t_ab window. The baseline rises ~34 eV, so a
# 2x improvement reads as ~17 eV against a ~0.3% noise floor -- no model in the way, unlike
# the fitted asymptote which needs a long clean window.
#
# NOT a substitute for theta_implicit_em, which removes the heating by construction (exact
# energy conservation at theta = 0.5) for ~2.1x the cost. This matrix exists to find out
# whether the cheap explicit levers get close enough to make that unnecessary.

set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
GPUS="${GPUS:-0,1}"
VARIANTS=(h0_baseline h1_filter8 h2_shape3 h3_filter8_shape3 h4_ppc400)
[[ $# -gt 0 ]] && VARIANTS=("$@")

echo "=== grid-heating pilot matrix: ${#VARIANTS[@]} variants, serial, GPUs $GPUS ==="
df -h "$ROOT" | tail -1 | awk '{print "  disk free: "$4" ("$5" used)"}'
echo

for v in "${VARIANTS[@]}"; do
    d="runs/heat_phase/$v"
    [[ -d $d ]] || { echo "  SKIP $v (no such run dir)"; continue; }
    # launch.sh refuses to overwrite a populated diags/, which is what we want on a re-run;
    # say so rather than letting the loop look like it silently skipped.
    if compgen -G "$d/diags/*" >/dev/null; then
        echo "  SKIP $v -- diags/ already populated (move it aside to redo)"
        continue
    fi
    echo "--- $v  $(date +%H:%M:%S) ---"
    scripts/launch.sh -g "$GPUS" --every-pct 25 --poll 60 "$d" || echo "  $v exited nonzero"
    echo "    done $(date +%H:%M:%S)  $(grep -oE 'mean [0-9.]+ s/step' "$d/progress.log" 2>/dev/null | tail -1)"
    df -h "$ROOT" | tail -1 | awk '{print "    disk free now: "$4}'
done

echo
echo "=== MATRIX DONE -- measure with ==="
echo "  for v in runs/heat_phase/*/; do python scripts/plot_upstream_beta.py \$v; done"
echo "  # or, for the single-number comparison:"
echo "  python scripts/grid_heating.py runs/heat_phase/h*/ --upstream-frac 0.95"
