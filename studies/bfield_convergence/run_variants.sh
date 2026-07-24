#!/usr/bin/env bash
# Numerical-convergence test for B-field fluctuations.
#
# Runs the SAME physics with ONE numerical knob changed per variant, long enough
# (default t*wci~0.56) for the far-upstream fluctuation to develop, so
# analyze.py / bfield_diagnostic.py can check whether it COLLAPSES (=> numerical,
# a grid instability from under-resolved lambda_D) or PERSISTS (=> physical).
#
# Usage:
#   studies/bfield_convergence/run_variants.sh <base_deck> [out_dir] [max_step] [threads]
#   WARPX=/path/to/warpx.1d studies/bfield_convergence/run_variants.sh runs/R1_core_half/inputs_kinshock_R1_core_half
#
# Scratch WarpX output goes to out_dir (default ./scratch, gitignored).
set -euo pipefail
WARPX="${WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
BASE_DECK="${1:?usage: run_variants.sh <base_inputs_deck> [out_dir] [max_step] [threads]}"
OUT="${2:-$(dirname "$0")/scratch}"
STEPS="${3:-25000}"
THREADS="${4:-8}"
DIAG=$((STEPS/5))
mkdir -p "$OUT"
run(){ local name=$1 steps=$2 diag=$3; shift 3
  local wd="$OUT/$name"; rm -rf "$wd"; mkdir -p "$wd"; cd "$wd"
  cp "$BASE_DECK" inputs
  { echo; echo "# convergence-variant override (ParmParse last-wins)";
    echo "max_step = $steps"; echo "diag1.intervals = $diag";
    echo "EP.intervals = 100000000"; echo "PN.intervals = 100000000";
    echo "warpx.use_filter = 1"; printf '%s\n' "$@"; } >> inputs
  echo "[$(date +%H:%M:%S)] $name (max_step=$steps) ..."
  OMP_NUM_THREADS=$THREADS OMP_PROC_BIND=spread OMP_PLACES=cores \
    taskset -c 0-$((THREADS-1)) "$WARPX" inputs > run.log 2>&1
  echo "[$(date +%H:%M:%S)] $name done: $(grep -oE 'STEP [0-9]+ ends' run.log | tail -1)"
}
NC=$(grep -oE 'amr\.n_cell = [0-9]+' "$BASE_DECK" | grep -oE '[0-9]+' | tail -1)
run baseline "$STEPS"      "$DIAG"                              # dz as-is, filter=1, shape=2
run filt8    "$STEPS"      "$DIAG" "warpx.filter_npass_each_dir = 8"   # heavy current filter
run shape3   "$STEPS"      "$DIAG" "algo.particle_shape = 3"           # cubic shape
run finer_dz $((STEPS*2))  $((DIAG*2)) "amr.n_cell = $((NC*2))"        # resolve Debye (dz/2, dt/2)
echo "ALL VARIANTS DONE -> $OUT  (analyze: studies/bfield_convergence/analyze.py $OUT)"
