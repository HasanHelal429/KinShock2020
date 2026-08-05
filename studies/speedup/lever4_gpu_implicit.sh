#!/bin/bash
# Levers 2 x 3: theta-implicit ON the GPU. The only combination that could give both
# 7.89x throughput AND immunity to the grid heating that is actually blocking
# R1_paper_470eV -- everything else so far trades one for the other.
#
# FEASIBILITY, checked before spending anything:
#   * No GPU guards or CPU-only aborts anywhere in Source/FieldSolver/ImplicitSolvers/
#     or Source/NonlinearSolvers/.
#   * Those two directories are UNCHANGED since the Jul 28 CUDA build (`git log --since`),
#     so unlike ParticleHeater the binary's implicit code is current. The Jul 28 build has
#     separately been validated against the CPU build on this deck: 30 diagnostic rows over
#     1479 steps, no diverging column, ambient electron count bit-identical.
#   * ImplicitSolver.cpp:679 -- mass matrices require Direct or Villasenor deposition. We
#     use villasenor.
#   * amr.max_grid_size=30000 is mandatory, not optional: the deck's CPU-tuned 235 boxes
#     cost 6.5x on the GPU (1.21x vs 7.89x for explicit).
#
# REFERENCES to beat, all measured on this deck:
#   explicit GPU 1 box   0.01415 s/step   -> 0.68 d full run
#   explicit CPU 20 thr  0.06143          -> 2.93 d
#   implicit CPU, no PC, cfl 0.75  0.43389 -> 20.71 d
#   implicit CPU, no PC, cfl 3.0   0.49834 ->  5.95 d
#
# THE PRIZE. At cfl 7.5 the run is 322,405 steps, so a per-step cost S gives
# S * 322405 * 1.279 / 86400 days. S = 0.06 would be 0.29 d -- faster than the explicit GPU
# run AND energy-conserving. S = 0.14 matches the explicit GPU run exactly.
#
# CONVERGENCE OVER SPEED. newton.require_convergence=0 means a non-converging solve at
# large dt still reports a plausible s/step, so the last point runs with full GMRES
# verbosity and the summary greps for max-iteration hits rather than trusting the timing.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="${1:-$HERE/../../runs/R1_phase/R1_paper_470eV}"
[[ -d "$RUN" ]] || { echo "no run dir: $RUN" >&2; exit 2; }
GPU="${KINSHOCK_WARPX_CUDA:-/home/hhelal/warpx-cda/build_cuda1d/bin/warpx.1d}"
[[ -x "$GPU" ]] || { echo "no CUDA binary: $GPU" >&2; exit 2; }

export CUDA_VISIBLE_DEVICES=0        # leave card 1 for other users

BASE=(algo.evolve_scheme=theta_implicit_em implicit_evolve.theta=0.5
      implicit_evolve.nonlinear_solver=newton
      algo.current_deposition=villasenor
      implicit_evolve.use_mass_matrices_jacobian=1
      implicit_evolve.skip_particle_picard_init=1
      newton.require_convergence=0 newton.verbose=1
      amr.max_grid_size=30000)

echo "=== theta-implicit on GPU (1 box) ==="
nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader
echo "beat: explicit GPU 0.01415 (0.68 d) | implicit CPU cfl3.0 0.49834 (5.95 d)"
echo

for CFL in 0.75 3.0 7.5; do
    echo "--- GPU implicit, cfl=$CFL, 200 steps (cap 900 s) ---"
    "$HERE/bench.sh" "l4_gpu_impl_cfl${CFL}" "$RUN" 1 "$GPU" 200 900 \
        "${BASE[@]}" gmres.verbose_int=0 warpx.cfl="$CFL"
done

echo "--- GPU implicit, cfl=7.5, 5 steps, FULL verbosity (cap 600 s) ---"
"$HERE/bench.sh" l4_gpu_impl_cfl7.5_verbose "$RUN" 1 "$GPU" 5 600 \
    "${BASE[@]}" gmres.verbose_int=2 warpx.cfl=7.5

echo
echo "=== results ==="
cat "$HERE"/out/l4_*/result.txt 2>/dev/null
echo
echo "=== convergence (a failed solve is SILENT in s/step) ==="
for d in "$HERE"/out/l4_*/; do
    n=$(basename "$d")
    steps=$(grep -c "This step" "$d/run.log" 2>/dev/null)
    maxit=$(grep -ciE "maximum iteration reached|did not converge|not converged" "$d/run.log" 2>/dev/null)
    conv=$(grep -ciE "satisfied (relative|absolute) tolerance" "$d/run.log" 2>/dev/null)
    printf "  %-34s steps=%-5s converged=%-6s hit_max_iter=%s\n" \
        "$n" "$steps" "$conv" "$maxit"
done
echo
echo "READ: full run at cfl 7.5 = 322,405 steps x S x 1.279. Under 0.142 beats the"
echo "explicit GPU run outright; anything with hit_max_iter > 0 is not a result."
