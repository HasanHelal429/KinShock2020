#!/bin/bash
# Lever 2 follow-up: does theta-implicit stay affordable as dt grows?
#
# WHY NOW. The retune found that the cost was the PRECONDITIONER, not the particle Picard
# loop I first blamed: pc_curl_curl_mlmg 23.74 s/step vs pc_jacobi 0.43928 vs NO
# preconditioner 0.43389. pc_jacobi ~= no-PC, so MLMG is uniquely pathological here -- 54x
# against using nothing. With that removed, Newton is 4.36x explicit at cfl 0.75 and 3.1x
# FASTER than Picard, which moves break-even from cfl > 13.5 down to cfl > 3.3 and makes
# large dt worth measuring for the first time.
#
# WHAT DECIDES IT. At cfl X the step count is 3,224,046 * 0.75/X, so the projected full run
# is (that) * S * 1.279 / 86400 days for a measured per-step cost S:
#
#   cfl 3.0 ->   806,012 steps    S=0.434 -> 5.18 d    S=0.87 -> 10.35 d
#   cfl 7.5 ->   322,405 steps    S=0.434 -> 2.07 d    S=0.87 ->  4.14 d
#
# Measured references: explicit 8 thr 5.33 d, explicit 20 thr 2.93 d, GPU 1 box 0.68 d.
# To match the GPU at cfl 7.5 needs S < 0.142, i.e. S would have to FALL 3x while dt grows
# 10x. That is not going to happen, so this is not a search for the fastest configuration --
# it is establishing whether the energy-conserving scheme is affordable at all, because it
# is the only lever that removes the grid heating rather than out-running it.
#
# CONVERGENCE IS THE REAL OUTPUT. newton.require_convergence=0 keeps a bad step from being
# fatal, which means a non-converged solve proceeds silently and the numbers look fine. The
# script greps the log for it explicitly rather than trusting the s/step.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="${1:-$HERE/../../runs/R1_phase/R1_paper_470eV}"
[[ -d "$RUN" ]] || { echo "no run dir: $RUN" >&2; exit 2; }
BIN="${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
THR="${THREADS:-8}"

# N4's configuration: mass-matrix Jacobian, cheap particle path, and NO preconditioner.
BASE=(algo.evolve_scheme=theta_implicit_em implicit_evolve.theta=0.5
      implicit_evolve.nonlinear_solver=newton
      algo.current_deposition=villasenor
      implicit_evolve.use_mass_matrices_jacobian=1
      implicit_evolve.skip_particle_picard_init=1
      newton.require_convergence=0 newton.verbose=1)

echo "=== Lever 2 follow-up: large dt, no preconditioner (threads=$THR) ==="
echo "explicit 0.09945 | Picard 1.33961 | Newton no-PC @cfl 0.75 0.43389"
echo

for CFL in 3.0 7.5; do
    echo "--- cfl=$CFL, 200 steps (cap 900 s) ---"
    "$HERE/bench.sh" "l2D_cfl${CFL}" "$RUN" "$THR" "$BIN" 200 900 \
        "${BASE[@]}" gmres.verbose_int=0 warpx.cfl="$CFL"
done

# Full GMRES verbosity, 5 steps only. Separate from the timing points on purpose: the
# per-iteration log writes are heavy enough to distort s/step, which is exactly the mistake
# that left the first Newton attempt undiagnosable.
echo "--- cfl=7.5, 5 steps, FULL verbosity (iteration counts, cap 600 s) ---"
"$HERE/bench.sh" l2D_cfl7.5_verbose "$RUN" "$THR" "$BIN" 5 600 \
    "${BASE[@]}" gmres.verbose_int=2 warpx.cfl=7.5

echo
echo "=== results ==="
cat "$HERE"/out/l2D_*/result.txt 2>/dev/null
echo
echo "=== convergence (require_convergence=0, so a failed solve is SILENT in s/step) ==="
for d in "$HERE"/out/l2D_*/; do
    n=$(basename "$d")
    steps=$(grep -c "This step" "$d/run.log" 2>/dev/null)
    maxit=$(grep -ciE "maximum iteration reached|did not converge|not converged" "$d/run.log" 2>/dev/null)
    conv=$(grep -ciE "satisfied (relative|absolute) tolerance" "$d/run.log" 2>/dev/null)
    gm=$(grep -cE "GMRES" "$d/run.log" 2>/dev/null)
    printf "  %-22s steps=%-5s converged=%-6s hit_max_iter=%-5s gmres_lines=%s\n" \
        "$n" "$steps" "$conv" "$maxit" "$gm"
done
echo
echo "READ: hit_max_iter > 0 means the solve is NOT converging at that dt and the s/step"
echo "is meaningless -- large dt bought nothing. Projected full run at cfl 7.5 is"
echo "322,405 steps x S x 1.279; compare against GPU 0.68 d and explicit-20thr 2.93 d."
