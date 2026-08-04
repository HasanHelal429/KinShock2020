#!/bin/bash
# Lever 2: algo.evolve_scheme = theta_implicit_em.
#
# WHY THIS IS THE ONE THAT MATTERS. The docs (parameters.rst:246-249) claim exactly the
# two properties that are blocking R1_paper_470eV:
#   "Robust to finite-grid instability (does not require cells that resolve the plasma
#    Debye length)"    <- the pilot failed on dz/lambda_D,amb = 6.07
#   "Numerically stable for large dt (does not require resolving the plasma period or
#    satisfying the CFL condition for light waves)"
# and theta = 0.5 is *exact* energy conservation, i.e. no grid heating by construction.
#
# COMPATIBILITY, checked in source rather than assumed:
#   - Custom operators: TargetInjector and ParticleHeater are applied in the OUTER Evolve
#     loop (WarpXEvolve.cpp:286,291), before the `if (m_implicit_solver)` branch at :299.
#     They are scheme-agnostic. (Contrast TargetInjector.cpp:225, which hard-aborts under
#     mesh refinement -- there is no such guard for implicit.)
#   - Collisions: mypc->doCollisions() is called inside the implicit branch
#     (WarpXEvolve.cpp:301), and the scheme's first reference is Angus et al.,
#     "...implicit particle-in-cell method coupled with a binary Monte-Carlo algorithm
#     for Coulomb collisions" -- this pairing is the published use case.
#   - Deposition: the deck uses Esirkepov (run.log:39, the FDTD default). Esirkepov IS
#     allowed for theta_implicit_em, but NOT with use_mass_matrices_jacobian=true, which
#     is the option that "can enable large speed ups for simulations with many particles"
#     -- and we have 6.0e6. Stage C therefore switches to villasenor, which is both
#     charge-conserving and mass-matrix compatible.
#   - 1D: Examples/Tests/implicit/inputs_test_1d_theta_implicit_picard exists and uses
#     particle_shape = 2, same as this deck.
#
# STAGED WITH GATES, because each stage's cost depends on the previous stage's answer.
# Stage A is the go/no-go and costs ~2 min. Pass a stage name to run just one.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$HERE/../../runs/R1_paper_470eV"
BIN="${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
THR="${THREADS:-8}"          # hold threads fixed so this measures the SCHEME, not lever 1
STAGE="${1:-all}"

IMPL=(algo.evolve_scheme=theta_implicit_em implicit_evolve.theta=0.5)

echo "=== Lever 2: theta_implicit_em (threads=$THR, held fixed) ==="
echo "baseline: 0.1113 s/step explicit. Expect implicit to be SLOWER PER STEP"
echo "(each Picard iteration redoes the push+deposit that is 71% of runtime) and to"
echo "win only if it buys a larger dt and/or removes the grid heating."
echo

# --- Stage A: does it run at all with our operators, collisions and BCs? -------------
# picard.verbose=1 so the log carries the per-step iteration count -- that number is
# what sets the per-step cost multiplier, so it is the real output of this stage.
# require_convergence=0: a non-convergent step should be VISIBLE, not fatal.
if [[ $STAGE == all || $STAGE == A ]]; then
    echo "--- Stage A: 100-step smoke, Picard, same dt (cap 300 s) ---"
    "$HERE/bench.sh" l2A_smoke_picard "$RUN" "$THR" "$BIN" 100 300 \
        "${IMPL[@]}" \
        implicit_evolve.nonlinear_solver=picard \
        picard.verbose=1 picard.require_convergence=0 \
        picard.max_iterations=100 picard.relative_tolerance=1.0e-6
    echo
    # PicardSolver.H:222 prints "Picard: exiting at iter = N" once per step -- that N is
    # the iteration count, and the per-step cost multiplier is roughly 1 + 0.71*(N-1)
    # since each iteration redoes the push+deposit (71% of explicit runtime).
    echo "Picard iterations per step (count, iters):"
    grep -oE "exiting at iter = *[0-9]+" "$HERE/out/l2A_smoke_picard/run.log" 2>/dev/null |
        grep -oE "[0-9]+$" | sort -n | uniq -c | tail -20 || echo "  (no Picard lines -- check run.log)"
    echo
    echo "GATE: if Stage A aborted, stop and read out/l2A_smoke_picard/run.log."
fi

# --- Stage B: per-step cost multiplier at the SAME dt --------------------------------
if [[ $STAGE == all || $STAGE == B ]]; then
    echo "--- Stage B: 1500-step timing, Picard, same dt (cap 1800 s) ---"
    "$HERE/bench.sh" l2B_time_picard "$RUN" "$THR" "$BIN" 1500 1800 \
        "${IMPL[@]}" \
        implicit_evolve.nonlinear_solver=picard \
        picard.verbose=0 picard.require_convergence=0 \
        picard.max_iterations=100 picard.relative_tolerance=1.0e-6
fi

# --- Stage C: the actual prize -- can dt grow? ---------------------------------------
# Picard is documented as "requires small time steps; often non-convergent for large
# time steps", so large dt needs Newton (PS-JFNK) with a preconditioner and the mass
# matrices. cfl 0.75 -> 7.5 is 10x dt, which would take the 3.22e6-step production run
# to 3.2e5 steps. Two cfl points, because 10x may simply not converge.
if [[ $STAGE == all || $STAGE == C ]]; then
    NEWTON=(implicit_evolve.nonlinear_solver=newton
            algo.current_deposition=villasenor
            implicit_evolve.use_mass_matrices_jacobian=1
            implicit_evolve.use_mass_matrices_pc=1
            jacobian.pc_type=pc_curl_curl_mlmg
            newton.verbose=0 newton.require_convergence=0
            pc_curl_curl_mlmg.verbose=0 gmres.verbose_int=0)
    for CFL in 3.0 7.5; do
        echo "--- Stage C: 1500 steps, Newton + mass matrices, cfl=$CFL (cap 1800 s) ---"
        "$HERE/bench.sh" "l2C_newton_cfl${CFL}" "$RUN" "$THR" "$BIN" 1500 1800 \
            "${IMPL[@]}" "${NEWTON[@]}" warpx.cfl="$CFL"
    done
fi

echo
echo "=== Lever 2 summary ==="
cat "$HERE"/out/l2*/result.txt 2>/dev/null
echo
echo "NOTE: none of these stages measures GRID HEATING -- they measure cost and"
echo "convergence only. Confirming the scheme actually kills the +3.9%/3.41 t_ab needs a"
echo "run matched to the pilot's physical time, and needs numerics.evolve_scheme plumbed"
echo "through scripts/make_inputs.py first (no passthrough block exists today, and the"
echo "repo rule is never to hand-edit a deck). See README.md."
