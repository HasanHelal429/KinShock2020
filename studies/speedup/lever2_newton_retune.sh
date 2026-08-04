#!/bin/bash
# Lever 2, retune: why was Newton (PS-JFNK) 568x explicit, and is it salvageable?
#
# THE BUG IN THE FIRST ATTEMPT. ImplicitSolver.cpp:811-817 gates the cheap particle path
# on BOTH flags:
#
#     if (m_use_mass_matrices_jacobian && m_skip_particle_picard_init) {
#         options.max_particle_iterations = 1;
#         options.particle_tolerance     = 0.0;
#     } else {
#         options.max_particle_iterations = m_max_particle_iterations;   // default 21
#     }
#
# Stage C set use_mass_matrices_jacobian=1 but left skip_particle_picard_init at its
# default false, so EVERY Newton iteration ran a full 21-iteration particle Picard update.
# At ~5 Newton iterations that is ~100 particle pushes per step, which is the right order
# for the observed 568x. The docs say the same thing in prose: skip_particle_picard_init
# "skips the full Picard update of the particles on the initial Newton step ... can enhance
# the overall efficiency of the Newton solver."
#
# The other Stage C mistake was diagnostic, not physical: newton.verbose=0 and
# gmres.verbose_int=0 meant the log carried no iteration counts, so the cost could not be
# attributed. Verbosity is ON here -- that is the entire point of the exercise.
#
# Runs SHORT: at even 5 s/step a 20-step point is 100 s, and the iteration counts (not the
# step count) are the measurement. N1 is deliberately the BROKEN config with verbosity on,
# so the ladder shows the mechanism rather than just asserting it.
#
# Reference numbers at cfl 0.75, 8 threads, this deck:
#   explicit                 0.09945 s/step (median)
#   theta-implicit, Picard   1.33961 s/step = 13.5x, 21-23 iterations, all converged
# Break-even for implicit against explicit needs dt to grow past ~13.5x, i.e. cfl > 10.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Resolved ONCE here, which is fine now that the run is settled at its phased path -- but
# it was not fine mid-reorganisation. On 2026-08-04 this script started seconds before
# runs/R1_paper_470eV was moved under R1_phase/; N1 captured the old path, the move landed,
# and N2-N4 all died on `cd: No such file or directory`. A start-time fallback protects
# against a move that already happened, not one that happens while you run. If the tree is
# in flux, pass the run dir explicitly.
RUN="${1:-$HERE/../../runs/R1_phase/R1_paper_470eV}"
[[ -d "$RUN" ]] || RUN="$HERE/../../runs/R1_paper_470eV"
[[ -d "$RUN" ]] || { echo "no run dir found (tried R1_phase/ and bare)" >&2; exit 2; }
BIN="${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
THR="${THREADS:-8}"

BASE=(algo.evolve_scheme=theta_implicit_em implicit_evolve.theta=0.5
      implicit_evolve.nonlinear_solver=newton
      algo.current_deposition=villasenor
      newton.require_convergence=0)
# Newton iteration counts are a handful per step, so verbose=1 is nearly free and is what
# lets each point be attributed. GMRES can be hundreds per Newton step -- verbose only in
# the diagnostic point, or the log write itself distorts the timing.
VERB=(newton.verbose=1 gmres.verbose_int=0 pc_curl_curl_mlmg.verbose=0)

echo "=== Lever 2 retune: Newton, verbosity on (threads=$THR) ==="
echo "explicit 0.09945 | Picard 1.33961 (13.5x) | broken Newton 56.49 (568x)"
echo

# --- N1: the BROKEN config, verbosity on. Shows the 21 particle iterations directly. ----
echo "--- N1: reproduce Stage C with verbosity ON (5 steps, cap 400 s) ---"
"$HERE/bench.sh" l2N1_broken_verbose "$RUN" "$THR" "$BIN" 5 400 \
    "${BASE[@]}" newton.verbose=1 gmres.verbose_int=2 \
    implicit_evolve.use_mass_matrices_jacobian=1 \
    implicit_evolve.use_mass_matrices_pc=1 \
    jacobian.pc_type=pc_curl_curl_mlmg

# --- N2: the hypothesised fix -- both flags, so max_particle_iterations collapses to 1 ---
echo "--- N2: + skip_particle_picard_init=1  (20 steps, cap 600 s) ---"
"$HERE/bench.sh" l2N2_skipinit "$RUN" "$THR" "$BIN" 20 600 \
    "${BASE[@]}" "${VERB[@]}" \
    implicit_evolve.use_mass_matrices_jacobian=1 \
    implicit_evolve.skip_particle_picard_init=1 \
    implicit_evolve.use_mass_matrices_pc=1 \
    jacobian.pc_type=pc_curl_curl_mlmg

# --- N3: cheaper preconditioner. pc_jacobi captures plasma response via the diagonal
# mass matrices only -- no multigrid solve per GMRES iteration.
echo "--- N3: N2 but pc_jacobi  (20 steps, cap 600 s) ---"
"$HERE/bench.sh" l2N3_pc_jacobi "$RUN" "$THR" "$BIN" 20 600 \
    "${BASE[@]}" "${VERB[@]}" \
    implicit_evolve.use_mass_matrices_jacobian=1 \
    implicit_evolve.skip_particle_picard_init=1 \
    implicit_evolve.use_mass_matrices_pc=1 \
    jacobian.pc_type=pc_jacobi pc_jacobi.verbose=0

# --- N4: no preconditioner at all. Isolates how much the PC is worth: without one, GMRES
# iteration count should climb sharply, so if N4 beats N2 the MLMG PC is the problem.
echo "--- N4: N2 with NO preconditioner  (20 steps, cap 600 s) ---"
"$HERE/bench.sh" l2N4_no_pc "$RUN" "$THR" "$BIN" 20 600 \
    "${BASE[@]}" "${VERB[@]}" \
    implicit_evolve.use_mass_matrices_jacobian=1 \
    implicit_evolve.skip_particle_picard_init=1

echo
echo "=== retune summary ==="
cat "$HERE"/out/l2N*/result.txt 2>/dev/null
echo
echo "--- Newton iterations per step, and GMRES per Newton where logged ---"
for d in "$HERE"/out/l2N*/; do
    n=$(basename "$d")
    it=$(grep -oE "Newton: iter = *[0-9]+" "$d/run.log" 2>/dev/null | grep -oE "[0-9]+$" | sort -n | tail -1)
    gm=$(grep -cE "GMRES" "$d/run.log" 2>/dev/null)
    pi=$(grep -oE "max particle iterations: *[0-9]+" "$d/run.log" 2>/dev/null | head -1 | grep -oE "[0-9]+$")
    printf "  %-24s max_newton_iter=%-5s gmres_lines=%-6s max_particle_iterations=%s\n" \
        "$n" "${it:-?}" "${gm:-0}" "${pi:-?}"
done
echo
echo "READ: if N2 collapses the cost, the 568x was the 21-iteration particle Picard update"
echo "and NOT the method. Net viability still needs dt: at cfl 0.75 anything above 13.5x"
echo "per step loses to explicit, and above ~1x loses to running explicit at 16 threads."
