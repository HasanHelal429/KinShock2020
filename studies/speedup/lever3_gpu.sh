#!/bin/bash
# Lever 3: GPU (build_cuda1d).
#
# This is the lever with the best cost/benefit and the worst logistics.
#
# WHY IT SHOULD WIN. 98% of this run is particle work (pilot TinyProfiler: GatherAndPush
# 49.9%, CurrentDeposition 20.9%, collisions 9.9%, Redistribute 8.5%) over 6.0e6
# macroparticles in 1D. That is close to ideal GPU shape, and the state is small:
# 6.0e6 particles is well inside one RTX 4070's 12 GiB.
#
# THE BINARY QUESTION. /home/hhelal/warpx-cda/build_cuda1d/bin/warpx.1d is dated
# Jul 28 18:05. It does contain the custom operators (17 `target_injector` strings; a
# test_1d_particle_heater target in bin/). But commit 9f981dea2 landed Jul 31 -- AFTER
# that build -- and its message reads:
#     "nvcc rejects an extended __device__ lambda inside a private or protected member
#      function ... which broke the CUDA build."
# So the binary predates a CUDA-specific fix to the very operator this run depends on.
# The change is behaviour-neutral by its own description ("The function needs no member
# state"), so the Jul 28 binary is probably fine -- but "probably" is not good enough for
# a physics result, which is why stage 3C below is an explicit agreement check against
# the CPU binary rather than a timing run.
#
# A REBUILD IS EXPENSIVE and should be a fallback, not the first move: ParticleHeater.H
# is included by Source/WarpX.H, which is included by nearly all 391 translation units,
# so the edit invalidates almost the whole CUDA build. ccache cannot save it (the
# preprocessed output of every TU changed). Estimate 30-50 min at -j24, up to ~90.
#
# GPU 1 is left alone throughout so the second card stays free for other users.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$HERE/../../runs/R1_paper_470eV"
GPUBIN="${KINSHOCK_WARPX_CUDA:-/home/hhelal/warpx-cda/build_cuda1d/bin/warpx.1d}"
CPUBIN="${KINSHOCK_WARPX:-/home/hhelal/warpx-cda/build/bin/warpx.1d}"
STEPS=1500
CAP=600

export CUDA_VISIBLE_DEVICES=0

echo "=== Lever 3: GPU ==="
nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv
echo "gpu binary $GPUBIN ($(date -r "$GPUBIN" -Is 2>/dev/null))"
echo "cpu binary $CPUBIN ($(date -r "$CPUBIN" -Is 2>/dev/null))"
echo "baseline: 0.1113 s/step on 8 CPU threads"
echo

# --- 3A: as-configured. 235 grids of 112-128 cells is a CPU-tuned decomposition and a
# bad one for a GPU -- each kernel launch gets ~25k particles. Measure it anyway, since
# it is the zero-change baseline.
echo "--- 3A: GPU, deck's default decomposition (235 grids) ---"
"$HERE/bench.sh" l3A_gpu_default "$RUN" 1 "$GPUBIN" "$STEPS" "$CAP"

# --- 3B: one big box. On a GPU you want few large boxes so kernels have enough work.
echo "--- 3B: GPU, single box (max_grid_size = whole domain) ---"
"$HERE/bench.sh" l3B_gpu_1box "$RUN" 1 "$GPUBIN" "$STEPS" "$CAP" amr.max_grid_size=30000
echo "--- 3B2: GPU, 8 boxes ---"
"$HERE/bench.sh" l3B2_gpu_8box "$RUN" 1 "$GPUBIN" "$STEPS" "$CAP" amr.max_grid_size=4096

# --- 3C: agreement check. This is the stage that decides whether the stale binary is
# usable. Bit-identity is NOT expected -- reduction order differs and the collision
# operator draws from a different RNG stream -- so the test is that the reduced
# diagnostics track to the level the collision noise allows, not that they match exactly.
echo
echo "--- 3C: CPU reference at the same step count, for the agreement check ---"
"$HERE/bench.sh" l3C_cpu_ref "$RUN" 8 "$CPUBIN" "$STEPS" "$CAP"

echo
echo "=== Lever 3 summary ==="
cat "$HERE"/out/l3*/result.txt 2>/dev/null
echo
echo "--- GPU vs CPU agreement (EP = particle energy, PN = particle number) ---"
python3 "$HERE/compare_diags.py" \
    "$HERE/out/l3B_gpu_1box" "$HERE/out/l3C_cpu_ref" 2>&1 || true
