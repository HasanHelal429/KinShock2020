#!/bin/bash -l
# Build the warpx-cda 1D CUDA binary on Perlmutter, at a named commit.
#
#   perlmutter/build_warpx.sh A     # the specular-wall binary (matches existing results)
#   perlmutter/build_warpx.sh B     # + reflect_symmetry_axis (the pi-rotation wall)
#   perlmutter/build_warpx.sh both
#
# WHY TWO BINARIES. `boundary.reflect_symmetry_axis` is fork-only and lived on an
# unmerged branch, so every result this project has produced silently used plain specular
# reflection instead of the pi-rotation symmetry wall its config asked for (RESULTS
# 2026-08-11). A is what produced those results; B is A plus the two cherry-picked commits
# that implement it. The A/B pair is the first measurement of the difference -- see
# scripts/ab_wall_test.py. If you only want the sweep, build whichever SWEEP_BUILD names.
#
# The heavy lifting is upstream's: Tools/machines/perlmutter-nersc/ carries the module
# set, AMREX_CUDA_ARCH=8.0 (A100) and an installer for boost/adios2/blaspp/lapackpp.
# Run install_gpu_dependencies.sh ONCE before the first build; it is idempotent but slow.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$HERE/site.conf" ]] || { echo "create $HERE/site.conf from site.conf.example" >&2; exit 1; }
set +u   # see _common.sh: site.conf may reference $PSCRATCH
# shellcheck disable=SC1090
source "$HERE/site.conf"
set -u

WHICH="${1:-both}"

[[ -n "${MY_PROFILE:-}" ]] || {
    echo "build: sourcing the WarpX Perlmutter profile"
    # shellcheck disable=SC1090
    source "$HOME/perlmutter_gpu_warpx.profile"
}
: "${AMREX_CUDA_ARCH:?the WarpX profile did not load - check \$HOME/perlmutter_gpu_warpx.profile}"

build_one() {
    local tag="$1" commit="$2" bdir="$3"
    echo "=== building $tag at $commit -> $bdir"
    cd "$WARPX_SRC"
    # Detached HEAD on purpose: these are build trees, not working branches, and it makes
    # the binary's provenance a commit rather than "whatever the branch was that day" --
    # which is exactly how the reflect_symmetry_axis gap went unnoticed for three weeks.
    git checkout --detach "$commit"
    cmake -S . -B "$bdir" \
        -DWarpX_DIMS=1 -DWarpX_COMPUTE=CUDA \
        -DAMReX_CUDA_ARCH="${AMREX_CUDA_ARCH}" \
        -DWarpX_MPI=ON -DWarpX_MPI_THREAD_MULTIPLE=ON \
        -DWarpX_PRECISION=DOUBLE -DWarpX_OPENPMD=ON \
        -DWarpX_EB=ON -DWarpX_QED=ON -DWarpX_FFT=OFF \
        -DWarpX_APP=ON -DWarpX_LIB=ON -DWarpX_PYTHON=OFF \
        -DCMAKE_BUILD_TYPE=Release
    cmake --build "$bdir" -j 16
    echo "--- $tag built: $(ls "$bdir"/bin/warpx.1d* 2>/dev/null | head -1)"
    # The check that would have caught the original bug: does the binary actually contain
    # the feature? A silently-missing fork-only input is invisible until a run ignores it.
    if strings "$bdir"/bin/warpx.1d* 2>/dev/null | grep -q reflect_symmetry_axis; then
        echo "--- $tag DOES implement reflect_symmetry_axis (pi-rotation wall)"
    else
        echo "--- $tag does NOT implement reflect_symmetry_axis (specular wall)"
    fi
}

case "$WHICH" in
    A)    build_one A "$WARPX_COMMIT_A" "$WARPX_BUILD_A" ;;
    B)    build_one B "$WARPX_COMMIT_B" "$WARPX_BUILD_B" ;;
    both) build_one A "$WARPX_COMMIT_A" "$WARPX_BUILD_A"
          build_one B "$WARPX_COMMIT_B" "$WARPX_BUILD_B" ;;
    *)    echo "usage: $0 [A|B|both]" >&2; exit 2 ;;
esac
echo
echo "Reminder: commit $WARPX_COMMIT_B must be REACHABLE FROM ORIGIN for this to work"
echo "on Perlmutter. As of 2026-08-11 the cherry-picks exist only in the chablis clone."
