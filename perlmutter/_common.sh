#!/bin/bash -l
# Shared environment + the one function that actually launches WarpX on Perlmutter.
# Sourced by job.sbatch (inside the allocation) and by submit.sh (outside it).

set -euo pipefail
KINSHOCK_PM="${KINSHOCK_PM:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export KINSHOCK_PM

[[ -f "$KINSHOCK_PM/site.conf" ]] || {
    echo "perlmutter: create $KINSHOCK_PM/site.conf from site.conf.example" >&2; exit 1; }
# site.conf is user-authored config, not code. Source it with `set -u` OFF: it legitimately
# refers to $PSCRATCH, which does not exist off Perlmutter, and an unbound-variable abort
# there makes every script in this directory impossible to dry-run anywhere else.
set +u
# shellcheck disable=SC1090
source "$KINSHOCK_PM/site.conf"
set -u

: "${NERSC_ACCOUNT:?set NERSC_ACCOUNT in site.conf (must end in _g for GPU)}"
: "${KINSHOCK_ROOT:?set KINSHOCK_ROOT in site.conf}"
[[ "$NERSC_ACCOUNT" == *_g ]] || echo "perlmutter: WARNING NERSC_ACCOUNT '$NERSC_ACCOUNT' does not end in _g; GPU jobs will be rejected" >&2

pm_profile() {
    [[ -n "${MY_PROFILE:-}" ]] && return 0
    # shellcheck disable=SC1090
    source "$HOME/perlmutter_gpu_warpx.profile"
}

# Resolve the binary for a build tag (A or B).
pm_binary() {
    local tag="${1:-$SWEEP_BUILD}" d
    case "$tag" in
        A) d="$WARPX_BUILD_A" ;;
        B) d="$WARPX_BUILD_B" ;;
        *) echo "pm_binary: unknown build tag '$tag'" >&2; return 2 ;;
    esac
    local bin
    bin="$(ls "$d"/bin/warpx.1d 2>/dev/null || ls "$d"/bin/warpx.1d.* 2>/dev/null | head -1)"
    [[ -x "$bin" ]] || { echo "pm_binary: no binary in $d/bin -- run perlmutter/build_warpx.sh $tag" >&2; return 1; }
    echo "$bin"
}

# run_warpx <run_dir_relative_to_repo> [label]
#
# THE INVARIANT THIS FUNCTION EXISTS TO PRESERVE: the generated deck sets no
# diag*.file_prefix, so WarpX writes plotfiles to diags/ RELATIVE TO THE LAUNCH CWD.
# Two runs started from a shared directory clobber each other's output -- that cost a
# rerun of the R2/R3 controls (RESULTS 2026-07-26) and is why scripts/launch.sh exists on
# chablis. Same rule here: always cd first.
#
# With a label, the run happens in a scratch working directory instead of the repo run
# dir -- that is how the A/B wall test gets a second replicate of an identical config
# without creating a duplicate config.yaml that would misrepresent what is being varied.
run_warpx() {
    local run_rel="$1" label="${2:-}"
    local run_dir="$KINSHOCK_ROOT/$run_rel"
    local build="${BINARY:-$SWEEP_BUILD}"
    local bin; bin="$(pm_binary "$build")"
    local deck; deck="$(ls "$run_dir"/inputs_* | head -1)"
    [[ -f "$deck" ]] || { echo "run_warpx: no deck in $run_dir" >&2; return 1; }

    local work="$run_dir"
    if [[ -n "$label" ]]; then
        work="${WORKROOT:?WORKROOT must be set when using a label}/$label"
        mkdir -p "$work"
        cp "$deck" "$work/"
        cp "$run_dir/config.yaml" "$work/"     # so the analysis scripts can read it
    fi

    if compgen -G "$work/diags/*" > /dev/null; then
        echo "run_warpx: $work/diags already has output -- refusing to overwrite it in place." >&2
        echo "           Move it aside or delete it, then resubmit." >&2
        return 1
    fi

    pm_profile
    export MPICH_OFI_NIC_POLICY=GPU
    export OMP_NUM_THREADS=16          # 16 physical cores per GPU; avoids hyperthreading

    cd "$work"                          # THE POINT: diags/ lands here
    echo "run_warpx: $(basename "$work")  build=$build  bin=$bin"
    echo "run_warpx: cwd=$PWD"

    # Progress logger: every run must leave a progress.log (project convention). Started
    # before srun so it sees run.log from the first poll; it exits when WarpX does.
    local logger_pid=""
    if [[ -f "$KINSHOCK_ROOT/scripts/run_progress_logger.py" ]]; then
        python3 "$KINSHOCK_ROOT/scripts/run_progress_logger.py" "$work" \
            > "$work/logger.out" 2>&1 &
        logger_pid=$!
    fi

    # `set -e` off across the run: a WarpX failure must NOT abort this function, or the
    # array task dies before the logger is reaped and before --verify reports why. A
    # non-zero exit is information to be surfaced, not a reason to skip the diagnostics.
    local rc=0
    set +e
    srun --cpu-bind=cores "$bin" "$(basename "$deck")" > run.log 2>&1
    rc=$?
    set -e
    if [[ -n "$logger_pid" ]]; then
        sleep 35                       # let the logger write its DONE line before it goes
        kill "$logger_pid" 2>/dev/null || true
    fi

    echo "run_warpx: exit $rc"
    if [[ $rc -ne 0 ]]; then
        echo "run_warpx: FAILED — last 15 lines of run.log:" >&2
        tail -15 run.log >&2 || true
    fi
    # The post-run check that closes the "config = what was simulated" loop, INCLUDING the
    # unused-input scan that would have caught reflect_symmetry_axis being ignored.
    if [[ -z "$label" ]]; then
        python3 "$KINSHOCK_ROOT/scripts/make_inputs.py" "$run_dir" --verify || true
    fi
    return $rc
}
