#!/bin/bash -l
# Submit a KinShock2020 job array on Perlmutter.
#
#   perlmutter/submit.sh sweep          # the six S_phase points, one GPU each
#   perlmutter/submit.sh eps            # the four E_phase eps-ladder rungs, one GPU each
#   perlmutter/submit.sh mag            # the four M_phase B0-scan rungs, one GPU each
#   perlmutter/submit.sh mu             # the two U_phase physical-mass-ratio twins
#   perlmutter/submit.sh wall           # the W_phase reflecting-wall smoke test
#   perlmutter/submit.sh ab             # the wall A/B: 2 runs per binary
#   perlmutter/submit.sh sweep --dry    # print the sbatch command, submit nothing
#
# Everything site-specific comes from perlmutter/site.conf. sbatch #SBATCH directives
# cannot read shell variables, which is why -A/-q/-t/--array are passed here on the
# command line rather than baked into job.sbatch.

set -euo pipefail
KINSHOCK_PM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export KINSHOCK_PM
# shellcheck disable=SC1090
source "$KINSHOCK_PM/_common.sh"

WHAT="${1:-}"; shift || true
DRY=0; QOS_OVERRIDE=""; TIME_OVERRIDE=""
# Run dirs given positionally REPLACE the target's default list, keeping its binary, QOS,
# WORKROOT and job name. This exists so a partial re-run is a subset of the same target
# rather than a hand-rolled sbatch line -- the three E_phase rungs that died on
# blocking_factor (2026-08-17) had to go back WITHOUT the one still running, and
# duplicating the launch invariants to do that is exactly how they drift.
SUBSET=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry)  DRY=1; shift ;;
        --qos)  QOS_OVERRIDE="$2"; shift 2 ;;
        --time) TIME_OVERRIDE="$2"; shift 2 ;;
        -*)     echo "unknown option '$1'" >&2; exit 2 ;;
        *)      SUBSET+=("$1"); shift ;;
    esac
done

case "$WHAT" in
  sweep)
    # Ordered cheapest-first so a mistake surfaces on the 8-minute run, not the 100-minute
    # one. ss_dz1_ppc100 is the domain control and should be read before the rest.
    RUNS=(runs/S_phase/ss_dz1_ppc100
          runs/S_phase/ss_dz2_ppc50
          runs/S_phase/ss_dz1_ppc400
          runs/S_phase/ss_dz2_ppc100
          runs/S_phase/ss_dz4_ppc25
          runs/S_phase/ss_dz4_ppc100)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-5"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="ksweep"
    ;;
  eps)
    # The eps = v_te,ab/c ladder (runs/E_phase). Four rungs between the 470 eV and 47 keV
    # endpoints, every Table I group AND every numerics group held fixed -- see
    # runs/E_phase/README.md. Ordered cheapest-first, as the sweep is: es_47keV is ~1 min
    # and es_1p5keV ~1.3 h, so a config mistake surfaces immediately.
    RUNS=(runs/E_phase/es_47keV
          runs/E_phase/es_15keV
          runs/E_phase/es_4p7keV
          runs/E_phase/es_1p5keV)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-3"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="keps"
    ;;
  mag)
    # The B0 scan at fixed eps (runs/M_phase). Four rungs above the anchor, which IS
    # runs/E_phase/es_47keV and is already run -- it is deliberately NOT in this list and
    # NOT duplicated as a config, per CLAUDE.md's A/B note. Every Table I group that E_phase
    # holds fixed is still fixed here EXCEPT the magnetization family; in particular
    # lambda_D,0/d_i0 and gamma are frozen while w_pe0/w_ce moves a decade, which is the
    # separation E_phase structurally cannot make. See runs/M_phase/README.md.
    #
    # Cheapest-first, as the others are: mb_r20 is ~2 min and mb_r100 ~10 min. Rungs 2-3
    # (mb_r20, mb_r32) are the ones to READ FIRST -- they reach w_pe0/w_ce = 31.6 while
    # beta_0 is still <= 2, where mb_r56/mb_r100 have a thermal-pressure-dominated upstream
    # that may not support a perpendicular shock at all.
    RUNS=(runs/M_phase/mb_r20
          runs/M_phase/mb_r32
          runs/M_phase/mb_r56
          runs/M_phase/mb_r100)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-3"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="kmag"
    ;;
  mu)
    # The physical-mass-ratio twins (runs/U_phase). Each is an existing run with
    # reference.mass_ratio as the ONLY choice changed, so grid heating, particle noise and
    # box-filling are common mode with the partner. They cross the foot-instability
    # threshold v_u/v_te = M_A/sqrt(beta_e mu) from 3.15 (Buneman) to 0.74 (MTSI) --
    # Umeda+2014. See runs/U_phase/README.md.
    #
    # 470 eV first: it is the one that carries the result. The 47 keV twin is the control
    # whose partner has no notch, so it is only interesting if it grows one.
    #
    # BOTH are ~4.5 GPU-h against SWEEP_TIME's 4 h default -- pass --time 10:00:00 or they
    # die at the wall. `shared` bills elapsed, not requested, so the headroom is free.
    RUNS=(runs/U_phase/um1836_470eV
          runs/U_phase/um1836_47keV)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-1"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="kmu"
    ;;
  wall_smoke)
    # The reflecting-wall go/no-go. PASSED 2026-09-05 (drift to 0.01%, reflection
    # 0 -> 0.17, both ramps off the wall at +11.5 d_i0/t*). Kept for re-running after any
    # deck change. ~6 GPU-min; use --time 00:20:00 so it backfills.
    RUNS=(runs/W_phase/ws_smoke_470eV)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-0"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="kwall"
    ;;
  wall47k)
    # The CLEAN-SIDE control, ~0.22 GPU-h. Submitted SEPARATELY from ws_470eV and with a
    # short wall on purpose: the two differ 30x in cost, and a 9 h request queues as badly
    # as the expensive one (2026-09-04: a 4 h request waited 4 h while 20 min backfilled in
    # minutes). Use --time 01:00:00.
    RUNS=(runs/W_phase/ws_47keV)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-0"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="kw47k"
    ;;
  wall470_t14)
    # The t* = 1.4 RERUN, ~24 GPU-h -- the project's most expensive single run. ws_470eV
    # stopped pre-shock at t* = 0.70 because its clock was taken from the piston driver;
    # the 47 keV control only steepens at t* ~ 1.4. Single GPU on gpu_shared (2-day cap, so
    # no checkpointing) and deliberately the SAME launch path ws_470eV completed on, not
    # the untested 4-GPU one. Use --time 36:00:00: shared bills elapsed, not requested, so
    # headroom is free except for queue position, and a timeout here costs 24 GPU-h.
    RUNS=(runs/W_phase/ws_470eV_t14)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-0"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="kw470t14"
    ;;
  wall470)
    # The comparison itself, ~6.7 GPU-h (2,783,015 steps at 65,592 cells). Give it real
    # headroom -- um1836_470eV timed out at 73% on a 6 h wall. Use --time 10:00:00.
    RUNS=(runs/W_phase/ws_470eV)
    BINARY="$SWEEP_BUILD"
    ARRAY="0-0"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_work"
    JOBNAME="kw470"
    ;;
  ab)
    # The wall A/B. Two runs per binary because WarpX on GPU is NOT reproducible --
    # ablastr's RandomSeed.H says so outright ("one should not expect to obtain the same
    # random numbers, even if a fixed random_seed is provided") and two runs of one deck
    # confirmed it on chablis. Without a same-binary replicate there is no noise floor and
    # |A-B| cannot be interpreted at all.
    #
    # A1/B1 land in the repo run dirs; A2/B2 are labelled replicates in $WORKROOT, because
    # a replicate is the SAME config and must not become a second config.yaml.
    RUNS=(runs/S_phase/ss_dz1_ppc100
          runs/S_phase/ss_dz1_ppc100:A2
          runs/S_phase/ss_dz1_ppc100_symwall
          runs/S_phase/ss_dz1_ppc100_symwall:B2)
    ARRAY="0-3"
    WORKROOT="${PSCRATCH:-$HOME}/kinshock_ab"
    JOBNAME="kab"
    # tasks 0,1 use binary A; tasks 2,3 use binary B -- job.sbatch reads BINARY, so the
    # array is split into two submissions rather than smuggling the mapping into the body.
    ;;
  *) echo "usage: $0 {sweep|eps|mag|mu|wall_smoke|wall47k|wall470|wall470_t14|ab} [run_dir ...] [--dry]" >&2; exit 2 ;;
esac

if [[ ${#SUBSET[@]} -gt 0 ]]; then
    # `ab` splits its array across two binaries by index, so a subset would silently
    # remap which binary each task gets. Refuse rather than mis-tag the runs.
    if [[ "$WHAT" == "ab" ]]; then
        echo "submit: run-dir subsets are not supported for 'ab' (the array index selects the binary)" >&2
        exit 2
    fi
    for r in "${SUBSET[@]}"; do
        [[ -f "$KINSHOCK_ROOT/${r%%:*}/config.yaml" ]] || {
            echo "submit: no config.yaml under '$r' -- refusing to submit a path that does not exist" >&2
            exit 2; }
    done
    RUNS=("${SUBSET[@]}")
    ARRAY="0-$(( ${#RUNS[@]} - 1 ))"
    echo "submit: run-dir subset given -- ${#RUNS[@]} of the '$WHAT' target"
fi

# The run list goes to a FILE, not into --export. sbatch's --export is a comma-separated
# list, so a value containing spaces (or commas) is silently mangled -- passing
# RUNS="a b c" produced `--export=ALL,RUNS=a b c` and would have broken every submission.
# One line per array index; job.sbatch reads the line matching SLURM_ARRAY_TASK_ID.
mkdir -p "$KINSHOCK_PM/.runlists"
RUNLIST="$KINSHOCK_PM/.runlists/${WHAT}-$(date +%Y%m%dT%H%M%S).txt"
printf '%s\n' "${RUNS[@]}" > "$RUNLIST"
echo "run list -> $RUNLIST"
cat -n "$RUNLIST" | sed 's/^/    /'

submit() {   # submit <array> <binary> <name>
    local arr="$1" bin="$2" name="$3"
    # Slurm writes -o/-e relative to the SUBMISSION cwd, so the #SBATCH defaults in
    # job.sbatch dropped 72 slurm-*.out/.err into the repo root. Pass an absolute path
    # here instead, for the same reason -A/-q/-t/--array are passed here: #SBATCH
    # directives cannot read shell variables. Slurm does NOT create the directory and
    # the job fails at launch if it is missing, hence the mkdir.
    local logdir="$KINSHOCK_ROOT/slurm"
    mkdir -p "$logdir"
    local qos="${QOS_OVERRIDE:-${SWEEP_QOS:-shared}}"
    local wall="${TIME_OVERRIDE:-${SWEEP_TIME:-03:00:00}}"
    # debug is capped at 30 min AND 5 jobs per user (measured 2026-08-11) -- catch the
    # walltime mistake here rather than after the queue wait.
    if [[ "$qos" == "debug" && "$wall" > "00:30:00" ]]; then
        echo "submit: -q debug caps walltime at 00:30:00 (asked $wall)" >&2; return 2
    fi
    local cmd=(sbatch -A "$NERSC_ACCOUNT" -q "$qos"
               -t "$wall" --array="$arr" -J "$name"
               -o "$logdir/slurm-%A_%a.out" -e "$logdir/slurm-%A_%a.err"
               --export="ALL,KINSHOCK_PM=$KINSHOCK_PM,RUNLIST=$RUNLIST,BINARY=$bin,WORKROOT=$WORKROOT"
               "$KINSHOCK_PM/job.sbatch")
    echo "+ ${cmd[*]}"
    [[ $DRY -eq 1 ]] || "${cmd[@]}"
}

if [[ "$WHAT" == "ab" ]]; then
    submit "0-1" A "${JOBNAME}A"
    submit "2-3" B "${JOBNAME}B"
else
    submit "$ARRAY" "$BINARY" "$JOBNAME"
fi

[[ $DRY -eq 1 ]] && echo "(--dry: nothing submitted)"
exit 0
