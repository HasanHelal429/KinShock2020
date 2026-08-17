#!/bin/bash -l
# Submit a KinShock2020 job array on Perlmutter.
#
#   perlmutter/submit.sh sweep          # the six S_phase points, one GPU each
#   perlmutter/submit.sh eps            # the four E_phase eps-ladder rungs, one GPU each
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
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry)  DRY=1; shift ;;
        --qos)  QOS_OVERRIDE="$2"; shift 2 ;;
        --time) TIME_OVERRIDE="$2"; shift 2 ;;
        *) echo "unknown option '$1'" >&2; exit 2 ;;
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
  *) echo "usage: $0 {sweep|eps|ab} [--dry]" >&2; exit 2 ;;
esac

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
    local qos="${QOS_OVERRIDE:-${SWEEP_QOS:-shared}}"
    local wall="${TIME_OVERRIDE:-${SWEEP_TIME:-03:00:00}}"
    # debug is capped at 30 min AND 5 jobs per user (measured 2026-08-11) -- catch the
    # walltime mistake here rather than after the queue wait.
    if [[ "$qos" == "debug" && "$wall" > "00:30:00" ]]; then
        echo "submit: -q debug caps walltime at 00:30:00 (asked $wall)" >&2; return 2
    fi
    local cmd=(sbatch -A "$NERSC_ACCOUNT" -q "$qos"
               -t "$wall" --array="$arr" -J "$name"
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
