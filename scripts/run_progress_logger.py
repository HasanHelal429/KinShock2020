#!/usr/bin/env python3
"""Sidecar progress logger for a WarpX run — real wall-clock checkpoints.

Watches a run's ``run.log`` and appends a checkpoint line to ``<run_dir>/progress.log``
every N percent of ``max_step``. Each line reports real (wall-clock) elapsed time and
ETA, the WarpX compute rate, a contention factor (wall-rate / compute-rate, >1 when the
machine is shared), and the system load — so compute cost is trackable after the fact
and runs can be paced/scheduled without babysitting them.

Launch it right after starting WarpX (it waits for run.log to appear):

    # foreground WarpX writing run.log, logger in the background:
    OMP_NUM_THREADS=8 warpx.1d inputs > runs/R1_phase/R1_core/run.log 2>&1 &
    python scripts/run_progress_logger.py runs/R1_phase/R1_core &

    python scripts/run_progress_logger.py runs/R1_phase/R1_core --every-pct 5 --poll 20
    python scripts/run_progress_logger.py runs/R1_phase/R1_core --total 125000 --log run.log

--total is auto-detected from warpx_used_inputs / the input deck (last ``max_step``
wins, matching ParmParse), so appended overrides are respected. The logger stops when
the run reaches max_step, prints its end marker, or goes stale.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import time

STEP_RE = re.compile(r"STEP (\d+) ends")
AVG_RE = re.compile(r"Avg\. per step = ([0-9.eE+-]+)")
EVOLVE_RE = re.compile(r"Evolve time = ([0-9.eE+-]+)")
MAXSTEP_RE = re.compile(r"^\s*max_step\s*=\s*(\d+)", re.MULTILINE)
END_MARKERS = ("AMReX finalized", "Total Time")


def detect_total(run_dir: str) -> int | None:
    """Last max_step across warpx_used_inputs and any input deck (ParmParse last-wins)."""
    candidates = [os.path.join(run_dir, "warpx_used_inputs")]
    candidates += sorted(glob.glob(os.path.join(run_dir, "inputs*")))
    last = None
    for path in candidates:
        try:
            with open(path) as fh:
                hits = MAXSTEP_RE.findall(fh.read())
            if hits:
                last = int(hits[-1])
        except OSError:
            continue
    return last


def read_state(log_path: str):
    """(current_step, warpx_avg_s_per_step, cumulative_evolve_s, ended?) from run.log."""
    try:
        with open(log_path, errors="replace") as fh:
            txt = fh.read()
    except OSError:
        return None
    steps = STEP_RE.findall(txt)
    step = int(steps[-1]) if steps else 0
    avg = float(AVG_RE.findall(txt)[-1]) if AVG_RE.search(txt) else float("nan")
    ev = float(EVOLVE_RE.findall(txt)[-1]) if EVOLVE_RE.search(txt) else float("nan")
    ended = any(m in txt for m in END_MARKERS)
    return step, avg, ev, ended


def fmt_dur(sec: float) -> str:
    sec = int(max(sec, 0))
    return f"{sec // 3600}h{(sec % 3600) // 60:02d}m"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="run directory containing run.log")
    ap.add_argument("--log", default="run.log", help="log filename inside run_dir (default run.log)")
    ap.add_argument("--out", default="progress.log", help="output filename inside run_dir")
    ap.add_argument("--total", type=int, default=None, help="total steps (default: auto-detect max_step)")
    ap.add_argument("--every-pct", type=float, default=10.0, help="checkpoint every this many %% (default 10)")
    ap.add_argument("--poll", type=float, default=30.0, help="seconds between polls (default 30)")
    ap.add_argument("--stale-min", type=float, default=15.0,
                    help="give up if run.log stops advancing this long (default 15 min)")
    args = ap.parse_args()

    log_path = os.path.join(args.run_dir, args.log)
    out_path = os.path.join(args.run_dir, args.out)
    total = args.total or detect_total(args.run_dir)
    if not total:
        raise SystemExit(f"could not determine total steps; pass --total (looked in {args.run_dir})")

    def loadavg() -> float:
        try:
            return float(open("/proc/loadavg").read().split()[0])
        except OSError:
            return float("nan")

    def emit(line: str):
        with open(out_path, "a") as fh:
            fh.write(line + "\n")

    # wait for run.log
    t_wait = time.time()
    while not os.path.isfile(log_path):
        if time.time() - t_wait > args.stale_min * 60:
            raise SystemExit(f"run.log never appeared at {log_path}")
        time.sleep(args.poll)

    header = (
        f"# WarpX progress — {os.path.abspath(args.run_dir)}\n"
        f"# logging_started {time.strftime('%Y-%m-%dT%H:%M:%S')}  total_steps={total}  "
        f"every={args.every_pct:g}%  poll={args.poll:g}s\n"
        f"# wall_time | step/total (pct) | wall_elapsed | wall_rate | warpx_rate | "
        f"ETA_wall | contention | load"
    )
    emit(header)

    t0 = time.time()
    last_pct_bucket = -1
    prev_wall, prev_step = t0, 0
    last_advance_wall, last_seen_step = t0, 0

    while True:
        time.sleep(args.poll)
        st = read_state(log_path)
        if st is None:
            continue
        step, avg, ev, ended = st
        now = time.time()

        if step > last_seen_step:
            last_seen_step, last_advance_wall = step, now
        stale = (now - last_advance_wall) > args.stale_min * 60

        pct = 100.0 * step / total
        bucket = int(pct // args.every_pct)
        done = step >= total or ended
        checkpoint = bucket > last_pct_bucket or done or stale

        if checkpoint and step > 0:
            d_wall = now - prev_wall
            d_step = step - prev_step
            wall_rate = d_wall / d_step if d_step > 0 else float("nan")
            eta = (total - step) * wall_rate if d_step > 0 else float("nan")
            cont = wall_rate / avg if (avg == avg and avg > 0) else float("nan")
            line = (
                f"{time.strftime('%Y-%m-%dT%H:%M:%S')} | {step}/{total} ({pct:4.1f}%) | "
                f"{fmt_dur(now - t0)} | {wall_rate:.4f} s/step(wall) | {avg:.4f} s/step(warpx) | "
                f"ETA ~{fmt_dur(eta)} | x{cont:.2f} | load {loadavg():.1f}"
            )
            emit(line)
            last_pct_bucket = bucket
            prev_wall, prev_step = now, step

        if done or stale:
            mean = (now - t0) / max(step, 1)
            status = ("reached max_step" if step >= total
                      else "end marker" if ended else "STALLED (no progress)")
            emit(f"DONE {time.strftime('%Y-%m-%dT%H:%M:%S')} | {step}/{total} "
                 f"({100.0*step/total:.1f}%) | total wall {fmt_dur(now - t0)} | "
                 f"mean {mean:.4f} s/step(wall) | {status}")
            break


if __name__ == "__main__":
    main()
