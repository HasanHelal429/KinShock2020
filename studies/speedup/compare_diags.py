#!/usr/bin/env python
"""Compare two benchmark runs' reduced diagnostics column by column.

Used for the GPU-vs-CPU agreement check in lever3_gpu.sh. Bit-identity is not the
bar and never could be: the reductions run in a different order and the pairwise
Coulomb operator draws from a different RNG stream on the device. What we need to
rule out is a *systematic* divergence -- a wrong heating rate, a lost species, a
mis-scaled weight -- which shows up as a growing relative difference rather than
noise at a fixed level.

So this prints, per column, the relative difference at the FIRST and LAST recorded
step. If the last-step difference is not much larger than the first, the two runs
agree to the level the stochastics allow.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def load(run_dir, name):
    p = Path(run_dir) / "diags" / "reducedfiles" / f"{name}.txt"
    if not p.is_file():
        return None, None
    header = p.open().readline().lstrip("#").split()
    return header, np.loadtxt(p, comments="#", ndmin=2)


def reldiff(a, b):
    scale = max(abs(a), abs(b))
    return 0.0 if scale == 0.0 else abs(a - b) / scale


MIN_ROWS = 3          # fewer than this is not a trajectory, it is an initial condition


def compare(dir_a, dir_b, name):
    """Returns True if a real comparison happened, False if there was nothing to compare.

    The distinction matters more than it looks. Twice in one session this script reported
    a non-comparison in a way that read like a pass: once with the runs' reduced-diag
    interval (5000) larger than the run length (1500), giving "1 common rows, steps 0..0";
    once with a stale path, giving "missing in A -- skipped". Both printed the encouraging
    footer and exited 0. A check that cannot fail is not a check, so "nothing compared" is
    now an explicit FAILURE rather than a quiet skip.
    """
    ha, da = load(dir_a, name)
    hb, db = load(dir_b, name)
    if da is None or db is None:
        which = "A" if da is None else "B"
        print(f"  {name}: *** NOT COMPARED *** — no {name}.txt under run {which}. "
              f"This is not a pass; the run probably failed or never wrote reduced diags.")
        return False
    n = min(len(da), len(db))
    if n < MIN_ROWS:
        print(f"  {name}: *** NOT COMPARED *** — only {n} common row(s), need >= {MIN_ROWS}. "
              f"The reduced-diag interval is coarser than the run length, so this would "
              f"compare the initial condition and nothing else.")
        return False
    if ha != hb:
        print(f"  {name}: WARNING column headers differ")
    print(f"  {name}: {n} common rows, steps {da[0,0]:.0f}..{da[n-1,0]:.0f}")
    print(f"    {'column':<46}{'first':>12}{'last':>12}")
    grew = []
    for i, col in enumerate(ha):
        if i < 2:                      # step, time
            continue
        first = reldiff(da[0, i], db[0, i])
        last = reldiff(da[n - 1, i], db[n - 1, i])
        flag = "  <-- GROWING" if last > max(10 * first, 1e-3) else ""
        grew.append(bool(flag))
        print(f"    {col[:46]:<46}{first:12.3e}{last:12.3e}{flag}")
    return not any(grew)


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        print("usage: compare_diags.py <run_dir_A> <run_dir_B>")
        return 2
    a, b = argv[1], argv[2]
    print(f"A = {a}\nB = {b}")
    results = [compare(a, b, name) for name in ("PN", "EP")]
    if not all(results):
        print("\n  RESULT: **NOT A PASS.** Either nothing was compared, or a column is "
              "GROWING.\n  Do not use the GPU binary for physics on the strength of this "
              "run.")
        return 1
    print("\n  RESULT: pass — both diagnostics compared over a trajectory and no column\n"
          "  diverges. Residual differences at this level are the independent RNG streams\n"
          "  (per-species energy differences scale as 1/sqrt(N)), not a behavioural change.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
