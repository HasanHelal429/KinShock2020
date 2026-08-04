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


def compare(dir_a, dir_b, name):
    ha, da = load(dir_a, name)
    hb, db = load(dir_b, name)
    if da is None or db is None:
        print(f"  {name}: missing in {'A' if da is None else 'B'} -- skipped")
        return
    n = min(len(da), len(db))
    if n == 0:
        print(f"  {name}: no rows")
        return
    if ha != hb:
        print(f"  {name}: WARNING column headers differ")
    print(f"  {name}: {n} common rows, steps {da[0,0]:.0f}..{da[n-1,0]:.0f}")
    print(f"    {'column':<46}{'first':>12}{'last':>12}")
    for i, col in enumerate(ha):
        if i < 2:                      # step, time
            continue
        first = reldiff(da[0, i], db[0, i])
        last = reldiff(da[n - 1, i], db[n - 1, i])
        flag = "  <-- GROWING" if last > max(10 * first, 1e-3) else ""
        print(f"    {col[:46]:<46}{first:12.3e}{last:12.3e}{flag}")


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        print("usage: compare_diags.py <run_dir_A> <run_dir_B>")
        return 2
    a, b = argv[1], argv[2]
    print(f"A = {a}\nB = {b}")
    for name in ("PN", "EP"):
        compare(a, b, name)
    print("\n  A 'first' column that is already ~1e-16 and a 'last' column still small\n"
          "  means the runs agree. Anything marked GROWING needs explaining before the\n"
          "  GPU binary is used for physics.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
