#!/usr/bin/env python
"""Is the pi-rotation symmetry wall different from specular reflection? Measured, not assumed.

    python scripts/ab_wall_test.py \
        --a runs/S_phase/ss_dz1_ppc100 <scratch>/A2 \
        --b runs/S_phase/ss_dz1_ppc100_symwall <scratch>/B2

WHY THIS EXISTS. `boundary: {lo: symmetry}` emits `boundary.reflect_symmetry_axis`, a
fork-only input that lived on an unmerged branch -- so every binary this project has used
parsed it, discarded it, and fell back to plain specular reflection. 27 runs asked for the
symmetry wall and none got it, and `symmetry` and `reflecting` were therefore the SAME
simulation (R0_half vs R0_half_sym differ in key_params by nothing). CLAUDE.md's "~5%
near-wall artifact" for the specular approximation has never been measured on this problem,
because the comparison that would have measured it was between two identical decks.

THE CONFIGS ARE IDENTICAL. The variable is the BINARY:
    A  build_cuda1d      feature/hybrid-laser @ acc2d6621
    B  build_cuda1d_sym  + cherry-picks d5f2e9917, 05d74af41

WHY EACH SIDE NEEDS TWO RUNS. WarpX on GPU is NOT reproducible -- ablastr's RandomSeed.H
says outright that "when GPU simulations are run, one should not expect to obtain the same
random numbers, even if a fixed random_seed is provided", and two runs of one deck confirm
it. So a single A-vs-B difference is meaningless on its own: it has to be read against the
run-to-run scatter of the SAME binary. That is the same measured-noise-floor method used to
validate the 2-GPU run against the 1-GPU run (RESULTS 2026-08-04).

READ IT AS: `|A-B| / floor`. Below ~1 the wall is doing nothing this run can see. Well
above ~1 it is real, and every one-sided result in the repo carries it as a systematic.

THE PISTON FRONT IS THE SENSITIVE METRIC. The wall sits at z = 0 where the piston is, and
the two BCs differ only in whether the transverse v_perp is flipped along with v_z. That is
a gyro-phase effect on exactly the particles the wall reflects, so if it shows up anywhere
it shows up in how the piston launches -- more than in the far upstream.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import kinshock                                       # noqa: E402
from kinshock import io as kio                        # noqa: E402
from check_domain_control import measure              # noqa: E402


def piston_front(run_dir, tmax, pct=99.5):
    """(t*wci0[], z_front/d_i0[]) from the piston ions -- the driver the wall acts on."""
    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    t, zf = [], []
    for p in kio.plotfiles(run_dir):
        fr = kio.load_frame(p)
        tt = fr.time * sc.wci0
        if tt > tmax:
            break
        try:
            z = np.asarray(fr.ds.all_data()[("piston_ions", "particle_position_x")])
        except Exception:
            continue
        if not len(z):
            continue
        t.append(tt)
        zf.append(float(np.percentile(z, pct)) / sc.di0)
    return np.array(t), np.array(zf)


def at(t, y, targets):
    """Linear interpolation onto the requested times; NaN outside the measured range."""
    out = []
    for tt in targets:
        out.append(float(np.interp(tt, t, y)) if len(t) > 1 and t.min() <= tt <= t.max()
                   else float("nan"))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", nargs="+", required=True, help="run dirs built with binary A")
    ap.add_argument("--b", nargs="+", required=True, help="run dirs built with binary B")
    ap.add_argument("--tmax", type=float, default=0.30)
    ap.add_argument("--zwindow", type=float, nargs=2, default=(8.0, 12.0))
    ap.add_argument("--edge", type=float, default=0.5)
    ap.add_argument("--at", type=float, nargs="*", default=[0.05, 0.10, 0.20, 0.29])
    args = ap.parse_args()

    if len(args.a) < 2:
        print("! --a needs at least TWO runs: without a replicate there is no noise floor "
              "and |A-B| cannot be interpreted. Refusing to report a bare difference.")
        return 1

    series = {}
    for tag, dirs in (("A", args.a), ("B", args.b)):
        for i, d in enumerate(dirs, 1):
            rid, t, ez, bp, sc = measure(d, args.tmax, *args.zwindow, args.edge, 0.2)
            tp, zf = piston_front(d, args.tmax)
            series[f"{tag}{i}"] = dict(dir=d, t=t, ez=ez, bp=bp, tp=tp, zf=zf)
            print(f"  {tag}{i}  {d}  ({len(t)} field frames, {len(tp)} particle frames)")

    metrics = [("rms Ez / vA B0", "t", "ez"),
               ("coherent |Bperp|/B0", "t", "bp"),
               ("piston front z/d_i0", "tp", "zf")]
    print(f"\nA/B wall test — |A-B| against the same-binary noise floor |A1-A2|\n")
    for name, tkey, ykey in metrics:
        vals = {k: at(v[tkey], v[ykey], args.at) for k, v in series.items()}
        A = np.nanmean([vals[k] for k in vals if k.startswith("A")], axis=0)
        B = np.nanmean([vals[k] for k in vals if k.startswith("B")], axis=0)
        floor = np.abs(vals["A1"] - vals["A2"])
        if "B2" in vals:
            floor = np.maximum(floor, np.abs(vals["B1"] - vals["B2"]))
        print(f"  {name}")
        print(f"    {'t*wci0':>8} {'A(mean)':>10} {'B(mean)':>10} {'|A-B|':>10} "
              f"{'floor':>10} {'ratio':>7}")
        for j, tt in enumerate(args.at):
            r = abs(A[j] - B[j]) / floor[j] if floor[j] > 0 else float("nan")
            print(f"    {tt:8.2f} {A[j]:10.4g} {B[j]:10.4g} {abs(A[j]-B[j]):10.4g} "
                  f"{floor[j]:10.4g} {r:7.2f}")
        print()
    print("  ratio < ~1: the wall is invisible at this configuration.")
    print("  ratio >> 1: it is real, and every one-sided run in the repo carries it.")
    print("  floor is the max over same-binary replicate pairs — GPU runs are not "
          "reproducible (ablastr RandomSeed.H), so this is irreducible.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
