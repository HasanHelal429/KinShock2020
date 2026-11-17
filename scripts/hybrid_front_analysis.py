#!/usr/bin/env python3
"""Scalar comparison of hybrid_phase runs: piston-front speed, electron temperature, and
the dimensionless identity of a 470 eV / 47 keV twin pair.

    python scripts/hybrid_front_analysis.py runs/hybrid_phase/hc_470eV_t2 \
        runs/hybrid_phase/hy_470eV_t2 [--pair runs/hybrid_phase/hc_470eV_t2 runs/hybrid_phase/hc_47keV_t2]

Piston front = the weight-weighted 99th percentile of piston-ion z (the definition RESULTS
2026-10-12 used for the kinetic 12.1-13.7 and hybrid 6.2-9.2 d_i0/t* numbers), fitted
linearly over t*w_ci0 windows 0.1-0.5, 0.3-1.0, 0.5-2.0. T_e is the hybrid fluid field,
reported in units of the run's own ambient T_0 so twins at different temperatures compare
directly: the slab maximum, and the value at the piston front.

--pair: two runs on the same dimensionless grid and clock (same steps); reports the
relative L2 difference of B_x/B0, n_i/n_amb and T_e/T_0 at matched steps. The hybrid has
no lambda_D or c, so for a pair built to be the same dimensionless problem this measures
numerics, not physics.

Shock features (false shock, notch) are judged BY EYE from the movies, not here.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

from kinshock import deck  # noqa: E402
from kinshock import io as kio  # noqa: E402
from hybrid_phase_movie import dumps, scales_from_used_inputs  # noqa: E402

QE = 1.602176634e-19
WINDOWS = ((0.1, 0.5), (0.3, 1.0), (0.5, 2.0))


def te0_ev(run_dir):
    d = deck.parse_inputs(os.path.join(run_dir, "warpx_used_inputs"))
    ns = deck.resolve_constants(d, allow_unresolved=True)
    return float(ns["Te0"])


def fronts(run_dir, s):
    out = []
    for st, p in sorted(dumps(run_dir, "diag_phase").items()):
        fr = kio.load_frame(p, add_applied_field=False)
        ad = fr.ds.all_data()
        z = np.asarray(ad[("piston_ions", "particle_position_x")])
        w = np.asarray(ad[("piston_ions", "particle_weight")])
        if z.size == 0:
            continue
        o = np.argsort(z)
        cw = np.cumsum(w[o]) / w.sum()
        zf = z[o][np.searchsorted(cw, 0.99)]
        t = fr.time * s["wci0"]
        out.append((st, t, zf / s["di"]))
    return np.array(out)


def field_frames(run_dir, s):
    fr = {}
    for st, p in sorted(dumps(run_dir, "diag_fields").items()):
        f = kio.load_frame(p, add_applied_field=False,
                           fields=("Te", "rho_piston_ions", "rho_amb_ions"))
        zero = np.zeros_like(f.Bx)
        fr[st] = dict(t=f.time * s["wci0"], z=f.z_centers / s["di"], bx=f.Bx / s["B0"],
                      te=f.comps.get("Te", zero),
                      ni=(f.comps.get("rho_piston_ions", zero) + f.comps.get("rho_amb_ions", zero))
                      / (QE * s["namb"]))
    return fr


def speeds(fr):
    rows = []
    for lo, hi in WINDOWS:
        m = (fr[:, 1] >= lo - 1e-9) & (fr[:, 1] <= hi + 1e-9)
        rows.append(np.polyfit(fr[m, 1], fr[m, 2], 1)[0] if m.sum() >= 3 else np.nan)
    return rows


def report(run_dir):
    s = scales_from_used_inputs(run_dir)
    t0 = te0_ev(run_dir)
    fr = fronts(run_dir, s)
    ff = field_frames(run_dir, s)
    v = speeds(fr)
    print(f"\n== {os.path.basename(run_dir)}  (T_0 = {t0:.4g} eV; front = 99th pct piston-ion z)")
    print("  front speed [d_i0/t*]  " + "  ".join(f"t* {lo}-{hi}: {x:6.2f}"
                                               for (lo, hi), x in zip(WINDOWS, v)))
    print("     t*   z_front/d_i0   max T_e/T_0   T_e(front)/T_0")
    fsteps = np.array(sorted(ff))
    for st, t, zf in fr[:: max(1, len(fr) // 10)]:
        f = ff[int(fsteps[max(0, np.searchsorted(fsteps, st, side='right') - 1)])]
        i = int(np.argmin(np.abs(f["z"] - zf)))
        print(f"  {t:5.2f}   {zf:10.3f}   {f['te'].max() / t0:11.1f}   {f['te'][i] / t0:12.1f}")
    return dict(s=s, t0=t0, ff=ff, speeds=v)


def pair(a, b):
    print(f"\n== identity check {os.path.basename(a['dir'])} vs {os.path.basename(b['dir'])} "
          "(relative L2 at matched steps)")
    common = sorted(set(a["ff"]) & set(b["ff"]))
    print("     t*     B_x/B0    n_i/n_amb   T_e/T_0")
    for st in common[:: max(1, len(common) // 10)]:
        fa, fb = a["ff"][st], b["ff"][st]
        r = []
        for k, na, nb in (("bx", 1, 1), ("ni", 1, 1), ("te", a["t0"], b["t0"])):
            xa, xb = fa[k] / na, fb[k] / nb
            r.append(np.linalg.norm(xa - xb) / max(np.linalg.norm(xa), 1e-300))
        print(f"  {fa['t']:5.2f}   " + "   ".join(f"{x:9.3e}" for x in r))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--pair", nargs=2, default=None)
    args = ap.parse_args()
    done = {}
    for r in args.runs + (args.pair or []):
        r = os.path.abspath(r)
        if r not in done:
            done[r] = report(r)
            done[r]["dir"] = r
    if args.pair:
        pair(done[os.path.abspath(args.pair[0])], done[os.path.abspath(args.pair[1])])


if __name__ == "__main__":
    main()
