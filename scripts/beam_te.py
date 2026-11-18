"""Test B add-on (requested by the test-A session): electron and ion temperatures in the BEAM
region, from the slab edge to the farthest ambient ion with v_z > 0.3 v_sh.

The region is cut into 0.1 d_i0 bins and each bin's variance is taken about its own mean flow,
so bulk-velocity gradients do not count as temperature. Reported per dump:
  - weighted mean T_e / T0 over the region (all electron species);
  - T_i,core / T0 (ambient ions with v_z < 0.3 v_sh);
  - their ratio, T_e / T_i,core;
  - T_e at the fast-ion front bin.
Usage: beam_te.py <run_dir> <out_dir> [t_max_wci0]
"""
import glob
import os
import sys

import numpy as np
import yaml
import yt

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
from kinshock import units  # noqa: E402

yt.set_log_level(40)
C, QE = 299792458.0, 1.602176634e-19
run, out = sys.argv[1], sys.argv[2]
tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
cfg = yaml.safe_load(open(os.path.join(run, "config.yaml")))
s = units.derive(cfg)
me, mi = s.me, s.mi
T0 = float(cfg["plasma"]["ambient"]["theta_0"]) * me * C**2 / QE
slab = float(cfg["geometry"]["slab_halfwidth_di"]) * s.di
dzb = 0.1 * s.di0


def load(ad, sp, m):
    z = np.asarray(ad[sp, "particle_position_x"])
    w = np.asarray(ad[sp, "particle_weight"])
    u = np.stack([np.asarray(ad[sp, f"particle_momentum_{c}"]) for c in "xyz"]) / (m * C)
    return z, w, u


def binned_T(z, w, u, m, lo, hi):
    """weighted-mean temperature [eV] over [lo, hi), each 0.1 d_i0 bin about its own mean."""
    edges = np.arange(lo, hi + dzb, dzb)
    idx = np.digitize(z, edges) - 1
    num = den = 0.0
    last = np.nan
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.sum() < 50:
            continue
        ws = w[sel]
        us = u[:, sel]
        mu = (us * ws).sum(1) / ws.sum()
        T = ((((us - mu[:, None]) ** 2) * ws).sum(1) / ws.sum()).mean() * m * C**2 / QE
        num += T * ws.sum()
        den += ws.sum()
        last = T
    return (num / den if den else np.nan), last


rows = []
for d in sorted(glob.glob(os.path.join(run, "diags", "diag1*")),
                key=lambda p: int(os.path.basename(p)[len("diag1"):])):
    ds = yt.load(d)
    t = float(ds.current_time)
    tw = t * s.wci0
    if tw > tmax + 1e-6 or tw < 0.05:
        continue
    ad = ds.all_data()
    zi, wi, ui = load(ad, "amb_ions", mi)
    vz = ui[2] * C / np.sqrt(1 + (ui**2).sum(0))
    fast = vz > 0.3 * s.vsh_model
    if fast.sum() < 50:
        continue
    zf = np.percentile(zi[fast], 99.5)
    lo, hi = slab, zf
    zs, ws, us = [], [], []
    for sp in ("amb_electrons", "piston_electrons"):
        z, w, u = load(ad, sp, me)
        zs.append(z); ws.append(w); us.append(u)
    ze, we, ue = np.concatenate(zs), np.concatenate(ws), np.concatenate(us, axis=1)
    Te, Te_front = binned_T(ze, we, ue, me, lo, hi)
    core = ~fast
    Ti, _ = binned_T(zi[core], wi[core], ui[:, core], mi, lo, hi)
    rows.append((tw, (zf - slab) / s.di0, Te / T0, Ti / T0, Te / Ti, Te_front / T0))
    print(os.path.basename(run), " ".join(f"{x:.4g}" for x in rows[-1]), flush=True)
os.makedirs(out, exist_ok=True)
np.savetxt(os.path.join(out, os.path.basename(run.rstrip("/")) + "_beam.csv"), np.array(rows),
           delimiter=",", header="t_wci0,beam_extent_di0,Te_beam/T0,Ti_core_beam/T0,Te/Ti,Te_front/T0")
