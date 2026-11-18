"""Test B (FALSE_SHOCK_AUDIT): heating ledger from particle plotfiles.

Per dump: T_e, T_i (3-component variance of u = p/(m c), times m c^2, about the window mean) in
  - far upstream: z in [L - 4, L - 2.5] d_i0, ambient species; ions restricted to v_z < 0.3 v_sh
    (the ambient core; the fraction excluded is reported as beam contamination);
  - slab: z in [0, slab half-width], piston species.
Writes <out>/<run>.csv.  Usage: heat_ledger.py <run_dir> <out_dir> [t_max_wci0]
"""
import glob
import math
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
me = s.me
mi = s.mi
T0 = float(cfg["plasma"]["ambient"]["theta_0"]) * me * C**2 / QE
Tab = float(cfg["plasma"]["piston"]["theta_e_heat"]) * me * C**2 / QE
slab = float(cfg["geometry"]["slab_halfwidth_di"]) * s.di
L = s.domain_halfwidth * (1 if cfg["geometry"].get("layout") == "one_sided" else 2)
zlo, zhi = L - 4 * s.di0, L - 2.5 * s.di0
wpe_ab = s.wpe


def temp(ds, sp, m, zmin, zmax, vcut=None):
    ad = ds.all_data()
    z = np.asarray(ad[sp, "particle_position_x"])
    sel = (z >= zmin) & (z < zmax)
    if not sel.any():
        return math.nan, 0.0
    w = np.asarray(ad[sp, "particle_weight"])[sel]
    u = np.stack([np.asarray(ad[sp, f"particle_momentum_{c}"])[sel] for c in "xyz"]) / (m * C)
    frac = 0.0
    if vcut is not None:
        g = np.sqrt(1 + (u**2).sum(0))
        core = u[2] / g * C < vcut
        frac = 1 - w[core].sum() / w.sum()
        u, w = u[:, core], w[core]
    mu = (u * w).sum(1) / w.sum()
    var = (((u - mu[:, None]) ** 2) * w).sum(1) / w.sum()
    return var.mean() * m * C**2 / QE, frac


rows = []
for d in sorted(glob.glob(os.path.join(run, "diags", "diag1*")), key=lambda p: int(os.path.basename(p)[len("diag1"):])):
    ds = yt.load(d)
    t = float(ds.current_time)
    if t * s.wci0 > tmax + 1e-6:
        continue
    Te_up, _ = temp(ds, "amb_electrons", me, zlo, zhi)
    Ti_up, fbeam = temp(ds, "amb_ions", mi, zlo, zhi, vcut=0.3 * s.vsh_model)
    Te_sl, _ = temp(ds, "piston_electrons", me, 0.0, slab)
    Ti_sl, _ = temp(ds, "piston_ions", mi, 0.0, slab)
    rows.append((t * s.wci0, t * wpe_ab, Te_up / T0, Ti_up / T0, fbeam, Te_sl / Tab, Ti_sl / Tab))
    print(os.path.basename(run), f"t*wci0={t*s.wci0:.3f}", " ".join(f"{x:.4g}" for x in rows[-1][2:]),
          flush=True)
os.makedirs(out, exist_ok=True)
np.savetxt(os.path.join(out, os.path.basename(run.rstrip("/")) + ".csv"), np.array(rows),
           delimiter=",", header="t_wci0,t_wpe_ab,Te_up/T0,Ti_up/T0,beam_frac_up,Te_slab/Tab,Ti_slab/Tab")
