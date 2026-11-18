"""Plot the heating ledger (test B) and print a summary table."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

D = "/global/homes/h/hhelal/.claude/jobs/32b8a20b/tmp/ledger"
OUT = "/pscratch/sd/h/hhelal/KinShock2020/media/F_phase/heating_ledger.png"
RUNS = [  # id, label, forms?
    ("R1_470eV_s3f8_t2", "470 eV prod (dz/λD 0.76)", True),
    ("fs_470eV_dz0379_t1", "470 eV fine (0.38)", True),
    ("fs_470eV_T100_t1", "470 eV, T0 100 eV", True),
    ("ft_15keV_t2", "15 keV", False),
    ("ft_30keV_t2", "30 keV", False),
    ("ek_47keV_t2", "47 keV", False),
    ("fs_47keV_dz1p5", "47 keV coarse (1.5)", False),
]
warm = iter(["#b2182b", "#ef8a62", "#d6604d"])
cool = iter(["#2166ac", "#4393c3", "#1b7837", "#5aae61"])
fig, ax = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
cols = [(2, "upstream T_e / T_0"), (3, "upstream T_i / T_0 (core v_z < 0.3 v_sh)"),
        (5, "slab T_e / T_e,ab"), (6, "slab T_i / T_e,ab")]
print("run | verdict | Te_up/T0 @0.3 | @last clean | Ti_up/T0 @0.3 | @last clean | last clean t (wpe t) | beams arrive | slabTe/Tab @end | "
      "slabTi/Tab @end | t_end | wpe_ab t_end")
for rid, lab, forms in RUNS:
    p = os.path.join(D, rid + ".csv")
    if not os.path.exists(p):
        continue
    a = np.atleast_2d(np.loadtxt(p, delimiter=","))
    c = next(warm) if forms else next(cool)
    ls = "-" if forms else "--"
    ok = a[:, 4] < 0.01   # upstream window free of fast ions (beam fraction < 1%)
    arr = a[~ok, 0].min() if (~ok).any() else None
    ok = a[:, 0] < arr if arr is not None else ok   # nothing after the first arrival
    for k, (j, ttl) in enumerate(cols):
        x = ax.flat[[0, 1, 3, 4][k]]
        m = ok if j in (2, 3) else np.ones(len(a), bool)
        x.plot(a[m, 0], a[m, j], ls, color=c, lw=1.6, label=f"{lab} ({'forms' if forms else 'clean'})"
               + (f", beams arrive {arr:.2f}" if arr is not None else ""))
        x.set_title(ttl)
        x.set_xlabel("t·ω_ci0")
    ax[0, 2].plot(a[ok, 1], a[ok, 2], ls, color=c, lw=1.6)
    ax[1, 2].plot(a[ok, 1], a[ok, 3], ls, color=c, lw=1.6)
    i3 = np.argmin(abs(a[:, 0] - 0.3))
    ie = np.where(ok)[0][-1]
    print(f"{rid} | {'forms' if forms else 'clean'} | {a[i3,2]:.3f} | {a[ie,2]:.3f} | {a[i3,3]:.3f} | "
          f"{a[ie,3]:.3f} | {a[ie,0]:.2f} ({a[ie,1]:.3g}) | {arr} | {a[-1,5]:.3f} | {a[-1,6]:.3f} | {a[-1,0]:.2f} | {a[-1,1]:.3g}")
for x in (ax[0, 2], ax[1, 2]):
    x.set_xscale("log")
    x.set_xlabel("t·ω_pe,ab (plasma periods ×2π)")
ax[0, 2].set_title("upstream T_e / T_0 vs plasma time")
ax[1, 2].set_title("upstream T_i / T_0 vs plasma time")
ax[0, 0].legend(fontsize=8)
fig.suptitle("Test B — heating ledger (particle dumps). Upstream curves stop when fast ions (v_z > 0.3 v_sh) exceed 1% of the window. Far upstream: z ∈ [L−4, L−2.5] d_i0; "
             "slab: z ∈ [0, 2 d_i,ab]", fontsize=11)
fig.savefig(OUT, dpi=130)
print("wrote", OUT)
