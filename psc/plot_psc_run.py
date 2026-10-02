"""Overview of a psc_schaeffer2020 run: ion (z, u_z) phase space with B_y overlaid at
selected t*w_ci0, ablation theta_e / n_e vs time, and an optional movie.

    python psc/plot_psc_run.py <run_dir> [--movie]

z >= 0 only (the paper's half), in d_i0 = 100 d_e,ab; time in t*w_ci0 (1/w_ci0 =
1e4 / w_pe). Reads ps.*.bin (read_ps.py) and pfd.*.h5 (fields).
"""
import argparse, glob, os, re, subprocess, sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_ps import load_header, load_frame

DT, WCI_INV, DI0, B0 = 0.182595, 1.0e4, 100.0, 0.01

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("--movie", action="store_true")
args = ap.parse_args()
run = args.run
hdr = load_header(run)
steps = sorted(int(re.search(r"ps\.(\d+)\.bin", f).group(1)) for f in glob.glob(f"{run}/ps.*.bin"))


def by_profile(step):
    fn = f"{run}/pfd.{step:06d}_p000000.h5"
    if not os.path.exists(fn):
        return None, None
    f = h5py.File(fn, "r")
    c = f[[k for k in f if k.startswith("jeh")][0]]["hy_fc"]
    pats = sorted(c.keys(), key=lambda s: int(c[s].attrs["global_patch"][0]))
    by = np.concatenate([c[p]["3d"][:, :, 0].mean(axis=1) for p in pats])
    z = np.linspace(hdr["zlo"], hdr["zhi"], len(by) + 1)[:-1] + 0.5 * (hdr["zhi"] - hdr["zlo"]) / len(by)
    return z, by


def panel(ax, step, vmax=None):
    fr = load_frame(run, step, hdr, masses=[1, 100])
    ion = hdr["kinds"].index(next(k for k in hdr["kinds"] if k["name"] == "i"))
    H = fr["hist"][ion, 0].T
    zpos = fr["z"] >= 0
    ki = hdr["kinds"][ion]
    Hp = H[:, zpos]
    vmax = vmax or Hp.max()
    ax.imshow(np.where(Hp > 0, Hp, np.nan), origin="lower", aspect="auto", cmap="viridis",
              extent=[0, hdr["zhi"] / DI0, ki["ulo"], ki["uhi"]],
              norm=LogNorm(vmin=vmax * 1e-4, vmax=vmax))
    z, by = by_profile(step)
    if z is not None:
        ax2 = ax.twinx()
        m = z >= 0
        ax2.plot(z[m] / DI0, by[m] / B0, color="w", lw=0.8)
        ax2.set_ylim(-1, 9)
        ax2.set_ylabel("B_y / B0", color="0.4", fontsize=8)
        ax2.tick_params(labelsize=7)
    ax.set_ylabel("ion u_z = γβ_z", fontsize=8)
    ax.text(0.01, 0.92, f"t·ω_ci0 = {step * DT / WCI_INV:.2f}", transform=ax.transAxes,
            color="w", fontsize=9)
    ax.tick_params(labelsize=7)
    return vmax


# --- overview figure at selected times
want = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
picks = [min(steps, key=lambda s: abs(s * DT / WCI_INV - t)) for t in want]
fig, axs = plt.subplots(len(picks), 1, figsize=(10, 2.1 * len(picks)), sharex=True)
for ax, s in zip(axs, picks):
    panel(ax, s)
axs[-1].set_xlabel("z / d_i0")
axs[0].set_title("PSC (Schaeffer 2020 setup), ion phase space, z ≥ 0, with B_y/B0 (white)")
fig.tight_layout()
out = f"{run}/overview_ion_phase.png"
fig.savefig(out, dpi=110)
print("wrote", out)

# --- ablation plasma check vs time
t, th, ne = [], [], []
for s in steps:
    fr = load_frame(run, s, hdr, masses=[1, 100])
    m = (np.abs(fr["z"]) > 25) & (np.abs(fr["z"]) < 75)
    n = fr["n"][0][m]
    t.append(s * DT / WCI_INV)
    th.append(np.nansum(fr["T"][0][m] * n) / np.nansum(n))
    ne.append(n.mean())
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(t, th, label="θ_e, 25 < |z| < 75 d_e (pressure moment)")
ax.axhline(0.092, ls="--", c="k", lw=0.8, label="Table I θ_e,ab = 0.092")
ax.set_xlabel("t·ω_ci0"); ax.set_ylabel("θ_e"); ax.legend(loc="lower right", fontsize=8)
ax2 = ax.twinx(); ax2.plot(t, ne, c="gray", ls=":"); ax2.axhline(1.25, c="gray", lw=0.6)
ax2.set_ylabel("n_e there (gray; line = 1.25)")
fig.tight_layout(); fig.savefig(f"{run}/ablation_check.png", dpi=110)
print("wrote", f"{run}/ablation_check.png")
print(f"late-time (t·ω_ci0 > 1) θ_e = {np.mean([a for a, b in zip(th, t) if b > 1]):.4f}, "
      f"n_e = {np.mean([a for a, b in zip(ne, t) if b > 1]):.3f}")

# --- movie
if args.movie:
    fdir = f"{run}/movie_frames"
    os.makedirs(fdir, exist_ok=True)
    vmax = None
    for i, s in enumerate(steps):
        fig, ax = plt.subplots(figsize=(10, 3.4))
        vmax = panel(ax, s, vmax=None)
        ax.set_xlabel("z / d_i0")
        fig.tight_layout(); fig.savefig(f"{fdir}/f{i:04d}.png", dpi=100); plt.close(fig)
    ff = os.environ.get("FFMPEG", "ffmpeg")
    subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", "12", "-i", f"{fdir}/f%04d.png",
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                    f"{run}/ion_phase.mp4"], check=True)
    print("wrote", f"{run}/ion_phase.mp4")
