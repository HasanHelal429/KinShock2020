"""Shared plotting / movie helpers and media paths for KinShock2020.

Centralizes figure styling, the media-folder layout (all output under
``KinShock2020/media/``, see REPLICATION_PLAN.md §6), and ffmpeg encoding, so the
driver scripts stay thin. Matplotlib uses the non-interactive Agg backend.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# consistent series colors (piston / ambient / reference), matching the reference script
C_PISTON = "#c7522a"
C_AMBIENT = "#1f6f8b"
C_REF = "#8a8a8a"

# bright emissive tints for the additive phase-space distribution (dark panel):
# piston = warm, ambient = cool; where both overlap the sum tends to white.
PHASE_COLORS = {"piston": (1.00, 0.52, 0.22), "ambient": (0.28, 0.72, 1.00)}
PHASE_BG = "#0b0b12"

# repo root = two levels up from this file (.../KinShock2020/src/kinshock/plotting.py)
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(_PKG_DIR))
MEDIA = os.path.join(ROOT, "media")

FFMPEG = (os.environ.get("FFMPEG")
          or shutil.which("ffmpeg")
          or os.path.join(sys.exec_prefix, "bin", "ffmpeg"))


def media_dir(run_id=None, testing=False):
    """Path under media/: ``media/testing`` (progress) or ``media/<run_id>`` (final)."""
    d = os.path.join(MEDIA, "testing") if testing else \
        (os.path.join(MEDIA, run_id) if run_id else MEDIA)
    os.makedirs(d, exist_ok=True)
    return d


def savefig(fig, name, run_id=None, testing=False, dpi=130):
    """Save ``fig`` to the appropriate media dir and return its path."""
    out = os.path.join(media_dir(run_id=run_id, testing=testing), name)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)
    return out


def style_axes(ax):
    ax.grid(True, alpha=0.25, linewidth=0.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _asinh_norm(H, softclip=0.03):
    """Map a 2D histogram to [0,1] with an asinh (log-like) stretch, so both the
    dense downstream core and the sparse reflected-ion beam are visible at once.
    ``softclip`` sets the fraction of the peak below which the scale is ~linear."""
    hmax = float(H.max()) if H.size else 0.0
    if hmax <= 0:
        return H
    a = softclip * hmax
    return np.arcsinh(H / a) / np.arcsinh(hmax / a)


def phase_distribution(ax, species_data, z_edges, v_edges, vline=1.0,
                       colors=None, bg=PHASE_BG, legend=False):
    """Draw a 2D phase-space *distribution* f(z, v) as an additive two-colour image.

    ``species_data`` maps a role key ('piston' / 'ambient') -> (z, v, weight) arrays
    already in plot units. Each species is weight-histogrammed on the shared
    (``z_edges``, ``v_edges``) grid and asinh-normalised **to its own peak** (so the
    ~250x-rarer ambient population is still visible), then tinted with its colour;
    where populations overlap the colours add toward white. Returns the RGB image.

    Replaces the old per-particle scatter: identical axes, but the local macro-
    particle density is now legible instead of a saturated dot cloud.
    """
    colors = colors or PHASE_COLORS
    nz, nv = len(z_edges) - 1, len(v_edges) - 1
    rgb = np.zeros((nv, nz, 3))
    present = []
    for key, (z, v, w) in species_data.items():
        if len(z) == 0:
            continue
        H, _, _ = np.histogram2d(z, v, bins=[z_edges, v_edges], weights=w)
        inten = _asinh_norm(H).T                      # (nv, nz)
        rgb += inten[..., None] * np.asarray(colors.get(key, (1, 1, 1)))
        present.append(key)
    rgb = np.clip(rgb, 0.0, 1.0)

    ax.imshow(rgb, origin="lower", aspect="auto", interpolation="nearest",
              extent=[z_edges[0], z_edges[-1], v_edges[0], v_edges[-1]])
    ax.set_facecolor(bg)
    ax.set_xlim(z_edges[0], z_edges[-1])
    ax.set_ylim(v_edges[0], v_edges[-1])
    if vline is not None:
        ax.axhline(vline, color="w", ls=":", lw=0.9, alpha=0.7)
        ax.text(0.985, vline, r"$v_z=v_{sh}$", transform=ax.get_yaxis_transform(),
                va="bottom", ha="right", fontsize=7, color="w", alpha=0.7)
    if legend:
        y = 0.965
        for key in present:
            ax.text(0.025, y, key, transform=ax.transAxes, va="top", ha="left",
                    fontsize=8, fontweight="bold", color=colors.get(key, (1, 1, 1)))
            y -= 0.085
    return rgb


def stamp(ax, cfg, scales, extra=None):
    """Annotate a figure with run_id + key parameters read from the config/scales,
    so every plot in media/ is self-describing (plan §6.0a)."""
    rid = cfg.get("meta", {}).get("run_id", "?")
    txt = (f"{rid}  |  M_A={scales.MA:.1f}  M_ms={scales.Mms:.1f}  "
           f"m_i/m_e={scales.mass_ratio:.0f}  n_t/n_e0={scales.nt/scales.namb:.0f}")
    if extra:
        txt += "  |  " + extra
    ax.text(0.0, 1.02, txt, transform=ax.transAxes, va="bottom", ha="left",
            fontsize=8, color=C_REF)


def encode(framedir, out, fps=8):
    """Encode frame_%03d.png in ``framedir`` to an mp4 at ``out`` (libx264)."""
    subprocess.run(
        [FFMPEG, "-y", "-loglevel", "error", "-framerate", str(fps),
         "-i", os.path.join(framedir, "frame_%03d.png"),
         "-c:v", "libx264", "-crf", "20", "-pix_fmt", "yuv420p", out],
        check=True,
    )
    print("wrote", out)
    return out
