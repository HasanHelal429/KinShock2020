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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# consistent series colors (piston / ambient / reference), matching the reference script
C_PISTON = "#c7522a"
C_AMBIENT = "#1f6f8b"
C_REF = "#8a8a8a"

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
