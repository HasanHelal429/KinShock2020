"""kinshock — reusable analysis library for the KinShock2020 project.

Reproduction and verification of the Schaeffer et al. 2020 piston-driven magnetized
collisionless shocks (Phys. Plasmas 27, 042901) using WarpX's ported laser-ablation
surrogate. See ../../REPLICATION_PLAN.md.

Design: per-run ``config.yaml`` holds the PRIMARY parameters (single source of truth);
:mod:`kinshock.units` derives all physical/normalized scales from them; the analysis
modules hold no hard-coded physical constants.

Import layers (lightest first):
  * :mod:`kinshock.units`   — pure Python (math only)
  * :mod:`kinshock.deck`    — pure Python; config -> WarpX deck generation + verification
  * :mod:`kinshock.config`  — + PyYAML
  * :mod:`kinshock.metrics` — + numpy
  * :mod:`kinshock.io`      — + yt (lazy)
  * :mod:`kinshock.plotting`— + matplotlib
"""

from __future__ import annotations

from . import units, deck, config, metrics  # light deps; safe to import eagerly

__all__ = ["units", "deck", "config", "metrics", "io", "plotting", "load",
           "find_runs", "unphased_runs"]
__version__ = "0.1.0"

# convenience re-exports
load = config.load
find_runs = config.find_runs            # runs/<phase>_phase/<run_id>/ discovery
unphased_runs = config.unphased_runs


def __getattr__(name):
    # lazily expose io (yt) and plotting (matplotlib) without importing them at package load
    if name in ("io", "plotting"):
        import importlib
        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
