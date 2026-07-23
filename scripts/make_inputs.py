#!/usr/bin/env python3
"""Generate a WarpX input deck from a run's config.yaml (the source of truth).

config.yaml holds only the intuitive primaries (densities/n0, theta=kT/m_e c^2,
lengths in d_e/d_i, speeds/c); this writes the corresponding WarpX deck so the
deck never has to be hand-edited — edit config.yaml and regenerate. See
REPLICATION_PLAN.md §6.0a and kinshock.deck.

Usage:
    python scripts/make_inputs.py runs/R1_core            # write runs/R1_core/inputs_kinshock_R1_core
    python scripts/make_inputs.py runs/R1_core --stdout   # print the deck, don't write
    python scripts/make_inputs.py runs/R1_core --check     # diff vs the existing deck, don't write
    python scripts/make_inputs.py runs/R1_core --verify    # after a run: check warpx_used_inputs vs config
    python scripts/make_inputs.py runs/R1_core -o path     # write to an explicit path
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import kinshock  # noqa: E402
from kinshock import deck  # noqa: E402


def _deck_path(cfg: dict, run_dir: str) -> str:
    name = cfg.get("meta", {}).get("deck") or f"inputs_kinshock_{cfg['meta']['run_id']}"
    return os.path.join(run_dir, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=os.path.join(ROOT, "runs", "R1"),
                    help="run directory containing config.yaml")
    ap.add_argument("-o", "--output", help="write the deck here instead of the default path")
    ap.add_argument("--stdout", action="store_true", help="print the deck to stdout, do not write")
    ap.add_argument("--check", action="store_true",
                    help="render and diff against the existing deck (parse+resolve); do not write")
    ap.add_argument("--verify", action="store_true",
                    help="compare warpx_used_inputs (post-run) against config; do not write")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    rid = cfg["meta"]["run_id"]

    if args.verify:
        used = os.path.join(args.run_dir, "warpx_used_inputs")
        if not os.path.isfile(used):
            sys.exit(f"no warpx_used_inputs in {args.run_dir} (run WarpX first)")
        warns = deck.verify(cfg, used)
        print(f"{rid}: warpx_used_inputs vs config —",
              "OK (WarpX ran exactly this config)" if not warns else "MISMATCH")
        for w in warns:
            print("  !", w)
        sys.exit(1 if warns else 0)

    text = deck.render(cfg)

    if args.stdout:
        sys.stdout.write(text)
        return

    if args.check:
        existing = args.output or _deck_path(cfg, args.run_dir)
        if not os.path.isfile(existing):
            sys.exit(f"no existing deck to check against: {existing}")
        warns = deck.verify(cfg, existing)
        print(f"{rid}: existing deck vs config —",
              "OK (physically equivalent)" if not warns else "DIFFERS")
        for w in warns:
            print("  !", w)
        sys.exit(1 if warns else 0)

    out = args.output or _deck_path(cfg, args.run_dir)
    with open(out, "w") as fh:
        fh.write(text)
    print(f"{rid}: wrote {out}")
    # self-check: the deck we just wrote must resolve back to the config
    warns = deck.verify(cfg, out)
    if warns:
        print("  ! WARNING: generated deck does not round-trip to config:")
        for w in warns:
            print("    -", w)
    else:
        print("  verified: deck resolves back to config primaries")


if __name__ == "__main__":
    main()
