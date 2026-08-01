#!/usr/bin/env python
"""Write ``runs/<ID>/README.md`` — the run's own page, with a source for every number.

    python scripts/make_run_readme.py runs/R1_warm
    python scripts/make_run_readme.py --all
    python scripts/make_run_readme.py --all --check    # CI: fail if any is stale

Hand-written prose between the ``<!-- prose:begin -->`` / ``<!-- prose:end -->`` markers
is preserved across regeneration; the tables are always rewritten from ``config.yaml``.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import kinshock  # noqa: E402
from kinshock import runreadme  # noqa: E402


def one(run_dir: str, check: bool = False) -> tuple[bool, str]:
    """Returns (changed, message)."""
    rid = os.path.basename(os.path.normpath(run_dir))
    cfg = kinshock.load(run_dir)
    text = runreadme.render(cfg, run_dir)
    path = os.path.join(run_dir, "README.md")
    old = open(path).read() if os.path.exists(path) else None
    if old == text:
        return False, f"{rid}: up to date"
    if check:
        return True, f"{rid}: STALE — run make_run_readme.py"
    with open(path, "w") as fh:
        fh.write(text)
    return True, f"{rid}: {'updated' if old else 'created'} {path}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="*")
    ap.add_argument("--all", action="store_true", help="every runs/*/ with a config.yaml")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if any README is stale")
    args = ap.parse_args()

    dirs = list(args.run_dir)
    if args.all:
        root = os.path.join(os.path.dirname(__file__), "..", "runs")
        dirs += sorted(os.path.dirname(p)
                       for p in glob.glob(os.path.join(root, "*", "config.yaml")))
    if not dirs:
        ap.error("give one or more run dirs, or --all")

    stale = 0
    for d in dirs:
        try:
            changed, msg = one(os.path.normpath(d), check=args.check)
        except Exception as exc:                      # keep going across runs
            changed, msg = True, f"{os.path.basename(d)}: ERROR {type(exc).__name__}: {exc}"
        stale += bool(changed and args.check)
        print(msg)
    return 1 if stale else 0


if __name__ == "__main__":
    raise SystemExit(main())
