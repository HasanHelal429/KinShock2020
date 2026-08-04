#!/usr/bin/env python
"""Migrate ``field.vA_over_c`` -> ``field.B0_tesla`` in one or more run configs.

B0 was previously derived as ``vA_over_c * c * sqrt(mu0 * namb * m_i)``. This applies
exactly that map and rewrites the key in place, so the resulting deck is numerically
identical and ``make_inputs.py --verify`` still passes against an existing
``warpx_used_inputs``.

    python scripts/migrate_field_b0.py runs/R1_phase/R1_warm [runs/... ]
    python scripts/migrate_field_b0.py --all
    python scripts/migrate_field_b0.py --all -n      # dry run

The edit is line-based rather than a YAML round-trip so that comments, key order and
formatting survive untouched.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import yaml  # noqa: E402

from kinshock import config, units  # noqa: E402

_KEY = re.compile(r"^(\s*)vA_over_c\s*:\s*([0-9eE.+-]+)\s*(#.*)?$")


def migrate(run_dir: str, dry: bool = False) -> str:
    path = os.path.join(run_dir, "config.yaml")
    if not os.path.exists(path):
        return f"{run_dir}: no config.yaml — skipped"
    with open(path) as fh:
        text = fh.read()
    cfg = yaml.safe_load(text)
    fld = cfg.get("field", {})
    if "B0_tesla" in fld:
        return f"{os.path.basename(run_dir)}: already migrated"
    if "vA_over_c" not in fld:
        return f"{os.path.basename(run_dir)}: no field.vA_over_c — skipped"

    # namb and m_i come from the same config, so the map is self-contained.
    n0 = float(cfg["reference"]["n0"])
    mass_ratio = float(cfg["reference"]["mass_ratio"])
    namb = float(cfg["plasma"]["ambient"]["density_over_n0"]) * n0
    mi = mass_ratio * units.ME
    vA_over_c = float(fld["vA_over_c"])
    B0 = units.b0_from_vA_over_c(vA_over_c, namb, mi)

    out, hit = [], False
    for line in text.splitlines(keepends=True):
        m = _KEY.match(line.rstrip("\n"))
        if m and not hit:
            hit = True
            indent, comment = m.group(1), (m.group(3) or "")
            nl = "\n" if line.endswith("\n") else ""
            out.append(f"{indent}B0_tesla: {B0!r}"
                       f"        # PRIMARY [T]. Migrated from vA_over_c = {vA_over_c!r}"
                       f" (identical B0).{nl}")
            out.append(f"{indent}# v_A is DERIVED: B0/sqrt(mu0*namb*m_i) ="
                       f" {B0 / (units.MU0 * namb * mi) ** 0.5 / units.C:.6g} c{nl}")
            if comment.strip():
                out.append(f"{indent}# (was: {comment.strip()}){nl}")
        else:
            out.append(line)
    if not hit:
        return f"{os.path.basename(run_dir)}: vA_over_c present in YAML but not matched"

    if dry:
        return f"{os.path.basename(run_dir)}: would set B0_tesla = {B0:.6g} T"
    with open(path, "w") as fh:
        fh.write("".join(out))
    return f"{os.path.basename(run_dir)}: B0_tesla = {B0:.6g} T  (was vA/c = {vA_over_c})"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="*")
    ap.add_argument("--all", action="store_true",
                    help="every run dir under runs/ (runs/<phase>_phase/<run_id>/)")
    ap.add_argument("-n", "--dry-run", action="store_true")
    args = ap.parse_args()

    dirs = list(args.run_dir)
    if args.all:
        root = os.path.join(os.path.dirname(__file__), "..", "runs")
        dirs += config.find_runs(root)
        stray = config.unphased_runs(root)
        if stray:
            print("note: not filed under a <phase>_phase/ folder: "
                  + ", ".join(os.path.basename(d) for d in stray), file=sys.stderr)
    if not dirs:
        ap.error("give one or more run dirs, or --all")
    for d in dirs:
        print(migrate(os.path.normpath(d), dry=args.dry_run))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
