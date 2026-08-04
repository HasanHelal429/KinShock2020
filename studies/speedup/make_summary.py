#!/usr/bin/env python
"""Generate ``runs/opt_phase/SUMMARY.md`` from the benchmark points, and file the
per-point artifacts alongside it.

Written as a generator rather than a hand-authored page for the same reason the run
READMEs are: the numbers should come from the measurements, not from whatever I
remembered while typing. Re-run it after any new point and the table updates.

    python studies/speedup/make_summary.py            # write SUMMARY.md
    python studies/speedup/make_summary.py --sync      # + copy artifacts into runs/opt_phase/

WHICH STATISTIC. The mean, not the median. bench.sh reports both, and for these runs the
mean is the honest one: the only outliers in a clean run are the collision supercycle
(`ndt_supercycle: 10`), so every 10th step is genuinely expensive and you pay for it. The
median hides that. The exception is a point contaminated by external I/O, where the mean is
inflated by work that is not yours -- those points are marked and superseded by a clean
re-run rather than silently kept.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "studies" / "speedup" / "out"
DEST = ROOT / "runs" / "opt_phase"

STEPS_FULL = 3_224_046          # R1_paper_470eV max_step = 220 t_ab
DRIFT = 1.279                   # particle-count growth over the full run

RESULT = re.compile(
    r"^(\S+)\s+thr=\s*(\d+)\s+n=\s*(\d+)\s+s/step=\s*([0-9.]+)\s+median=\s*([0-9.]+)"
    r"\s+status=(\S+)")

# label -> (lever, human description, trust)
NOTES = {
    "l1_thr4":            ("1", "4 OMP threads", "ok"),
    "l1_thr8":            ("1", "8 OMP threads (pilot baseline)", "ok"),
    "l1_thr12":           ("1", "12 OMP threads", "ok"),
    "l1_thr16":           ("1", "16 OMP threads", "ok"),
    "l1_thr20":           ("1", "20 OMP threads", "CONTAMINATED"),
    "l1_thr24":           ("1", "24 OMP threads", "CONTAMINATED"),
    "l1_thr20_clean":     ("1", "20 OMP threads (clean re-run)", "ok"),
    "l1_thr24_clean":     ("1", "24 OMP threads (clean re-run)", "ok"),
    "l1_thr16_mgs64":     ("1", "16 threads, max_grid_size=64", "ok"),
    "l1_thr24_mgs64":     ("1", "24 threads, max_grid_size=64", "ok"),
    "l2A_smoke_picard":   ("2", "theta-implicit, Picard, cfl 0.75 (smoke)", "ok"),
    "l2B_time_picard":    ("2", "theta-implicit, Picard, cfl 0.75", "ok"),
    "l2C_newton_cfl3.0":  ("2", "Newton, cfl 3.0 -- MISCONFIGURED", "BROKEN"),
    "l2N1_broken_verbose": ("2", "Newton, reproduce the misconfig with verbosity", "diagnostic"),
    "l2N2_skipinit":      ("2", "Newton + skip_particle_picard_init", "ok"),
    "l2N3_pc_jacobi":     ("2", "Newton + skip_init, pc_jacobi", "ok"),
    "l2N4_no_pc":         ("2", "Newton + skip_init, no preconditioner", "ok"),
    "l3A_gpu_default":    ("3", "GPU, deck default (235 boxes)", "ok"),
    "l3B_gpu_1box":       ("3", "GPU, 1 box (max_grid_size=30000)", "ok"),
    "l3B2_gpu_8box":      ("3", "GPU, 8 boxes (max_grid_size=4096)", "ok"),
    "l3C_cpu_ref":        ("3", "CPU reference, 8 threads", "ok"),
    "l3D_gpu_agree":      ("3", "GPU, agreement-check run", "ok"),
    "l3D_cpu_agree":      ("3", "CPU, agreement-check run", "ok"),
}


def read_points():
    pts = {}
    for d in sorted(OUT.iterdir()) if OUT.is_dir() else []:
        f = d / "result.txt"
        if not f.is_file():
            continue
        m = RESULT.match(f.read_text().strip())
        if not m:
            continue
        label, thr, n, mean, med, status = m.groups()
        pts[label] = dict(label=label, thr=int(thr), n=int(n), mean=float(mean),
                          median=float(med), status=status,
                          lever=NOTES.get(label, ("?", label, "ok"))[0],
                          desc=NOTES.get(label, ("?", label, "ok"))[1],
                          trust=NOTES.get(label, ("?", label, "ok"))[2])
    return pts


def table(pts, lever, base):
    hdr = (f"| point | config | thr | steps | s/step (mean) | median | vs base | full run |\n"
           f"|---|---|---|---|---|---|---|---|\n")
    rows = ""
    for p in sorted(pts.values(), key=lambda x: (x["lever"], x["mean"])):
        if p["lever"] != lever:
            continue
        spd = base / p["mean"] if p["mean"] else 0.0
        days = STEPS_FULL * p["mean"] * DRIFT / 86400.0
        flag = "" if p["trust"] == "ok" else f" **[{p['trust']}]**"
        dur = f"{days:.2f} d" if days < 400 else f"{days/365:.1f} yr"
        rows += (f"| `{p['label']}`{flag} | {p['desc']} | {p['thr']} | {p['n']} | "
                 f"**{p['mean']:.5f}** | {p['median']:.5f} | {spd:.2f}x | {dur} |\n")
    return hdr + rows if rows else "_no points yet_\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sync", action="store_true",
                    help="copy per-point artifacts into runs/opt_phase/")
    args = ap.parse_args()

    pts = read_points()
    base = pts["l1_thr8"]["mean"] if "l1_thr8" in pts else 0.11169

    DEST.mkdir(parents=True, exist_ok=True)
    body = f"""# Optimization sweep — `R1_paper_470eV`

Generated by `studies/speedup/make_summary.py` from `studies/speedup/out/*/result.txt`.
Do not hand-edit the tables. Harness and rationale live in `studies/speedup/README.md`;
findings and their interpretation are in `RESULTS.md`.

These are **performance points, not physics runs** — they carry no `config.yaml`, which is
why `kinshock.find_runs()` does not enumerate them.

Baseline: **{base:.5f} s/step** at 8 OMP threads. Full run is {STEPS_FULL:,} steps; the
"full run" column applies the measured {DRIFT} particle-growth drift, so it is a projected
wall-clock, not steps x s/step.

**Statistic: the mean.** In a clean run the only outliers are the collision supercycle
(`ndt_supercycle: 10`) — every 10th step is genuinely expensive and you pay for it, so the
median understates real throughput. Points marked CONTAMINATED had their means inflated by
*my own* filesystem work (a `find` over 130 GB, `git add -A` across the tree) during the
runs//media/ regrouping; their medians were clean, which is the tell. They are superseded
by `_clean` re-runs and kept only so the artefact is on the record.

## Lever 1 — OMP threads

{table(pts, "1", base)}
Monotonic to 24 threads with returns flattening after 20 (20 -> 24 buys ~5%). **No cliff**,
which is conditional: the earlier ~20x collapse above 12 threads was oversubscription
against other users on a busy box with a 5x smaller deck. `max_grid_size=64` is
consistently worse — the deck's default 235 grids is already right for CPU.

## Lever 2 — `algo.evolve_scheme = theta_implicit_em`

{table(pts, "2", base)}
Picard converges cleanly (21-23 iterations, every step exiting on relative tolerance) but
costs ~13.5x per step, so break-even against explicit needs dt to grow past cfl ~10. It is
**not a speedup**. What it does buy is exact energy conservation at theta = 0.5, i.e. no
grid heating by construction — the thing actually blocking this run.

The `l2C` point is a misconfiguration of mine, kept as evidence: `ImplicitSolver.cpp:811`
gates the cheap particle path on `use_mass_matrices_jacobian && skip_particle_picard_init`,
and I set only the first, so every Newton iteration ran a full 21-iteration particle Picard
update. The `l2N*` points retune it with verbosity on.

## Lever 3 — GPU

{table(pts, "3", base)}
**The decomposition matters more than the hardware.** The deck's CPU-tuned 235 boxes give
~1.2x; one box gives ~7.9x — a 6.5x swing from `amr.max_grid_size=30000` alone. A GPU
benchmark left at the default decomposition would have read as "not worth it".

Speed is contingent on the agreement check (`l3D_*`): the CUDA binary is dated Jul 28 and
commit `9f981dea2` landed Jul 31 fixing an nvcc rejection in `ParticleHeater`, so it
predates a fix to an operator this deck needs. It runs without aborting, which is necessary
and not sufficient.

## Reproducing

```bash
studies/speedup/lever1_threads.sh
studies/speedup/lever2_implicit.sh all
studies/speedup/lever2_newton_retune.sh
studies/speedup/lever3_gpu.sh
python studies/speedup/make_summary.py --sync
```
"""
    (DEST / "SUMMARY.md").write_text(body)
    print(f"wrote {DEST / 'SUMMARY.md'} ({len(pts)} points)")

    if args.sync:
        n = 0
        for label in pts:
            src, dst = OUT / label, DEST / label
            dst.mkdir(parents=True, exist_ok=True)
            for name in ("result.txt", "meta.txt"):
                if (src / name).is_file():
                    shutil.copy2(src / name, dst / name)
            rf = src / "diags" / "reducedfiles"
            if rf.is_dir():
                shutil.copytree(rf, dst / "reducedfiles", dirs_exist_ok=True)
            # run.log can carry hundreds of MB of solver verbosity. Keep only the two parts
            # with provenance value: the opening banner (Grids Summary, deposition, pusher,
            # gathering, shape factor) and the tail (TinyProfiler table). 150 lines each --
            # 800 was mostly per-step noise repeated identically across 17 near-identical
            # points, which is volume without information.
            log = src / "run.log"
            if log.is_file():
                lines = log.read_text(errors="replace").splitlines()
                head, tail = lines[:150], lines[-150:]
                (dst / "run.log.excerpt").write_text(
                    "\n".join(head)
                    + f"\n\n... [{max(0, len(lines) - 300)} lines elided] ...\n\n"
                    + "\n".join(tail) + "\n")
            n += 1
        print(f"synced {n} point directories into {DEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
