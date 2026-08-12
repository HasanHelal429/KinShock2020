# Running KinShock2020 on Perlmutter (NERSC)

> **Executed 2026-08-11.** This was written blind on chablis; it has since been run
> end-to-end on Perlmutter and needed no changes beyond one bug in `build_warpx.sh`'s
> self-check (fixed). Both binaries built, the wall A/B ran 4/4 clean, and the S_phase
> sweep followed. The ⚠ marks below were the author's unverifiable guesses — each one is
> now annotated with what actually happened. See RESULTS 2026-08-11 (Perlmutter).

## What Perlmutter actually buys here

Not per-GPU speed — and **measurement (2026-08-11) was harsher than the guess above it.**
The A100 turned in **479 s for the cost-1 point** (`ss_dz1_ppc100`: WarpX `Total Time`
483.1 s and 475.6 s over two runs), against the ~459 s/point the chablis estimate implied.
So the A100 was **not faster at all** here, let alone the 1.5–2.5× guessed — these 1D
decks are latency-bound, and one box on one GPU leaves an A100 badly under-occupied.

**The win is entirely concurrency.** All six S_phase points run at once, so wall-clock is
the *longest single run* — `ss_dz4_ppc100` at 16 cost units ≈ **2 h 09 m** — instead of the
31-unit, ~4.1 h serial total. Note the original "~1 h" here was too optimistic by 2×: it
overlooked that the longest point is 16× the baseline, not ~8×.

Scaling that holds: cost = `n_cell × ppc × max_step`, i.e. the README table's 1/2/4/4/4/16,
at **479 s/unit** plus ~43 s/job of Slurm + logger overhead. One caveat — `ss_dz4_ppc25`
has 16× the reference's *cell*-steps but only 4× its *particle*-steps, so it is
under-costed by this model and may run ~40–50 min rather than 33.

## ~~Blocker: the cherry-picks are local-only~~ — RESOLVED 2026-08-11

The cherry-picks were pushed to the shared lab branch, and `feature/hybrid-laser`
fast-forwarded `acc2d6621 → fcb48c9fe`. **That tip is `WARPX_COMMIT_B` verbatim** — the
rebuilt commit reproduced the chablis SHA exactly, so `site.conf.example` needed no
change and `build_warpx.sh both` works as written. Nothing to do here.

⚠ **Do not build B from `origin/feature/reflect-symmetry-axis`.** `d5f2e9917` and
`05d74af41` are on that branch too, but it forked *before* the heater merge and carries
**zero** `ParticleHeater`/`TargetInjector` files (15 commits behind `development`). B must
come from `feature/hybrid-laser`. Checked out and confirmed 2026-08-11.

## One-time setup

```bash
# 1. profile (upstream's, edit line 2: proj must end in _g)
cp $PSCRATCH/warpx-cda/Tools/machines/perlmutter-nersc/perlmutter_gpu_warpx.profile.example \
   $HOME/perlmutter_gpu_warpx.profile
vi $HOME/perlmutter_gpu_warpx.profile      # set proj
source $HOME/perlmutter_gpu_warpx.profile

# 2. dependencies (boost, c-blosc, adios2, blaspp, lapackpp). Slow, idempotent, once.
bash $PSCRATCH/warpx-cda/Tools/machines/perlmutter-nersc/install_gpu_dependencies.sh

# 3. both repos on $PSCRATCH — NOT $HOME. Run dirs receive the plotfiles, and $HOME is
#    quota'd and not for parallel I/O.
cd $PSCRATCH
git clone git@github.com:Schaeffer-Lab/warpx-cda.git
git clone https://github.com/HasanHelal429/KinShock2020.git

# 4. site config
cd $PSCRATCH/KinShock2020
cp perlmutter/site.conf.example perlmutter/site.conf
vi perlmutter/site.conf                    # NERSC_ACCOUNT, paths, SWEEP_BUILD
```

⚠ Lustre striping: WarpX plotfiles are many small files. If write time shows up in
`progress.log`, `stripe_medium $PSCRATCH/KinShock2020/runs` before the first run.

## Build

```bash
perlmutter/build_warpx.sh both     # A and B; or just the one SWEEP_BUILD names
```

The script checks out a **detached commit**, not a branch, so the binary's provenance is a
SHA. It also greps the binary for `reflect_symmetry_axis` and says which wall you got —
the check whose absence let 27 runs silently use the wrong boundary.

## Submit

```bash
perlmutter/submit.sh sweep --dry   # print the sbatch line, submit nothing
perlmutter/submit.sh sweep         # six S_phase points, one GPU each
perlmutter/submit.sh ab            # the wall A/B, two runs per binary
```

⚠ **QOS.** `site.conf` defaults to `SWEEP_QOS=shared`, which bills per-GPU and suits
these single-GPU runs. If `shared` is not available for GPU jobs on your allocation, set
`SWEEP_QOS=regular` — that takes whole 4-GPU nodes per task and wastes three of them, so
prefer `shared` if it works. Confirm with `sacctmgr show assoc user=$USER format=qos%40`.

Each array task is one run, ordered cheapest-first so a mistake surfaces on the 8-minute
run rather than the 100-minute one.

**Why single-GPU per run and not one big job:** every S_phase config sets
`numerics.max_grid_size = n_cell`, i.e. one box, and one box cannot be split across ranks.
`scripts/launch.sh` refuses that on chablis for the same reason. More GPUs would idle.

## Read the results

```bash
# the domain control -- must pass before the other five mean anything
python scripts/check_domain_control.py runs/S_phase/ss_dz1_ppc100 \
       runs/R1_phase/R1_paper_470eV --tmax 0.30      # needs the reference diags present

# the wall A/B
python scripts/ab_wall_test.py \
   --a runs/S_phase/ss_dz1_ppc100          $PSCRATCH/kinshock_ab/A2 \
   --b runs/S_phase/ss_dz1_ppc100_symwall  $PSCRATCH/kinshock_ab/B2

# the sweep itself
python scripts/plot_ez.py runs/S_phase/<ID> --tmax 0.3
python scripts/compare_phase.py runs/S_phase/ss_dz1_ppc100 runs/S_phase/ss_dz4_ppc25 \
       --at 0.11 0.21 0.30
```

⚠ `check_domain_control.py` compares against `R1_paper_470eV`, whose ~18 GB of diags live
on chablis (and `/mnt/cellar`). Either copy the `diag_fields*` directories for
steps ≤ 149090 (~60 frames, a few hundred MB) to Perlmutter, or run that one comparison on
chablis. The other analyses need only Perlmutter-local data.

## Two things to keep straight

**The A/B needs replicates, and that is not optional.** WarpX on GPU is not reproducible —
`Source/ablastr/math/RandomSeed.H` states that "when GPU simulations are run, one should
not expect to obtain the same random numbers, even if a fixed random_seed is provided",
and two runs of one deck on chablis confirmed it. So `submit.sh ab` runs each binary
**twice**; `ab_wall_test.py` refuses to report a bare `|A−B|` without a same-binary
replicate to measure the floor against. `A2`/`B2` are labelled scratch runs under
`$PSCRATCH/kinshock_ab`, deliberately *not* second copies of the config — the config is
identical and duplicating it would misrepresent what is being varied.

**Which binary the sweep should use.** If `SWEEP_BUILD=B`, the sweep uses the correct
π-rotation wall, and its `ss_dz1_ppc100` no longer matches `R1_paper_470eV` (specular) —
so the domain control has to be re-read in light of whatever the A/B says the wall is
worth. If `SWEEP_BUILD=A`, the sweep stays directly comparable to every existing result
and the A/B stands on its own. Both are defensible; the A/B is what tells you how much
the choice matters, so run `ab` first — it is four 8-minute runs.
