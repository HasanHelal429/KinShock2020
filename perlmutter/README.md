# Running KinShock2020 on Perlmutter (NERSC)

> **Executed 2026-08-11.** Written blind on chablis, then verified and run on the real
> machine. Verified from chablis first: NERSC reachable with the existing ssh cert, account
> `m5032_g`, `$PSCRATCH=/pscratch/sd/h/hhelal`, and the QOS table below via
> `sbatch --test-only`. Then executed end-to-end on Perlmutter — both binaries built, the
> wall A/B ran 4/4 clean, and the S_phase sweep followed — needing no changes beyond one bug
> in `build_warpx.sh`'s self-check (fixed). The ⚠ marks below were the author's unverifiable
> guesses; each is now annotated with what actually happened. The build environment is not
> invented; it comes from upstream `Tools/machines/perlmutter-nersc/` in the `warpx-cda`
> tree. See RESULTS 2026-08-11 (Perlmutter).

## What Perlmutter actually buys here

Both, as it turns out — the guess above was about right, and **the whole sweep was measured
on 2026-08-11** (job `56715249`). Concurrency is still the bigger win.

## ~~Blocker: the cherry-picks are local-only~~ — RESOLVED 2026-08-11

`feature/hybrid-laser` was advanced `acc2d6621..fcb48c9fe` on `Schaeffer-Lab/warpx-cda` — a
clean fast-forward of two commits, 8 files, +232/−0, all additive. The rebuilt commit
**reproduced the chablis SHA exactly**, so `site.conf.example` needed no change and
`build_warpx.sh both` works as written. A fresh Perlmutter clone can build **both**:

| tag | commit | wall |
|---|---|---|
| A | `acc2d6621` | specular — what every existing result used |
| B | `fcb48c9fe` | π-rotation — what the configs have always asked for |

Both were built on Perlmutter 2026-08-11 and confirmed by `strings`: A carries 0 hits of
`reflect_symmetry_axis`, B carries 12.

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

**QOS — verified on Perlmutter 2026-08-11** with `sbatch --test-only` (validates and
estimates, submits nothing). `shared` works for GPU jobs and is by far the right choice:

| QOS | partition | resources | max wall | max jobs | would start |
|---|---|---|---|---|---|
| `shared` | `shared_gpu_ss11` | 1 GPU / 32 cores | 2 d | 5000 | **+11 h** |
| `debug` | `gpu_ss11` | whole node | **30 min** | **5** | **+5 h** |
| `regular` | `gpu_ss11` | whole node | 2 d | 5000 | **+6 days** |

`regular` is six days deep and wastes three of four GPUs per task — do not use it for this.
`--test-only` reports an upper bound, so real starts may be sooner, but the ordering is
the point.

⚠ **The multi-GPU production decks have no choice but `regular`.** `gpu_shared` is capped
at `gres/gpu=2` (`sacctmgr show qos`), so anything above two GPUs — e.g.
`R1_paper_470eV_ppc400` at `max_grid_size = n_cell/4`, four boxes on four GPUs — must go to
`regular` and be billed a whole node. Weigh that against the six-day queue depth above: for
a ~30 h run, going 2 → 4 GPUs saves order 13 h of *runtime* and can cost days of *queue*.
**Measured 2026-08-12, and it is not close:**

| config | QOS | would start | wait |
|---|---|---|---|
| 4 GPU | `regular` | 2026-08-17 | **~5.8 days** |
| 2 GPU | `shared` | 2026-08-12 | **~6.5 h** |

Time-to-result was ~36 h on two GPUs against ~6.5 days on four, so `R1_paper_470eV_ppc400`
runs on **two**. Shortening the walltime ask does not help `regular` — 24 h and 48 h return
the same date. **Measure the queue with `--test-only` before choosing a resource shape**;
it is free, submits nothing, and the runtime win from more GPUs is easily swamped by it.
The estimate moves with load, so re-measure rather than trusting this table.

**Put the A/B in `debug`.** It is four ~8-minute runs, which is exactly within debug's
5-job and 30-minute limits, and it starts about 6 h sooner:

```bash
perlmutter/submit.sh ab --qos debug --time 00:30:00
```

`submit.sh` refuses `--qos debug` with a walltime over 30 min rather than letting you
discover it after the wait. The sweep cannot use debug — `ss_dz4_ppc25` (~33 min) and
`ss_dz4_ppc100` (~1 h 42 m) exceed the cap, and six tasks exceed the job limit — so it
goes to `shared`.

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
