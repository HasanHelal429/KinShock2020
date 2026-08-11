# Running KinShock2020 on Perlmutter (NERSC)

> **None of this has been executed.** It was written on chablis, which has no Perlmutter
> access, so treat the first submission as a shakeout. Everything marked ⚠ below is a
> specific thing I could not verify from here. The build environment is not invented —
> it comes from upstream `Tools/machines/perlmutter-nersc/` in the `warpx-cda` tree.

## What Perlmutter actually buys here

Not per-GPU speed. An A100 is maybe 1.5–2.5× an RTX 4070 on these decks, and the small
points are latency-bound rather than bandwidth-bound, so the low end is likelier. **The
win is concurrency**: all six S_phase points run at once, so wall-clock is the *longest
single run* (~1 h) instead of the 3 h 37 m serial total measured on chablis.

## Blocker: the cherry-picks are local-only

`fcb48c9fe` (= `acc2d6621` + `d5f2e9917` + `05d74af41`, the `reflect_symmetry_axis`
fix) exists **only in the chablis clone**. Perlmutter clones from
`git@github.com:Schaeffer-Lab/warpx-cda.git` and will not find it. Before building
binary B:

```bash
# on chablis, after you're satisfied with the cherry-pick
cd ~/warpx-cda && git push origin feature/hybrid-laser
```

That pushes to a **shared lab branch**, so it is your call, not a mechanical step. If you
would rather not move `feature/hybrid-laser`, push it as its own branch
(`git push origin feature/hybrid-laser:build/symwall`) and point `WARPX_COMMIT_B` at that.

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
