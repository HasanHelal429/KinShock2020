# `scripts/` — deck generation, bring-up & analysis drivers

Command-line entry points for the KinShock2020 replication of Schaeffer 2020.
Each script is config-driven: it reads `runs/<RUN_ID>/config.yaml`, derives all
physical scales through `kinshock.units`, and either generates the WarpX deck or
writes outputs under `media/<run_id>/` (or `media/testing/` for bring-up checks).

**`config.yaml` is the single source of truth.** You author the intuitive
primaries there (densities as fractions of n0, temperatures as θ = kT/m_e c²,
lengths in d_e/d_i, speeds as fractions of c); `make_inputs.py` generates the
WarpX input deck from it. The deck is a build artifact — never hand-edit it; edit
`config.yaml` and regenerate.

All scripts are run from the repository root:

```bash
python scripts/<script>.py [run_dir] [options]
```

The `run_dir` positional argument defaults to `runs/R1` for every script.

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `make_inputs.py` | Generate the WarpX input deck from the config | `config.yaml` | `inputs_kinshock_<id>` |
| `make_run_readme.py` | Write `runs/<ID>/README.md` — the run's own page, with a source for every number | `config.yaml` | `<run_dir>/README.md` |
| `migrate_field_b0.py` | One-shot `field.vA_over_c` → `field.B0_tesla` (exact map; decks keep the same numeric B0) | `config.yaml` | `config.yaml` (in place) |
| `launch.sh` | **The** way to start a run — cd's into the run dir so `diags/` lands there | `config.yaml`, deck | `<run_dir>/{run.log,diags/}` |
| `run_checks.py` | Bring-up / progress checks (works before any sim output exists) | `config.yaml`, plotfiles, reduced diags | `media/testing/*.png` |
| `make_figures.py` | Reproduce the paper's shock diagnostics (analyses A–D) | `config.yaml`, plotfiles | `media/<run_id>/*.png`, `criteria.json` |
| `make_movies.py` | Animated density + phase-space movies | `config.yaml`, plotfiles | `media/<run_id>/*.mp4` |
| `make_thomson.py` | Synthetic Thomson spectra: EPW + IAW spectrograms | `config.yaml`, plotfiles | `media/<run_id>/thomson_{epw,iaw}.png`, `thomson_spectra.npz` |
| `run_progress_logger.py` | Sidecar: real wall-clock progress/ETA log at %-checkpoints | `run.log`, input deck | `<run_dir>/progress.log` |
| `bfield_diagnostic.py` | B-field fluctuation: physical vs numerical (spectra, polarization, particle-response) | `config.yaml`, plotfiles | `media/<run_id>/bfield_diagnostic.png` |

Typical workflow for a run:

```bash
python scripts/make_inputs.py  runs/R1          # config.yaml -> inputs_kinshock_R1
scripts/launch.sh -b -L        runs/R1          # start WarpX + the progress logger
python scripts/make_inputs.py  runs/R1 --verify # confirm warpx_used_inputs == config
python scripts/run_checks.py   runs/R1          # sanity: scales vs Table I, conservation
python scripts/make_figures.py runs/R1          # A–D diagnostics + criteria table
python scripts/make_movies.py  runs/R1          # optional animations
```

---

## `launch.sh`

Starts a WarpX run. **Use this rather than invoking the binary yourself**: the generated deck sets
no `diag*.file_prefix`, so WarpX writes plotfiles to `diags/` *relative to the launch CWD*, and two
runs launched from the repo root will share `./diags/` and clobber each other (this cost a rerun of
the R2/R3 controls — see RESULTS 2026-07-26). `launch.sh` always `cd`s into the run dir first.

It also applies the benchmarked thread settings (`OMP_NUM_THREADS=8 OMP_PROC_BIND=spread
OMP_PLACES=cores` — near-linear to 8 cores, ~1.8× vs 4), resolves the run's single `inputs_*` deck,
sends stdout/stderr to `<run_dir>/run.log`, and **refuses to start when `diags/` already holds
output**, since relaunching overwrites a finished run's plotfiles in place.

```bash
scripts/launch.sh runs/R1                    # foreground, logs to runs/R1/run.log
scripts/launch.sh -b -L runs/R1              # detach + start the progress logger
scripts/launch.sh -n runs/R1                 # dry run: print the command, change nothing
scripts/launch.sh runs/R0 -- max_step=20     # ParmParse overrides (smoke tests)
```

| Flag | Default | Description |
|---|---|---|
| `-j, --threads N` | `8` | `OMP_NUM_THREADS`. The benchmarked sweet spot is 8; beyond that the run is memory-bandwidth-bound. |
| `-w, --warpx PATH` | `$KINSHOCK_WARPX`, else `/home/hhelal/warpx-cda/build/bin/warpx.1d` | WarpX binary. |
| `-b, --background` | off | Detach via `nohup` and return immediately, printing the PID. |
| `-L, --logger` | **on** | Start `run_progress_logger.py` in the background (stdout → `<run_dir>/logger.out`, gitignored). On by default — every run should leave a `progress.log`. |
| `--no-logger` | off | Suppress the progress logger. Rarely wanted: the run then leaves no wall-clock record. |
| `--every-pct P` | `10` | Logger checkpoint every P percent of `max_step`. On a 12 h run the default 10% is ~1.2 h per line — use `0.5`–`1` for a run you want to watch. |
| `--poll S` | `30` | Logger poll interval, seconds. |
| `-f, --force` | off | Launch even though `diags/` already has output, overwriting it. |
| `-n, --dry-run` | off | Print what would run and exit. Warns instead of failing on a populated `diags/`. |
| `-- <args>` | — | Everything after `--` is appended as ParmParse overrides. **Not** reflected in `config.yaml`, so `make_inputs.py --verify` will flag them afterwards — smoke tests only, never physics. |

---

## `make_inputs.py`

Generate a WarpX input deck from a run's `config.yaml` — the forward direction of
the config↔deck relationship (REPLICATION_PLAN.md §6.0a).

`config.yaml` holds only the primary, physically-intuitive parameters;
`kinshock.deck.render` maps them onto a WarpX deck whose `my_constants` are
written *symbolically* (`nt = 2.5*n0`, `slab = 2.0*di`,
`B0 = vA*sqrt(mu0*namb*Mi)`) so the deck stays readable and WarpX still records
the fully-resolved values in `warpx_used_inputs`. After writing, the script
self-checks that the deck resolves back to the config primaries.

**Positional argument**

| Argument | Required | Default | Description |
|---|---|---|---|
| `run_dir` | no | `runs/R1` | Run directory containing `config.yaml`. The deck is written to `<run_dir>/<meta.deck>` (default `inputs_kinshock_<run_id>`). |

**Options**

| Flag | Type | Default | Description |
|---|---|---|---|
| `-o, --output <path>` | path | `<run_dir>/<meta.deck>` | Write the deck to an explicit path. |
| `--stdout` | flag | off | Print the deck to stdout instead of writing it. |
| `--check` | flag | off | Render and diff against the existing deck (parse + resolve both, compare), without writing. Reports `OK (physically equivalent)` or lists differences. |
| `--verify` | flag | off | **Post-run:** parse `warpx_used_inputs` and confirm the numbers WarpX actually used match the config. Reports `OK (WarpX ran exactly this config)` or lists mismatches. Exits non-zero on mismatch. |

The comparison is numeric (resolves every expression), so it is immune to
formatting, comments, or whether a value was written as `20.*de` or `2.0*di`. It
value-checks only `my_constants` present in both decks (WarpX prunes unused ones
from `warpx_used_inputs`) plus the scalar settings (`max_step`, `n_cell`, `cfl`,
`tau`, per-species ppc, diagnostic intervals, …).

**Diagnostics (config `diagnostics:` block).** `plotfile_intervals` and
`reduced_intervals` set the cadence of the `Full` plotfiles (`diag1`: fields +
particles, for phase space) and the `EP`/`PN` reduced diags. Optionally,
`field_intervals` adds a second, **field-only** diagnostic (`diag_fields`,
`write_species = 0`) at that high cadence — full grid resolution, all of
`Ex Ey Ez Bx By Bz jx jy jz rho` plus per-species `rho`. Fields are ~1–2 MB/frame
(a few % of a particle frame), so this buys high-fidelity streaks / field
structure cheaply. Omit `field_intervals` to keep the single-diagnostic deck.
`make_figures.py` uses the dense `diag_fields` series for the B_perp streak when
present and falls back to `diag1` otherwise.

**Examples**

```bash
python scripts/make_inputs.py runs/R1_core            # write the deck (+ self-check)
python scripts/make_inputs.py runs/R1_core --stdout   # preview without writing
python scripts/make_inputs.py runs/R1_core --check    # is the existing deck still in sync?
python scripts/make_inputs.py runs/R1_core --verify   # after a run: did WarpX run this config?
```

---

## `run_checks.py`

Bring-up / progress checks for a run — writes to `media/testing/`.

Always emits a **config-summary figure** (derived scales vs Table I targets) so
progress is visible even before any WarpX output exists. If the run has produced
plotfiles / reduced diagnostics, it also emits:

- **loaded-state sanity** — initial per-species density and B_perp profiles, with
  the measured C_s,ab / v_A / M_A / M_ms / β_0 stamped against their Table I targets;
- **operator balance** — energy-conservation and piston-inventory histories
  (heater ↔ injector balance), from the `EP`/`PN` reduced diagnostics.

**Positional argument**

| Argument | Required | Default | Description |
|---|---|---|---|
| `run_dir` | no | `runs/R1` | Run directory to check (parsed positionally via `sys.argv`; no flags). |

**Example**

```bash
python scripts/run_checks.py runs/R1
```

Prints the derived-scales table and validation status to stdout. If no plotfiles
exist yet, it writes only the config-summary figure and exits.

---

## `make_figures.py`

Reproduce the Schaeffer 2020 shock diagnostics for a run. Config-driven; writes to
`media/<run_id>/`. Implements analyses A–D of REPLICATION_PLAN.md §6:

- **A.** B_perp(z,t) streak plot with piston/shock speed lines → `shock_streak.png`;
  shock-front trajectory + measured v_sh, M_A, M_ms → `shock_trajectory.png`
- **B.** ion (z, u_z) phase space at several times → `shock_phase.png`;
  density / B_perp line-outs → `shock_lineouts.png`;
  species-resolved 3-row phase space (ambient-ion / piston-ion / electron, Fig. 7)
  → `shock_fig7.png` (and `shock_fig7_rho_i0.png` with `--fig7-xunits rho_i0`)
- **C.** seven shock-formation criteria per frame → `criteria.json` (+ stdout summary)
- **D.** reflected-ambient-ion fraction G(t), F(z), t\*₁, z\*₁ → `shock_reflected.png`

**Positional argument**

| Argument | Required | Default | Description |
|---|---|---|---|
| `run_dir` | no | `runs/R1` | Run directory to analyze. |

**Options**

| Flag | Type | Default | Description |
|---|---|---|---|
| `--nframes <N>` | int | `5` | Number of time frames sampled for the multi-frame figures (line-outs and phase space). Frames are chosen evenly across the available plotfiles. |

**Example**

```bash
python scripts/make_figures.py runs/R1 --nframes 5
```

Prints derived scales, the frame count and time span, the measured v_sh / M_A /
M_ms, and the first-precursor / first-shock times from the criteria table.

**Fig. 7 velocity axes** are pinnable per row, in units of v_z/v_sh:

```bash
python scripts/make_figures.py runs/R1_paper --only fig7 \
    --fig7-xunits d_i0 rho_i0 --phase-times 0.15 0.49 0.73 0.98 1.25 \
    --v-ambient -0.5 2.0 --v-piston -0.2 1.7 --v-electron -6 12
```

Without them the two ion rows share `fig_phase`'s band (−1 … 3) and the electron row is
auto-sized **symmetrically** from a percentile of |v_z|/v_sh — which means it differs
between runs, so pin it whenever panels are to be compared run-to-run. The electron row's
twin axis (B_x/B₀, n_e/n_e0) keeps its zero aligned with v_z = 0 at any range, including
asymmetric ones.

---

## `make_movies.py`

Animated diagnostics for a run → `media/<run_id>/` (config-driven, yt + ffmpeg):

- `shock_ni.mp4` — ion density line-out n_i/n_e0(z) vs time
- `shock_phase.mp4` — ambient + piston ion (z, v_z/v_sh) phase space vs time
  (the v_z > v_sh dotted line marks the reflected-ion threshold)

Frames are rendered to `$MOVIE_SCRATCH` (default: session scratchpad) then encoded.
Uses the same shock speed the figures use — the model value from config, unless a
measured one is passed with `--vsh-c`.

**Positional argument**

| Argument | Required | Default | Description |
|---|---|---|---|
| `run_dir` | no | `runs/R1` | Run directory to animate. |

**Options**

| Flag | Type | Default | Description |
|---|---|---|---|
| `--fps <N>` | int | `8` | Frames per second of the encoded movies. |
| `--vsh-c <v>` | float | model value from config | Shock speed in units of c, used to normalize the phase-space v_z axis. Pass the value measured by `make_figures.py` for a self-consistent movie. |

**Environment**

| Variable | Default | Description |
|---|---|---|
| `MOVIE_SCRATCH` | session scratchpad dir | Directory where per-frame PNGs are rendered before ffmpeg encoding. |

**Example**

```bash
python scripts/make_movies.py runs/R1 --fps 8 --vsh-c 0.14
```


---

## `make_thomson.py`

Synthetic Thomson scattering spectra from a run's particle plotfiles → `media/<run_id>/`:

- `thomson_epw.png` — electron feature, Doppler-broadened by the electron thermal speed
- `thomson_iaw.png` — ion-acoustic feature, the ±k·C_s doublet
- `thomson_spectra.npz` — `t`, both spectrograms and wavelength axes, `alpha_*`, `n_e`

Each figure carries three panels: the absolute spectrogram, one normalised per timestep
(the scattered power climbs ~2 decades as the piston arrives, which otherwise saturates
the late frames), and α vs time.

Forward-modelling uses the Schaeffer PlasmaPy fork (branch `feature/pic-thomson-pipeline`),
found at `$KINSHOCK_PLASMAPY` or `~/Schaeffer_PlasmaPy/src`.

**Two stages, because no single env here has both dependencies** — the reader needs `yt`
(in `physics`) and the forward model needs `torch` (in `tsnn`). Binned phase spaces are
cached to `runs/<ID>/thomson_cache/` (gitignored) as the handoff:

```bash
python scripts/make_thomson.py runs/R1_paper                                  # reads  (physics)
/opt/anaconda3/envs/tsnn/bin/python scripts/make_thomson.py runs/R1_paper     # models (tsnn)
```

`--stage auto` (the default) does whichever half the current interpreter supports and
prints the command for the other. Delete `thomson_cache/` to force a re-read.

**Read the reported α before interpreting the IAW panel.** α = 1/(k λ_D) ∝ √n/T, so
collectivity is set by the run's absolute density: the same deck gives α ~ 1e-5 at
n0 = 1e18 m⁻³ and ~0.5 at Table I's 6e26. Below α ~ 1 the electron susceptibility vanishes
and there is no real ion feature — the IAW panel then shows the broad electron feature and
any bulk drift, not an acoustic doublet.

The IAW window is sized from C_s,ab (via `units.derive`), *not* from the ion thermal
spread, which the piston drift dominates — using the latter widens the window ~6× and
leaves the doublet spanning 2–3 pixels. Adjust with `--iaw-halfwidths`.

**`--velocity-scale-factor R`** divides all velocities by √R, putting the spectra on
real-ion-mass scales rather than the simulation's reduced-µ_p ones. `physical` derives it
from the config as (m_p/m_e)/mass_ratio — 18.36 at µ_p = 100, so v/4.285. Omitted by
default, leaving velocities as the simulation reports them. It matters more than a
relabelling: α ∝ 1/v_te, so on R1_paper it lifts α from 0.22–0.89 to **1.04–4.25**,
crossing into the collective regime where EPW satellites appear (RESULTS 2026-08-03).
Treat it as a *partial* physical correction — it does not undo the reduced-c temperature
offset (T_e,ab = 47 keV here against Table I's 470 eV). Outputs are suffixed `_scaled`
(or `--tag`) so the two variants never overwrite each other, and R is recorded in the npz.

**`--notch LO HI`** blanks a stray-light band in the EPW window, as a real diagnostic does
— without it the unshifted probe light and the ion feature swamp the far weaker EPW
satellites. Defaults to exactly the IAW window, so the two figures are complementary: what
the EPW panel blanks is what the IAW panel resolves. `--no-notch` disables it.

Useful options: `--probe-wavelength` (nm, default 532), `--angle` (deg, default 90),
`--position` (m, default domain centre), `--smoothing-window` / `--smoothing-iterations`
(defaults 31/2, well above the pipeline's own 9/1: above α ~ 1 the spectrum carries
|1 − χ_e/ε|², which diverges at the EPW resonance and turns VDF shot noise into speckle).

---

## `run_progress_logger.py`

A sidecar that watches a run's `run.log` and appends a checkpoint line to
`<run_dir>/progress.log` every N percent of `max_step`. It records **real
(wall-clock) elapsed time and ETA**, the WarpX compute rate, a **contention
factor** (wall-rate ÷ compute-rate — >1 when the machine is shared), and the
system load, so compute cost is trackable after the fact and runs can be paced
without babysitting. `--total` is auto-detected from `warpx_used_inputs` / the
deck (last `max_step` wins, matching ParmParse), so appended overrides are honored.

Launch it right after starting WarpX (it waits for `run.log` to appear):

```bash
OMP_NUM_THREADS=8 warpx.1d inputs > runs/R1_core/run.log 2>&1 &
python scripts/run_progress_logger.py runs/R1_core &          # every 10%, 30s poll
python scripts/run_progress_logger.py runs/R1_core --every-pct 5 --poll 20
```

Example `progress.log` line:

```
2026-07-23T17:02:31 | 25000/125000 (20.0%) | 0h13m | 0.0331 s/step(wall) | 0.0299 s/step(warpx) | ETA ~55m | x1.11 | load 18.4
```

The `x1.11` means wall time ran 11% over pure compute due to sharing the node —
a direct readout of the contention cost the CPU benchmark quantified.


---

## `make_run_readme.py`

One README per run directory. Every row carries its provenance, so "where did this number come
from?" is answerable without re-deriving anything:

* **`config.yaml:<key>`** — a primary. Editing it is the only way to change the run.
* **`derived: <formula>`** — computed by `units.derive`; never stored, never hand-copied.
* **Table I** — Schaeffer 2020's value, with the ratio and a `**OFF**` flag past 20%, plus the
  known cause when the repo already understands the deviation.

| Argument | Default | Meaning |
|---|---|---|
| `run_dir ...` | — | One or more run directories. |
| `--all` | off | Every `runs/*/` that has a `config.yaml`. |
| `--check` | off | Do not write; exit 1 if any README is stale. Use in CI / before a commit. |

**Prose is preserved.** Anything between `<!-- prose:begin -->` and `<!-- prose:end -->` is
hand-written and survives regeneration — that is where a run's story lives (why it exists, what
it showed, what to distrust). Everything outside those markers is rewritten from the config, so
never hand-edit the tables.

## `migrate_field_b0.py`

`field.B0_tesla` (tesla) is the primary background field; `v_A` is derived from it and the
ambient density. This script converts configs that still carry the old `field.vA_over_c`
primary, applying exactly the map the old code used
(`B0 = vA_over_c * c * sqrt(mu0 * namb * m_i)`), so the regenerated deck holds the *same*
numeric `B0` and `make_inputs.py --verify` still passes against an existing
`warpx_used_inputs`. Edits are line-based, so comments and key order survive.

`units.derive` refuses a config still carrying `vA_over_c` rather than guessing — the two
parameterizations disagree the moment the ambient density changes, which is the failure this
inversion exists to prevent (RESULTS 2026-07-31).

| Argument | Default | Meaning |
|---|---|---|
| `run_dir ...` | — | One or more run directories. |
| `--all` | off | Every `runs/*/` that has a `config.yaml`. |
| `-n, --dry-run` | off | Report the B0 each config would get; write nothing. |

## `table1.py`

Renders Schaeffer 2020 Table I in **three unit systems side by side** — PSC code units, the
physical plasma, and our WarpX SI deck — with the paper's own printed value in a fourth
column for checking. `--deck` prints the `config.yaml` values implied; `--show-work` prints
the Coulomb-logarithm algebra.

```bash
python scripts/table1.py                       # the standard set (below)
python scripts/table1.py --deck --show-work
python scripts/table1.py --Te-ab-eV 300 --mu 1836 --beta-ab 30
python scripts/table1.py --beta-0 0.02         # checks an over-determined value
```

**The input is a choice of real plasma.** A run of this problem is a set of dimensionless
numbers and corresponds to a whole *family* of real plasmas. You pick one member by choosing,
in physical units: `n_e,ab` (6e26 m⁻³), `T_e,ab` (470 eV), `lambda_ab` (20 d_e), `mu = m_i/m_e`
(100), `beta_ab` (1150), `T_0` (10 eV) and `n_e0/n_e,ab` (0.008). Everything else follows, and
**`beta_ab` sets B₀** → 7.026 T.

**`mu = 100` is a physical choice here, not just a code convenience** — the represented plasma
really is a light-ion plasma. That is what makes the set self-consistent, and it is why the ion
rows do **not** reproduce Table I's own SI column: that column is real hydrogen (µ = 1836), and
the paper's caption calls it "one possible set of experimentally-relevant physical values",
i.e. an illustration rather than a unit map. At µ = 100, C_s,ab = 909 km/s (not 210),
d_i,ab = 2.17 µm (not 9.31), 1/ω_ci0 = 80.9 ps (not 1.5 ns). Everything *dimensionless* still
matches: β_ab = 1150, β_0 = 0.196, M_A = 13.95, M_ms = 12.76, λ_ab = 20, d_i0 = 11.18 d_i,ab.
A bonus: at µ = 100 the reduced c is self-consistent at **c_sim/c_phys = 0.100** across both
the temperature and the velocity rows, whereas Table I's printed 0.02 only fits its velocity
rows — the 4.3× = √(1836/100) discrepancy is an artifact of its µ = 1836 physical column.

**`beta_0` is NOT free.** Given `n_e0`, `T_0` and B₀ (itself set by `beta_ab`), it is
determined: µ₀n_e0kT_0/B₀² = 0.196. Pass `--beta-0` to check a value; the script reports what
a different one would require.

**The one degree of freedom left is the speed of light.** §II sets
`c = sqrt(mu_p/T_e,ab) C_s,ab`, so the code's `theta_e,ab` *is* its speed of light. PSC picks
0.092 (O(0.1), non-relativistic), which at 470 eV means c_sim/c_phys = 0.100. WarpX has no
reduced-c option, so it must use the physical `theta_e = 9.2e-4`. Both then represent the same
plasma and agree on every dimensionless row **except the one that is the speed of light**,
c/C_s,ab: 33.0 for PSC against 329.7 for WarpX. The script prints a "reduced-c ledger"
isolating that, because it has two unavoidable costs:

1. **10× the timesteps** for the same 220 t_ab (dt is CFL-locked to dz/c while t_ab ∝ c):
   322,364 → 3,224,046.
2. **dz/λ_D,ab goes 0.99 → 9.89.** λ_D is c-independent but d_e ∝ c, so PSC's d_e,ab is only
   3.3 λ_D while at real c it is 33 λ_D — the same 0.3 d_e,ab cell is then 9.9 λ_D wide and
   will grid-heat. Resolving it needs dz ~10× smaller *and* dt with it, ~100× on top of the 10×.

**The payoff.** At T_e,ab = 470 eV the deck needs **no Coulomb-log dial**: `quantity:
lambda_ab, value: 20` resolves to lnΛ = 12.2 on its own, because ν_ei ∝ T^-3/2 and dropping
T_e,ab from 47 keV (θ_e = 0.092 at real c) to 470 eV raises ν_ei ~1000×. That is what turns
the deck's lnΛ = 1.22e5 into a physical number. Do *not* switch to `value: physical`, which
uses `24 − ln(√n/T)` = 6.2 — a different quantity that would not give λ_ab = 20.
