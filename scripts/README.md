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
| `-L, --logger` | off | Also start `run_progress_logger.py` in the background (stdout → `<run_dir>/logger.out`, gitignored). |
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
