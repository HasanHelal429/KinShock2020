# `scripts/` — analysis & bring-up drivers

Command-line entry points for the KinShock2020 replication of Schaeffer 2020.
Each script is config-driven: it reads `runs/<RUN_ID>/config.yaml`, derives all
physical scales through `kinshock.units`, and writes outputs under
`media/<run_id>/` (or `media/testing/` for bring-up checks).

All scripts are run from the repository root:

```bash
python scripts/<script>.py [run_dir] [options]
```

The `run_dir` positional argument defaults to `runs/R1` for every script except
`make_config.py`, where it is required.

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `make_config.py` | Regenerate/verify `config.yaml` from the deck WarpX actually ran | `warpx_used_inputs` | `config.generated.yaml` or `config.yaml` |
| `run_checks.py` | Bring-up / progress checks (works before any sim output exists) | `config.yaml`, plotfiles, reduced diags | `media/testing/*.png` |
| `make_figures.py` | Reproduce the paper's shock diagnostics (analyses A–D) | `config.yaml`, plotfiles | `media/<run_id>/*.png`, `criteria.json` |
| `make_movies.py` | Animated density + phase-space movies | `config.yaml`, plotfiles | `media/<run_id>/*.mp4` |

Typical workflow for a finished run:

```bash
python scripts/make_config.py runs/R1 --write   # lock config to what was simulated
python scripts/run_checks.py    runs/R1         # sanity: scales vs Table I, conservation
python scripts/make_figures.py  runs/R1         # A–D diagnostics + criteria table
python scripts/make_movies.py   runs/R1         # optional animations
```

---

## `make_config.py`

Regenerate / verify a run's `config.yaml` from the WarpX inputs it actually ran.

WarpX writes `warpx_used_inputs` (the fully-resolved input deck) into each run
directory. This script parses it, resolves the `my_constants` expressions
numerically, maps them to the KinShock2020 config primaries, and either writes
`config.generated.yaml` or overwrites `config.yaml`. It always diffs against the
existing `config.yaml` so the config provably matches what was simulated
(REPLICATION_PLAN.md §6.0a).

**Positional argument**

| Argument | Required | Description |
|---|---|---|
| `run_dir` | yes | Run directory containing `warpx_used_inputs` and `config.yaml`. |

**Options**

| Flag | Type | Default | Description |
|---|---|---|---|
| `--inputs <path>` | path | `<run_dir>/warpx_used_inputs` | Parse this inputs file instead of the run's `warpx_used_inputs`. |
| `--write` | flag | off | Overwrite `config.yaml` in place. If omitted, writes `config.generated.yaml` and leaves `config.yaml` untouched. |

**Examples**

```bash
python scripts/make_config.py runs/R1                # verify + write config.generated.yaml
python scripts/make_config.py runs/R1 --write        # overwrite config.yaml
python scripts/make_config.py runs/R1 --inputs deck  # parse an arbitrary deck instead
```

Verification prints `OK` when the config primaries (n0, mass ratio, `vA_over_c`,
`max_step`, grid resolution, heater θ, densities, …) match the simulated deck
within a `1e-3` relative tolerance, or lists each mismatch.

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
  density / B_perp line-outs → `shock_lineouts.png`
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
