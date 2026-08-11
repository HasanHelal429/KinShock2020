# S_phase — why does the 470 eV run make a shock too early?

**The question.** At t·ω_ci0 = 0.39 `R1_paper_470eV` already shows ambient ions
accelerated into a shock-like structure. `R1_paper` — the *same dimensionless problem*,
identical n_amb, dz, dt, ρ_i0/dz, β₀ and M_ms — has no shock structure at all that early
(RESULTS 2026-08-05, `scripts/compare_phase.py`). Something in the 470 eV run is
accelerating ambient ions before any sweeping-up could have happened. **Is it a
resolution artifact, and if so which one?**

**What is already established, and what it does not settle.** Normalized by the run's own
units, the upstream E_z is 13× larger in the 470 eV run at t·ω_ci0 = 0.02 — *before any
shock exists* (30.2 vs 2.25 in v_A·B₀). That is not a mystery on its own: the exact
identity for the thermal-fluctuation field,

&nbsp;&nbsp;&nbsp;&nbsp;E_th = √(nT/ε₀) ⟹ **E_th/(v_A B₀) = β · c/v_ti**

puts both runs at the same 5–7 % of their own thermal noise scale, and β is identical by
construction, so the whole 13× is c/v_ti — the reduced-c trick being undone (RESULTS
2026-08-05). What that argument does *not* say is whether an upstream E_z at 5 % of the
thermal scale is enough to *start a shock early*, or whether the early acceleration has a
different cause and the noise is incidental. Only a resolution scan answers that, and it
has to be a scan of the actual shock, not of a quiet box — which is what separates this
phase from `H_phase`.

## The reduced problem

`R1_paper_470eV`'s physics **exactly** — same piston, heater, injector, collisions, B₀,
densities, temperatures, cfl, particle_shape — with two truncations:

| | production | here |
|---|---|---|
| domain | 9000 d_e,ab = 80.5 d_i0 | 1344 d_e,ab = **12.02 d_i0** |
| window | 2 784 400 steps = 5.60/ω_ci0 | 150 000 steps = **0.3018/ω_ci0** |
| cost | 7.9 h on 2 GPUs | ~1/124 of that per baseline point |

By t·ω_ci0 = 0.30 the piston has travelled 3.12 d_i0 and the model shock 4.19 d_i0, so
the domain gives 2.9× head-room ahead of the front. The window brackets the anomaly:
whatever starts at or before 0.39 starts inside it.

### The one thing to check before trusting the truncated domain

`field_hi = open` is **pec**, which *reflects* fields — only particles are absorbed
(`kinshock.deck._BC_MAP`; Silver-Mueller is not available because the B-field divergence
cleaner runs whenever a background B is set). Light crosses 12.0 d_i0 in ~6000 steps, so
the piston turn-on precursor makes ~25 round trips inside this window against ~4 in the
full domain. That is a real difference.

**Which is why `ss_dz1_ppc100` is a control, not just the cheapest point.** It is the
production parameter point (dz/λ_D = 6.07, N_D = 16.5) on the truncated domain, and it
must reproduce the full `R1_paper_470eV` over t·ω_ci0 ≤ 0.30 — upstream E_z level,
|B_perp| growth, and the ambient phase space — before any other point here means
anything. The full run's `diag_fields` cadence (2500 steps) gives ~60 frames in that
window, so **the control costs no new compute**; it is a comparison against data already
on disk. If it fails, the domain has to grow until it passes.

## The grid, and why these six points

Same two-parameter logic as `H_phase` (read that README's grid section — the argument is
identical and is not repeated here): **dz/λ_D** controls aliasing, **N_D = ppc·λ_D/dz**
controls the thermal-noise amplitude, refining dz improves both, raising ppc improves only
N_D, and at fixed N_D ppc is 4× cheaper.

| run | dz/λ_D | ppc | N_D | cells | steps | cost |
|---|---|---|---|---|---|---|
| `ss_dz1_ppc100` | 6.07 | 100 | **16.5** | 4 480 | 150 000 | 1 |
| `ss_dz1_ppc400` | 6.07 | 400 | 66.0 | 4 480 | 150 000 | 4 |
| `ss_dz2_ppc50` | 3.03 | 50 | **16.5** | 8 960 | 300 000 | 2 |
| `ss_dz2_ppc100` | 3.03 | 100 | 33.0 | 8 960 | 300 000 | 4 |
| `ss_dz4_ppc25` | 1.52 | 25 | **16.5** | 17 920 | 600 000 | 4 |
| `ss_dz4_ppc100` | 1.52 | 100 | 66.0 | 17 920 | 600 000 | 16 |

- **Aliasing at fixed noise:** `dz1_ppc100 → dz2_ppc50 → dz4_ppc25` (N_D = 16.5 throughout).
- **Noise at fixed aliasing:** `dz1_ppc100 → dz1_ppc400` (dz/λ_D = 6.07 throughout).
- **N_D = 66 at two resolutions:** `dz1_ppc400` vs `dz4_ppc100`.

No dz/8 point is staged. It would be 64 units — more than everything above combined — and
it is only worth buying if dz/4 has not converged.

## What to measure

The diagnostics are set up for the *field* question, not the late-time shock: 150
field-only frames (`diag_fields`, every 1000 steps at the coarsest dz) and 30 full frames
with particles.

1. **Upstream E_z**, `scripts/plot_ez.py --tmax 0.3` — the level far ahead of the piston,
   in v_A·B₀ and against the β·c/v_ti thermal prediction. Does it fall as 1/√N_D (noise),
   faster (aliasing), or not at all?
2. **Coherent |B_perp|/B₀**, `scripts/plot_bperp_pileup.py` — a real shock *builds* its
   barrier from B₀ upward. ⚠ this script still takes a domain-wide max and needs boundary
   exclusion added before it is used here; at 12 d_i0 the wall spike is a much larger
   fraction of the domain than it was at 80.5.
3. **Ambient phase space**, `scripts/compare_phase.py --at 0.11 0.21 0.30` — the direct
   read: how far ahead of the piston does the accelerated population reach, and does it
   retreat as N_D or dz/λ_D improves?

Read 1 and 3 together. The result that would settle it is the one where the accelerated
ambient population disappears along one of the two lines above and not the other.

## Running

```bash
python scripts/make_inputs.py runs/S_phase/<ID>
scripts/launch.sh -b -g 0 runs/S_phase/<ID>              # one GPU, one box
python scripts/make_inputs.py runs/S_phase/<ID> --verify
```

## Measured cost

Benchmarked 2026-08-11, same harness as `H_phase` (1200 steps, `s/step` from WarpX's
cumulative `Evolve time` between steps 200 and 1200; one RTX 4070; arena capped). Wall
times carry a ×1.10 allowance for the injector growing the piston population over the
window — the production run grows 1.279× over 220 t_ab and this window is 10.2 t_ab.

| run | s/step | **1 GPU** | CPU, 8 thr¹ | output |
|---|---|---|---|---|
| `ss_dz1_ppc100` | 0.00315 | **8 m** | 5 h 15 m | 1.4 GB |
| `ss_dz2_ppc50` | 0.00303 | **16 m** | 10 h 14 m | 1.4 GB |
| `ss_dz1_ppc400` | 0.01013 | **27 m** | 22 h | 5.2 GB |
| `ss_dz2_ppc100` | 0.00517 | **28 m** | 21 h | 2.7 GB |
| `ss_dz4_ppc25` | 0.00300 | **33 m** | 27 h | 1.6 GB |
| `ss_dz4_ppc100` | 0.00932 | **1 h 42 m** | 84 h | 5.4 GB |
| **all six** | | **3 h 37 m** serial · **1 h 48 m** on both cards | | 17.7 GB |

`ss_dz1_ppc100` — the control that has to pass before the rest mean anything — is **8
minutes**. Run it and compare against the full run's first 60 field frames before
launching anything else here.

¹ **The CPU column is not the best CPU number** — see the same footnote in
`runs/H_phase/README.md`. `max_grid_size = n_cell` starves OpenMP by a measured 7.3× and
must be deleted from `config.yaml` for a CPU run.

## What the answer costs: the full production run

The decision this phase and `H_phase` exist to make is *which resolution the production
`R1_paper_470eV` rerun should use*. Scaled from the measured anchor (0.01533 s/step on 1
GPU for the production deck), over the same 9000 d_e,ab domain and 2 784 400-step window,
with the ×1.279 particle-growth factor `opt_phase` uses. Cost goes as **k²·(ppc/100)** —
refining dz by k buys k× the cells *and* k× the steps, because dt is CFL-locked to dz.

| dz | dz/λ_D | ppc | N_D | cells | steps | 1 GPU | **2 GPUs** |
|---|---|---|---|---|---|---|---|
| dz/1 | 6.07 | 100 | 16.5 | 30 000 | 2.78 M | 15.2 h | **8.6 h** ← what already ran |
| dz/2 | 3.03 | 50 | 16.5 | 60 000 | 5.57 M | 30.3 h | **17.1 h** |
| dz/1 | 6.07 | 400 | 65.9 | 30 000 | 2.78 M | 60.7 h | **34.3 h** |
| dz/2 | 3.03 | 100 | 33.0 | 60 000 | 5.57 M | 60.7 h | **34.3 h** |
| dz/4 | 1.52 | 25 | 16.5 | 120 000 | 11.14 M | 60.7 h | **34.3 h** |
| dz/4 | 1.52 | 100 | 65.9 | 120 000 | 11.14 M | 243 h | **5.7 d** |
| dz/8 | 0.76 | 50 | 65.9 | 240 000 | 22.28 M | 485 h | **11.4 d** |
| dz/8 | 0.76 | 100 | 131.9 | 240 000 | 22.28 M | 971 h | **22.9 d** |

**Read the three 34.3 h rows together — that is the whole point of the sweep.** The same
budget buys either N_D 16.5 → 66 with the aliasing untouched (ppc 400), or N_D → 33 with
aliasing 2× better (dz/2, ppc 100), or aliasing 4× better with the noise untouched (dz/4,
ppc 25). They are not interchangeable, and nothing measured so far says which one is the
right purchase. ~9 h of sweep decides how ~34 h (or 5.7 d) of production is spent.

Device memory is never the constraint: 331 MB for the production deck's 6×10⁶ particles
(55 B/particle), so even dz/8 at ppc 100 is 3.4 GB. Cap `amrex.the_arena_init_size`
anyway when sharing a card.
