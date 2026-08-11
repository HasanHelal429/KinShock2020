# H_phase — where does the upstream heat come from?

> ## ⛔ SUPERSEDED 2026-08-11 — do not run the remaining 7 points
>
> The control point `hs_dz1_ppc100` ran, and **it heats by 0.00 eV** (10.00 → 9.98 eV over
> 30 t_ab) where the full deck `h0_baseline` heats 10.0 → 40.7 eV over the same window at
> the same dz/λ_D and N_D. The box is not inert — it carries the production upstream's
> E_z noise exactly (29–30 v_A B₀ against 30.2) — so the noise and the heating are
> **different problems**, and the quiescent discretization does not heat this plasma.
>
> `h0`'s own data agrees: its far upstream is flat at 9.98 eV until t = 2.50 t_ab and
> rises only after the **light-crossing time of 2.59 t_ab**, with no piston particle ever
> arriving. `R1_paper` shows the identical behaviour at *its* light-crossing time of
> 25.9 t_ab. The upstream is heated by an EM precursor whose arrival is set by **c/v_ti**,
> which is 10× too large here — and **no dz and no ppc changes c/v_ti**.
>
> Full argument and numbers: RESULTS.md 2026-08-11 (later). The design below is kept
> because the *box* is sound and correctly measured a real null; it is the premise that
> was wrong. Anything reusing it should ask a question about a quiescent plasma.


**The question.** `R1_paper_470eV` warms its far upstream by ~+30 eV over 30 t_ab
(`h0_baseline`, RESULTS 2026-08-05) from a nominal T₀ = 10 eV. The production run is
190 t_ab. If that rate persists the upstream is no longer the 10 eV plasma the Table I
regime is defined against, and M_ms follows it down. **Which discretization knob buys the
heating back, and at what price?**

**Why a separate box rather than more full runs.** Three full-deck attempts have now
failed to move it — `h1_filter8` (12.5%), `h2_shape3` (1.7% on the E_z noise),
`i0/i1` (θ-implicit and Villasenor, indistinguishable). Each cost hours and each varied
*one* knob at the production resolution. The two knobs that were never varied are the two
that set the physics of the finite-grid instability, and they are the expensive ones: `dz`
and `ppc`. In the full deck a 4× refinement is a 5-day run. In a box with no piston it is
minutes.

## The box

A uniform, periodic, magnetized ambient plasma and **nothing else** — no piston, no
heater, no injector, no shock:

| | |
|---|---|
| n | 4.8×10²⁴ m⁻³ (= 0.008 n_e,ab, as production) |
| T_e = T_i | 10 eV |
| B₀ | 7.0264 T along x, perpendicular |
| box | 230.4 d_e,ab = **2.06 d_i0** = 4.66 thermal ion gyroradii = 4659 λ_D, fixed in metres across the sweep |
| boundaries | periodic, both faces |
| collisions | the production block, unchanged (lnΛ = 12.23) |
| window | 440 000 steps at the coarsest dz = **30.0 t_ab**, the same window `h0_baseline` measured |

Everything else — n₀, mass_ratio, B₀, θ₀, cfl, particle_shape, the collision target — is
byte-identical to `R1_paper_470eV`. The resolution decision has to be made *with* the
operator set the production run actually carries, not against a cleaner surrogate.

**The measurement is free.** With no piston and no injector, the domain-wide mean kinetic
energy per particle *is* the temperature, so `EP`/`PN` (the reduced diagnostics, ~1000
rows per run) carry the whole result. This is exactly what `grid_heating.py`'s docstring
says the reduced diagnostic *cannot* do in the full deck, where EP is dominated by real
piston heating and the measurement needs a spatial window far upstream. Here there is no
upstream and no downstream.

**The box keeps the piston block in `config.yaml` anyway** — `theta_e_heat` defines
C_s,ab and hence t_ab, the time unit every normalization in this repo is drawn against.
No species carries `role: piston`, so none of those numbers reach a particle, and
`operators:` is absent so no heater/injector is emitted (`kinshock.deck.render`,
2026-08-11). That is not merely tidier: `particle_heater`'s rate is H ∼ 1/foil_width
(`ParticleHeater.cpp:207`), and `slab_halfwidth_di: 0` — which is what makes the ambient
profile uniform — would be a division by zero, not a no-op.

## The grid, and why these eight points

Two *different* things go wrong when a PIC plasma is under-resolved, and the production
run is bad at both. They are separated by two parameters, not one:

- **dz/λ_D** — aliasing. Modes above the grid Nyquist fold back and drive the
  finite-grid instability. Production: **6.07**. `R1_paper`, which behaves: **0.60**.
- **N_D = ppc·λ_D/dz** — the number of particles per Debye length, which sets the
  *amplitude* of the thermal-fluctuation noise (∼1/√N_D). Production: **16.5**.
  `R1_paper`: **167**.

Refining dz improves *both*; raising ppc improves *only* N_D. So a one-dimensional dz
scan cannot tell you which one matters — and the answer decides the price, because at
fixed N_D **ppc is 4× cheaper than dz** (cost ∝ ppc, but ∝ dz⁻² once the timestep follows
the CFL). The grid below is built so that lines of constant N_D cross lines of constant
dz/λ_D:

| run | dz/λ_D | ppc | N_D | cells | steps | cost |
|---|---|---|---|---|---|---|
| `hs_dz1_ppc100` | 6.07 | 100 | **16.5** | 768 | 440 000 | 1 |
| `hs_dz1_ppc400` | 6.07 | 400 | 66.0 | 768 | 440 000 | 4 |
| `hs_dz2_ppc50` | 3.03 | 50 | **16.5** | 1 536 | 880 000 | 2 |
| `hs_dz2_ppc100` | 3.03 | 100 | 33.0 | 1 536 | 880 000 | 4 |
| `hs_dz4_ppc25` | 1.52 | 25 | **16.5** | 3 072 | 1 760 000 | 4 |
| `hs_dz4_ppc100` | 1.52 | 100 | 66.0 | 3 072 | 1 760 000 | 16 |
| `hs_dz8_ppc25` | 0.76 | 25 | 33.0 | 6 144 | 3 520 000 | 16 |
| `hs_dz8_ppc50` | 0.76 | 50 | 66.0 | 6 144 | 3 520 000 | 32 |

`hs_dz1_ppc100` **is the production parameter point** — same dz/λ_D, same N_D. It is the
control, and it must reproduce `h0_baseline`'s +29.7 eV/30 t_ab (or explain why a
piston-free box does not).

The three readings the grid supports:

- **Aliasing, at fixed noise.** `dz1_ppc100 → dz2_ppc50 → dz4_ppc25`, all at N_D = 16.5,
  dz/λ_D falling 6.07 → 3.03 → 1.52. If the heating collapses along this line the cause is
  the finite-grid instability and only dz will fix it.
- **Noise, at fixed aliasing.** `dz1_ppc100 → dz1_ppc400`, dz/λ_D pinned at 6.07,
  N_D 16.5 → 66. If the heating collapses *here*, ppc is the lever and the production run
  costs 4× rather than 16×.
- **Aliasing again at high N_D.** `dz1_ppc400 → dz4_ppc100 → dz8_ppc50`, all N_D = 66.
  Confirms whichever of the two above is real, and `dz8_ppc50` (dz/λ_D = 0.76, N_D = 66)
  is the converged anchor: it is inside `R1_paper`'s regime on both axes, so if *it*
  still heats, the cause is neither.

`hs_dz8_ppc25` (N_D = 33 at dz/λ_D = 0.76) pairs with `hs_dz2_ppc100` (N_D = 33 at 3.03)
as the second matched-N_D pair, so the aliasing reading has two independent lines.

## Running and reading

```bash
python scripts/make_inputs.py runs/H_phase/<ID>          # deck from config
scripts/launch.sh -b -g 0 runs/H_phase/<ID>              # one GPU, one box
python scripts/make_inputs.py runs/H_phase/<ID> --verify  # after the run
python scripts/heating_rate.py runs/H_phase/*            # the sweep table
```

`heating_rate.py` reads `EP`/`PN` and reports T_e(t), T_i(t) in eV, the fitted dT/dt in
eV per t_ab, and the projection to the 190 t_ab production window. Quote the **rate
ratios** against `hs_dz1_ppc100`, not the absolute projections — the ratios cancel the
extrapolation risk, which is the same trap the `h0_baseline` fitted asymptote fell into
(RESULTS 2026-08-05: 132 → 87 eV extrapolated 2–3× past its data and is not quotable).

## Measured cost

Every point benchmarked 2026-08-11 with its own deck, 1200 steps, `s/step` taken from
WarpX's cumulative `Evolve time` between step 200 and 1200 so init and the first sort drop
out. One RTX 4070, `amrex.the_arena_init_size` capped. The harness was validated against a
known number: `R1_paper_470eV` measured 0.01533 s/step here against 0.01415 recorded
2026-08-04 (+8 %, and this deck carries the 2-box decomposition that measurement noted
costs 3.3 %).

| run | s/step | **1 GPU** | CPU, 8 thr¹ | output |
|---|---|---|---|---|
| `hs_dz1_ppc100` | 0.00102 | **7 m** | 2 h 23 m | 0.16 GB |
| `hs_dz2_ppc50` | 0.00097 | **14 m** | 4 h 44 m | 0.18 GB |
| `hs_dz1_ppc400` | 0.00232 | **17 m** | 9 h 50 m | 0.60 GB |
| `hs_dz2_ppc100` | 0.00135 | **19 m** | 9 h 40 m | 0.32 GB |
| `hs_dz4_ppc25` | 0.00094 | **27 m** | 9 h 32 m | 0.21 GB |
| `hs_dz4_ppc100` | 0.00207 | **1 h 00 m** | 38 h | 0.65 GB |
| `hs_dz8_ppc25` | 0.00131 | **1 h 16 m** | 37 h | 0.41 GB |
| `hs_dz8_ppc50` | 0.00204 | **1 h 59 m** | 76 h | 0.71 GB |
| **all eight** | | **5 h 43 m** serial · **2 h 51 m** on both cards | | 3.2 GB |

The three cheapest points (`dz1_ppc100`, `dz2_ppc50`, `dz4_ppc25`) are the entire
constant-N_D aliasing line and cost **48 minutes together** — run those first.

Peak device memory is 331 MB for the 6×10⁶-particle production deck (55 B/particle), so
nothing here is close to memory-bound. The arena still has to be capped when a card is
shared: AMReX allocates 3/4 of *total* device memory at init regardless of need.

¹ **The CPU column is not the best CPU number.** These configs set
`numerics.max_grid_size = n_cell` — one box, which is right for a single GPU and is what
`launch.sh -g` requires. On CPU that is a trap: AMReX does not tile in 1D, so one box is
one tile and only one OpenMP thread gets work. Measured directly on `hs_dz1_ppc400`,
8 threads: **0.1073 s/step with one box against 0.0148 with `max_grid_size = 64`, a 7.3×
penalty.** To run these on CPU, delete the key from `config.yaml` and regenerate the deck
— do not override it on the command line, which `--verify` would then flag. Even then the
GPU is ~6× faster on this deck, so use it unless both cards are taken.

⚠ **The same trap is live in `R1_paper_470eV` itself.** That config gained
`max_grid_size: 15000` on 2026-08-04 *for GPU*, which is two boxes — and two boxes on 8
threads measures **0.3924 s/step against the 0.11169 in `runs/opt_phase/SUMMARY.md`, 3.5×**
(2026-08-11, this harness). The opt_phase CPU rows were measured before that key existed,
so **every CPU projection in that table is stale for the config as it now stands**: the
5.33 d at 8 threads is really ~18 d. The GPU rows are unaffected.
