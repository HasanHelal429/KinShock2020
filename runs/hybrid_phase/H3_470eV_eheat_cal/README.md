# H3_470eV_eheat_cal — one calibration point for the electron-fluid heater amplitude

**Phase.** H3, `HYBRID_PLAN.md` §6.2 + §8
**Question.** `theta_e_equiv` was calibrated for the **ion** target; on the electron fluid the
same amplitude gives a 2945 eV steady-state slab `T_e` against the 470 eV being modelled.
**What `theta_e_equiv` lands the fluid on 470 eV?**
**Hypothesis.** `T_ss ∝ θ^p` with `p` between 1.0 and 1.5 — see below; the two ends are
physical limits, not fit slop.
**Expected.** **297 eV if `p` = 1.5, 638 eV if `p` = 1.0.** A factor 2.1, which is why this
point was placed where it was.
**Falsified by.** `T_ss` outside [250, 750] eV. That would mean the steady state is not a
power law in the deposition rate at all — most likely because the slab is not in balance by
`t*` = 5.6, which the parent's plateau (`dT_e/dt*` → 0 from `t*` ≈ 3.3) says it is.
**Kinetic reference.** **NONE — this run is UNVALIDATED and is not a physics run.** Its only
output is `T_ss(θ)`. A calibrated run cannot be cited for the quantity it was calibrated on
(`CLAUDE.md`), and here that is the whole run.

## Geometry

```
1D  |  shock normal z  |  lengths in d_i0 (ambient) = 24.26 um  |  solver: hybrid

      ###...............................................................
      ^                                                                ^
      symmetry                                                      open
      z = +0                                                z = +80

  #  piston   : 250 n_amb, +/-2 d_i0 about +0;  closure T_e = 396.9 eV
  .  ambient  : 1 n_amb, theta_i = 1.95695e-07;  closure T_e = 10 eV
     closure  : P_e = n0 T_e0 (n/n0)^1.66667, n0_ref = 1 n_amb  =>  T_e contrast 39.7x (adiabatic limit)
  B  field    : B0 = 7.026 T along x (perp. to z), 1/wci0 = 80.92 ps, v_A = 0.001 c
     grid     : 3200 cells, dz = 0.025 d_i0, dt = 0.05057 ps = 0.000625/wci0, 8960 steps = 5.6/wci0
```

Identical to `H3_470eV_eheat`. The box is **not** grown here even though the parent clamps at
`t*` ≈ 5.28 — this run is slower, not faster, because it deposits 10× less power, so the
clamp does not arise. Growing it would also break the one-variable comparison.

## Setup

Parent is **`H3_470eV_eheat`**. **One number changed:**

```yaml
operators.heater.amplitude.theta_e_equiv:  9.19767056e-4  ->  1.993066e-4
```

Rendered amplitude `H`: **6.197268695e24 → 6.251238548e23 m²/s³**, a ratio of 0.100871 —
which is `(θ/θ₀)^1.5` to six figures, confirming the mapping numerically rather than assuming
it.

Everything else byte-identical: `dz`, `dt`, 8960 steps, 38 substeps, ppc 3200/1600,
`random_seed`, box, injector, `target: electron_fluid`, `mode: advected`, and **collisions
off** — the parent's 2945 eV point is collisionless, so every point in this fit must be.

Requires the **`warpx-cda` fork, branch `feature/particle-heater`**.

### Why this is a run and not an edit

The source term in `dT_e/dt` is `m_e H` **exactly** — `S = (3/2) n_e m_e H` and
`dT_e/dt = (2/3)S/n_e`, so the density cancels — and `H = 8 θ^1.5 c³/(√M · width)` is exactly
`θ^1.5` (`deck.py:704`). So `T_ss = m_e H τ` and everything turns on the residence time `τ`:

| assumption on `τ` | law | `θ` for 470 eV |
|---|---|---|
| independent of `T` | `T_ss ∝ θ^1.5` | 2.706e-4 |
| `τ ∝ 1/√T` | `T_ss ∝ θ` | 1.468e-4 |

The second is the physical case if the heat leaves at the sound speed the heating itself sets
— and it is the **electron** pressure that drives the expansion carrying it out, so it is not
obviously wrong. The two are a factor **1.84** apart in `θ`. Both are defensible a priori; the
answer is an exponent, and an exponent is measured.

**θ = 1.993066e-4 is the geometric mean of the two predictions**, placed where they disagree
most usefully. Combined with the parent's `(9.19767056e-4, 2945 eV)` it pins `p` exactly, and
the production run at the fitted `θ` then **tests** the law rather than assuming it — the
pattern `H2_ionheat_cal2` used, where the refitted `M_A² = M_A0² + k·θ` predicted a held-out
point to 2.1 %.

### Reading it

The number wanted is the **plateau** slab `T_e`, not a final-frame value: the parent's approach
is strongly sub-linear (4002, 1799, 1155, 810, 474, 286, 162, 80, 4 eV per `t*`) and I once
called it a runaway by extrapolating its first two points. Take `T_ss` as the mean over
`t*` ∈ [4.0, 5.6] and quote `dT_e/dt*` over that window alongside it, so "it plateaued" is a
measurement rather than an impression.

## Cost

3200 cells × 8960 steps × 38 substeps, ppc 3200/1600 — identical to the parent, which measured
**13 min** at 8 threads. The heater's cost is a per-5-step pass over the slab and does not
depend on its amplitude. **Expect ~13 min**, `diags/` ~5 GB. Take wall time from a quiet box.

## Gates

`run_checks.py`: **6 pass, 0 warn, 0 fail, 2 info, 1 post-run** — all inherited unchanged from
the parent, as a one-number delta should give.

| Gate | Value | Pass? |
|---|---|---|
| GH1 `dt*wci0` / particle CFL | 6.25e-4 / 0.26 | pass |
| GH2 substeps vs requirement | 38 vs 3.49 (11× margin) | pass |
| GH3 `dz/d_i` (ambient / piston) | 0.025 / 0.395 | pass |
| GH4 `n_floor` clamp fraction | parent measured exactly 0 | post-run |
| GH5 ambipolar-front `dz` convergence | inherited from the dz ladder | pass |
| GH6 energy closure incl. the fluid | `ElectronFluidEnergy` reduced diag | post-run |
| GH7 agreement with the kinetic reference | **N/A — calibration run** | info |
| GH9 no-drive control | `H3_470eV_advected` | pass |

The **advective CFL** is enforced by WarpX rather than `run_checks.py` and is the gate most
likely to bite when heating raises flow speeds — but this run heats *less* than the parent,
which cleared it, so it has more margin, not less.

## Result

<pending — run not yet launched>

## Retracted

Nothing yet.
