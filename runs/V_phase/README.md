# V_phase — is the notch caused by an unphysically fast shock?

**The hypothesis.** At µ = 100 the ablation sound speed is √(1836/100) = **4.2849× too high**,
so the shock runs at **4182 km/s** against Table I's physical **980**. Every dimensionless
number is still right — M_A 14, M_ms 13, β₀ 0.2 — because absolute speed is the intended
casualty of a reduced mass ratio. But if the notch depends on the *absolute* speed, or on a
ratio Table I's set does not pin, the µ = 100 runs are in the wrong place.

**Why the ablation temperature cannot simply fix it.** v_sh = 4.6·√(T_ab/m_i) and
v_te,ab = √(T_ab/m_e), so

> **v_sh / v_te,ab = 4.6/√µ — exactly, with T_ab cancelling.**

Verified across 25.6 eV, 470 eV and 47 keV: 0.4600 at µ = 100 and 0.1074 at µ = 1836,
independent of T_ab. The shock is 4.28× too fast *relative to the ablation electrons* at any
ablation temperature. Only µ changes that. What T_ab **can** restore is the absolute speed:
T_ab = 470 × 100/1836 = **25.599 eV** gives v_sh = 976 km/s at µ = 100.

## The trilemma

At µ = 100 you cannot have physical v_sh **and** M_A = 14 **and** β₀ = 0.2 **and** the paper's
T₀ = 10 eV. This phase runs all three corners, so which constraint is load-bearing is
measured rather than argued.

| run | change | v_sh | M_A | M_ms | β₀ | ε₀ | cells | steps | cost |
|---|---|---|---|---|---|---|---|---|---|
| reference `R1_470eV_s3f8_t2` | — | 4182 | 13.95 | 12.76 | 0.196 | 0.00442 | 101,984 | 7.95 M | 29.8 GPU-h |
| **`vA_lowTab`** | T_ab → 25.6 eV | **976** | **3.26** | 2.98 | 0.196 | 0.00442 | 59,624 | 10.16 M | **~22 GPU-h** |
| **`vB_lowB0`** | + B₀ → 1.640 T | **976** | 13.95 | 6.51 | **3.59** | 0.00442 | 59,624 | 10.16 M | **~22 GPU-h** |
| **`vC_scaled`** | + T₀ → 0.545 eV | **976** | 13.95 | 12.76 | 0.196 | **0.00103** | 255,496 | 43.5 M | **~409 GPU-h** |

All three hit **976.1 km/s** against Table I's 980 (0.4%).

**vA and vB cost the same, and not by accident:** vA needs 4.28× more t\* to reach the same
piston travel, while vB needs 4.28× more steps per t\* because ω_ci0 ∝ B₀ falls. The two
factors cancel exactly.

**vC is 18.4× dearer** because the lower B₀ costs 4.28× in clock *and* the lower T₀ costs
4.28× in grid (λ_D ∝ √T₀).

## What each isolates

- **`vA_lowTab`** — only T_ab moves. Physical v_sh, but M_A collapses to 3.26, barely above
  the ~2.7 first critical Mach number for a perpendicular shock. This is a **marginally
  supercritical** shock and may reflect few ions, so a clean result is ambiguous between
  "physical speed fixed it" and "M_A 3.3 is too weak to do anything". vB breaks that tie.
- **`vB_lowB0`** — the paper's shock strength at the paper's speed. The cost is β₀ = 3.59,
  and **β₀ was already excluded as a notch discriminator by crossing pair** (RESULTS
  2026-09-04), which makes vB the cleanest of the three: the one quantity it breaks is one
  already shown not to matter. ⚠ M_ms does drop to 6.51, which is *not* excluded.
- **`vC_scaled`** — every Table I number restored, but only via a 0.545 eV upstream, which
  drops ε₀ 4.28× **below** the notching value. The ε₀ discriminator therefore predicts vC
  notches *harder*: it confirms rather than discriminates. Run it only if vA and vB come back
  clean and the speed hypothesis needs its strongest test.

## What the existing discriminator predicts

ε₀ separated notching from clean across six runs at two mass ratios (RESULTS 2026-09-05).
vA and vB hold ε₀ at 0.00442, the **notching** value, so **ε₀ predicts they still notch
despite a physical v_sh**. If they come back clean, ε₀ is not the cause and the absolute
speed is. That is the whole point of running them.

## Held fixed, and one deliberate departure

Grid, ppc, filter and shape are `R1_470eV_s3f8_t2`'s unchanged — dz/λ_D = 0.758, ppc 100,
filter_npass 8, particle_shape 3 — so nothing turns on resolution or N_D. Box 20 d_i0. The
clock is set by **piston travel of 6.2 d_i0**, where the reference notch has saturated, *not*
by t\*: t\* is not comparable across these because ω_ci0 changes with B₀.

**Collisions are OFF**, the single departure from the reference. At 25.6 eV neither
collisional option was defensible: Table I's λ_ab = 20 inverts to lnΛ = 0.036, far below the
Coulomb logarithm's domain of validity, while the physical NRL value there (3.321, against
6.231 at 470 eV) makes ν_ei ∝ n T^−3/2 lnΛ about **42× stronger** than the reference — a
large confound in a test about speed. Off removes the nuisance parameter and is consistent
across all three. Justified by measurement: R1_coll reproduces R1_warm on six of seven
criteria. If V_phase shows a difference, a collisions-on control is the follow-up.

## Running them

```bash
perlmutter/submit.sh vphase_ab  --time 30:00:00   # vA + vB, ~45 GPU-h total
perlmutter/submit.sh vphase_c   --time 48:00:00   # vC alone, ~409 GPU-h -- NOT recommended yet
```
