"""Derived scales for a KinShock2020 run.

Everything the analysis needs beyond the raw primaries in ``config.yaml`` is
computed HERE from those primaries, so there is exactly one source of truth and
no risk of a script drifting out of sync with the deck (the design smell in the
reference script ``warpx-cda/heating_operator/scripts/make_shock_figures.py``,
which hard-codes deck constants at the top).

Reproduction convention (see REPLICATION_PLAN.md §1): the dimensionless Schaeffer
2020 run is reproduced with the REAL speed of light, so every "x c" in Table I is
a fraction of the real c. The absolute reference density ``n0`` is a free choice;
only ratios matter for the shock physics.

Pure Python (``math`` only) so it imports and runs without yt / numpy — enabling
the config↔Table-I unit test in the plan's checklist.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict

# --- physical constants (SI, CODATA-ish; match make_shock_figures.py) ---
ME = 9.1093837015e-31      # electron mass [kg]
QE = 1.602176634e-19       # elementary charge [C]
C = 299792458.0            # speed of light [m/s]
EPS0 = 8.8541878128e-12    # vacuum permittivity [F/m]
MU0 = 1.25663706212e-6     # vacuum permeability [H/m]
ME_C2_J = ME * C * C       # electron rest energy [J]
ME_C2_EV = ME_C2_J / QE    # electron rest energy [eV] (~511 keV)

# NRL Plasma Formulary (Braginskii) electron-ion collision rate, the standard
# definition behind both "mean free path" and the paper's collisionality:
#     nu_ei = 2.91e-6 * n_e[cm^-3] * lnLambda * T_e[eV]^{-3/2}   [s^-1]
NU_EI_NRL = 2.91e-6


@dataclass
class Scales:
    """Derived physical + normalized scales for one run.

    All quantities are SI unless the name carries an explicit unit suffix.
    Speeds are also reported as fractions of c; temperatures as theta = kT/(m_e c^2).
    """

    # reference / masses
    n0: float                 # ablation reference density n_e,ab [m^-3]
    mass_ratio: float
    mi: float                 # ion mass [kg]

    # ablation (piston) scales
    wpe: float                # omega_pe at n0 [rad/s]  -> sets d_e, t normalization
    de: float                 # electron skin depth c/omega_pe [m]
    di: float                 # ion inertial length d_e*sqrt(mass_ratio) [m]
    t_ab: float               # ablation time d_i / C_s,ab [s]
    Cs_ab: float              # ablation ion-acoustic speed sqrt(theta_e/mass_ratio) c [m/s]
    nt: float                 # piston/target density [m^-3]

    # upstream (ambient) scales
    namb: float               # ambient density n_e0 [m^-3]
    de0: float                # ambient electron skin depth [m]
    di0: float                # ambient ion inertial length [m]
    Cs0: float                # upstream sound speed [m/s]
    B0: float                 # background field [T]
    vA: float                 # Alfven speed [m/s]
    wci0: float               # ion cyclotron frequency eB0/mi [rad/s]
    wci0_inv: float           # 1/wci0 [s]
    rho_i0: float             # vp-gyroradius vp/wci0 [m] (z* normalization)

    # model speeds (Table I; measured values override post-run)
    vp_model: float           # piston speed [m/s]
    vsh_model: float          # shock speed [m/s]

    # dimensionless
    MA: float                 # Alfvenic Mach vsh/vA
    Mms: float                # magnetosonic Mach vsh/sqrt(vA^2+Cs0^2)
    beta_ab: float            # ablation electron beta
    beta_0: float             # upstream electron beta

    # numerics
    dz: float                 # cell size [m]
    dt: float                 # time step [s]
    dt_wpe: float             # dt * omega_pe (dimensionless)
    n_cell: int
    domain_halfwidth: float   # [m]
    steps_per_wci0: float     # timesteps per upstream ion gyroperiod (1/wci0)

    # --- collisions (ablation-unit electron-ion; None when the config has no
    #     `collisions` block, i.e. the run is collisionless) ---
    Te_ab_eV: float = 0.0        # ablation electron temperature [eV]
    vte_ab: float = 0.0          # ablation electron thermal speed sqrt(kT_e/m_e) [m/s]
    wce_ab: float = 0.0          # electron cyclotron frequency eB0/m_e [rad/s]
    coulomb_log_nrl: float = 0.0  # the PHYSICAL lnLambda at (n_e,ab, T_e,ab)
    coulomb_log: float | None = None   # lnLambda actually used in the deck
    nu_ei_ab: float | None = None      # electron-ion collision rate [s^-1]
    mfp_ei_ab: float | None = None     # electron mean free path v_te/nu_ei [m]
    lambda_ab: float | None = None     # paper's collisionality wce_ab/nu_ei_ab
    nu_ei_dt: float | None = None      # nu_ei * dt (must be << 1 for Takizuka-Abe)

    # normalization helpers below use these; keep a copy of a few config bits
    _meta: dict = field(default_factory=dict)

    # ---- normalization helpers (SI -> paper units) ----
    def z_over_de(self, z):      return _div(z, self.de)
    def z_over_di0(self, z):     return _div(z, self.di0)
    def z_over_rho_i0(self, z):  return _div(z, self.rho_i0)
    def t_wci0(self, t):         return _mul(t, self.wci0)
    def t_wpe(self, t):          return _mul(t, self.wpe)
    def v_over_Csab(self, v):    return _div(v, self.Cs_ab)
    def v_over_vsh(self, v):     return _div(v, self.vsh_model)
    def n_over_n0(self, n):      return _div(n, self.namb)   # upstream-normalized (paper's n/n_e0)
    def b_over_b0(self, b):      return _div(b, self.B0)

    def report(self) -> dict:
        d = asdict(self)
        d.pop("_meta", None)
        # convenient dimensionless echoes
        d["Cs_ab_over_c"] = self.Cs_ab / C
        d["vA_over_c"] = self.vA / C
        d["vp_model_over_c"] = self.vp_model / C
        d["vsh_model_over_c"] = self.vsh_model / C
        d["di_over_de"] = self.di / self.de
        d["wci0_inv_over_wpe_inv"] = self.wci0_inv * self.wpe
        d["rho_i0_over_de"] = self.rho_i0 / self.de
        if self.mfp_ei_ab is not None:
            d["mfp_ei_ab_over_de"] = self.mfp_ei_ab / self.de
            d["mfp_ei_ab_over_di0"] = self.mfp_ei_ab / self.di0
            d["rho_e_ab_over_de"] = (self.vte_ab / self.wce_ab) / self.de
        return d

    def pretty(self) -> str:
        r = self.report()
        lines = [f"Scales for run {self._meta.get('run_id', '?')} "
                 f"(tier {self._meta.get('tier', '?')}):"]
        order = ["n0", "mass_ratio", "wpe", "de", "di", "Cs_ab_over_c", "vA_over_c",
                 "B0", "wci0_inv", "rho_i0", "vp_model_over_c", "vsh_model_over_c",
                 "MA", "Mms", "beta_ab", "beta_0", "dt", "dt_wpe", "n_cell",
                 "domain_halfwidth", "steps_per_wci0"]
        if self.coulomb_log is not None:
            order += ["Te_ab_eV", "coulomb_log_nrl", "coulomb_log", "nu_ei_ab",
                      "nu_ei_dt", "mfp_ei_ab_over_de", "mfp_ei_ab_over_di0", "lambda_ab"]
        for k in order:
            lines.append(f"  {k:24s} = {r[k]:.5g}")
        return "\n".join(lines)


def _div(x, y):
    try:
        return [xi / y for xi in x]
    except TypeError:
        return x / y


def _mul(x, y):
    try:
        return [xi * y for xi in x]
    except TypeError:
        return x * y


def derive(cfg: dict) -> Scales:
    """Compute all derived scales from a loaded config dict (see runs/R1/config.yaml)."""
    ref = cfg["reference"]
    pis = cfg["plasma"]["piston"]
    amb = cfg["plasma"]["ambient"]
    geo = cfg["geometry"]
    num = cfg["numerics"]
    mdl = cfg["model"]

    n0 = float(ref["n0"])
    mass_ratio = float(ref["mass_ratio"])
    mi = mass_ratio * ME

    # ablation scales (built on n0 = n_e,ab)
    wpe = math.sqrt(n0 * QE * QE / (EPS0 * ME))
    de = C / wpe
    di = de * math.sqrt(mass_ratio)

    theta_e = float(pis["theta_e_heat"])
    theta_0 = float(amb["theta_0"])
    Cs_ab = math.sqrt(theta_e / mass_ratio) * C           # sqrt(Z k_B T_e,ab / m_i), Z=1
    Cs0 = math.sqrt(theta_0 / mass_ratio) * C
    t_ab = di / Cs_ab

    nt = float(pis["density_over_n0"]) * n0
    namb = float(amb["density_over_n0"]) * n0

    # upstream scales
    wpe0 = math.sqrt(namb * QE * QE / (EPS0 * ME))
    de0 = C / wpe0
    di0 = de0 * math.sqrt(mass_ratio)

    vA = float(cfg["field"]["vA_over_c"]) * C
    B0 = vA * math.sqrt(MU0 * namb * mi)
    wci0 = QE * B0 / mi
    wci0_inv = 1.0 / wci0

    vp_model = float(mdl["vp_over_c"]) * C
    vsh_model = float(mdl["vsh_over_Csab"]) * Cs_ab
    rho_i0 = vp_model / wci0

    MA = vsh_model / vA
    Mms = vsh_model / math.sqrt(vA * vA + Cs0 * Cs0)
    beta_ab = 2.0 * MU0 * n0 * theta_e * ME_C2_J / (B0 * B0)
    beta_0 = 2.0 * MU0 * namb * theta_0 * ME_C2_J / (B0 * B0)

    # ablation electron kinetics (used by the collision block below)
    Te_ab_eV = theta_e * ME_C2_EV
    vte_ab = math.sqrt(theta_e) * C
    wce_ab = QE * B0 / ME

    dz = float(geo["dz_over_de"]) * de
    dt = float(num["cfl"]) * dz / C                       # 1D CFL
    dt_wpe = dt * wpe
    halfwidth_de = float(geo["domain_halfwidth_de"])
    domain_halfwidth = halfwidth_de * de
    # one-sided [0, half] has half the cells of symmetric [-half, +half] at same dz.
    _span = 1.0 if geo.get("layout", "symmetric") == "one_sided" else 2.0
    n_cell = int(round(_span * halfwidth_de / float(geo["dz_over_de"])))
    steps_per_wci0 = wci0_inv / dt

    coulomb_log_nrl = _lnL_nrl(n0, Te_ab_eV)
    coulomb_log = nu_ei_ab = mfp_ei_ab = lambda_ab = nu_ei_dt = None
    if cfg.get("collisions"):
        coulomb_log = coulomb_log_for(cfg["collisions"], n0, Te_ab_eV, vte_ab, de, wce_ab)
        nu_ei_ab = nu_ei(n0, Te_ab_eV, coulomb_log)
        mfp_ei_ab = vte_ab / nu_ei_ab
        lambda_ab = wce_ab / nu_ei_ab
        nu_ei_dt = nu_ei_ab * dt

    return Scales(
        n0=n0, mass_ratio=mass_ratio, mi=mi,
        wpe=wpe, de=de, di=di, t_ab=t_ab, Cs_ab=Cs_ab, nt=nt,
        namb=namb, de0=de0, di0=di0, Cs0=Cs0, B0=B0, vA=vA,
        wci0=wci0, wci0_inv=wci0_inv, rho_i0=rho_i0,
        vp_model=vp_model, vsh_model=vsh_model,
        MA=MA, Mms=Mms, beta_ab=beta_ab, beta_0=beta_0,
        dz=dz, dt=dt, dt_wpe=dt_wpe, n_cell=n_cell,
        domain_halfwidth=domain_halfwidth, steps_per_wci0=steps_per_wci0,
        Te_ab_eV=Te_ab_eV, vte_ab=vte_ab, wce_ab=wce_ab,
        coulomb_log_nrl=coulomb_log_nrl, coulomb_log=coulomb_log,
        nu_ei_ab=nu_ei_ab, mfp_ei_ab=mfp_ei_ab, lambda_ab=lambda_ab, nu_ei_dt=nu_ei_dt,
        _meta=dict(cfg.get("meta", {})),
    )


# --------------------------------------------------------------------------- #
# Coulomb collisions
# --------------------------------------------------------------------------- #
def nu_ei(n_e: float, Te_eV: float, coulomb_log: float) -> float:
    """NRL electron-ion collision rate [s^-1] for ``n_e`` [m^-3], ``Te_eV`` [eV]."""
    return NU_EI_NRL * (n_e / 1e6) * coulomb_log * Te_eV ** -1.5


def _lnL_nrl(n_e: float, Te_eV: float) -> float:
    """The PHYSICAL Coulomb logarithm (NRL e-i, T_e > 10 eV): 24 - ln(sqrt(n)/T)."""
    return 24.0 - math.log(math.sqrt(n_e / 1e6) / Te_eV)


def coulomb_log_for(coll: dict, n0: float, Te_ab_eV: float, vte_ab: float,
                    de: float, wce_ab: float) -> float:
    """Resolve the ``collisions.target`` block to the lnLambda the deck must use.

    Because :mod:`kinshock` runs the paper's dimensionless problem at the REAL
    speed of light, theta_e = 0.078 means T_e,ab ~ 40 keV rather than the paper's
    ~470 eV, and nu_ei ~ T^-3/2 makes the plasma collisionless by ~4-9 orders of
    magnitude at any attainable density (see RESULTS / R1_warm header). lnLambda
    is therefore used as the single knob that dials in a chosen collisionality:
    nu_ei is exactly linear in it, so requesting a target mean free path (or the
    paper's magnetization ratio) inverts to one number.

    ``target.quantity`` is one of

    ``mfp_over_de``
        electron mean free path v_te,ab/nu_ei,ab in ablation skin depths.
    ``lambda_ab``
        the paper's collisionality lambda_ab = omega_ce,ab/nu_ei,ab (Schaeffer
        2020 Table I uses 20; note this equals mfp/rho_e,ab, NOT mfp/d_e).
    ``coulomb_log``
        set lnLambda directly (``value: physical`` uses the NRL formula).
    """
    tgt = coll.get("target") or {}
    quantity = tgt.get("quantity", "coulomb_log")
    value = tgt.get("value", "physical")
    # rate per unit lnLambda: nu_ei = rate1 * lnLambda
    rate1 = nu_ei(n0, Te_ab_eV, 1.0)
    if quantity == "coulomb_log":
        return _lnL_nrl(n0, Te_ab_eV) if value == "physical" else float(value)
    if quantity == "mfp_over_de":
        return (vte_ab / (float(value) * de)) / rate1
    if quantity == "lambda_ab":
        return (wce_ab / float(value)) / rate1
    raise ValueError(f"unknown collisions.target.quantity {quantity!r}; expected one of "
                     "'mfp_over_de', 'lambda_ab', 'coulomb_log'")


def collision_pairs(cfg: dict) -> list[tuple[str, str]]:
    """The species pairs to collide, from ``collisions.pairs``.

    ``all`` (default) = every unordered pair of the run's species INCLUDING the
    intra-species self-pairs, i.e. the full Takizuka-Abe treatment the paper
    applies. An explicit list of ``[a, b]`` pairs overrides it.
    """
    coll = cfg.get("collisions") or {}
    names = list(cfg["species"])
    spec = coll.get("pairs", "all")
    if spec == "all":
        return [(a, b) for i, a in enumerate(names) for b in names[i:]]
    pairs = []
    for p in spec:
        a, b = (p.split() if isinstance(p, str) else list(p))
        for s in (a, b):
            if s not in names:
                raise ValueError(f"collisions.pairs references unknown species {s!r}")
        pairs.append((a, b))
    return pairs
