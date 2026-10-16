//! Quasi-static quasi-star envelope (port of `envelope.jl`): hydrostatic, constant luminosity L
//! from the central black hole, local mixing-length convection (Böhm-Vitense; Kippenhahn &
//! Weigert cubic), point mass M_BH at the centre. The independent variable is x = ln P and the
//! envelope is integrated inward from the outer boundary with an adaptive Dormand–Prince 5(4)
//! method.
use crate::constants::*;
use crate::eos::{eos_state, rho_from_p_t};
use crate::opacity::KapTable;
use std::f64::consts::PI;

/// Model parameters (cgs). See [`Params::quasi_star`] for the defaults of the text.
#[derive(Clone, Copy, Debug)]
pub struct Params<'a> {
    /// Total mass M (black hole + envelope), g.
    pub m: f64,
    /// Black-hole (core) mass M_BH, g.
    pub m_bh: f64,
    /// Luminosity through the envelope L, erg s⁻¹.
    pub l: f64,
    /// Mixing-length parameter α_MLT (ℓ = α H_P).
    pub alpha: f64,
    /// Inner radius in units of R, ε_in, at which m must equal M_BH.
    pub eps_in: f64,
    /// Opacity table.
    pub tab: &'a KapTable,
}

impl<'a> Params<'a> {
    /// The parameters of the text for total mass `m_msun` (M_sun) and black-hole fraction
    /// `f_bh = M_BH/M`: L = L_Edd(κ_es) = 4π c G M/(0.2(1+X)) (see [`crate::constants::l_edd`]),
    /// α_MLT = 0.8, ε_in = 10⁻³.
    pub fn quasi_star(tab: &'a KapTable, m_msun: f64, f_bh: f64) -> Self {
        let m = m_msun * MSUN;
        Params { m, m_bh: f_bh * m, l: l_edd(m), alpha: 0.8, eps_in: 1e-3, tab }
    }
}

/// Local structure at one point of the envelope, see [`mlt`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mlt {
    /// Actual temperature gradient ∇ = d ln T/d ln P.
    pub nabla: f64,
    /// Radiative gradient ∇_rad.
    pub nabla_rad: f64,
    /// Adiabatic gradient ∇_ad.
    pub nabla_ad: f64,
    /// Convective velocity v_c (cm s⁻¹; 0 if radiative).
    pub v: f64,
    /// Adiabatic sound speed c_s.
    pub cs: f64,
    /// Density ρ.
    pub rho: f64,
    /// Rosseland mean opacity κ.
    pub kappa: f64,
    /// Pressure scale height H_P = P/(ρ g).
    pub hp: f64,
    /// Local Eddington factor Γ = κ L/(4π c G m).
    pub gamma: f64,
    /// β = P_gas/P.
    pub beta: f64,
    /// First adiabatic exponent Γ₁.
    pub gamma1: f64,
}

/// Equations (eq-mlt-cubic) and (eq-mlt-grad) of the text: the MLT temperature gradient and
/// convective velocity at (P, T, r, m).
///
/// `∇_rad = 3 κ L P / (16π a c G m T⁴)`, `Γ = κ L/(4π c G m)`, `g = G m/r²`, `H_P = P/(ρ g)`.
/// If `∇_rad ≤ ∇_ad` the point is radiative and ∇ = ∇_rad. Otherwise, with ℓ = α H_P,
/// `U = 3 a c T³/(c_P ρ² κ ℓ²) · (8 H_P/(g δ))^{1/2}` and `W = ∇_rad − ∇_ad`, the cubic
/// (eq-mlt-cubic) `(ξ − U)³ + (8U/9)(ξ² − U² − W) = 0` is solved for ξ by bisection on
/// [U, (W + U²)^{1/2}] (to 10⁻¹⁴ relative, as in the Julia code), and (eq-mlt-grad)
/// `∇ = ξ² − U² + ∇_ad`, `v_c² = g δ (ξ − U)² ℓ²/(8 H_P)`.
///
/// Returns `None` if there is no gas state (P ≤ P_rad(T)).
pub fn mlt(p: f64, t: f64, r: f64, m: f64, l: f64, par: &Params) -> Option<Mlt> {
    let rho = rho_from_p_t(p, t);
    if rho.is_nan() {
        return None;
    }
    let th = eos_state(rho, t);
    let kap = par.tab.kappa_r(rho, t);
    let g = G * m / (r * r);
    let hp = p / (rho * g);
    let nabla_rad = 3.0 * kap * l * p / (16.0 * PI * ARAD * C * G * m * t.powi(4));
    let gam = kap * l / (4.0 * PI * C * G * m);
    let mut out = Mlt {
        nabla: nabla_rad,
        nabla_rad,
        nabla_ad: th.nabla_ad,
        v: 0.0,
        cs: th.c_s,
        rho,
        kappa: kap,
        hp,
        gamma: gam,
        beta: th.beta,
        gamma1: th.gamma1,
    };
    if nabla_rad <= th.nabla_ad {
        return Some(out);
    }
    let ell = par.alpha * hp;
    let u = 3.0 * ARAD * C * t.powi(3) / (th.c_p * rho * rho * kap * ell * ell)
        * (8.0 * hp / (g * th.delta)).sqrt();
    let w = nabla_rad - th.nabla_ad;
    let f = |xi: f64| (xi - u).powi(3) + (8.0 * u / 9.0) * (xi * xi - u * u - w);
    let (mut a, mut b) = (u, (w + u * u).sqrt());
    for _ in 0..200 {
        let mid = 0.5 * (a + b);
        if f(mid) > 0.0 {
            b = mid;
        } else {
            a = mid;
        }
        if b - a < 1e-14 * b {
            break;
        }
    }
    let xi = 0.5 * (a + b);
    out.nabla = xi * xi - u * u + th.nabla_ad;
    out.v = (g * th.delta * (xi - u).powi(2) * ell * ell / (8.0 * hp)).max(0.0).sqrt();
    Some(out)
}

/// Equation (eq-structure) of the text: right-hand side of the structure equations in x = ln P,
/// `dr/dx = −P r²/(G m ρ)`, `dm/dx = 4π r² ρ dr/dx`, `d ln T/dx = ∇` (∇ from [`mlt`]), for the
/// state y = (r, m, ln T). Returns `None` where the EOS has no state.
pub fn structure_rhs(x: f64, y: [f64; 3], par: &Params) -> Option<([f64; 3], Mlt)> {
    let [r, m, lt] = y;
    let p = x.exp();
    let t = lt.exp();
    let q = mlt(p, t, r, m, par.l, par)?;
    if q.nabla.is_nan() {
        return None;
    }
    let drdx = -p * r * r / (G * m * q.rho);
    Some(([drdx, 4.0 * PI * r * r * q.rho * drdx, q.nabla], q))
}

/// How an inward integration ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// r reached r_in = ε_in R with m > M_BH (envelope mass left over).
    ReachedRin,
    /// m fell to M_BH before r_in (envelope mass exhausted).
    MassExhausted,
    /// No EOS state, step size underflow, or too many steps.
    Failed,
}

impl Status {
    /// Numeric code used by the WebAssembly exports: 1 reached r_in, 2 mass exhausted, 0 failed.
    pub fn code(self) -> f64 {
        match self {
            Status::ReachedRin => 1.0,
            Status::MassExhausted => 2.0,
            Status::Failed => 0.0,
        }
    }
}

/// Envelope profile saved at the start of every accepted step (outside in).
#[derive(Clone, Debug, Default)]
pub struct Profile {
    /// Radius r (cm).
    pub r: Vec<f64>,
    /// Enclosed mass m (g), including M_BH.
    pub m: Vec<f64>,
    /// Pressure P (dyn cm⁻²).
    pub p: Vec<f64>,
    /// Temperature T (K).
    pub t: Vec<f64>,
    /// Density ρ (g cm⁻³).
    pub rho: Vec<f64>,
    /// ∇ = d ln T/d ln P.
    pub nabla: Vec<f64>,
    /// ∇_ad.
    pub nabla_ad: Vec<f64>,
    /// ∇_rad.
    pub nabla_rad: Vec<f64>,
    /// Convective velocity v_c (cm s⁻¹).
    pub v_conv: Vec<f64>,
    /// Sound speed c_s (cm s⁻¹).
    pub cs: Vec<f64>,
    /// Rosseland opacity κ (cm² g⁻¹).
    pub kappa: Vec<f64>,
    /// Local Eddington factor Γ = κ L/(4π c G m).
    pub gamma: Vec<f64>,
    /// β = P_gas/P.
    pub beta: Vec<f64>,
    /// Pressure scale height H_P (cm).
    pub hp: Vec<f64>,
    /// First adiabatic exponent Γ₁.
    pub gamma1: Vec<f64>,
}

impl Profile {
    /// Number of saved points.
    pub fn len(&self) -> usize {
        self.r.len()
    }
    /// True if no point was saved.
    pub fn is_empty(&self) -> bool {
        self.r.is_empty()
    }
    /// The fields in the order used by the WebAssembly export `envelope_profile`:
    /// r, m, P, T, ρ, ∇, ∇_ad, ∇_rad, v_conv, c_s, κ, Γ, β, H_P, Γ₁.
    pub fn fields(&self) -> [&Vec<f64>; 15] {
        [
            &self.r, &self.m, &self.p, &self.t, &self.rho, &self.nabla, &self.nabla_ad,
            &self.nabla_rad, &self.v_conv, &self.cs, &self.kappa, &self.gamma, &self.beta,
            &self.hp, &self.gamma1,
        ]
    }
    fn push(&mut self, y: &[f64; 3], x: f64, q: &Mlt) {
        self.r.push(y[0]);
        self.m.push(y[1]);
        self.p.push(x.exp());
        self.t.push(y[2].exp());
        self.rho.push(q.rho);
        self.nabla.push(q.nabla);
        self.nabla_ad.push(q.nabla_ad);
        self.nabla_rad.push(q.nabla_rad);
        self.v_conv.push(q.v);
        self.cs.push(q.cs);
        self.kappa.push(q.kappa);
        self.gamma.push(q.gamma);
        self.beta.push(q.beta);
        self.hp.push(q.hp);
        self.gamma1.push(q.gamma1);
    }
}

/// Result of [`integrate_inward`].
#[derive(Clone, Debug)]
pub struct Integration {
    /// How the integration ended.
    pub status: Status,
    /// m(r_in) (linear interpolation within the last step) if r_in was reached; M_BH if the mass
    /// was exhausted; NaN on failure.
    pub m_rin: f64,
    /// Saved profile (if requested).
    pub profile: Option<Profile>,
}

// Dormand–Prince 5(4) tableau
const DP_C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
const DP_A: [[f64; 6]; 7] = [
    [0.0; 6],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
];
const DP_B: [f64; 7] = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0];
const DP_BS: [f64; 7] = [
    5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0,
];

#[inline]
fn tadd(y: &[f64; 3], ks: &[[f64; 3]; 7], coeffs: &[f64], h: f64) -> [f64; 3] {
    let mut o = [0.0; 3];
    for i in 0..3 {
        let mut s = 0.0;
        for (j, cj) in coeffs.iter().enumerate() {
            s += cj * ks[j][i];
        }
        o[i] = y[i] + h * s;
    }
    o
}

/// Relative tolerance of the inward integration (as in the Julia code).
pub const RTOL: f64 = 1e-7;

/// Equation (eq-structure) of the text integrated inward in x = ln P from (P_b, T_b) at r = R,
/// m = M, with adaptive Dormand–Prince 5(4) (rtol 10⁻⁷; error norm
/// `max_i |y5_i − y4_i| / (rtol (|y5_i| + 10⁻³⁰) + [i = ln T] rtol)`, step factor
/// `min(4, 0.9 err^{−1/5} + 10⁻³)` on acceptance and `max(0.1, 0.9 err^{−1/5})` on rejection,
/// first step 0.01; exactly the Julia `integrate_inward`).
///
/// Stops when r crosses r_in = ε_in R (status [`Status::ReachedRin`], m(r_in) linearly interpolated
/// within the step) or when m ≤ M_BH ([`Status::MassExhausted`]). The shooting condition of the
/// text is m(r_in) = M_BH. If `save`, the [`Profile`] at the start of each accepted step is returned.
pub fn integrate_inward(rr: f64, pb: f64, tb: f64, par: &Params, save: bool) -> Integration {
    let maxsteps = 200_000;
    let mut x = pb.ln();
    let mut y = [rr, par.m, tb.ln()];
    let mut h = 0.01;
    let rin = par.eps_in * rr;
    let mut prof = if save { Some(Profile::default()) } else { None };
    let fail = |prof| Integration { status: Status::Failed, m_rin: f64::NAN, profile: prof };
    for _ in 0..maxsteps {
        let (k1, q) = match structure_rhs(x, y, par) {
            Some(v) => v,
            None => return fail(prof),
        };
        if let Some(p) = prof.as_mut() {
            p.push(&y, x, &q);
        }
        loop {
            let mut ks = [[0.0; 3]; 7];
            ks[0] = k1;
            let mut bad = false;
            for s in 1..7 {
                let ys = tadd(&y, &ks, &DP_A[s][..s], h);
                if ys[0] <= 0.0 {
                    bad = true;
                    break;
                }
                match structure_rhs(x + DP_C[s] * h, ys, par) {
                    Some((k, _)) if !k[0].is_nan() => ks[s] = k,
                    _ => {
                        bad = true;
                        break;
                    }
                }
            }
            if bad {
                h /= 4.0;
                if h < 1e-12 {
                    return fail(prof);
                }
                continue;
            }
            let y5 = tadd(&y, &ks, &DP_B, h);
            let y4 = tadd(&y, &ks, &DP_BS, h);
            let mut err: f64 = 0.0;
            for i in 0..3 {
                let e = (y5[i] - y4[i]).abs()
                    / (RTOL * (y5[i].abs() + 1e-30) + if i == 2 { RTOL } else { 0.0 });
                if e.is_nan() || e > err {
                    err = e; // NaN propagates and rejects the step
                }
            }
            if err <= 1.0 && y5[0] > 0.0 {
                if y5[0] <= rin {
                    let w = (y[0] - rin) / (y[0] - y5[0]);
                    return Integration {
                        status: Status::ReachedRin,
                        m_rin: y[1] + w * (y5[1] - y[1]),
                        profile: prof,
                    };
                }
                if y5[1] <= par.m_bh {
                    return Integration { status: Status::MassExhausted, m_rin: par.m_bh, profile: prof };
                }
                x += h;
                y = y5;
                h *= (0.9 * err.powf(-0.2) + 1e-3).min(4.0);
                break;
            } else {
                h *= (0.9 * err.powf(-0.2)).max(0.1);
                if h < 1e-12 {
                    return fail(prof);
                }
            }
        }
    }
    fail(prof)
}
