//! Equation of state of the Z = 0 envelope: an ideal gas of H, H⁺, He, He⁺, He²⁺ and e⁻ in Saha
//! (LTE) ionization equilibrium plus blackbody radiation (port of `eos.jl`). No degeneracy,
//! no pressure ionization and no molecules.
//!
//! The thermodynamic derivatives are exact: the code is generic over a forward-mode dual number
//! seeded with ∂/∂ln ρ and ∂/∂ln T, and the implicit charge-balance solve is differentiated by
//! Newton steps taken in the dual type from the converged `f64` root (as the Julia code does
//! with ForwardDiff).
use crate::constants::*;
use crate::dual::{Dual, Num};
use std::f64::consts::PI;

/// Number densities (cm⁻³) of the Saha equilibrium at given (ρ, T).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Saha {
    /// Free electrons n_e.
    pub ne: f64,
    /// Neutral hydrogen n(H⁰).
    pub n_h0: f64,
    /// Protons n(H⁺).
    pub n_hp: f64,
    /// Neutral helium n(He⁰).
    pub n_he0: f64,
    /// Singly ionized helium n(He⁺).
    pub n_he1: f64,
    /// Doubly ionized helium n(He²⁺).
    pub n_he2: f64,
    /// Total hydrogen nuclei n_H = Xρ/m_H.
    pub n_h: f64,
    /// Total helium nuclei n_He = Yρ/m_He.
    pub n_he: f64,
    /// Ionized fraction of hydrogen x_H = n(H⁺)/n_H.
    pub x_h: f64,
}

impl Saha {
    /// Fraction of helium that is singly ionized, y_He1 = n(He⁺)/n_He.
    pub fn y_he1(&self) -> f64 {
        self.n_he1 / self.n_he
    }
    /// Fraction of helium that is doubly ionized, y_He2 = n(He²⁺)/n_He.
    pub fn y_he2(&self) -> f64 {
        self.n_he2 / self.n_he
    }
}

#[derive(Clone, Copy)]
pub(crate) struct SahaG<T> {
    pub ne: T,
    pub n_h0: T,
    pub n_hp: T,
    pub n_he0: T,
    pub n_he1: T,
    pub n_he2: T,
    pub n_h: T,
    pub n_he: T,
    pub x_h: T,
}

/// Right-hand sides of the three Saha equations, λ_e⁻³ e^{−χ/kT} (with statistical weight 4 for He⁰→He⁺).
#[inline]
fn saha_rates<T: Num>(t: T) -> (T, T, T) {
    let kt = t * KB;
    let l = (t * (2.0 * PI * ME * KB) / (H * H)).powf(1.5);
    let sh = l * (-(kt.recip() * CHI_H)).exp();
    let s1 = l * 4.0 * (-(kt.recip() * CHI_HE1)).exp();
    let s2 = l * (-(kt.recip() * CHI_HE2)).exp();
    (sh, s1, s2)
}

/// Ionization fractions for a trial n_e, and the charge sum n_H x_H + n_He (y1 + 2 y2).
#[inline]
fn fractions<T: Num>(ne: T, nh: T, nhe: T, sh: T, s1: T, s2: T) -> (T, T, T, T, T) {
    let xh = sh / (ne + sh);
    let r1 = s1 / ne;
    let r2 = s2 / ne;
    let f0 = (r1 + r1 * r2 + 1.0).recip();
    let y1 = r1 * f0;
    let y2 = r1 * r2 * f0;
    let q = nh * xh + nhe * (y1 + y2 * 2.0);
    (xh, f0, y1, y2, q)
}

/// Charge-balance residual in ln n_e: g = ln(charge sum) − ln n_e (strictly decreasing in ln n_e).
#[inline]
fn charge_residual<T: Num>(lne: T, nh: T, nhe: T, sh: T, s1: T, s2: T) -> T {
    fractions(lne.exp(), nh, nhe, sh, s1, s2).4.ln() - lne
}

/// Float64 charge-balance solve: safeguarded Newton in ln n_e inside the bisection bracket of the
/// Julia code, [ln(10⁻³⁰ n_e,max + 10⁻³⁰⁰), ln n_e,max], n_e,max = n_H + 2 n_He.
fn ln_ne_f64(nh: f64, nhe: f64, sh: f64, s1: f64, s2: f64) -> f64 {
    let nemax = nh + 2.0 * nhe;
    let mut lo = (nemax * 1e-30 + 1e-300).ln();
    let mut hi = nemax.ln();
    let c = |v: f64| Dual::<1>::new(v, [0.0]);
    let (dnh, dnhe, dsh, ds1, ds2) = (c(nh), c(nhe), c(sh), c(s1), c(s2));
    let mut x = hi; // start from full ionization
    for _ in 0..200 {
        let g = charge_residual(Dual::new(x, [1.0]), dnh, dnhe, dsh, ds1, ds2);
        if g.v == 0.0 {
            return x;
        }
        // g decreasing: g > 0 means the root lies above x
        if g.v > 0.0 {
            lo = lo.max(x);
        } else {
            hi = hi.min(x);
        }
        let mut xn = x - g.v / g.d[0];
        if !(xn > lo && xn < hi) || !g.d[0].is_finite() {
            xn = 0.5 * (lo + hi);
        }
        let dx = (xn - x).abs();
        x = xn;
        if dx < 1e-13 * (1.0 + x.abs()) || hi - lo < 1e-14 {
            break;
        }
    }
    x
}

/// Saha state in any number type: f64 solve, then (for duals) two Newton polish steps of the
/// charge balance in the dual type, which make the first derivatives exact.
pub(crate) fn saha_g<T: Num>(rho: T, t: T) -> SahaG<T> {
    let nh = rho * (X / MH);
    let nhe = rho * (Y / MHE);
    let (sh, s1, s2) = saha_rates(t);
    let x0 = ln_ne_f64(nh.val(), nhe.val(), sh.val(), s1.val(), s2.val());
    let mut lne = T::cst(x0);
    if T::IS_DUAL {
        let c = |v: f64| Dual::<1>::new(v, [0.0]);
        for _ in 0..2 {
            let g = charge_residual(lne, nh, nhe, sh, s1, s2);
            let gp = charge_residual(
                Dual::new(lne.val(), [1.0]),
                c(nh.val()),
                c(nhe.val()),
                c(sh.val()),
                c(s1.val()),
                c(s2.val()),
            )
            .d[0];
            lne = lne - g / gp;
        }
    }
    let ne = lne.exp();
    let (xh, f0, y1, y2, _) = fractions(ne, nh, nhe, sh, s1, s2);
    SahaG {
        ne,
        n_h0: nh * (-xh + 1.0),
        n_hp: nh * xh,
        n_he0: nhe * f0,
        n_he1: nhe * y1,
        n_he2: nhe * y2,
        n_h: nh,
        n_he: nhe,
        x_h: xh,
    }
}

/// Equation (eq-saha) of the text: the Saha equations
/// `n(H⁺) n_e / n(H⁰) = λ_e⁻³ exp(−χ_H/kT)`, `n(He⁺) n_e / n(He⁰) = 4 λ_e⁻³ exp(−χ_HeI/kT)`,
/// `n(He²⁺) n_e / n(He⁺) = λ_e⁻³ exp(−χ_HeII/kT)`, with `λ_e⁻³ = (2π m_e k T / h²)^{3/2}`,
/// closed by charge neutrality `n_e = n(H⁺) + n(He⁺) + 2 n(He²⁺)`.
///
/// The charge balance is solved by safeguarded Newton iteration in ln n_e inside the bisection
/// bracket of the Julia code (to ~10⁻¹³ in ln n_e).
pub fn saha(rho: f64, t: f64) -> Saha {
    let s = saha_g(rho, t);
    Saha {
        ne: s.ne,
        n_h0: s.n_h0,
        n_hp: s.n_hp,
        n_he0: s.n_he0,
        n_he1: s.n_he1,
        n_he2: s.n_he2,
        n_h: s.n_h,
        n_he: s.n_he,
        x_h: s.x_h,
    }
}

#[inline]
fn pgas_g<T: Num>(s: &SahaG<T>, t: T) -> T {
    (s.n_h + s.n_he + s.ne) * t * KB
}
#[inline]
fn prad_g<T: Num>(t: T) -> T {
    let t2 = t * t;
    t2 * t2 * (ARAD / 3.0)
}
#[inline]
fn eint_g<T: Num>(s: &SahaG<T>, rho: T, t: T) -> T {
    let t2 = t * t;
    ((s.n_h + s.n_he + s.ne) * t * (1.5 * KB)
        + s.n_hp * CHI_H
        + s.n_he1 * CHI_HE1
        + s.n_he2 * (CHI_HE1 + CHI_HE2))
        / rho
        + t2 * t2 * ARAD / rho
}
#[inline]
fn entropy_g<T: Num>(s: &SahaG<T>, rho: T, t: T) -> T {
    let kt = t * KB;
    let q = |m: f64| (kt * (2.0 * PI * m) / (H * H)).powf(1.5);
    let term = |n: T, g: f64, m: f64| {
        if n.val() > 1e-280 {
            n * ((q(m) * g / n).ln() + 2.5)
        } else {
            T::cst(0.0)
        }
    };
    let sum = term(s.n_h0, 2.0, MH)
        + term(s.n_hp, 1.0, MH)
        + term(s.n_he0, 1.0, MHE)
        + term(s.n_he1, 2.0, MHE)
        + term(s.n_he2, 1.0, MHE)
        + term(s.ne, 2.0, ME);
    sum * KB / rho + t * t * t * (4.0 * ARAD) / (rho * 3.0)
}

/// Radiation pressure `P_rad = a T⁴ / 3` (part of Equation (eq-pressure) of the text).
pub fn prad(t: f64) -> f64 {
    prad_g(t)
}

/// Gas pressure `P_gas = (n_H + n_He + n_e) k T` (part of Equation (eq-pressure) of the text).
pub fn pgas(rho: f64, t: f64) -> f64 {
    pgas_g(&saha_g(rho, t), t)
}

/// Equation (eq-pressure) of the text: total pressure `P = (n_H + n_He + n_e) k T + a T⁴ / 3`.
pub fn pressure(rho: f64, t: f64) -> f64 {
    pgas(rho, t) + prad(t)
}

/// Specific internal energy (erg g⁻¹): translational `(3/2)(n_H + n_He + n_e) k T / ρ`, plus
/// ionization `[n(H⁺) χ_H + n(He⁺) χ_HeI + n(He²⁺)(χ_HeI + χ_HeII)] / ρ`, plus radiation `a T⁴ / ρ`.
pub fn internal_energy(rho: f64, t: f64) -> f64 {
    eint_g(&saha_g(rho, t), rho, t)
}

/// Equation (eq-entropy) of the text: specific entropy (erg g⁻¹ K⁻¹), the Sackur–Tetrode sum over
/// the six species plus radiation,
/// `s = (k/ρ) Σ_i n_i [5/2 + ln( (g_i/n_i) (2π m_i k T / h²)^{3/2} )] + 4 a T³ / (3ρ)`,
/// with g = 2, 1, 1, 2, 1, 2 for H⁰, H⁺, He⁰, He⁺, He²⁺, e⁻ (species with n < 10⁻²⁸⁰ are skipped).
pub fn entropy(rho: f64, t: f64) -> f64 {
    entropy_g(&saha_g(rho, t), rho, t)
}

/// Full thermodynamic state at (ρ, T), see [`eos_state`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EosState {
    /// Density ρ (g cm⁻³).
    pub rho: f64,
    /// Temperature T (K).
    pub t: f64,
    /// Total pressure P (dyn cm⁻²).
    pub p: f64,
    /// Gas pressure P_gas.
    pub pgas: f64,
    /// Radiation pressure a T⁴/3.
    pub prad: f64,
    /// Specific internal energy e (erg g⁻¹).
    pub e: f64,
    /// Specific entropy s (erg g⁻¹ K⁻¹).
    pub s: f64,
    /// χ_ρ = ∂ln P/∂ln ρ|_T.
    pub chi_rho: f64,
    /// χ_T = ∂ln P/∂ln T|_ρ.
    pub chi_t: f64,
    /// s_ρ = ∂s/∂ln ρ|_T.
    pub s_rho: f64,
    /// s_T = ∂s/∂ln T|_ρ.
    pub s_t: f64,
    /// c_V = ∂e/∂T|_ρ.
    pub c_v: f64,
    /// Adiabatic gradient ∇_ad.
    pub nabla_ad: f64,
    /// First adiabatic exponent Γ₁.
    pub gamma1: f64,
    /// δ = −∂ln ρ/∂ln T|_P = χ_T/χ_ρ.
    pub delta: f64,
    /// Specific heat at constant pressure c_P.
    pub c_p: f64,
    /// Adiabatic sound speed c_s = (Γ₁ P/ρ)^{1/2}.
    pub c_s: f64,
    /// β = P_gas/P = 1 − a T⁴/(3P).
    pub beta: f64,
    /// The underlying Saha state.
    pub saha: Saha,
}

/// Equation (eq-thermo) of the text: the full EOS state at (ρ, T) with
/// `∇_ad = 1/(χ_T − χ_ρ s_T/s_ρ)`, `Γ₁ = χ_ρ − χ_T s_ρ/s_T`, `δ = χ_T/χ_ρ`, `c_P = s_T − s_ρ δ`,
/// where `χ_ρ = ∂ln P/∂ln ρ|_T`, `χ_T = ∂ln P/∂ln T|_ρ`, `s_ρ = ∂s/∂ln ρ|_T`, `s_T = ∂s/∂ln T|_ρ`;
/// also `c_s = (Γ₁ P/ρ)^{1/2}` and `β = 1 − a T⁴/(3P)`.
///
/// P, e, s are Equations (eq-pressure) and (eq-entropy); the derivatives are exact (forward-mode
/// dual numbers in (ln ρ, ln T) carried through one Saha solve, Equation (eq-saha)).
pub fn eos_state(rho: f64, t: f64) -> EosState {
    let rd = Dual::<2>::new(rho, [rho, 0.0]); // ∂/∂ln ρ
    let td = Dual::<2>::new(t, [0.0, t]); // ∂/∂ln T
    let sd = saha_g(rd, td);
    let pd = pgas_g(&sd, td) + prad_g(td);
    let ed = eint_g(&sd, rd, td);
    let s_d = entropy_g(&sd, rd, td);
    let p = pd.v;
    let chi_rho = pd.d[0] / p;
    let chi_t = pd.d[1] / p;
    let (s_rho, s_t) = (s_d.d[0], s_d.d[1]);
    let nabla_ad = 1.0 / (chi_t - chi_rho * s_t / s_rho);
    let gamma1 = chi_rho - chi_t * s_rho / s_t;
    let delta = chi_t / chi_rho;
    let c_p = s_t - s_rho * delta;
    let pr = prad(t);
    EosState {
        rho,
        t,
        p,
        pgas: p - pr,
        prad: pr,
        e: ed.v,
        s: s_d.v,
        chi_rho,
        chi_t,
        s_rho,
        s_t,
        c_v: ed.d[1] / t,
        nabla_ad,
        gamma1,
        delta,
        c_p,
        c_s: (gamma1 * p / rho).sqrt(),
        beta: 1.0 - ARAD * t.powi(4) / (3.0 * p),
        saha: Saha {
            ne: sd.ne.v,
            n_h0: sd.n_h0.v,
            n_hp: sd.n_hp.v,
            n_he0: sd.n_he0.v,
            n_he1: sd.n_he1.v,
            n_he2: sd.n_he2.v,
            n_h: sd.n_h.v,
            n_he: sd.n_he.v,
            x_h: sd.x_h.v,
        },
    }
}

/// Density from pressure and temperature: solves Equation (eq-pressure) of the text,
/// `P_gas(ρ, T) = P − a T⁴/3`, for ρ by safeguarded Newton iteration in ln ρ (bracket and
/// stopping rule |Δ ln P_gas| < 10⁻¹⁴ as in the Julia `rho_PT`). Returns NaN if P ≤ P_rad(T).
pub fn rho_from_p_t(p: f64, t: f64) -> f64 {
    let pg = p - prad(t);
    if pg <= 0.0 {
        return f64::NAN;
    }
    let lpg = pg.ln();
    let mut lo = (pg * 0.5 * MH / (KB * t)).ln();
    let mut hi = (pg * 1.4 * MH / (KB * t)).ln();
    let mut lr = 0.5 * (lo + hi);
    let td = Dual::<1>::new(t, [0.0]);
    for _ in 0..200 {
        let rd = Dual::<1>::new(lr, [1.0]).exp();
        let s = saha_g(rd, td);
        let f = pgas_g(&s, td).ln() - lpg;
        if f.v.abs() < 1e-14 {
            break;
        }
        if f.v > 0.0 {
            hi = lr;
        } else {
            lo = lr;
        }
        let lrn = lr - f.v / f.d[0];
        lr = if lo < lrn && lrn < hi { lrn } else { 0.5 * (lo + hi) };
    }
    lr.exp()
}
