//! Outer boundary conditions and the shooting on R (port of `boundary.jl` and of the root search
//! of `scripts/branches.jl`).
//!
//! (A) No infall: Eddington grey photosphere at τ = 2/3 (Equation (eq-bc-phot)).
//! (B) Spherical free-fall infall at rate Ṁ onto a radiative accretion shock at r = R
//! (Equation (eq-bc-infall)), with the infall optical depth computed in LTE at the Lucy
//! temperature (Equation (eq-lucy)). For Ṁ → 0, (B) reduces exactly to (A).
use crate::constants::*;
use crate::envelope::{integrate_inward, Params, Profile, Status};
use crate::eos::{prad, rho_from_p_t};
use crate::opacity::KapTable;
use std::f64::consts::PI;

/// Free-fall speed `v_ff = (2 G M / r)^{1/2}`.
pub fn v_ff(m: f64, r: f64) -> f64 {
    (2.0 * G * m / r).sqrt()
}

/// Density of the free-fall infall `ρ_in(r) = Ṁ / (4π r² v_ff(r))` (∝ r^{−3/2}).
pub fn rho_infall(mdot: f64, m: f64, r: f64) -> f64 {
    mdot / (4.0 * PI * r * r * v_ff(m, r))
}

/// Optical depth of the infall and the radius of the τ = 2/3 surface within it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InfallTau {
    /// τ_in = ∫_R^{10⁴R} κ_R ρ_in dr.
    pub tau: f64,
    /// Radius where τ = 2/3 inside the infall (R if τ_in < 2/3).
    pub r_ph: f64,
}

/// Equation (eq-lucy) of the text: the LTE optical depth of the infall,
/// `τ_in = ∫_R^∞ κ_R(ρ_in, T) ρ_in dr`, with the Lucy (1971) temperature
/// `T⁴(r) = T_L⁴ [W(r) + (3/4) τ(r)]`, dilution `W = ½[1 − (1 − R²/r²)^{1/2}]`,
/// `T_L⁴ = L_tot/(4π σ R²)`; integrated inward from 10⁴ R on 4000 logarithmic radii (midpoint
/// rule, as in the Julia code). Also returns the radius r_ph where τ = 2/3 (R if never reached).
pub fn infall_tau(mdot: f64, m: f64, rr: f64, l_tot: f64, tab: &KapTable) -> InfallTau {
    if mdot <= 0.0 {
        return InfallTau { tau: 0.0, r_ph: rr };
    }
    let n = 4000usize;
    let rmax_fac = 1e4;
    let tl4 = l_tot / (4.0 * PI * SIGMA_SB * rr * rr);
    let (la, lb) = ((rmax_fac * rr).ln(), rr.ln());
    let rs = |i: usize| (la + (lb - la) * (i as f64 / (n - 1) as f64)).exp();
    let mut tau = 0.0;
    let mut r_ph = rr;
    let mut found = false;
    let mut r1 = rs(0);
    for i in 0..n - 1 {
        let r2 = rs(i + 1);
        let rm = (r1 * r2).sqrt();
        let rho = rho_infall(mdot, m, rm);
        let w = 0.5 * (1.0 - (1.0 - (rr / rm).powi(2)).max(0.0).sqrt());
        let t = (tl4 * (w + 0.75 * tau)).powf(0.25);
        let kap = tab.kappa_r(rho, t);
        let dtau = kap * rho * (r1 - r2);
        if !found && tau + dtau >= 2.0 / 3.0 {
            r_ph = r1 - (2.0 / 3.0 - tau) / (kap * rho);
            found = true;
        }
        tau += dtau;
        r1 = r2;
    }
    InfallTau { tau, r_ph: if found { r_ph } else { rr } }
}

/// Outer boundary state, see [`outer_bc`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OuterBc {
    /// Total pressure at the base of the boundary layer P_b.
    pub pb: f64,
    /// Temperature T_b.
    pub tb: f64,
    /// Infall optical depth τ_in (0 without infall).
    pub tau_in: f64,
    /// Photosphere radius in the infall (R if τ_in < 2/3).
    pub r_ph: f64,
    /// L_tot = L + G M Ṁ/R.
    pub l_tot: f64,
    /// Ram pressure ρ_in v_ff² at R.
    pub ram: f64,
    /// T_L = [L_tot/(4π σ R²)]^{1/4}.
    pub teff_l: f64,
    /// Γ_κ = κ_b F/(c g) at the boundary.
    pub gamma_kappa: f64,
    /// Gas pressure P_b − a T_b⁴/3.
    pub pgas_b: f64,
}

/// Equations (eq-bc-phot) and (eq-bc-infall) of the text: (P_b, T_b) at radius R for infall rate
/// `mdot` (g s⁻¹; 0 gives the photosphere).
///
/// (eq-bc-infall): `P_b = ρ_in v_ff² + (2/3) g/κ_b + (F_tot/c)(τ_in + 2/3)`,
/// `T_b⁴ = (3/4) T_L⁴ (τ_in + 4/3)`, `T_L⁴ = L_tot/(4π σ R²)`, `L_tot = L + G M Ṁ/R`,
/// `F_tot = L_tot/(4π R²)`, `g = G M/R²`, κ_b = κ_R(ρ(P_b, T_b), T_b), τ_in from [`infall_tau`]
/// (Equation (eq-lucy)). For Ṁ = 0 this is (eq-bc-phot): `T_b = T_eff`,
/// `P_b = (2/3)(g/κ_b + F/c)` ⇔ `P_gas,b = (2/3)(g/κ_b)(1 − Γ_κ)`, `Γ_κ = κ_b F/(c g)`.
///
/// The implicit equation for P_b is solved by bisection in ln P on
/// [ln(P_rad(T_b)(1 + 10⁻¹⁰)), ln(P_rad(T_b) + 10⁸)] to 10⁻¹²; `None` if not bracketed.
pub fn outer_bc(rr: f64, par: &Params, mdot: f64) -> Option<OuterBc> {
    let l_acc = if mdot > 0.0 { G * par.m * mdot / rr } else { 0.0 };
    let l_tot = par.l + l_acc;
    let inf = infall_tau(mdot, par.m, rr, l_tot, par.tab);
    let tl4 = l_tot / (4.0 * PI * SIGMA_SB * rr * rr);
    let tb = (0.75 * tl4 * (inf.tau + 4.0 / 3.0)).powf(0.25);
    let g = G * par.m / (rr * rr);
    let f = l_tot / (4.0 * PI * rr * rr);
    let ram = if mdot > 0.0 { rho_infall(mdot, par.m, rr) * v_ff(par.m, rr).powi(2) } else { 0.0 };
    let prad_b = prad(tb);
    let fp = |p: f64| {
        let rho = rho_from_p_t(p, tb);
        if rho.is_nan() {
            f64::NEG_INFINITY
        } else {
            p - (ram + (2.0 / 3.0) * g / par.tab.kappa_r(rho, tb) + f / C * (inf.tau + 2.0 / 3.0))
        }
    };
    let mut lo = (prad_b * (1.0 + 1e-10)).ln();
    let mut hi = (prad_b + 1e8).ln();
    if fp(hi.exp()) < 0.0 || fp(lo.exp()) > 0.0 {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if fp(mid.exp()) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo < 1e-12 {
            break;
        }
    }
    let pb = (0.5 * (lo + hi)).exp();
    let kb = par.tab.kappa_r(rho_from_p_t(pb, tb), tb);
    Some(OuterBc {
        pb,
        tb,
        tau_in: inf.tau,
        r_ph: inf.r_ph,
        l_tot,
        ram,
        teff_l: tl4.powf(0.25),
        gamma_kappa: kb * f / (C * g),
        pgas_b: pb - prad_b,
    })
}

/// One shot of the shooting method, see [`shoot_sign`].
#[derive(Clone, Copy, Debug)]
pub struct Shot {
    /// +1 if m(r_in) > M_BH (envelope mass left over), −1 if m(r_in) < M_BH or the mass was
    /// exhausted before r_in, 0 if no boundary state or the integration failed.
    pub sign: i32,
    /// m(r_in) (M_BH if exhausted, NaN if failed).
    pub m_rin: f64,
    /// Integration status (Failed also when there is no boundary state).
    pub status: Status,
    /// The boundary state.
    pub bc: Option<OuterBc>,
}

/// Mismatch of the shooting condition m(ε_in R) = M_BH at trial radius R: boundary from
/// [`outer_bc`] (Equations (eq-bc-phot)/(eq-bc-infall)), then Equation (eq-structure) integrated
/// inward ([`integrate_inward`]). See [`Shot`] for the sign convention.
pub fn shoot_sign(rr: f64, par: &Params, mdot: f64) -> Shot {
    let bc = match outer_bc(rr, par, mdot) {
        Some(b) => b,
        None => return Shot { sign: 0, m_rin: f64::NAN, status: Status::Failed, bc: None },
    };
    let it = integrate_inward(rr, bc.pb, bc.tb, par, false);
    let sign = match it.status {
        Status::Failed => 0,
        Status::MassExhausted => -1,
        Status::ReachedRin => {
            if it.m_rin > par.m_bh {
                1
            } else {
                -1
            }
        }
    };
    Shot { sign, m_rin: it.m_rin, status: it.status, bc: Some(bc) }
}

/// The scan grid of log10 R: `lr_min, lr_min + dlr, …` up to `lr_max` (inclusive, with a 10⁻⁹
/// tolerance, like the Julia range `lr_min:dlr:lr_max`; 14.0:0.025:16.6 has 105 points).
pub fn scan_grid(lr_min: f64, lr_max: f64, dlr: f64) -> Vec<f64> {
    if !(dlr > 0.0) || !(lr_max >= lr_min) {
        return vec![lr_min];
    }
    let n = ((lr_max - lr_min) / dlr + 1e-9).floor() as usize + 1;
    (0..n).map(|i| lr_min + i as f64 * dlr).collect()
}

/// A steady solution found by [`solve_all`] or [`solve_first`].
#[derive(Clone, Debug)]
pub struct Root {
    /// Outer radius R (cm).
    pub r: f64,
    /// Boundary state at R.
    pub bc: OuterBc,
    /// Observed temperature `T_obs = [L_tot/(4π σ r_ph²)]^{1/4}` (= T_L without an infall photosphere).
    pub t_obs: f64,
}

fn root_at(lr: f64, par: &Params, mdot: f64) -> Option<Root> {
    let r = 10f64.powf(lr);
    let bc = outer_bc(r, par, mdot)?;
    let t_obs = (bc.l_tot / (4.0 * PI * SIGMA_SB * bc.r_ph * bc.r_ph)).powf(0.25);
    Some(Root { r, bc, t_obs })
}

/// Bisection of the shooting sign on [a, b] in log10 R (Julia: `s == 0` moves b).
fn bisect(mut a: f64, mut b: f64, sa: i32, par: &Params, mdot: f64, iters: usize, tol: f64) -> f64 {
    for _ in 0..iters {
        let mid = 0.5 * (a + b);
        let s = shoot_sign(10f64.powf(mid), par, mdot).sign;
        if s == 0 {
            b = mid;
        } else if s == sa {
            a = mid;
        } else {
            b = mid;
        }
        if b - a < tol {
            break;
        }
    }
    0.5 * (a + b)
}

/// All steady solutions R with m(ε_in R) = M_BH (the shooting condition with Equations
/// (eq-structure), (eq-bc-phot)/(eq-bc-infall)), found exactly as `scripts/branches.jl`: the sign
/// of [`shoot_sign`] on the log10 R grid [`scan_grid`]`(lr_min, lr_max, dlr)`, then 40 bisections of
/// every adjacent pair with opposite non-zero signs. `mdot` in g s⁻¹. Also returns the scan signs.
pub fn solve_all(par: &Params, mdot: f64, lr_min: f64, lr_max: f64, dlr: f64) -> (Vec<Root>, Vec<i32>) {
    let grid = scan_grid(lr_min, lr_max, dlr);
    let sg: Vec<i32> = grid.iter().map(|&lr| shoot_sign(10f64.powf(lr), par, mdot).sign).collect();
    let mut roots = Vec::new();
    for i in 0..grid.len().saturating_sub(1) {
        if sg[i] == 0 || sg[i + 1] == 0 || sg[i] == sg[i + 1] {
            continue;
        }
        let lr = bisect(grid[i], grid[i + 1], sg[i], par, mdot, 40, 0.0);
        if let Some(rt) = root_at(lr, par, mdot) {
            roots.push(rt);
        }
    }
    (roots, sg)
}

/// A single solution with its profile, see [`solve_first`].
#[derive(Clone, Debug)]
pub struct Solution {
    /// The root.
    pub root: Root,
    /// Envelope profile at the root (from [`integrate_inward`] with `save`).
    pub profile: Profile,
    /// Status of that integration.
    pub status: Status,
    /// Scan signs.
    pub signs: Vec<i32>,
}

/// The first steady solution (smallest R) on a scan of `nscan` equally spaced log10 R in
/// [lr_min, lr_max], bisected up to 60 times to |Δ log R| < 10⁻⁷, with its envelope profile — the
/// Julia `solve_R` (defaults there: 14.5, 17.5, 31). `mdot` in g s⁻¹.
pub fn solve_first(par: &Params, mdot: f64, lr_min: f64, lr_max: f64, nscan: usize) -> Option<Solution> {
    let nscan = nscan.max(2);
    let grid: Vec<f64> =
        (0..nscan).map(|i| lr_min + (lr_max - lr_min) * i as f64 / (nscan - 1) as f64).collect();
    let sg: Vec<i32> = grid.iter().map(|&lr| shoot_sign(10f64.powf(lr), par, mdot).sign).collect();
    let idx = (0..nscan - 1).find(|&i| sg[i] != 0 && sg[i + 1] != 0 && sg[i] != sg[i + 1])?;
    let lr = bisect(grid[idx], grid[idx + 1], sg[idx], par, mdot, 60, 1e-7);
    let root = root_at(lr, par, mdot)?;
    let it = integrate_inward(root.r, root.bc.pb, root.bc.tb, par, true);
    Some(Solution { root, profile: it.profile.unwrap_or_default(), status: it.status, signs: sg })
}
