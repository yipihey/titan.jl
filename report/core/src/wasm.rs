//! WebAssembly exports used by the laboratory (built with
//! `wasm-pack build core --target web --release --out-dir ../app/pkg`).
//!
//! All model functions take the total mass in M_sun, the black-hole fraction f_BH = M_BH/M, the infall
//! rate Ṁ in M_sun/yr, α_MLT and ε_in; the luminosity is always L = L_Edd(κ_es) of that mass
//! ([`crate::constants::l_edd`]). Results are flat `Float64Array`s whose layout is documented per
//! function. The opacity table is not compiled in: the page fetches `kaptable.bin` and passes its
//! bytes to [`load_table`] once, before calling anything that needs an opacity.
use crate::boundary;
use crate::constants::{l_edd, MSUN, YR};
use crate::envelope::{integrate_inward, Params, Profile, Status};
use crate::eos;
use crate::opacity::KapTable;
use std::cell::RefCell;
use std::rc::Rc;
use wasm_bindgen::prelude::*;

thread_local! {
    static TABLE: RefCell<Option<Rc<KapTable>>> = const { RefCell::new(None) };
}

fn table() -> Result<Rc<KapTable>, JsError> {
    TABLE.with(|t| t.borrow().clone()).ok_or_else(|| JsError::new("opacity table not loaded: call load_table(bytes) first"))
}

fn params(tab: &KapTable, m_msun: f64, f_bh: f64, alpha: f64, eps_in: f64) -> Params<'_> {
    let mut p = Params::quasi_star(tab, m_msun, f_bh);
    p.alpha = alpha;
    p.eps_in = eps_in;
    p
}

/// Parse and install the opacity table (the bytes of `kaptable.bin`; see [`crate::opacity`]).
/// Replaces any previously loaded table. Throws on a malformed file.
#[wasm_bindgen]
pub fn load_table(bytes: &[u8]) -> Result<(), JsError> {
    let t = KapTable::from_bytes(bytes).map_err(|e| JsError::new(&e))?;
    TABLE.with(|c| *c.borrow_mut() = Some(Rc::new(t)));
    Ok(())
}

/// True once [`load_table`] has succeeded.
#[wasm_bindgen]
pub fn table_loaded() -> bool {
    TABLE.with(|c| c.borrow().is_some())
}

/// Equation of state at (ρ in g cm⁻³, T in K) ([`crate::eos::eos_state`]). Returns 14 numbers:
/// `[P, P_gas, e, s, ∇_ad, Γ₁, δ, c_P, c_s, β, n_e, x_H, y_He1, y_He2]` (cgs; x_H = n(H⁺)/n_H,
/// y_He1 = n(He⁺)/n_He, y_He2 = n(He²⁺)/n_He). Needs no table.
#[wasm_bindgen]
pub fn eos_state(rho: f64, t: f64) -> Vec<f64> {
    let s = eos::eos_state(rho, t);
    vec![
        s.p, s.pgas, s.e, s.s, s.nabla_ad, s.gamma1, s.delta, s.c_p, s.c_s, s.beta, s.saha.ne, s.saha.x_h,
        s.saha.y_he1(), s.saha.y_he2(),
    ]
}

/// Density ρ(P, T) (NaN if P ≤ a T⁴/3), [`crate::eos::rho_from_p_t`]. Needs no table.
#[wasm_bindgen]
pub fn rho_from_p_t(p: f64, t: f64) -> f64 {
    eos::rho_from_p_t(p, t)
}

/// Opacities at (ρ, T): `[κ_R, κ_P, in_table]` (cm² g⁻¹; in_table = 1 inside the tabulated range, 0 if the
/// arguments were clamped to its edge).
#[wasm_bindgen]
pub fn kappa(rho: f64, t: f64) -> Result<Vec<f64>, JsError> {
    let tab = table()?;
    Ok(vec![tab.kappa_r(rho, t), tab.kappa_p(rho, t), if tab.in_table(rho, t) { 1.0 } else { 0.0 }])
}

/// κ_R and κ_P along log T at fixed log R (R = ρ/T₆³): `n` equally spaced log T from `log_t_min` to
/// `log_t_max` (inclusive). Returns `2n` numbers interleaved `[κ_R(T₀), κ_P(T₀), κ_R(T₁), κ_P(T₁), …]`.
#[wasm_bindgen]
pub fn kappa_vs_t(log_r: f64, log_t_min: f64, log_t_max: f64, n: u32) -> Result<Vec<f64>, JsError> {
    let tab = table()?;
    let n = n.max(1) as usize;
    let mut out = Vec::with_capacity(2 * n);
    for i in 0..n {
        let lt = if n == 1 { log_t_min } else { log_t_min + (log_t_max - log_t_min) * i as f64 / (n - 1) as f64 };
        let t = 10f64.powf(lt);
        let rho = 10f64.powf(log_r) * (t / 1e6).powi(3);
        out.push(tab.kappa_r(rho, t));
        out.push(tab.kappa_p(rho, t));
    }
    Ok(out)
}

/// Saha ionization fractions at fixed ρ along T: `n` logarithmically spaced temperatures from `t_min` to
/// `t_max` (inclusive). Returns `3n` numbers interleaved `[x_H, y_He1, y_He2]` per temperature
/// ([`crate::eos::saha`]). Needs no table.
#[wasm_bindgen]
pub fn ionization_vs_t(rho: f64, t_min: f64, t_max: f64, n: u32) -> Vec<f64> {
    let n = n.max(1) as usize;
    let (a, b) = (t_min.ln(), t_max.ln());
    let mut out = Vec::with_capacity(3 * n);
    for i in 0..n {
        let t = if n == 1 { t_min } else { (a + (b - a) * i as f64 / (n - 1) as f64).exp() };
        let s = eos::saha(rho, t);
        out.push(s.x_h);
        out.push(s.y_he1());
        out.push(s.y_he2());
    }
    out
}

/// Shooting scan ([`crate::boundary::shoot_sign`]) on the grid [`crate::boundary::scan_grid`]`(log_r_min,
/// log_r_max, dlog_r)` of log10 R (R in cm). Returns `2 n_scan` numbers interleaved `[sign, m(r_in)/M_BH]` per scan
/// point: sign = +1 (envelope mass left at r_in), −1 (mass exhausted or m(r_in) < M_BH), 0 (no boundary state or
/// failed integration); m(r_in)/M_BH is 1 when exhausted and NaN when failed.
#[wasm_bindgen]
pub fn shoot_scan(
    m_msun: f64, f_bh: f64, mdot_msun_yr: f64, alpha: f64, eps_in: f64, log_r_min: f64, log_r_max: f64, dlog_r: f64,
) -> Result<Vec<f64>, JsError> {
    let tab = table()?;
    let p = params(&tab, m_msun, f_bh, alpha, eps_in);
    let mdot = mdot_msun_yr * MSUN / YR;
    let grid = boundary::scan_grid(log_r_min, log_r_max, dlog_r);
    let mut out = Vec::with_capacity(2 * grid.len());
    for lr in grid {
        let s = boundary::shoot_sign(10f64.powf(lr), &p, mdot);
        out.push(s.sign as f64);
        out.push(s.m_rin / p.m_bh);
    }
    Ok(out)
}

/// All steady solutions ([`crate::boundary::solve_all`]: scan of log10 R on `log_r_min:dlog_r:log_r_max`, then 40
/// bisections per sign change; the report uses 14, 16.6, 0.025). Returns 10 numbers per root, roots in increasing R:
/// `[R, T_b, P_b, P_gas,b, ram, τ_in, r_ph, T_obs, Γ_κ,b, L_tot]` (cgs: cm, K, dyn cm⁻², erg s⁻¹; ram = ρ_in v_ff² at
/// R; r_ph = R when τ_in < 2/3; T_obs = [L_tot/(4π σ r_ph²)]^{1/4}). Empty if there is no solution.
#[wasm_bindgen]
pub fn solve_all(
    m_msun: f64, f_bh: f64, mdot_msun_yr: f64, alpha: f64, eps_in: f64, log_r_min: f64, log_r_max: f64, dlog_r: f64,
) -> Result<Vec<f64>, JsError> {
    let tab = table()?;
    let p = params(&tab, m_msun, f_bh, alpha, eps_in);
    let (roots, _) = boundary::solve_all(&p, mdot_msun_yr * MSUN / YR, log_r_min, log_r_max, dlog_r);
    let mut out = Vec::with_capacity(10 * roots.len());
    for r in roots {
        let b = r.bc;
        out.extend_from_slice(&[r.r, b.tb, b.pb, b.pgas_b, b.ram, b.tau_in, b.r_ph, r.t_obs, b.gamma_kappa, b.l_tot]);
    }
    Ok(out)
}

/// Number of fields per profile row of [`envelope_profile`].
pub const PROFILE_FIELDS: usize = 15;

/// Envelope profile for outer radius `r` (cm) ([`crate::boundary::outer_bc`] then
/// [`crate::envelope::integrate_inward`] with `save`). Layout: a 4-number header
/// `[n_points, n_fields (= 15), status, m(r_in)]` (status 1 = reached r_in, 2 = mass exhausted before r_in,
/// 0 = failed or no boundary state; m(r_in) in g), followed by `n_points` rows (outside in, one per accepted
/// step) of `n_fields` numbers:
/// `[r, m, P, T, ρ, ∇, ∇_ad, ∇_rad, v_conv, c_s, κ, Γ, β, H_P, Γ₁]` (cgs; Γ = κL/(4πcGm) the local Eddington
/// factor, β = P_gas/P). For a steady solution pass an R from [`solve_all`].
#[wasm_bindgen]
pub fn envelope_profile(m_msun: f64, f_bh: f64, mdot_msun_yr: f64, alpha: f64, eps_in: f64, r: f64) -> Result<Vec<f64>, JsError> {
    let tab = table()?;
    let p = params(&tab, m_msun, f_bh, alpha, eps_in);
    let (status, m_rin, prof) = match boundary::outer_bc(r, &p, mdot_msun_yr * MSUN / YR) {
        Some(bc) => {
            let it = integrate_inward(r, bc.pb, bc.tb, &p, true);
            (it.status, it.m_rin, it.profile.unwrap_or_default())
        }
        None => (Status::Failed, f64::NAN, Profile::default()),
    };
    let n = prof.len();
    let mut out = Vec::with_capacity(4 + n * PROFILE_FIELDS);
    out.extend_from_slice(&[n as f64, PROFILE_FIELDS as f64, status.code(), m_rin]);
    let cols = prof.fields();
    for i in 0..n {
        for c in cols.iter() {
            out.push(c[i]);
        }
    }
    Ok(out)
}

/// Eddington luminosity for electron scattering `L_Edd = 4π c G M/(0.2(1 + X))` in erg s⁻¹ for M in M_sun
/// (the envelope luminosity of the model).
#[wasm_bindgen]
pub fn luminosity_eddington(m_msun: f64) -> f64 {
    l_edd(m_msun * MSUN)
}
