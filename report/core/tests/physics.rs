//! Physical limits and sanity checks of the model (independent of the Julia reference values).
//! Run with `cargo test --release -- --nocapture` to see the measured numbers.
mod common;
use common::*;
use qs_core::boundary::{outer_bc, solve_all};
use qs_core::constants::*;
use qs_core::envelope::{integrate_inward, Params, Status};
use qs_core::eos::{eos_state, prad, rho_from_p_t};

/// Fully ionized gas plus radiation (eq-thermo): Γ₁(β) = (32 − 24β − 3β²)/(24 − 21β) is reproduced to 1e-5 relative
/// for β = 0.5, 0.1, 0.01, both at T = 4×10⁵ K (where Saha leaves the smallest neutral/He⁺ fractions at these β,
/// ~1e-6) and at T = 10⁷ K (the check of the Julia code; there residual He⁺ ≈ 1e-4 at β = 0.5 gives the largest
/// deviation, 4e-6). ρ is found from P = P_rad/(1 − β) with `rho_from_p_t`.
#[test]
fn gamma1_fully_ionized_radiation() {
    let mut w = Worst::new("Gamma1(beta)", 1e-5);
    for t in [4e5, 1e7] {
        for beta in [0.5, 0.1, 0.01] {
            let rho = rho_from_p_t(prad(t) / (1.0 - beta), t);
            let st = eos_state(rho, t);
            assert!((st.beta - beta).abs() < 1e-9 * beta.max(1e-3), "beta {} vs {}", st.beta, beta);
            let want = (32.0 - 24.0 * beta - 3.0 * beta * beta) / (24.0 - 21.0 * beta);
            println!("  T={t:e} beta={beta}: Gamma1={:.9} analytic={want:.9} rel {:.1e}", st.gamma1, rel(st.gamma1, want));
            w.check(st.gamma1, want, &format!("T={t:e} beta={beta}"));
        }
    }
}

/// ∇_ad limits (eq-thermo): → 2/5 for a cold, neutral, gas-dominated monatomic gas (T = 3000 K, ρ = 10⁻⁶ g cm⁻³,
/// β ≈ 1; the EOS has no molecules) and → 1/4 when radiation dominates (T = 10⁷ K, ρ = 10⁻¹⁰, β ≈ 6e-9), to 1e-3
/// absolute; Γ₁ → 5/3 and 4/3 likewise.
#[test]
fn nabla_ad_limits() {
    let gas = eos_state(1e-6, 3000.0);
    let rad = eos_state(1e-10, 1e7);
    println!("  neutral gas: nabla_ad={:.7} Gamma1={:.7} beta={:.7}", gas.nabla_ad, gas.gamma1, gas.beta);
    println!("  radiation:   nabla_ad={:.7} Gamma1={:.7} beta={:.2e}", rad.nabla_ad, rad.gamma1, rad.beta);
    assert!(gas.beta > 0.999 && gas.saha.x_h < 1e-3);
    assert!((gas.nabla_ad - 0.4).abs() < 1e-3);
    assert!((gas.gamma1 - 5.0 / 3.0).abs() < 1e-3);
    assert!(rad.beta < 1e-6);
    assert!((rad.nabla_ad - 0.25).abs() < 1e-3);
    assert!((rad.gamma1 - 4.0 / 3.0).abs() < 1e-3);
}

/// High-temperature limit of the Rosseland mean (eq-rosseland): at log T = 7.9 and low density (log R = −10, −8)
/// the table gives pure Thomson scattering on the free electrons of fully ionized gas,
/// κ_R = σ_T (X/m_H + 2Y/m_He) = 0.35015 cm² g⁻¹, to 1e-4 relative; this is the documented 0.350 (within 1%), below
/// 0.2(1 + X) = 0.352 because m_H rather than m_p enters n_e.
#[test]
fn kappa_high_t_limit() {
    let tab = table();
    let thomson = SIGMA_T * (X / MH + 2.0 * Y / MHE);
    let t = 10f64.powf(7.9);
    for lr in [-10.0, -8.0] {
        let rho = 10f64.powf(lr) * (t / 1e6).powi(3);
        let k = tab.kappa_r(rho, t);
        println!("  log T=7.9 log R={lr}: kappa_R={k:.6}  sigma_T n_e/rho={thomson:.6}  rel {:.1e}", rel(k, thomson));
        assert!(tab.in_table(rho, t));
        assert!(rel(k, thomson) < 1e-4);
        assert!(rel(k, 0.350) < 0.01);
    }
}

fn reference_envelope() -> (f64, qs_core::envelope::Profile, Params<'static>) {
    let p = Params::quasi_star(table(), 1e6, 0.1);
    let (roots, _) = solve_all(&p, 0.0, 15.8, 16.0, 0.025);
    assert_eq!(roots.len(), 1);
    let r = roots[0].r;
    let it = integrate_inward(r, roots[0].bc.pb, roots[0].bc.tb, &p, true);
    assert_eq!(it.status, Status::ReachedRin);
    (r, it.profile.unwrap(), p)
}

/// In the f_BH = 0.1, Ṁ = 0 steady envelope (M = 10⁶ M_sun, L = L_Edd) β = P_gas/P lies in [0.0055, 0.0080] for
/// 0.05 R ≤ r ≤ 0.9 R (measured: 0.00594–0.00720, decreasing inward); the Eddington standard-model quartic
/// 1 − β = 0.003 μ⁴ β⁴ (M/M_sun)² gives β = 0.00725 for μ = 0.588 (fully ionized X = 0.76), within 1% of the
/// envelope's largest bulk value and 22% above its smallest.
#[test]
fn eddington_quartic_beta() {
    let (r, pf, _) = reference_envelope();
    let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
    for i in 0..pf.len() {
        let x = pf.r[i] / r;
        if (0.05..=0.9).contains(&x) {
            lo = lo.min(pf.beta[i]);
            hi = hi.max(pf.beta[i]);
        }
    }
    println!("  beta in 0.05-0.9 R: {lo:.5} .. {hi:.5}");
    assert!(lo >= 0.0055 && hi <= 0.0080, "beta range {lo}..{hi}");
    let mu: f64 = 1.0 / (2.0 * X + 0.75 * Y);
    let f = |b: f64| 1.0 - b - 0.003 * mu.powi(4) * b.powi(4) * 1e12;
    let (mut a, mut b) = (0.0, 1.0);
    for _ in 0..100 {
        let c = 0.5 * (a + b);
        if f(c) > 0.0 {
            a = c
        } else {
            b = c
        }
    }
    println!("  Eddington quartic (mu={mu:.4}): beta = {a:.5}");
    assert!((mu - 0.588).abs() < 1e-3);
    assert!((a - 0.00725).abs() < 5e-5);
    assert!(a > lo && rel(a, hi) < 0.01);
}

/// Sanity of a steady profile (eq-structure): going inward, r and m decrease and P and T increase strictly,
/// m stays above M_BH, the luminosity used in ∇_rad is
/// the constant L to 1e-12 (Γ m = κ L/(4π c G) at every point), and every accepted step satisfies hydrostatic
/// equilibrium dP/dr = −G m ρ/r² in the form Δln P = ∫ −G m ρ/(r P) dln r to 3e-2 relative (measured 2.1e-2; this
/// is the error of the trapezoid rule over the ~100 large Dormand–Prince steps, so it is only a coarse check).
const TOL_HSE: f64 = 3e-2;
#[test]
fn luminosity_constant_and_hydrostatic() {
    let (_, pf, p) = reference_envelope();
    let n = pf.len();
    assert!(n > 50);
    for i in 1..n {
        assert!(pf.r[i] < pf.r[i - 1] && pf.m[i] < pf.m[i - 1], "r, m monotone at {i}");
        assert!(pf.p[i] > pf.p[i - 1] && pf.t[i] > pf.t[i - 1], "P, T monotone at {i}");
        assert!(pf.m[i] > p.m_bh);
    }
    let mut wl = Worst::new("L const", 1e-12);
    for i in 0..n {
        let l = pf.gamma[i] * 4.0 * std::f64::consts::PI * C * G * pf.m[i] / pf.kappa[i];
        wl.check(l, p.l, &format!("point {i}"));
    }
    // hydrostatic equilibrium per accepted step: Δln P = ∫ dln P/dln r dln r with dln P/dln r = −G m ρ/(r P),
    // trapezoid rule in ln r between consecutive saved points
    let g = |i: usize| -G * pf.m[i] * pf.rho[i] / (pf.r[i] * pf.p[i]);
    let mut wh = Worst::new("hydrostatic", TOL_HSE);
    for i in 0..n - 1 {
        let dlnp = pf.p[i + 1].ln() - pf.p[i].ln();
        let trap = 0.5 * (g(i) + g(i + 1)) * (pf.r[i + 1].ln() - pf.r[i].ln());
        wh.check(trap, dlnp, &format!("step {i}"));
    }
    // the photospheric boundary is consistent: P_b at the root equals the profile's first point
    let r0 = pf.r[0];
    let bc = outer_bc(r0, &p, 0.0).unwrap();
    assert_eq!(bc.pb, pf.p[0]);
}
