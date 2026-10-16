//! Port validation: every function of the Rust core reproduces the Julia implementation
//! (`QuasiStarInfall.jl`, the specification) on the values exported to `data/reference.json`.
//! Run with `cargo test --release -- --nocapture` to see the measured deviations.
mod common;
use common::*;
use qs_core::boundary::{outer_bc, shoot_sign, solve_all, solve_first};
use qs_core::constants::{self, MSUN, YR};
use qs_core::envelope::{integrate_inward, mlt, Params, Status};
use qs_core::eos::{eos_state, rho_from_p_t, saha};

fn params(f_bh: f64) -> Params<'static> {
    Params::quasi_star(table(), 1e6, f_bh)
}

/// The constants (and the derived radiation constant, Eddington luminosity, defaults) are those of the Julia code.
#[test]
fn constants_match_reference() {
    let c = &reference()["constants"];
    let mut w = Worst::new("constants", 1e-15);
    for (k, v) in [
        ("G", constants::G), ("c", constants::C), ("kB", constants::KB), ("h", constants::H), ("me", constants::ME),
        ("mH", constants::MH), ("mHe", constants::MHE), ("sigmaSB", constants::SIGMA_SB), ("arad", constants::ARAD),
        ("Msun", MSUN), ("yr", YR), ("chiH", constants::CHI_H), ("chiHe1", constants::CHI_HE1),
        ("chiHe2", constants::CHI_HE2), ("X", constants::X), ("Y", constants::Y),
    ] {
        w.check(v, f(c, k), k);
    }
    let p = params(0.1);
    w.check(p.m, f(c, "M"), "M");
    w.check(p.l, f(c, "L"), "L");
    w.check(p.alpha, f(c, "alpha"), "alpha");
    w.check(p.eps_in, f(c, "eps_in"), "eps_in");
}

/// Saha ionization equilibrium (eq-saha): n_e, x_H, n(He⁺), n(He²⁺) agree with Julia to 1e-9 relative
/// on nine (ρ, T) points from neutral to fully ionized.
#[test]
fn saha_matches_reference() {
    let mut wne = Worst::new("n_e", 1e-9);
    let mut wx = Worst::new("x_H", 1e-9);
    let mut w1 = Worst::new("n_He+", 1e-9);
    let mut w2 = Worst::new("n_He++", 1e-9);
    for v in reference()["saha"].as_array().unwrap() {
        let (rho, t) = (f(v, "rho"), f(v, "T"));
        let s = saha(rho, t);
        let ctx = format!("rho={rho:e} T={t}");
        wne.check(s.ne, f(v, "ne"), &ctx);
        wx.check(s.x_h, f(v, "xH"), &ctx);
        w1.check(s.n_he1, f(v, "nHe1"), &ctx);
        w2.check(s.n_he2, f(v, "nHe2"), &ctx);
    }
}

/// Equation of state (eq-pressure, eq-entropy, eq-thermo): P, e, s agree with Julia to 1e-10 relative and
/// ∇_ad, Γ₁, δ, c_P (exact derivatives) to 1e-6 relative.
#[test]
fn eos_matches_reference() {
    let mut ws: Vec<Worst> = [("P", 1e-10), ("e", 1e-10), ("s", 1e-10), ("nabla_ad", 1e-6), ("Gamma1", 1e-6),
        ("delta", 1e-6), ("cP", 1e-6), ("cs", 1e-6), ("beta", 1e-6)]
        .iter().map(|&(n, t)| Worst::new(n, t)).collect();
    for v in reference()["eos"].as_array().unwrap() {
        let (rho, t) = (f(v, "rho"), f(v, "T"));
        let st = eos_state(rho, t);
        let ctx = format!("rho={rho:e} T={t}");
        let got = [st.p, st.e, st.s, st.nabla_ad, st.gamma1, st.delta, st.c_p, st.c_s, st.beta];
        for (w, g) in ws.iter_mut().zip(got) {
            let want = f(v, w.name);
            w.check(g, want, &ctx);
        }
    }
}

/// ρ(P, T), the inversion of Equation (eq-pressure), agrees with Julia to 1e-9 relative (and is NaN where P ≤ P_rad).
#[test]
fn rho_pt_matches_reference() {
    let mut w = Worst::new("rho(P,T)", 1e-9);
    for v in reference()["rho_PT"].as_array().unwrap() {
        let (p, t) = (f(v, "P"), f(v, "T"));
        let r = rho_from_p_t(p, t);
        let want = f(v, "rho");
        if want.is_nan() {
            assert!(r.is_nan(), "P={p} T={t}: expected NaN, got {r}");
        } else {
            w.check(r, want, &format!("P={p:e} T={t}"));
        }
    }
}

/// Opacity table interpolation (eq-rosseland): κ_R and κ_P from the shipped (float32) table agree with the
/// Julia bicubic interpolation of the float64 table to 1e-5 relative, including clamped points outside it.
#[test]
fn kappa_matches_reference() {
    let tab = table();
    let mut wr = Worst::new("kappa_R", 1e-5);
    let mut wp = Worst::new("kappa_P", 1e-5);
    for v in reference()["kap"].as_array().unwrap() {
        let (rho, t) = (f(v, "rho"), f(v, "T"));
        let ctx = format!("rho={rho:e} T={t}");
        wr.check(tab.kappa_r(rho, t), f(v, "kR"), &ctx);
        wp.check(tab.kappa_p(rho, t), f(v, "kP"), &ctx);
    }
}

/// Mixing-length theory (eq-mlt-cubic, eq-mlt-grad) at eight points of the f_BH = 0.1, Ṁ = 0 steady
/// envelope: ∇ and v_c agree with Julia to 1e-6 relative (∇_rad, ∇_ad, ρ, κ, Γ to 1e-6 as well).
#[test]
fn mlt_matches_reference() {
    let p = params(0.1);
    let mut ws: Vec<Worst> = ["nabla", "v", "nabla_rad", "nabla_ad", "rho", "kappa", "Gamma"]
        .iter().map(|&n| Worst::new(n, 1e-6)).collect();
    for v in reference()["mlt"].as_array().unwrap() {
        let q = mlt(f(v, "P"), f(v, "T"), f(v, "r"), f(v, "m"), p.l, &p).expect("EOS state");
        let ctx = format!("P={:e} T={}", f(v, "P"), f(v, "T"));
        let got = [q.nabla, q.v, q.nabla_rad, q.nabla_ad, q.rho, q.kappa, q.gamma];
        for (w, g) in ws.iter_mut().zip(got) {
            let want = f(v, w.name);
            if want == 0.0 {
                assert_eq!(g, 0.0, "{} at {}", w.name, ctx);
            } else {
                w.check(g, want, &ctx);
            }
        }
    }
}

/// Outer boundary (eq-bc-phot, eq-bc-infall, eq-lucy): P_b, T_b, τ_in, r_ph and Γ_κ agree with Julia to
/// 1e-6 relative, without infall and for Ṁ = 0.3, 1, 3 M_sun/yr.
#[test]
fn outer_bc_matches_reference() {
    let mut ws: Vec<Worst> = ["Pb", "Tb", "tau_in", "r_ph", "Gamma_kappa"].iter().map(|&n| Worst::new(n, 1e-6)).collect();
    for v in reference()["outer_bc"].as_array().unwrap() {
        let p = params(f(v, "fBH"));
        let md = f(v, "Mdot_Msun_yr");
        let bc = outer_bc(f(v, "R"), &p, md * MSUN / YR).expect("boundary state");
        let ctx = format!("fBH={} Mdot={} R={:e}", f(v, "fBH"), md, f(v, "R"));
        let got = [bc.pb, bc.tb, bc.tau_in, bc.r_ph, bc.gamma_kappa];
        for (w, g) in ws.iter_mut().zip(got) {
            let want = f(v, w.name);
            if want == 0.0 {
                assert_eq!(g, 0.0, "{} at {}", w.name, ctx);
            } else {
                w.check(g, want, &ctx);
            }
        }
    }
}

/// Inward integration (eq-structure) with Dormand–Prince 5(4): the integration ends the same way as in Julia
/// (reached r_in / mass exhausted) and m(r_in) agrees to 1e-4 relative.
#[test]
fn shoot_matches_reference() {
    let mut w = Worst::new("m(r_in)", 1e-4);
    for v in reference()["shoot"].as_array().unwrap() {
        let p = params(f(v, "fBH"));
        let md = f(v, "Mdot_Msun_yr");
        let r = f(v, "R");
        let bc = outer_bc(r, &p, md * MSUN / YR).unwrap();
        let it = integrate_inward(r, bc.pb, bc.tb, &p, false);
        let want = match v["status"].as_str().unwrap() {
            "reached_rin" => Status::ReachedRin,
            "mass_exhausted" => Status::MassExhausted,
            _ => Status::Failed,
        };
        let ctx = format!("fBH={} Mdot={} R={:e}", f(v, "fBH"), md, r);
        assert_eq!(it.status, want, "status at {ctx}");
        w.check(it.m_rin, f(v, "m_rin"), &ctx);
        let s = shoot_sign(r, &p, md * MSUN / YR);
        assert_eq!(s.status, want);
    }
}

/// Run `solve_all` with the scan of `scripts/branches.jl` for every (f_BH, Ṁ) of branches.jl on two threads.
fn all_branches() -> &'static Vec<(f64, f64, Vec<qs_core::boundary::Root>, Vec<i32>)> {
    use std::sync::OnceLock;
    static B: OnceLock<Vec<(f64, f64, Vec<qs_core::boundary::Root>, Vec<i32>)>> = OnceLock::new();
    B.get_or_init(|| {
        let t0 = std::time::Instant::now();
        let mdots = [0.0, 0.3, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 30.0, 100.0];
        let jobs: Vec<(f64, f64)> = [0.1, 0.01].iter().flat_map(|&fb| mdots.iter().map(move |&m| (fb, m))).collect();
        // two worker threads (each solve is single-threaded, as in the browser)
        let mut out: Vec<_> = std::thread::scope(|sc| {
            let hs: Vec<_> = (0..2)
                .map(|k| {
                    let ch: Vec<(f64, f64)> = jobs.iter().copied().skip(k).step_by(2).collect();
                    sc.spawn(move || {
                        ch.into_iter()
                            .map(|(fbh, md)| {
                                let (roots, sg) = solve_all(&params(fbh), md * MSUN / YR, 14.0, 16.6, 0.025);
                                (fbh, md, roots, sg)
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
        });
        out.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
        println!("  solve_all for 22 (f_BH, Mdot): {:.1} s wall on 2 threads", t0.elapsed().as_secs_f64());
        out
    })
}

/// The central test of the port: `solve_all` with the scan of `scripts/branches.jl` (log R = 14:0.025:16.6,
/// 40 bisections per sign change) reproduces the steady solutions of the Julia code (`branches_current` in
/// reference.json: the same procedure run with the current Julia source) for all 22 (f_BH, Ṁ) of branches.jl
/// (Ṁ = 0 … 100 M_sun/yr): the same number of roots everywhere (none for Ṁ ≥ 10 at f_BH = 0.1 and Ṁ ≥ 30 at
/// f_BH = 0.01), every R to 2e-3 relative, T_b and T_obs to 1e-2 relative.
#[test]
fn branches_match_reference() {
    let refs = reference()["branches_current"].as_array().unwrap();
    let mut wr = Worst::new("R", 2e-3);
    let mut wt = Worst::new("T_b", 1e-2);
    let mut wo = Worst::new("T_obs", 1e-2);
    let mut nsign = 0;
    for (fbh, md, roots, sg) in all_branches() {
        let v = refs.iter().find(|v| f(v, "fBH") == *fbh && f(v, "Mdot_Msun_yr") == *md).expect("reference case");
        println!("  fBH={fbh:<4} Mdot={md:<5}: R = {:?}", roots.iter().map(|r| format!("{:.6e}", r.r)).collect::<Vec<_>>());
        assert_eq!(roots.len() as f64, f(v, "nroots"), "number of roots at fBH={fbh} Mdot={md}");
        let jsg: Vec<i64> = v["signs"].as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect();
        nsign += sg.iter().zip(&jsg).filter(|(a, b)| **a as i64 != **b).count();
        for (k, rt) in roots.iter().enumerate() {
            let ctx = format!("fBH={fbh} Mdot={md} root {k}");
            wr.check(rt.r, v["R"][k].as_f64().unwrap(), &ctx);
            wt.check(rt.bc.tb, v["Tb"][k].as_f64().unwrap(), &ctx);
            wo.check(rt.t_obs, v["T_obs"][k].as_f64().unwrap(), &ctx);
        }
    }
    println!("  scan signs differing from Julia: {nsign} of {}", 22 * 105);
}

/// With the same (Float32-rounded) opacity data as the Rust core, the Julia code's steady solutions
/// (`branches_current_f32table`) are reproduced with identical scan signs at all 22 × 105 scan points, the same
/// number of roots, and every R, T_b, T_obs to 1e-6 relative (measured: ≤ 4e-7; the residual is floating-point
/// roundoff — the Dormand–Prince error estimate y5 − y4 cancels ~7 digits, so step sizes carry ~1e-9 noise that the
/// 40 bisections resolve near a root). The 2.7e-5 differences in `branches_match_reference` therefore come from
/// storing the table in Float32, not from the port.
#[test]
fn branches_match_julia_f32_table() {
    let refs = reference()["branches_current_f32table"].as_array().unwrap();
    let mut wr = Worst::new("R f32", 1e-6);
    let mut wt = Worst::new("T_b f32", 1e-6);
    let mut wo = Worst::new("T_obs f32", 1e-6);
    for (fbh, md, roots, sg) in all_branches() {
        let v = refs.iter().find(|v| f(v, "fBH") == *fbh && f(v, "Mdot_Msun_yr") == *md).expect("reference case");
        assert_eq!(roots.len() as f64, f(v, "nroots"), "number of roots at fBH={fbh} Mdot={md}");
        let jsg: Vec<i32> = v["signs"].as_array().unwrap().iter().map(|x| x.as_i64().unwrap() as i32).collect();
        assert_eq!(sg, &jsg, "scan signs at fBH={fbh} Mdot={md}");
        for (k, rt) in roots.iter().enumerate() {
            let ctx = format!("fBH={fbh} Mdot={md} root {k}");
            wr.check(rt.r, v["R"][k].as_f64().unwrap(), &ctx);
            wt.check(rt.bc.tb, v["Tb"][k].as_f64().unwrap(), &ctx);
            wo.check(rt.t_obs, v["T_obs"][k].as_f64().unwrap(), &ctx);
        }
    }
}

/// The table printed in the book (`steady/output/branches.csv`, Table tbl-branches), regenerated with the current
/// Julia code on 2026-10-16 (the 15 October version came from an older code and differed by up to 4.6e-3 in R).
/// Pins the same number of roots for every (f_BH, Ṁ), the compact-branch R to 2e-3, every R to 5e-3, T_b and T_obs
/// to 1e-2 — loose bounds kept from when the table was stale; against the regenerated table the agreement is ~3e-5.
#[test]
fn published_branch_table() {
    let refs = reference()["branches"].as_array().unwrap();
    let mut wc = Worst::new("R compact", 2e-3);
    let mut wr = Worst::new("R all", 5e-3);
    let mut wt = Worst::new("T_b", 1e-2);
    let mut wo = Worst::new("T_obs", 1e-2);
    for (fbh, md, roots, _) in all_branches() {
        let want: Vec<_> = refs.iter().filter(|v| f(v, "fBH") == *fbh && f(v, "Mdot_Msun_yr") == *md).collect();
        assert_eq!(roots.len(), want.len(), "number of roots at fBH={fbh} Mdot={md}");
        for (k, (rt, v)) in roots.iter().zip(&want).enumerate() {
            let ctx = format!("fBH={fbh} Mdot={md} R={:e}", f(v, "R"));
            if roots.len() == 2 && k == 0 {
                wc.check(rt.r, f(v, "R"), &ctx);
            }
            wr.check(rt.r, f(v, "R"), &ctx);
            wt.check(rt.bc.tb, f(v, "Tb"), &ctx);
            wo.check(rt.t_obs, f(v, "T_obs"), &ctx);
        }
    }
}

/// The single-solution search (`solve_first`, the Julia `solve_R` with scan 15:0.05:16.8, 60 bisections) gives
/// the radius of the exported reference profile (f_BH = 0.1, Ṁ = 0) to 1e-4 relative. Integrating inward from that
/// radius with the same (Float32-rounded) opacity data as Julia's `profile_f01_md0_f32table` takes the same accepted
/// steps and reproduces r, m, P, T, ρ, ∇, v_c, Γ, β at every step, and m(r_in), to 1e-6 relative (measured ≤ 3.4e-8,
/// roundoff in the step-size control).
///
/// (Against the Float64-table Julia profile the step sequences drift apart after ~25 steps and r(P), m(P) differ by
/// up to ~1e-3 in the deep interior while T(P) agrees to ~1e-6: at the root radius the shooting is ill-conditioned —
/// m(r_in)/M_BH jumps from 0.999 to "exhausted" within Δlog R ≈ 2e-4 — so table rounding of ~1e-7 is amplified.)
#[test]
fn profile_matches_reference() {
    let p = params(0.1);
    let pr = &reference()["profile_f01_md0"];
    let s = solve_first(&p, 0.0, 15.0, 16.8, 37).expect("a solution");
    let mut wr = Worst::new("R (solve_R)", 1e-4);
    wr.check(s.root.r, f(pr, "R"), "f_BH=0.1 Mdot=0");
    assert_eq!(s.status, Status::ReachedRin);

    let j = &reference()["profile_f01_md0_f32table"];
    let r0 = f(j, "R");
    let bc = outer_bc(r0, &p, 0.0).unwrap();
    let it = integrate_inward(r0, bc.pb, bc.tb, &p, true);
    assert_eq!(it.status, Status::ReachedRin);
    let mut wm = Worst::new("m(r_in) f32", 1e-6);
    wm.check(it.m_rin, f(j, "m_rin"), "f_BH=0.1 Mdot=0");
    let pf = it.profile.unwrap();
    assert_eq!(pf.len(), j["r"].as_array().unwrap().len(), "number of accepted steps");
    let names = ["r", "m", "P", "T", "rho", "nabla", "v", "Gamma", "beta"];
    let ours = [&pf.r, &pf.m, &pf.p, &pf.t, &pf.rho, &pf.nabla, &pf.v_conv, &pf.gamma, &pf.beta];
    for (name, col) in names.iter().zip(ours) {
        let mut w = Worst::new(name, 1e-6);
        for (k, g) in col.iter().enumerate() {
            let want = j[*name][k].as_f64().unwrap();
            if want == 0.0 {
                assert!(g.abs() < 1e-300, "{name} at step {k}: {g} vs 0");
            } else {
                w.check(*g, want, &format!("step {k}"));
            }
        }
    }
}
