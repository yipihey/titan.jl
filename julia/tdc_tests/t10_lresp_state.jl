# M6: state tests of the network linear response on the material's CK network (needs the loaded material: run through the LRESPTEST hook of popiii_m2c.jl, restart rp_k650;
# `m = LRESP_MAT`, the modules PM = PopIIIMaterial, TX = TitanResidFast in Main). Log: ENV["LRESP_LOG"] (docs/logs/subconv/lresp/t10_lresp_state.log).
# (a-i) limits against an independent equilibrium solve and a frozen finite-difference EOS; (a-iii) AD of Γ₁(ω), δ(ω), c_P(ω) w.r.t. (ρ, T, x_H⁺, x_H₂) vs central FD;
# (a-iv) the layer's interfaces: M6 against M2's equilibrium/frozen blend.
using ForwardDiff, Printf, LinearAlgebra
m = LRESP_MAT; a_ = Int(TF.getscalar(:ngrs)); b_ = Int(TF.getscalar(:ngre))
const LR = TX.LResp
LOGIO = open(get(ENV, "LRESP_LOG", "/tmp/t10_lresp_state.log"), "w")
lp(args...) = (println(args...); println(LOGIO, args...); flush(LOGIO); flush(stdout))
NFAIL = Ref(0)
check(name, ok::Bool, detail = "") = (lp(ok ? "PASS " : "FAIL ", name, "  ", detail); ok || (NFAIL[] += 1); ok)
rel(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300)
lp("t10_lresp_state: restart rp_k650 (3.2 Msun), TitanCore of the branch (not aace2a0); zones numbered from the inner boundary; HM88=1, SLIMP=2 stack, local variables (ρ, T, x_H⁺, x_H₂) at the start-of-step star-light terms")
zstate(k) = (PM.TM.v(:dn)[k], PM.TM.v(:tn)[k], TX.XSN[][1][k], TX.XSN[][2][k], TX.XSN[][3][k])
zs = [z for z in 33:50 if TX.cplz(a_ + z - 1) && !m.eqzone[a_ + z - 1]]
lp("coupled network zones tested: ", zs)
arad = LR.ARAD_LR
function thermo_derivs(Pρ, PT, eρ, eT, ρ, T, P)       # independent re-implementation (real derivatives)
    δ = T * PT / (ρ * Pρ); cP = eT - (eρ - P / ρ^2) * PT / Pρ; α = P / (ρ * Pρ); nad = P * δ / (T * ρ * cP)
    (; gam1 = 1 / (α - δ * nad), del = δ, cp = cP)
end
d5(f, x, h) = (8 * (f(x + h) - f(x - h)) - (f(x + 2h) - f(x - 2h))) / (12h)      # 4th-order central difference
# --- (a-i) ---------------------------------------------------------------------------------------------------------------------------------------------------------------
lp("\n(a-i) limits. Frozen: lresp_frozen of the AD data vs 4th-order FD of the CK EOS at fixed (x_H⁺, x_H₂, x_HI). Equilibrium: base state := the network's own equilibrium S(x) = 0 at the zone's (ρ, T)")
function eq_solve(k, ρ, T, x0, xn0)
    u = [log(x0[1] / xn0), log(x0[2])]
    function R(u)
        x2 = exp(u[2]); x1, xn_ = TX.hlogit(u[1], x2); S = PM.cpl_g(m, k, [ρ + zero(x1), T + zero(x1), x1, x2]; xn = ForwardDiff.value(xn_))
        [S[3], S[4]]
    end
    for it in 1:80
        r0 = R(u); J = ForwardDiff.jacobian(R, u); du = -(J \ r0)
        s = 1.0; while s > 1e-8 && (any(!isfinite, R(u .+ s .* du)) || norm(R(u .+ s .* du) ./ (abs.(J) * abs.(u) .+ 1e-300)) > 2 * norm(r0 ./ (abs.(J) * abs.(u) .+ 1e-300)) + 1e-30); s /= 2; end
        u .+= s .* du
        maximum(abs.(s .* du)) < 1e-14 && break
    end
    x2 = exp(u[2]); x1, xn_ = TX.hlogit(u[1], x2); (ForwardDiff.value(x1), x2, ForwardDiff.value(xn_), maximum(abs.(R(u))))
end
for z in zs[1:min(end, 8)]
    k = a_ + z - 1; ρ, T, x1, x2, xn = zstate(k)
    # frozen
    zf = PM.lr_zone(m, k, ρ, T, (x1, x2, xn), nothing, false)
    pu(r, t) = (y = PM.cpl_g(m, k, [r, t, x1, x2]; xn)[1:2]; (y[1] + arad * t^4 / 3, y[2] + arad * t^4 / r))
    Pρ = d5(r -> pu(r, T)[1], ρ, 1e-3 * ρ); PT = d5(t -> pu(ρ, t)[1], T, 1e-3 * T); eρ = d5(r -> pu(r, T)[2], ρ, 1e-3 * ρ); eT = d5(t -> pu(ρ, t)[2], T, 1e-3 * T)
    ff = thermo_derivs(Pρ, PT, eρ, eT, ρ, T, pu(ρ, T)[1])
    ef = max(rel(zf.fr.gam1, ff.gam1), rel(zf.fr.del, ff.del), rel(zf.fr.cp, ff.cp))
    check(@sprintf("zone %d frozen: Γ₁ %.9g δ %.7g c_P %.5e vs FD, max rel diff %.2e", z, zf.fr.gam1, zf.fr.del, zf.fr.cp, ef), ef < 1e-8)
    # equilibrium
    xe = try eq_solve(k, ρ, T, (x1, x2), xn) catch e; nothing end
    if xe === nothing || !(xe[4] < 1e-6 * abs(PM.cpl_g(m, k, [ρ, T, xe[1], xe[2]]; xn = xe[3])[3] + 1e-300) + 1e-3)
        lp(@sprintf("      zone %d: no network equilibrium found at its (ρ, T) (residual %s): skipped", z, xe === nothing ? "-" : @sprintf("%.2e", xe[4]))); continue
    end
    x1e, x2e, xne, _ = xe
    zeq = PM.lr_zone(m, k, ρ, T, (x1e, x2e, xne), nothing, false)
    function eqpu(r, t)
        x1_, x2_, xn_, res = eq_solve(k, r, t, (x1e, x2e), xne)
        y = PM.cpl_g(m, k, [r, t, x1_, x2_]; xn = xn_)[1:2]; (y[1] + arad * t^4 / 3, y[2] + arad * t^4 / r)
    end
    Pρ = d5(r -> eqpu(r, T)[1], ρ, 1e-3 * ρ); PT = d5(t -> eqpu(ρ, t)[1], T, 1e-3 * T); eρ = d5(r -> eqpu(r, T)[2], ρ, 1e-3 * ρ); eT = d5(t -> eqpu(ρ, t)[2], T, 1e-3 * T)
    fe = thermo_derivs(Pρ, PT, eρ, eT, ρ, T, eqpu(ρ, T)[1])
    l0 = LR.lresp(zeq.d, 0.0); fz = zeq.fr
    ee = max(rel(l0.gam1, fe.gam1), rel(l0.del, fe.del), rel(l0.cp, fe.cp))
    check(@sprintf("zone %d equilibrium base (x_H⁺ %.6g, x_H₂ %.3e): lresp(ω=0) Γ₁ %.9g δ %.7g c_P %.5e vs equilibrium solves, max rel diff %.2e (frozen Γ₁ %.9g; τ_max %.2e s)", z, x1e, x2e, l0.gam1, l0.del, l0.cp, ee, fz.gam1, zeq.tau), ee < 1e-7)
end
# --- (a-iii) -------------------------------------------------------------------------------------------------------------------------------------------------------------
lp("\n(a-iii) AD (outer Duals through lr_data's Hessian) vs central FD of the Float64 chain, Γ₁(ω), δ(ω), c_P(ω) in (ρ, T, x_H⁺, x_H₂); ω = 1/τ_max × (0.1, 1, 10)")
for z in filter(z -> z in (33, 36, 38, 40, 41, 42, 43, 44, 46, 48), zs)
    k = a_ + z - 1; ρ, T, x1, x2, xn = zstate(k)
    zz = PM.lr_zone(m, k, ρ, T, (x1, x2, xn), nothing, false); τ = zz.tau
    for fω in (0.1, 1.0, 10.0)
        ω = fω / τ
        F(v; dual = true) = (z_ = PM.lr_zone(m, k, v[1], v[2], (v[3], v[4], dual ? xn : xn - (v[3] - x1) - 2 * (v[4] - x2)), nothing, dual); r = LR.lresp(z_.d, ω); [r.gam1, r.del, r.cp])
        v0 = [ρ, T, x1, x2]
        Jad = ForwardDiff.jacobian(v -> F(v), v0)
        worst = 0.0; wl = ""; nres = 0; nun = 0; f0 = F(v0; dual = false)
        for l in 1:4
            sx = l == 3 ? min(x1, xn) : v0[l]                       # the scale of the variable (x_H⁺ near 1: the neutral fraction is the scale of its relative change)
            h = 1e-3 * sx; vp = copy(v0); vm = copy(v0); vp[l] += h; vm[l] -= h; vp2 = copy(v0); vp2[l] += 2h; vm2 = copy(v0); vm2[l] -= 2h
            fd = (8 * (F(vp; dual = false) - F(vm; dual = false)) - (F(vp2; dual = false) - F(vm2; dual = false))) / (12h)      # 4th-order stencil
            for i in 1:3
                if max(abs(fd[i]), abs(Jad[i, l])) * sx / abs(f0[i]) < 1e-9; nun += 1; continue; end      # sensitivity below 1e-9: below the rounding of the FD, unresolved
                e = abs(Jad[i, l] - fd[i]) / max(abs(fd[i]), abs(Jad[i, l])); nres += 1; e > worst && (worst = e; wl = "var $l row $i")
            end
        end
        check(@sprintf("zone %d ω τ = %g: d(Γ₁, δ, c_P)/d(ρ, T, x_H⁺, x_H₂) AD vs FD (%d resolved of 12, sensitivity ≥ 1e-9), worst %.2e (%s)", z, fω, nres, worst, wl), worst < 1e-5)
    end
end
# --- (a-iv) --------------------------------------------------------------------------------------------------------------------------------------------------------------
lp("\n(a-iv) interfaces of the layer (zones j−1, j both coupled): M6 Re Γ₁,b, δ_b, c_P,b against M2's blend f = ℓ/(ℓ + t_chem a) of the TABLE equilibrium and the CK-frozen values (conv_thermo, t_chem = 1/max|∂S_i/∂x_i|, as M2 coded it).")
lp("      g = G M/r² with M = 3.2 Msun (approximate), H_P from the zones' p_gas + E_r/3, ω = √(g δ_fr/(8 H_P)); columns: interface j, T, ωτ_max, f (M2), Γ₁,fr, Γ₁,b (M6), Γ₁,b (M2), Γ₁(ω→0) network, Γ₁,eq (table), δ_b M6 / M2, c_P,b M6 / M2")
Gc = 6.6743e-8; MsunG = 1.98847e33; rn_ = PM.TM.v(:rn); pg_ = PM.TM.v(:pgn); er_ = PM.TM.v(:ern)
function m2_blob(k, ρ, T, x1, x2, xn, ω0)
    th = PM.conv_thermo(m, k, ρ, T; x = (x1, x2, xn), tablefr = true)
    Jc = ForwardDiff.jacobian(z -> PM.cpl_g(m, k, z; xn)[3:4], [ρ, T, x1, x2])
    tch = 1 / max(abs(Jc[1, 3]), abs(Jc[2, 4])); f = 1 / (1 + tch * ω0)
    (; Γ = 1 / (f / th.gam1_eq + (1 - f) / th.gam1_fr), δ = f * th.del_eq + (1 - f) * th.del_fr, cp = f * th.cp_eq + (1 - f) * th.cp_fr, f, Γeq = th.gam1_eq, Γfr = th.gam1_fr)
end
nrows = 0; dmax_small = 0.0
for z in zs[2:end]
    (z - 1) in zs || continue
    kA = a_ + z - 2; kB = a_ + z - 1
    A = zstate(kA); B = zstate(kB); r = rn_[kB]; g = Gc * 3.2MsunG / r^2
    zA = PM.lr_zone(m, kA, A[1], A[2], A[3:5], nothing, false); zB = PM.lr_zone(m, kB, B[1], B[2], B[3:5], nothing, false)
    ρi = 0.5 * (A[1] + B[1]); Ti = 0.5 * (A[2] + B[2]); Pi = 0.5 * (pg_[kA] + pg_[kB] + (er_[kA] + er_[kB]) / 3); HP = Pi / (ρi * g)
    ω0 = sqrt(g * 0.5 * (zA.fr.del + zB.fr.del) / (8HP))
    tA = zA.cpl ? LR.lresp(zA.d, ω0) : zA.fr; tB = zB.cpl ? LR.lresp(zB.d, ω0) : zB.fr
    Γb = 0.5 * (tA.gam1 + tB.gam1); δb = 0.5 * (tA.del + tB.del); cpb = 0.5 * (tA.cp + tB.cp)
    a2 = m2_blob(kA, A[1], A[2], A[3], A[4], A[5], ω0); b2 = m2_blob(kB, B[1], B[2], B[3], B[4], B[5], ω0)
    Γ2 = 0.5 * (a2.Γ + b2.Γ); δ2 = 0.5 * (a2.δ + b2.δ); cp2 = 0.5 * (a2.cp + b2.cp); f2 = 0.5 * (a2.f + b2.f)
    Γ0 = 0.5 * (LR.lresp(zA.d, 0.0).gam1 + LR.lresp(zB.d, 0.0).gam1); Γeq = 0.5 * (a2.Γeq + b2.Γeq)
    wτ = ω0 * max(zA.tau, zB.tau)
    lp(@sprintf("      %3d  T %8.1f  ωτ %9.3e  f %.3f  Γ₁fr %.6f  Γ₁b M6 %.6f M2 %.6f  Γ₁(0) %.6f  Γ₁eq(table) %.6f   δ %.5f / %.5f   cP %.4e / %.4e", z, Ti, wτ, f2, 0.5 * (zA.fr.gam1 + zB.fr.gam1), Γb, Γ2, Γ0, Γeq, δb, δ2, cpb, cp2))
    global nrows += 1
    wτ < 1e-2 && (global dmax_small = max(dmax_small, abs(Γ0 - Γeq) / Γeq))
end
lp(@sprintf("      interfaces compared: %d; in the ωτ_max < 1e-2 interfaces |Γ₁(0)_network − Γ₁,eq_table|/Γ₁ ≤ %.2e (the table-vs-network equilibrium difference)", nrows, dmax_small))
# --- cost ----------------------------------------------------------------------------------------------------------------------------------------------------------------
let k = a_ + 40 - 1; ρ, T, x1, x2, xn = zstate(k); v0 = [ρ, T, x1, x2]
    chain(v) = (z_ = PM.lr_zone(m, k, v[1], v[2], (v[3], v[4], xn), nothing, true); r = LR.lresp(z_.d, 1e-5); [r.gam1, r.del, r.cp])
    PM.LRCACHE[][k] = nothing; PM.lr_zone(m, k, ρ, T, (x1, x2, xn), nothing, true); ForwardDiff.jacobian(chain, v0)      # (warm-up: compilation)
    tm = minimum(@elapsed((PM.LRCACHE[][k] = nothing; PM.lr_zone(m, k, ρ, T, (x1, x2, xn), nothing, true))) for _ in 1:5)
    th = minimum(@elapsed(PM.lr_zone(m, k, ρ, T, (x1, x2, xn), nothing, true)) for _ in 1:5)
    tj = minimum(@elapsed(ForwardDiff.jacobian(chain, v0)) for _ in 1:5)
    tc = minimum(@elapsed(ForwardDiff.jacobian(v -> PM.cpl_g(m, k, v; xn)[1:4], v0)) for _ in 1:5)
    lp(@sprintf("\ncost, zone 40, n = 4 local variables [s, min of 5]: Hessian + Jacobian (cache miss) %.2e; cache hit %.2e; Dual chain (4 partials) lr_zone + lresp, cache hit %.2e; for comparison one ForwardDiff Jacobian of cpl_g %.2e", tm, th, tj, tc))
end
lp("\nt10_lresp_state failures: ", NFAIL[]); close(LOGIO)
