# M6: stand-alone tests of julia/LResp.jl (the network's linear response of the blob thermodynamics); no solver, no material.
# (the state tests on the material's CK network are julia/tdc_tests/t10_lresp_state.jl, run through the LRESPTEST hook of popiii_m2c.jl)
# usage: julia --project=<env with ForwardDiff> t10_lresp.jl      log: docs/logs/subconv/lresp/t10_lresp.log
using ForwardDiff, LinearAlgebra, Printf
include(joinpath(@__DIR__, "..", "LResp.jl")); using .LResp

const LOGDIR = normpath(joinpath(@__DIR__, "..", "..", "docs", "logs", "subconv", "lresp")); mkpath(LOGDIR)
const LOGIO = open(joinpath(LOGDIR, "t10_lresp.log"), "w")
lp(args...) = (println(args...); println(LOGIO, args...); flush(LOGIO))
const NFAIL = Ref(0)
check(name, ok::Bool, detail = "") = (lp(ok ? "PASS " : "FAIL ", name, "  ", detail); ok || (NFAIL[] += 1); ok)
rel(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300)

# a synthetic zone with physical scales: an ideal gas (R_g = k_B/(μ m_u), μ = 1.2) whose degree x1 of ionisation (energy E per unit x1) relaxes to x_eq(ρ, T) = ½(1 + tanh((T − T0)/ΔT) − ...) with rate 1/τ,
# a second species x2 that is coupled weakly to x1; S_i = dx_i/dt; every derivative is analytic.
const Rg = 1.380649e-16 / (1.2 * 1.66053906660e-24); const Eion = 2.18e-11 / (1.66053906660e-24)
function mk(ρ, T; τ1 = 3.0e4, τ2 = 4.0e5, T0 = 1.0e4, ΔT = 2.0e3)
    x1 = 0.5 * (1 + tanh((T - T0) / ΔT)); sech2 = 1 - tanh((T - T0) / ΔT)^2; x2 = 0.1 * x1
    xeq1_T = 0.5 * sech2 / ΔT; xeq1_ρ = -0.05 * x1 / ρ                  # (x_eq(ρ, T) = x1(T) (1 − 0.05 ln ρ/ρ0): a weak density dependence)
    (; ρ, T, p = Rg * ρ * T * (1 + x1 + x2), u = 1.5 * Rg * T * (1 + x1 + x2) + Eion * x1,
       dp = (Rg * T * (1 + x1 + x2), Rg * ρ * (1 + x1 + x2), Rg * ρ * T, Rg * ρ * T), du = (0.0, 1.5 * Rg * (1 + x1 + x2), 1.5 * Rg * T + Eion, 1.5 * Rg * T),
       dS1 = (xeq1_ρ / τ1, xeq1_T / τ1, -1 / τ1, 0.0), dS2 = (0.1 * xeq1_ρ / τ2, 0.1 * xeq1_T / τ2, 0.1 / τ2 * 0.5, -1 / τ2))
end

# (i) independent complex implementation: Julia's Complex matrix solve
function ref_resp(d, ω)
    J = [d.dS1[3] d.dS1[4]; d.dS2[3] d.dS2[4]]; A = im * ω * I(2) - J
    Xρ = A \ [d.dS1[1], d.dS2[1]]; XT = A \ [d.dS1[2], d.dS2[2]]
    ρ, T = d.ρ, d.T; a = LResp.ARAD_LR; P = d.p + a * T^4 / 3
    Pρ = d.dp[1] + d.dp[3] * Xρ[1] + d.dp[4] * Xρ[2]; PT = d.dp[2] + 4a * T^3 / 3 + d.dp[3] * XT[1] + d.dp[4] * XT[2]
    eρ = d.du[1] - a * T^4 / ρ^2 + d.du[3] * Xρ[1] + d.du[4] * Xρ[2]; eT = d.du[2] + 4a * T^3 / ρ + d.du[3] * XT[1] + d.du[4] * XT[2]
    δ = T * PT / (ρ * Pρ); cP = eT - (eρ - P / ρ^2) * PT / Pρ; α = P / (ρ * Pρ); nad = P * δ / (T * ρ * cP)
    (; gam1 = real(1 / (α - δ * nad)), del = real(δ), cp = real(cP), nad = real(nad), alpha = real(α))
end
d0 = mk(1e-9, 1.1e4)
lp("T10 (a-i) LResp against an independent complex-matrix implementation (synthetic ionising gas, ρ = 1e-9, T = 1.1e4, τ = 3e4 and 4e5 s)")
for ω in (0.0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-1, 10.0)
    r1 = lresp(d0, ω); r2 = ref_resp(d0, ω)
    e = maximum(rel(getfield(r1, f), getfield(r2, f)) for f in (:gam1, :del, :cp, :nad, :alpha))
    check(@sprintf("ω = %g: Γ₁ %.9g δ %.6g cP %.6g, max rel diff %.2e", ω, r1.gam1, r1.del, r1.cp, e), e < 1e-12)
end
# limits: ω → ∞ is the frozen EOS, ω → 0 the equilibrium of the network: dx = −J⁻¹(S_ρ δρ + S_T δT) solved independently
fr = lresp_frozen(d0); big = lresp(d0, 1e12)
check("ω → ∞ equals the frozen composition", maximum(rel(getfield(big, f), getfield(fr, f)) for f in (:gam1, :del, :cp, :nad)) < 1e-9, @sprintf("Γ₁ frozen %.9g, ω = 1e12: %.9g", fr.gam1, big.gam1))
function eq_resp(d)      # the equilibrium derivatives by implicit differentiation done explicitly (a different route: solve for dx/dρ, dx/dT, then substitute)
    J = [d.dS1[3] d.dS1[4]; d.dS2[3] d.dS2[4]]
    dxρ = -(J \ [d.dS1[1], d.dS2[1]]); dxT = -(J \ [d.dS1[2], d.dS2[2]])
    (; dxρ, dxT)
end
let e = eq_resp(d0), r0 = lresp(d0, 0.0), a = LResp.ARAD_LR, ρ = d0.ρ, T = d0.T
    P = d0.p + a * T^4 / 3
    Pρ = d0.dp[1] + d0.dp[3] * e.dxρ[1] + d0.dp[4] * e.dxρ[2]; PT = d0.dp[2] + 4a * T^3 / 3 + d0.dp[3] * e.dxT[1] + d0.dp[4] * e.dxT[2]
    eρ = d0.du[1] - a * T^4 / ρ^2 + d0.du[3] * e.dxρ[1] + d0.du[4] * e.dxρ[2]; eT = d0.du[2] + 4a * T^3 / ρ + d0.du[3] * e.dxT[1] + d0.du[4] * e.dxT[2]
    δ = T * PT / (ρ * Pρ); cP = eT - (eρ - P / ρ^2) * PT / Pρ; α = P / (ρ * Pρ); Γ = 1 / (α - δ * P * δ / (T * ρ * cP))
    check("ω = 0 equals the equilibrium response", maximum((rel(r0.gam1, Γ), rel(r0.del, δ), rel(r0.cp, cP))) < 1e-12, @sprintf("Γ₁(0) = %.9g (frozen %.9g)", r0.gam1, fr.gam1))
end

# (ii) single-species Debye form: one relaxing species with a constant-T EOS coupling (S_T = 0, p_T = e_ρ-terms absent), the in-phase part of the dynamic derivatives
#      Re P_ρ^ω = P_ρ^fr + (P_ρ^eq − P_ρ^fr)/(1 + (ωτ)²)
lp("\nT10 (a-ii) Debye form of the dynamic derivative (one species: dS1 = (S_ρ, 0, −1/τ, 0), no second species, base P_ρ small so that the increment is not lost in rounding)")
let τ = 3.7e4, ρ = 1.0, T = 3000.0
    d = (; ρ, T, p = 10.0, u = 1.0e3, dp = (2.0, 0.5, 3.0, 0.0), du = (-1.0, 5.0, -4.0, 0.0), dS1 = (4.0e-2, 0.0, -1 / τ, 0.0), dS2 = (0.0, 0.0, 0.0, -1.0))
    Pfr = lresp_derivs(d, nothing)[1].re; Peq = lresp_derivs(d, 0.0)[1].re
    emax = 0.0
    for ωτ in (1e-3, 0.1, 0.5, 1.0, 2.0, 10.0, 1e3)
        Pω = lresp_derivs(d, ωτ / τ)[1].re; deb = (Peq - Pfr) / (1 + ωτ^2)
        e = rel(Pω - Pfr, deb); emax = max(emax, e); lp(@sprintf("    ωτ = %-7g  Re P_ρ − P_ρ,fr = %.10e   Debye (P_ρ,eq − P_ρ,fr)/(1+ωτ²) = %.10e   rel %.1e", ωτ, Pω - Pfr, deb, e))
    end
    check("Re P_ρ^ω − P_ρ,fr follows the Debye form (P_ρ,eq − P_ρ,fr)/(1 + (ωτ)²)", emax < 1e-12, @sprintf("max rel %.2e", emax))
    # the same for the energy derivative (an exactly linear function of X), and the imaginary part −ωτ/(1+ωτ²)
    emax = 0.0
    for ωτ in (0.1, 1.0, 10.0)
        r = lresp_derivs(d, ωτ / τ)[1]; Δeq = Peq - Pfr
        emax = max(emax, rel(r.im, -Δeq * ωτ / (1 + ωτ^2)))
    end
    check("Im P_ρ^ω = −(P_ρ,eq − P_ρ,fr) ωτ/(1 + (ωτ)²)", emax < 1e-12, @sprintf("max rel %.2e", emax))
    g0 = lresp(d, 1e-7 / τ).gam1; gi = lresp(d, 1e7 / τ).gam1
    check("Re Γ₁(ω): ωτ = 1e-7 is the equilibrium value, ωτ = 1e7 the frozen one", rel(g0, lresp(d, 0.0).gam1) < 1e-12 && rel(gi, lresp_frozen(d).gam1) < 1e-6, @sprintf("Γ₁: %.8g (eq %.8g) … %.8g (frozen %.8g)", g0, lresp(d, 0.0).gam1, gi, lresp_frozen(d).gam1))
end

# (iii) AD of Γ₁(ω), δ(ω), c_P(ω) with respect to the inputs (ρ, T, p, u and all the derivatives, ω), against central differences; and nested Duals (an outer ForwardDiff over an inner one)
lp("\nT10 (a-iii) derivatives of Γ₁(ω), δ(ω), c_P(ω): ForwardDiff against central FD, and second-order (nested Dual) consistency")
function pack(v)      # v = (ρ, T, p, u, dp1..4, du1..4, dS1_1..4, dS2_1..4, ω) -> (Γ, δ, cP)
    d = (; ρ = v[1], T = v[2], p = v[3], u = v[4], dp = (v[5], v[6], v[7], v[8]), du = (v[9], v[10], v[11], v[12]), dS1 = (v[13], v[14], v[15], v[16]), dS2 = (v[17], v[18], v[19], v[20]))
    r = lresp(d, v[21]); [r.gam1, r.del, r.cp]
end
let d = mk(1e-9, 1.1e4); ω = 3e-5
    v0 = Float64[d.ρ, d.T, d.p, d.u, d.dp..., d.du..., d.dS1..., d.dS2..., ω]
    J = ForwardDiff.jacobian(pack, v0); worst = 0.0
    for j in eachindex(v0)
        h = 1e-3 * abs(v0[j]); vp = copy(v0); vm = copy(v0); vp[j] += h; vm[j] -= h; vp2 = copy(v0); vp2[j] += 2h; vm2 = copy(v0); vm2[j] -= 2h
        fd = (8 * (pack(vp) - pack(vm)) - (pack(vp2) - pack(vm2))) / (12h)
        for i in 1:3
            scale = max(abs(fd[i]), abs(J[i, j]))
            scale > 1e-12 * abs(pack(v0)[i]) / abs(v0[j]) && (worst = max(worst, abs(J[i, j] - fd[i]) / scale))
        end
    end
    check("AD Jacobian of (Γ₁, δ, c_P) in the 21 inputs vs 4th-order central FD (h = 1e-3 relative)", worst < 1e-6, @sprintf("worst relative difference %.2e", worst))
    # nested: the gradient of Γ₁ (inner) differentiated again (outer), against FD of the inner gradient
    g(v) = ForwardDiff.gradient(w -> pack(w)[1], v)
    H = ForwardDiff.jacobian(g, v0); worst = 0.0
    for j in (1, 2, 21)
        h = 1e-5 * abs(v0[j]); vp = copy(v0); vm = copy(v0); vp[j] += h; vm[j] -= h
        fd = (g(vp) - g(vm)) / 2h
        for i in eachindex(v0); s = max(abs(fd[i]), abs(H[i, j])); s > 1e-12 * abs(g(v0)[i]) / abs(v0[j]) && (worst = max(worst, abs(H[i, j] - fd[i]) / s)); end
    end
    check("nested ForwardDiff (Hessian of Γ₁ in ρ, T, ω against FD of the AD gradient)", worst < 1e-6, @sprintf("worst relative difference %.2e", worst))
end
lp("\nT10 (LResp) failures: ", NFAIL[]); close(LOGIO); exit(NFAIL[] == 0 ? 0 : 1)
