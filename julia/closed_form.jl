# Closed-form ionization of the infall (the analytic limit of analytic_infall.jl), checked against the full non-LTE
# solution. Assumptions: steady free fall; the Lyman continuum in detailed balance (case B); Lyα and Lyβ optically
# thick and in radiative balance with the star's diluted far-UV continuum, so n₂/n₁ = (g₂/g₁) J(Lyα)/(2hν³/c²) and
# n₃/n₁ = (g₃/g₁) J(Lyβ)/(2hν³/c²) (radiative "Boltzmann" populations); photoionization of n = 2 and n = 3 by the star's
# diluted, Thomson-trapped Balmer and Paschen continua, Γ_i = 4π ∫ J σ_i/(hν) dν; case-B recombination α_B(T).
# Balance per volume:  n_e n_p α_B = n₁ (n₂/n₁ Γ₂ + n₃/n₁ Γ₃),  n_e = n_p = x n_H,  n₁ ≈ (1 − x) n_H:
#     x² / (1 − x) = Q(r) / (n_H α_B),   Q = (n₂/n₁) Γ₂ + (n₃/n₁) Γ₃ .
# Star: the emergent spectrum of the TITAN star without infall (star_t0), J = (W + 3τ*R²/4r²) I_ν with τ* from x itself
# (iterated). Temperatures: those of the W30 non-LTE solution (for the comparison) — they enter only through α_B.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia --project=env closed_form.jl
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
# CF_NLTE=NAME (default W30_late_ext): the non-LTE solution compared with (its radii and temperatures); CF_STAR=NAME
# (default star_t0): the stellar spectrum; CF_TEFF, CF_L: the star's T_eff and L (they set its radius)
cfn = get(ENV, "CF_NLTE", "W30_late_ext")
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$cfn.jls"))
d0 = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", get(ENV, "CF_STAR", "star_t0") * ".jls"))
R0 = 8.4352e15; M = 1e6 * Z0.Msun; Mdot = 3Z0.Msun / Z0.yr; Teff = parse(Float64, get(ENV, "CF_TEFF", "7310")); L = parse(Float64, get(ENV, "CF_L", "1.4204e44")); X = 0.76
Rst = sqrt(L / (4π * Z0.σSB * Teff^4)); W(r) = 0.5 * (1 - sqrt(max(1 - (Rst / r)^2, 0.0)))
Ltot0 = sum(0.5 * (d0.L[i] + d0.L[i+1]) * (d0.ν[i+1] - d0.ν[i]) for i in 1:length(d0.ν)-1)
Iν(ν) = (i = clamp(searchsortedlast(d0.ν, ν), 1, length(d0.ν) - 1); t = (log(ν) - log(d0.ν[i])) / (log(d0.ν[i+1]) - log(d0.ν[i]));
         exp(log(d0.L[i]) * (1 - t) + log(d0.L[i+1]) * t) * L / Ltot0 / (4π^2 * Rst^2))
h, c, kB = Z0.h, Z0.c, Z0.kB
νL(l, u) = (NH.En[u] - NH.En[l]) / h
# Γ_i per unit dilution (J = I_ν): 4π ∫ I_ν σ_i/(hν) dν over the edge to 13.6 eV
function Γunit(i)
    ν0 = NH.χn[i] / h; νs = exp.(range(log(ν0 * 1.0001), log(Z0.χH / h * 0.9999), length = 400))
    sum(0.5 * (f(νs[k]) + f(νs[k+1])) * (νs[k+1] - νs[k]) for k in 1:length(νs)-1 for f in (ν -> 4π * Iν(ν) * NH.σbf(i, ν) / (h * ν),))
end
G2, G3 = Γunit(2), Γunit(3)
fuv = parse(Float64, get(ENV, "CF_FUV", "1"))      # CF_FUV=f: far UV (Lyα, Lyβ pumping) × f, a sensitivity test
p2 = fuv * 4 * Iν(νL(1, 2)) / (2h * νL(1, 2)^3 / c^2); p3 = fuv * 9 * Iν(νL(1, 3)) / (2h * νL(1, 3)^3 / c^2)   # (g_u/g_1) I/(2hν³/c²)
αB(T) = 2.59e-13 * (T / 1e4)^(-0.8)          # case-B recombination coefficient (Osterbrock & Ferland; power-law fit)
vff(r) = sqrt(2Z0.G * M / r)
rs = dn.r; Ts = dn.T; nH = [X * Mdot / (4π * r^2 * vff(r)) / Z0.mH for r in rs]
xfull = [st.np / (st.np + sum(st.n)) for st in dn.nl]
τs = zeros(length(rs)); x = zeros(length(rs))
for it in 1:60
    global τs, x
    for k in eachindex(rs)
        f = W(rs[k]) + 3τs[k] * Rst^2 / (4rs[k]^2)                         # dilution × trapping
        Q = f^2 * (p2 * G2 + p3 * G3)                                      # pumping ∝ f and ionizing field ∝ f
        a = Q / (nH[k] * αB(Ts[k]))                                         # x²/(1−x) = a
        x[k] = (-a + sqrt(a^2 + 4a)) / 2
    end
    ne = x .* nH; acc = 0.0; τn = zeros(length(rs))
    for k in length(rs)-1:-1:1; acc += 0.5 * (ne[k] / rs[k]^2 + ne[k+1] / rs[k+1]^2) * 6.6524e-25 * (rs[k+1] - rs[k]); τn[k] = rs[k]^2 * acc; end
    maximum(abs.(τn .- τs)) < 1e-4 && break
    τs = 0.5 .* (τs .+ τn)
end
col(v) = sum(0.5 * (v[k] + v[k+1]) * (rs[k+1] - rs[k]) for k in 1:length(rs)-1)
@printf("star (star_t0) at the surface: J(Lyα)/(2hν³/c²) ×g₂/g₁ = %.3e, Γ₂ = %.3e s⁻¹, Γ₃ = %.3e s⁻¹ (undiluted)\n", p2, G2, G3)
@printf("%8s %7s %9s %9s\n", "r/R0", "T", "x closed", "x full")
for xx in (1.0, 1.2, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 50.0, 80.0, 115.0)
    k = argmin(abs.(rs ./ R0 .- xx)); @printf("%8.2f %7.0f %9.4f %9.4f\n", rs[k] / R0, Ts[k], x[k], xfull[k])
end
@printf("neutral column: closed %.3e, full %.3e cm⁻²; Thomson depth: closed %.3f, full %.3f\n",
        col((1 .- x) .* nH), col((1 .- xfull) .* nH), col(x .* nH .* 6.6524e-25), col([st.ne for st in dn.nl] .* 6.6524e-25))
fuv == 1 && serialize(joinpath(@__DIR__, "..", "lib", "spectra", cfn == "W30_late_ext" ? "closedform_W30.jls" : "closedform_$cfn.jls"), (r = rs, T = Ts, x = x, xfull = xfull, p2 = p2, p3 = p3, G2 = G2, G3 = G3))
