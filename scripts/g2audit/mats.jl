# The two materials at (ρ, T), the network one with its populations in LTE (Saha with the network's own species data).
const SP = CN.SPECIES
lnK(R, P, T) = CN.lnKeq(R, P, T)

"Network species in LTE (ChemNet partition functions and E0): number densities [cm⁻³]."
function net_lte(ρ, T)
    NH = ρ / CN.μH(); NHe = CN.yHe * NH; NLi = CN.yLi * NH
    K1 = exp(lnK([:H], [:Hp, :e], T)); Km = exp(lnK([:H, :e], [:Hm], T)); K2 = exp(lnK([:H, :H], [:H2], T))
    K2p = exp(lnK([:H, :Hp], [:H2p], T))
    S1 = exp(lnK([:He], [:Hep, :e], T)); S2 = exp(lnK([:Hep], [:Hepp, :e], T)); SLi = exp(lnK([:Li], [:Lip, :e], T))
    function parts(ne)
        b = 1 + K1 / ne + Km * ne; a = 2K2 + 2K2p * K1 / ne
        nH = 2NH / (b + sqrt(b^2 + 4a * NH))
        nHp = K1 * nH / ne; nHm = Km * nH * ne; nH2 = K2 * nH^2; nH2p = K2p * nH * nHp
        r1 = S1 / ne; r2 = S2 / ne; he0 = NHe / (1 + r1 + r1 * r2)
        li0 = NLi / (1 + SLi / ne)
        (; nH, nHp, nHm, nH2, nH2p, nHe = he0, nHep = he0 * r1, nHepp = he0 * r1 * r2, nLi = li0, nLip = NLi - li0)
    end
    g(l) = (ne = exp(l); q = parts(ne); q.nHp + q.nH2p + q.nHep + 2q.nHepp + q.nLip - q.nHm - ne)
    lo, hi = log(1e-280), log(NH * (1 + 2CN.yHe) * 1.01)
    for _ in 1:300
        m = 0.5(lo + hi); g(m) > 0 ? (lo = m) : (hi = m); hi - lo < 1e-14 && break
    end
    ne = exp(0.5(lo + hi)); merge(parts(ne), (; ne, NH))
end
"Network EOS (ChemNet.eos: ideal, u = Σ n (c_v kT + E0)/ρ) of the composition s."
function net_eos(ρ, T, s)
    kT = CN.kB * T
    n = (s.nH, s.nHp, s.nHm, s.nH2, s.nH2p, s.nHe, s.nHep, s.nHepp, s.nLi, s.nLip, s.ne)
    nm = (:H, :Hp, :Hm, :H2, :H2p, :He, :Hep, :Hepp, :Li, :Lip, :e)
    u = 0.0; N = 0.0
    for (x, k) in zip(n, nm)
        sp = SP[k]; cv = k === :H2 ? CN.u_H2(T) : sp.cv(T)
        u += x * (cv * kT + sp.E0); N += x
    end
    (u / ρ, N * kT)
end
net_pu(ρ, T) = (s = net_lte(ρ, T); e = net_eos(ρ, T, s); (e[2], e[1]))
tab_pu(ρ, T) = (r = EQ.eos(ρ, T); (r[2], r[1]))

"Thermodynamic response from p(ρ,T), u(ρ,T) by centred differences in ln ρ, ln T (equilibrium composition): χρ, χT, c_V, Γ1, ∇ad, c_P."
function thermo(f, ρ, T; h = 1e-4)
    p0, u0 = f(ρ, T)
    pr1, ur1 = f(ρ * exp(h), T); pr0, ur0 = f(ρ * exp(-h), T)
    pt1, ut1 = f(ρ, T * exp(h)); pt0, ut0 = f(ρ, T * exp(-h))
    χρ = (log(pr1) - log(pr0)) / 2h; χT = (log(pt1) - log(pt0)) / 2h; cV = (ut1 - ut0) / (2h * T)
    G3m1 = p0 * χT / (ρ * T * cV); Γ1 = χT * G3m1 + χρ; ∇ad = G3m1 / Γ1; cP = cV * Γ1 / χρ
    (; p = p0, u = u0, χρ, χT, cV, Γ1, ∇ad, cP)
end

# ---------------------------------------------------------------- network continuum means (PopIIIMaterial.nlte_means arithmetic)
const NLF = Ref{Any}(nothing)
function nlf()
    NLF[] === nothing || return NLF[]
    νs = COp.grid(; n = 1000); n = length(νs) - 1
    ν = [0.5 * (νs[i] + νs[i+1]) for i in 1:n]; dν = [νs[i+1] - νs[i] for i in 1:n]; σ = ν ./ COp.c; λμ = 1e4 ./ σ
    NLF[] = (; ν, dν, σ, λμ, ν3 = ν .^ 3, σhm = COp.σ_Hm_bf.(λμ), ray2 = COp.σ_ray_H2.(λμ .* 1e3), ray1 = COp.σ_ray_H.(λμ .* 1e4))
end
hbf_levels(T, n2) = ntuple(i -> (n = i + 1; n2 * (n^2 / 4) * exp(-13.598 * EVERG * (0.25 - 1 / n^2) / (CN.kB * T))), 29)
function hbf_pre(ν, n1, nn)
    χ = 13.598 * EVERG; σ0 = 7.906e-18; hν = COp.h * ν; α = 0.0
    hν >= χ && (α += n1 * σ0 * (χ / hν)^3)
    for n in 2:30
        hν >= χ / n^2 || continue
        α += nn[n-1] * σ0 * n * (χ / n^2 / hν)^3
    end
    α
end
"""(κ_E, κ_P, κ_R) [cm² g⁻¹] of the network's continuum code (nlte_means, hbfP = true: LTE) and the parts of κ_P and of 1/κ_R
contributions (Planck-weighted absorption by process). `extra(ν, T)` adds an absorption coefficient [cm⁻¹] (prototype terms)."""
function net_means(ρ, T, s; n2 = 4s.nH * exp(-E21 / (CN.kB * T)), hbfP = true, extra = nothing, scaextra = nothing)
    F = nlf(); t = COp.tables(); sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0
    n1 = s.nH; nH2 = s.nH2; nHe = s.nHe; nHp = s.nHp; ne = s.ne; nHm = s.nHm
    cia = nH2 > 1e-3 * n1; acia = 0.0; nnl = hbf_levels(T, n2)
    part = zeros(8)   # Planck-mean parts: cia, Hm bf, Hm ff, H+ ff, H bf, extra
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; σ = F.σ[i]; λμ = F.λμ[i]
        stim = -expm1(-COp.h * ν / (COp.kB * T))
        cia && i % 8 == 1 && (acia = COp.k_h2h2(t, σ, T) * nH2^2 + COp.kcia(t.H2He, σ, T) * nH2 * nHe + COp.kcia(t.H2H, σ, T) * nH2 * n1)
        aff = COp.k_Hm_ff(λμ, T) * n1 * ne * COp.kB * T; apff = 3.692e8 * ne * nHp / (sqrt(T) * F.ν3[i]) * stim
        ahm = F.σhm[i] * stim * nHm
        ab = hbf_pre(ν, n1, nnl) * stim
        ex = extra === nothing ? 0.0 : extra(ν, T)
        a0 = acia + aff + apff + ex
        sca = F.ray2[i] * nH2 + F.ray1[i] * n1 + COp.σT * ne + (scaextra === nothing ? 0.0 : scaextra(ν, T))
        B = COp.Bν(ν, T); D = COp.dBdT(ν, T); hb = hbfP ? ab : 0.0
        sP += (a0 + ahm + hb) * B * dν; sB += B * dν
        part[1] += acia * B * dν; part[2] += ahm * B * dν; part[3] += aff * B * dν; part[4] += apff * B * dν; part[5] += hb * B * dν; part[6] += ex * B * dν
        sR += D / max(a0 + ahm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (; κP = sP / sB / ρ, κR = sD / sR / ρ, parts = part ./ sB ./ ρ)
end
tab_kap(ρ, T) = (SO.kappa_R(ρ, T, XH), SO.kappa_P(ρ, T, XH))
z0_kap(ρ, T; lines = true, lines_R = true) = Z0.kappa_means(ρ, T; comp = Z0.Comp(XH, 1 - XH), lines, lines_R)
