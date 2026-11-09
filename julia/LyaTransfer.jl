# Included in PopIIIMaterial (S6). Lyman-alpha transfer for the three-level H atom (LYART=1, with STARLIGHT=1 SL3=1).
"""
LYART=1: Lyα transfer on the zones from continuum depth τ ≥ LYATAU [20] (inside the photosphere) to the outer edge, replacing the Sobolev
escape and the diluted-starlight pumping of the three-level atom (tl3) by the transferred field. Without it the Lyα photons of the hot,
ionized envelope are part of the gray field, and the cold neutral gas just outside it, which is opaque in the Lyα damping wings
(N_HI ~ 1e24 cm⁻²) but has no Lyα opacity in κ_E, saw only W B_ν(T_c ≈ 4800 K) at Lyα.

- Line: Voigt profile H(a, x)/(√π v_D) (a = A21 λ/(4π v_D), radiative damping; approximation H ≈ e^{−x²} + a (1 − e^{−x²})/(√π x²)), integrated
  opacity (πe²/m_e c) f n_1 (stimulated emission neglected), static atmosphere (no Doppler shifts: the core thermalizes inside a zone, the
  transport across zones is in wings 10³–10⁵ km s⁻¹ wide; the flow is ≲ 200 km s⁻¹). Symmetric profile: v ≥ 0 only, bins log-spaced in
  v from 0 to LYAVMAX [2e5 km/s] (beyond it the Lorentz wing is neither correct nor needed).
- Continuum: TITAN's gray extinction χρ (absorption, source B_ν(T) at ν_Lyα): it bounds the line formation in the hot zones.
- Redistribution: complete in the Doppler core, coherent in the wings (the R_II approximation of Hubeny & Lites 1995 with the elastic
  collision rate ≪ A21: S_ν = S_L + γ_ν (J_ν − J̄), γ_ν = smoothstep of x = v/v_D from 2.5 to 4; LYAGAMMA=0: complete redistribution).
- Atom: the three-level statistical equilibrium of tl3 with the Lyα rates explicit: n_2s (3A21 + Λ_2s + 4 n_e q21 + β) = n_1 (n_e q12 +
  Λ_2s n̄_2γ + 3A21 n̄(J̄)) + α_B n_e n_p, S_L = (2hν³/c²) n_2s/n_1 (2s and 2p in statistical equilibrium, g_2p/g_1 = 3).
- Formal solution: spherical rays (impact parameters at the interfaces plus LYANCORE [6] through the inner sphere, which emits at the inner
  zone's source), zone-constant source; per zone the direction mean of the segment-averaged intensity. ALI with the local operator
  (Λ* = 1 − ⟨(1 − e^{−Δτ})/Δτ⟩); per zone the bins are eliminated exactly: J̄ = J_ext + Λ_eff S_L, so the atom gets
  R_Lyα = 3A21 (1 − Λ_eff) and P12 = 3A21 n̄(J_ext) + Λ_2s n̄_2γ — the Sobolev form of tl3 with the transferred field (the chemistry keeps
  the local, stiff part implicit). The incident part is accumulated as ⟨g I_in⟩ (no cancellation in thick bins).
Energy: as tl3/LYABB (RAYQ ↔ the gray field, where the hot zones' Lyα emission is): the net Lyα absorption E21 (n_1 P12 − n_2s (R+Λ)) of
a network zone is taken from the gray field there. Left out: H₂ absorption at Lyα, H⁻ at 1216 Å, Stark/resonance broadening, the
red–blue wing asymmetry, Doppler shifts.
"""
const LYART = Ref(get(ENV, "LYART", "0") == "1")
const LYATAU = Ref(parse(Float64, get(ENV, "LYATAU", "20")))
const LYAGAMMA = Ref(get(ENV, "LYAGAMMA", "1") != "0")
const LYAOVR = Ref{Vector{Any}}(Any[])
const LYADIAG = Ref{Any}(nothing); const LYAPREV = Ref{Any}(nothing)
const LYA_NU = CNn.c / LYA_LAM; const LYA_B0 = 2CNn.hP * LYA_NU^3 / CNn.c^2
const LYA_KV = 0.0265400 * 0.41641 * LYA_LAM          # (πe²/m_e c) f λ: ∫ χ dv = LYA_KV n_1 [cm⁻¹ · cm s⁻¹]
const LYA_AV = A21LYA * LYA_LAM / (4π)                # a = LYA_AV / v_D
voigtH(a, x) = (x2 = x * x; exp(-x2) + a / sqrt(π) * (x2 < 1e-10 ? 1.0 : -expm1(-x2) / x2))
function lya_vgrid()
    vmax = parse(Float64, get(ENV, "LYAVMAX", "2e5")) * 1e5; h = parse(Float64, get(ENV, "LYAH", "0.2")); v0 = 1e5
    nb = ceil(Int, log1p(vmax / v0) / h); e = [v0 * expm1(h * i) for i in 0:nb]
    vc = [0.5(e[i] + e[i+1]) for i in 1:nb]; dv = diff(e)
    vc, dv
end
smooth01(x) = (x <= 0 ? 0.0 : x >= 1 ? 1.0 : x * x * (3 - 2x))

"""Formal solution of one bin on shells r[1..n+1] (zone-constant source S, extinction χ [cm⁻¹], inner sphere at r[1] emitting Ib).
Fills, per zone, esc = ⟨(1 − e^{−Δτ})/Δτ⟩ and Ain = ⟨g I_in⟩ (direction means over the zone's segments)."""
function lya_formal!(esc, Ain, r, χ, S, Ib, P, W)
    n = length(χ); fill!(esc, 0.0); fill!(Ain, 0.0)
    @inbounds for (j, p) in enumerate(P)
        kmin = 1; while kmin <= n && r[kmin+1] <= p; kmin += 1; end
        kmin > n && continue
        hitcore = p < r[1]
        I = 0.0
        for k in n:-1:kmin                                            # inward
            z1 = sqrt(max(r[k+1]^2 - p^2, 0.0)); z0 = sqrt(max(r[k]^2 - p^2, 0.0)); Δτ = χ[k] * (z1 - z0)
            g = Δτ < 1e-6 ? 1 - 0.5Δτ : -expm1(-Δτ) / Δτ
            w = 0.5 * W[j, k]; esc[k] += w * g; Ain[k] += w * g * I
            I = S[k] + (I - S[k]) * (Δτ > 700 ? 0.0 : exp(-Δτ))
        end
        hitcore && (I = Ib)
        for k in kmin:n                                               # outward
            z1 = sqrt(max(r[k+1]^2 - p^2, 0.0)); z0 = sqrt(max(r[k]^2 - p^2, 0.0)); Δτ = χ[k] * (z1 - z0)
            g = Δτ < 1e-6 ? 1 - 0.5Δτ : -expm1(-Δτ) / Δτ
            w = 0.5 * W[j, k]; esc[k] += w * g; Ain[k] += w * g * I
            I = S[k] + (I - S[k]) * (Δτ > 700 ? 0.0 : exp(-Δτ))
        end
    end
    nothing
end
"Ray impact parameters and the μ weight of each ray in each shell (bands between neighbouring rays, at the shell's mid radius)."
function lya_rays(r; nc = parse(Int, get(ENV, "LYANCORE", "6")))
    n = length(r) - 1
    P = vcat([r[1] * sqrt((i - 0.5) / nc) for i in 1:nc], r[1:n])
    edges = vcat(0.0, [0.5(P[i] + P[i+1]) for i in 1:length(P)-1], Inf)
    W = zeros(length(P), n)
    for k in 1:n
        rm = 0.5(r[k] + r[k+1]); μ(p) = sqrt(max(1 - (min(p, rm) / rm)^2, 0.0))
        for j in eachindex(P)
            edges[j] >= rm && break
            P[j] < r[k+1] || continue
            W[j, k] = μ(edges[j]) - μ(edges[j+1])
        end
        s = sum(view(W, :, k)); s > 0 && (W[:, k] ./= s)
    end
    P, W
end

"""Lyα transfer of this step (see LYART): sets LYAOVR[][k] = (; RL, P12, ...) for the network zones and LYADIAG[]."""
function lya_transfer!(m::NetGas, a, b, tauz, kph, rates, diag)
    t0 = time()
    ov = LYAOVR[]; length(ov) == length(TM.v(:dn)) || (LYAOVR[] = ov = Any[nothing for _ in TM.v(:dn)]); fill!(ov, nothing)   # (rates(k) below: Sobolev values)
    dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn); chi = TM.v(:chifn); fr = TM.v(:frn)
    kin = a; for k in kph:-1:a; tauz[k] >= LYATAU[] && (kin = k; break); end
    K = kin:b
    vc, dv = lya_vgrid(); nb = length(vc)
    # the transfer's layers: one per zone, except thick zones near the surface (τ_zone > 1, τ < LYASUBTOP [30] at their outer edge),
    # whose unresolved photosphere is rebuilt as LYANSUB [24] sublayers (geometric in τ from 0.01 at the top) with the gray Eddington
    # T⁴ = ¾ T_eff⁴ (τ + 2/3) of the local flux F = σ T_eff⁴ at the zone's outer interface, capped at the zone's T, and table
    # composition at that T: the zone-constant 16 kK of a τ = 600 zone radiated 25× the star's flux in the Lyα band
    subtop = parse(Float64, get(ENV, "LYASUBTOP", "30")); nsub = parse(Int, get(ENV, "LYANSUB", "24"))
    L_k = Int[]; L_T = Float64[]; L_sub = Bool[]; r = [rn[kin]]
    for k in K
        τz = chi[k] * dn[k] * (rn[k+1] - rn[k]); τtop = tauz[k] - τz
        Teff = (abs(fr[k+1]) / 5.670374419e-5)^0.25
        if nsub > 0 && m.eqzone[k] && τz > 1 && τtop < subtop && Teff > 100      # (no flux yet, e.g. right after ENVSPLIT: no rebuild)
            τ0 = 0.01
            f(q) = τ0 * (q^nsub - 1) / (q - 1) - τz; lo, hi = 1.0 + 1e-9, 100.0
            for _ in 1:200; mid = 0.5(lo + hi); f(mid) > 0 ? (hi = mid) : (lo = mid); end
            q = 0.5(lo + hi); te = [τ0 * (q^j - 1) / (q - 1) for j in 0:nsub]; te[end] = τz       # τ below the zone's top
            for j in nsub:-1:1                                                                    # inner → outer
                tm = τtop + 0.5(te[j] + te[j+1])
                push!(L_k, k); push!(L_T, min(tn[k], max((0.75 * Teff^4 * (max(tm, 0.0) + 2 / 3))^0.25, 100.0))); push!(L_sub, true)
                push!(r, j == 1 ? rn[k+1] : rn[k+1] - te[j] / (chi[k] * dn[k]))
            end
        else
            push!(L_k, k); push!(L_T, tn[k]); push!(L_sub, false); push!(r, rn[k+1])
        end
    end
    n = length(L_k)
    n1 = zeros(n); Q0 = zeros(n); Dc = zeros(n); nb2 = zeros(n); χc = zeros(n); Bc = zeros(n); SL = zeros(n); Psob = zeros(n); Rsob = zeros(n)
    χL = zeros(n, nb); wq = zeros(n, nb); γ = zeros(n, nb); rr = zeros(n, nb)
    prev = LYAPREV[]
    for i in 1:n
        k = L_k[i]; T = L_T[i]
        p = nothing
        if !L_sub[i]; rates(k); p = diag.tl3[k]; end
        x = m.eqzone[k] ? x_from_table(m, dn[k], T) : m.x[k]; d = CNn.full(x, m.net)
        nH = dn[k] / CNn.μH(); n1[i] = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
        if p === nothing                               # (sublayer, or no star light in this zone: LTE radiation at the gas T)
            nbT = 1 / expm1(min(E21 / (CNn.kB * T), 700.0)); p = (; β = beta2s_bb(T), RL = 3A21LYA, ne, P12 = (3A21LYA + LAM2S) * nbT)
        end
        nb2[i] = p.P12 / (p.RL + LAM2S); Psob[i] = p.P12; Rsob[i] = p.RL
        Q0[i] = n1[i] * (ne * q12_H(T) + LAM2S * nb2[i]) + alphaB_recfast(T) * ne * np
        Dc[i] = LAM2S + 4ne * q21_H(T) + p.β
        χc[i] = chi[k] * dn[k]; Bc[i] = LYA_B0 / expm1(min(E21 / (CNn.kB * T), 700.0))
        n2s0 = n2_tl3(T, n1[i], ne, np, p) / 4; SL[i] = n1[i] > 0 ? LYA_B0 * n2s0 / n1[i] : Bc[i]
        (!L_sub[i] && prev !== nothing && haskey(prev, k)) && (SL[i] = prev[k])        # warm start (outer scattering halo)
        vD = sqrt(2CNn.kB * T / 1.6735575e-24); av = LYA_AV / vD
        for ib in 1:nb
            φ = voigtH(av, vc[ib] / vD) / (sqrt(π) * vD)
            χL[i, ib] = LYA_KV * n1[i] * φ; wq[i, ib] = 2φ * dv[ib]
            γ[i, ib] = LYAGAMMA[] ? smooth01((vc[ib] / vD - 2.5) / 1.5) : 0.0
            rr[i, ib] = χL[i, ib] / (χL[i, ib] + χc[i] + 1e-300)
        end
        s = sum(view(wq, i, :)); s > 0 && (wq[i, :] ./= s)
    end
    oml, Jext, Jbar, SL, n2s, it, δ = lya_solve(r, n1, Q0, Dc, χc, Bc, SL, χL, wq, γ, rr)
    LYAPREV[] = Dict(L_k[i] => SL[i] for i in 1:n if !L_sub[i])
    zi = Dict(L_k[i] => i for i in 1:n if !L_sub[i])
    for k in K
        haskey(zi, k) || continue; i = zi[k]
        (m.eqzone[k] || n1[i] <= 0) && continue
        ov[k] = (; RL = 3A21LYA * oml[i], P12 = 3A21LYA * Jext[i] / LYA_B0 + LAM2S * nb2[i], nJ = Jbar[i] / LYA_B0, nS = SL[i] / LYA_B0, Psob = Psob[i], Rsob = Rsob[i])
    end
    # outward Lyα luminosity at the outer edge (one more formal pass for the flux is not needed: diagnostics only from J̄)
    sel(v) = [haskey(zi, k) ? v[zi[k]] : NaN for k in K]
    nsub_ = count(L_sub); ksub = unique(L_k[L_sub])
    LYADIAG[] = (; kin, n = length(K), nlayer = n, ksub, Tsurf = nsub_ > 0 ? L_T[findlast(L_sub)] : NaN, nb, it, δ, time = time() - t0, K,
                 nJ = sel(Jbar ./ LYA_B0), nS = sel(SL ./ LYA_B0), oml = sel(oml), Jext = sel(Jext ./ LYA_B0), Psob = sel(Psob), Rsob = sel(Rsob))
    nothing
end

"""Lyα ALI on shells r (see LYART): per zone n_1, Q0 (n_2s creation without Lyα), Dc (n_2s destruction without Lyα), continuum χc, Bc,
initial S_L, and per zone and bin line extinction χL, profile weights wq, coherence γ, line share rr. Returns (1 − Λ_eff, J_ext, J̄, S_L, n_2s, iterations, max |Δ ln S_L|)."""
function lya_solve(r, n1, Q0, Dc, χc, Bc, SL, χL, wq, γ, rr)
    n, nb = size(χL); SL = copy(SL)
    P, W = lya_rays(r)
    Sν = repeat(SL, 1, nb); esc = zeros(n, nb); Ain = zeros(n, nb)
    Jext = zeros(n); oml = ones(n); Jbar = zeros(n); n2s = zeros(n)
    maxit = parse(Int, get(ENV, "LYAMAXIT", "40")); tol = parse(Float64, get(ENV, "LYATOL", "1e-3")); it = 0; δ = Inf
    ng = parse(Int, get(ENV, "LYANG", "4")); hist = Vector{Vector{Float64}}()
    while it < maxit && δ > tol
        it += 1
        Threads.@threads for ib in 1:nb
            χ = [χL[i, ib] + χc[i] for i in 1:n]; S = [rr[i, ib] * Sν[i, ib] + (1 - rr[i, ib]) * Bc[i] for i in 1:n]
            e_ = zeros(n); A_ = zeros(n)
            lya_formal!(e_, A_, r, χ, S, S[1], P, W)
            esc[:, ib] .= e_; Ain[:, ib] .= A_
        end
        δ = 0.0; SLold = copy(SL)
        for i in 1:n
            G0 = 0.0; G2 = 0.0; Ge = 0.0
            for ib in 1:nb
                ρ_ = rr[i, ib]; E = 1 - ρ_ + ρ_ * esc[i, ib]; L = ρ_ * (1 - esc[i, ib]); g = γ[i, ib]; dd = 1 - g + g * E
                A = Ain[i, ib] + (1 - esc[i, ib]) * (1 - ρ_) * Bc[i]
                w = wq[i, ib]; G0 += w * A / dd; G2 += w * g * L / dd; Ge += w * E / dd
            end
            Jext[i] = max(G0 / (1 + G2), 0.0); oml[i] = Ge / (1 + G2)
            if n1[i] > 0
                n2s[i] = (Q0[i] + n1[i] * 3A21LYA * Jext[i] / LYA_B0) / (3A21LYA * oml[i] + Dc[i])
                s_new = LYA_B0 * n2s[i] / n1[i]
            else
                s_new = Bc[i]
            end
            Dc[i] > 1e-3 * 3A21LYA * oml[i] && (δ = max(δ, abs(log(max(s_new, 1e-300) / max(SL[i], 1e-300)))))   # (scattering halo: warm-started over steps)
            SL[i] = s_new; Jbar[i] = Jext[i] + (1 - oml[i]) * SL[i]
            for ib in 1:nb
                ρ_ = rr[i, ib]; E = 1 - ρ_ + ρ_ * esc[i, ib]; g = γ[i, ib]; dd = 1 - g + g * E
                A = Ain[i, ib] + (1 - esc[i, ib]) * (1 - ρ_) * Bc[i]
                Sν[i, ib] = max((SL[i] - g * Jbar[i] + g * A) / dd, 0.0)
            end
        end
        if ng > 0                                       # Ng (1974) acceleration of order 2 on ln S_L; the bin sources follow S_L
            push!(hist, log.(max.(SL, 1e-300))); length(hist) > 4 && popfirst!(hist)
            if length(hist) == 4 && it % ng == 0
                x0, x1, x2, x3 = hist[4], hist[3], hist[2], hist[1]
                d0 = x0 .- x1; d1 = x1 .- x2; d2 = x2 .- x3; u1 = d0 .- d1; u2 = d0 .- d2
                A11 = sum(u1 .* u1); A12 = sum(u1 .* u2); A22 = sum(u2 .* u2); b1 = sum(d0 .* u1); b2 = sum(d0 .* u2)
                det = A11 * A22 - A12^2
                if det > 1e-14 * A11 * A22 && isfinite(det)
                    α = (b1 * A22 - b2 * A12) / det; β_ = (b2 * A11 - b1 * A12) / det
                    xn = (1 - α - β_) .* x0 .+ α .* x1 .+ β_ .* x2
                    for i in 1:n
                        sn = exp(xn[i]); f = sn / max(SL[i], 1e-300)
                        for ib in 1:nb; Sν[i, ib] *= f; end
                        SL[i] = sn
                    end
                    empty!(hist)
                end
            end
        end
    end
    oml, Jext, Jbar, SL, n2s, it, δ
end
