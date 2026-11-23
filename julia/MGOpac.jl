# ------------------------------------------------------------------ MGROUPS: material side (group opacity records, group photo-rates)
# Included into PopIIIMaterial. Design: docs/POPIII_MULTIGROUP.md. Nothing here runs unless MGROUPS is set (MGON false: no call).
const MGON = haskey(ENV, "MGROUPS") && ENV["MGROUPS"] != "0"
"Group edges of the code's own thresholds [eV]: H⁻ photodetachment (λ₀ = 1.6419 μm, ContinuumOpacity.σ_Hm_bf), H n = 3, 2 (χ_H/n²), the
lower end of the Lyman–Werner band (11.2 eV, as RADGROUPS; Draine & Bertoldi 1996), H, He I, He II ionisation (Z0Physics χ_H, χ_He1, χ_He2:
the thresholds of the G2 continuum that supplies the group opacities)."
mg_thresholds() = (Hm = COp.h * COp.c * (1e4 / 1.6419) / EVERG, Pa = Z0Physics.χH / 9 / EVERG, Ba = Z0Physics.χH / 4 / EVERG, LW = 11.2,
                   H = Z0Physics.χH * (1 - 1 / 1600) / EVERG, He1 = Z0Physics.χHe1 / EVERG, He2 = Z0Physics.χHe2 / EVERG)
# (the H edge at the onset of the dissolved Lyman pseudo-continuum, q = 1/n² − ν/ν_H ≤ 1/1600 as G2Opacity's DIS window (Hubeny, Hummer &
#  Lanz 1994): that 6e-4-wide sliver of ground-state absorption belongs with the ionising group; in the LW band's Planck mean it raised κ_E,LW
#  ×1e5 in the LTE zones below the photosphere and thermalised the band there, one zone deeper and 2500 K hotter than the front)
function mg_edges(s)
    t = mg_thresholds()
    (occursin(',', s) || occursin('.', s)) && return parse.(Float64, split(s, ','))
    n = parse(Int, s)
    n == 4 ? [t.Hm, t.LW, t.H] : n == 6 ? [t.Hm, t.LW, t.H, t.He1, t.He2] :
    error("MGROUPS=$n: edge sets for 4 (IR | optical | LW | ionising) and 6 (+ He I, He II) groups, or give the edges in eV")
end
"Group containing photon energy e [eV]."
mg_gof(e) = searchsortedlast(TX.MGEDGE[], e) + 1
"Roles: the group of the H⁻ band (just above its threshold), of the Balmer continuum, of Lyα (E21), of the LW rate (12.87 eV), the first ionizing group."
const MGROLE = Ref{Any}(nothing)
function mg_roles!()
    t = mg_thresholds()
    MGROLE[] = (; Hm = mg_gof(t.Hm * 1.0001), Ba = mg_gof(t.Ba * 1.0001), Lya = mg_gof(E21 / EVERG), LW = mg_gof(12.87), H = mg_gof(t.H * 1.0001))
    r = MGROLE[]      # the LW closure and the H(1s) pair need the LW band as a group of its own, a group below it for the fluorescence, an ionising group above it
    (1 < r.LW < r.H <= length(TX.MGEDGE[]) + 1 && mg_gof(11.2 * 1.0001) == r.LW) || error("MGROUPS: the edges must contain the LW band's lower end (11.2 eV) and the H edge ($(t.H) eV), with a group below the LW band")
    # Raman pairs (TX.MGRAMP): every group between Lyα and the Lyman limit (G2O.lwband) with the groups of its Stokes partners hν − E_21
    e = TX.MGEDGE[]; G = length(e) + 1; e21 = E21 / EVERG; hz = EVERG / CNn.hP; pairs = Tuple{Int,Int,Float64,Float64}[]
    for gL in 2:G-1
        G2O.lwband(e[gL-1] * hz, e[gL] * hz) || continue
        for gp in 1:G
            e1 = max(e[gL-1], (gp == 1 ? 0.0 : e[gp-1]) + e21); e2 = min(e[gL], (gp == G ? Inf : e[gp]) + e21)
            e2 > e1 && push!(pairs, (gL, gp, e1, e2))
        end
    end
    TX.mg_rampairs!(pairs)
    MGROLE[]
end
const MGGIDX = Ref{Vector{Int}}(Int[])
const MGTMIN = 3000.0                                         # in-group shape temperature floor of the group means [K] (docs/POPIII_MULTIGROUP.md §4, §19)

const MGEHZ = Ref(Float64[])                                  # interior group edges [Hz]
"The record of a zone without group data (MGROUPS off, or the group means failed): `nlo_it` then writes the default record."
const MGNOV = Float64[]
"""
Set-up of the material side (popiii_m2c mg_install!, after TX.mg_setup!): roles and Raman pairs, the groups of the G2 nodes, the line
table, per-zone diagnostics, the LW group of the H₂ closure (`m_free`: the network's free species).
"""
function mg_init!(m_free = Symbol[])
    n = length(TM.v(:dn))
    MGEHZ[] = TX.MGEDGE[] .* (EVERG / CNn.hP); MGGIDX[] = G2O.group_index(MGEHZ[])
    mg_roles!()
    MGWREC[] = Any[nothing for _ in 1:n]
    # stage 2c: H₂ Lyman–Werner absorption in the LW group (mg_h2lw!), its fluorescence into the group below
    ro = MGROLE[]; e = TX.MGEDGE[]
    if (:H2 in m_free) && ro.LW > 1 && abs(e[ro.LW-1] - 11.2) < 1e-6 && ro.H == ro.LW + 1
        TX.MGH2G[] = ro.LW; TX.MGH2T[] = ro.LW - 1; MGLWDNU[] = MGEHZ[][ro.LW] - MGEHZ[][ro.LW-1]
    else
        TX.MGH2G[] = 0
    end
    mg_linetab!()
    nothing
end

"""
The lines' share of the grey Planck-type means and its distribution over the groups (stage 2b).

The network zones' grey κ_P (and κ_E) is the continuum code's mean of the zone's composition times (P_t/P_q)^w (NLOIT: the opacity
table's Planck mean over the continuum code's in equilibrium, w = nlo_w(T)). The table's Planck mean (stellaropacity_build.jl) is
log-blended, weight b(T) = smoothstep((log₁₀T − 3.8)/0.3), from the low-temperature continuum mean (no lines) to
Z0Physics.kappa_means(lines = true, lines_R = false): continuum plus the Planck-weighted strengths of the H, He I and He II lines. So
the grey mean carries the factor e^{w b λ}, λ(ρ, T) = ln(P_L/P_C) of Z0Physics with and without the lines, which is the lines' part,
φ = 1 − e^{−w b λ} of the mean; the rest of P_t/P_q is the continuum's calibration. λ and the lines' group shares s_g (each line
in the group of its hν₀, weight u₀³/(e^{u₀} − 1)(1 − e^{−u₀}) S as in kappa_means) are tabulated here on (log T, log ρ) from the same
functions and composition as the table. Lines inside a group whose H(1s) line wings are carried exactly (G2O.lwband: Lyβ and above,
Rayleigh + Raman of HScatter, docs §16) are left out of the shares, which are renormalised: the Lyman lines of that band scatter
(or Raman-convert), they do not thermalise, and their LTE Planck mean in κ_E,LW would absorb the band below the photosphere.
"""
const MGLT = Ref{Any}(nothing)
const MGLTB = (3.8, 0.3)                                      # the table's blend (stellaropacity_build.jl: wlo = 1 − smoothstep((lT − 3.8)/0.3))
function mg_linetab!()
    G = TX.MGN[]; e = MGEHZ[]; Z = Z0Physics
    X = 1 / (1 + 4PC.yHe); comp = Z.Comp(X, 1 - X)
    lT = collect(MGLTB[1]:0.05:5.0); lr = collect(-18.0:1.0:-4.0); nt = length(lT); nr = length(lr)
    excl = [g > 1 && g < G && G2O.lwband(e[g-1], e[g]) for g in 1:G]
    λ = zeros(nt, nr); S = zeros(G, nt, nr); t0 = time()
    Threads.@threads for idx in 1:nt*nr
        i = (idx - 1) % nt + 1; j = (idx - 1) ÷ nt + 1; T = 10.0^lT[i]; ρ = 10.0^lr[j]
        _, pc = Z.kappa_means(ρ, T; comp, lines = false); _, pl = Z.kappa_means(ρ, T; comp, lines = true, lines_R = false)
        λ[i, j] = (pc > 0 && pl > pc) ? log(pl / pc) : 0.0
        sh = zeros(G)
        for (ν0, S_) in Z.line_strengths(Z.opacity_state(ρ, T; comp), T)
            u0 = Z.h * ν0 / (Z.kB * T); u0 > 300 && continue
            g = mg_gof(Z.h * ν0 / EVERG); excl[g] && continue
            sh[g] += u0^3 / expm1(u0) * (-expm1(-u0)) * S_
        end
        tot = sum(sh)
        for g in 1:G; S[g, i, j] = (isfinite(tot) && tot > 0) ? sh[g] / tot : TX.mg_bfrac(g, T); end
    end
    MGLT[] = (; t0 = lT[1], dt = 0.05, nt, r0 = lr[1], dr = 1.0, nr, λ, S)
    @printf("MGROUPS line table: %d × %d points (log T %.2f–%.2f, log ρ %.0f–%.0f), λ = ln(P_L/P_C) %.2f–%.2f, %.1f s\n", nt, nr, lT[1], lT[end], lr[1], lr[end],
            minimum(λ), maximum(λ), time() - t0)
    nothing
end
"(φ, s_g) at (ρ, T): the lines' part of the grey Planck-type means and the lines' group shares (bilinear in log T, log ρ; clamped). φ = 0 below the table's blend."
function mg_lines(ρ, T)
    L = MGLT[]; lt = log10(T)
    (L === nothing || lt <= MGLTB[1]) && return 0.0, MGNOV
    w = nlo_w(T); w > 0 || return 0.0, MGNOV
    x = clamp((lt - L.t0) / L.dt, 0.0, L.nt - 1 - 1e-9); i = floor(Int, _val(x)) + 1; a = x - (i - 1)
    y = clamp((log10(ρ) - L.r0) / L.dr, 0.0, L.nr - 1 - 1e-9); j = floor(Int, _val(y)) + 1; b = y - (j - 1)
    w00 = (1 - a) * (1 - b); w10 = a * (1 - b); w01 = (1 - a) * b; w11 = a * b
    λ = w00 * L.λ[i, j] + w10 * L.λ[i+1, j] + w01 * L.λ[i, j+1] + w11 * L.λ[i+1, j+1]
    bl = clamp((lt - MGLTB[1]) / MGLTB[2], 0.0, 1.0); bl = bl * bl * (3 - 2bl)
    φ = -expm1(-bl * w * λ)
    φ > 0 || return 0.0, MGNOV
    φ, [w00 * L.S[g, i, j] + w10 * L.S[g, i+1, j] + w01 * L.S[g, i, j+1] + w11 * L.S[g, i+1, j+1] for g in 1:TX.MGN[]]
end

"""
The group means of a network composition at (ρ, T): G2Opacity.group_fast (absorption weighted with B_ν(T_w), T_w the colour temperature
as NLOPACJ; H bound–free in κ_E/κ_P only in deep zones, NLOHBFDEEP). Its 5th–7th values are the grey means of `means_fast` on the same
nodes: `nlo_ratio` takes the zone's grey means from this call, so the group ratios and the grey means belong to one evaluation.
`nothing` if the evaluation throws.
"""
mg_group(ρ, T; kw...) = try
    G2O.group_fast(ρ, T, MGGIDX[], TX.MGN[]; Tmin = MGTMIN, edges = MGEHZ[], kw...)
catch
    nothing
end
"""
The zone's group record (TX.mg_nv() numbers, layout: MGRad.jl) from the group means y of `mg_group` at (ρ, T): the ratios of the group
means to the grey means of the same call, continuum and lines (`mg_lines`):
  r_P,g = (1 − φ) κ_P,g/κ_P + φ ℓ_g,  ℓ_g = s_g/b_g(T)        (Σ_g b_g(T) r_P,g = 1),
  r_E,g = (1 − φ) κ_E,g/κ_E + φ ℓ_g/Σ_g′ b_g′(T_w) ℓ_g′       (Σ_g b_g(T_w) r_E,g = 1),   r_Ec,g = (1 − φ) κ_E,g/κ_E,
r_R,g = χ_g/χ (the table's Rosseland mean has no lines); the absolute H(1s) (and n ≥ 2) bound–free absorption K_1 of the
groups at and above χ_H (LW band) where H bound–free is not in κ_E (case-B zones); the Milne emission coefficient
A_1 = σ̄_1s,g n_e n_p λ_e³/ρ of the same groups and zones (emission 4π A_1 e^{χ_H/kT} b_g(T) σT⁴/π: spontaneous + induced recombination
to the ground state in detailed balance with K_1; in LTE, n_1 = n_e n_p λ_e³ e^{χ_H/kT}, emission = absorption of a Planck field);
κ_Ram and the Raman pairs' κ_Ram C, κ_Ram ω (MGRad.jl MGRAMP; n_2s = n_2/4 of the three-level atom). A ratio that is not finite and
positive is 1 (grey).
"""
function mg_pack(ρ, T, Tw, n1, ne, np, n2, hbf, y)
    G = TX.MGN[]; r = MGROLE[]
    κE, κP, χ, k1, gE, gP, gχ, k2, kram, s1 = y
    o = zeros(typeof(one(ρ) * one(T) * one(Tw) * one(n1) * one(ne) * one(np) * one(n2) * one(gE)), TX.mg_nv())        # (generic in the number type: Duals give the record's exact slopes)
    ok(v) = isfinite(v) && v > 0
    φ, sh = mg_lines(ρ, T); ℓ = MGNOV; nrm = 0.0
    if φ > 0
        ℓ = [(b = TX.mg_bfrac(g, T); b > 1e-150 ? sh[g] / b : 0.0) for g in 1:G]
        for g in 1:G; nrm += TX.mg_bfrac(g, Tw) * ℓ[g]; end
        ok(nrm) || (φ = 0.0)
    end
    λ3 = (CNn.hP^2 / (2π * CNn.me * CNn.kB * T))^1.5
    for g in 1:G
        rE = gE > 0 ? κE[g] / gE : 1.0; ok(rE) || (rE = 1.0)
        rP = gP > 0 ? κP[g] / gP : 1.0; ok(rP) || (rP = 1.0)
        rR = (gχ > 0 && χ[g] > 0) ? χ[g] / gχ : 1.0; ok(rR) || (rR = 1.0)
        j = TX.mgv(g, 1)
        o[j] = φ > 0 ? (1 - φ) * rE + φ * ℓ[g] / nrm : rE
        o[j+1] = φ > 0 ? (1 - φ) * rP + φ * ℓ[g] : rP
        o[j+2] = rR; o[j+3] = (1 - φ) * rE
        # H n ≥ 2 bound–free (three-level n = 2, Boltzmann × HM88 above it) of the groups at and above the LW band: NLOPAC leaves all H
        # bound–free out of κ_E in case-B zones because the three-level atom carries the n = 2 photo-exchange as a chemistry heat with the
        # optical group's dilution — right for the Balmer continuum, but the LW and ionising photons absorbed by n ≥ 2 were then not removed
        o[j+4] = ((!hbf && g >= r.H && isfinite(k1[g])) ? k1[g] : 0.0) + ((!hbf && g >= r.LW && isfinite(k2[g])) ? k2[g] : 0.0)
        o[j+5] = (!hbf && g >= r.H && ok(s1[g])) ? s1[g] * ne * np * λ3 / ρ : 0.0
        o[j+6] = ok(kram[g]) ? kram[g] : 0.0
    end
    b2 = n1 > 0 ? n2 / (4n1) * exp(min(E21 / (CNn.kB * Tw), 600.0)) : 0.0; isfinite(b2) || (b2 = 0.0)
    P = TX.MGRAMP[]; i = 1
    while i <= length(P)
        gL = P[i][1]; i1 = i; while i1 < length(P) && P[i1+1][1] == gL; i1 += 1; end       # (the pairs of one gL are consecutive, mg_roles!)
        KR = o[TX.mgv(gL, 7)]; tot = 0.0
        for q in i:i1; tot += TX.mg_band(P[q][3], P[q][4], Tw); end
        for q in i:i1
            bs = TX.mg_band(P[q][3], P[q][4], Tw); jp = TX.mgvp(q, 1)
            o[jp] = KR * (bs > 1e-150 ? b2 * bs / TX.mg_bfrac(P[q][2], Tw) : 0.0)
            # (the shares where the band's Planck weight underflows, T_w ≲ 250 K: by width; no light there)
            o[jp+1] = KR * (tot > 1e-250 ? bs / tot : (P[q][4] - P[q][3]) / (TX.MGEDGE[][gL] - TX.MGEDGE[][gL-1]))
        end
        i = i1 + 1
    end
    o[end-1] = n2; o[end] = Tw
    o
end
"The new level's records of all zones (MGVN) from the present state (install, restart): the sweep of fill_opac! without touching TITAN's arrays."
function mg_fill!(m::NetGas)
    OPAC[] === nothing && return nothing
    dn, tn = TM.v(:dn), TM.v(:tn); V = TX.MGVN[]
    TX.zone_sweep!(TM.krange(), k -> @inbounds(dn[k] > 0 && tn[k] > 0)) do k
        TX.material_opacg!(view(V, :, k), m, dn[k], tn[k])
    end
    nothing
end
"(diagnostic) per zone of the last chemistry: (grey W, group W_g, T_c) of the star light (MGROUPS only; not in the checkpoint)."
const MGWREC = Ref(Any[])
# (sized by mg_init!, serially: the chemistry's zones may run on CHEMTHREADS tasks)
mg_wrec!(k, W, Wg, Tc) = (k <= length(MGWREC[]) && (MGWREC[][k] = (W, collect(Wg), Tc)); nothing)
"""
H(1s) bound–free sums of group g over the fast grid's nodes, with the stimulated factor of the gas temperature T and the in-group shape
B_ν(Ts): (S1, SN, SB) = (Σ σ_1s (1 − e^{−hν/kT}) B_ν(Ts) w, the same with 1/(hν), Σ B_ν(Ts) w), all three times e^{xL}, xL = hν_L/kTs of
the group's lower edge (returned fourth), so that the sums and the derivatives of their ratios stay finite where B_ν(Ts) underflows.
S1/SB is the group's energy absorption cross-section per ground-state atom (K_1 of the record, G2Opacity.group_fast), S1/SN the mean
energy ε̄_g of the photons it absorbs.
"""
function mg_h1s(g, T, Ts)
    F = G2O.fastgrid(); gi_ = MGGIDX[]                 # (the grid's nodes ascend in ν: a group's nodes are one range)
    _mg_h1s(F.ν, F.w, F.SH, searchsortedfirst(gi_, g):searchsortedlast(gi_, g), CNn.hP * MGEHZ[][g-1], T, Ts)
end
function _mg_h1s(νs, ws, SH, rng, hνL, T, Ts)          # (function barrier: the grid is held untyped)
    kT = CNn.kB * T; kTs = CNn.kB * Ts; xL = hνL / kTs
    S1 = zero(T) * zero(Ts); SN = S1; SB = zero(Ts)
    @inbounds for i in rng
        ν = νs[i]; x = CNn.hP * ν / kTs; ForwardDiff.value(x - xL) > 700 && continue
        B = 2CNn.hP * ν^3 / CNn.c^2 * exp(xL - x) / (-expm1(-x)) * ws[i]; SB += B
        σ = SH[1, i]; σ > 0 || continue
        a = σ * (-expm1(-CNn.hP * ν / kT)) * B; S1 += a; SN += a / (CNn.hP * ν)
    end
    S1, SN, SB, xL
end
"In-group shape temperature of the H(1s) absorption mean K_1 (the Tw, Tmin of the record's `mg_group` call: T_c with NLOPACJ, floored at MGTMIN)."
mg_t1s(T, Tc) = max(NLOPACJ[] ? Tc + zero(T) : T + zero(Tc), MGTMIN)
"""
The H(1s) photoionisation / ground-state recombination pair of the ionising groups, photon-conserving (docs/POPIII_MULTIGROUP.md §18).
The group rows exchange ENERGY with H(1s): absorption c E_g n_1 S1/SB (shape T1s = `mg_t1s`), Milne emission 4π n_e n_p λ_e³ s1_g(T) (σT⁴/π) X_g(T)
(mg_theta!). The network needs photon NUMBERS. Both are the row's energy term divided by the same mean photon energy ε̄_g = S1/SN of the
group's absorbed light:
    k_ion [s⁻¹ per H]      = Σ_g c u_g SN/SB                       (`mg_kion`)
    α_1   [cm³ s⁻¹]        = Σ_g 4π λ_e³ s1_g(T) (σT⁴/π) X_g(T) SN/S1   (`mg_arec`)
so that (ionisations − ground-state recombinations) × ε̄_g = the group's net energy exchange with H(1s), identically: where the group is
opaque and its row is in on-the-spot balance (absorption = emission) the pair cancels exactly, whatever T and T_c are, and every net
ionisation costs the transported field the energy ε̄_g. Exact for star light of colour T_c; the gross recombination coefficient differs from
the Milne value by ε̄_em(T)/ε̄_g(T1s) (a few % for |T − T_c| ~ 10³ K; it cancels on the spot). The earlier forms counted the two sides
with different weights (photon count with B(T_c) on the grid over the analytic band fraction, emission at the iterate's T against a
recombination coefficient frozen at the start-of-step T): the on-the-spot loop then had a gain ≠ 1 and created or destroyed ionisations at
a rate ∝ the gross loop rate (t_mb09–t_mb13).
"""
function mg_kion(cu, T, Tc)
    Ts = mg_t1s(T, Tc); s = zero(cu[1]) * zero(T)
    for g in (MGROLE[].H::Int):TX.MGN[]
        _, SN, SB = mg_h1s(g, T, Ts)
        ForwardDiff.value(SB) > 0 && (s += cu[g] * SN / SB)
    end
    s
end
function mg_arec(T, Tc)
    Ts = mg_t1s(T, Tc); λ3 = (CNn.hP^2 / (2π * CNn.me * CNn.kB * T))^1.5; plf = 5.670374419e-5 * T^4 / π; s = zero(T)
    for g in (MGROLE[].H::Int):TX.MGN[]
        S1, SN, _ = mg_h1s(g, T, Ts); ForwardDiff.value(S1) > 0 || continue
        E1, _, EB, xL = mg_h1s(g, T, T)
        # (group_fast takes the emission mean with the weight B_ν(max(T, MGTMIN)) where the group's Σ B_ν(T) w underflows, < 1e-280)
        (ForwardDiff.value(EB) > 0 && log(ForwardDiff.value(EB)) - ForwardDiff.value(xL) > log(1e-280)) || ((E1, _, EB) = mg_h1s(g, T, max(T, MGTMIN)))
        ForwardDiff.value(EB) > 0 || continue
        s += 4π * λ3 * (E1 / EB) * plf * TX.mg_xfrac(g, T, TX.MGXH[]) * SN / S1
    end
    s
end
"The two reactions (zero unless the zone's env carries the group fields `mgq` with `h1`: star-lit case-B zones)."
function mg_h1s_reactions()
    on(env) = hasproperty(env, :mgq) && env.mgq.h1
    [CNn.Reaction("H + γ → H+ + e (ionising groups, k_ion)", [:H], [:Hp, :e], (T, env) -> on(env) ? mg_kion(env.mgq.cu, T, env.mgq.Tc) : zero(T)),
     CNn.Reaction("H+ + e → H(1s) + γ (ionising groups, α_1)", [:Hp, :e], [:H], (T, env) -> on(env) ? mg_arec(T, env.mgq.Tc) : zero(T))]
end
"""
Stage 2c: the star light of a zone IS its group fields. c u_g = c(E_g − E_cmb b_g(T_cmb)) [erg cm⁻² s⁻¹] of every group (floored at MGCUFL so that
ratios and logarithms exist where a group has underflowed or lies below its CMB part), and the dilution of group g against a blackbody of the colour temperature T_c,
W_g = c u_g/(4σT_c⁴ b_g(T_c)). No free-streaming floor (SLWFLUX), no cap at 1, no W = 1 rule for thermalised zones: the group energies carry
all of that (docs/POPIII_MULTIGROUP.md §18).
"""
const MGCUFL = 1e-60       # [erg cm⁻² s⁻¹] (photo-rates below 1e-66 s⁻¹; products of four such factors in the rate partials stay representable: 1e-240 gave 0·∞)
mg_cu(k, Ecmb) = (E = TX.MGEN[]; ntuple(g -> max(CNn.c * (E[g][k] - Ecmb * TX.MGBCMB_[g]), MGCUFL), TX.MGN[]))
mg_wg(cu, Tc) = (nrm = 4 * 5.670374419e-5 * Tc^4; ntuple(g -> (b = TX.mg_bfrac(g, Tc); b > 1e-180 ? G2O._qd(cu[g], nrm * b) : zero(cu[g])), length(cu)))

"""
(TX.MGSLQ, the callback `mgqt` of `sl_sweep!`) Zone k at the Newton iterate: q_g = c u_g/c u_g(start of step) of every group into M.qs[1:G, k] (M.egn: the
groups' E_g of the iterate), the LW factor of the record's step start into M.qs[G+1, k] (`mg_h2lw!` replaces it where the H₂ closure acts); returns q of
group MGSLG, which the sweep stores as the zone's grey q. The colour temperature is the sweep's own (T_src, w_q and T_r of the iterate, as the grey path).
"""
function sl_qt_mg(M, k, Tc)
    o = one(Tc); G = TX.MGN[]; qs = M.qs
    if !(k <= length(SLW0[]) && (r = SLW0[][k]) !== nothing && hasproperty(r, :cu0))
        @inbounds for g in 1:G+1; qs[g, k] = o; end
        return o
    end
    @inbounds for g in 1:G
        qs[g, k] = max(CNn.c * (M.egn[g][k] - r.Ecmb * TX.MGBCMB_[g]), MGCUFL) / r.cu0[g]
    end
    @inbounds qs[G+1, k] = r.s0 + zero(o)
    qs[TX.MGSLG[], k]
end

# ---------------------------------------------------------------- H₂ Lyman–Werner band: photon-conserving absorption (stage 2c)
"""
H₂ photodissociation and the LW group's energy, one closure (docs/POPIII_MULTIGROUP.md §18).

Rate per H₂ (Abel et al. 1997, k27): k = 1.1e8 F_ν(12.87 eV) f_sh s⁻¹, F_ν = 4πJ_ν in erg s⁻¹ cm⁻² Hz⁻¹; f_sh(N, T) of Wolcott-Green, Haiman &
Bryan (2011, eq. 10: Draine & Bertoldi 1996 eq. 37 with exponent 1.1), `fsh_h2`. Behind a column N the photons that dissociate are gone: the
unshielded band flux that an infinite column can still use is ∝ J(N) = ∫_N^∞ f_sh dN′ (closed form below), so the LW group's energy flux obeys
d ln(c u)/dN = −f_sh/J ≡ −σ_eff: every dissociation removes e_d = Δν_LW/(1.1e8 J(0)) of band energy (59 eV at b = 3 km s⁻¹; ⟨hν⟩/p_diss =
12.4 eV/0.15 = 83 eV from DB96's dissociation probability) and nothing else does. The group value c u(N) = c u(0) J(N)/J(0) is then the
attenuated field, and the rate from it is k = 1.1e8 (c u φ) f_sh J(0)/J, φ = B_ν(ν_LW, T_c)/∫_LW B_ν dν the in-group shape.

Returns (f_sh, J, J(0)) at column N [cm⁻²] and temperature T; J in cm⁻².
"""
function mg_fshJ(N, T)
    b5 = sqrt(CNn.kB * T / 1.6735575e-24) / 1e5; x = N / 5e14; a = 8.5e-4; u = sqrt(1 + x); E = exp(-a * u)
    f = 0.965 / (1 + x / b5)^1.1 + 0.035 / u * E
    J = 0.965 * 5e15 * b5 * (1 + x / b5)^-0.1 + 0.035 * 5e14 * (2 / a) * E
    J0 = 0.965 * 5e15 * b5 + 0.035 * 5e14 * (2 / a) * exp(-a)
    f, J, J0
end
const MGLWDNU = Ref(0.0)                                      # width of the LW group [Hz]
const MGKPH = Ref(0)                                          # the photosphere zone of the step's star-light setup
const MGX2 = Ref(Float64[])                                   # x_H₂ of the zones at the star-light setup (for the zones whose species are not unknowns)
const MGLWNU = 12.87 * EVERG / CNn.hP
"φ(T_c) Δν_LW: F_ν(12.87 eV) Δν_LW / (c u_LW) for the in-group shape B_ν(T_c); 0 where the band's Planck fraction underflows."
mg_lwshape(Tc) = (b = TX.mg_bfrac(MGROLE[].LW, Tc); (Tc > 0 && b > 1e-180) ? G2O._qd(π * _bnu(MGLWNU, Tc) * MGLWDNU[], 5.670374419e-5 * Tc^4 * b) : 0.0)
"Unshielded H₂ dissociation rate per molecule [s⁻¹] in the LW group field c u [erg cm⁻² s⁻¹] (AAZN97 k27 with F_ν = c u φ)."
mg_klw(cu, Tc) = 1.1e8 * cu * mg_lwshape(Tc) / MGLWDNU[]
"Effective H₂ column for the shielding: the static column from the photosphere, or the Sobolev column n_H₂ b/|du/dr| where that is smaller."
mg_ncol(Ns, n2, T, dudr) = dudr > 0 ? min(Ns, n2 * sqrt(CNn.kB * T / 1.6735575e-24) / dudr) : Ns
"Share of the absorbed LW energy that stays in the gas: 4.48 eV binding + 0.4 eV kinetic (AAZN97) of e_d; the rest is fluorescence (group MGH2T)."
mg_fgas(J0) = min(1.0, 4.88 * EVERG * 1.1e8 * J0 / MGLWDNU[])
"""
κ_H₂ [cm² g⁻¹] of the LW group in zones a..b from the state (dn, tn, rn, un, x_H₂ = x2): κ[k] the photon-conserving absorption (see `mg_fshJ`),
κg[k] its share that heats the gas, qs[G+1, k] the factor s = f_sh J(0)/J of the rate. Columns from the step's photosphere zone MGKPH outward at
the given state; only the local x2[k] carries derivatives (the column is non-local: an inexact Newton column, the residual is exact). `tcf(k)`:
the zone's colour temperature (0: no star light there).
"""
function mg_h2lw!(κ, κg, qs, a, b, dn, tn, rn, un, x2, tcf)
    T = eltype(κ); V = _val; G = TX.MGN[]; μH = CNn.μH(); X0 = MGX2[]
    fill!(κ, zero(T)); fill!(κg, zero(T))
    TX.MGH2G[] > 0 || return nothing
    isD = T <: ForwardDiff.Dual; isD1 = T <: ForwardDiff.Dual{<:Any,<:Any,1}
    k1 = max(MGKPH[] + 1, a)
    xk(j) = TX.cplz(j) ? x2[j] : (j <= length(X0) ? X0[j] : 0.0)
    # column of zone j seen from zone k: Duals of the zones in the band window of k (`winv`, as sl_sweep!'s Lyman-Werner columns), values elsewhere
    cj(j, k) = winv(dn[j], j, k) / μH * max(winv(xk(j), j, k), 0) * (winv(rn[j+1], j + 1, k) - winv(rn[j], j, k))
    P = get!(() -> zeros(b + 2), TX.SLQBUF, (:mg_h2P, b + 2))::Vector{Float64}; fill!(P, 0.0)       # value prefix sums of the columns of the zones below k
    @inbounds for j in k1:b; P[j+1] = P[j] + V(dn[j]) / μH * max(V(xk(j)), 0.0) * (V(rn[j+1]) - V(rn[j])); end
    @inbounds for k in k1:b
        Tc = tcf(k); Tc > 0 || continue
        ck = cj(k, k); N = P[k] + 0.5 * ck
        if isD                                       # the Duals inside the window (all zones for the 1-partial product)
            for j in (isD1 ? k1 : max(k1, k - 3)):k-1
                c_ = cj(j, k); N += c_ - V(c_)
            end
        end
        n2 = dn[k] / μH * max(xk(k), 0)
        f, J, J0 = mg_fshJ(mg_ncol(N, n2, tn[k], abs(un[k+1] - un[k]) / (rn[k+1] - rn[k])), tn[k])
        qs[G+1, k] = f * J0 / J
        κ[k] = max(xk(k), zero(xk(k))) / μH * (mg_lwshape(Tc) * f / J); κg[k] = κ[k] * mg_fgas(J0)
    end
    nothing
end
"(TX.MGH2FN) `mg_h2lw!` at the Newton iterate into the residual's group workspace M (slt: T_c/T_c,start of the iterate)."
function mg_h2iter!(M, dn, tn, rn, un, x2, slt, a, b)
    R = SLW0[]; V = ForwardDiff.value
    tcf = k -> (SLSTAR[] && k <= length(R) && (r = R[k]) !== nothing && hasproperty(r, :cu0) && r.lw) ? r.Tc0 * slt[k] : 0.0
    mg_h2lw!(M.kh2, M.kh2g, M.qs, a, b, dn, tn, rn, un, x2, tcf)
end
"`mg_h2lw!` at the colour temperatures slt (T_c/T_c,start) of a Float state: the LW factor s into qs[G+1, :] (zones with a star-light record), κ and κg per zone."
function mg_h2fill!(slt, qs, κ, κg)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); R = SLW0[]
    (TX.MGH2G[] > 0 && SLSTAR[] && length(TX.XSN[][2]) >= b) || return nothing
    tcf = k -> (SLSTAR[] && k <= length(R) && (r = R[k]) !== nothing && hasproperty(r, :cu0) && r.lw) ? r.Tc0 * Float64(slt[k]) : 0.0
    mg_h2lw!(κ, κg, qs, a, b, TM.v(:dn), TM.v(:tn), TM.v(:rn), TM.v(:un), TX.XSN[][2], tcf)
    nothing
end
"`mg_h2lw!` at the current Float state (formal solution, diagnostics): (κ, κg, s) per zone; the star-light sweep of the state supplies T_c."
function mg_h2now()
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); dn = TM.v(:dn); n = length(dn); G = TX.MGN[]
    κ = zeros(n); κg = zeros(n); qs = ones(G + 1, n)
    (TX.MGH2G[] > 0 && SLSTAR[] && SLIMP && length(TX.XSN[][2]) >= b) || return κ, κg, qs[G+1, :]
    S = (; q = ones(n), t = ones(n), fs = zeros(n), g = ones(n), ts = zeros(n), wq = zeros(n), dp = zeros(n)); Mq = (; egn = TX.MGEN[], qs)
    sl_sweep!(S, a, b, TM.v(:rn), dn, TM.v(:tn), TM.v(:un), TM.v(:ern), TM.v(:frn), TX.XSN[][1], TX.XSN[][2], TX.cplz; mgqt = (k, Tc) -> sl_qt_mg(Mq, k, Tc))
    mg_h2fill!(S.t, qs, κ, κg)
    κ, κg, qs[G+1, :]
end

"The photo-rate environment of a zone in the group fields cu at colour temperature Tc (ChemNet.photorates_mgq reads it)."
mg_env(cu, W, Tc, h1) = (; Tc, cu, W, ehz = MGEHZ[], h1)
"""
(stage 2c) The coupled zone's chemistry terms at the Newton iterate's group fields: qs[g] = c u_g/c u_g(record), s the LW factor f_sh J(0)/J,
tc = T_c/T_c,record, g the velocity gradient of the Lyα escape, dp the deep weight (as `cpl_tq`). Every photo-process of the network takes its rate from
the group that holds its photons: the photo list (H⁻, He, He⁺, Li, …) through ChemNet.photorates_mgq (the static field plus the group fields), the
three-level atom's Balmer-continuum and Lyα terms from the groups of those bands, H(1s) from the ionising groups, H₂ from the LW group with the
photon-conserving closure (`mg_fshJ`). The colour temperature enters in the number type of the iterate (the photo-rate tables are interpolated in ln T_c: exact slopes).
"""
function cpl_tqg(t, ρ, T, x1, x2, qs, s, tc, g = nothing, dp = nothing; need_cool::Bool = true, xn = nothing)
    s_ = t.slq; p3_0 = s_.p3; r = s_.rec; G = TX.MGN[]; ro = MGROLE[]; Tc = r.Tc0 * tc
    tenv0 = t.env
    tenv = if CNn.HM88 && hasproperty(tenv0, :hm) && s_.hm0 !== nothing       # HM88: the perturber densities of the iterate (as cpl_tq)
        nHt = ρ / CNn.μH(); merge(tenv0, (; hm = (nHt * (max(xneut(x1, x2, xn), 0) + s_.hm0[1] + x2), nHt * (x1 + s_.hm0[2]))))
    else
        tenv0
    end
    p3_1 = dp === nothing ? p3_0 : merge(p3_0, (; deep = dp))          # TL3DEEP: the iterate's deep weight
    RL = g === nothing ? p3_1.RL : lya_escape(ρ, x1, x2, g; xn)
    cu = ntuple(i -> qs[i] * r.cu0[i], G); W = mg_wg(cu, Tc)
    zf(x) = iszero(ForwardDiff.value(x)) ? 0.0 : x              # (an exact 0 stays a Float 0: tl3_heat has 0/0 partials otherwise, see cpl_tq)
    mgq = mg_env(cu, W, Tc, r.kHon)
    p3q = merge(p3_1, (; β = zf(W[ro.Ba] * beta2s_bb(Tc)), P12 = zf(_div_expm1((RL + LAM2S) * W[ro.Lya], E21 / (CNn.kB * Tc))), hν = balmer_mean(Tc), Tc, RL, mgq))
    env = merge(tenv, (; mgq, tl3 = p3q))
    ph_ = s_.ph::Vector{CNn.Photo}; pc = Ref{Any}(nothing)      # (the steady-state solves ask for the same T repeatedly: one evaluation per T)
    phq = Tv -> (c_ = pc[]; (c_ !== nothing && c_[1] === Tv) ? c_[2] : (v = _phvec(ph_, Tv, env); pc[] = (Tv, v); v))
    kLW = r.lw ? mg_klw(cu[ro.LW], Tc) * s : zero(s)
    extra = starlight_terms(s_.mm, 0.0, 0.0, 0.0, kLW)      # (H⁻ and the other photo-processes: the photo list, ChemNet.photorates_mgq; H(1s): mg_h1s_reactions)
    c0 = s_.c0
    cooling = need_cool ? ((ρ_, T_, d) -> (c0 === nothing ? zero(T_) : c0(ρ_, T_, d)) - tl3_heat(ρ_, T_, d, merge(p3q, (; ne = ρ_ / CNn.μH() * max(d[:Hp], zero(d[:Hp])))))[1]) : nothing
    (; env, phc = phq, cooling, extra, kLW)
end
"Has zone term set t a group-field star-light record (stage 2c)?"
mg_hasrec(t) = t !== nothing && hasproperty(t, :slq) && (r = t.slq.rec) !== nothing && hasproperty(r, :cu0)
@inline cpl_closeqg(vJ, v) = cpl_close(vJ, v) && (TX.AD_JV[] || all(abs(log(v[i] / vJ[i])) < CPLJLAG[1] for i in 5:length(v)-4) && abs(v[end-3] - vJ[end-3]) <= CPLJLAG[1] * max(v[end-3], vJ[end-3], 1e-3) &&
                             abs(log(v[end-2] / vJ[end-2])) < CPLJLAG[1] && abs(v[end-1] - vJ[end-1]) <= CPLJLAG[1] && abs(log(v[end] / vJ[end])) < CPLJLAG[1])
"(MGROUPS) local values and Jacobian at v = (ρ, T, x1, x2, q_1 … q_G, s, g, dp, t): AD in all variables (T_c enters the photo-rate tables and the record as a Dual)."
function cpl_vjqg(m, k, v, xn = nothing)
    y, J = cpl_vj(m, k, v, xn)
    if TX.SLIMPDBG[] > 0 && !all(isfinite, J)
        TX.SLIMPDBG[] -= 1; bad = [(i, j) for i in 1:8, j in 1:size(J, 2) if !isfinite(J[i, j])]
        println("SLIMPDBG cpl_vjqg zone $k v $v: non-finite J at $bad; y $y")
    end
    y, J
end
"(stage 2c) cpl_zone with the star-light state (q_1 … q_G, s) = qv, g, dp and t = T_c/T_c,start as local variables (see cpl_zone_q, cpl_tqg)."
function cpl_zone_qg(m::NetGas, k, ρ, T, a, b, qv, tc, g, dp, xn)
    V(x) = Float64(ForwardDiff.value(x)); n = length(qv)
    v = (V(ρ), V(T), V(a), V(b), ntuple(i -> V(qv[i]), n)..., V(g), V(dp), V(tc))
    needJ = any(x -> x isa ForwardDiff.Dual, (ρ, T, a, b, tc, g, dp)) || eltype(qv) <: ForwardDiff.Dual
    length(CPLCACHEQ[]) < k && (CPLCACHEQ[] = vcat(CPLCACHEQ[], Any[nothing for _ in length(CPLCACHEQ[])+1:k]))
    c = CPLCACHEQ[][k]
    if c === nothing || c[1] != v || c[5] != xn
        if needJ && !(c !== nothing && cpl_closex(cpl_closeqg, c, v, xn))
            y, J = cpl_vjqg(m, k, v, xn); c = (v, y, v, J, xn, xn)
        else
            y = cpl_g(m, k, [v...]; xn); c = c === nothing ? (v, y, nothing, nothing, xn, nothing) : (v, y, c[3], c[4], xn, c[6])
        end
        CPLCACHEQ[][k] = c
    elseif needJ && !cpl_closex(cpl_closeqg, c, v, xn)
        y, J = cpl_vjqg(m, k, v, xn); c = (v, y, v, J, xn, xn); CPLCACHEQ[][k] = c
    end
    needJ ? cpl_liftq(c[2], c[4], ρ, T, a, b, ntuple(i -> qv[i], n)..., g, dp, tc) : Tuple(c[2])
end
