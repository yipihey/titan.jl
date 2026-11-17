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
    n == 3 ? [t.Hm, t.LW] : n == 4 ? [t.Hm, t.LW, t.H] : n == 6 ? [t.Hm, t.Ba, t.LW, t.H, t.He1] :
    n == 7 ? [t.Hm, t.Ba, t.LW, t.H, t.He1, t.He2] : n == 8 ? [t.Hm, t.Pa, t.Ba, t.LW, t.H, t.He1, t.He2] :
    error("MGROUPS=$n: default edge sets for 3, 4, 6, 7, 8 groups (or give the edges in eV)")
end
"Group containing photon energy e [eV]."
mg_gof(e) = searchsortedlast(TX.MGEDGE[], e) + 1
"Roles: the group of the H⁻ band (just above its threshold), of the Balmer continuum, of Lyα (E21), of the LW rate (12.87 eV), the first ionizing group."
const MGROLE = Ref{Any}(nothing)
function mg_roles!()
    t = mg_thresholds()
    MGROLE[] = (; Hm = mg_gof(t.Hm * 1.0001), Ba = mg_gof(t.Ba * 1.0001), Lya = mg_gof(E21 / EVERG), LW = mg_gof(12.87), H = mg_gof(t.H * 1.0001))
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
const MGTMIN = parse(Float64, get(ENV, "MGTMIN", "3000"))     # in-group shape temperature floor of the group means [K]

const MGEHZ = Ref(Float64[])                                  # interior group edges [Hz]
"The record of a zone without group data (MGROUPS off, or the group means failed): `nlo_it` then writes the default record."
const MGNOV = Float64[]
"""
Set-up of the material side (popiii_m2c mg_install!, after TX.mg_setup!): roles and Raman pairs, the groups of the G2 nodes, the line
table, per-zone diagnostics.
"""
function mg_init!()
    n = length(TM.v(:dn))
    MGEHZ[] = TX.MGEDGE[] .* (EVERG / CNn.hP); MGGIDX[] = G2O.group_index(MGEHZ[])
    mg_roles!()
    MGWREC[] = Any[nothing for _ in 1:n]; MGA1[] = zeros(n)
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
    x = clamp((lt - L.t0) / L.dt, 0.0, L.nt - 1 - 1e-9); i = floor(Int, x) + 1; a = x - (i - 1)
    y = clamp((log10(ρ) - L.r0) / L.dr, 0.0, L.nr - 1 - 1e-9); j = floor(Int, y) + 1; b = y - (j - 1)
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
r_R,g = χ_g/χ (the table's Rosseland mean has no lines); the absolute H(1s) (and, MGN2BF, n ≥ 2) bound–free absorption K_1 of the
groups at and above χ_H (LW band) where H bound–free is not in κ_E (case-B zones); the Milne emission coefficient
A_1 = σ̄_1s,g n_e n_p λ_e³/ρ of the same groups and zones (emission 4π A_1 e^{χ_H/kT} b_g(T) σT⁴/π: spontaneous + induced recombination
to the ground state in detailed balance with K_1; in LTE, n_1 = n_e n_p λ_e³ e^{χ_H/kT}, emission = absorption of a Planck field);
κ_Ram and the Raman pairs' κ_Ram C, κ_Ram ω (MGRad.jl MGRAMP; n_2s = n_2/4 of the three-level atom). A ratio that is not finite and
positive is 1 (grey).
"""
function mg_pack(ρ, T, Tw, n1, ne, np, n2, hbf, y)
    G = TX.MGN[]; o = zeros(TX.mg_nv()); r = MGROLE[]
    κE, κP, χ, k1, gE, gP, gχ, k2, kram, s1 = y
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
        o[j+4] = ((!hbf && g >= r.H && isfinite(k1[g])) ? k1[g] : 0.0) + ((MGN2BF && !hbf && g >= r.LW && isfinite(k2[g])) ? k2[g] : 0.0)
        o[j+5] = (MGMILNE && !hbf && g >= r.H && ok(s1[g])) ? s1[g] * ne * np * λ3 / ρ : 0.0
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
"Elementwise logarithmic slope of two records (0 where an entry is not positive in both, or the records differ in length)."
mg_slope(a, b, h) = length(a) == length(b) ? [(a[i] > 0 && b[i] > 0) ? log(a[i] / b[i]) / h : 0.0 for i in eachindex(a)] : zeros(length(b))
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
mg_wrec!(k, W, Wg, Tc) = (k <= length(MGWREC[]) && (MGWREC[][k] = (W, copy(Wg), Tc)); nothing)
"""
(per zone, set in starlight_setup where the group H(1s) photoionisation is on) the matching ground-state recombination of the network:
dx_H⁺/dt −= a1 x_e x_H⁺, a1 = n_H λ_e³ e^{χ_H/kT} × k_1s(B_ν(T)) (spontaneous Milne recombination: the photoionisation rate in J = B_ν(T) times the
Saha population) — the network is case B (no recombination to 1s), the groups now carry those photons, so with k_1s from the groups the
balance needs them (detailed balance: k_1s n_1s = a1 n_e n_p/n_H in LTE).
"""
const MGA1 = Ref(Float64[])
"MGN2BF=0: no H n ≥ 2 bound–free in the LW and ionising groups' κ_E of case-B zones (before 2a-fix)."
const MGN2BF = get(ENV, "MGN2BF", "1") == "1"
"MGMILNE=0: no H(1s) recombination emission in the ionising groups and no ground-state recombination in the network (stage-1 case B)."
const MGMILNE = get(ENV, "MGMILNE", "1") == "1"
"H(1s) photoionisation rate per atom of group g in J_ν = B_ν(T_c) (W = 1): 4π Σ σ_1s B_ν/(hν) w over the group's G2 nodes."
function mg_kH1s(g, Tc)
    F = G2O.fastgrid(); s = 0.0; gi_ = MGGIDX[]
    @inbounds for i in eachindex(F.ν)
        gi_[i] == g || continue
        σ = F.SH[1, i]; σ > 0 || continue
        s += σ * _bnu(F.ν[i], Tc) / (CNn.hP * F.ν[i]) * F.w[i]
    end
    4π * s
end
"""
Group dilutions W_g of zone k for colour temperature Tc (the grey star-light rule per group): c u_g = c(E_g − E_cmb b_g(T_cmb)) (SLWFLUX:
at least the group's outward flux), W_g = c u_g/(4σT_c⁴ b_g(T_c)), capped at 1 where the grey W is (`cap`). Thermalised zones (wq ≥ 1): 1.
"""
function mg_wg(k, Tc, Ecmb, cap::Bool, therm::Bool)
    G = TX.MGN[]; therm && return ones(G)
    E = TX.MGEN[]; F = TX.MGFN[]; W = zeros(G)
    for g in 1:G
        Ec = Ecmb * TX.MGBCMB_[g]
        cu = SLWFLUX ? max(CNn.c * (E[g][k] - Ec), 0.5 * (max(F[g][k], 0.0) + max(F[g][k+1], 0.0))) : CNn.c * (E[g][k] - Ec)
        w = max(cu, 0.0) / (4 * 5.670374419e-5 * Tc^4 * TX.mg_bfrac(g, Tc))
        W[g] = cap ? min(w, 1.0) : w
    end
    W
end
"J_ν = W_g B_ν(T_c) piecewise in the groups."
mg_J(Wg, Tc) = ν -> Wg[searchsortedlast(TX.MGEDGE[], CNn.hP * ν / EVERG) + 1] * _bnu(ν, Tc)
"H(1s) photoionisation rate of the ionising groups for dilutions Wg."
mg_kH(Wg, Tc) = (gH = MGROLE[].H; s = 0.0; for g in gH:TX.MGN[]; Wg[g] > 0 && (s += Wg[g] * mg_kH1s(g, Tc)); end; s)

"""
(SLIMP with MGROUPS) (q, t) of zone k at the iterate: as `sl_qt` with the dilution of group MGSLG (the H⁻ band) from its E_g, F_g, F_g+1,
and T_r (the local part of T_c, SLIMP=2) from the grey E_r.
"""
function sl_qt_mg(k, E, Eg, Fk, Fk1)
    o = one(E) + zero(Eg) + zero(Fk) + zero(Fk1)
    (k <= length(SLW0[]) && (r = SLW0[][k]) !== nothing && hasproperty(r, :Wg)) || return (o, o)
    arad = 7.565723e-15; g = TX.MGSLG[]; Ec = r.Ecmb * TX.MGBCMB_[g]
    cu = SLWFLUX ? max(CNn.c * (Eg - Ec), 0.5 * (max(Fk, zero(Fk)) + max(Fk1, zero(Fk1)))) : CNn.c * (Eg - Ec) + zero(Fk) + zero(Fk1)
    Tr = (max(E, r.Ecmb + zero(E)) / arad)^0.25 + zero(Fk) + zero(Fk1) + zero(Eg)
    if r.mode === :blend
        Tc = SLIMPN >= 2 ? r.Tsrc^(1 - r.wq) * Tr^r.wq : r.Tc0 + zero(Tr)
        W = r.wq >= 1 ? o : min(max(cu, zero(cu)) / (4 * 5.670374419e-5 * Tc^4 * TX.mg_bfrac(g, Tc)), one(cu))
    elseif r.mode === :te
        Tc = r.Tc0 + zero(Tr); W = max(cu, zero(cu)) / (4 * 5.670374419e-5 * r.Tc0^4 * TX.mg_bfrac(g, r.Tc0))
    else
        Tc = SLIMPN >= 2 ? Tr : r.Tc0 + zero(Tr); W = o
    end
    (max(W / r.W0, 1e-6 + zero(W)), Tc / r.Tc0)
end
