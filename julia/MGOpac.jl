# ------------------------------------------------------------------ MGROUPS: material side (group opacity ratios, group photo-rates)
# Included into PopIIIMaterial. Design: docs/POPIII_MULTIGROUP.md. Nothing here runs unless MGROUPS is set (MGON false: no call).
const MGON = haskey(ENV, "MGROUPS") && ENV["MGROUPS"] != "0"
"Group edges of the code's own thresholds [eV]: H⁻ photodetachment (λ₀ = 1.6419 μm, ContinuumOpacity.σ_Hm_bf), H n = 3, 2 (χ_H/n²), the
lower end of the Lyman–Werner band (11.2 eV, as RADGROUPS; Draine & Bertoldi 1996), H, He I, He II ionisation (Z0Physics χ_H, χ_He1, χ_He2:
the thresholds of the G2 continuum that supplies the group opacities)."
mg_thresholds() = (Hm = COp.h * COp.c * (1e4 / 1.6419) / EVERG, Pa = Z0Physics.χH / 9 / EVERG, Ba = Z0Physics.χH / 4 / EVERG, LW = 11.2,
                   H = Z0Physics.χH / EVERG, He1 = Z0Physics.χHe1 / EVERG, He2 = Z0Physics.χHe2 / EVERG)
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
end
const MGGIDX = Ref{Vector{Int}}(Int[])
const MGTMIN = parse(Float64, get(ENV, "MGTMIN", "3000"))     # in-group shape temperature floor of the group means [K]

"""
Per-zone group ratios at the start of a hydro step (`TX.MGRE/MGRP/MGRR`, `TX.MGK1`): network zones with a composition record (NLOST:
the last chemistry's species and three-level parameters) take G2Opacity.group_fast on that composition (κ_E weighted with B_ν(T_c) as
NLOPACJ, H bound–free in κ_E/κ_P only in deep zones, NLOHBFDEEP), ratio = group mean / grey mean of the same call; the absolute H(1s)
absorption of the groups at and above χ_H where H bound–free is not in κ_E (case B zones). Table zones and zones without a record: 1, 0
(the grey field split by Planck fractions: the diffusion limit of the optically thick interior).
"""
function mg_ratios!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); dn, tn = TM.v(:dn), TM.v(:tn); G = TX.MGN[]
    isempty(MGGIDX[]) && (MGGIDX[] = G2O.group_index(TX.MGEDGE[] .* (EVERG / CNn.hP)))
    MGROLE[] === nothing && mg_roles!()
    gH = MGROLE[].H
    for g in 1:G; fill!(TX.MGRE[][g], 1.0); fill!(TX.MGRP[][g], 1.0); fill!(TX.MGRR[][g], 1.0); fill!(TX.MGK1[][g], 0.0); end
    st = NLOST[]
    for k in a:b
        (m.eqzone[k] || k > length(st) || st[k] === nothing) && continue
        d, p = st[k].d, st[k].p
        ρ = dn[k]; T = tn[k]; nH = ρ / CNn.μH(); nHe = PC.yHe * nH
        n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
        CNn.HECHEM && (nHe = nH * max(d[:He], 0.0))
        nHmS = n1 * ne / CNn.Keq([:Hm], [:H, :e], T); nHm = haskey(d, :Hm) ? nH * max(d[:Hm], 0.0) : nHmS
        Tw = NLOPACJ[] && hasproperty(p, :Tc) ? p.Tc : T; hbf = hbf_deep(p)
        y = try
            G2O.group_fast(ρ, T, MGGIDX[], G; n1, n2 = n2_tl3(T, n1, ne, np, p), nH2 = nH * max(d[:H2], 0.0), nHe, nHp = np, ne, nHm, nHmS,
                           hbfP = hbf, Tw, nHep = nH * max(get(d, :Hep, 0.0), 0.0), nHepp = nH * max(get(d, :Hepp, 0.0), 0.0), Tmin = MGTMIN)
        catch
            nothing
        end
        y === nothing && continue
        κE, κP, χ, k1, gE, gP, gχ = y
        for g in 1:G
            rE = gE > 0 ? κE[g] / gE : 1.0; rP = gP > 0 ? κP[g] / gP : 1.0; rR = gχ > 0 && χ[g] > 0 ? χ[g] / gχ : 1.0
            TX.MGRE[][g][k] = isfinite(rE) && rE > 0 ? rE : 1.0
            TX.MGRP[][g][k] = isfinite(rP) && rP > 0 ? rP : 1.0
            TX.MGRR[][g][k] = isfinite(rR) && rR > 0 ? rR : 1.0
            (!hbf && g >= gH && isfinite(k1[g])) && (TX.MGK1[][g][k] = k1[g])
        end
    end
    nothing
end

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
