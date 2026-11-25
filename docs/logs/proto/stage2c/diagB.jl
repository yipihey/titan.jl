# Task B: old-MLT NaN at the core interfaces (Stage-1 code 700e60e), included by a hook after one step of the S16 restart (Float64 only).
# Prints per core interface: ∇ (raw), ∇_ad, c_P, δ, κ, U, ζ, the MLT_HS blend (lp, lph, sw), and the τ gate (convcoef's cvw) with its first NaN.
using Printf
# (included in PopIIIRun's module, where TX, PM, TF are bound)
dm = TX.MATERIAL[]; dh = PM.HOT[]
dc = TX.context(); dbase = dc.base; ngrs = dc.ngrs; ngre = dc.ngre
dn = dbase[:, 3]; tn = dbase[:, 5]; dr = dbase[:, 1]
println("ngrs ", ngrs, " ngre ", ngre, " T1 ", dh.T1, " T2 ", dh.T2, " ZONE0[] (stale, left by the last sweep) = ", TX.ZONE0[], " eqzone[ZONE0] = ", dm.eqzone[TX.ZONE0[]],
        " COLDMLT ", PM.COLDMLT[], " MLT_NOINFALL ", TX.MLT_NOINFALL[], " TAU_MLT ", TX.TAU_MLT[], " MLT_AD ", TX.MLT_AD[], " MLT_HS ", TX.MLT_HS[], " MLT_RCAP ", TX.MLT_RCAP[], " MLT_HSDISC ", TX.MLT_HSDISC[])
z0 = TX.ZONE0[]
println("\n--- zone opacity at (rho_k, T_k) with ZONE[] = k (correct) and with ZONE[] = ZONE0 (what convcoef's pre-sweep and its serial tail use)")
let
for k in ngrs+1:min(ngre, ngrs + 28)
    (dn[k] > 0 && tn[k] > 0) || continue
    TX.ZONE[] = k;  a = try TX.material_opac(dm, dn[k], tn[k])[1] catch e; NaN end
    TX.ZONE[] = z0; b = try TX.material_opac(dm, dn[k], tn[k])[1] catch e; NaN end
    flag = (isnan(a) ? " NaN(own zone)" : "") * (isnan(b) ? " NaN(ZONE0)" : "")
    @printf("k %3d  eq %d  rho %.4g T %.4g  chi(ZONE=k) %.5g  chi(ZONE=ZONE0) %.5g%s\n", k, dm.eqzone[k], dn[k], tn[k], a, b, flag)
end
end
TX.ZONE[] = z0
println("\n--- tau accumulation, outside in (convcoef's loop, ZONE = ZONE0 as in convcoef), first NaN")
let τ = 0.0, firstnan = 0
for k in ngre:-1:ngrs+1
    (dn[k] > 0 && tn[k] > 0) || continue
    χ = TX.material_opac(dm, dn[k], tn[k])[1]
    τ += χ * dn[k] * (dbase[k+1, 1] - dbase[k, 1])
    if isnan(τ) && firstnan == 0; firstnan = k; @printf("first NaN tau (outside in) at k = %d: chi = %g rho = %g T = %g dr = %g (ZONE=ZONE0=%d)\n", k, χ, dn[k], tn[k], dbase[k+1,1]-dbase[k,1], z0); end
end
@printf("tau at the centre %g\n", τ)
end
println("\n--- NaN origin: network zone ZONE0 = ", z0, " populations (NLOST) applied at the hot core state (rho_k, T_k): material_opac -> nlo_it -> nlo_ratio")
let s0 = PM.NLOST[][z0]
for k in (22, 24, 25, 26)
    ρ = dn[k]; T = tn[k]; CN = PM.CNn
    TX.ZONE[] = z0
    d, p = PM.nlo_comp(dm, z0, s0, ρ, T, 1.0, 1.0)
    nH = ρ / CN.μH()
    @printf("k %d rho %.4g T %.4g | zone %d d: H %.3g e %.3g Hp %.3g He %.3g H2 %.3g Hm %s | Keq(H- <-> H+e) %.4g\n", k, ρ, T, z0, d[:H], d[:e], d[:Hp], d[:He], d[:H2], haskey(d, :Hm) ? string(d[:Hm]) : "-", CN.Keq([:Hm], [:H, :e], T))
    n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
    n2 = PM.n2_tl3(T, n1, ne, np, p)
    bPR = PM.nlo_base(ρ, T)
    R = try PM.nlo_ratio(ρ, T, d, p) catch e; sprint(showerror, e) end
    @printf("   n1 %.4g ne %.4g np %.4g n2_tl3 %.4g | nlo_base (bP,bR) %s | nlo_ratio (R_E,R_P,R_R),ne = %s\n", n1, ne, np, n2, string(bPR), string(R))
    HE = (PM.G2O.ON && true)
    for hb in (false, true)
        y = try PM.nlte_means_fast(ρ, T; n1, n2 = n2, nH2 = nH * max(d[:H2], 0.0), nHe = nH * max(d[:He], 0.0), nHp = np, ne, nHm = n1 * ne / CN.Keq([:Hm], [:H, :e], T), hbfP = hb, Tw = T, nHep = nH * max(get(d, :Hep, 0.0), 0.0), nHepp = nH * max(get(d, :Hepp, 0.0), 0.0)) catch e; sprint(showerror, e) end
        println("   nlte_means_fast hbfP=", hb, " -> ", y)
    end
end
end

println("\n--- nlo_base(rho, T) (called by nlo_ratio, hot state, w = nlo_w(T) = 1): pieces at (31.32, 1.142e8) and the NaN boundary")
let CN = PM.CNn, PC_ = PM.PC, EQ_ = PM.EQ, G2 = PM.G2O
  for (ρ, T) in ((31.32, 1.142e8), (7.451, 1.385e8))
    nH = ρ / CN.μH(); w = PM.nlo_w(T); Rt, Pt = PM.OPAC[](ρ, T)
    y = EQ_.composition(ρ, min(T, 1e6)); xp = y[1]; xh2 = y[3]; n1e = nH * max(1 - xp - 2xh2, 0.0); nee = nH * xp; nHee = PC_.yHe * nH
    z1, z2 = EQ_.composition_he(ρ, min(T, 1e6)); nee2 = nH * (xp + z1 + 2z2); nHee2 = nH * max(PC_.yHe - z1 - z2, 0.0)
    Keq = CN.Keq([:Hm], [:H, :e], T); nHme = n1e * nee2 / Keq
    n2g = G2.n2_lte(T, n1e, nee2, n1e + nHee2 + nH * xh2)
    @printf("rho %.4g T %.4g w %.3g | OPAC (Rt,Pt) %.4g %.4g | EQ.composition xp %.6g xh2 %.4g | he z1 %.4g z2 %.4g | n1e %.4g nee %.4g nHee %.4g | Keq %.4g nHm %.4g | G2O.n2_lte %.4g\n", ρ, T, w, Rt, Pt, xp, xh2, z1, z2, n1e, nee2, nHee2, Keq, nHme, n2g)
    y2 = try PM.nlte_means_fast(ρ, T; n1 = n1e, n2 = n2g, nH2 = nH * xh2, nHe = nHee2, nHp = nH * xp, ne = nee2, nHm = nHme, hbfP = true, nHep = nH * z1, nHepp = nH * z2) catch e; sprint(showerror, e) end
    println("   nlte_means_fast(n1 = ", n1e, ", ...) = ", y2)
  end
  println("nlo_base NaN scan, T = 1.142e8: ", [(ρ, any(isnan, PM.nlo_base(ρ, 1.142e8))) for ρ in (2.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 31.32)])
  println("nlo_base NaN scan, rho = 31.32: ", [(T, any(isnan, PM.nlo_base(31.32, T))) for T in (1e6, 1e7, 3e7, 1e8, 1.142e8, 1.4e8)])
end
println("\n--- core interfaces: the recording pass (conv_diag_pass) + dh.thermo re-evaluated at (rho_i, T_i), ZONE = k")
dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
dcd = TX.conv_diag_pass(dt)
arad = 4dc.p.csigr / dc.p.cc
for k in ngrs+1:min(ngrs+22, ngre)
    ρi = dcd.rhoi[k]; Ti = dcd.Ti[k]; l = dcd.lmix[k]; g = dcd.g[k]; HP = dcd.HP[k]
    TX.ZONE[] = k; th = dh.thermo(ρi, Ti, k); κ = TX.material_opac(dm, ρi, Ti)[1]
    sup = dcd.nab[k] - dcd.nad_used[k]; sup < 0 && continue      # (the first run stopped at k = 26: sqrt of the negative sup of the stable layer)
    U = 3arad * dc.p.cc * Ti^3 / (dcd.cp_used[k] * ρi^2 * dcd.kap[k] * l^2) * sqrt(8HP / (g * dcd.del_used[k])); ζ = sqrt(U^2 + sup) - U
    @printf("k %2d rho_i %.4g T_i %.4g | lp_raw %.4g lph %.4g sw %.3g | nab_raw %.4g nab %.4g nad %.4g (dh.thermo %.4g) cp %.4g (%.4g) del %.4g (%.4g) kap %.4g (%.4g) | sup %.4g U %.3g zeta %.4g | fpot %.4g wtau(cvw) %.3g tau %.4g | fconv %.4g\n",
            k, ρi, Ti, dcd.lp_raw[k], dcd.lph[k], dcd.hs_w[k], dcd.nab_raw[k], dcd.nab[k], dcd.nad_used[k], th[1], dcd.cp_used[k], th[2], dcd.del_used[k], th[3], dcd.kap[k], κ, sup, U, ζ, dcd.fpot[k], dcd.cvw[k], dcd.tau[k], dcd.fconv[k])
end
println("fconv NaN interfaces: ", findall(isnan, dcd.fconv[ngrs+1:ngre]) .+ ngrs, "   cvw NaN interfaces (first/last): ", let v = findall(isnan, dcd.cvw[ngrs+1:ngre]) .+ ngrs; isempty(v) ? v : (first(v), last(v)) end)
