# The zdiagB cliff zones: both materials at the zones' (ρ, T); the network also with the zones' own (non-LTE) composition.
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
include(joinpath(JD, "G2Opacity.jl")); const G2O = G2Opacity
# frames exported by python (frames.py dump): columns zone logrho logT logTrad xHp xH2 logxHm logtau
rows = [parse.(Float64, split(l)) for l in readlines(ARGS[1]) if !isempty(strip(l))]
g2(ρ, T, s; nHm = s.nHm, nHmS = s.nHm, nH = s.nH, ne = s.ne, nHp = s.nHp, nH2 = s.nH2) =
    G2O.means(ρ, T; n1 = nH, n2 = G2O.n2_lte(T, nH, ne, nH + s.nHe + nH2), nH2, nHe = s.nHe, nHep = s.nHep, nHepp = s.nHepp, nHp, ne, nHm, nHmS, hbfP = true)
@printf("%4s %7s %6s %6s %6s | %9s %9s %9s | %9s %9s | %9s %9s %9s | %9s %9s %9s %9s | %9s %9s %9s\n", "zone", "lrho", "lT", "lTr", "ltau",
        "xe_tab", "xe_tabXE", "xe_LTE", "cP_tab", "cP_net", "kR_tab", "kR_net", "kR_G2", "kP_tab", "kP_net", "kP_G2", "kP_G2noLy", "xe_zone", "kR_zone", "kP_zone")
for r in rows
    k, lρ, lT, lTr, xp, x2, lxm, lτ = r; ρ = 10^lρ; T = 10^lT
    s = net_lte(ρ, T); tt = thermo(tab_pu, ρ, T); tn = thermo(net_pu, ρ, T)
    y = EQ.composition(ρ, T); z = EQ.composition_he(ρ, T); xet = y[1] + z[1] + 2z[2]
    ys, z1, z2, _ = EQ.saha_he(ρ, T); w = EQ.node(lρ, lT)[5]; xfix = ys[1] + z1 + 2z2
    if w > 0
        sc = SCVH.state(ρ, T); X2, X1 = sc.xH2, sc.xH; Xp = max(1 - X2 - X1, 0.0) / 2
        xfix = (1 - w) * ys[1] + w * Xp / (2X2 + X1 + Xp) + z1 + 2z2
    end
    κt = tab_kap(ρ, T); m = net_means(ρ, T, s); gm = g2(ρ, T, s)
    gnl = G2O.means(ρ, T; n1 = s.nH, n2 = G2O.n2_lte(T, s.nH, s.ne, s.nH + s.nHe + s.nH2), nH2 = s.nH2, nHe = s.nHe, nHep = s.nHep, nHepp = s.nHepp, nHp = s.nHp,
                    ne = s.ne, nHm = s.nHm, nHmS = s.nHm, hbfP = false)
    # the zone's own composition: H⁺, H₂, H⁻ from the frame; He and Li as LTE; H atoms by conservation
    NH = s.NH; nHp = xp * NH; nH2 = x2 * NH; nHm = 10^lxm * NH; nH = max(NH - nHp - 2nH2 - nHm, 0.0)
    ne = nHp + s.nHep + 2s.nHepp + s.nLip - nHm
    sz = merge(s, (; nH, nHp, nH2, nHm, ne))
    nHmS = nH * ne * exp(lnK([:H, :e], [:Hm], T))
    mz = g2(ρ, T, sz; nHm, nHmS, nH, ne, nHp, nH2)
    @printf("%4d %7.3f %6.3f %6.3f %6.2f | %9.3e %9.3e %9.3e | %9.3e %9.3e | %9.3e %9.3e %9.3e | %9.3e %9.3e %9.3e %9.3e | %9.3e %9.3e %9.3e\n", k, lρ, lT, lTr, lτ,
            xet, xfix, s.ne / NH, tt.cP, tn.cP, κt[1], m.κR, gm[3], κt[2], m.κP, gm[2], gnl[2], ne / NH, mz[3], mz[2])
end
