# HM88 test: the network's LTE x_e with K_ion × Z_el/Z_w (what ChemNet's kinetics relax to with HM88=1) vs Z0 (HM88 Saha) and SCVH(XE);
# eos derivatives by ForwardDiff vs finite differences.
ENV["HM88"] = "1"
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
function net_lte_hm(ρ, T; it = 30)
    s = net_lte(ρ, T)
    for _ in 1:it
        Nneu = s.nH + s.nHe + s.nH2; Nion = s.ne
        f = CN.ZH_el(T) / (2CN.ZH_hm(T, Nneu, Nion))
        NH = ρ / CN.μH(); NHe = CN.yHe * NH; NLi = CN.yLi * NH
        K1 = exp(lnK([:H], [:Hp, :e], T)) * f; Km = exp(lnK([:H, :e], [:Hm], T)); K2 = exp(lnK([:H, :H], [:H2], T)); K2p = exp(lnK([:H, :Hp], [:H2p], T))
        S1 = exp(lnK([:He], [:Hep, :e], T)); S2 = exp(lnK([:Hep], [:Hepp, :e], T)); SLi = exp(lnK([:Li], [:Lip, :e], T))
        function parts(ne)
            b = 1 + K1 / ne + Km * ne; a = 2K2 + 2K2p * K1 / ne
            nH = 2NH / (b + sqrt(b^2 + 4a * NH)); nHp = K1 * nH / ne; r1 = S1 / ne; r2 = S2 / ne; he0 = NHe / (1 + r1 + r1 * r2); li0 = NLi / (1 + SLi / ne)
            (; nH, nHp, nHm = Km * nH * ne, nH2 = K2 * nH^2, nH2p = K2p * nH * nHp, nHe = he0, nHep = he0 * r1, nHepp = he0 * r1 * r2, nLi = li0, nLip = NLi - li0)
        end
        g(l) = (ne = exp(l); q = parts(ne); q.nHp + q.nH2p + q.nHep + 2q.nHepp + q.nLip - q.nHm - ne)
        lo, hi = log(1e-280), log(NH * 1.2)
        for _ in 1:300; m = 0.5(lo + hi); g(m) > 0 ? (lo = m) : (hi = m); hi - lo < 1e-14 && break; end
        ne = exp(0.5(lo + hi)); s = merge(parts(ne), (; ne, NH))
    end
    s
end
@printf("%6s %6s | %9s %9s %9s %9s | %8s %8s\n", "lrho", "lT", "xe_net", "xe_netHM", "xe_Z0HM", "xe_SCVHXE", "Zel/2", "Zw/2")
for (lρ, lT) in ((-6.42, 3.98), (-5.33, 4.26), (-4.67, 4.38), (-4.0, 4.4), (-3.08, 4.78), (-2.0, 4.9), (-1.0, 5.0), (-8.0, 4.5), (-10.0, 5.0))
    ρ = 10^lρ; T = 10^lT; s0 = net_lte(ρ, T); s = net_lte_hm(ρ, T)
    o = Z0.opacity_state(ρ, T; comp = Z0.Comp(XH, 1 - XH))
    xfix = NaN
    w = EQ.node(lρ, lT)[5]; ys, z1, z2, _ = EQ.saha_he(ρ, T); xfix = ys[1] + z1 + 2z2
    if w > 0
        sc = SCVH.state(ρ, T); X2, X1 = sc.xH2, sc.xH; Xp = max(1 - X2 - X1, 0.0) / 2; xfix = (1 - w) * ys[1] + w * 2Xp / (1 + 3X2 + X1) + z1 + 2z2
    end
    @printf("%6.2f %6.2f | %9.3e %9.3e %9.3e %9.3e | %8.2f %8.2f\n", lρ, lT, s0.ne / s0.NH, s.ne / s.NH, o.ne / o.nH, xfix, CN.ZH_el(T) / 2, CN.ZH_hm(T, s.nH + s.nHe, s.ne))
end
# eos AD vs FD with HM88 (a composition vector of a compiled toy network: use the species of `net_lte` through ChemNet.eos needs a Network)
net = CN.Network([:H, :Hp, :e, :He, :Hep, :Hepp, :H2], [:Hp, :H2, :Hep, :Hepp], CN.Reaction[])
ρ = 1e-5; T = 2e4; s = net_lte_hm(ρ, T); x = [s.nHp, s.nH2, s.nHep, s.nHepp] ./ s.NH
u(r, t) = CN.eos(r, t, x, net)[1]; p(r, t) = CN.eos(r, t, x, net)[2]
using ForwardDiff
for (nm, f) in (("u", u), ("p", p))
    dT = ForwardDiff.derivative(t -> f(ρ, t), T); dr = ForwardDiff.derivative(r -> f(r, T), ρ)
    fT = (f(ρ, T * (1 + 1e-6)) - f(ρ, T * (1 - 1e-6))) / (2e-6 * T); fr = (f(ρ * (1 + 1e-6), T) - f(ρ * (1 - 1e-6), T)) / (2e-6 * ρ)
    @printf("%s: d/dT AD %.8e FD %.8e | d/dρ AD %.8e FD %.8e\n", nm, dT, fT, dr, fr)
end
