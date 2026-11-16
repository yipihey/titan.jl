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
"Network EOS with HM88 (as ChemNet.eos with HM88=1) of the HM88-LTE composition: (p, u)."
function net_pu_hm(ρ, T)
    s = net_lte_hm(ρ, T); u, p = net_eos(ρ, T, s)
    eex, pex = CN.hm_exc(T, s.nH + s.nHe + s.nH2, s.ne)
    (p + s.nH * pex * CN.kB * T, u + s.nH * eex / ρ)
end
const QE = 4.80320471e-10
"Non-ideality parameters at (ρ, T), LTE (table) composition: Γ_e = e²/(a_e kT), a_e = (3/(4π n_e))^(1/3); D = n_e λ_e³/2; P_HM = −ln w_1 (HM88 4.71)."
function nonideal(ρ, T)
    y = EQ.composition(ρ, T); z = EQ.composition_he(ρ, T); nH = ρ / CN.μH()
    ne = nH * (y[1] + z[1] + 2z[2]); nneu = nH * (max(1 - y[1] - 2y[3], 0.0) + y[3] + max(CN.yHe - z[1] - z[2], 0.0))
    kT = CN.kB * T
    Γ = ne > 0 ? QE^2 / (kT * (3 / (4π * ne))^(1 / 3)) : 0.0
    λ = CN.hP / sqrt(2π * CN.me * kT); D = ne * λ^3 / 2
    PHM = -log(CN.w_hm(1, 1, nneu, ne))
    Pneu = (4π / 3) * nneu * (2 * 5.29177211e-9)^3      # HM88 neutral (hard-sphere) term of the ground state
    (; Γ, D, PHM, Pneu, ne, nneu)
end
