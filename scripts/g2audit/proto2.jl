# PROTO2: G2Opacity.means on populations from a Saha with Hummer–Mihalas occupation-probability partition functions
# (Z0Physics.opacity_state: the H/He ionisation the network would have with HM88 partition functions) for T ≥ 10^3.7 K;
# below, the network's own LTE populations (molecules matter there, occupation probabilities do not).
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
include(joinpath(JD, "G2Opacity.jl")); const G2O = G2Opacity
const PROTO = (ρ, T, s) -> begin
    if T < 10^3.7
        n2 = G2O.n2_lte(T, s.nH, s.ne, s.nH + s.nHe + s.nH2)
        E, P, R = G2O.means(ρ, T; n1 = s.nH, n2, nH2 = s.nH2, nHe = s.nHe, nHep = s.nHep, nHepp = s.nHepp, nHp = s.nHp, ne = s.ne, nHm = s.nHm, nHmS = s.nHm, hbfP = true)
    else
        o = Z0.opacity_state(ρ, T; comp = Z0.Comp(XH, 1 - XH))
        nHm = o.pH[1] * o.ne * exp(lnK([:H, :e], [:Hm], T)) * CN.ZH_el(T) / 2      # H⁻ on the ground state (ChemNet K is per H atom with Z_H)
        E, P, R = G2O.means(ρ, T; n1 = o.nH0, n2 = o.pH[2], nH2 = 0.0, nHe = o.nHe0, nHep = o.nHe1, nHepp = o.nHe2, nHp = o.nHp, ne = o.ne, nHm, nHmS = nHm, hbfP = true)
    end
    (; κR = R, κP = P)
end
G2O.freq()
include(joinpath(@__DIR__, "grid.jl"))
