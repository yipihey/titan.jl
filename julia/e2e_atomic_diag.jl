# E2E-A diagnostic: contributions of the photo pairs to the atomic `rad` cooling term in the dense core of the z = 19 collapse (netatomicli, J = CMB only).
#   TITAN_NOLIB=1 julia --project=env e2e_atomic_diag.jl
using Printf, ForwardDiff
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics", "TitanMaterials", "ProblemSetup", "PrimordialChem", "SCVH", "PopIIIEqEOS",
          "ContinuumOpacity", "ChemNet", "ChemRates", "AtomicHHe")
    include(joinpath(@__DIR__, f * ".jl"))
end
const PC = PrimordialChem; const CN = ChemNet
netw, phot = AtomicHHe.atomic(; li = true)
tcmb = 54.5
Bν_(ν, T) = 2CN.hP * ν^3 / CN.c^2 / expm1(CN.hP * ν / (CN.kB * T))
env = (Trad = tcmb, J = ν -> Bν_(ν, tcmb), static = true)
eqd(ρ, T) = (y = PopIIIEqEOS.composition(ρ, T); Dict(:H => 1 - y[1] - 2y[3], :Hp => y[1], :H2 => y[3], :He => CN.yHe, :Hep => 0.0, :Hepp => 0.0, :e => y[1], :Hm => 0.0, :H2p => 0.0))
for (nH, T) in ((1e14, 1900.0), (1e16, 2200.0), (3e16, 2310.0), (1e17, 2500.0), (1e18, 2500.0), (1e19, 3000.0))
    ρ = nH * CN.μH(); d = eqd(ρ, T)
    # full network composition at the Saha/equilibrium values: H⁻ from Saha with the electrons
    nHe = nH; ne = d[:e] * nH
    println("nH = $nH T = $T  x_e = $(d[:e]) f_H2 = $(2d[:H2])")
    tot = 0.0
    thermal = 1.5 * CN.kB * T * nH * (1 + d[:H2] * 0 + CN.yHe)
    for p in phot
        RE, AE, AX = CN.photoenergy_t(p, T, env)
        xion = get(d, p.P[1], 0.0); xX = get(d, p.R, 0.0)
        c = nH^2 * xion * get(d, :e, 0.0) * 1.0 * AE - nH * xX * RE
        tot += c
        abs(c) > 1e-12 * thermal && @printf("   %-22s emission-absorption %+.3e erg/cm3/s (AE %.2e RE %.2e, x_ion %.2e x_R %.2e)\n", p.name, c, AE, RE, xion, xX)
    end
    @printf("   sum %.3e, thermal energy density %.3e erg/cm3 -> t_cool = %.3e s\n", tot, thermal, thermal / max(abs(tot), 1e-300))
end
