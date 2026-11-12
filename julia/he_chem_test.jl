using Printf, ForwardDiff
dir = ARGS[1]; mode = ARGS[2]
for f in ("PrimordialChem", "ContinuumOpacity", "ChemNet", "ChemRates"); include(joinpath(dir, f * ".jl")); end
const CN = ChemNet
net, photos = ChemRates.modern()
println("free: ", net.free, "  HECHEM ", isdefined(CN, :HECHEM) ? CN.HECHEM : "n/a")
function main(mode)
if mode == "bits"
    for (nH, T, xp) in ((1e10, 3000.0, 1e-6), (1e13, 6000.0, 1e-3), (1e11, 9000.0, 1e-2), (1e9, 4000.0, 1e-4), (1e15, 2500.0, 1e-9))
        ρ = nH * CN.μH(); x0 = zeros(length(net.free)); x0[1] = xp; x0[3] = 0.1
        u = CN.eos(ρ, T, x0, net)[1]
        for env in ((Trad = T, J = ν -> 2 * 6.62607015e-27 * ν^3 / 2.99792458e10^2 / expm1(6.62607015e-27 * ν / (1.380649e-16 * T)), static = true), (Trad = 0.0, J = ν -> 0.0, static = true, caseB = true))
            x, T1 = CN.step(x0, ρ, u, T, 1e6, net; env, photos)
            println(join([bitstring(v)[1:64] |> s -> string(parse(UInt64, s; base = 2), base = 16) for v in vcat(x, T1)], " "))
        end
    end
else
    iHp = findfirst(==(:Hep), net.free); iHpp = findfirst(==(:Hepp), net.free); yHe = CN.yHe
    me = 9.1093837e-28; h = 6.62607015e-27; kB = 1.380649e-16
    for (nH, T) in ((1e12, 2e4), (1e12, 4e4), (1e14, 6e4), (1e10, 1e5), (1e13, 1.5e5))
        ρ = nH * CN.μH()
        for (lab, env) in (("J=B(T)", (Trad = T, J = ν -> 2 * 6.62607015e-27 * ν^3 / 2.99792458e10^2 / expm1(6.62607015e-27 * ν / (1.380649e-16 * T)), static = true)), ("no light", (Trad = 0.0, J = ν -> 0.0, static = true, caseB = false)))
            x0 = zeros(length(net.free)); x0[1] = 0.9
            u = CN.eos(ρ, T, x0, net)[1]
            x, T1 = CN.step(x0, ρ, u, T, 1e12, net; env, photos, Tfix = T)
            d = CN.full(x, net); ne = d[:e] * nH
            # Saha with the network's own partition functions (checks the network reaches its detailed-balance limit)
            f = (2π * me * kB * T / h^2)^1.5
            s1 = 2f * CN.ZHep_el(T) / CN.ZHe_el(T) * exp(-CN.χHe1 / (kB * T)) / ne
            s2 = 2f * 1.0 / CN.ZHep_el(T) * exp(-CN.χHe2 / (kB * T)) / ne
            @printf("nH %.0e T %.1e %-8s: He %.4e He+ %.4e He++ %.4e  ne/nH %.4e | (He+/He)/Saha %.6f (He++/He+)/Saha %.6f  sum %.3e\n",
                    nH, T, lab, d[:He], d[:Hep], d[:Hepp], d[:e], d[:Hep] / d[:He] / s1, d[:Hepp] / d[:Hep] / s2, d[:He] + d[:Hep] + d[:Hepp] - yHe)
        end
    end
    # energy: fixed u, start hot and neutral
    nH = 1e12; ρ = nH * CN.μH(); T0 = 1.5e5; x0 = zeros(length(net.free)); x0[1] = 1.0
    u0 = CN.eos(ρ, T0, x0, net)[1]
    x, T1 = CN.step(x0, ρ, u0, T0, 1e4, net; env = (Trad = 0.0, J = ν -> 0.0, static = true, caseB = false), photos)
    d = CN.full(x, net)
    @printf("energy: T %.4e -> %.4e, He+ %.3e He++ %.3e, u0 %.12e u1 %.12e rel %.2e\n", T0, T1, d[:Hep], d[:Hepp], u0, CN.eos(ρ, T1, x, net)[1], CN.eos(ρ, T1, x, net)[1] / u0 - 1)
end
end
main(mode)
