# Population III project, M2: one-zone free-fall tracks with the RHFC02 chemistry (PrimordialChem).
#   julia --project=env popiii_m2_onezone.jl [cooling] [k5]     cooling ∈ thin | ra04lines | ra04 (lines + CIE); k5 ∈ pss | ra04 | db (detailed balance with k4)
# Initial state: RHFC02 run P100 centre (Table 1): ρ/m_H = 2.6e6, T = 270 K, f_H2 = 5e-4 (fraction of H nuclei in H₂),
# f_e = 1e-10 (taken as x_e = x_H⁺ = n/n_H). Prints the track at every 0.5 dex of ρ/m_H. Log: lib/popiii/m2_onezone_*.log
using Printf
include(joinpath(@__DIR__, "PrimordialChem.jl")); const PC = PrimordialChem
cool = isempty(ARGS) ? "ra04" : ARGS[1]; k5s = length(ARGS) > 1 ? ARGS[2] : "pss"
cf = cool == "thin" ? PC.cool_H2_thin : cool == "ra04lines" ? PC.cool_H2_ra04 :
     (ρ, T, y) -> PC.cool_H2_ra04(ρ, T, y) + PC.cool_CIE_ra04(ρ, T, y)
rates = PC.rates_rhfc02(; k5 = k5s == "ra04" ? PC.k5_ra04 : k5s == "db" ? (T -> PC.k5_db(T)) : PC.k5_pss)
nmax = parse(Float64, get(ENV, "NMAX", "1e18"))
io = open(joinpath(@__DIR__, "..", "lib", "popiii", "m2_onezone_$(cool)_k5$(k5s).log"), "w")
pr(f, a...) = (s = Printf.format(Printf.Format(f), a...); print(s); print(io, s))
pr("one-zone free fall, RHFC02 rates (k5 = %s), cooling = %s\n", k5s, cool)
pr("%9s %9s %9s %10s %10s %10s %10s %10s\n", "log rho/mH", "log nH", "T [K]", "f_H2", "x_e", "x_H-", "Lambda/(p/t_ff)", "t [yr]")
nH0 = 2.6e6 * PC.mH / PC.μH; y0 = (1e-10, 0.0, 2.5e-4, 0.0)
next = Ref(log10(2.6e6))
PC.onezone(; nH0, T0 = 270.0, y0, rates, cooling = cf, nHmax = nmax, dlnρ = parse(Float64, get(ENV, "DLNRHO", "0.01")),
    out = function (ρ, T, y, t)
        lr = log10(ρ / PC.mH)
        if lr >= next[]
            p = PC.eos_rhfc02(ρ, T, y)[2]; tff = sqrt(3π / (32PC.G * ρ))
            pr("%9.2f %9.2f %9.1f %10.3e %10.3e %10.3e %10.3e %10.4e\n", lr, log10(ρ / PC.μH), T, 2y[3], PC.xe(y), y[2],
               cf(ρ, T, y) / (p / tff), t / 3.156e7)
            next[] += 0.5
        end
    end)
close(io)
