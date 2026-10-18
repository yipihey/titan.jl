# Tests of NLTEH.jl.  julia --project=env nlte_test.jl
# (1) atomic data against standard values; (2) with J = B(T) the statistical equilibrium must return LTE (b = 1) for
# any escape probabilities and collision scaling; (3) the non-LTE opacity with LTE populations has S_th = B(T).
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Printf
const NH = NLTEH; const SP = Spectrum; const Z0 = Z0Physics
@printf("f(Lyα) %.4f (0.4162)  f(Lyβ) %.4f (0.0791)  f(Hα) %.4f (0.6407)  f(Hβ) %.4f (0.1193)  f(Paα) %.4f (0.8421)\n",
        NH.fosc(1, 2), NH.fosc(1, 3), NH.fosc(2, 3), NH.fosc(2, 4), NH.fosc(3, 4))
@printf("σ(n=1, threshold) %.3e (6.30e-18)  σ(n=2, threshold) %.3e (1.4e-17)\n", NH.σbf(1, Z0.χH / Z0.h * 1.0001), NH.σbf(2, Z0.χH / 4Z0.h * 1.0001))
@printf("A(Lyα) %.3e (4.70e8, l-averaged)  A(Hα) %.3e (4.41e7)\n",
        let νl = (NH.En[2] - NH.En[1]) / Z0.h; 2Z0.h * νl^3 / Z0.c^2 * NH.gn[1] / NH.gn[2] * 4π * NH.πe2mc * NH.fosc(1, 2) / (Z0.h * νl) end,
        let νl = (NH.En[3] - NH.En[2]) / Z0.h; 2Z0.h * νl^3 / Z0.c^2 * NH.gn[2] / NH.gn[3] * 4π * NH.πe2mc * NH.fosc(2, 3) / (Z0.h * νl) end)
ν = SP.grid(; λmax_μm = 30.0, n = 250, nedge = 12)
w = [i == 1 ? 0.5(ν[2] - ν[1]) : i == length(ν) ? 0.5(ν[end] - ν[end-1]) : 0.5(ν[i+1] - ν[i-1]) for i in eachindex(ν)]
for (ρ, T) in ((1e-15, 6600.0), (1e-12, 9000.0), (1e-16, 4000.0))
    s = Z0.saha(ρ, T); pH, _ = Z0.level_pops(T, s)
    J = SP.Bν.(ν, T)
    worst = 0.0
    for (β0, cs) in ((1.0, 1.0), (1e-3, 0.01), (0.3, 100.0))
        n, np, ne, _, _ = NH.solve_point(T, s.nH, s.ne - s.nHp, 1.6e8, 8.4e15, -0.5, ν, w, J, 0.0, pH .* 1.3, s.nHp * 0.7; cscale = cs)
        b = [n[i] / (np * ne * NH.Φ(i, T)) for i in 1:NH.NL]
        worst = max(worst, maximum(abs.(b .- 1)), abs(np / s.nHp - 1))
    end
    @printf("ρ=%.0e T=%.0f: J = B ⇒ max |b − 1|, |n_p/n_p,LTE − 1| = %.2e\n", ρ, T, worst)
end
let ρ = 1e-15, T = 6600.0, s = Z0.saha(ρ, T), pH = Z0.level_pops(T, s)[1]
    atm = (r = [1e15], ρ = [ρ], T = [T]); op = SP.Opac(atm)
    no = NH.NOpac(op, Any[(n = pH, np = s.nHp, ne = s.ne)])
    worst = maximum(abs(NH.SP.coeffs3(no, 1, x)[3] / SP.Bν(x, T) - 1) for x in ν if SP.Bν(x, T) > 1e-250)
    @printf("LTE populations in NOpac: max |S_th/B − 1| = %.2e\n", worst)
end
