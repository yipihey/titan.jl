# atomic_cooling_test.jl — checks of julia/AtomicCooling.jl: the CIE cooling curve of
# primordial gas (shape and magnitude), the He fractions, and the ChemNet adapter.
#   julia --project=env atomic_cooling_test.jl
using Printf
include(joinpath(@__DIR__, "AtomicCooling.jl"))
const AC = AtomicCooling
pass = Ref(0); fail = Ref(0)
chk(name, ok) = (ok ? (pass[] += 1) : (fail[] += 1); @printf("%-70s %s\n", name, ok ? "PASS" : "FAIL"))

# CIE composition of H (k1/k2) and He (k3/k4, k5/k6) at n_H = 1
function cie(T)
    r = AC.k1(T) / AC.k2(T); xHII = r / (1 + r)
    f1, f2, f3 = AC.he_fractions(T, 1.0)
    nHe = AC.yHe
    (nHI = 1 - xHII, nHII = xHII, nHeI = f1 * nHe, nHeII = f2 * nHe, nHeIII = f3 * nHe, ne = xHII + nHe * (f2 + 2 * f3))
end
println("CIE cooling curve, Λ/n_H² [erg cm³ s⁻¹] (thermal-energy convention):")
println("  log T   H_exc    H_ci     H_rec    He_exc   He_ci    He_rec   brem     total    He/total  x_HII   x_HeII  x_HeIII  ceHI(T)")
Ts = 10 .^ (4.0:0.25:7.0)
tot = Float64[]; Hs = Float64[]; Hes = Float64[]
for T in Ts
    c = cie(T); ch = AC.channels(T; c...)
    push!(tot, ch.total); push!(Hs, ch.H); push!(Hes, ch.He)
    @printf("  %5.2f  %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e  %6.3f   %.3f  %.3f  %.3f  %8.1e\n", log10(T), ch.H_exc, ch.H_ci, ch.H_rec, ch.He_exc, ch.He_ci, ch.He_rec, ch.brem, ch.total, ch.He / ch.total, c.nHII, c.nHeII / AC.yHe, c.nHeIII / AC.yHe, AC.ceHI(T))
end
iH = argmax(Hs .* (Ts .< 6e4))
chk("H (Lyα) peak between 1.5e4 and 3e4 K", 1.5e4 < Ts[iH] < 3e4)
chk("H peak Λ/n_H² within [3e-23, 3e-22] (Dalgarno & McCray 1972 range)", 3e-23 < Hs[iH] < 3e-22)
share = Hes ./ tot
s5 = share[findfirst(isapprox.(log10.(Ts), 5.0; atol = 1e-6))]
chk("He share at 1e5 K > 0.3 and above its values at 3e4 and 1e6 K (He⁺/He²⁺ excitation)", s5 > 0.3 && s5 > share[findfirst(isapprox.(log10.(Ts), 4.5; atol = 1e-6))] && s5 > share[findfirst(isapprox.(log10.(Ts), 6.0; atol = 1e-6))])
chk("bremsstrahlung dominates at 1e7 K", (AC.brem(1e7) * (cie(1e7).nHII + cie(1e7).nHeII) + AC.brem_Z(1e7, 2) * cie(1e7).nHeIII) * cie(1e7).ne / tot[end] > 0.8)
f = AC.he_fractions(5e4, 1.0); chk("He at 5e4 K (CIE): mostly He⁺", f[2] > 0.5)
f = AC.he_fractions(2e5, 1.0); chk("He at 2e5 K (CIE): mostly He²⁺", f[3] > 0.9)
f = AC.he_fractions(1e4, 1.0); chk("He at 1e4 K (CIE): neutral", f[1] > 0.99)
fs = AC.he_fractions(2e4, 1e10; mode = :saha); chk("Saha at 2e4 K, n_e = 1e10: He partly ionised (> CIE)", fs[2] > AC.he_fractions(2e4, 1e10)[2])
# the ChemNet adapter: ionisation energy in u → no ciHI, recombination removes χ + kinetic
d = Dict(:H => 0.5, :Hp => 0.5, :e => 0.5, :He => AC.yHe)
T = 2e4; nH = 1.0
fh = AC.he_fractions(T, 0.5)
c = AC.channels(T; nHI = 0.5, nHII = 0.5, nHeI = fh[1] * AC.yHe, nHeII = fh[2] * AC.yHe, nHeIII = fh[3] * AC.yHe, ne = 0.5)
Λu = AC.cooling_atomic(nH, T, d; ionization_in_u = true); Λt = AC.cooling_atomic(nH, T, d; ionization_in_u = false)
chk("adapter, thermal convention = channels total", isapprox(Λt, c.total; rtol = 1e-12))
chk("adapter, χ-in-u convention = total − ciHI + χ_H k2 n_e n_H⁺", isapprox(Λu, c.total - c.H_ci + AC.χH * AC.k2(T) * 0.5 * 0.5; rtol = 1e-12))
chk("ρ form consistent", isapprox(AC.cooling_atomic_rho(AC.μH * nH, T, d), Λu; rtol = 1e-12))
@printf("\n%d passed, %d failed\n", pass[], fail[])
