# Photo-heating audit check (docs/POPIII_MULTIGROUP.md §31). For the star light J = W B_ν(T_c) of a coupled zone:
#  (1) H⁻ photodetachment: the energy rate of the network's rate (4π∫σ J (1 − e^{−hν/kT}) dν per H⁻, the stimulated factor κ_E has) against the
#      absorbed power of G2Opacity's κ_E (H⁻ bound–free alone, T_w = T_c, NLOPACJ): ρ κ_E c E_star, c E_star = 4σ W T_c⁴. Also the network's rate
#      R (ChemNet.photorate) and the mean photon energy ⟨hν⟩ = ∫4πσJ dν / R.
#  (2) Li photoionisation (in the network's photo list, not in κ_E): its power at photo-ionisation equilibrium against that of H⁻.
# Run: G2OPAC=2 julia --project=<env> photoheat_check.jl
using Printf
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "G2Opacity.jl"))
include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl")); include(joinpath(@__DIR__, "ChemRates.jl"))
const CN = ChemNet; const COp = ContinuumOpacity; const G2O = G2Opacity
G2O.ON || error("set G2OPAC=2")
net, photos = ChemRates.modern(li = true, he = true)
pHm = photos[findfirst(p -> p.R === :Hm, photos)]; pLi = photos[findfirst(p -> p.R === :Li, photos)]
σSB = 5.670374419e-5; c = CN.c; kB = CN.kB; h = CN.hP; eV = 1.602176634e-12
B(ν, T) = 2h * ν^3 / c^2 / expm1(h * ν / (kB * T))
"4π ∫ σ(ν) J(ν) s(ν) dν over a log grid from ν0 [erg/s per absorber]; s = 1 or the stimulated factor of gas T."
function eint(σ, J, ν0, νmax, stim; n = 6000)
    νs = [ν0 * (νmax / ν0)^(i / n) for i in 0:n]; s = 0.0
    for i in 1:n
        ν = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]
        s += 4π * σ(ν) * J(ν) * stim(ν) * dν
    end
    s
end
println("H⁻ photodetachment: network rate with the star light J = W B(T_c) against the κ_E exchange (per H⁻; n_H, x_H⁻ cancel)")
for (T, Tc, W, nH) in ((4514.0, 5871.0, 0.345, 3e12), (5205.0, 5916.0, 0.597, 2e12), (4094.0, 5851.0, 0.225, 4e12), (10332.0, 10332.0, 1.0, 1e12), (3000.0, 5000.0, 0.01, 1e12))
    ρ = nH * CN.μH(); nHm = 1e-9 * nH
    J = ν -> W * B(ν, Tc); env = (Trad = Tc, J = J)
    R = CN.photorate(pHm, T, env)
    σ = ν -> COp.σ_Hm_bf(1e4 / (ν / c))
    νmax = 1e6 * c; ν0 = c * 1.0
    E0 = eint(σ, J, ν0, νmax, ν -> 1.0); Es = eint(σ, J, ν0, νmax, ν -> -expm1(-h * ν / (kB * T)))
    κE = G2O.means(ρ, T; n1 = 0.0, n2 = 0.0, nH2 = 0.0, nHe = 0.0, nHp = 0.0, ne = 0.0, nHm = nHm, hbfP = false, Tw = Tc)[1]
    Pκ = ρ * κE * 4σSB * W * Tc^4 / nHm
    @printf("  T %6.0f T_c %6.0f W %.3f: R = %.4e s⁻¹, ⟨hν⟩ = %.3f eV (thermal gain hν − 0.754 = %.3f eV per event) | rate-side %.5e erg/s, κ_E-side %.5e erg/s, ratio %.5f\n",
            T, Tc, W, R, E0 / R / eV, E0 / R / eV - 0.754, Es, Pκ, Pκ / Es)
end
println("Li photoionisation (not in κ_E) at photo-ionisation equilibrium, against the H⁻ power at x_H⁻ = 1e-9 (per H nucleus)")
for (T, Tc, W, nH, xe) in ((4514.0, 5871.0, 0.345, 3e12, 1e-4), (5205.0, 5916.0, 0.597, 2e12, 1e-3), (10332.0, 10332.0, 1.0, 1e12, 1e-1))
    J = ν -> W * B(ν, Tc); env = (Trad = Tc, J = J)
    RLi, ALi = CN.photorates(pLi, T, env)                    # photons/s per Li, recombination coefficient [cm³ s⁻¹]
    f0 = ALi * xe * nH / (RLi + ALi * xe * nH)               # neutral share at photo-ionisation equilibrium
    Pli = CN.yLi * f0 * RLi * 6 * eV                          # erg/s per H nucleus, ⟨hν⟩ taken as 6 eV (threshold 5.39 eV)
    σ = ν -> COp.σ_Hm_bf(1e4 / (ν / c)); Phm = 1e-9 * eint(σ, J, c, 1e6 * c, ν -> 1.0)
    @printf("  T %6.0f T_c %6.0f W %.3f x_e %.0e: R_Li = %.3e s⁻¹, Li⁰ share %.2e, Li power %.3e erg/s/H, H⁻ power %.3e erg/s/H, ratio %.1e\n", T, Tc, W, xe, RLi, f0, Pli, Phm, Pli / Phm)
end
