# Tests of Spectrum.jl. (1) Isothermal, purely absorbing, optically thick sphere: L_ν = 4π²R²B_ν.
# (2) Isothermal sphere with scattering (ρ = 1e-9): L_ν ≤ 4π²R²B_ν (scattering lowers the emergent flux of an
#     isothermal medium; at 3e15 Hz below the Lyman edge ε is smallest).
#   julia -t 2 --project=env spectrum_test.jl
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl"))
using Printf
const SP = Spectrum
R = 1e14; n = 120
r = R .* (1 .- 10 .^ range(0, -6, length = n)) .+ 1e-6R          # zone centres crowding towards R
r = sort(unique(vcat(range(0.5R, 0.999R, length = 60), r[r .> 0.999R])))
T = 1e4
for ρ in (1e-6, 1e-9)
    atm = (r = r, ρ = fill(ρ, length(r)), T = fill(T, length(r)))
    ν = [3e14, 1e15, 3e15]
    L, its = SP.spectrum(atm, ν)
    for i in eachindex(ν)
        @printf("ρ=%.0e ν=%.1e: L/(4π²R²B) = %.4f (ALI its %d)\n", ρ, ν[i], L[i] / (4π^2 * r[end]^2 * SP.Bν(ν[i], T)), its[i])
    end
end
# (3) The direct Λ-matrix solution (used where the diagonal ALI converges too slowly) against ALI iterated to
#     convergence, on the scattering-dominated case at the three frequencies; and the Verner cross sections at threshold.
atm = (r = r, ρ = fill(1e-9, length(r)), T = fill(T, length(r)))
op = SP.Opac(atm); _, Rr = SP.rays(atm.r; ncore = 12); Rr = SP.weights!(atm, Rr)
for ν in (3e14, 1e15, 3e15)
    La, _, ia = SP.monochromatic(atm, op, Rr, ν; maxit = 5000, maxit_direct = false)
    Ld, _, _ = SP.monochromatic(atm, op, Rr, ν; maxit = 1)
    @printf("ν=%.1e: direct / converged ALI (%d its) = %.6f\n", ν, ia, Ld / La)
end
@printf("Verner σ at threshold: H I %.3e (6.30e-18), He I %.3e (7.4e-18) cm²\n", SP.σH1(13.6), SP.σHe1(24.59))
