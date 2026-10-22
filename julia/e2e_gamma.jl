# E2E-A diagnostic: Eddington factor profile Gamma(m) = kappa_R (L_rad + L_conv-free total flux) / (4 pi c G m) of an exported star (e2e_handoff.jl file)
#   TITAN_NOLIB=1 julia --project=env e2e_gamma.jl <file.star>
include(joinpath(@__DIR__, "s4_common.jl"))
using Serialization
h = deserialize(ARGS[1]); n = length(h.rho)
G = Z0.G; c = Z0.c
@printf("M = %.2f Msun, R = %.2f Rsun, L_surf = %.3e Lsun, Gamma_e(0.35) = %.3f\n", h.m[end] / Msun, h.r[end] / Rsun, 4π * h.r[end]^2 * h.Fr[end] / Lsun, 0.35 * 4π * h.r[end]^2 * h.Fr[end] / (4π * c * G * h.m[end]))
@printf("%5s %10s %10s %10s %10s %8s %8s %8s\n", "zone", "q_ext", "r/Rsun", "rho", "T", "kappa", "Gamma", "L/Lsun")
gmax = 0.0
for j in [1:10:n-30; n-29:n]
    kap = StellarOpacity.kappa_R(h.rho[j], h.T[j], 0.75)
    L = 4π * h.r[j+1]^2 * h.Fr[j+1]; Γ = kap * L / (4π * c * G * h.m[j+1])
    @printf("%5d %10.3e %10.3f %10.3e %10.3e %8.3f %8.3f %10.3e\n", j, 1 - h.m[j+1] / h.m[end], h.r[j+1] / Rsun, h.rho[j], h.T[j], kap, Γ, L / Lsun)
end
