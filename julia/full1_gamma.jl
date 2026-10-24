# FULL-1: Eddington factor Gamma = (L* + L_acc) / L_Edd,es (kappa_es = 0.34) along the protostar track.  julia --project=env full1_gamma.jl [star csv]
using DelimitedFiles, Printf
f = isempty(ARGS) ? joinpath(@__DIR__, "..", "lib", "popiii", "s6b_e2e_abn02z19_fb_star.csv") : ARGS[1]
A, h = readdlm(f, ','; header = true); h = vec(String.(h)); c(n) = Float64.(A[:, findfirst(==(n), h)])
M = c("M"); LE = 4π * 2.99792458e10 * 6.674e-8 * 1.98847e33 .* M ./ 0.34 / 3.828e33
G = (c("Lstar_Lsun") .+ c("Lacc_Lsun")) ./ LE; G0 = c("Lstar_Lsun") ./ LE
for m in (10, 30, 60, 80, 100, 105, 110, 113, 115)
    i = argmin(abs.(M .- m)); abs(M[i] - m) < 0.02m && @printf("M=%.1f  Gamma(L*+Lacc) = %.3f  Gamma(L*) = %.3f  R = %.1f\n", M[i], G[i], G0[i], c("R_Rsun")[i])
end
i = argmax(G); @printf("max Gamma(L*+Lacc) = %.3f at M = %.1f, R = %.1f; max Gamma(L*) = %.3f; final Gamma = %.3f / %.3f\n", G[i], M[i], c("R_Rsun")[i], maximum(G0), G[end], G0[end])
