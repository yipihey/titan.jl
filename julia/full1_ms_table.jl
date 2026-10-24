# FULL-1: main-sequence table (X_c milestones) and Schaerer 2002 / Marigo 2001 interpolation.  julia --project=env full1_ms_table.jl [csv]
using DelimitedFiles, Printf
f = isempty(ARGS) ? joinpath(@__DIR__, "..", "lib", "popiii", "s5_e2e_abn02z19_fb_ms.csv") : ARGS[1]
A, h = readdlm(f, ','; header = true); h = vec(String.(h)); c(n) = Float64.(A[:, findfirst(==(n), h)])
t = c("t_yr"); Xc = c("Xc")
println("X_c  t[yr]  R  logL  logTeff  Tc  rho_c  Mcc/M")
for x in (0.7, 0.5, 0.3, 0.1, 0.01, 1e-3)
    i = findfirst(<=(x), Xc); i === nothing && (i = length(Xc))
    @printf("%.4g  %.4e  %.2f  %.3f  %.3f  %.4e  %.2f  %.2f\n", Xc[i], t[i], c("R_Rsun")[i], c("logL")[i], c("logTeff")[i], c("Tc")[i], c("rho_c")[i], c("Mcc_over_M")[i])
end
i = argmin(c("R_Rsun")[1:min(end, 500)]); @printf("R minimum (first 500 rows): t=%.4e R=%.3f logL=%.3f\n", t[i], c("R_Rsun")[i], c("logL")[i])
# Schaerer 2002 Table 3/4 (validation.md), log-log interpolation in M between 80 and 120 and (for 200) extrapolation not used
M = 115.7; Ms = (80.0, 120.0); logL = (5.947, 6.243); lT = (4.970, 4.981); tau = (3.012e6, 2.521e6); Q = (5.938e49, 1.069e50)
w = log(M / Ms[1]) / log(Ms[2] / Ms[1])
@printf("Schaerer 2002 at %.1f Msun (interp 80-120 in log M): log L = %.3f, log Teff = %.3f, tau_H = %.3e yr, Q = %.3e\n", M, logL[1] + w * (logL[2] - logL[1]), lT[1] + w * (lT[2] - lT[1]), exp(log(tau[1]) + w * log(tau[2] / tau[1])), exp(log(Q[1]) + w * log(Q[2] / Q[1])))
# Marigo 2001 Table 1: 100 Msun 2.9509e6 (the table ends there); slope from Schaerer tau ~ M^-0.35 (80-120) -> 115.7:
@printf("Marigo 100 Msun 2.951e6 yr; scaled to %.1f Msun with the Schaerer 80-120 slope: %.3e\n", M, 2.9509e6 * (M / 100)^(log(tau[2] / tau[1]) / log(Ms[2] / Ms[1])))
