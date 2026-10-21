# S6a: validation table and curve features from lib/popiii/s6b_<tag>.csv (columns of s6b_run.jl).
#   julia s6_analyze.jl TAG [Mdot_Msun_per_yr]
using Printf, DelimitedFiles
tag = ARGS[1]; Mdot = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.1e-3
path = joinpath(@__DIR__, "..", "lib", "popiii", "s6b_$(tag).csv")
D, hdr = readdlm(path, ','; header = true); hdr = vec(hdr)
col(n) = D[:, findfirst(==(n), hdr)]
t, M, R, Ls, La, Teff, Tc, ρc, Tmax, mTmax, LD, Lpp, Lcno, fD, tKH, tacc, Sion, rph = (col(n) for n in ("t_yr", "M", "R_Rsun", "Lstar_Lsun", "Lacc_Lsun", "Teff", "Tc", "rho_c", "Tmax", "mTmax", "LD", "Lpp", "Lcno", "fD", "tKH_yr", "tacc_yr", "Sion", "rph_Rsun"))
Lnuc = LD .+ Lpp .+ Lcno
at(m, v) = (i = findfirst(>=(m), M); i === nothing ? NaN : v[i])
@printf("%s: %d rows, M = %.3f … %.3f Msun, t = %.4g yr, Mdot = %.3e Msun/yr\n", tag, length(M), M[1], M[end], t[end], Mdot)
@printf("%7s %9s %10s %10s %9s %9s %10s %9s %8s %8s %8s %9s\n", "M", "R*[Rsun]", "L*[Lsun]", "Lacc", "Teff_ph", "Tc[K]", "rho_c", "Tmax[K]", "m(Tmax)", "L_D/L*", "f_D", "tKH/tacc")
for m in (0.3, 0.5, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60)
    M[end] >= m || break
    i = findfirst(>=(m), M)
    @printf("%7.2f %9.2f %10.3e %10.3e %9.0f %9.3e %10.3e %9.3e %8.3f %8.2e %8.4f %9.2f\n", M[i], R[i], Ls[i], La[i], Teff[i], Tc[i], ρc[i], Tmax[i], mTmax[i], LD[i] / max(Ls[i], 1e-30), fD[i], tKH[i] / tacc[i])
end
# phase I: R* ∝ M^0.27 Mdot^0.41 (SPS86 eq. 21): local slope of ln R vs ln M between 0.5 and 3 Msun
sel = findall(m -> 0.5 <= m <= 3.0, M)
if length(sel) > 5
    X = log.(M[sel]); Y = log.(R[sel]); A = [ones(length(X)) X]; c = A \ Y
    @printf("phase I: d ln R*/d ln M (0.5–3 Msun) = %.3f (SPS86/OP03: 0.27)\n", c[2])
end
# R peak and its mass; luminosity wave: largest jump of ln L* between neighbouring masses; D-burning episodes: ε_D L*/M* criterion ≈ L_D/L* > 0.1 (L*/M* average)
ip = argmax(R); @printf("R* maximum %.1f Rsun at M = %.2f Msun\n", R[ip], M[ip])
for Lth in (30.0, 100.0, 1e3, 1e4)      # luminosity wave: mass at which L* (smoothed over 9 rows) first exceeds a threshold, after 1 Msun
    sm = [sum(Ls[max(i - 4, 1):min(i + 4, end)]) / length(max(i - 4, 1):min(i + 4, length(Ls))) for i in eachindex(Ls)]
    i = findfirst(j -> M[j] > 1.0 && sm[j] > Lth, eachindex(Ls)); i === nothing || @printf("L* > %.0f Lsun first at M = %.2f Msun (R* = %.1f Rsun)\n", Lth, M[i], R[i])
end
iD = findall(x -> x > 0.1, LD ./ max.(Ls, 1e-30)); isempty(iD) || @printf("L_D > 0.1 L*: M = %.2f … %.2f Msun; max L_D/L* = %.2f at M = %.2f; f_D(end) = %.4f\n", M[iD[1]], M[iD[end]], maximum(LD ./ Ls), M[argmax(LD ./ Ls)], fD[end])
# phase II: L* ∝ M^3 and t_KH ≈ t_acc
sel = findall(m -> 15 <= m <= 40, M)
if length(sel) > 5
    X = log.(M[sel]); c = [ones(length(X)) X] \ log.(Ls[sel]); c2 = [ones(length(X)) X] \ log.(R[sel])
    @printf("phase II (15–40 Msun): d ln L*/d ln M = %.2f (L* ∝ M³), d ln R*/d ln M = %.2f (R* ∝ Mdot/M²: −2); t_KH/t_acc = %.2f … %.2f\n", c[2], c2[2], tKH[sel[1]] / tacc[sel[1]], tKH[sel[end]] / tacc[sel[end]])
end
# ZAMS: Lnuc ≈ L* with H burning: first row with L_pp + L_CNO > 0.5 L_total
iz = findfirst(i -> (Lpp[i] + Lcno[i]) > 0.5 * (Ls[i] + 1e-30) && Lcno[i] > 0.1 * Lpp[i], eachindex(M)); iz === nothing || @printf("H burning (pp+CNO > L*/2, CNO > 0.1 pp): M = %.1f Msun, R* = %.2f Rsun, L* = %.3g Lsun, Tc = %.3e K\n", M[iz], R[iz], Ls[iz], Tc[iz])
@printf("S_ion(end) = %.3e s^-1 (T_eff,ph = %.0f K, r_ph = %.1f Rsun)\n", Sion[end], Teff[end], rph[end])
