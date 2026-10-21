# S6b analysis of an s6b_run.jl CSV (lib/popiii/s6b_<OUT>.csv): milestones and R*(M), L*(M), T_eff(M), S_ion(M).
#   julia s6c_analyze.jl <OUT> [<OUT of a reference run to compare with>]
using Printf
const Lsun = 3.828e33
function readcsv(tag)
    L = readlines(joinpath(@__DIR__, "..", "lib", "popiii", "s6b_$(tag).csv"))
    hdr = split(L[1], ","); rows = [parse.(Float64, split(l, ",")) for l in L[2:end] if !isempty(strip(l))]
    col(n) = [r[findfirst(==(n), hdr)] for r in rows]
    (; hdr, rows, col)
end
"Interpolate y(x) (x increasing) at xq, linear in ln x."
function at(x, y, xq)
    (xq < x[1] || xq > x[end]) && return NaN
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1); w = log(xq / x[i]) / log(x[i+1] / x[i]); (1 - w) * y[i] + w * y[i+1]
end
function milestones(d)
    M = d.col("M"); R = d.col("R_Rsun"); Ls = d.col("Lstar_Lsun"); Teff = d.col("Teff"); fD = d.col("fD"); Lcno = d.col("Lcno"); Lpp = d.col("Lpp"); Tc = d.col("Tc"); t = d.col("t_yr")
    first_over(v, thr) = (i = findfirst(>(thr), v); i === nothing ? NaN : M[i])
    ip = argmax(R)
    println(@sprintf("  M range %.4g – %.4g Msun, t = %.4g yr", M[1], M[end], t[end] - t[1]))
    println(@sprintf("  luminosity wave: L* > 1e2 / 1e3 / 1e4 Lsun at M = %.2f / %.2f / %.2f", first_over(Ls, 1e2), first_over(Ls, 1e3), first_over(Ls, 1e4)))
    println(@sprintf("  R* peak %.1f Rsun at M = %.2f Msun (L* %.3g Lsun, T_eff,ph %.0f K)", R[ip], M[ip], Ls[ip], Teff[ip]))
    ifD = findfirst(<(0.99), fD); ifD5 = findfirst(<(0.5), fD); ifD1 = findfirst(<(0.1), fD)
    println(@sprintf("  D: f_D < 0.99 / 0.5 / 0.1 at M = %.2f / %.2f / %.2f", ifD === nothing ? NaN : M[ifD], ifD5 === nothing ? NaN : M[ifD5], ifD1 === nothing ? NaN : M[ifD1]))
    Lnuc = Lcno .+ Lpp
    iz = findfirst(i -> Lcno[i] > 0.9 * Ls[i] && Ls[i] > 0, eachindex(M)); ic = findfirst(i -> Lcno[i] > 0.01 * Ls[i] && Ls[i] > 0, eachindex(M))
    println(@sprintf("  CNO ignition (L_CNO > 0.01 L*) at M = %.2f; ZAMS (L_CNO > 0.9 L*) at M = %.2f (R = %.2f Rsun)", ic === nothing ? NaN : M[ic], iz === nothing ? NaN : M[iz], iz === nothing ? NaN : R[iz]))
    println("  R*(M): ", join([@sprintf("%g:%.1f", m, at(M, R, m)) for m in (0.5, 1, 3, 5, 8, 10, 15, 20, 30, 40, 50, 70, 100, 150) if !isnan(at(M, R, m))], "  "))
    println("  L*(M): ", join([@sprintf("%g:%.3g", m, at(M, Ls, m)) for m in (1, 3, 5, 8, 10, 15, 20, 30, 40, 50, 70, 100, 150) if !isnan(at(M, Ls, m))], "  "))
    println("  T_eff,ph(M): ", join([@sprintf("%g:%.0f", m, at(M, Teff, m)) for m in (1, 3, 5, 8, 10, 20, 40, 70, 100, 150) if !isnan(at(M, Teff, m))], "  "))
    if "Sion" in d.hdr
        Si = d.col("Sion"); println("  S_ion(M) [1/s, black body of T_eff,ph from r_ph]: ", join([@sprintf("%g:%.2e", m, at(M, Si, m)) for m in (1, 3, 5, 8, 10, 20, 40, 70, 100, 150) if !isnan(at(M, Si, m))], "  "))
    end
end
"Black-body ionising photon rate [1/s] of a sphere of radius R [cm] and temperature T [K] (photons above 13.6 eV)."
function sion_bb(R, T)
    x0 = 13.6 * 1.602177e-12 / (1.380649e-16 * T)
    I = sum(exp(-n * x0) * (x0^2 / n + 2x0 / n^2 + 2 / n^3) for n in 1:200)
    4π * R^2 * 2π * (1.380649e-16 * T)^3 / (6.62607e-27^3 * 2.99792458e10^2) * I
end
"Table of the protostar at the given masses: t, Mdot, R*, L*, L_acc, T_eff,ph, T_c, f_D, S_ion, L_nuc/L*."
function table(d, Ms = (0.3, 1, 3, 5, 10, 15, 20, 30, 40, 60, 100, 150))
    M = d.col("M"); g(n) = d.col(n)
    haskey_(n) = n in d.hdr
    println(@sprintf("  %7s %9s %10s %9s %10s %10s %7s %9s %7s %10s %8s %9s %10s", "M", "t [yr]", "Mdot", "R* [Rsun]", "L*", "L_acc", "Teff,ph", "T_c", "f_D", "S_ion,ph", "L_nuc/L*", "Teff,*", "S_ion,*"))
    for m in Ms
        isnan(at(M, g("R_Rsun"), m)) && continue
        Md = haskey_("Mdot") ? at(M, g("Mdot"), m) : NaN
        Ln = (at(M, g("Lcno"), m) + at(M, g("Lpp"), m) + at(M, g("LD"), m)) / max(at(M, g("Lstar_Lsun"), m), 1e-30)
        Rs = at(M, g("R_Rsun"), m); Lstar = at(M, g("Lstar_Lsun"), m); Tst = (Lstar * Lsun / (4π * 5.670374e-5 * (Rs * 6.957e10)^2))^0.25
        println(@sprintf("  %7.3g %9.4g %10.3e %9.1f %10.3e %10.3e %7.0f %9.3e %7.3f %10.2e %8.3f %9.0f %10.2e", m, at(M, g("t_yr"), m), Md, Rs, Lstar, at(M, g("Lacc_Lsun"), m), at(M, g("Teff"), m), at(M, g("Tc"), m), at(M, g("fD"), m), at(M, g("Sion"), m), Ln, Tst, sion_bb(Rs * 6.957e10, Tst)))
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    d = readcsv(ARGS[1]); println(ARGS[1], ":"); milestones(d); table(d)
    if length(ARGS) > 1
        d2 = readcsv(ARGS[2]); println(ARGS[2], " (reference):"); milestones(d2)
        M = d.col("M"); M2 = d2.col("M")
        println("  ratio ", ARGS[1], " / ", ARGS[2], " at equal M for R*, L*, T_c:")
        for m in (0.5, 1, 3, 5, 8, 10, 15, 20, 30, 40, 50, 60)
            a = (at(M, d.col("R_Rsun"), m), at(M, d.col("Lstar_Lsun"), m), at(M, d.col("Tc"), m)); b = (at(M2, d2.col("R_Rsun"), m), at(M2, d2.col("Lstar_Lsun"), m), at(M2, d2.col("Tc"), m))
            any(isnan, a) || any(isnan, b) || println(@sprintf("   M = %5.1f: R %.1f / %.1f = %.3f   L* %.3g / %.3g = %.3f   T_c %.3g / %.3g = %.3f", m, a[1], b[1], a[1] / b[1], a[2], b[2], a[2] / b[2], a[3], b[3], a[3] / b[3]))
        end
    end
end
