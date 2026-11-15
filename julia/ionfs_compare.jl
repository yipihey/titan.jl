# IONFS test against the reference run from the same checkpoint: time series and profiles from the RUNREC snapshot streams.
#   julia --project=<env> ionfs_compare.jl <ion run dir> <reference snaps .jls> [out prefix]
# Prints per snapshot: t, dt, L(edge), Ṁ at 1e2/1e3/1e4/1e5 R☉, the ionization front (ray), Q_H at r_b and at 1e3/1e4 R☉; and, at nearest
# times, the x_e and T profiles of both runs.
using Serialization, Printf
include(joinpath(@__DIR__, "runrec.jl"))
const μH = 1.6735e-24 * (1 + 4 / 12)
ion = lineage(ARGS[1]); ref = readstream(ARGS[2])
iy(x, s) = findfirst(==(s), x.spn)
function mdot(x, R)
    i = findfirst(>=(R * Rs), x.r); (i === nothing || i < 2 || i > length(x.ρ)) && return NaN
    -4π * x.r[i]^2 * x.u[i] * 0.5 * (x.ρ[i-1] + x.ρ[i]) * yr / Ms
end
Ledge(x) = 4π * x.r[end-2]^2 * x.fr[end-2] / Ls
"Radius where x_H crosses 0.5 (outermost crossing from inside: the H I front of the gas itself)."
function hfront(x)
    j = iy(x, :H); xh = x.Y[j, :]
    for k in 2:length(xh); (xh[k-1] < 0.5 <= xh[k]) && return 0.5 * (x.r[k] + x.r[k+1]) / Rs; end
    NaN
end
function series(L, lab)
    println("\n== $lab: step t[yr] dt[s] L_edge[L☉] Mdot(1e2,1e3,1e4,1e5 R☉)[M☉/yr] r(x_H=0.5)[R☉] | ray: k0 r_b Q_src kfront r_front Q(1e3) Q(1e4) eps")
    for x in L
        @printf("%d %.6f %.3e %.4e %.3e %.3e %.3e %.3e %.4g", x.step, x.t / yr - 2.2305e7, x.dt, Ledge(x), mdot(x, 1e2), mdot(x, 1e3), mdot(x, 1e4), mdot(x, 1e5), hfront(x))
        if hasproperty(x, :ionfs)
            p = x.ionfs; q(R) = (i = findfirst(>=(R * Rs), x.r); i === nothing ? NaN : p.Q[i])
            @printf(" | %d %.1f %.3e %d %.4g %.3e %.3e %.2e", p.k0 - x.a + 1, p.rb / Rs, p.Qsrc[1], p.kf == 0 ? 0 : p.kf - x.a + 1, p.rf / Rs, q(1e3), q(1e4), p.eps)
        end
        println()
    end
end
series(ion, "IONFS run"); series(ref, "reference")
# profiles at three times of the IONFS run, against the reference snapshot nearest in time
nearest(t) = ref[argmin([abs(y.t - t) for y in ref])]
for x in (isempty(ion) ? Any[] : unique([ion[1], ion[length(ion)÷2+1], ion[end]]))
    y = nearest(x.t)
    println("\n== profiles at step $(x.step) (t $(round(x.t / yr, digits = 4)) yr; reference t $(round(y.t / yr, digits = 4)) yr): zone r[R☉] | x_e T u[km/s] (IONFS) | x_e T u (reference)" * (hasproperty(x, :ionfs) ? " | Q_H in" : ""))
    for k in 60:min(length(x.ρ), length(y.ρ))
        @printf("  %3d %10.4e | %7.4f %9.0f %8.1f | %7.4f %9.0f %8.1f", k, x.r[k] / Rs, x.Y[iy(x, :e), k], x.T[k], x.u[k] / 1e5, y.Y[iy(y, :e), k], y.T[k], y.u[k] / 1e5)
        hasproperty(x, :ionfs) && @printf(" | %.3e", x.ionfs.Q[k])
        println()
    end
end
