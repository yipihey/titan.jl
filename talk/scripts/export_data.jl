# Exports the Pop III numbers the talk plots to talk/data/*.csv, from the run snapshots in lib/popiii (not in git).
# From talk/: julia --project=<env with the repository's packages> scripts/export_data.jl
using Serialization, Printf
const L = joinpath(@__DIR__, "..", "..", "lib", "popiii"); const D = joinpath(@__DIR__, "..", "data")
const Rs = 6.957e10; const Ms = 1.989e33; const Ls = 3.83e33; const yr = 3.156e7; const G = 6.674e-8; const mH = 1.6726e-24
snaps(t) = deserialize(joinpath(L, "snaps_abn02z19_lines_n200_$t.jls"))
# the hot end of the central track
open(joinpath(D, "collapse_hot.csv"), "w") do io
    println(io, "run,log_n,Tc")
    for t in ("rad", "n3", "c2"), x in snaps(t); @printf(io, "%s,%.4f,%.1f\n", t, log10(x.ρ[1] / mH), x.T[1]); end
end
# accretion history on the fixed mesh; mass and radius at the point where the infall has slowed to 10 km/s
open(joinpath(D, "accretion.csv"), "w") do io
    println(io, "run,step,t_yr,M,R,Mdot,L500,Tc,umin")
    for t in ("c1", "c2", "c3")
        A = [x for x in snaps(t) if x.step >= 23220]
        for x in A
            j = findfirst(<(-1e6), x.u); jm = min(j + 6, length(x.u)); n = findfirst(>(500Rs), x.r)
            @printf(io, "%s,%d,%.5f,%.5f,%.3f,%.4f,%.4e,%.1f,%.2f\n", t, x.step, (x.t - A[1].t) / yr, x.m[j] / Ms, x.r[j] / Rs,
                    4π * x.r[jm]^2 * x.ρ[jm] * abs(x.u[jm]) * yr / Ms, 4π * x.r[n]^2 * x.fr[n] / Ls, x.T[1], minimum(x.u) / 1e5)
        end
    end
end
# profiles: the same mass in the runs without and with convection, and the mesh-independent last state of each
function prof(name, x)
    open(joinpath(D, "profile_$name.csv"), "w") do io
        println(io, "r,m,rho,T,u,s,L")
        for j in 1:length(x.ρ)
            @printf(io, "%.5e,%.6f,%.5e,%.1f,%.4f,%.5e,%.4e\n", x.r[j] / Rs, x.m[j] / Ms, x.ρ[j], x.T[j], x.u[j] / 1e5, x.pg[j] / x.ρ[j]^(5 / 3), 4π * x.r[j]^2 * x.fr[j] / Ls)
        end
    end
    j = findfirst(<(-1e6), x.u); @printf("%s: step %d, M %.3f, R %.1f\n", name, x.step, x.m[j] / Ms, x.r[j] / Rs)
end
c1 = snaps("c1"); c2 = snaps("c2")
mshock(x) = x.m[findfirst(<(-1e6), x.u)] / Ms
target = 1.0
for (n, A) in (("noconv", c1), ("conv", c2))
    i = findfirst(x -> x.step > 26000 && mshock(x) >= target, A); prof(n, A[i])
end
# opacity against temperature
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "julia"))
include(joinpath(@__DIR__, "..", "..", "julia", "StellarOpacity.jl"))
open(joinpath(D, "opacity.csv"), "w") do io
    println(io, "rho,logT,kR")
    for ρ in (1e-12, 1e-8, 1e-4, 1.0), lT in 3.0:0.02:8.0
        @printf(io, "%.0e,%.2f,%.5e\n", ρ, lT, StellarOpacity.kappa_R(ρ, 10.0^lT, 0.75))
    end
end
