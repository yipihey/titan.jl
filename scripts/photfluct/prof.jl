# profile around shock/photosphere for given steps:  julia prof.jl <dir> step1,step2,...
using Serialization, Printf
isdefined(Main, :NucNet) || (@eval Main module NucNet end)
include(joinpath(@__DIR__, "..", "..", "julia", "runrec.jl"))
const arad_ = 7.5657e-15
dir = ARGS[1]; steps = parse.(Int, split(ARGS[2], ","))
L = lineage(dir)
for x in L
    x.step in steps || continue
    n = length(x.ρ); r = x.r
    τ = zeros(n + 1); for k in n:-1:1; τ[k] = τ[k+1] + x.chi[k] * x.ρ[k] * (r[k+1] - r[k]); end
    Lf = 4π .* r .^ 2 .* x.fr ./ Ls
    kph = findlast(>(2 / 3), τ); js = findfirst(<(-2e5), x.u)
    @printf("== step %d t %.6e dt %.3e kph %d js %d\n", x.step, x.t, x.dt, kph, js)
    println("  k  r[Rsun](iface)  u[km/s]   L(iface)    tau(iface) | zone: dr[Rs] rho  T  Trad  tau_z  F/cE(iface)  kP*rho*c [1/s]")
    for k in max(js - 4, 1):min(kph + 6, n)
        Tr = (max(x.er[k], 0) / arad_)^0.25
        fce = x.fr[k] / (cl * 0.5 * (x.er[max(k - 1, 1)] + x.er[k]))
        @printf("%3d %10.4f %9.3f %11.4g %11.4g | %8.4f %10.3e %8.1f %8.1f %10.3e %10.3e %10.3e\n", k, r[k] / Rs, x.u[k] / 1e5, Lf[k], τ[k],
                (r[k+1] - r[k]) / Rs, x.ρ[k], x.T[k], Tr, x.chi[k] * x.ρ[k] * (r[k+1] - r[k]), fce, x.kp[k] * x.ρ[k] * cl)
    end
end
