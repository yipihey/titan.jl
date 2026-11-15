# Zone table of a RUNREC snapshot for the ionizing-precursor analysis (post-processing only).
#   julia --project=<env> ionprec_zones.jl <snaps_*.jls> [index|last]
using Serialization, Printf
include(joinpath(@__DIR__, "runrec.jl"))
v = readstream(ARGS[1]); x = length(ARGS) >= 2 && ARGS[2] != "last" ? v[parse(Int, ARGS[2])] : v[end]
println(length(v), " snaps; step ", x.step, " t ", x.t / yr, " yr; species ", x.spn)
n = length(x.ρ); τ = zeros(n + 1)
for k in n:-1:1; τ[k] = τ[k+1] + x.chi[k] * x.ρ[k] * (x.r[k+1] - x.r[k]); end
@printf("nz %d kph %d Teff %.0f M %.2f Msun\n", n, x.kph, x.Teff, x.m[end] / Ms)
iy(s) = findfirst(==(s), x.spn)
μH = 1.6735e-24 * (1 + 4 / 12)
for k in 1:n
    nH = x.ρ[k] / μH
    @printf("%3d r %10.4e Rs ρ %9.2e nH %9.2e T %9.0f τ %9.3e eq %d xe %9.2e xH %9.2e xHe+ %9.2e xHe++ %9.2e L %9.2e Tr %8.0f u %8.1f\n",
            k, x.r[k] / Rs, x.ρ[k], nH, x.T[k], τ[k], x.eq[k], x.Y[iy(:e), k], x.Y[iy(:H), k], x.Y[iy(:Hep), k], x.Y[iy(:Hepp), k],
            4π * x.r[k]^2 * x.fr[k] / Ls, (max(x.er[k], 0) / arad)^0.25, x.u[k] / 1e5)
end
