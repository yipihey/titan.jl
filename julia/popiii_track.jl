# Checks of a fixed-mesh accretion run from its snapshots (popiii_m2c.jl, SNAPS=1), for docs/POPIII_GOAL.md milestone M1.
#   julia popiii_track.jl <tag> [first step for the budget]
# 1. R(M) against the adiabatic-phase scaling of Stahler et al. 1986, as written in Hosokawa & Omukai 2009 eq. 9,
#    R = 26 Rsun (M/Msun)^0.27 (Mdot/1e-3 Msun/yr)^0.41, averaged over 21 snapshots because the shock moves in whole zones.
# 2. Energy budget of a fixed sphere around the star: change of (gas + kinetic + radiation + gravitational) energy against
#    what crosses the face (gas enthalpy, kinetic and potential energy, advected radiation, radiative flux) and H2 line losses.
#    The whole-domain budget says little, since the CMB field in the far zones dominates it.
#    Inflow is the mass change of the sphere, because the snapshots do not hold the code's own face flux.
using Serialization, Printf
const Rs = 6.957e10; const Ms = 1.989e33; const yr = 3.156e7; const G = 6.674e-8; const Lsun = 3.83e33
tag = ARGS[1]; st0 = parse(Int, get(ARGS, 2, "0"))
ALL = deserialize(joinpath(@__DIR__, "..", "lib", "popiii", "snaps_abn02z19_lines_n200_$tag.jls"))
shock(x) = findfirst(<(-1e6), x.u)

println("M, R(shock), R(eq. 9), ratio, Mdot, L(500 Rsun), GM Mdot/R, ratio, T_max, m(T_max)")
for Mt in (0.7, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 70, 100)
    i = findfirst(x -> (j = shock(x); j !== nothing && x.m[j] / Ms >= Mt), ALL); i === nothing && continue
    r = map(ALL[max(1, i - 10):min(length(ALL), i + 10)]) do x
        j = shock(x); jm = j + 6; n = findfirst(>(500Rs), x.r); Tm, im = findmax(x.T[1:j])
        Md = 4π * x.r[jm]^2 * x.ρ[jm] * abs(x.u[jm])
        [x.m[j] / Ms, x.r[j] / Rs, Md * yr / Ms, 4π * x.r[n]^2 * x.fr[n] / Lsun, G * x.m[j] * Md / x.r[j] / Lsun, Tm, x.m[im] / Ms]
    end
    M, R, Md, Lo, La, Tm, mm = sum(r) / length(r)
    Req = 26 * M^0.27 * (Md / 1e-3)^0.41
    @printf("%7.2f %7.1f %7.1f %6.2f %8.4f %9.2e %9.2e %5.2f %9.2e %6.2f\n", M, R, Req, R / Req, Md, Lo, La, Lo / La, Tm, mm)
end

A = [x for x in ALL if x.step >= st0]
function sphere(x, Rb)
    n = findfirst(>=(Rb), x.r) - 1
    V = 4π / 3 .* (x.r[2:n+1] .^ 3 .- x.r[1:n] .^ 3); dm = x.ρ[1:n] .* V
    mc = (x.m[1:n] .+ x.m[2:n+1]) ./ 2; rc = ((x.r[1:n] .^ 3 .+ x.r[2:n+1] .^ 3) ./ 2) .^ (1 / 3)
    E = sum(dm .* x.e[1:n]) + 0.5sum(dm .* (x.u[1:n] .^ 2 .+ x.u[2:n+1] .^ 2) ./ 2) + sum(V .* x.er[1:n]) - G * sum(mc .* dm ./ rc)
    u = x.u[n+1]; r = x.r[n+1]; up = u < 0 ? n + 1 : n
    (; E, M = sum(dm), q = x.e[up] + x.pg[up] / x.ρ[up] + u^2 / 2 - G * x.m[n+1] / r, Fra = -4π * r^2 * u * (4 / 3) * x.er[up],
       Lr = 4π * r^2 * x.fr[n+1], Ll = sum(V .* x.Λext[1:n]))
end
@printf("\nbudget, steps %d..%d, %.1f yr\nRb[Rsun], dE, predicted, gas in, radiated, lines, (dE - predicted)/radiated\n", A[1].step, A[end].step, (A[end].t - A[1].t) / yr)
for Rb in (250, 300, 500, 1000, 3000)
    C = [sphere(x, Rb * Rs) for x in A]; gas = 0.0; rad = 0.0; line = 0.0; adv = 0.0
    for i in 2:length(A)
        dt = A[i].t - A[i-1].t
        gas += (C[i].M - C[i-1].M) * (C[i].q + C[i-1].q) / 2; adv += (C[i].Fra + C[i-1].Fra) / 2 * dt
        rad += (C[i].Lr + C[i-1].Lr) / 2 * dt; line += (C[i].Ll + C[i-1].Ll) / 2 * dt
    end
    dE = C[end].E - C[1].E; pred = gas + adv - rad - line
    @printf("%6d %+.4e %+.4e %+.3e %+.3e %.2e %+.3f\n", Rb, dE, pred, gas, rad, line, (dE - pred) / rad)
end
