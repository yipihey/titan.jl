# SUBZONE unit test (TitanResidFast.sz_tau / sz_avchi!): derivatives of the sub-zone interface optical depth carried by
# ForwardDiff against central finite differences, limits away from a front, and the front values of zdiagB's zones 63/64.
#   TITAN_NOLIB=1 julia --project=<env> subzone_test.jl
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl"))
const TX = TitanResidFast
nbad = Ref(0)
check(ok, msg) = (ok || (nbad[] += 1); println(ok ? "  ok   " : "  FAIL ", msg))

# 1. sz_avchi! on a 4-zone strip (interfaces 2..4): d ln avchi_k / d ln (χ, ρ) of all zones by AD vs central FD
TX.SUBZONE[] = true
function strip(p)        # p = (χ1..χ4, ρ1..ρ4): avchi at the three interior interfaces of zones 1..4 (r = 1..5 R_sun)
    T = eltype(p); χ = p[1:4]; ρ = p[5:8]
    r = [1.0, 1.08, 1.17, 1.26, 1.36] .* 9e11; dvol = [(r[k+1]^3 - r[k]^3) / 3 for k in 1:4]; rmu = r .^ 2
    av = T[2 / (1 / χ[k-1] + 1 / χ[k]) for k in 2:4]; avf = vcat(zero(T), av, zero(T))
    TX.sz_avchi!(avf, χ, ρ, dvol, rmu, 1, 4)
    avf[2:4]
end
cases = (("front (zdiagB 63|64: hot thick, cold thin)", [543.0, 21.7, 3.1e-5, 2.2e-5, 6.0e-7, 2.9e-7, 1.73e-7, 1.2e-7]),
         ("handover (cold zone at tau ~ 1)",            [543.0, 21.7, 3.0e-4, 2.2e-5, 6.0e-7, 2.9e-7, 1.73e-7, 1.2e-7]),
         ("thick interior",                             [2e4, 1.4e4, 1.0e4, 7.0e3, 1.8e-5, 1.2e-5, 8.3e-6, 5.6e-6]),
         ("thin envelope",                              [2.6e-5, 2.1e-5, 1.9e-5, 1.8e-5, 1.2e-7, 9.4e-8, 8.1e-8, 7.4e-8]))
for (name, p) in cases
    # elasticities ∂ln avchi_i / ∂ln p_j (dimensionless, O(1)): AD vs central FD in ln p (step 1e-6)
    f(x) = log.(strip(exp.(x))); x0 = log.(p)
    J = ForwardDiff.jacobian(f, x0); worst = 0.0; h = 1e-6
    for j in eachindex(x0)
        xp = copy(x0); xm = copy(x0); xp[j] += h; xm[j] -= h
        fd = (f(xp) .- f(xm)) ./ 2h
        worst = max(worst, maximum(abs.(fd .- J[:, j])))
    end
    check(worst < 1e-6, @sprintf("%-46s AD vs FD elasticities of avchi: worst absolute difference %.2e (max |elasticity| %.3f)", name, worst, maximum(abs, J)))
end

# 2. limits: away from a front the interface depth is TITAN's (ratio Δτ_SZ/Δτ_T), at the front τ_f + b-like and less sensitive to b
τf = TX.SZTAUF[]
for (a, b, tol) in ((1e7, 2e6, 1e-12), (1e3, 1.2e3, 1e-10), (0.05, 0.06, 1e-5), (0.1, 0.12, 2e-4))
    T = 4a * b / (a + b); s, w = TX.sz_tau(a, b, T, τf)
    check(abs(s / T - 1) < tol, @sprintf("no front: a %.3g b %.3g  dtau_SZ/dtau_T - 1 = %.2e (w %.2e)", a, b, s / T - 1, w))
end
Mm, Mp = 1.08e4, 6.98e3; a = 21.7Mm
dT(χ) = (Mm + Mp) * 2 / (1 / 21.7 + 1 / χ)
for χ in (1e-5, 3.1e-5, 6.2e-5, 1e-4, 3e-4, 1e-3)
    s, w = TX.sz_tau(a, χ * Mp, dT(χ), τf)
    g = ForwardDiff.derivative(x -> log(TX.sz_tau(a, exp(x) * Mp, dT(exp(x)), τf)[1]), log(χ))
    @printf("  front: chi+ %.1e  b %.3f  dtau_T %.3f  dtau_SZ %.3f  w %.3f  dln(dtau_SZ)/dln(chi+) %.3f (TITAN 1.0)\n", χ, χ * Mp, dT(χ), s, w, g)
end
println(nbad[] == 0 ? "SUBZONE unit test: all passed" : "SUBZONE unit test: $(nbad[]) FAILED")
