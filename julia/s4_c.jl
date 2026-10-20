# S4 test c1: mixing operator on a synthetic Lagrangian grid (no TITAN state): conservation, uniformity, untouched radiative zones.
#   (c2, mixing of the real convective core of the relaxed star, is part of s4_d.jl, MIXTEST=1.)
using Printf
include(joinpath(@__DIR__, "s4_common.jl"))
n = 60; NS = NN.NS
r = cumsum(fill(1e10, n + 1)); ρ = fill(30.0, n)
dm = [4π * (r[k+1]^3 - r[k]^3) / 3 * ρ[k] for k in 1:n]       # g, consistent with ρ = 30 and the radii
area = [(k == 1 || k == n + 1) ? 0.0 : 4π * r[k]^2 * ρ[k] for k in 1:n+1]
conv = 15:35                                                  # zones; interfaces k in 16:35 convective
D = zeros(n + 1); for k in 16:35; D[k] = 3e16; end             # cm²/s: D = vℓ/3 with ℓ = 1e11 cm, v = 9e5 cm/s (core of a 100 Msun star)
Y = zeros(NS, n); Y[1, :] .= 0.7; Y[4, :] .= 0.075
Y[1, 15:25] .= 0.5; Y[1, 26:35] .= 0.7; Y[4, 15:25] .= 0.1; Y[4, 26:35] .= 0.075       # a step inside the convective region
Y0 = copy(Y)
tot0 = [sum(Y[i, :] .* dm) for i in 1:NS]
# turnover time ℓ/v: with D = vℓ/3 and ℓ = 1e11 cm: v = 3D/ℓ
ℓ = 1e11; v = 3 * 3e16 / ℓ; tt = ℓ / v
# diffusion time across the convective region of mass M_c: t_D ~ (Δm)²/(area² D)
@printf("turnover time l/v = %.3e s (region width %.1f l)\n", tt, (r[36] - r[15]) / ℓ)
res = Dict()
function evolve(T, nsub)
    Ym = copy(Y0); for _ in 1:nsub; SM.mix_species!(Ym, D, dm, area, T / nsub); end; Ym
end
for (lab, T) in (("1 t_turn", tt), ("3 t_turn", 3tt), ("10 t_turn", 10tt), ("30 t_turn", 30tt), ("100 t_turn", 100tt))
    for nsub in (1, 200)
        Ym = evolve(T, nsub)
        tot1 = [sum(Ym[i, :] .* dm) for i in 1:NS]
        cons = maximum(abs(tot1[i] - tot0[i]) / tot0[i] for i in 1:NS if tot0[i] > 0)
        xc = Ym[1, 15:35]; spread = (maximum(xc) - minimum(xc)) / abs(Y0[1, 15] - Y0[1, 30])
        out = maximum(abs.(Ym[:, [1:14; 36:n]] .- Y0[:, [1:14; 36:n]]))
        res[(lab, nsub)] = (cons, spread, out)
        @printf("  T = %-10s %3d step(s): conservation %.2e; spread of X_H in the convective region / initial step = %.3e; max change in radiative zones = %.1e\n", lab, nsub, cons, spread, out)
    end
end
ok = all(x -> x[1] < 1e-14 && x[3] == 0.0, values(res)) && res[("30 t_turn", 200)][2] < 1e-3 && res[("10 t_turn", 200)][2] < 0.05 &&
     res[("3 t_turn", 200)][2] < res[("1 t_turn", 200)][2]
Ym = evolve(1e4tt, 1)
xm0 = sum(Y0[1, 15:35] .* dm[15:35]) / sum(dm[15:35])
@printf("fully mixed X_H = %.10f vs mass-weighted initial mean over the region %.10f\n", Ym[1, 20], xm0)
println(ok ? "PASS" : "FAIL")
