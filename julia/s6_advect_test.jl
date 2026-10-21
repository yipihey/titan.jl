# S6a: tests of the composition advection through a moving mass grid (StarMaterial.advect_composition!).
#   TITAN_NOLIB=1 julia --project=env s6_advect_test.jl
# (1) random grid motion with inflow: Σ M Y_i conserved (to round-off, apart from the inflow), Σ A_i Y_i = 1, 0 ≤ Y; (2) a step profile of the matter
# (X_H: 0.75 → 0 at mass coordinate m_s) seen by a grid whose interfaces drift upward in mass by a fraction c of a zone per step: the step must stay at
# m_s; first order vs second order: width of the transition and L1 error after 200 steps; (3) a stationary grid is the identity.
include(joinpath(@__DIR__, "s4_common.jl"))
using Random, Printf
ns = SM.NS; AA = NN.AA
function randY(n, rng)
    Y = zeros(ns, n)
    for j in 1:n
        x = rand(rng) * 0.75; d = rand(rng) * 3e-5 * x; he = (1 - x - d) / 4
        Y[SM.IH, j] = x; Y[SM.ID, j] = d / 2; Y[NN.IDX[:he4], j] = he
        Y[:, j] ./= sum(AA[i] * Y[i, j] for i in 1:ns)
    end
    Y
end
rng = MersenneTwister(1)
pass = true
# ---- (1) conservation under random grid motion, with inflow
n = 120; Y = randY(n, rng); Yin = zeros(ns); Yin[SM.IH] = 0.75; Yin[SM.ID] = 1.25e-5; Yin[NN.IDX[:he4]] = 0.25 / 4
Yin ./= sum(AA[i] * Yin[i] for i in 1:ns)
M = 1 .+ rand(rng, n)
cmax = 0.0; ncons = 0.0; nsum = 0.0; ymin = 0.0
sum0 = [sum(AA[i] * Y[i, j] for i in 1:ns) for j in 1:n]; sum0[:] .= sum0       # (the random zones are normalised to Σ A Y = 1 only up to the isotopic masses: compare with the initial sums)
for step in 1:300
    δ = 0.3 .* (rand(rng, n + 1) .- 0.5) .* [min(M[min(j, n)], M[max(j - 1, 1)]) for j in 1:n+1]
    δ[1] = 0.0; δ[n+1] = abs(δ[n+1])                    # regular centre; inflow at the outer boundary
    global M, cmax, ncons, nsum, ymin
    tot0 = [sum(M[j] * Y[i, j] for j in 1:n) for i in 1:ns]
    infl = δ[n+1] .* Yin
    r = SM.advect_composition!(Y, copy(M), δ, Yin; order = 2)
    M = r.M
    tot1 = [sum(M[j] * Y[i, j] for j in 1:n) for i in 1:ns]
    ncons = max(ncons, maximum(abs.(tot1 .- tot0 .- infl) ./ max.(tot0, 1e-300) .* (tot0 .> 0)))
    nsum = max(nsum, maximum(abs(sum(AA[i] * Y[i, j] for i in 1:ns) - sum0[1]) / sum0[1] for j in 1:n))
    ymin = min(ymin, minimum(Y)); cmax = max(cmax, r.cfl)
end
@printf("(1) random motion, inflow, 300 steps: max CFL %.2f, max |Δ(Σ M Y_i) − inflow|/Σ M Y_i = %.2e, max |Σ A Y − 1| = %.2e, min Y = %.2e\n", cmax, ncons, nsum, ymin)
ok1 = ncons < 1e-13 && nsum < 1e-13 && ymin >= 0
global pass &= ok1
# ---- (2) step profile, grid drifting upward in mass
function stepprofile(order; n = 200, c = 0.15, nsteps = 200, ms = 120.0)
    Mz = ones(n); m = [0.0; cumsum(Mz)]                  # interface masses
    Y = zeros(ns, n)
    # cell-averaged step in the matter: X_H = 0.75 for m < ms (below), 0 above (in the mass coordinate of the matter)
    for j in 1:n
        lo, hi = m[j], m[j+1]; f = clamp((ms - lo) / (hi - lo), 0.0, 1.0)
        Y[SM.IH, j] = 0.75 * f; Y[NN.IDX[:he4], j] = 0.25 / 4 * f + 0.25 * (1 - f)       # (pure He above)
    end
    Mc = copy(Mz); cflmax = 0.0
    for s in 1:nsteps
        δ = fill(c, n + 1); δ[1] = 0.0                    # interior mass of every interface (but the centre) grows by c per step
        r = SM.advect_composition!(Y, Mc, δ, Y[:, n]; order); Mc = r.M; cflmax = max(cflmax, r.cfl)
        # boundary inflow carries the abundances of the top zone: irrelevant for the step far below
    end
    # exact: grid interface j (j≥2) now at mass-coordinate-in-matter m_j' = m_j + c·nsteps (the matter below ms stays X_H = 0.75).
    # The matter between interface j and j+1 (new zone) was "m_j' .. m_j'+M_j"; here positions are measured in the matter's original mass coordinate
    mnew = [0.0; cumsum(Mc)]
    xs = Y[SM.IH, :] ./ 0.75
    ex = [clamp((ms - mnew[j]) / (mnew[j+1] - mnew[j]), 0.0, 1.0) for j in 1:n]
    l1 = sum(abs.(xs .- ex) .* Mc) / 1.0                  # L1 error in units of mass of the step (mass-weighted)
    width = sum(x * (1 - x) for x in xs)                   # number of zones between 10 % and 90 %, loosely: Σ x(1−x)
    tot = sum(Mc .* Y[SM.IH, :]); tot0 = 0.75 * ms
    (; l1, width, cfl = cflmax, tot, tot0)
end
for order in (1, 2)
    r = stepprofile(order)
    @printf("(2) step profile, order %d: L1 error %.3f (units of zone mass), Σ x(1−x) = %.2f zones, CFL %.2f, ΣM X_H = %.6f (inflow-corrected exact %.6f)\n", order, r.l1, r.width, r.cfl, r.tot, r.tot0)
end
r1 = stepprofile(1); r2 = stepprofile(2)
ok2 = r2.l1 < 0.5 * r1.l1 && r2.width < 3.0
global pass &= ok2
# ---- (3) identity
Y = randY(30, rng); Y0 = copy(Y); M0 = 1 .+ rand(rng, 30)
r = SM.advect_composition!(Y, copy(M0), zeros(31), Yin; order = 2)
ok3 = maximum(abs.(Y .- Y0)) < 1e-15 && r.M == M0
global pass &= ok3
@printf("(3) stationary grid: identity %s\n", ok3)
# ---- (4) large grid motion (several zone masses per step): substepping keeps conservation and positivity
n = 150; Y = randY(n, rng); M = 1 .+ rand(rng, n); tot0 = [sum(M[j] * Y[i, j] for j in 1:n) for i in 1:ns]
δ = [3.0 * sin(2π * j / 200 + 0.3) for j in 1:n+1]; δ[1] = 0.0; δ[n+1] = 0.0
r = SM.advect_composition!(Y, copy(M), δ, Yin; order = 2)
tot1 = [sum(r.M[j] * Y[i, j] for j in 1:n) for i in 1:ns]
c4 = maximum(abs.(tot1 .- tot0) ./ max.(tot0, 1e-300) .* (tot0 .> 0)); ymin4 = minimum(Y); s4 = maximum(abs(sum(AA[i] * Y[i, j] for i in 1:ns) - 1) for j in 1:n)
ok4 = c4 < 1e-13 && ymin4 >= 0 && s4 < 1e-13 && r.nsub > 1
@printf("(4) large motion: CFL (last substep) %.2f, %d substeps, conservation %.2e, min Y %.1e, |ΣAY−1| %.1e\n", r.cfl, r.nsub, c4, ymin4, s4)
global pass &= ok4
# ---- (5) self-similar growth (matter flows down through every interface), huge CFL in the thin surface zones; step in the matter
let n = 120
    qs = vcat(range(0.0, 0.9; length = 61), 1 .- 10 .^ range(log10(0.1), -7; length = 61)[2:end])       # interior fractions: linear to 0.9, then log-spaced exterior mass fraction down to 1e-7
    mI = qs .* 1.0; mI[end] = 1.0; Mz = diff(mI)
    ms = 0.5                                                            # step in the matter at mass coordinate 0.5 of the initial mass
    Y = zeros(ns, n)
    for j in 1:n
        f = clamp((ms - mI[j]) / (mI[j+1] - mI[j]), 0.0, 1.0)         # fraction of zone j below ms
        Y[SM.IH, j] = 0.75 * f; Y[NN.IDX[:he4], j] = (0.25 / 4) * f + 0.25 * (1 - f) * 1.0
    end
    Y .*= 1; for j in 1:n; Y[:, j] ./= sum(AA[i] * Y[i, j] for i in 1:ns); end
    Yacc = zeros(ns); Yacc[SM.IH] = 0.75; Yacc[NN.IDX[:he4]] = 0.25 / 4; Yacc ./= sum(AA[i] * Yacc[i] for i in 1:ns)
    Mc = copy(Mz); Mtot0 = sum(Mz); tot0 = [sum(Mc[j] * Y[i, j] for j in 1:n) for i in 1:ns]; infl = zeros(ns); cflmax = 0.0
    nst = 100; f = 0.8 / nst                                           # M grows by a factor 1.8 in 100 steps
    for s in 1:nst
        Mt = sum(Mc); m = [0.0; cumsum(Mc)]
        δ = f .* m; δ[1] = 0.0                                          # interior mass of every interface grows by the fraction f (the outer one by the inflow)
        infl .+= δ[n+1] .* Yacc
        r = SM.advect_composition!(Y, Mc, δ, Yacc; order = 2); Mc = r.M; cflmax = max(cflmax, r.cfl)
    end
    tot1 = [sum(Mc[j] * Y[i, j] for j in 1:n) for i in 1:ns]
    cons5 = maximum(abs.(tot1 .- tot0 .- infl) ./ max.(tot0 .+ infl, 1e-300))
    mnew = [0.0; cumsum(Mc)]; xs = Y[SM.IH, :] ./ 0.75; ex = [clamp((ms - mnew[j]) / (mnew[j+1] - mnew[j]), 0.0, 1.0) + clamp((mnew[j+1] - Mtot0) / (mnew[j+1] - mnew[j]), 0.0, 1.0) for j in 1:n]
    ok5 = cons5 < 1e-13 && minimum(Y) >= 0 && cflmax > 100 && sum(abs.(xs .- ex) .* Mc) < 0.05 * 2
    @printf("(5) self-similar growth ×%.2f in %d steps: max donor CFL %.1e, conservation %.2e, min Y %.1e, L1 error of the matter step %.3e (mass units; zone mass at the step %.3e)\n",
            sum(Mc) / Mtot0, nst, cflmax, cons5, minimum(Y), sum(abs.(xs .- ex) .* Mc), Mz[findfirst(>(ms), mI) - 1])
    global pass &= ok5
end
println(pass ? "PASS" : "FAIL")
