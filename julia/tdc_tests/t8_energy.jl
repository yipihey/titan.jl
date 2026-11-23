include(joinpath(@__DIR__, "harness.jl"))
open_log("t8_energy")
lp("T8: discrete turbulent-energy budget  Σ_k Q_k Δt = Δ Σ_k m e + Δt Σ_j Pt_j dVdt_j  (e = w², closed domain, w = 0 boundaries)")
using Random
rng = MersenneTwister(20261122)
maxrel = 0.0; maxgas = 0.0
for (μ, nz, cen) in ((0, 40, 0.0), (2, 57, 0.0), (2, 57, 0.4), (2, 31, 1.0), (0, 25, 0.7))
    r = sort(unique(vcat(1.0, 1.0 .+ cumsum(0.2 .+ rand(rng, nz)))))     # non-uniform interfaces
    nz = length(r) - 1
    rho = 0.5 .+ rand(rng, nz)
    G0 = make_grid(r; μ = μ, rho = rho)
    # break continuity on purpose (mzn ≠ mzo), random mass fluxes of both signs, random volume rates
    G = TDCGeom(nz, G0.r, G0.Az, G0.mzn .* (1 .+ 0.1 .* randn(rng, nz)), G0.mzo, 3 .* randn(rng, nz + 1) .* (μ == 2 ? 5.0 : 1.0), randn(rng, nz))
    B = make_back(G; a = 1.0 .+ rand(rng, nz + 1), U = rand(rng, nz + 1), sup = randn(rng, nz + 1), ell = 0.5 .+ rand(rng, nz + 1), rho = rho, ellz = 0.5 .+ rand(rng, nz))
    par = TDCPar(α_t = 0.3, α_p = 2 / 3)
    w = abs.(randn(rng, nz + 1)); w[1] = 0; w[end] = 0
    wo = abs.(randn(rng, nz + 1)); wo[1] = 0; wo[end] = 0
    dt = 0.07
    Q, parts = tdc_exchange(w, wo, G, B, dt, par; cen = cen)
    rs, sc = tdc_budget_check(Q, w, wo, G, B, dt, par; cen = cen)
    gas = tdc_gas_exchange(Q)
    gres = abs(sum(gas) + sum(Q)) / sum(abs, Q)
    global maxrel = max(maxrel, abs(rs) / sc); global maxgas = max(maxgas, gres)
    lp(@sprintf("    μ=%d nz=%d cen=%.1f: |resid|/scale = %.2e  (scale %.3g),  |Σ gas + Σ Q|/Σ|Q| = %.2e", μ, nz, cen, abs(rs) / sc, sc, gres))
end
check("random non-uniform grids, continuity-violating masses, both flow signs, μ = 0, 2, cen > 0: budget closes to round-off", maxrel < 1e-14 && maxgas < 1e-14, @sprintf("max %.2e / gas %.2e", maxrel, maxgas))
# the steady local state has Q = 0 (production = dissipation) and a growing one Q > 0
let nz = 12; G = make_grid(collect(range(0.0, 12.0, length = nz + 1))); par = TDCPar(α_t = 0.0)
    B = make_back(G; a = 1.0, U = 0.4, sup = 0.8, ell = 1.0); wref = tdc_mlt_w(1.0, 0.4, 0.8)
    w = fill(wref, nz + 1); w[1] = 0; w[end] = 0
    Q, _ = tdc_exchange(w, w, G, B, 0.1, par)
    lp(@sprintf("    steady MLT state: max |Q_k| interior = %.2e (should vanish: no storage, no flux divergence)", maximum(abs.(Q[3:end-2]))))
    Qw, _ = tdc_exchange(w .* 1.1, w, G, B, 0.1, par)
    check("steady uniform MLT state: Q = 0 in the interior; growing w: Q > 0", maximum(abs.(Q[3:end-2])) < 1e-13 && Qw[6] > 0)
end
# the tridiagonal-Jacobian assumption (3-colour ForwardDiff) vs the dense ForwardDiff Jacobian
let nz = 30, rr = collect(range(1.0, 4.0, length = 31)); G = make_grid(rr; μ = 2, mdot = 2.0); par = TDCPar(α_t = 0.2)
    B = make_back(G; a = 1.0, U = 0.1, sup = 0.5, ell = 0.3); wo = 0.1 .* rand(rng, nz + 1); wo[1] = 0; wo[end] = 0
    w = 0.1 .* rand(rng, nz + 1)
    f!(φ, x) = tdc_rows!(φ, x, wo, G, B, 0.2, par)
    Jd = ForwardDiff.jacobian(x -> (φ = zeros(eltype(x), nz + 1); f!(φ, x); φ), w)
    Jt = Matrix(jac_tridiag(f!, w))
    check("3-colour tridiagonal Jacobian == dense ForwardDiff Jacobian (rows are 3-point)", maximum(abs.(Jd .- Jt)) <= 1e-12 * maximum(abs.(Jd)), @sprintf("max diff %.2e", maximum(abs.(Jd .- Jt))))
end
lp("T8 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
