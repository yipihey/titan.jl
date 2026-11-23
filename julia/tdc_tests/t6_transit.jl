include(joinpath(@__DIR__, "harness.jl"))
open_log("t6_transit")
lp("T6: unstable layer of thickness L advected by a uniform settling flow u (Eulerian grid, mass flux ṁ = ρ u through the domain, w = 0 inflow)")
lp("    plane, ℓ = a = 1, sup_u = 1 (w_s = 1, N = √8), stable sup = -1 outside, α_t = 0, U = 0; TITAN sign: dmdt > 0 = inward (towards lower index)")
lp("    analytic transit (parcel enters the layer with w = 0): w(d) = w_s tanh(w_s d/(u ℓ)), d = distance travelled in the layer;")
lp("    measured at the downstream-most unstable interface, d = L − h (the last interface in the layer is d = (n_unstable) h from the stable edge interface)")
N = sqrt(8.0); L = 10.0; xa = 10.0; Xd = 30.0
par = TDCPar(α_t = 0.0)
function transit(ncells, Xs; dir = +1, cen = 0.0)
    h = L / ncells; nz = round(Int, Xd / h); xs = collect(range(0.0, Xd, length = nz + 1))
    u = L * N / Xs
    md = dir * 1.0 * u
    G = make_grid(xs; mdot = md)
    sup(x) = (xa < x < xa + L) ? 1.0 : -1.0
    B = make_back(G; a = 1.0, U = 0.0, sup = sup, ell = 1.0)
    dt = min(Xd / u, 5h / u); bud = 0.0      # CFL <= 5: with CFL >> 1 the semismooth Newton needs ~ncells iterations (the active set moves one cell per iteration), see doc
    cb = (s, t, w, wo, info) -> begin
        Q, _ = tdc_exchange(w, wo, G, B, dt, par; cen = cen); rs, sc = tdc_budget_check(Q, w, wo, G, B, dt, par; cen = cen)
        bud = max(bud, abs(rs) / max(sc, 1e-300))
    end
    w, info = evolve(zeros(nz + 1), G, B, dt, 4000, par; stop = 1e-13, cb = cb, cen = cen)
    ids = findall(k -> xa < xs[k] < xa + L, 1:nz+1)
    k = dir > 0 ? ids[1] : ids[end]           # downstream-most unstable interface
    d = dir > 0 ? (xa + L - xs[k]) : (xs[k] - xa)
    wan = tanh(d / u)
    return w[k], wan, info, bud, w, xs
end
res = Dict(); rows = []
lp(@sprintf("\n    %-8s %-8s %-9s | %s", "(L/u)N", "u/w_s", "w_an", "rel. deviation |w−w_an|/w_an for n cells across L = 10 | 30 | 100   (Newton failures, max budget residual)"))
let maxdev = Dict(10 => 0.0, 30 => 0.0, 100 => 0.0), maxbud = 0.0, nbad = 0
    for Xs in (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
        line = @sprintf("    %-8g %-8.3g ", Xs, L * N / Xs)
        devs = Float64[]; wan_ = NaN
        for nc in (10, 30, 100)
            wn, wan, info, bud, _, _ = transit(nc, Xs); wan_ = wan
            dev = abs(wn - wan) / wan; push!(devs, dev); maxdev[nc] = max(maxdev[nc], dev); maxbud = max(maxbud, bud); nbad += info.nbad + !info.steady
        end
        line *= @sprintf("%-9.5f | %.3e  %.3e  %.3e", wan_, devs...)
        lp(line); push!(rows, (Xs, L * N / Xs, wan_, devs...))
    end
    write_csv("t6_transit.csv", "tN_settle,u_over_ws,w_analytic,dev_n10,dev_n30,dev_n100", rows)
    lp(@sprintf("    max deviation: n=10 %.3e, n=30 %.3e, n=100 %.3e;  convergence order 10→100: %.2f", maxdev[10], maxdev[30], maxdev[100], log10(maxdev[10] / maxdev[100])))
    check("transit amplitude converges to w_s tanh(w_s L/(uℓ)) first order in h", maxdev[100] < 0.03 && maxdev[10] / maxdev[100] > 5 && nbad == 0)
    check("T8 energy budget (advective + diffusive e fluxes telescope) closes in every step of every T6 run", maxbud < 1e-13, @sprintf("max |resid|/scale = %.2e", maxbud))
end
# mirrored flow (outward, sign check of the donor logic): grid reversed and mirrored background must give the same amplitude
let
    wn1, wan, _, _, _, _ = transit(30, 3.0; dir = +1)
    h = L / 30; nz = round(Int, Xd / h); xs = collect(range(0.0, Xd, length = nz + 1)); u = L * N / 3.0
    G = make_grid(xs; mdot = -u)                       # outward flow, inflow at x = 0
    sup(x) = (Xd - xa - L < x < Xd - xa) ? 1.0 : -1.0   # mirrored layer
    B = make_back(G; a = 1.0, U = 0.0, sup = sup, ell = 1.0)
    w, info = evolve(zeros(nz + 1), G, B, 5h / u, 4000, par; stop = 1e-13)
    ids = findall(k -> Xd - xa - L < xs[k] < Xd - xa, 1:nz+1); wn2 = w[ids[end]]
    check("outward flow, mirrored layer: same amplitude as inward flow", abs(wn1 - wn2) < 1e-12, @sprintf("%.15f vs %.15f", wn1, wn2))
end
# profile through the layer, one case
let
    wn, wan, info, _, w, xs = transit(100, 1.0)
    lp("\n    profile along the flow, (L/u)N = 1, n = 100: x, w_num, w_an (tanh(d/u))")
    u = L * N / 1.0
    for x in 20.0:-2.0:8.0
        k = argmin(abs.(xs .- x)); d = x < xa + L && x > xa ? xa + L - xs[k] : NaN
        lp(@sprintf("      x = %5.1f  w = %.6f  w_an = %s", xs[k], w[k], isnan(d) ? "—" : @sprintf("%.6f", tanh(d / u))))
    end
end
lp("T6 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
