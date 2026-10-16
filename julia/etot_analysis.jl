# Phase 4 item 2 analysis: ledgers (std vs etot), SOD wall reflection vs the exact reflected
# state, early shock formation (SUBC/SUPC).  usage: julia etot_analysis.jl [tags...]
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const ROOT = normpath(joinpath(@__DIR__, ".."))
tags = isempty(ARGS) ? ["std", "etot"] : ARGS
ld(name, tag, f) = deserialize(joinpath(ROOT, "lib", "p2_runs", "$(name)_$tag", f))
out = open(joinpath(ROOT, "lib", "p4_runs", "etot_analysis.txt"), "w")
say(x...) = (println(x...); println(out, x...))
for name in ("SOD", "SUBC", "SUPC"), tag in tags
    isfile(joinpath(ROOT, "lib", "p2_runs", "$(name)_$tag", "ledger.jls")) || continue
    L = ld(name, tag, "ledger.jls")
    cum = cumsum([l.imb for l in L]); W = cumsum([l.fout - l.fin for l in L]); E0 = L[1].Enew - L[1].dE
    sc = max(abs(E0), maximum(abs, W))
    cuts = haskey(L[1], :cuts) ? sum(l.cuts for l in L) : -1
    say(@sprintf("%-5s %-8s steps %4d  t_end %.6e  Newton it/step %.2f  dt-cut retries %3d  ΣΔE_imbalance/scale %+.2e  max|step|/scale %.1e  max dual-switch zones %d",
                 name, tag, length(L), L[end].t, sum(l.nit for l in L) / length(L), cuts, cum[end] / sc,
                 maximum(abs(l.imb) for l in L) / sc, maximum(l.nsw for l in L)))
end
# ---------------- SOD: exact state behind the shock reflected from the right wall
γ = 1.4; ρ3, u3, p3 = 0.265574, 0.927453, 0.30313
A = 2 / ((γ + 1) * ρ3); B = (γ - 1) / (γ + 1) * p3
f(p) = u3 - (p - p3) * sqrt(A / (p + B))
pw = let lo = p3, hi = 10p3
    for _ in 1:200; m = (lo + hi) / 2; f(m) > 0 ? (lo = m) : (hi = m); end
    (lo + hi) / 2
end; r = pw / p3; ρw = ρ3 * (r + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * r + 1)
Tw_rel = (pw / ρw) / (p3 / ρ3)
say(@sprintf("\nSOD reflected-shock exact state (wall at rest): p_w = %.5f  ρ_w = %.5f  (p/ρ)_w / (p/ρ)_3 = %.5f", pw, ρw, Tw_rel))
for tag in vcat(["fortran"], tags)
    p = joinpath(ROOT, "lib", "p2_runs", "SOD_$tag", "snaps.jls"); isfile(p) || continue
    S = deserialize(p)
    for t0 in (0.32, 0.36)
        i = argmin([abs(s.time - t0) for s in S]); s = S[i]; n = length(s.rho)
        rc = 0.5 .* (s.r[1:end-1] .+ s.r[2:end])
        # plateau = zones within 0.02 of the wall excluding the 2 wall zones
        pl = [k for k in 1:n-2 if rc[k] > s.r[end] - 0.03]
        ρp = isempty(pl) ? NaN : sum(s.rho[pl]) / length(pl); pp = isempty(pl) ? NaN : sum(s.p[pl]) / length(pl)
        say(@sprintf("  %-8s t=%.4f  wall zone: ρ/ρ_w %.4f  p/p_w %.4f  T-ratio/exact %.4f | next zone ρ/ρ_w %.4f | plateau(0.03) ρ/ρ_w %.4f p/p_w %.4f  (wall zone Δr %.4f)",
                     tag, s.time, s.rho[n] / ρw, s.p[n] / pw, (s.p[n] / s.rho[n]) / (pw / ρw), s.rho[n-1] / ρw, ρp / ρw, pp / pw,
                     s.r[end] - s.r[end-1]))
    end
    s = S[end]; n = length(s.rho)
    say(@sprintf("  %-8s t_end=%.4f  wall-zone T/T(5 zones in) %.4f   max|u| %.3e", tag, s.time,
                 (s.p[n] / s.rho[n]) / (s.p[n-5] / s.rho[n-5]), maximum(abs, s.u)))
end
# ---------------- shock formation (first steps)
for name in ("SUBC", "SUPC")
    say("\n$name early phase (shock formation):")
    for tag in tags
        p = joinpath(ROOT, "lib", "p2_runs", "$(name)_$tag", "ledger.jls"); isfile(p) || continue
        L = ld(name, tag, "ledger.jls")
        for nst in (10, 50, 100, 200)
            say(@sprintf("  %-8s after %3d steps: t = %.4e   Newton iterations %4d", tag, nst, L[nst].t, sum(l.nit for l in L[1:nst])))
        end
    end
end
close(out)
