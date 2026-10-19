# Regression check for the Population III project: tier-(c) observables of the pure-Julia runs
# (lib/p2_runs/<NAME>_pure/snaps.jls, written by `phase6_run.jl pure NAME`) against the fixed Fortran reference
# reference_runs/metrics_reference_fixed.csv (PORTING_NOTES P4.1), with tierc.jl's tolerances.
# usage: julia --project=env popiii_regress.jl [SOD SUBC SUPC]   → prints a table; exit code 1 on any FAIL
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const ROOT = normpath(joinpath(@__DIR__, ".."))
const POS_FRAC = 5e-3; const TEMP_REL = 1e-2; const STEP_REL = 3e-2
# acceptance times as in tierc.jl; other reference times (SOD 0.5: after the wall reflection) are printed as "info"
const ACCEPT_T = Dict("SOD" => (0.1, 0.2), "SUBC" => (5.4e3, 1.7e4, 3.8e4, 6.0e4), "SUPC" => (4.0e3, 7.5e3, 1.3e4))
const REFSTEPS = Dict("SOD" => 452, "SUBC" => 587, "SUPC" => 713)       # fixed Fortran reference (P4.1, P6)
const OBS = Dict("SOD" => (Metrics.sod_metrics, (:r_shock, :r_contact, :r_head, :r_tail)),
                 "SUBC" => (Metrics.shock_metrics, (:r_shock, :r_spike, :precursor, :width)),
                 "SUPC" => (Metrics.shock_metrics, (:r_shock, :r_spike, :precursor, :width)))
ref = Dict{Tuple{String,Float64,Symbol},Float64}()
for l in readlines(joinpath(ROOT, "reference_runs", "metrics_reference_fixed.csv"))[2:end]
    a = split(l, ","); ref[(a[1], parse(Float64, a[2]), Symbol(a[3]))] = parse(Float64, a[4])
end
nfail = 0
for name in (isempty(ARGS) ? ["SOD", "SUBC", "SUPC"] : ARGS)
    global nfail
    S = deserialize(joinpath(ROOT, "lib", "p2_runs", "$(name)_pure", "snaps.jls"))
    f, posk = OBS[name]; ms = [f(s) for s in S]
    ns = length(S) - 1; dn = ns / REFSTEPS[name] - 1
    okn = abs(dn) <= STEP_REL; okn || (nfail += 1)
    @printf("=== %s: steps %d vs reference %d (%+.2f%%) %s\n", name, ns, REFSTEPS[name], 100dn, okn ? "ok" : "FAIL")
    for t in sort(unique(k[2] for k in keys(ref) if k[1] == name))
        b = Metrics.at_time(ms, t); b === nothing && continue
        L = name == "SOD" ? 1.0 : ref[(name, t, :width)]
        for k in keys(b)
            (k == :time || !haskey(ref, (name, t, k))) && continue
            va = ref[(name, t, k)]; vb = b[k]
            d, crit, tol = k in posk ? (abs(vb - va) / L, "|Δ|/L", POS_FRAC) : (abs(vb - va) / max(abs(va), 1e-300), "rel", TEMP_REL)
            ok = d <= tol; acc = t in ACCEPT_T[name]; (ok || !acc) || (nfail += 1)
            @printf("  t = %-8.4g %-10s ref %-14.8g run %-14.8g %s = %.2e %s\n", t, k, va, vb, crit, d,
                    !acc ? (ok ? "info (within tol)" : "info (outside tol)") : ok ? "ok" : "FAIL")
        end
    end
end
println(nfail == 0 ? "ACCEPT" : "$nfail FAIL")
exit(nfail == 0 ? 0 : 1)
