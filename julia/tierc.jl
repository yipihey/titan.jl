# Tier-(c) acceptance: physical observables of a candidate run vs the Fortran
# reference run (both recorded per step as Metrics.Snap by phase2_run.jl).
# usage: julia tierc.jl [candidate_mode=julia]  -> prints tables, writes
#        reference_runs/metrics_reference.csv and lib/p2_runs/tierc_<mode>.txt
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const ROOT = normpath(joinpath(@__DIR__, ".."))
cand = isempty(ARGS) ? "julia" : ARGS[1]
# reference run mode (env TIERC_REF): "fortran" = TITAN 1.2 bitwise (original bookkeeping),
# "std" = Julia driver, all-Fortran routines, clock/energy fixes ON (the Phase-4 reference)
const REF = get(ENV, "TIERC_REF", "fortran")
out = open(joinpath(ROOT, "lib", "p2_runs", "tierc_$(cand)_vs_$(REF).txt"), "w")
refcsv = open(joinpath(ROOT, "reference_runs", REF == "fortran" ? "metrics_reference.csv" : "metrics_reference_$(REF).csv"), "w")
println(refcsv, "problem,time,metric,value")
say(x...) = (println(x...); println(out, x...))

load(name, mode) = deserialize(joinpath(ROOT, "lib", "p2_runs", "$(name)_$mode", "snaps.jls"))

# acceptance thresholds (see PORTING_NOTES §P2.3)
const POS_FRAC = 5e-3      # positions: |Δr| ≤ 0.5% of the length scale L
const TEMP_REL = 1e-2      # temperatures, plateau values: ≤ 1% relative
const EERR_ABS = 1e-3      # energy-error metric: |Δ| ≤ 1e-3 (absolute)
const STEP_REL = 3e-2      # step count: ≤ 3%

function series(snaps, f)
    s0 = snaps[1]
    [merge(f(s), (E_err = Metrics.energy_err(s, s0),)) for s in snaps]
end

function compare(name, f, times, posk, Lfun; exact = nothing)
    R = load(name, REF); C = load(name, cand)
    mr = series(R, f); mc = series(C, f)
    say("\n=== $name: steps ref = $(length(R)), $cand = $(length(C)) ",
        @sprintf("(Δ = %+.2f%%)  %s", 100 * (length(C) / length(R) - 1),
                 abs(length(C) / length(R) - 1) <= STEP_REL ? "PASS" : "FAIL"))
    say(@sprintf("    final time ref = %.6e, %s = %.6e", R[end].time, cand, C[end].time))
    nfail = abs(length(C) / length(R) - 1) <= STEP_REL ? 0 : 1
    for t in times
        a = Metrics.at_time(mr, t); b = Metrics.at_time(mc, t)
        (a === nothing || b === nothing) && continue
        L = Lfun(a)
        say(@sprintf("  t = %.4g   (length scale L = %.4g)", t, L))
        for k in keys(a)
            k == :time && continue
            va, vb = a[k], b[k]
            println(refcsv, "$name,$t,$k,$(va)")
            if k in posk
                d = abs(vb - va) / L; ok = d <= POS_FRAC; crit = "|Δ|/L"
            elseif k == :E_err
                d = abs(vb - va); ok = d <= EERR_ABS; crit = "|Δ|"
            elseif k in (:width,)
                d = abs(vb - va) / L; ok = d <= POS_FRAC; crit = "|Δ|/L"
            else
                d = abs(vb - va) / max(abs(va), 1e-300); ok = d <= TEMP_REL; crit = "rel"
            end
            ok || (nfail += 1)
            ex = exact !== nothing && haskey(exact(t), k) ? @sprintf("  exact %.6g", exact(t)[k]) : ""
            say(@sprintf("    %-11s ref %.8g   %s %.8g   %s = %.2e  %s%s", k, va, cand, vb, crit, d,
                         ok ? "ok" : "FAIL", ex))
        end
    end
    say("  => $name: ", nfail == 0 ? "ACCEPT" : "$nfail metric(s) outside tolerance")
end

compare("SOD", Metrics.sod_metrics, (0.1, 0.2), (:r_shock, :r_contact, :r_head, :r_tail),
        a -> 1.0; exact = t -> Metrics.sod_exact(t))
compare("SUBC", Metrics.shock_metrics, (5.4e3, 1.7e4, 2.8e4, 3.8e4, 6.0e4),
        (:r_shock, :r_spike, :precursor), a -> a.width)
compare("SUPC", Metrics.shock_metrics, (4.0e3, 7.5e3, 1.3e4),
        (:r_shock, :r_spike, :precursor), a -> a.width)
close(out); close(refcsv)
