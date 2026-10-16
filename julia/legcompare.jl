# Phase 4, item 1: observables for all 9 benchmark legs, orig (TITAN 1.2 bitwise)
# vs fixed (clock + energy-booking fixes).  Writes
#   reference_runs/metrics_reference_fixed.csv  (new 'fixed' reference)
#   reference_runs/metrics_reference_orig_all.csv (all legs, original bookkeeping)
#   lib/p4_legs/legcompare.txt
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const ROOT = normpath(joinpath(@__DIR__, ".."))
const D = joinpath(ROOT, "lib", "p4_legs")
load(name, mode) = deserialize(joinpath(D, "$(name)_$mode", "snaps.jls"))
tosnap(s) = Metrics.Snap(s.jstep, s.t, s.r, s.rho, s.T, s.Er, s.p, s.u, s.tetot, s.tee, s.tew)
rc(s) = 0.5 .* (s.r[1:end-1] .+ s.r[2:end])
"steepest drop of y per unit length -> refined crossing position"
function steep(x, y)
    k = argmax([(y[i] - y[i+1]) / (x[i+1] - x[i]) for i in 1:length(y)-1])
    lev = 0.5 * (maximum(y[max(1, k-2):k]) + minimum(y[k+1:min(end, k+3)]))
    for i in max(1, k-3):min(length(y)-1, k+3)
        (y[i] - lev) * (y[i+1] - lev) <= 0 && y[i] != y[i+1] &&
            return x[i] + (x[i+1] - x[i]) * (lev - y[i]) / (y[i+1] - y[i])
    end
    x[k]
end
fourpi(s) = s.mu == 2 ? 4π : 1.0
Lout(s) = fourpi(s) * (s.mu == 2 ? s.r[end]^2 : 1.0) * s.F[end]
obs = Dict(
    "SOD"  => s -> Metrics.sod_metrics(tosnap(s)),
    "SUBC" => s -> Metrics.shock_metrics(tosnap(s)),
    "SUPC" => s -> Metrics.shock_metrics(tosnap(s)),
    "BLST" => s -> (k = argmax(s.rho); kp = argmax(s.p);
                    (; time = s.t, rho_max = s.rho[k], r_rhomax = rc(s)[k], p_max = s.p[kp], r_pmax = rc(s)[kp])),
    "SDOV" => s -> (; time = s.t, r_shock = steep(rc(s), s.rho), rho_max = maximum(s.rho)),
    "HEAT" => s -> (; time = s.t, L_out = Lout(s), T_in = s.T[1], T_mid = s.T[end ÷ 2]),
    "COOL" => s -> (; time = s.t, L_out = Lout(s), T_in = s.T[1], T_mid = s.T[end ÷ 2]),
    "RLST" => s -> (; time = s.t, r_shock = steep(rc(s), s.rho), T_max = maximum(s.T), L_out = Lout(s)),
    "RLST_leg2" => s -> (; time = s.t, r_shock = steep(rc(s), s.rho), T_max = maximum(s.T), L_out = Lout(s)))
times = Dict("SOD" => (0.1, 0.2, 0.5), "SUBC" => (5.4e3, 1.7e4, 3.8e4, 6.0e4), "SUPC" => (4.0e3, 7.5e3, 1.3e4),
             "BLST" => (0.016, 0.028, 0.038), "SDOV" => (5.0e4, 1.0e5, 1.9e5),
             "HEAT" => (1e13, 1e14, 1e15, 2.5e15), "COOL" => (1e15, 1e16, 1e17, 2e17),
             "RLST" => (1e-21, 5e-20), "RLST_leg2" => (1e2, 1e3, 5e3))
order = ["SOD", "SUBC", "SUPC", "BLST", "SDOV", "HEAT", "COOL", "RLST", "RLST_leg2"]
function at(S, f, t0)
    ts = [s.t for s in S]
    i = findlast(<=(t0), ts)
    (i === nothing || i >= length(S)) && return nothing
    j = i + 1; while j <= length(S) && S[j].t == S[i].t; j += 1; end
    j > length(S) && return nothing
    a, b = f(S[i]), f(S[j]); w = (t0 - S[i].t) / (S[j].t - S[i].t)
    NamedTuple{keys(a)}(Tuple(k == :time ? t0 : a[k] + w * (b[k] - a[k]) for k in keys(a)))
end
"TITAN's energy diagnostic: tetot drift since step 1, relative to the largest ledger term"
ediag(S) = (s = S[end]; s1 = S[2]; (s.tetot - s1.tetot) / max(abs(s.tee), abs(s.tew), abs(s.tel), abs(s1.tee), floatmin()))
out = open(joinpath(D, "legcompare.txt"), "w")
cf = open(joinpath(ROOT, "reference_runs", "metrics_reference_fixed.csv"), "w")
co = open(joinpath(ROOT, "reference_runs", "metrics_reference_orig_all.csv"), "w")
for io in (cf, co); println(io, "leg,time,metric,value"); end
say(x...) = (println(x...); println(out, x...))
for name in order
    O = load(name, "orig"); F = load(name, "fixed")
    say(@sprintf("\n=== %s: steps orig %d / fixed %d; t_end orig %.6e / fixed %.6e; Σdt %.6e; retried steps %d; TITAN energy diagnostic orig %+.3e / fixed %+.3e",
                 name, length(O) - 1, length(F) - 1, O[end].t, F[end].t, sum(s.dt for s in O[2:end]),
                 count(i -> O[i].t == O[i-1].t, 2:length(O)), ediag(O), ediag(F)))
    for t in times[name]
        a = at(O, obs[name], t); b = at(F, obs[name], t)
        (a === nothing || b === nothing) && (say("  t=$t: outside run"); continue)
        for k in keys(a)
            k in (:time, :E_err) && continue
            println(co, "$name,$t,$k,$(a[k])"); println(cf, "$name,$t,$k,$(b[k])")
            say(@sprintf("  t=%-9.3g %-11s orig %+.8e  fixed %+.8e  rel shift %+.2e", t, k, a[k], b[k],
                         (b[k] - a[k]) / max(abs(a[k]), floatmin())))
        end
    end
end
close(out); close(cf); close(co)
