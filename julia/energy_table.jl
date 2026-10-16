# Fixed-time comparison of the independent energy residual across the sweep.
# usage: julia energy_table.jl  -> lib/p3_runs/energy_table.txt
using Serialization, Printf
const D = joinpath(@__DIR__, "..", "lib", "p3_runs")
function resid_series(S; w)                 # θ-weighted (w) work, +q, KE_zone; true step dt
    n = length(S); W = zeros(n); L = zeros(n)
    for i in 2:n
        x(f) = w * getfield(S[i], f) + (1 - w) * getfield(S[i-1], f)
        W[i] = W[i-1] + S[i].dt * ((x(:Win_p) + x(:Win_r) + x(:Win_q)) - (x(:Wout_p) + x(:Wout_r) + x(:Wout_q)))
        L[i] = L[i-1] + S[i].dt * (x(:Lout) - x(:Lin))
    end
    E(s) = s.Eg + s.Er + s.KEz
    ttrue = cumsum([0.0; [S[i].dt for i in 2:n]])            # physical time incl. retried steps
    res = [E(S[i]) - E(S[1]) - W[i] + L[i] for i in 1:n]
    ttrue, res, W, [E(s) for s in S]
end
at(t, x, y) = (i = searchsortedlast(t, x); i < 1 || i >= length(t) ? NaN : y[i] + (y[i+1] - y[i]) * (x - t[i]) / (t[i+1] - t[i]))
cases = Dict("SUBC" => ((5.4e3, 1.7e4, 3.8e4, 8.0e4), 1.0), "SUPC" => ((4.0e3, 7.5e3, 1.3e4), 1.0),
             "SOD" => ((0.2, 0.5, 1.1), 0.55))
open(joinpath(D, "energy_table.txt"), "w") do io
    for (name, (times, θ)) in sort(collect(cases))
        files = filter(f -> startswith(f, "ledger_$(name)_") && endswith(f, ".jls"), readdir(D))
        println(io, "\n$name: residual (ΔE − W + L) normalised by |W(t)| (piston problems) or E0 (SOD), at physical time t")
        @printf(io, "  %-14s %6s %10s", "run", "steps", "t_end")
        for t in times; @printf(io, "   t=%-9.3g", t); end
        println(io)
        for f in sort(files)
            tag = f[length("ledger_$(name)_")+1:end-4]
            S = deserialize(joinpath(D, f))
            th = occursin("thet", tag) ? parse(Float64, replace(split(tag, "thet")[2], "075" => "0.75", "1" => "1.0")) : θ
            tt, res, W, E = resid_series(S; w = th)
            @printf(io, "  %-14s %6d %10.4g", tag, length(S) - 1, tt[end])
            for t in times
                sc = name == "SOD" ? E[1] : abs(at(tt, t, W))
                @printf(io, "   %+.3e", at(tt, t, res) / sc)
            end
            println(io)
        end
    end
end
print(read(joinpath(D, "energy_table.txt"), String))
