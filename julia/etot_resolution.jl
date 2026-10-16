# Item 2: is the std-vs-etot difference a discretisation difference? Compare both at 100
# zones against the 200-zone runs of each formulation (SUBC: shock metrics; SOD: vs exact).
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const R = joinpath(@__DIR__, "..", "lib", "p2_runs")
ld(t) = deserialize(joinpath(R, t, "snaps.jls"))
out = open(joinpath(@__DIR__, "..", "lib", "p4_runs", "etot_resolution.txt"), "w")
say(x...) = (println(x...); println(out, x...))
runs = ["SUBC_std", "SUBC_etot", "SUBC_std_n200", "SUBC_etot_n200"]
M = Dict(r => [Metrics.shock_metrics(s) for s in ld(r)] for r in runs)
for t0 in (5.4e3, 1.7e4, 3.8e4, 6.0e4)
    m = Dict(r => Metrics.at_time(M[r], t0) for r in runs)
    L = m["SUBC_std_n200"].width
    say(@sprintf("SUBC t=%.2g (L=width n200=%.3e)", t0, L))
    for k in (:r_shock, :rho_post, :rho_pre, :T_post, :T_pre, :Trad_post, :precursor)
        ref = 0.5 * (m["SUBC_std_n200"][k] + m["SUBC_etot_n200"][k])
        nrm = k in (:r_shock, :precursor) ? L : abs(ref)
        say(@sprintf("   %-10s (x − mean n200)/%s:  std100 %+.2e  etot100 %+.2e  std200 %+.2e  etot200 %+.2e", k,
                     k in (:r_shock, :precursor) ? "L" : "|x|", (m["SUBC_std"][k] - ref) / nrm, (m["SUBC_etot"][k] - ref) / nrm,
                     (m["SUBC_std_n200"][k] - ref) / nrm, (m["SUBC_etot_n200"][k] - ref) / nrm))
    end
end
for r in ("SOD_std", "SOD_etot", "SOD_std_n200", "SOD_etot_n200")
    ms = [Metrics.sod_metrics(s) for s in ld(r)]
    for t0 in (0.1, 0.2)
        m = Metrics.at_time(ms, t0); ex = Metrics.sod_exact(t0)
        say(@sprintf("%-14s t=%.1f  Δ vs exact: r_shock %+.2e r_contact %+.2e rho_2 %+.2e rho_3 %+.2e p_star %+.2e", r, t0,
                     m.r_shock - ex.r_shock, m.r_contact - ex.r_contact, m.rho_2 / ex.rho_2 - 1, m.rho_3 / ex.rho_3 - 1, m.p_star / ex.p_star - 1))
    end
end
close(out)
