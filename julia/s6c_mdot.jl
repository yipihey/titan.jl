# S6b Part 1 analysis: accretion history through the sink of popiii_m2c.jl (SINK=1) from its events file lib/popiii/s6c_<TAG>_events.csv.
#   julia s6c_mdot.jl <TAG> [<snapshot file of the collapse start>] [<TAG2> ...]       (plain Julia, no TITAN)
# Mdot(M) = Δm/Δt between successive mass shells crossing the sink radius (t shifted by the free-fall time from r_s to R_* = 0.05 AU, i.e. the time at which the shell
# reaches the protostar), compared with c_s³/G (T of the shell in the collapse snapshot and T at the sink) and with Omukai, Hosokawa & Yoshida 2010, Z = 0 (Fig. B16, read off).
# Output: table on stdout and lib/popiii/s6c_<TAG>_mdot.csv (M [Msun], Mdot [Msun/yr], t_arr [yr]) used as MDOTFILE of s6b_run.jl.
using Printf, Serialization
const G = 6.6743e-8; const Msun = 1.98847e33; const AU = 1.495978707e13; const YR = 3.15576e7; const kB = 1.380649e-16; const mH = 1.6726e-24
# OHY10 Fig. B16, Z = 0 (read off the figure, ±30 %): (M*/Msun, Mdot [Msun/yr])
const OHY = [(0.003, 0.5), (0.01, 0.2), (0.03, 0.12), (0.1, 0.085), (0.3, 0.068), (1.0, 0.055), (3.0, 0.045), (10.0, 0.035), (30.0, 0.03), (100.0, 0.02), (300.0, 0.011), (1000.0, 0.004)]
function ohy(M)
    M <= OHY[1][1] && return OHY[1][2]; M >= OHY[end][1] && return NaN
    i = findlast(p -> p[1] <= M, OHY); w = log(M / OHY[i][1]) / log(OHY[i+1][1] / OHY[i][1])
    exp((1 - w) * log(OHY[i][2]) + w * log(OHY[i+1][2]))
end
mu_of(fH2, xe) = 1.333 / (1.0833 - 0.5fH2 + xe)
cs3G(T, mu) = (kB * T / (mu * mH))^1.5 / G * YR / Msun         # Msun/yr
"free-fall time [yr] from r_s (inward speed u_s [cm/s] > 0) to r_end onto the point mass m [g]"
function tfall(rs, us, m, rend)
    n = 400; xs = exp.(range(log(rs), log(rend); length = n)); t = 0.0
    for i in 1:n-1
        v(r) = sqrt(max(us^2 + 2G * m * (1 / r - 1 / rs), 1e-10))
        t += (xs[i] - xs[i+1]) / (0.5 * (v(xs[i]) + v(xs[i+1])))
    end
    t / YR
end
function analyse(tag; snap = nothing, rend = 0.05 * AU, quiet = false)
    L = readlines(joinpath(@__DIR__, "..", "lib", "popiii", "s6c_$(tag)_events.csv"))[2:end]
    ev = [parse.(Float64, split(l, ",")) for l in L if !isempty(strip(l))]
    sort!(ev, by = e -> e[3])
    # columns: step, t_cross_yr, m_interface, r_s_AU, u_kms, rho, T, fH2, xe, e, p, cs_kms, rho_up, T_up, u_up_kms
    m = [e[3] for e in ev]; tc = [e[2] for e in ev]; rs = [e[4] * AU for e in ev]; us = [-e[5] * 1e5 for e in ev]
    ta = [tc[i] + tfall(rs[i], max(us[i], 1e5), m[i] * Msun, rend) for i in eachindex(ev)]
    Ts = [e[7] for e in ev]; mu = [mu_of(e[8], e[9]) for e in ev]
    # pre-collapse temperature of the shell (snapshot)
    Tsnap = fill(NaN, length(m)); musnap = fill(NaN, length(m))
    if snap !== nothing
        s = deserialize(snap)[end]; ms = s.m ./ Msun
        for i in eachindex(m)
            j = clamp(searchsortedfirst(ms, m[i]) - 1, 1, length(s.T)); Tsnap[i] = s.T[j]; musnap[i] = mu_of(2s.y[j][3], s.y[j][1])
        end
    end
    # Mdot = ΔM/Δt over windows of events spanning Δ ln M ≥ WIN (default 0.5): the shells cross the sink radius one zone at a time and the crossing times of single shells scatter
    # (alternating Mdot of ×4 between neighbours at 5 % zone masses); the mean rate is the ratio of the sums, not the mean of the single-pair rates
    WIN = 0.5
    rows = []
    for i in 1:length(m)-1
        j1 = i; j2 = i + 1
        while log(m[j2] / m[j1]) < WIN && (j1 > 1 || j2 < length(m))
            j1 > 1 && (j1 -= 1); log(m[j2] / m[j1]) >= WIN && break
            j2 < length(m) && (j2 += 1)
        end
        dt = ta[j2] - ta[j1]; dt > 0 || continue
        push!(rows, (; M = sqrt(m[j1] * m[j2]), Mdot = (m[j2] - m[j1]) / dt, t = 0.5 * (ta[j1] + ta[j2]), Tsink = Ts[i], Tsnap = Tsnap[i], musnap = musnap[i], rs = rs[i] / AU, lo = m[j1], hi = m[j2]))
    end
    open(joinpath(@__DIR__, "..", "lib", "popiii", "s6c_$(tag)_mdot.csv"), "w") do io
        println(io, "M_Msun,Mdot_Msun_yr,t_arr_yr,T_sink,T_snapshot,r_s_AU,M_lo,M_hi")
        for r in rows; println(io, r.M, ",", r.Mdot, ",", r.t, ",", r.Tsink, ",", r.Tsnap, ",", r.rs, ",", r.lo, ",", r.hi); end
    end
    open(joinpath(@__DIR__, "..", "lib", "popiii", "s6c_$(tag)_mdot_smooth.csv"), "w") do io      # input of s6b_run.jl (MDOTFILE); one row per ~0.15 in ln M
        println(io, "M_Msun,Mdot_Msun_yr")
        last = -Inf
        for r in rows; log(r.M) - last >= 0.15 && (println(io, r.M, ",", r.Mdot); last = log(r.M)); end
    end
    if !quiet
        @printf("%s: %d events, M = %.4g–%.4g Msun, shells arrive over %.4g yr\n", tag, length(m), m[1], m[end], ta[end] - ta[1])
        @printf("%10s %12s %10s %9s %9s %12s %8s %10s %8s\n", "M [Msun]", "Mdot [Msun/yr]", "t_arr [yr]", "T_snap", "T_sink", "cs3G(Tsnap)", "phi", "OHY10 Z=0", "ratio")
        step = max(1, length(rows) ÷ 28)
        for r in rows[1:step:end]
            c3 = cs3G(r.Tsnap, r.musnap)
            @printf("%10.4g %12.4e %10.4g %9.0f %9.0f %12.3e %8.1f %10.3e %8.2f\n", r.M, r.Mdot, r.t - rows[1].t, r.Tsnap, r.Tsink, c3, r.Mdot / c3, ohy(r.M), r.Mdot / ohy(r.M))
        end
    end
    rows
end
if abspath(PROGRAM_FILE) == @__FILE__
    tag = ARGS[1]; snap = length(ARGS) > 1 ? ARGS[2] : nothing
    analyse(tag; snap)
end
