# pilot_snap.jl <rundir> <prefix> : physics observables of a pilot from its snaps_*.jls (SNAPEVERY=5, CONVDIAG=1): prints "<prefix>.<key> <value>" lines.
# snaps: n, step span, t_sim_s (first to last snapshot), dt_min / dt_median over snapshots (as the M3b/M4 tables), R* / Mdot / M* at the end,
# fc_ftot_max (max F_conv/F_tot over zones and snapshots), pt_p_max & pt_p_ps_max (TDC only: P_t/P, max over snapshots / its post-shock-zone scalar),
# w_max (TDC), mach_max. Keys with no data are not printed.
using Serialization, Statistics, Printf
d, pre = ARGS[1], ARGS[2]
fs = sort(filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(d)))
S = Any[]
for f in fs
    open(joinpath(d, f)) do io
        while !eof(io); push!(S, deserialize(io)); end
    end
end
out(k, v) = (isfinite(v) && @printf("%s.%s %.10g\n", pre, k, v))
nz(v) = [isnan(q) ? 0.0 : q for q in v]
out("snaps", length(S))
isempty(S) && exit(0)
dts = [Float64(s.dt) for s in S]
out("step_first", S[1].step); out("step_last", S[end].step)
out("t_sim_s", S[end].t - S[1].t); out("dt_min", minimum(dts)); out("dt_median", median(dts)); out("dt_last", dts[end])
fc = -Inf; ptm = -Inf; ptps = -Inf; wmx = -Inf; mach = -Inf
for s in S
    c = s.cd; c === nothing && continue
    global fc = max(fc, maximum(nz(c.Fconv_over_Ftot)))
    global mach = max(mach, maximum(nz(c.mach_i)))
    hasproperty(c, :tdc_PtP) && (global ptm = max(ptm, maximum(nz(c.tdc_PtP))))
    hasproperty(c, :tdc_scal) && (global ptps = max(ptps, c.tdc_scal[2]))
    hasproperty(c, :tdc_w) && (global wmx = max(wmx, maximum(nz(c.tdc_w))))
end
out("fc_ftot_max", fc); out("pt_p_max", ptm); out("pt_p_ps_max", ptps); out("w_max", wmx); out("mach_max", mach)
s = S[end]; js = findfirst(<(-2e5), s.u)
if js !== nothing
    out("rstar_rsun", s.r[js] / 6.957e10)
    jj = js + 3
    jj <= length(s.ρ) && out("mdot_msun_yr", 4π * s.r[jj]^2 * s.ρ[jj] * abs(s.u[jj]) * 3.156e7 / 1.989e33)
end
out("mstar_msun", s.m[min(length(s.m), js === nothing ? length(s.m) : js)] / 1.989e33)   # (mass coordinate at R*: approximate; only a regression key)
