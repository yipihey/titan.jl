# stage2c_table.jl <rundir> [out.tsv]: the every-10th-step record of Task A (Stage 2c) from <rundir>/snaps_*.jls (CONVDIAG snapshots) and the DTERR lines of <rundir>/run_s*.out.
# core = interfaces with r < RCORE (1e11 cm); r_cb = outer end of the contiguous w > 0 run from the centre; layer = the sub-shock unstable layer of CONVDIAG (n_interfaces, max Fc/Ft).
using Serialization, Printf
dir = ARGS[1]; out = length(ARGS) > 1 ? ARGS[2] : joinpath(dir, "table.tsv"); RCORE = 1e11
snaps = Dict{Int,Any}()
for f in sort(filter(x -> occursin(r"^snaps_.*jls$", x), readdir(dir)))
    open(joinpath(dir, f)) do io
        while !eof(io); s = try deserialize(io) catch; break end; snaps[s.step] = s; end
    end
end
rx = r"DTERR jstep (\d+) t \S+ yr: dt (\S+) err (\S+) .*?rej (\d+) fail (\d+)"; rxe = r"ETOT (\d+) solves: G\(last \w+\) (\S+), sum (\S+), max\|G\| (\S+)"
dl = Dict{Int,Any}()
for f in sort(filter(x -> occursin(r"^run_s.*out$", x), readdir(dir)))
    for l in eachline(joinpath(dir, f))
        m = match(rx, l); m === nothing && continue; e = match(rxe, l)
        dl[parse(Int, m[1])] = (rej = parse(Int, m[4]), fail = parse(Int, m[5]), err = parse(Float64, m[3]), ns = e === nothing ? NaN : parse(Int, e[1]), sum = e === nothing ? NaN : parse(Float64, e[3]), maxG = e === nothing ? NaN : parse(Float64, e[4]))
    end
end
t0 = 8.752646409286999e12
open(out, "w") do io
    println(io, "step\tt_s\tdt_s\trej_err\tfail_newton\trho_c\tT_c\tw_max\tk_wmax\tr_wmax\tPtP_max\tFcFt_core_min\tFcFt_core_mean\tu_core_max\tmach_core\tr_cb\tcs_core\tt_sound\tt_turn_core\tlw_core_med\tlayer_n\tlayer_FcFt\tetot_sum\tetot_maxG")
    for st in sort(collect(keys(snaps)))
        s = snaps[st]; c = s.cd; n = length(s.ρ); r = s.r
        core = [j for j in 2:n if r[j] < RCORE]; isempty(core) && continue
        w = [isfinite(c.tdc_w[j]) ? c.tdc_w[j] : 0.0 for j in 1:n+1]
        jw = argmax(w[core]); kw = core[jw]
        PtP = maximum(filter(isfinite, c.tdc_PtP[core]); init = -Inf)
        F = filter(isfinite, c.Fconv_over_Ftot[core])
        ucore = maximum(abs.(s.u[2:maximum(core)])); cs = [sqrt(5 / 3 * s.pg[j] / s.ρ[j]) for j in 1:n]
        jcb = 2; while jcb + 1 <= n && w[jcb+1] > 0; jcb += 1; end
        csc = sum(cs[1:maximum(core)]) / maximum(core); rcb = r[jcb]
        tt = filter(isfinite, c.t_turn[core]); lay = c.layer
        d = get(dl, st, (rej = -1, fail = -1, err = NaN, ns = NaN, sum = NaN, maxG = NaN))
        @printf(io, "%d\t%.6f\t%.5e\t%d\t%d\t%.5g\t%.5g\t%.4e\t%d\t%.4e\t%.4e\t%.4g\t%.4g\t%.4e\t%.4g\t%.4e\t%.4e\t%.4g\t%.4g\t%.4g\t%.0f\t%.4g\t%.3e\t%.3e\n", st, s.t - t0, s.dt, d.rej, d.fail, s.ρ[1], s.T[1], w[kw], kw, r[kw], PtP,
                isempty(F) ? NaN : minimum(F), isempty(F) ? NaN : sum(F) / length(F), ucore, ucore / csc, rcb, csc, rcb / csc, isempty(tt) ? NaN : sum(tt) / length(tt), (lw = [c.lmix[j] / w[j] for j in core if w[j] > 1e3 && isfinite(c.lmix[j])]; isempty(lw) ? NaN : sort(lw)[(end + 1) ÷ 2]), lay[7], lay[5], d.sum, d.maxG)
    end
end
println("wrote ", out, " (", length(snaps), " snapshots)")
