# S7b protostar side of the coupled protostellar feedback: hooks for s6b_run.jl (included when S7B_STAR=1). Counterpart of s7b_env.jl (see the header there).
#
# Phase :pre — the protostar runs as in S6b Part 3 (Mdot from the MDOTFILE table of the collapse) until its mass reaches S7B_MSYNC [Msun]. Then (barrier 0) it writes its record star_00000.txt and waits for
# env_00000.txt. Phase :coupled — Mdot is constant over each interval between barriers: the envelope's rate over its last window (env record) plus a correction that closes the difference between the
# mass accreted by the envelope's sink (Macc) and by the star (M − M_sync) over S7B_NSLAVE [15] intervals (limited to ±20 % of the rate, and the rate changes by at most S7B_SLEW [5 %] per interval); the time step is cut so that the star lands on the barrier of the common clock,
# where it writes star_<j>.txt (M, R*, L*, L_acc, L_∞, T_ph, r_ph, T_s, τ, T_c, … in cgs) and reads env_<j>.txt. The run ends when the envelope reports done = 1.
module S7bStar
using Printf
const S = Main.S7bX
const SR = 6.957e10
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
const Msun = 1.98847e33; const YR = 3.15576e7
const TAG = get(ENV, "S7B_TAG", "run"); const DIR = get(ENV, "S7B_DIR", joinpath(Main.TF.ROOT, "lib", "popiii", "s7b_x_$(TAG)")); mkpath(DIR)
const MSYNC = envf("S7B_MSYNC", 10.0) * Msun; const NSLAVE = envf("S7B_NSLAVE", 15.0); const SLEW = envf("S7B_SLEW", 0.05); const CLIP = envf("S7B_CLIP", 0.2)
const SAVEPREFIX = get(ENV, "S7B_CKPREFIX", "s7b_$(TAG)")
const PERFD = Dict{Symbol,Tuple{Int,Float64}}()
function perf(name, f)
    t0 = time_ns(); v = f(); cs = get(PERFD, name, (0, 0.0)); PERFD[name] = (cs[1] + 1, cs[2] + (time_ns() - t0) * 1e-9); v
end
mutable struct SS
    phase::Symbol; j::Int; t0::Float64; Mstar0::Float64; tb::Float64; dtau::Float64; Mdot::Float64; dtfree::Float64; texch::Float64; Macc_env::Float64
end
const ST = SS(:pre, 0, 0.0, 0.0, Inf, 0.0, 0.0, 0.0, 0.0, 0.0)
fname(kind, j) = joinpath(DIR, @sprintf("%s_%05d.txt", kind, j))

"Record of the star at the end of the step (d: the NamedTuple of s6b_run.jl `diag`), cgs."
record(d) = Dict("t" => d.t, "M" => d.M, "R" => d.R, "Lstar" => d.Lstar, "Lacc" => d.Lacc, "Linf" => d.Linf, "Tph" => d.Teff, "rph" => d.rph, "Ts" => d.Ts, "tau" => d.tau, "Tc" => d.Tc,
                 "rhoc" => d.ρc, "fD" => d.fD, "tKH" => d.tKH, "Mdot" => d.Mdot * Msun / YR, "Sion_ph" => d.Sion)

function exchange!(d)
    t_w0 = time(); j = ST.j
    perf(:exchange, () -> begin
        S.write_rec(fname("star", j), record(d))
        env = S.wait_rec(fname("env", j))
        ST.Macc_env = env["Macc"]
        Macc_star = d.M - ST.Mstar0
        Mdp = env["Mdot"]
        corr = (env["Macc"] - Macc_star) / (NSLAVE * env["dtau"])
        tgt = clamp(Mdp + corr, (1 - CLIP) * Mdp, (1 + CLIP) * Mdp)
        ST.Mdot = ST.j == 0 ? tgt : clamp(tgt, (1 - SLEW) * ST.Mdot, (1 + SLEW) * ST.Mdot)       # slew limit: the precursor boundary condition stalls on abrupt changes of Mdot
        ST.dtau = env["dtau"]; ST.tb = ST.t0 + env["tc"] + env["dtau"]
        env["save"] > 0.5 && Main.checkpoint(SAVEPREFIX * @sprintf("_%05d", j))
        if env["done"] > 0.5
            ST.phase = :done; Core.eval(Main, :(mend = 0.0))
        end
        @printf("S7B-star j=%d t_c=%.4f yr M=%.4f (env M_acc %.4f, star %.4f) Mdot_env=%.3e -> Mdot_star=%.3e Msun/yr, next barrier in %.4f yr\n", j, (d.t - ST.t0) / YR, d.M / Msun, env["Macc"] / Msun + ST.Mstar0 / Msun,
                d.M / Msun, Mdp * YR / Msun, ST.Mdot * YR / Msun, ST.dtau / YR)
    end)
    ST.j += 1; ST.texch += time() - t_w0
    nothing
end

"Called after every accepted step with the diagnostics d."
function on_step(d)
    if ST.phase === :pre
        if d.M >= MSYNC
            ST.phase = :coupled; ST.t0 = d.t; ST.Mstar0 = d.M; ST.Mdot = d.Mdot * Msun / YR
            @printf("S7B-star: coupling starts at t = %.6e s, M = %.5f Msun, R = %.1f Rsun, Mdot = %.3e Msun/yr\n", d.t, d.M / Msun, d.R / SR, d.Mdot)
            get(ENV, "S7B_SAVESYNC", "0") == "1" && Main.checkpoint(SAVEPREFIX * "_sync")
            exchange!(d)
        end
    elseif ST.phase === :coupled
        if d.t >= ST.tb * (1 - 1e-9)
            exchange!(d)
            ST.phase === :coupled && Main.TF.setscalar!(:dtime, max(Main.TF.getscalar(:dtime), 0.7 * ST.dtfree))
        end
    end
    nothing
end
"Cap of the next time step: land on the barrier (a step that would end within 2 dt of it is split in two)."
function dt_cap(dtn)
    ST.phase === :coupled || return dtn
    rem = ST.tb - Main.TF.getscalar(:timen)
    if rem <= dtn; return max(rem, 1e-3 * dtn)
    elseif rem < 2dtn; return 0.5 * rem
    else ST.dtfree = dtn; return dtn; end
end
mdot_now() = ST.phase === :coupled || ST.phase === :done ? ST.Mdot : nothing
function finish()
    mkpath(joinpath(Main.TF.ROOT, "lib", "popiii", "perf"))
    open(joinpath(Main.TF.ROOT, "lib", "popiii", "perf", "s7b_$(TAG)_star_sections.csv"), "w") do io
        println(io, "section,calls,total_s,mean_ms")
        for (kk, (c, s)) in sort(collect(PERFD); by = x -> -x[2][2]); @printf(io, "%s,%d,%.4f,%.4f\n", kk, c, s, 1e3 * s / max(c, 1)); end
        @printf(io, "wall_total,1,%.4f,%.4f\n", time() - Main.t0wall, 1e3 * (time() - Main.t0wall))
    end
end
Main.S7BH[] = on_step; Main.S7BDT[] = dt_cap; Main.S7BMDOT[] = mdot_now
end # module
