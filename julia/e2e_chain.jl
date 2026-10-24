# E2E-A: the chained ABN02 z = 19 run, from the initial cloud to the end of the main sequence of the star that forms (here without radiative feedback).
#   julia --project=env e2e_chain.jl [stage ...] [--force] [--status]        stages (default: all, in order, each resumes from the saved output of the previous one):
#     collapse    popiii_m2c.jl lines 200   cloud -> n_c = 1e21 (final2 settings: IC=abn02z19 TCMB=54.5, TR-BDF2, AVZONE=2, NLTE=hyb, CHEM=netatomicli); final state saved
#     formation   popiii_m2c.jl lines 200   RESTART collapse snapshot, past protostar formation (3000 steps, n_c ~ 7e22; LINEEVERY=5); final state saved
#     sink        popiii_m2c.jl lines 200   RESTART formation snapshot, SINK=1 (S6b Part 1 method) until M_pt >= E2E_MCAP or the mass flux < E2E_MDOTSTOP; then s6c_mdot.jl -> Mdot(M*)
#     protostar   s6b_run.jl                accreting protostar fed by Mdot(M*) (S6b Part 3) to M_pt of the sink run; HANDOFF_OUT = the star (grid, profiles, Y_i(m))
#     ms          s5_run.jl                 HANDOFF = that star, no accretion, to central hydrogen exhaustion (X_c < 1e-3)
#   Every stage writes its log (lib/popiii/*.log), its performance record lib/popiii/perf/<E2E_PERF>_<stage>_* and a marker <E2E_WORK>/<stage>.done; a stage with a marker is
#   skipped (--force reruns it and everything after it). protostar and ms restart from their last checkpoint when they stop early (up to E2E_RETRIES times).
#   FULL-1 (E2E_FB=1): the chain WITH feedback, stages  collapse -> formation -> coupled -> ms:
#     coupled     popiii_m2c.jl (S7B=1, SINK=1, envelope) + s6b_run.jl (S7B_STAR=1, protostar) as two processes exchanging at barriers (S7b two-way scheme, S7c-1 precursor fixes PROTGAMMA 0.25,
#                 S8a RT_LYMAN=1 RT_LYA=mc and the Mdot starvation estimator, HI21 inline, EDDCAP 0.95); ends when the envelope reports done (Mdot_win < E2E_MDOTSTOP [1e-3] or t_c > E2E_TCMAX [40000 yr]);
#                 HANDOFF_OUT = the star at that time.  A dead star process is restarted from its last checkpoint with S7B_RESUME (replay of the recorded exchange) while the envelope waits at the barrier;
#                 a dead envelope fails the stage (rerun it with --force: the envelope has no resume).
#     ms          as below, from the hand-off of `coupled`.
#   E2E_REUSE=<work dir of an earlier chain> imports its collapse / formation snapshots (the stages are then marked done with `reused=<dir>`).
#   S7b (feedback, no-FB chain): set E2E_TAG / E2E_PERF to new names and put the feedback switches of the drivers into E2E_EXTRA_<STAGE> ("KEY=VAL KEY=VAL"), e.g. E2E_EXTRA_PROTOSTAR.
# Environment: E2E_TAG [abn02z19_nofb], E2E_WORK [lib/popiii/e2e_work/<tag>] (snapshots, checkpoints, hand-off: not committed), E2E_THREADS [4], E2E_MCAP [250] Msun,
#   E2E_MDOTSTOP [1e-5] Msun/yr, E2E_PERF [e2e], E2E_RETRIES [4].
using Printf
const J = @__DIR__; const ROOT = normpath(J, ".."); const LIB = joinpath(ROOT, "lib", "popiii")
const FB = get(ENV, "E2E_FB", "0") == "1"
const TAG = get(ENV, "E2E_TAG", FB ? "abn02z19_fb" : "abn02z19_nofb")
const WORK = get(ENV, "E2E_WORK", joinpath(LIB, "e2e_work", TAG)); mkpath(WORK); mkpath(joinpath(WORK, "ck"))
const THREADS = get(ENV, "E2E_THREADS", "4")
const MCAP = parse(Float64, get(ENV, "E2E_MCAP", "250")); const MDOTSTOP = get(ENV, "E2E_MDOTSTOP", "1e-5")
const PERF = get(ENV, "E2E_PERF", "e2e"); const RETRIES = parse(Int, get(ENV, "E2E_RETRIES", "4"))
const STAGES = FB ? ["collapse", "formation", "coupled", "ms"] : ["collapse", "formation", "sink", "protostar", "ms"]
marker(st) = joinpath(WORK, st * ".done")
done(st) = isfile(marker(st))
readdone(st) = Dict(split(l, "=", limit = 2) |> p -> String(p[1]) => String(p[2]) for l in readlines(marker(st)) if occursin('=', l))
function writedone(st, d::Pair...)
    open(marker(st), "w") do io
        println(io, "stage=", st); println(io, "finished=", Base.Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
        for (k, v) in d; println(io, k, "=", v); end
    end
end
extra(st) = [String(p[1]) => String(p[2]) for p in (split(w, "=", limit = 2) for w in split(get(ENV, "E2E_EXTRA_" * uppercase(st), ""))) if length(p) == 2]
"The command of `julia -t THREADS --project=env script args` in J with the environment `env` (+ E2E_EXTRA_<STAGE>), log WORK/<stage>[suffix].out."
function mkcmd(stage, script, args, env; suffix = "", threads = THREADS)
    full = Dict{String,String}(env); for (k, v) in extra(stage); full[k] = v; end
    full["TITAN_NOLIB"] = "1"
    out = joinpath(WORK, stage * suffix * ".out")
    # the child is started as `e2ejl ... -e include(ENV["E2E_SCRIPT"])`: a `pkill julia` / `pkill -f <driver>.jl` aimed at other jobs on this shared machine does not match it
    full["E2E_SCRIPT"] = joinpath(J, script)
    cmd = setenv(`bash -c 'exec -a e2ejl julia -t "$0" --project=env -e "include(ENV[\"E2E_SCRIPT\"])" "$@"' $threads $args`, merge(Dict(ENV), full); dir = J)
    @printf("[e2e] %s: %s %s  (log %s)\n", stage, script, join(args, " "), out)
    pipeline(ignorestatus(cmd); stdout = out, stderr = out)
end
function launch(stage, script, args, env; suffix = "", threads = THREADS)
    t0 = time()
    p = run(mkcmd(stage, script, args, env; suffix, threads))
    p.termsignal != 0 && (@printf("[e2e] %s: killed by signal %d after %.0f s\n", stage, p.termsignal, time() - t0); return -p.termsignal)
    @printf("[e2e] %s: exit %d after %.0f s\n", stage, p.exitcode, time() - t0)
    p.exitcode
end
logtail(file, n = 40) = isfile(file) ? join(last(readlines(file), n), "\n") : ""
# ----------------------------------------------------------------------------------------------------------------- stages 1-3: popiii_m2c.jl
"Environment of the ABN02 z = 19 collapse driver (final2 settings of the notes, S6c chemistry)."
collapse_env() = Dict("IC" => "abn02z19", "TCMB" => "54.5", "AVZONE" => "2", "INTEGRATOR" => "trbdf2", "CHEM" => "netatomicli", "NLTE" => "hyb", "SNAPS" => "1", "SNAPDIR" => WORK,
                                 "XFLOOR" => "1e-12", "CHEMFALLBACK" => "1", "THICKCONT" => "3", "TEQMAX" => "1e30")   # (atomic net: He+, He++, Li+ start at 0; ChemNet's Newton needs XFLOOR, see notes E2E-A)
"Snapshot file of a popiii_m2c.jl run (its name carries the zone number after a sink removed zones: n200 -> n<N>)."
function snapfile(tag)
    c = filter(f -> startswith(f, "snaps_abn02z19_lines_n") && endswith(f, tag * ".jls") && occursin(r"^snaps_abn02z19_lines_n\d+" * tag * r"\.jls$", f), readdir(WORK))
    isempty(c) ? joinpath(WORK, "snaps_abn02z19_lines_n200$(tag).jls") : joinpath(WORK, first(c))
end
"popiii_m2c.jl catches its errors and exits 0: success means no `stopped:` line in its log."
m2c_ok(tag) = (l = joinpath(LIB, "m2c_abn02z19_lines_n200$(tag).log"); isfile(l) && (L = readlines(l); any(startswith(x, "wall time:") for x in L) && !any(startswith(x, "stopped:") for x in L)))
"""
Run popiii_m2c.jl for `tag0`; a run that was killed or failed is continued from the last snapshot of the previous attempt (`<tag0>_r1`, `_r2`, ...; the collapse drivers write one every 0.25 dex of n_c),
a finished attempt is reused. Returns the tag of the finishing attempt.
"""
function run_m2c(stage, tag0, env; resumable = true, nsh = "200", resume_env = Dict{String,String}())
    tags = [tag0; ["$(tag0)_r$i" for i in 1:RETRIES]]; prev = nothing
    for t in tags
        m2c_ok(t) && return t
        if isfile(joinpath(LIB, "m2c_abn02z19_lines_n200$(t).log"))      # an earlier, unfinished attempt: continue from its last snapshot
            isfile(snapfile(t)) && (prev = snapfile(t)); continue
        end
        e = copy(env); e["TAG"] = t
        resumable && prev !== nothing && (e["RESTART"] = prev; e["RESIDX"] = "0"; haskey(env, "PERF") && (e["PERF"] = env["PERF"] * "_" * last(split(t, "_")));  e["TMAXYR"] = get(env, "TMAXYR", "1e8"); merge!(e, resume_env))
        rc = launch(stage, "popiii_m2c.jl", ["lines", nsh], e; suffix = t == tag0 ? "" : last(split(t, "_")) )
        m2c_ok(t) && rc == 0 && return t
        isfile(snapfile(t)) && (prev = snapfile(t))
        @printf("[e2e] %s: attempt %s did not finish (rc %d)\n", stage, t, rc)
    end
    error("$stage failed after $(RETRIES + 1) attempts, see $(WORK)/$stage*.out")
end
"E2E_REUSE=<work dir>: take the stage's snapshot from an earlier chain run (its <stage>.done), copy it here and mark the stage done."
function reuse_stage(st)
    src = get(ENV, "E2E_REUSE", ""); isempty(src) && return false
    m = joinpath(src, st * ".done"); isfile(m) || return false
    d = Dict(split(l, "=", limit = 2) |> p -> String(p[1]) => String(p[2]) for l in readlines(m) if occursin('=', l))
    snap = d["snapshot"]; isfile(snap) || (snap = joinpath(src, basename(snap)))
    isfile(snap) || return false
    dst = joinpath(WORK, basename(snap)); isfile(dst) || cp(snap, dst)
    writedone(st, "snapshot" => dst, "log" => get(d, "log", ""), "reused" => src); @printf("[e2e] %s: reused from %s\n", st, src); true
end
function stage_collapse()
    reuse_stage("collapse") && return
    tag = "_e2e_$(TAG)_c1"
    env = merge(collapse_env(), Dict("TAG" => tag, "NSTOP" => "1e21", "SNAPEND" => "1", "PERF" => "$(PERF)_collapse", "PERFSTAGE" => "collapse", "MAXSTEPS" => "60000", "TSTOP" => "1e30"))
    tag = run_m2c("collapse", tag, env)
    writedone("collapse", "snapshot" => snapfile(tag), "log" => joinpath(LIB, "m2c_abn02z19_lines_n200$(tag).log"))
end
function stage_formation()
    reuse_stage("formation") && return
    c = readdone("collapse"); tag = "_e2e_$(TAG)_c2"
    env = merge(collapse_env(), Dict("TAG" => tag, "RESTART" => c["snapshot"], "RESIDX" => "0", "LINEEVERY" => "5", "NSTOP" => "1e30", "TSTOP" => "1e30", "TMAXYR" => "10",
                                     "MAXSTEPS" => "3000", "SNAPEND" => "1", "PERF" => "$(PERF)_formation", "PERFSTAGE" => "formation", "SNAPS" => "1"))
    tag = run_m2c("formation", tag, env)
    writedone("formation", "snapshot" => snapfile(tag), "log" => joinpath(LIB, "m2c_abn02z19_lines_n200$(tag).log"))
end
"""
Smooth the windowed Mdot(M*) rows of s6c_mdot.jl (`s6c_<tag>_mdot.csv`: M, Mdot) with a Gaussian in ln M (sigma = SIGMA, weights on ln Mdot) at nodes every 0.15 in ln M.
The windowed rates zig-zag by +-20 % from node to node (the shells cross the sink one zone at a time); the protostar's luminosity-wave onset at 0.6 M* depends on those wiggles (S6b stall,
notes E2E-A), the physical rate does not.
"""
function smooth_mdot(infile, outfile; sigma = 0.4)
    rows = [parse.(Float64, split(l, ",")[1:2]) for l in readlines(infile)[2:end] if !isempty(strip(l))]
    lm = log.([r[1] for r in rows]); ly = log.([r[2] for r in rows])
    nodes = collect(range(lm[1], lm[end]; step = 0.15))
    open(outfile, "w") do io
        println(io, "M_Msun,Mdot_Msun_yr")
        for x in nodes
            w = exp.(-0.5 .* ((lm .- x) ./ sigma) .^ 2); println(io, exp(x), ",", exp(sum(w .* ly) / sum(w)))
        end
    end
    outfile
end
function stage_sink()
    f = readdone("formation"); tag = "_e2e_$(TAG)_sink"
    env = merge(collapse_env(), Dict("TAG" => tag, "RESTART" => f["snapshot"], "RESIDX" => "0", "SINK" => "1", "SPLIT" => "2", "RSINKAU" => "5", "RSINKP" => "1", "RSMREF" => "1",
                                     "RSMFLOOR" => "0.01", "CDYN" => "0.03", "NLTE" => "ep", "XFLOOR" => "1e-12", "CHEMFALLBACK" => "1", "TMAXYR" => "1e5", "NSTOP" => "1e30", "TSTOP" => "1e30",
                                     "MAXSTEPS" => "400000", "SNAPS" => "1", "SNAPEVERY" => "1000", "MSTOPMSUN" => string(MCAP), "MDOTSTOP" => MDOTSTOP, "PERF" => "$(PERF)_sink", "PERFSTAGE" => "sink"))
    delete!(env, "SNAPEND")
    tag = run_m2c("sink", tag, env; resume_env = Dict("SPLIT" => "1"))      # a restart from a sink snapshot (already split, zones below r_s removed again by the driver)
    # accretion history -> Mdot(M*) tables (lib/popiii/s6c_<tag>_mdot{,_smooth}.csv)
    run(pipeline(setenv(`julia --project=env $(joinpath(J, "s6c_mdot.jl")) $tag $(f["snapshot"])`, merge(Dict(ENV), Dict("TITAN_NOLIB" => "1")); dir = J); stdout = joinpath(LIB, "s6c_$(tag)_mdot.txt")))
    L = readlines(joinpath(LIB, "s6c_$(tag)_mdot_smooth.csv"))[2:end]; Mend = parse(Float64, split(last(L), ",")[1])
    S = readlines(joinpath(LIB, "s6c_$(tag).csv")); Mpt = parse(Float64, split(last(S), ",")[6]); tyr = parse(Float64, split(last(S), ",")[2])
    writedone("sink", "tag" => tag, "mdot_table" => smooth_mdot(joinpath(LIB, "s6c_$(tag)_mdot.csv"), joinpath(LIB, "s6c_$(tag)_mdot_smooth2.csv")), "Mpt_final" => string(Mpt), "t_yr_final" => string(tyr), "Mtable_end" => string(Mend))
end
# ----------------------------------------------------------------------------------------------------------------- stage 4: protostar fed by Mdot(M*) (S6b Part 3)
"Newest checkpoint `<prefix>_auto{0,1}` in WORK/ck (by modification time of the .aux file), or nothing."
function latest_ck(drv, prefix)
    c = [(mtime(joinpath(WORK, "ck", "$(drv)_$(prefix)_auto$(i).aux")), "$(prefix)_auto$(i)") for i in 0:1 if isfile(joinpath(WORK, "ck", "$(drv)_$(prefix)_auto$(i).aux"))]
    isempty(c) ? nothing : last(sort(c))[2]
end
finished(log, key) = isfile(log) && any(l -> occursin(key, l), readlines(log))
function stage_protostar()
    s = readdone("sink"); Mend = parse(Float64, s["Mpt_final"]); out0 = "e2e_$(TAG)_protostar"
    base = Dict("OUT" => out0, "PRECTHICK" => "1", "BCSMOOTH" => "0.03", "BCSLOPE" => "8", "BCFALLBACK" => "1", "ESLIM" => "1.05", "TOL" => "0.05", "NIN" => "80", "QC" => "1.083",
                "QSPLIT" => "0.3", "NENV" => "70", "PSI" => "0", "MDOTFILE" => s["mdot_table"], "MDOTRAMP" => "4.41e-3", "M0GUESS" => "0.22", "MDOTRAMPP" => "2.5",
                "MEND" => string(Mend), "MAXSTEPS" => "400000", "EVERY" => "50", "CSVEVERY" => "200", "CKEVERY" => "300", "CKDIR" => joinpath(WORK, "ck"), "SAVE" => "e2e_$(TAG)_p3",
                "PROFM" => "1,3,10,30,60,100,150,200", "EDDCAP" => "0.95", "PERFSTAGE" => "protostar", "HANDOFF_OUT" => joinpath(WORK, "handoff_$(TAG).star"))
    ck0 = latest_ck("s6b", "e2e_$(TAG)_p3")
    attempt = ck0 === nothing ? 0 : 1; restart = ck0
    while true
        env = copy(base); env["PERF"] = attempt == 0 ? "$(PERF)_protostar" : "$(PERF)_protostar_r$(attempt)"
        attempt > 0 && (env["OUT"] = "$(out0)_r$(attempt)"; env["RESTART"] = restart)
        launch("protostar", "s6b_run.jl", String[], env; suffix = attempt == 0 ? "" : "_r$(attempt)")
        log = joinpath(LIB, "s6b_$(env["OUT"]).log")
        if isfile(base["HANDOFF_OUT"]) && finished(log, "mass accretion check") && !finished(log, "STOPPED at step") && mtime(base["HANDOFF_OUT"]) >= mtime(log) - 600
            writedone("protostar", "handoff" => base["HANDOFF_OUT"], "log" => log, "attempts" => string(attempt + 1), "Mend" => string(Mend)); return
        end
        attempt += 1; ck = latest_ck("s6b", "e2e_$(TAG)_p3")
        (attempt > RETRIES || ck === nothing) && error("protostar stage stopped early (see $(log)): $(logtail(log, 8))")
        restart = ck; @printf("[e2e] protostar stopped early; restart %d from checkpoint %s\n", attempt, ck)
        rm(base["HANDOFF_OUT"]; force = true)
    end
end
# ----------------------------------------------------------------------------------------------------------------- FULL-1 stage 3+4: coupled envelope + protostar with feedback
"Two-way coupled sink run (envelope, popiii_m2c.jl S7B=1) and protostar (s6b_run.jl S7B_STAR=1) of the S7b scheme, two processes of E2E_THREADS threads; see the header."
function stage_coupled()
    f = readdone("formation"); etag = "_e2e_$(TAG)_cpl"; X = joinpath(WORK, "x_cpl"); mkpath(X)
    tcmax = get(ENV, "E2E_TCMAX", "40000"); mdotstop = get(ENV, "E2E_MDOTSTOP_FB", "1e-3"); lya = get(ENV, "E2E_LYA", "mc")
    pre = get(ENV, "E2E_MDOTPRE", joinpath(LIB, "s6c__z19a2_mdot_smooth.csv"))      # Mdot(M*) of the collapse for the star until the coupling starts (M = 3 Msun); the S7b/S7c/S8a validations use this table
    shared = Dict("S7B_TAG" => TAG, "S7B_DIR" => X, "S7B_MSYNC" => "3", "S7B_SLEWDN" => "0.2", "S7B_DTMAX" => "60", "S7B_PERFPREFIX" => PERF, "RT_LYMAN" => "1", "RT_LYA" => lya, "S7B_STARVECUT" => get(ENV, "E2E_STARVECUT", "3"))
    envp = merge(collapse_env(), shared, Dict("TAG" => etag, "RESTART" => f["snapshot"], "RESIDX" => "0", "SINK" => "1", "SPLIT" => "2", "RSINKAU" => "5", "RSINKP" => "1", "RSMREF" => "1", "RSMFLOOR" => "0.01",
        "CDYN" => "0.03", "NLTE" => "ep", "TMAXYR" => "1e5", "NSTOP" => "1e30", "TSTOP" => "1e30", "MAXSTEPS" => "600000", "SNAPS" => "0", "S7B" => "1", "S7B_TCMAX" => tcmax, "S7B_MDOTSTOP" => mdotstop,
        "S7B_NPROF" => "12", "S7B_SAVEAT" => get(ENV, "E2E_SAVEAT", "10,30,50,75,100,120,130,140"), "HI21" => "1", "HI21EVERY" => get(ENV, "E2E_HI21EVERY", "200"), "HI21RUN" => "e2e_$(TAG)_coupled",
        "PERF" => "$(PERF)_coupled_env", "PERFSTAGE" => "coupled_env"))
    elog = joinpath(LIB, "m2c_abn02z19_lines_n200$(etag).log"); efinal = joinpath(X, "env_final.txt")
    out0 = "e2e_$(TAG)_star"; ho = joinpath(WORK, "handoff_$(TAG).star"); slog(o) = joinpath(LIB, "s6b_$(o).log")
    starbase = Dict("OUT" => out0, "PRECTHICK" => "1", "BCSMOOTH" => "0.03", "BCSLOPE" => "8", "BCFALLBACK" => "1", "ESLIM" => "1.05", "TOL" => "0.05", "NIN" => "80", "QC" => "1.083", "QSPLIT" => "0.3",
        "NENV" => "70", "PSI" => "0", "MDOTFILE" => pre, "MDOTRAMP" => "4.41e-3", "M0GUESS" => "0.22", "MDOTRAMPP" => "2.5", "MAXSTEPS" => "600000", "EVERY" => "50", "CSVEVERY" => "200", "CKEVERY" => "100",
        "CKDIR" => joinpath(WORK, "ck"), "SAVE" => "e2e_$(TAG)_cp", "S7B_STAR" => "1", "PROTGAMMA" => "0.25", "EDDCAP" => "0.95", "PROFM" => "1,3,10,30,60,100,150,200", "HANDOFF_OUT" => ho, "PERFSTAGE" => "coupled_star")
    merge!(starbase, shared)
    star_ok(log) = isfile(ho) && finished(log, "mass accretion check") && !finished(log, "STOPPED at step") && mtime(ho) >= mtime(log) - 600
    env_ok() = m2c_ok(etag) && isfile(efinal)
    t0 = time(); pe = nothing; ps = nothing; attempt = 0; slogf = slog(out0)
    if !env_ok()
        isfile(elog) && error("coupled: an unfinished envelope run exists ($elog): the envelope has no resume; rerun the stage with --force")
        rm(X; recursive = true, force = true); mkpath(X)
        pe = run(mkcmd("coupled_env", "popiii_m2c.jl", ["lines", "200"], envp; threads = THREADS); wait = false); @printf("[e2e] coupled: envelope pid %d\n", getpid(pe)); sleep(20)
    end
    if !star_ok(slogf)
        ck0 = latest_ck("s6b", "e2e_$(TAG)_cp")
        env = copy(starbase); env["PERF"] = "$(PERF)_coupled_star"
        if ck0 !== nothing && env_ok(); attempt = 1; env["RESTART"] = ck0; env["S7B_RESUME"] = "1"; env["S7B_ENVDIR"] = X; env["OUT"] = "$(out0)_r1"; env["PERF"] *= "_r1"; slogf = slog(env["OUT"]); end
        ps = run(mkcmd("coupled_star", "s6b_run.jl", String[], env; suffix = attempt == 0 ? "" : "_r$(attempt)"); wait = false); @printf("[e2e] coupled: star pid %d\n", getpid(ps))
    end
    while true
        sleep(30)
        e_run = pe !== nothing && process_running(pe); s_run = ps !== nothing && process_running(ps)
        if pe !== nothing && !e_run && !env_ok()
            ps !== nothing && s_run && kill(ps)
            error("coupled: the envelope process ended without finishing (see $(WORK)/coupled_env.out, $elog)")
        end
        if ps !== nothing && !s_run
            if star_ok(slogf)
                ps = nothing
            else
                attempt += 1; ck = latest_ck("s6b", "e2e_$(TAG)_cp")
                (attempt > RETRIES || ck === nothing) && (pe !== nothing && e_run && kill(pe); error("coupled: the star stopped early (see $(slogf)): $(logtail(slogf, 8))"))
                env = copy(starbase); env["RESTART"] = ck; env["S7B_RESUME"] = "1"; env["S7B_ENVDIR"] = X; env["OUT"] = "$(out0)_r$(attempt)"; env["PERF"] = "$(PERF)_coupled_star_r$(attempt)"
                slogf = slog(env["OUT"]); rm(ho; force = true)
                @printf("[e2e] coupled: star stopped early; restart %d from checkpoint %s (envelope waits at its barrier)\n", attempt, ck)
                ps = run(mkcmd("coupled_star", "s6b_run.jl", String[], env; suffix = "_r$(attempt)"); wait = false)
            end
        end
        (pe === nothing || !e_run) && ps === nothing && break
    end
    env_ok() || error("coupled: envelope did not finish cleanly ($elog)")
    star_ok(slogf) || error("coupled: no hand-off star (see $slogf)")
    writedone("coupled", "handoff" => ho, "log" => slogf, "env_log" => elog, "env_csv" => joinpath(LIB, "s7b_$(TAG)_env.csv"), "attempts" => string(attempt + 1), "wall_s" => string(round(Int, time() - t0)))
end
# ----------------------------------------------------------------------------------------------------------------- stage 5: main sequence
function stage_ms()
    p = readdone(FB ? "coupled" : "protostar"); out0 = "e2e_$(TAG)_ms"
    base = Dict("HANDOFF" => p["handoff"], "OUT" => out0, "INT" => "theta", "DXMAX" => "4e-3", "EFRAC" => "5", "CONV" => "1e-5", "XSTOP" => "1e-3", "MAXSTEPS" => "100000", "EVERY" => "50",
                "CKEVERY" => "100", "CKDIR" => joinpath(WORK, "ck"), "SAVE" => "e2e_$(TAG)_ms", "PROFILE_X" => "0.7,0.5,0.3,0.1,0.01", "EDDCAP" => "0.95", "PERFSTAGE" => "ms", "TOL" => "0.01", "NITER" => "40")      # (TOL 0.03 / NITER 25: Newton fails at t = 1.9e4 yr, R = 7.3 Rsun, see notes E2E-A)
    ck0 = latest_ck("s5", "e2e_$(TAG)_ms")      # a run that was killed leaves its checkpoints: continue from the newest one
    attempt = ck0 === nothing ? 0 : 1; restart = ck0
    while true
        env = copy(base); env["PERF"] = attempt == 0 ? "$(PERF)_ms" : "$(PERF)_ms_r$(attempt)"
        attempt > 0 && (env["OUT"] = "$(out0)_r$(attempt)"; env["RESTART"] = restart)      # (HANDOFF stays: it builds the grid, the checkpoint then supplies the state)
        launch("ms", "s5_run.jl", String[], env; suffix = attempt == 0 ? "" : "_r$(attempt)")
        log = joinpath(LIB, "s5_$(env["OUT"]).log")
        if finished(log, "central hydrogen exhausted") || finished(log, "TEND reached")
            writedone("ms", "log" => log, "attempts" => string(attempt + 1)); return
        end
        attempt += 1; ck = latest_ck("s5", "e2e_$(TAG)_ms")
        if attempt > RETRIES || ck === nothing
            writedone("ms_stopped", "log" => log); error("ms stage stopped early (see $(log)): $(logtail(log, 8))")
        end
        restart = ck; @printf("[e2e] ms stopped early; restart %d from checkpoint %s\n", attempt, ck)
    end
end
const FN = Dict("coupled" => stage_coupled, "collapse" => stage_collapse, "formation" => stage_formation, "sink" => stage_sink, "protostar" => stage_protostar, "ms" => stage_ms)
function main()
    args = filter(a -> !startswith(a, "--"), ARGS); force = "--force" in ARGS
    if "--status" in ARGS
        for st in STAGES; @printf("%-10s %s\n", st, done(st) ? "done  " * join(["$k=$v" for (k, v) in readdone(st) if k in ("snapshot", "handoff", "Mpt_final", "mdot_table")], " ") : "-"); end
        return
    end
    todo = isempty(args) ? STAGES : args
    for st in todo
        st in STAGES || error("unknown stage $st")
        if force
            for s2 in STAGES[findfirst(==(st), STAGES):end]; rm(marker(s2); force = true); end
        end
        if done(st); println("[e2e] $st: done (", marker(st), "), skipped"); continue; end
        i = findfirst(==(st), STAGES)
        i > 1 && !done(STAGES[i-1]) && error("stage $st needs the output of $(STAGES[i-1])")
        FN[st]()
    end
end
(abspath(PROGRAM_FILE) == @__FILE__()) && main()
