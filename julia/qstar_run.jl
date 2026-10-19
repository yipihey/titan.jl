# Phase 7d/e: quasi-star outer shell in pure-Julia TITAN.
#   TITAN_NOLIB=1 julia --project=env qstar_run.jl setup PROFILE.csv [rin nenv natm]   → reference_runs/init/QSTAR.tstate.gz
#   TITAN_NOLIB=1 julia --project=env qstar_run.jl relax [tmax_s maxsteps]           run (hydrostatic test, or infall if set up)
#   relax options (env): QS_NAME, QS_INIT (set-up to start from), QS_FROM (saved state to continue), QS_SET (scalars),
#     QS_GASROW / QS_COMPTON (=1: gas-only energy row / Compton exchange, see TitanResidFast), QS_RAMP_FAC/QS_RAMP_T, QS_EVERY,
#     QS_SINK_TAU (inner mass sink relaxation time [s]; with QS_FROM give QS_SINK_M0 [g] = the original shell mass),
#     QS_TRBDF2=1 (second-order time stepping), QS_CFL (advective Courant cap on Δt), QS_DTMAX (cap on Δt [s])
# Run convention (koine docs/plan-run-convention.md): every QS_* variable is a declared parameter (`--schema` prints them with
# units, groups and mutability) and a run writes a conforming directory lib/qs_runs/NAME/ (run.json, params.json, events.jsonl,
# outputs.json, control/) next to the old flat NAME_relaxed.tstate.gz / NAME_relax_snaps.jls, which are still written.
#   julia qstar_run.jl --schema                              the parameter schema
#   julia qstar_run.jl relax --show-params                   dry run: every resolved parameter, value and source
#   julia qstar_run.jl relax 1e30 5 --sink.tau=3e6           name=value on the command line (CLI > QS_* env > params file > default)
#   QS_FROM=lib/qs_runs/NAME/checkpoint-00000050.tstate.gz   inside a run directory: resumed as a child run (lineage recorded)
using Printf, Serialization
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl"))
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "QuasiStarSetup.jl")); include(joinpath(@__DIR__, "Trbdf2.jl"))
using .TitanDriver
using RunConvention
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const BS = BlockSolver
const TC = TitanCore; const TM = TitanMaterials; const QS = QuasiStarSetup; const Z0 = Z0Physics

initfile(name) = joinpath(TF.ROOT, "reference_runs", "init", "$name.tstate.gz")
outdir = joinpath(TF.ROOT, "lib", "qs_runs"); mkpath(outdir)

function install_julia!()
    TD.ROUTINES[:matgen] = () -> TX.matgen!()
    TD.ROUTINES[:matslv] = BS.Penta.matslv!
    for (s, f) in ((:peq, TC.peq!), (:reset, TC.reset!), (:update, TC.update!), (:maxdel, TC.maxdel!),
                   (:timstp, TC.timstp!), (:totnrg, TC.totnrg!))
        TD.ROUTINES[s] = f
    end
    TC.EDDFAC[] = TC.eddfac!
end

"Snapshot of the structure: interface r, u, L_rad; zone ρ, T, p_g; convective luminosity from the residual workspace."
function snap()
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    v(n) = copy(TF.var(n)::Vector{Float64})
    rn, un, frn, dn, tn, pgn, ern = v(:rn), v(:un), v(:frn), v(:dn), v(:tn), v(:pgn), v(:ern)
    fcv = isempty(TX.WSCACHE) ? zeros(length(rn)) : TX.ForwardDiff.value.(first(values(TX.WSCACHE)).a.fcv)
    (; t = TF.getscalar(:timen), dt = TF.getscalar(:dtime), a, b,
       r = rn[a:b+1], u = un[a:b+1], Lrad = [4π * rn[k]^2 * frn[k] for k in a:b+1],
       Lconv = [4π * rn[k]^2 * fcv[k] for k in a:b+1],
       d = dn[a:b], T = tn[a:b], pg = pgn[a:b], er = ern[a:b])
end


"Quantities `report` prints (and the run records as metrics)."
function qstats(s, L)
    cs = [sqrt(4 / 3 * (s.pg[i] + s.er[i] / 3) / s.d[i]) for i in eachindex(s.d)]
    mach = maximum(abs(s.u[i]) / cs[min(i, end)] for i in eachindex(s.u))
    Ltot = s.Lrad .+ s.Lconv
    (; mach, lmin = minimum(Ltot[2:end-1]) / L, lmax = maximum(Ltot[2:end-1]) / L, lconv = maximum(s.Lconv) / L)
end

function report(s, L, label)
    q = qstats(s, L)
    @printf("%s t = %.4e s  dt = %.3e  max|u|/cs = %.3e  L/L0: min %.4f max %.4f (interior)  Lconv/L max %.3f\n",
            label, s.t, s.dt, q.mach, q.lmin, q.lmax, q.lconv)
end

# ---------------------------------------------------------------------------------------------------------------------
# Parameters: each QS_* variable is declared once (name, type, default, unit, group, mutability, env). `p.group.name` reads.
# All declarations sit here, ahead of the mode dispatch, so that `--schema` lists them without running anything.
# ---------------------------------------------------------------------------------------------------------------------
@param run.tmax::Float64 = 3e9 unit="s" group="io" mutable=:restart doc="relax: stop at this model time (positional arg 2)"
@param run.maxsteps::Int = 400 group="io" mutable=:restart doc="relax: stop after this many steps (positional arg 3)"
@param io.name::String = "QSTAR" group="io" mutable=:restart env="QS_NAME" doc="run name: set-up file reference_runs/init/NAME.tstate.gz and the run directory lib/qs_runs/NAME"
@param io.init::String = "" group="io" mutable=:restart env="QS_INIT" doc="relax: set-up to start from (default: NAME)"
@param io.from::String = "" group="io" mutable=:restart env="QS_FROM" doc="relax: saved state to continue (the ramp keeps its clock); inside a run directory it is a resume"
@param io.every::Int = 100 group="io" mutable=:always env="QS_EVERY" doc="print/record every this many steps (and the first 3); checkpoint every 10x"
@param io.msg::Bool = false group="io" mutable=:always env="QS_MSG" doc="TITAN step messages to stdout"
@param io.bt::Bool = false group="io" mutable=:never env="QS_BT" doc="print the backtrace when a step fails"
@param setup.rin::Float64 = 0.4 unit="R" group="setup" mutable=:never doc="setup: inner radius in units of the star's radius (positional arg 3)"
@param setup.nenv::Int = 140 group="setup" mutable=:never doc="setup: envelope zones (positional arg 4)"
@param setup.natm::Int = 60 group="setup" mutable=:never doc="setup: atmosphere zones (positional arg 5)"
@param setup.grey::Bool = false group="setup" mutable=:never env="QS_GREY" doc="grey (rather than Z0 table) gas opacity; used by setup and relax"
@param setup.mdot::Float64 = 0.0 unit="Msun/yr" group="setup" mutable=:never env="QS_MDOT" doc="setup: infall rate at r_out"
@param setup.rout::Float64 = 2.0 unit="R" group="setup" mutable=:never env="QS_ROUT" doc="setup: outer radius of the infall region"
@param setup.dr0::Float64 = 5e-5 unit="R" group="setup" mutable=:never env="QS_DR0" doc="setup: first zone width"
@param setup.lgrid::Int = 2 group="setup" mutable=:never env="QS_LGRID" doc="setup: grid layout selector (QuasiStarSetup.setup! lgrid)"
@param setup.wcontact::Float64 = 0.0 group="setup" mutable=:never env="QS_WCONTACT" doc="setup: contact-layer weight (QuasiStarSetup.setup! wcontact)"
@param setup.equi::Bool = false group="setup" mutable=:never env="QS_EQUI" doc="setup: equilibrium initial state (QuasiStarSetup.setup! equi)"
@param setup.lady::String = "2,2,2" group="setup" mutable=:never env="QS_LADY" doc="setup: comma-separated triple of integers (QuasiStarSetup.setup! lady)"
@param physics.gasrow::Bool = false group="physics" mutable=:never env="QS_GASROW" doc="gas-only energy row (TitanResidFast.GASROW)"
@param physics.compton::Bool = false group="physics" mutable=:never env="QS_COMPTON" doc="Compton gas-radiation exchange (TitanResidFast.COMPTON)"
@param physics.mlt_noinfall::Bool = false group="physics" mutable=:never env="QS_MLT_NOINFALL" doc="no mixing-length convection in the inflow (TitanResidFast.MLT_NOINFALL)"
@param ramp.fac::Float64 = 1.0 group="ramp" mutable=:restart env="QS_RAMP_FAC" doc="inflow ramp: target/initial Mdot (density and fluxes at r_out grow exponentially by this factor); setup: Mdot_max"
@param ramp.t::Float64 = 1e6 unit="s" group="ramp" mutable=:restart env="QS_RAMP_T" doc="inflow ramp: e-folding time scale (the ramp lasts t*ln(fac))"
@param sink.tau::Float64 = 0.0 unit="s" group="sink" mutable=:never env="QS_SINK_TAU" requires="physics.gasrow" doc="Inner mass sink relaxation time; 0 = closed wall; needs physics.gasrow. The shell mass below the shock is held at sink.m0"
@param sink.m0::Float64 = 0.0 unit="g" group="sink" mutable=:never env="QS_SINK_M0" doc="shell mass the sink holds; 0 = measured from the loaded state (the value used is recorded in params.json)"
@param integrator.trbdf2::Bool = false group="integrator" mutable=:never env="QS_TRBDF2" doc="second-order TR-BDF2 time stepping (Trbdf2.jl)"
@param integrator.cfl::Float64 = 0.0 group="integrator" mutable=:restart env="QS_CFL" doc="advective Courant cap on the next step, dt <= cfl*min(dr/|u|); 0 = off"
@param integrator.dtmax::Float64 = 0.0 unit="s" group="integrator" mutable=:restart env="QS_DTMAX" doc="upper limit on dt; 0 = no limit (the old default was Inf, which canonical JSON cannot hold)"
@metric mach unit="" doc="max |u|/c_s over the interfaces"
@metric l_over_l0_min unit="" doc="min total luminosity / L0 over interior interfaces"
@metric l_over_l0_max unit="" doc="max total luminosity / L0 over interior interfaces"
@metric lconv_over_l_max unit="" doc="max convective luminosity / L0"
@metric sink_mdot unit="Msun/yr" doc="inner mass sink rate (when the sink is active)"
@metric sink_point_mass unit="Msun" doc="mass gained by the central point mass since the start (when the sink is active)"

# QS_SET=k=v,... wrote arbitrary Fortran COMMON scalars; now the open group `common` (--common.ql1=1e-4), generated from the
# COMMON tables of titan.com. Default 0 means "not set": a scalar is applied to the state only when set from env/CLI or nonzero.
function common_scalars()
    out = Tuple{String,Any,Any,String}[]
    for (blk, vs) in TF.BLOCKS, v in vs
        (isempty(v.dims) && v.T in (Float64, Int32)) || continue
        push!(out, (String(v.name), v.T, zero(v.T), "Fortran COMMON /$blk/ scalar"))
    end
    sort!(out; by = first)
end
@open_group common common_scalars

const SETUP_ONLY = ["setup.rin", "setup.nenv", "setup.natm", "setup.mdot", "setup.rout", "setup.dr0", "setup.lgrid",
                    "setup.wcontact", "setup.equi", "setup.lady"]
const RELAX_ONLY = ["run.tmax", "run.maxsteps", "io.init", "io.from", "io.every", "io.msg", "io.bt", "physics.gasrow",
                    "physics.compton", "physics.mlt_noinfall", "ramp.t", "sink.tau", "sink.m0", "integrator.trbdf2",
                    "integrator.cfl", "integrator.dtmax"]

maybe_print_schema() && exit()

# Positional arguments (mode, then tmax maxsteps | profile rin nenv natm) become --name=value, ahead of ARGS so that an explicit
# --name=value still wins; QS_SET entries go in first of all (they are the weakest of the three, as env is below the CLI).
pos = filter(a -> !startswith(a, "--"), ARGS)
mode = isempty(pos) ? error("mode: setup | relax (or --schema)") : pos[1]
posargv = String[]
if mode == "relax"
    length(pos) >= 2 && push!(posargv, "--run.tmax=" * pos[2])
    length(pos) >= 3 && push!(posargv, "--run.maxsteps=" * pos[3])
elseif mode == "setup"
    length(pos) >= 3 && push!(posargv, "--setup.rin=" * pos[3])
    length(pos) >= 4 && push!(posargv, "--setup.nenv=" * pos[4])
    length(pos) >= 5 && push!(posargv, "--setup.natm=" * pos[5])
end
qsset = Pair{String,String}[]
for kv in split(get(ENV, "QS_SET", ""), ","; keepempty = false)      # e.g. QS_SET=ql1=1e-4,cq2=1
    k, val = split(kv, "="; limit = 2); push!(qsset, String(k) => String(val))
end
qssetargv = ["--common.$k=$v" for (k, v) in qsset]
base_argv = vcat(qssetargv, posargv, ARGS)
# where QS_SET (not the command line) supplied a COMMON scalar, say so in params.json
set_from_env(n) = any(kv -> "common." * kv.first == n, qsset) && !any(a -> startswith(a, "--$n="), ARGS)
if "--show-params" in ARGS                      # dry run: the resolved parameters and where each value came from
    p = params(; argv = base_argv)
    vals, srcs = RunConvention.values_of(p), RunConvention.sources_of(p)
    for n in keys(p); println(n, "\t", vals[n], "\t", srcs[n] == "cli" && set_from_env(n) ? "env:QS_SET" : srcs[n]); end
    exit()
end

"Relabel sources the library cannot know (QS_SET is one variable feeding many parameters; derived values) and rewrite params.json."
function relabel!(rc, derived = String[])
    P = rc.params
    for n in P.order
        P.sources[n] == "cli" && set_from_env(n) && (P.sources[n] = "env:QS_SET")
    end
    for n in derived; P.sources[n] = "derived:state"; end
    RunConvention.write_json_atomic(joinpath(rc.dir, "params.json"), RunConvention.params_document(P))
end

"First free run directory lib/qs_runs/NAME, NAME.2, NAME.3, ... (an existing run is never touched)."
function fresh_dir(name)
    d = joinpath(outdir, name); k = 1
    while ispath(d); k += 1; d = joinpath(outdir, "$name.$k"); end
    d
end

"The run directory above a checkpoint path (a conforming parent), searched up to lib/qs_runs; nothing for a flat file."
function parent_run_dir(path)
    d = dirname(abspath(path))
    while startswith(d, normpath(outdir) * "/")
        isfile(joinpath(d, "run.json")) && return d
        d = dirname(d)
    end
    nothing
end

function read_common!(P)       # COMMON scalars set from env/CLI or nonzero (the former QS_SET)
    for n in P.order
        startswith(n, "common.") || continue
        k = Symbol(n[8:end]); v = P[n]
        (v != P.specs[n].default || P.sources[n] == "cli" || startswith(P.sources[n], "env:")) || continue
        TF.setscalar!(k, eltype(TF.var(k))(v))
    end
end

if mode == "setup"
    prof = pos[2]
    p0 = params(; argv = base_argv)
    name = p0.io.name
    rc = start_run(fresh_dir("$name-setup"); problem = (name = "qstar-setup", version = "1", params = Dict("profile" => prof)),
                   name = name, params = p0)
    relabel!(rc)
    for n in RELAX_ONLY; p0[n]; end; p0["io.name"]
    try
        rin, nenv, natm = p0.setup.rin, p0.setup.nenv, p0.setup.natm
        install_julia!()
        grey = p0.setup.grey
        Mdot = p0.setup.mdot * Z0.Msun / Z0.yr                         # QS_MDOT in Msun/yr
        rout = p0.setup.rout; dr0 = p0.setup.dr0
        lgrid = p0.setup.lgrid
        Mdot_max = Mdot * p0.ramp.fac
        wcontact = p0.setup.wcontact; equi = p0.setup.equi
        lady = Tuple(parse.(Int, split(p0.setup.lady, ",")))
        q = QS.setup!(prof; rin, nenv, natm, grey, Mdot, rout, dr0, lgrid, Mdot_max, wcontact, equi, lady)
        @printf("quasi-star shell: M = %.4e g, L = %.4e erg/s, R = %.4e cm, r_in = %.4e (%.2f R), m_in = %.4e g, shell mass %.4e g\n",
                q.M, q.L, q.R, q.r_in, q.r_in / q.R, q.m_in, q.shell)
        Mdot > 0 && @printf("inflow: Ṁ = %.3e g/s, r_out = %.2f R, ρ_out = %.3e, T_out = %.0f, u_out = %.3e cm/s, xmext = %.3e\n",
                            Mdot, rout, TF.getscalar(:dr), TF.getscalar(:tr), TF.getscalar(:ur), TF.getscalar(:xmext))
        @printf("outer region: top r/R = %.5f, T_top = %.0f K, ρ_top = %.3e, τ_top = %.3e; %d zones (%d envelope, %d atmosphere)\n",
                q.atm.r[end] / q.R, q.atm.T[end], q.atm.ρ[end], q.atm.τtop, nenv + natm, nenv, natm)
        TF.save_state(initfile(name))
        # setup_meta.jls: the last set-up (read by qstar_analyze.jl by default); setup_meta_NAME.jls: this one, kept
        meta = (M = q.M, L = q.L, R = q.R, r_in = q.r_in)
        serialize(joinpath(outdir, "setup_meta.jls"), meta)
        serialize(joinpath(outdir, "setup_meta_$name.jls"), meta)
        println("saved ", initfile(name))
        output!(rc, initfile(name); kind = "snapshot", format = "tstate.gz")
        output!(rc, joinpath(outdir, "setup_meta_$name.jls"); kind = "other", format = "jls")
        finish!(rc; code = 0, label = "ok")
    catch e
        finish!(rc; code = 300, label = "failed", message = sprint(showerror, e)); rethrow()
    end
elseif mode == "relax"
    p0 = params(; argv = base_argv)                                  # preliminary: name, init, from (all restart-class)
    name = p0.io.name; init = isempty(p0.io.init) ? name : p0.io.init; from = p0.io.from
    parent = isempty(from) ? nothing : parent_run_dir(from)
    if parent !== nothing && isempty(p0.io.init)                     # a resume starts from its parent's set-up unless told otherwise
        pinit = get(get(RunConvention.parse_json(read(joinpath(parent, "run.json"), String)), "problem", Dict()), "params", Dict())
        init = String(get(pinit, "init", init))
    end
    TF.load_state!(initfile(init))                                   # io.init: start a new run NAME from another set-up
    ρ0x, φ0x, φ2x = TF.getscalar(:dr), TF.getscalar(:phir0), TF.getscalar(:phir2)   # ramp base = set-up inflow
    isempty(from) || TF.load_state!(from)                            # continue from a saved state (ramp keeps its clock)
    function shellmass()
        a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); rn = TF.var(:rn); dn = TF.var(:dn); un = TF.var(:un)
        ks = argmin([un[k+1] - un[k] for k in a:b]) + a - 1          # shock: largest outward velocity drop
        sum((4π / 3 * (rn[k+1]^3 - rn[k]^3) * dn[k] for k in a:ks-1); init = 0.0)   # 0 without a shock (hydrostatic)
    end
    # sink.m0 unset with a sink asked for: measure it from the loaded state, and record the number that was used (it is
    # what a later resume needs, and a run's parameters must be its own).
    derived = String[]; argv = base_argv
    if p0.sink.tau > 0 && p0.sink.m0 == 0
        argv = vcat(["--sink.m0=" * repr(shellmass())], base_argv); push!(derived, "sink.m0")
    end
    rundir = fresh_dir(name)
    if parent === nothing
        rc = start_run(rundir; problem = (name = "qstar-relax", version = "1", params = Dict("init" => init, "from" => from)),
                       name = name, params = params(; argv))
    else
        rc = resume(rundir; from, argv, name = name)
        println("resume: child of run ", basename(parent), " (", from, ")")
    end
    p = rc.params; relabel!(rc, derived)
    for n in SETUP_ONLY; p[n]; end                                   # parameters of the other mode, not read here
    p["io.name"]; p["io.init"]; p["io.from"]                         # read through the preliminary resolution above
    iobt = p.io.bt
    tmax = p.run.tmax; maxsteps = p.run.maxsteps
    install_julia!()
    TX.CONVECTION[] = true; TX.ALPHA_MLT[] = 0.8
    TX.GASROW[] = p.physics.gasrow                                   # gas-only energy row (see TitanResidFast.GASROW)
    TX.COMPTON[] = p.physics.compton                                 # Compton gas-radiation exchange (TitanResidFast.COMPTON)
    TX.MLT_NOINFALL[] = p.physics.mlt_noinfall                       # no MLT in the inflow (TitanResidFast.MLT_NOINFALL)
    TM.use!(TM.Z0Gas(; grey = p.setup.grey))
    pio = p.io
    TD.MSGIO[] = pio.msg ? stdout : devnull
    read_common!(p)
    TF.setscalar!(:tmax, tmax); TF.setscalar!(:jstep, Int32(1)); TF.setscalar!(:jstepe, Int32(10^7))
    L = TF.getscalar(:xlum)
    snaps = Any[snap()]
    t0 = time()
    # Inflow ramp (ramp.fac = target/initial Ṁ, ramp.t = e-folding time scale in s): the external density, and with it
    # the mass and energy fluxes through r_out, grow exponentially from the set-up value by ramp.fac over ramp.t·ln(fac).
    rampfac = p.ramp.fac; rampT = p.ramp.t
    function ramp!()
        rampfac > 1 || return
        f = min(rampfac, exp(TF.getscalar(:timen) / rampT))
        TF.setscalar!(:dr, ρ0x * f); TF.setscalar!(:phir0, φ0x * f); TF.setscalar!(:phir2, φ2x * f)
    end
    # Inner mass sink (sink.tau = relaxation time in s; 0: closed wall). The mass below the shock is held at its
    # initial value M_s0: each step the inner wall passes Ṁ_sink = max(0, M_s − M_s0)/τ onto the central point mass,
    # carrying the gas and radiation enthalpy of the innermost zone (TitanResidFast.INNER_SINK).
    sinkτ = p.sink.tau
    Ms0 = p.sink.m0
    sinkMdot = Ref(0.0)
    Mc0 = 4π * (TF.getscalar(:xmtot) - TF.var(:xmen)[Int(TF.getscalar(:ngrs))]) / Z0.Msun
    function sink!()
        sinkτ > 0 || return
        a = Int(TF.getscalar(:ngrs)); r = TF.var(:rn)[a]; ρ = TF.var(:dn)[a]
        Md = max(0.0, shellmass() - Ms0) / sinkτ; sinkMdot[] = Md
        φ0 = -Md / (4π * r^2)
        TF.setscalar!(:phil0, φ0); TF.setscalar!(:phil2, φ0 * (TF.var(:egn)[a] + TF.var(:pgn)[a] / ρ))
    end
    if sinkτ > 0
        TX.INNER_SINK[] = true; TX.GASROW[] || error("sink.tau (QS_SINK_TAU) needs physics.gasrow (QS_GASROW=1)")
        @printf("inner sink: τ = %.3e s, M_shell0 = %.6e g (below the shock)\n", sinkτ, Ms0); sink!()
    end
    # Time integrator (integrator.trbdf2: second-order TR-BDF2, Trbdf2.jl) and an optional advective Courant limit on the
    # next step (integrator.cfl = C: Δt ≤ C·min_k Δr_k/|u_k|), for time accuracy where TITAN's relative-change control allows more.
    p.integrator.trbdf2 && (TD.STEPFN[] = Trbdf2.step!; println("time integrator: TR-BDF2"))
    cfl = p.integrator.cfl
    dtmax = p.integrator.dtmax > 0 ? p.integrator.dtmax : Inf        # plain upper limit on Δt [s]
    function cfl!()
        TF.getscalar(:dtime) > dtmax && TF.setscalar!(:dtime, dtmax)
        cfl > 0 || return
        a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); rn = TF.var(:rn); un = TF.var(:un)
        dtc = minimum((rn[k+1] - rn[k]) / max(abs(un[k]), abs(un[k+1]), 1.0) for k in a:b) * cfl
        TF.getscalar(:dtime) > dtc && TF.setscalar!(:dtime, dtc)
    end
    # Checkpoint: the state into the run directory (indexed), and the old flat NAME_relaxed.tstate.gz beside it.
    curstep = Ref(0)
    flat_state = joinpath(outdir, "$(name)_relaxed.tstate.gz"); flat_snaps = joinpath(outdir, "$(name)_relax_snaps.jls")
    function checkpoint_now!(r)
        ck = joinpath(r.dir, @sprintf("checkpoint-%08d.tstate.gz", curstep[]))
        TF.save_state(ck); cp(ck, flat_state; force = true)
        checkpoint!(r, ck; step = curstep[], time = TF.getscalar(:timen))
        sn = joinpath(r.dir, @sprintf("snaps-%08d.jls", curstep[]))
        serialize(sn, snaps); cp(sn, flat_snaps; force = true)
        output!(r, sn; kind = "snapshot", format = "jls", step = curstep[], time = TF.getscalar(:timen))
    end
    rc.on_checkpoint = checkpoint_now!
    stopped = Ref(false)
    n, why = try TD.run_julia!(maxsteps = maxsteps, output = function (k)
        curstep[] = k
        TD.MSGIO[] = pio.msg ? stdout : devnull                      # io.msg is steerable
        ramp!(); sink!(); cfl!()
        every = pio.every                                            # io.every is steerable
        if k <= 3 || k % every == 0
            s = snap(); push!(snaps, s); report(s, L, @sprintf("step %5d", k))
            Mc = 4π * (TF.getscalar(:xmtot) - TF.var(:xmen)[Int(TF.getscalar(:ngrs))]) / Z0.Msun - Mc0
            sinkτ > 0 && @printf("          sink: Ṁ = %.4e Msun/yr, point mass gained %.4e Msun\n", sinkMdot[] * Z0.yr / Z0.Msun, Mc)
            flush(stdout)
            st = qstats(s, L)
            progress!(rc; step = k, time = s.t, dt = s.dt, fraction = clamp(s.t / tmax, 0.0, 1.0))
            vals = Dict{String,Any}("mach" => st.mach, "l_over_l0_min" => st.lmin, "l_over_l0_max" => st.lmax,
                                    "lconv_over_l_max" => st.lconv)
            sinkτ > 0 && (vals["sink_mdot"] = sinkMdot[] * Z0.yr / Z0.Msun; vals["sink_point_mass"] = Mc)
            metric!(rc; step = k, values = vals)
            k % (10every) == 0 && checkpoint_now!(rc)                # checkpoint, so long runs can be inspected and continued
        end
        if poll_control!(rc; step = k) == :stop                      # control/*.json: set, checkpoint, stop, pause
            stopped[] = true; TF.setscalar!(:tmax, 0.0)              # the driver's own end-of-run test ends the loop
        end
    end)
    catch e
        println("STOPPED: ", sprint(showerror, e)); iobt && Base.show_backtrace(stdout, catch_backtrace())
        push!(snaps, snap()); (length(snaps) - 1, e)
    end
    failed = why isa Exception
    stopped[] && (why = :stopped)
    @printf("relax: %d steps (%s) in %.1f s, t = %.4e s\n", n, failed ? :failed : why, time() - t0, TF.getscalar(:timen))
    TD.STEPFN[] === Trbdf2.step! && println("TR-BDF2 statistics: ", Trbdf2.STATS)
    serialize(flat_snaps, snaps)
    TF.save_state(flat_state)
    try
        if !failed
            curstep[] = n; checkpoint_now!(rc)                       # the end state, resumable
        end
    catch e
        println("final checkpoint failed: ", sprint(showerror, e)); global failed = true
        global why = e
    end
    finish!(rc; code = failed ? 300 : stopped[] ? 100 : 0, label = failed ? "failed" : stopped[] ? "stopped" : "ok",
            message = failed ? sprint(showerror, why) : nothing)
else
    error("mode: setup | relax")
end
