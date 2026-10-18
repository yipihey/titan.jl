# Phase 7d/e: quasi-star outer shell in pure-Julia TITAN.
#   TITAN_NOLIB=1 julia --project=env qstar_run.jl setup PROFILE.csv [rin nenv natm]   → reference_runs/init/QSTAR.tstate.gz
#   TITAN_NOLIB=1 julia --project=env qstar_run.jl relax [tmax_s maxsteps]           run (hydrostatic test, or infall if set up)
#   relax options (env): QS_NAME, QS_INIT (set-up to start from), QS_FROM (saved state to continue), QS_SET (scalars),
#     QS_GASROW / QS_COMPTON (=1: gas-only energy row / Compton exchange, see TitanResidFast), QS_RAMP_FAC/QS_RAMP_T, QS_EVERY,
#     QS_SINK_TAU (inner mass sink relaxation time [s]; with QS_FROM give QS_SINK_M0 [g] = the original shell mass),
#     QS_TRBDF2=1 (second-order time stepping), QS_CFL (advective Courant cap on Δt), QS_DTMAX (cap on Δt [s])
using Printf, Serialization
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl"))
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "QuasiStarSetup.jl")); include(joinpath(@__DIR__, "Trbdf2.jl"))
using .TitanDriver
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

function report(s, L, label)
    a = s.a; cs = [sqrt(4 / 3 * (s.pg[i] + s.er[i] / 3) / s.d[i]) for i in eachindex(s.d)]
    mach = maximum(abs(s.u[i]) / cs[min(i, end)] for i in eachindex(s.u))
    Ltot = s.Lrad .+ s.Lconv
    @printf("%s t = %.4e s  dt = %.3e  max|u|/cs = %.3e  L/L0: min %.4f max %.4f (interior)  Lconv/L max %.3f\n",
            label, s.t, s.dt, mach, minimum(Ltot[2:end-1]) / L, maximum(Ltot[2:end-1]) / L, maximum(s.Lconv) / L)
end

mode = ARGS[1]
if mode == "setup"
    prof = ARGS[2]
    rin = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.4
    nenv = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 140
    natm = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 60
    install_julia!()
    grey = get(ENV, "QS_GREY", "0") == "1"
    Mdot = parse(Float64, get(ENV, "QS_MDOT", "0")) * Z0.Msun / Z0.yr      # QS_MDOT in Msun/yr
    rout = parse(Float64, get(ENV, "QS_ROUT", "2.0")); dr0 = parse(Float64, get(ENV, "QS_DR0", "5e-5"))
    lgrid = parse(Int, get(ENV, "QS_LGRID", "2"))
    Mdot_max = Mdot * parse(Float64, get(ENV, "QS_RAMP_FAC", "1"))
    wcontact = parse(Float64, get(ENV, "QS_WCONTACT", "0")); equi = get(ENV, "QS_EQUI", "0") == "1"
    lady = Tuple(parse.(Int, split(get(ENV, "QS_LADY", "2,2,2"), ",")))
    q = QS.setup!(prof; rin, nenv, natm, grey, Mdot, rout, dr0, lgrid, Mdot_max, wcontact, equi, lady)
    @printf("quasi-star shell: M = %.4e g, L = %.4e erg/s, R = %.4e cm, r_in = %.4e (%.2f R), m_in = %.4e g, shell mass %.4e g\n",
            q.M, q.L, q.R, q.r_in, q.r_in / q.R, q.m_in, q.shell)
    Mdot > 0 && @printf("inflow: Ṁ = %.3e g/s, r_out = %.2f R, ρ_out = %.3e, T_out = %.0f, u_out = %.3e cm/s, xmext = %.3e\n",
                        Mdot, rout, TF.getscalar(:dr), TF.getscalar(:tr), TF.getscalar(:ur), TF.getscalar(:xmext))
    @printf("outer region: top r/R = %.5f, T_top = %.0f K, ρ_top = %.3e, τ_top = %.3e; %d zones (%d envelope, %d atmosphere)\n",
            q.atm.r[end] / q.R, q.atm.T[end], q.atm.ρ[end], q.atm.τtop, nenv + natm, nenv, natm)
    TF.save_state(initfile(get(ENV, "QS_NAME", "QSTAR")))
    serialize(joinpath(outdir, "setup_meta.jls"), (M = q.M, L = q.L, R = q.R, r_in = q.r_in))
    println("saved ", initfile(get(ENV, "QS_NAME", "QSTAR")))
elseif mode == "relax"
    tmax = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 3e9
    maxsteps = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 400
    name = get(ENV, "QS_NAME", "QSTAR")
    TF.load_state!(initfile(get(ENV, "QS_INIT", name)))           # QS_INIT: start a new run NAME from another set-up
    ρ0x, φ0x, φ2x = TF.getscalar(:dr), TF.getscalar(:phir0), TF.getscalar(:phir2)   # ramp base = set-up inflow
    haskey(ENV, "QS_FROM") && TF.load_state!(ENV["QS_FROM"])     # continue from a saved state (ramp keeps its clock)
    install_julia!()
    TX.CONVECTION[] = true; TX.ALPHA_MLT[] = 0.8
    TX.GASROW[] = get(ENV, "QS_GASROW", "0") == "1"           # gas-only energy row (see TitanResidFast.GASROW)
    TX.COMPTON[] = get(ENV, "QS_COMPTON", "0") == "1"         # Compton gas-radiation exchange (TitanResidFast.COMPTON)
    TX.MLT_NOINFALL[] = get(ENV, "QS_MLT_NOINFALL", "0") == "1"   # no MLT in the inflow (TitanResidFast.MLT_NOINFALL)
    TM.use!(TM.Z0Gas(; grey = get(ENV, "QS_GREY", "0") == "1"))
    TD.MSGIO[] = get(ENV, "QS_MSG", "0") == "1" ? stdout : devnull
    for kv in split(get(ENV, "QS_SET", ""), ","; keepempty = false)      # e.g. QS_SET=ql1=1e-4,cq2=1
        k, val = split(kv, "="); T = eltype(TF.var(Symbol(k)))
        TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
    end
    TF.setscalar!(:tmax, tmax); TF.setscalar!(:jstep, Int32(1)); TF.setscalar!(:jstepe, Int32(10^7))
    L = TF.getscalar(:xlum)
    snaps = Any[snap()]
    t0 = time()
    # Inflow ramp (QS_RAMP_FAC = target/initial Ṁ, QS_RAMP_T = e-folding time scale in s): the external density, and with it
    # the mass and energy fluxes through r_out, grow exponentially from the set-up value by QS_RAMP_FAC over QS_RAMP_T·ln(fac).
    rampfac = parse(Float64, get(ENV, "QS_RAMP_FAC", "1")); rampT = parse(Float64, get(ENV, "QS_RAMP_T", "1e6"))
    function ramp!()
        rampfac > 1 || return
        f = min(rampfac, exp(TF.getscalar(:timen) / rampT))
        TF.setscalar!(:dr, ρ0x * f); TF.setscalar!(:phir0, φ0x * f); TF.setscalar!(:phir2, φ2x * f)
    end
    # Inner mass sink (QS_SINK_TAU = relaxation time in s; unset: closed wall). The mass below the shock is held at its
    # initial value M_s0: each step the inner wall passes Ṁ_sink = max(0, M_s − M_s0)/τ onto the central point mass,
    # carrying the gas and radiation enthalpy of the innermost zone (TitanResidFast.INNER_SINK).
    sinkτ = parse(Float64, get(ENV, "QS_SINK_TAU", "0"))
    function shellmass()
        a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); rn = TF.var(:rn); dn = TF.var(:dn); un = TF.var(:un)
        ks = argmin([un[k+1] - un[k] for k in a:b]) + a - 1          # shock: largest outward velocity drop
        sum(4π / 3 * (rn[k+1]^3 - rn[k]^3) * dn[k] for k in a:ks-1)
    end
    Ms0 = get(ENV, "QS_SINK_M0", "") == "" ? shellmass() : parse(Float64, ENV["QS_SINK_M0"])
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
        TX.INNER_SINK[] = true; TX.GASROW[] || error("QS_SINK_TAU needs QS_GASROW=1")
        @printf("inner sink: τ = %.3e s, M_shell0 = %.6e g (below the shock)\n", sinkτ, Ms0); sink!()
    end
    # Time integrator (QS_TRBDF2=1: second-order TR-BDF2, Trbdf2.jl) and an optional advective Courant limit on the
    # next step (QS_CFL = C: Δt ≤ C·min_k Δr_k/|u_k|), for time accuracy where TITAN's relative-change control allows more.
    get(ENV, "QS_TRBDF2", "0") == "1" && (TD.STEPFN[] = Trbdf2.step!; println("time integrator: TR-BDF2"))
    cfl = parse(Float64, get(ENV, "QS_CFL", "0"))
    dtmax = parse(Float64, get(ENV, "QS_DTMAX", "Inf"))          # plain upper limit on Δt [s]
    function cfl!()
        TF.getscalar(:dtime) > dtmax && TF.setscalar!(:dtime, dtmax)
        cfl > 0 || return
        a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); rn = TF.var(:rn); un = TF.var(:un)
        dtc = minimum((rn[k+1] - rn[k]) / max(abs(un[k]), abs(un[k+1]), 1.0) for k in a:b) * cfl
        TF.getscalar(:dtime) > dtc && TF.setscalar!(:dtime, dtc)
    end
    n, why = try TD.run_julia!(maxsteps = maxsteps, output = function (k)
        ramp!(); sink!(); cfl!()
        every = parse(Int, get(ENV, "QS_EVERY", "100"))
        (k <= 3 || k % every == 0) || return
        s = snap(); push!(snaps, s); report(s, L, @sprintf("step %5d", k))
        sinkτ > 0 && @printf("          sink: Ṁ = %.4e Msun/yr, point mass gained %.4e Msun\n", sinkMdot[] * Z0.yr / Z0.Msun,
                              4π * (TF.getscalar(:xmtot) - TF.var(:xmen)[Int(TF.getscalar(:ngrs))]) / Z0.Msun - Mc0)
        flush(stdout)
        if k % (10every) == 0          # checkpoint, so long runs can be inspected and continued (QS_FROM)
            TF.save_state(joinpath(outdir, "$(name)_relaxed.tstate.gz")); serialize(joinpath(outdir, "$(name)_relax_snaps.jls"), snaps)
        end
    end)
    catch e
        println("STOPPED: ", sprint(showerror, e)); get(ENV, "QS_BT", "0") == "1" && Base.show_backtrace(stdout, catch_backtrace())
        push!(snaps, snap()); (length(snaps) - 1, :failed)
    end
    @printf("relax: %d steps (%s) in %.1f s, t = %.4e s\n", n, why, time() - t0, TF.getscalar(:timen))
    TD.STEPFN[] === Trbdf2.step! && println("TR-BDF2 statistics: ", Trbdf2.STATS)
    serialize(joinpath(outdir, "$(name)_relax_snaps.jls"), snaps)
    TF.save_state(joinpath(outdir, "$(name)_relaxed.tstate.gz"))
else
    error("mode: setup | relax")
end
