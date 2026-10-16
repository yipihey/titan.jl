# Phase 7d/e: quasi-star outer shell in pure-Julia TITAN.
#   TITAN_NOLIB=1 julia --project=env qstar_run.jl setup PROFILE.csv [rin nenv natm]   → reference_runs/init/QSTAR.tstate.gz
#   TITAN_NOLIB=1 julia --project=env qstar_run.jl relax [tmax_s maxsteps]           hydrostatic test (walls, no infall)
using Printf, Serialization
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl"))
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "QuasiStarSetup.jl"))
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
    q = QS.setup!(prof; rin, nenv, natm, grey, Mdot, rout, dr0, lgrid)
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
    TF.load_state!(initfile(name))
    install_julia!()
    TX.CONVECTION[] = true; TX.ALPHA_MLT[] = 0.8
    TM.use!(TM.Z0Gas(; grey = get(ENV, "QS_GREY", "0") == "1"))
    TD.MSGIO[] = get(ENV, "QS_MSG", "0") == "1" ? stdout : devnull
    TF.setscalar!(:tmax, tmax); TF.setscalar!(:jstep, Int32(1)); TF.setscalar!(:jstepe, Int32(10^7))
    L = TF.getscalar(:xlum)
    snaps = Any[snap()]
    t0 = time()
    n, why = try TD.run_julia!(maxsteps = maxsteps, output = function (k)
        every = parse(Int, get(ENV, "QS_EVERY", "100"))
        (k <= 3 || k % every == 0) || return
        s = snap(); push!(snaps, s); report(s, L, @sprintf("step %5d", k)); flush(stdout)
    end)
    catch e
        println("STOPPED: ", sprint(showerror, e)); get(ENV, "QS_BT", "0") == "1" && Base.show_backtrace(stdout, catch_backtrace())
        push!(snaps, snap()); (length(snaps) - 1, :failed)
    end
    @printf("relax: %d steps (%s) in %.1f s, t = %.4e s\n", n, why, time() - t0, TF.getscalar(:timen))
    serialize(joinpath(outdir, "$(name)_relax_snaps.jls"), snaps)
    TF.save_state(joinpath(outdir, "$(name)_relaxed.tstate.gz"))
else
    error("mode: setup | relax")
end
