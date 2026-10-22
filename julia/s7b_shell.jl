# S7b test c: radiation-pressure shell in TITAN. Included by s7b_test.jl (needs CN, RT, photos, check, TITAN not yet loaded).
# A cold shell of fully ionised hydrogen (Thomson + free-free opacity of ContinuumOpacity) around an empty hot cavity, driven by a central monochromatic source L through
# RT.transfer (direct beam, diffuse = false) and the momentum-row hook TitanResidFast.AEXT (interface_accel of the zone accelerations, updated after every step, so the force of a step is
# that of the previous state). Checks: (1) total momentum gained by the shell = ∫ Σ a_j M_j dt (the hook delivers what RadTransfer deposits), (2) = (L_in − L_out)/c × t (momentum balance of the
# beam), thin (τ ≈ 0.01) and thick (τ ≈ 60) shells, (3) acceleration of the thin shell = κ L/(4π r² c).
include(joinpath(J, "TitanDriver.jl")); include(joinpath(J, "JuliaPhysics.jl")); include(joinpath(J, "TitanResidFast.jl")); include(joinpath(J, "BlockSolver.jl"))
include(joinpath(J, "TitanCore.jl")); include(joinpath(J, "Z0Physics.jl")); include(joinpath(J, "TitanMaterials.jl")); include(joinpath(J, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TM = TitanMaterials
mutable struct IdealIon <: TM.Material; μ::Float64; end
@inline function TX.material_gas(m::IdealIon, ρ, T)
    p = ρ * CN.kB * T / (m.μ * CN.mH); (p, 1.5 * p / ρ, sqrt(5.0 / 3.0 * p / ρ))
end
function TM.fill_eos!(m::IdealIon)
    dn, tn, pgn, egn, asn = TM.v(:dn), TM.v(:tn), TM.v(:pgn), TM.v(:egn), TM.v(:asn)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        pg, eg, as_ = TX.material_gas(m, dn[k], tn[k]); pgn[k] = pg; egn[k] = eg; asn[k] = as_
    end
    nothing
end
TM.fill_opac!(::IdealIon) = nothing

function shell_run(; ρsh, L, tend, hook = true, nvoid = 20, nsh = 60, r1 = 1e15, r2 = 2e15, Eph = 3.0)
    g = RT.edge_grid([2.0, 3.0, 4.0]); ab = RT.absorbers(photos, g)
    src = RT.Source([0.0, L], 0.0, L)      # all power in the bin 3-4 eV
    Tsh = 100.0; μ = 0.5; psh = ρsh * CN.kB * Tsh / (μ * CN.mH); Tv = 1e5; ρv = psh * μ * CN.mH / (CN.kB * Tv)
    rI = vcat(collect(range(0, r1, length = nvoid + 1))[1:end-1], collect(range(r1, r2, length = nsh + 1)))
    rho = vcat(fill(ρv, nvoid), fill(ρsh, nsh)); Tz = vcat(fill(Tv, nvoid), fill(Tsh, nsh)); u0 = zeros(nvoid + nsh + 1)
    mat = IdealIon(μ)
    av = (ql0 = 0.0, ql1 = 0.0, cq1 = 0.0, cq2 = 1.0)
    PS.build!(; rI, rho, T = Tz, u = u0, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = psh, selfgrav = false, material = mat, av)
    TX.AV_ZONE[] = 1.5; TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); TX.PEXT[] = TF.var(:pgn)[b]
    TF.setscalar!(:thet, 0.55); TC = TitanCore; TC.DT_DYN[] = 0.05
    TF.setscalar!(:dtime, 1e3); TF.setscalar!(:tmax, tend)
    n = b - a + 1; Wgain = Ref(0.0); Pgain = Ref(0.0); Pexp_beam = Ref(0.0); aprev = Ref(zeros(n + 1)); Mif = Ref(zeros(n + 1)); tlast = Ref(0.0); amax = Ref(0.0); nsteps = Ref(0)
    P_now() = (un = TF.var(:un); dn = TF.var(:dn); dv = TF.var(:dvoln); sum(4π * dn[k] * dv[k] * 0.5 * (un[k] + un[k+1]) for k in a:b))
    TX.AEXT[] = Float64[]
    Lfrac = Ref(0.0); uprev = Ref(zeros(n + 1))
    KE_now() = (un = TF.var(:un); dn = TF.var(:dn); dv = TF.var(:dvoln); sum(4π * dn[k] * dv[k] * 0.5 * (0.5 * un[k]^2 + 0.5 * un[k+1]^2) for k in a:b))
    function outcb(k)
        t = TF.getscalar(:timen); dt = t - tlast[]; tlast[] = t
        # force applied during this step (the AEXT of the step that just ran)
        Pgain[] += sum(aprev[] .* Mif[]) * dt; Wgain[] += sum(aprev[] .* Mif[] .* (0.5 .* (TF.var(:un)[a:b+1] .+ uprev[]))) * dt; uprev[] = copy(TF.var(:un)[a:b+1]); Pexp_beam[] += Lfrac[] * L / c * dt; nsteps[] += 1
        rn = TF.var(:rn); dn = TF.var(:dn); tn = TF.var(:tn); dv = TF.var(:dvoln)
        med = RT.Medium(rn[a:b+1], tn[a:b], dn[a:b] ./ CN.mH, Dict{Symbol,Vector{Float64}}(:Hp => fill(1.0, n), :e => fill(1.0, n)))
        mz = [4π * dn[j] * dv[j] for j in a:b]
        rf = RT.transfer(g, ab, src, med, RT.Opts(diffuse = false); rho = dn[a:b])
        Lfrac[] = (rf.ledger.Esrc - rf.ledger.Edir_out) / L
        if hook
            acc = RT.interface_accel(rf.arad, mz); full = zeros(length(TF.var(:un))); full[a:b+1] .= acc; TX.AEXT[] = full
            aprev[] = acc
            # masses on the interfaces in the units of the momentum row: Σ a_j M_j with M_j = ½(m_{j-1} + m_j)  [m = 4π ρ dvol]
            Mi = zeros(n + 1); Mi[1] = 0.5 * mz[1]; Mi[n+1] = 0.5 * mz[n]; for j in 2:n; Mi[j] = 0.5 * (mz[j-1] + mz[j]); end
            Mif[] = Mi
            amax[] = acc[nvoid+2]
        end
    end
    TD.run_julia!(maxsteps = 20000, output = outcb)
    (; P = P_now(), KE = KE_now(), W = Wgain[], Pgain = Pgain[], Pbeam = Pexp_beam[], steps = nsteps[], a = amax[], t = TF.getscalar(:timen), Lfrac = Lfrac[], rn = copy(TF.var(:rn)[a:b+1]), un = copy(TF.var(:un)[a:b+1]), dn = copy(TF.var(:dn)[a:b]))
end

function test_c()
    L = 1e43
    for (name, ρsh) in (("thin", 1e-17), ("thick", 2e-13))
        tend = name == "thin" ? 1e6 : 3e5
        o0 = shell_run(; ρsh, L, tend, hook = false)
        o1 = shell_run(; ρsh, L, tend, hook = true)
        ΔP = o1.P - o0.P
        # Thomson: κ = 0.4 σ... the shell optical depth from the transfer's own opacity (Lfrac = 1 − e^{-τ} for the beam)
        @printf("c %-5s: steps %d, t = %.2e s, τ_eff: absorbed fraction %.4f, control drift P0 = %.3e, ΔP = %.4e, hook ∫ΣaM dt = %.4e, beam (L/c) Σ f dt = %.4e, a_max = %.3e cm/s²\n", name, o1.steps, o1.t, o1.Lfrac, o0.P, ΔP, o1.Pgain, o1.Pbeam, o1.a)
        check("c1 hook momentum ($name)", abs(ΔP / o1.Pgain - 1) < 5e-3, @sprintf("Σ m u gained %.4e vs ∫ Σ a_j M_j dt %.4e (ratio %.5f)", ΔP, o1.Pgain, ΔP / o1.Pgain))
        check("c2 beam balance ($name)", abs(o1.Pgain / o1.Pbeam - 1) < 2e-2, @sprintf("∫ Σ a M dt vs ∫ (L_in − L_out)/c dt: ratio %.5f", o1.Pgain / o1.Pbeam))
        # energy: the work of the radiation force on the gas (kinetic energy gained) is a fraction ~v/c of the absorbed radiation; heating/ionisation of the absorbed power is accounted by the transfer ledger
        @printf("   %s: kinetic energy gained %.4e erg, ∫ Σ a M ū dt = %.4e erg (ratio %.4f), absorbed radiation energy %.4e erg: work/absorbed = %.2e\n", name, o1.KE - o0.KE, o1.W, (o1.KE - o0.KE) / o1.W, o1.Lfrac * L * o1.t, o1.W / (o1.Lfrac * L * o1.t))
        check("c4 work of the force ($name)", 0.8 < (o1.KE - o0.KE) / o1.W < 1.02, @sprintf("ΔKE / ∫ΣaMū dt = %.4f (≤ 1: the rest is dissipated by the artificial viscosity in the accelerated skin of the thick shell)", (o1.KE - o0.KE) / o1.W))
        if name == "thin"
            κ = 0.2 * (1 + 0.0)     # electron scattering of fully ionised H: σ_T/m_p = 0.4006 cm²/g
            aexp = 0.4006 * L / (4π * (1.0e15 + 1.5 * 1e15 / 60)^2 * c)       # at the first interface inside the shell
            check("c3 thin acceleration", abs(o1.a / aexp - 1) < 0.15, @sprintf("a_max = %.3e vs σ_T L/(4π m_p r1² c) = %.3e (ratio %.3f; the free-free adds ~%.0f%%)", o1.a, aexp, o1.a / aexp, 100 * (o1.a / aexp - 1)))
        end
    end
end
