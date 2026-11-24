# Noh problem (spherical, γ = 5/3): cold uniform gas flowing into a reflecting wall, the high-Mach limit of the accretion shock (K/e ≫ 1 ahead of the shock).
#   TITAN_NOLIB=1 [ECOMPAT=1] julia --project=env noh_test.jl [n] [T0]
# Exact solution inside the shock (r < |u0| t/3): u = 0, ρ = ((γ+1)/(γ−1))³ ρ0 = 64 ρ0, e = u0²/2 (all the kinetic energy thermalised), p = (γ−1) ρ e;
# ahead of it ρ = ρ0 (1 + |u0| t/r)², u = u0. Hydrodynamics only (neqn = 5), Eulerian fixed mesh, inner wall at r = RIN, outer inflow (ρ0, T0, u0).
# The run ends at t = 0.6 R/|u0|: gas at r < R − |u0| t (0.4 R) is initial gas, so the exact pre-shock state holds there (outside it the boundary's constant ρ0 flows in). Output: shock position and plateau values against the exact ones.
# ECOMPAT=1: the compatible internal-energy rows (TitanResidFast.ECOMPAT); unset: TITAN's own rows.
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TM = TitanMaterials
# INTEGRATOR=trbdf2 (needs DTERR=rtol in the environment: error-controlled TR-BDF2, as the production stack) | theta (default: TITAN's theta step)
const INTEG = get(ENV, "INTEGRATOR", "theta")
INTEG == "trbdf2" && include(joinpath(@__DIR__, "Trbdf2.jl"))
# AV=prod: production artificial viscosity (AVZONE=2, AVTENSION=0, ql0 = 0, ql1 = 0.15, cq1 = 0, cq2 = 1); AV=off: cq1 = cq2 = 0; unset: the state's defaults (effectively off)
avm = get(ENV, "AV", "")
av = avm == "prod" ? (ql0 = 0.0, ql1 = 0.15, cq1 = parse(Float64, get(ENV, "AVCQ1", "0")), cq2 = parse(Float64, get(ENV, "AVCQ2", "1"))) :
     avm == "off" ? (ql0 = 0.0, ql1 = 0.0, cq1 = 0.0, cq2 = 0.0) : (;)
n = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 200; T0 = length(ARGS) > 1 ? parse(Float64, ARGS[2]) : 1.0
γ = 5 / 3; μg = 1.0; R = 1e10; RIN = 0.01R; ρ0 = 1e-10; u0 = -1e6
rI = collect(range(RIN, R, length = n + 1)); uI = fill(u0, n + 1); uI[1] = 0.0         # the wall interface is at rest
PS.build!(; rI, rho = fill(ρ0, n), T = fill(T0, n), u = uI, geom = 2, lgrid = 2, inner = :wall, outer = :inflow, inflow = (ρ0, T0, u0), gam = γ, gmu = μg, av)
avm == "prod" && (TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "2")); TX.AV_TENSION[] = get(ENV, "AVTENSION", "0") == "1")
INTEG == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
PS.install_julia!(BlockSolver); TD.MSGIO[] = get(ENV, "NOHMSG", "0") == "1" ? stdout : devnull
etot = get(ENV, "ECOMPAT", "0")
etot == "1" && (TX.ECOMPAT[] = true; TX.GASROW[] = true)
mon = Any[]; get(ENV, "ETOTMON", "0") != "0" && (TD.STAGE_HOOK[] = (kind, dt) -> (push!(mon, TX.etot_monitor(dt)); nothing))
ck, cm0 = TF.getscalar(:ck), TF.getscalar(:cm0)
Tpost = (γ - 1) * μg * cm0 * 0.5u0^2 / ck                       # e = (k/(m0 μ)) T/(γ−1) = u0²/2
tend = 0.6R / abs(u0)
TF.setscalar!(:dtime, 1e-5R / abs(u0)); TF.setscalar!(:tmax, tend)
haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
nst, why = TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")))
t = TF.getscalar(:timen); r = TF.var(:rn)[a:b+1]; rc = 0.5 .* (r[1:end-1] .+ r[2:end]); d = TF.var(:dn)[a:b]; T = TF.var(:tn)[a:b]; u = TF.var(:un)[a:b+1]
rs = abs(u0) * t / 3
# shock position: first zone from the outside where ρ exceeds 32 ρ0 (halfway between the pre-shock 16 ρ0 and the post-shock 64 ρ0, scaled)
is = findlast(>(32ρ0), d); rsh = rc[is]
sel = findall(x -> 0.3rs < x < 0.8rs, rc)
@printf("Noh n=%d T0=%g K (Mach %.0f, K/e %.0f) ECOMPAT=%s: %d steps (%s), t = %.4e s, tend %.4e\n", n, T0, abs(u0) / sqrt(γ * ck * T0 / (cm0 * μg)),
        0.5u0^2 / (ck * T0 / (cm0 * μg * (γ - 1))), etot == "" ? "unset" : etot, nst, why, t, tend)
@printf("shock radius (ρ > 32 ρ0): %.4e, exact |u0|t/3 = %.4e, ratio %.4f\n", rsh, rs, rsh / rs)
@printf("plateau 0.3-0.8 r_s (%d zones): T/T_RH median %.4f min %.4f max %.4f; ρ/(64ρ0) median %.4f min %.4f max %.4f; max|u|/|u0| %.3e\n", length(sel),
        sort(T[sel] ./ Tpost)[end÷2+1], minimum(T[sel] ./ Tpost), maximum(T[sel] ./ Tpost), sort(d[sel] ./ (64ρ0))[end÷2+1], minimum(d[sel] ./ (64ρ0)), maximum(d[sel] ./ (64ρ0)),
        maximum(abs.(u[sel])) / abs(u0))
sp = findall(x -> 1.3rs < x < R - abs(u0) * t, rc)
@printf("pre-shock 1.3 r_s - (R - |u0| t) (%d zones): T min %.3e max %.3e K (T0 %.3g); ρ/ρ_exact median %.4f min %.4f max %.4f\n", length(sp), minimum(T[sp]), maximum(T[sp]), T0,
        sort(d[sp] ./ (ρ0 .* (1 .+ abs(u0) * t ./ rc[sp]) .^ 2))[end÷2+1], minimum(d[sp] ./ (ρ0 .* (1 .+ abs(u0) * t ./ rc[sp]) .^ 2)), maximum(d[sp] ./ (ρ0 .* (1 .+ abs(u0) * t ./ rc[sp]) .^ 2)))
@printf("T_RH = %.4e K; wall zone T/T_RH %.3f, ρ/(64ρ0) %.3f\n", Tpost, T[1] / Tpost, d[1] / (64ρ0))
if !isempty(mon)
    Gs = sort(abs.([m.G for m in mon])); @printf("ETOTMON %d solves: median|G| %.3e max|G| %.3e; sum imb/sum thr %.3e\n", length(mon), Gs[(length(Gs)+1)÷2], Gs[end], sum(m.imb for m in mon) / sum(m.thr for m in mon))
end
