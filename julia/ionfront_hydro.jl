# M3 test: Iliev et al. 2009 test 5, an H II region expanding with hydrodynamics, on TITAN's fixed mesh with the same pieces as the accretion run (IonRay, ChemNet, implicit upwind species).
# Pure H (CHEM_YHE=0), n_H = 1e-3, T = 100 K, u = 0, 5e48 photons/s of a 1e5 K blackbody, L = 15 kpc. No gray radiation: H line, recombination and free-free cooling explicit below (Iliev's list).
#   TITAN_NOLIB=1 CHEM_YHE=0 TITAN_MGR=600 julia --project=env ionfront_hydro.jl [nzones=256] [dtmax_Myr=0.5]
# Bars: paper's Fig. of r_I(t), v_I(t) and the profiles; the Spitzer law r_S (1 + 7 c_s t/(4 r_S))^(4/7) is only rough, as the paper says.
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl"))
include(joinpath(@__DIR__, "SCVH.jl")); include(joinpath(@__DIR__, "PopIIIEqEOS.jl"))
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl"))
include(joinpath(@__DIR__, "PopIIIMaterial.jl")); include(joinpath(@__DIR__, "IonRay.jl"))
using ForwardDiff
include(joinpath(@__DIR__, "Trbdf2.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const PM = PopIIIMaterial; const CN = ChemNet; const IR = IonRay
TD.STEPFN[] = Trbdf2.step!
nz = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 256; dtmaxMyr = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.5
const kpc = 3.0857e21; const Myr = 3.15576e13; const nH0 = 1e-3; const Ndot = 5e48; const χ = CN.χH
net, photos = ChemRates.modern()
envA = (Trad = 0.0, J = ν -> 0.0, static = true); envB = merge(envA, (; caseB = true))
const αB_tab = (T -> CN.photorates_t(photos[2], T, envB)[2])
# erg cm-3 s-1. The network's u holds χ per H⁺, so a recombination returns χ to the gas: Λ_rec carries it away as well as the electron's kT. Excitation: Cen 1992; free-free g = 1.3; thermal recombination: Cen 1992 case A shape.
function Λ(ρ, T, d)
    n = ρ / CN.μH(); ne = d[:e] * n; np = d[:Hp] * n; nh = d[:H] * n
    Tv = ForwardDiff.value(T)
    (ne * nh * 7.5e-19 * exp(-118348 / T) / (1 + sqrt(T / 1e5)) + ne * np * (αB_tab(Tv) * χ + 8.7e-27 * sqrt(T) * (T / 1e3)^-0.2 / (1 + (T / 1e6)^0.7)) + 1.3 * 1.42e-27 * sqrt(T) * ne * np) *
        (T > 50 ? one(T) : zero(T))
end
x0 = zeros(length(net.free)); x0[findfirst(==(:Hp), net.free)] = 1e-6
mat = PM.NetGas([copy(x0) for _ in 1:2000], net, photos, Λ; env = envB, pswitch = 1e30)
rI = collect(range(0.0, 15kpc, length = nz + 1)); ρ0 = nH0 * CN.μH()
av = (ql0 = 0.0, ql1 = 0.15, cq1 = 0.0, cq2 = 1.0)
PS.build!(; rI, rho = fill(ρ0, nz), T = fill(100.0, nz), u = zeros(nz + 1), geom = 2, lgrid = 2, inner = :wall, outer = :wall, selfgrav = false, material = mat, av, xmext = nothing)
TX.AV_ZONE[] = 2.0; TX.REGULAR_CENTRE[] = true; empty!(TX.XMINT[])
PS.install_julia!(BlockSolver)
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
for k in eachindex(mat.x); mat.x[k] = copy(x0); end
G = IR.groups(; species = [:H])
N0 = IR.blackbody(G, 1.0, 1e5); N0 .*= Ndot / sum(N0)
Q = zeros(length(mat.x))
PM.RAY[] = (; start = (a_, b_) -> (a_ - 1, N0), source = (ρ, r1, r2, Nin) -> IR.source(ρ, r1, r2, Nin, G, net), env = envB, Ec = G.Ec, Q, xmin = 0.0)
XP = Ref(copy(TF.var(:xmen)::Vector{Float64}))
function advect!()
    xme = TF.var(:xmen); old = XP[]; n = b - a + 1
    F = [xme[a+j-1] - old[a+j-1] for j in 1:n+1]; M = [xme[a+j-1] - xme[a+j] for j in 1:n]; M0 = [old[a+j-1] - old[a+j] for j in 1:n]
    lo = zeros(n); up = zeros(n); di = zeros(n)
    for j in 1:n
        di[j] = M[j] + max(F[j+1], 0.0) + max(-F[j], 0.0)
        j > 1 ? (lo[j] = -max(F[j], 0.0)) : (di[j] -= max(F[j], 0.0))
        j < n ? (up[j] = -max(-F[j+1], 0.0)) : (di[j] -= max(-F[j+1], 0.0))
    end
    cp_ = zeros(n); den = zeros(n); den[1] = di[1]; cp_[1] = up[1] / den[1]
    for j in 2:n; den[j] = di[j] - lo[j] * cp_[j-1]; cp_[j] = up[j] / den[j]; end
    xn = [M0[j] .* mat.x[a+j-1] for j in 1:n]; xn[1] = xn[1] ./ den[1]
    for j in 2:n; xn[j] = (xn[j] .- lo[j] .* xn[j-1]) ./ den[j]; end
    for j in n-1:-1:1; xn[j] = xn[j] .- cp_[j] .* xn[j+1]; end
    for j in 1:n; mat.x[a+j-1] = xn[j]; end
    copyto!(old, xme)
end
function front()
    r = TF.var(:rn); rc = [0.5 * (r[k] + r[k+1]) for k in a:b]; xp = [CN.full(mat.x[k], net)[:Hp] for k in a:b]
    j = findfirst(<(0.5), xp); j === nothing && return r[b+1]; j == 1 && return rc[1]
    rc[j-1] + (0.5 - xp[j-1]) / (xp[j] - xp[j-1]) * (rc[j] - rc[j-1])
end
const dtmax = dtmaxMyr * Myr; const tend = 500Myr
TF.setscalar!(:timen, 0.0); TF.setscalar!(:timeo, 0.0); TF.setscalar!(:dtime, 1e11); TF.setscalar!(:tmax, tend); TF.setscalar!(:jstepe, 10^7)
outs = [10.0, 30.0, 100.0, 200.0, 500.0]; io = Ref(1); tw = time()
rS = (3Ndot / (4π * 2.59e-13 * nH0^2))^(1 / 3)
@printf("  %5s %7s %8s %8s %8s %8s %8s %8s  %s\n", "step", "t[Myr]", "r_I[kpc]", "T(0.5r_I)", "x(0.5r_I)", "umax km/s", "Mass err", "dt[Myr]", "wall s")
M00 = sum(4π / 3 * ρ0 * (rI[j+1]^3 - rI[j]^3) for j in 1:nz)
res = TD.run_julia!(maxsteps = 200000, output = function (n)
    advect!(); dT = PM.chem!(mat)
    TF.setscalar!(:dtime, min(TF.getscalar(:dtime), dtmax, dT > 0.3 ? TF.getscalar(:dtime) * 0.7 : Inf))
    t = TF.getscalar(:timen)
    if io[] <= length(outs) && t >= outs[io[]] * Myr * (1 - 1e-6)
        r = front(); rn = TF.var(:rn); un = TF.var(:un); tn = TF.var(:tn); dn = TF.var(:dn)
        j = a + searchsortedfirst(rn[a:b], 0.5r) - 2
        Mt = sum(4π / 3 * dn[k] * (rn[k+1]^3 - rn[k]^3) for k in a:b)
        @printf("  %5d %7.1f %8.3f %8.0f %8.4f %8.2f %8.1e %8.2f  %.0f\n", n, t / Myr, r / kpc, tn[j], CN.full(mat.x[j], net)[:Hp], maximum(un[a:b+1]) / 1e5, Mt / M00 - 1, TF.getscalar(:dtime) / Myr, time() - tw)
        io[] += 1
    end
end)
open(joinpath(get(ENV, "OUTDIR", "."), "ionfront_hydro_profile_$(nz).txt"), "w") do f
    rn = TF.var(:rn); un = TF.var(:un); tn = TF.var(:tn); dn = TF.var(:dn)
    println(f, "# r/L n_H T u[km/s] x_HII at t = ", TF.getscalar(:timen) / Myr, " Myr")
    for k in a:b; @printf(f, "%.4f %.4e %.1f %.3f %.5f\n", 0.5 * (rn[k] + rn[k+1]) / 15kpc, dn[k] / CN.μH(), tn[k], un[k] / 1e5, CN.full(mat.x[k], net)[:Hp]); end
end
println("ended ", res, "; last t = ", TF.getscalar(:timen) / Myr, " Myr, rS(2.59e-13) = ", rS / kpc, " kpc")
