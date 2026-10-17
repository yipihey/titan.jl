# Switch a saved Eulerian quasi-star state to TITAN's adaptive grid (lgrid = 1) for continuation.
#   TITAN_NOLIB=1 julia --project=env to_adaptive.jl IN.tstate.gz OUT.tstate.gz [gridtau [lady [weights]]]
# The grid parameters are those of QuasiStarSetup (ladx = 1 with xscale = R, alph = 1.5, ibet = 1, resolution weights
# lady = "2,2,2": exterior mass, ln ρ, ln T). TITAN's grid equation forces the first two and the last two zones to equal
# width; the two boundary pairs are made equal here by a conservative re-split (mass, gas and radiation energy,
# velocities unchanged), so the first step does not have to move them.
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
using Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const PS = ProblemSetup
const TC = TitanCore; const Z0 = Z0Physics
TF.load_state!(ARGS[1]); PS.install_julia!(BlockSolver); TM.use!(TM.Z0Gas())
gridtau = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e3
lady = length(ARGS) >= 4 ? parse.(Int, split(ARGS[4], ",")) : [2, 2, 2]
wts = length(ARGS) >= 5 ? parse.(Float64, split(ARGS[5], ",")) : [1.0 for _ in lady]     # weights of the lady terms
v(s) = TF.var(s)::Vector{Float64}; s!(k, x) = (a = TF.var(k); a[] = eltype(a)(x))
a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
rn, dn, tn, un, ern, frn = v(:rn), v(:dn), v(:tn), v(:un), v(:ern), v(:frn)
vol(r1, r2) = (r2^3 - r1^3) / 3
"Re-split zones k, k+1 into two of equal width, conserving mass, gas energy and radiation energy (u unchanged)."
function resplit!(k)
    r1, r2, r3 = rn[k], rn[k+1], rn[k+2]; V1, V2 = vol(r1, r2), vol(r2, r3)
    m = dn[k] * V1 + dn[k+1] * V2
    e = dn[k] * V1 * Z0.gas_e(dn[k], tn[k]) + dn[k+1] * V2 * Z0.gas_e(dn[k+1], tn[k+1])
    E = ern[k] * V1 + ern[k+1] * V2
    rm = 0.5(r1 + r3); W1, W2 = vol(r1, rm), vol(rm, r3)
    @printf("zones %d,%d: interface %.6e -> %.6e (%.2e relative)\n", k, k + 1, r2, rm, rm / r2 - 1)
    # keep the profile shape: scale both zones' ρ by one factor (mass), both T by one factor (gas energy), both E by
    # one factor (radiation energy). Equal values in the two zones would make ∇ = d ln T/d ln P there 0/0 (MLT).
    fρ = m / (dn[k] * W1 + dn[k+1] * W2); ρ1, ρ2 = dn[k] * fρ, dn[k+1] * fρ
    eof(f) = ρ1 * W1 * Z0.gas_e(ρ1, tn[k] * f) + ρ2 * W2 * Z0.gas_e(ρ2, tn[k+1] * f)
    lo, hi = log(0.5), log(2.0)
    for _ in 1:100; mid = 0.5(lo + hi); eof(exp(mid)) > e ? (hi = mid) : (lo = mid); end
    fT = exp(0.5(lo + hi)); fE = E / (ern[k] * W1 + ern[k+1] * W2)
    rn[k+1] = rm; dn[k], dn[k+1] = ρ1, ρ2; tn[k], tn[k+1] = tn[k] * fT, tn[k+1] * fT; ern[k], ern[k+1] = ern[k] * fE, ern[k+1] * fE
end
resplit!(a); resplit!(b - 1)
s!(:lgrid, 1); s!(:ladx, 1); s!(:xscale, 8.4352e15); s!(:alph, 1.5); s!(:ibet, 1); s!(:tau, gridtau)
la = TF.var(:lady); wt = TF.var(:wt); ys = TF.var(:yscl); la .= 0; wt .= 0.0
for (i, l) in enumerate(lady); la[i] = l; wt[i] = l > 0 ? wts[i] : 0.0; end
ys[1] = TF.getscalar(:xmtot)
# derived quantities, masses, old = new (as ProblemSetup.build!)
mu = 2; xmn = v(:xmn); xmen = v(:xmen); xmtot = TF.getscalar(:xmtot)
for k in a+1:b+1; xmn[k] = xmn[k-1] + dn[k-1] * vol(rn[k-1], rn[k]); end
for k in a:b+1; xmen[k] = xmtot - xmn[k]; end
for k in a-1:b+2; v(:rmun)[k] = rn[k]^2; v(:rmup1n)[k] = rn[k]^3; v(:rmum1n)[k] = rn[k]; end
for k in a-1:b+1; v(:dvoln)[k] = vol(rn[k], rn[k+1]); end
s!(:delrl, rn[a+1] - rn[a]); s!(:delrr, rn[b+2] - rn[b+1])
TD.EOS[](); TD.OPAC[]()
rng = a-1:b+2
for (o, n) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(n), rng)); end
for (o, n) in TC.PEQ_OLD2_RH; copyto!(view(v(o), a-1:b+1), view(v(n), a-1:b+1)); end
copyto!(view(v(:ero), rng), view(v(:ern), rng)); copyto!(view(v(:fro), rng), view(v(:frn), rng))
for k in a:b+1; v(:xm0)[k] = xmn[k]; v(:r0)[k] = rn[k]; v(:u0)[k] = un[k]; end
TF.save_state(ARGS[2]); println("saved ", ARGS[2], " (lgrid = 1, tau = ", gridtau, ", lady = ", lady, ", weights = ", wts, ")")
