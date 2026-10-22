# S7a test c: Iliev et al. (2009, arXiv 0905.2920) Test 5, the expansion of an H II region in an initially uniform gas, with TITAN's spherical Lagrangian
# hydrodynamics and the photon-conserving transfer + H chemistry/energy of RadTransfer/HII (s7a_hii.jl).
#   TITAN_NOLIB=1 julia --project=env s7a_hydro.jl            env: NZ (default 300), TEND [Myr] (500), DTMAX [Myr] (0.5), AVZONE (1.5), TAG, DL (bin width 0.05)
# n_H = 1e-3 cm⁻³, T = 100 K, x = 0, Ṅ = 5e48 s⁻¹ (blackbody 1e5 K), box 15 kpc, pure H (case B, Cen cooling without Compton), source at r = 0.
# Material: pure H, p = (1 + x) ρ kT/m_H, e = (3/2) p/ρ (thermal; the ionisation energy is not part of the hydro energy).
# Output: lib/popiii/s7a_hydro_<TAG>.csv  t[Myr], r_I [kpc] (x = 0.5 outermost), r_I/r_S0, T_hii, r_shock [kpc] (density maximum), n_max, and the log.
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl"))
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl")); include(joinpath(@__DIR__, "AtomicHHe.jl"))
include(joinpath(@__DIR__, "RadTransfer.jl")); include(joinpath(@__DIR__, "s7a_hii.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore; const TM = TitanMaterials
const RT = RadTransfer
include(joinpath(@__DIR__, "Trbdf2.jl"))
get(ENV, "INTEGRATOR", "trbdf2") == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)

const kB = 1.380649e-16; const mH = 1.6735575e-24; const Myr = 3.15576e13; const kpc = 3.0857e21
mutable struct HIIGas <: TM.Material
    x::Vector{Float64}
end
@inline function TX.material_gas(m::HIIGas, ρ, T)
    x = m.x[TX.ZONE[]]; p = (1 + x) * ρ * kB * T / mH
    (p, 1.5 * p / ρ, sqrt(5.0 / 3.0 * p / ρ))
end
function TM.fill_eos!(m::HIIGas)
    dn, tn, pgn, egn, asn = TM.v(:dn), TM.v(:tn), TM.v(:pgn), TM.v(:egn), TM.v(:asn)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        pg, eg, as_ = TX.material_gas(m, dn[k], tn[k]); pgn[k] = pg; egn[k] = eg; asn[k] = as_
    end
    nothing
end
TM.fill_opac!(::HIIGas) = nothing

nz = parse(Int, get(ENV, "NZ", "300")); tend = parse(Float64, get(ENV, "TEND", "500")) * Myr; dtmax = parse(Float64, get(ENV, "DTMAX", "0.5")) * Myr
tag = get(ENV, "TAG", "n$(nz)"); dl = parse(Float64, get(ENV, "DL", "0.05"))
nH0 = 1e-3; L = 15kpc; Ndot = 5e48; T0 = 100.0
E = 13.598 .* exp.(range(0, log(100 / 13.598), step = dl)); E[1] = 13.5985
g = RT.edge_grid(E); src = RT.blackbody_photons(g, 1e5, Ndot)
σ = [RT.binavg(ν -> ChemRates.σ_verner(ν, 13.598, ChemRates.VERNER_HI), g, i) for i in 1:RT.nfreq(g)]
HII.AB_FIXED[] = 0.0
HII.init_ek([(ChemNet.χH * (1 - 1 / n^2), 2.0n^2, ChemNet.χH / n^2, ν -> AtomicHHe.σ_hyd(n, 1, ν)) for n in 2:15])
trec = 1 / (2.59e-13 * nH0); rS0 = cbrt(3Ndot / (4π * 2.59e-13 * nH0^2))
rI = collect(range(0, L, length = nz + 1)); rho = fill(nH0 * mH, nz); Tz = fill(T0, nz); u0 = zeros(nz + 1)
mat = HIIGas(fill(1e-8, 2000))
const AVZ = parse(Float64, get(ENV, "AVZONE", "1.5"))
av = (ql0 = 0.0, ql1 = 0.0, cq1 = 0.0, cq2 = 1.0)
PS.build!(; rI, rho, T = Tz, u = u0, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = nH0 * kB * T0 * (1 + 1e-8), selfgrav = false, material = mat, av)
TX.AV_ZONE[] = AVZ
TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); TX.PEXT[] = TF.var(:pgn)[b]
TF.setscalar!(:thet, parse(Float64, get(ENV, "THETA", "0.55")))
TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.05"))
TF.setscalar!(:dtime, 1e-3 * Myr); TF.setscalar!(:tmax, tend)
io = open(joinpath(@__DIR__, "..", "lib", "popiii", "s7a_hydro_$(tag).csv"), "w")
println(io, "t_Myr,r_I_kpc,r_I_over_rS0,T_hii,r_shock_kpc,n_max,xmean_in,step")
st = HII.State(fill(1e-8, nz), fill(T0, nz))
outs = [1.0, 3.0, 10.0, 30.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0] .* Myr; iout = Ref(1)
t0wall = time(); Nstep = Ref(0)
function chem_step(k)
    rn = TF.var(:rn); dn = TF.var(:dn); tn = TF.var(:tn); t = TF.getscalar(:timen); dt = t - TF.getscalar(:timeo)
    r = rn[a:b+1]; nH = dn[a:b] ./ mH
    st.T .= tn[a:b]
    HII.sweep!(st, g, σ, src, r, nH, dt; energy = true, collisions = true)
    for j in 1:nz; mat.x[a+j-1] = st.x[j]; tn[a+j-1] = st.T[j]; end
    for kk in a-2:a-1; kk >= 1 && (mat.x[kk] = mat.x[a]); end
    for kk in b+1:b+2; mat.x[kk] = mat.x[b]; end
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    TM.fill_eos!(mat)
    # time-step caps: ≤ DTMAX and ≤ 10 % of the elapsed time (R-type phase)
    TF.setscalar!(:dtime, min(TF.getscalar(:dtime), dtmax, max(0.1 * t, 1e-3 * Myr)))
    Nstep[] += 1
    if iout[] <= length(outs) && t >= outs[iout[]]
        xm = [st.x[j] for j in 1:nz]; kf = findlast(j -> xm[j] >= 0.5, 1:nz); rr = [0.5 * (r[j] + r[j+1]) for j in 1:nz]
        rIf = kf === nothing ? 0.0 : kf == nz ? rr[nz] : rr[kf] + (xm[kf] - 0.5) / (xm[kf] - xm[kf+1]) * (rr[kf+1] - rr[kf])
        # mass-weighted temperature of the ionised gas
        w = [dn[a+j-1] * (r[j+1]^3 - r[j]^3) * (xm[j] > 0.5) for j in 1:nz]; Th = sum(w .* st.T) / max(sum(w), 1e-300)
        ks = argmax(nH); println(io, @sprintf("%.3f,%.5f,%.5f,%.1f,%.5f,%.4e,%.4f,%d", t / Myr, rIf / kpc, rIf / rS0, Th, rr[ks] / kpc, nH[ks], sum(xm .* diff(r .^ 3)) / (r[end]^3 - r[1]^3), Nstep[]))
        flush(io)
        @printf("t = %7.2f Myr  step %5d  r_I = %.4f kpc  (%.3f rS0)  T_hii = %.0f  r(nmax) = %.3f kpc  n_max = %.3e  dt = %.3f Myr  wall %.0f s\n", t / Myr, Nstep[], rIf / kpc, rIf / rS0, Th, rr[ks] / kpc, nH[ks], dt / Myr, time() - t0wall)
        iout[] += 1
    end
end
try
    TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = chem_step)
catch e
    println("stopped: ", sprint(showerror, e)[1:min(end, 400)])
end
close(io)
@printf("rS0 = %.4f kpc, trec = %.2f Myr, steps %d, wall %.0f s\n", rS0 / kpc, trec / Myr, Nstep[], time() - t0wall)
