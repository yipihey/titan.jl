# Extend the outer (inflow) boundary of a saved quasi-star infall state to a larger radius.
# STATUS [tested, does not work]: the first Newton step after the graft fails even for NADD = 1 (huge corrections at the
# photosphere and above the shock; only the flux-equation residual differs from the ungrafted state, cause not found).
# The wide-domain run W30 was therefore set up from scratch (QS_ROUT=30, see PORTING_NOTES P7.10). Kept for reference.
#   TITAN_NOLIB=1 Z0_TABLE=ext julia --project=env extend_grid.jl IN.tstate.gz OUT.tstate.gz RMAX_over_R NADD [XSCALE_over_R]
# The NADD new zones continue the grid geometrically from the old outer boundary (first width = the last zone's width)
# to RMAX·R and are filled with the steady free fall at the inflow rate through the old boundary:
#   ρ = Ṁ/(4πr²v_ff), u = −v_ff, T from dilute radiative equilibrium continued from the old boundary,
#   T⁴ = T_out⁴·W(r)/W(r_out) with W the geometric dilution factor of the photosphere R, E = aT⁴, F = F_out (r_out/r)².
# The external (ghost) inflow state moves to the new boundary. Masses, volumes, old = new and the derived quantities are
# rebuilt as in QuasiStarSetup.setup!. XSCALE (default: unchanged) rescales the radius term of the adaptive-grid
# resolution, which otherwise would pull zones into the long new region.
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
using Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const PS = ProblemSetup
const TC = TitanCore; const Z0 = Z0Physics
TF.load_state!(ARGS[1]); PS.install_julia!(BlockSolver)
TX.CONVECTION[] = true; TX.ALPHA_MLT[] = 0.8; TM.use!(TM.Z0Gas()); TC.EDDFAC[] = TC.eddfac!
rmaxR = parse(Float64, ARGS[3]); nadd = parse(Int, ARGS[4])
v(s) = TF.var(s)::Vector{Float64}; s!(k, x) = (a = TF.var(k); a[] = eltype(a)(x))
a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
length(v(:rn)) >= b + nadd + 3 || error("arrays too short: set TITAN_MGR ≥ $(b + nadd + 3)")
# Every zone-indexed array (first dimension = the zone index, like rn) gets the old last zone's values in the new zones
# and ghosts first, so auxiliary quantities (artificial viscosity, Eddington factors, flux-limiter and bookkeeping
# arrays) start from sensible values; the physical variables are overwritten below.
let nz = length(TF.var(:rn)), src = b
    for (n, arr) in TF.JSTATE
        (ndims(arr) >= 1 && size(arr, 1) == nz && eltype(arr) === Float64) || continue
        n in (:rn,) && continue
        for k in b+1:b+nadd+2
            selectdim(arr, 1, k) .= selectdim(arr, 1, src)
        end
    end
end
R = 8.4352e15; M = 1e6 * Z0.Msun          # photospheric radius of the set-up and total mass (as Spectrum/qstar scripts)
rn, dn, tn, un, ern, frn = v(:rn), v(:dn), v(:tn), v(:un), v(:ern), v(:frn)
r0 = rn[b+1]; w0 = rn[b+1] - rn[b]
Mdot = -4π * r0^2 * TF.getscalar(:phir0)          # phir0 = −ρ_x v_ff (outward positive)
T0 = tn[b]; F0 = frn[b+1]
W(r) = 0.5 * (1 - sqrt(max(1 - (R / r)^2, 0.0)))
# geometric ratio q with w0 (q + q² + … + q^nadd) = RMAX·R − r0
L = rmaxR * R - r0
q = nadd == 0 ? 1.0 : let lo = 1.0 + 1e-9, hi = 2.0
    for _ in 1:200
        qq = 0.5(lo + hi); s = w0 * qq * (qq^nadd - 1) / (qq - 1); s > L ? (hi = qq) : (lo = qq)
    end
    0.5(lo + hi)
end
@printf("old: %d zones, r_out = %.5f R, last width %.3e R, Ṁ = %.3f Msun/yr; new: %d zones to %.1f R, ratio %.4f\n",
        b - a + 1, r0 / R, w0 / R, Mdot * Z0.yr / Z0.Msun, nadd, rmaxR, q)
vff(r) = sqrt(2Z0.G * M / r)
w = w0
for k in b+1:b+nadd
    global w *= q
    rn[k+1] = rn[k] + w
end
nadd > 0 && (rn[b+nadd+1] = rmaxR * R)
for k in b+1:b+nadd
    rc = ((rn[k]^3 + rn[k+1]^3) / 2)^(1 / 3)
    dn[k] = Mdot / (4π * rc^2 * vff(rc)); tn[k] = T0 * (W(rc) / W(0.5(rn[b] + rn[b+1])))^0.25
    ern[k] = Z0.arad * tn[k]^4
end
for k in b+1:b+nadd+1
    un[k] = -vff(rn[k]); frn[k] = F0 * (r0 / rn[k])^2
end
b += nadd; s!(:ngre, b)
# external inflow state (leobc = 2) at the new boundary
k = b + 1
rext = rn[k]; ρx = Mdot / (4π * rext^2 * vff(rext)); Tx = tn[b]
dn[k] = ρx; dn[k+1] = ρx; tn[k] = Tx; tn[k+1] = Tx; rn[k+1] = 2rn[k] - rn[k-1]; un[k] = -vff(rext); un[k+1] = -vff(rext)
s!(:ur, -vff(rext)); s!(:dr, ρx); s!(:tr, Tx); s!(:phir0, -ρx * vff(rext)); s!(:phir2, -ρx * vff(rext) * Z0.gas_e(ρx, Tx))
ern[k] = ern[k-1]; ern[k+1] = ern[k]; frn[k+1] = frn[k]
# masses, volumes, metric
xmup1 = 3.0; dvol(k) = (rn[k+1]^3 - rn[k]^3) / xmup1
xmn = v(:xmn)
for k in a+1:b+2; xmn[k] = xmn[k-1] + dn[k-1] * dvol(k - 1); end
xmext = TF.getscalar(:xmext); xmass = xmn[b+2]; xmtot = xmass + xmext
s!(:xmass, xmass); s!(:xmtot, xmtot); TF.var(:yscl)[1] = xmtot
xmen = v(:xmen); for k in a-1:b+2; xmen[k] = xmtot - xmn[k]; end
s!(:delmr, dn[b] * dvol(b)); s!(:delrr, rn[b+2] - rn[b+1])
for k in a-1:b+2; v(:rmun)[k] = rn[k]^2; v(:rmup1n)[k] = rn[k]^3; v(:rmum1n)[k] = rn[k]; end
for k in a-1:b+1; v(:dvoln)[k] = dvol(k); end
length(ARGS) >= 5 && (s!(:xscale, parse(Float64, ARGS[5]) * R); @printf("xscale = %.1f R\n", parse(Float64, ARGS[5])))
# bookkeeping of the new zones, old = new, derived quantities
newr = b-nadd+1:b+2
for n in (:dsn, :egsn, :etn, :etsn, :ersn, :egtn, :egtsn, :egrsn, :egrtn, :egrtsn, :usn, :frsn, :qx, :xnto, :xnt, :tescr, :tescr0)
    view(v(n), newr) .= 0.0
end
rng = a-1:b+2
for (o, n) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(n), rng)); end
copyto!(view(v(:ero), rng), view(ern, rng)); copyto!(view(v(:fro), rng), view(frn, rng))
for k in newr; v(:frnom)[k] = abs(frn[k]) + 1e-30; v(:unom)[k] = abs(un[k]) + 1e-30; end
TD.EOS[](); TD.OPAC[]()
for (o, n) in TC.PEQ_OLD2_RH; copyto!(view(v(o), a-1:b+1), view(v(n), a-1:b+1)); end
TC._ghost_and_derived!(:peq, () -> TD.EOS[](), () -> TD.OPAC[]())
for (o, n) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(n), rng)); end
for (o, n) in TC.PEQ_OLD2_RH; copyto!(view(v(o), a-1:b+1), view(v(n), a-1:b+1)); end
copyto!(view(v(:avchio), a:b+1), view(v(:avchin), a:b+1))
for k in a:b+1; v(:xm0)[k] = xmn[k]; v(:r0)[k] = rn[k]; v(:u0)[k] = un[k]; end
s!(:tfac, 1.0); s!(:dtime, min(TF.getscalar(:dtime), 1e3))
TF.save_state(ARGS[2])
@printf("saved %s: %d zones, r_out = %.2f R, ρ_out = %.3e, T_out = %.0f K, u_out = %.0f km/s, xmtot %.4e\n",
        ARGS[2], b - a + 1, rn[b+1] / R, ρx, Tx, -vff(rext) / 1e5, xmtot)
