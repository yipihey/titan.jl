"""
    TitanResidFast

Allocation-free version of TitanResid (same equations, verified against it and against
Fortran matgen): typed context, per-element-type cached workspace, one stencil-coloured
ForwardDiff sweep with 5 colours (35 partials for rad-hydro), side outputs taken from the
Dual values (no second Float64 pass).

Colouring note: the residual reaches ±3 blocks in exactly one place (radmom: ρ(k−3)
through the F/ρ interface interpolation, ~1e-8 of the row scale).  With 5 colours that
dependence aliases into the +2 entry of the same row (error ~1e-8 of the row scale);
TITAN drops it altogether.  `ncol = 7` (slower) separates it cleanly.
"""
module TitanResidFast

using ForwardDiff
using Printf
import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

val(x) = ForwardDiff.value(x)
gi(s) = Int(TD.g(s)); gf(s) = Float64(TD.g(s))

struct Par
    dt::Float64; θ::Float64; cadv::Float64; eps::Float64; sigd::Float64; sige::Float64
    phil0::Float64; phil2::Float64; phir0::Float64; phir2::Float64; cgrav::Float64; g::Float64
    xmtot::Float64; alph::Float64; tau::Float64; xscale::Float64; ql0::Float64; ql1::Float64
    cq1::Float64; cq2::Float64; cqvis::Float64; gam::Float64; gmu::Float64; ck::Float64; cm0::Float64
    csige::Float64; chif0::Float64; ratio::Float64; cc::Float64; cpi::Float64; csigr::Float64
    geddr::Float64; ximr::Float64; delml::Float64; delmr::Float64; delrl::Float64; delrr::Float64
    dl::Float64; dr::Float64
    ibet::Int; ladx::Int; lgrid::Int; leibc::Int; leobc::Int; llibc::Int; llobc::Int; lribc::Int
    lrobc::Int; lopac::Int
end

const NV = 7   # variable codes: 1 rn 2 xmen 3 dn 4 un 5 tn 6 ern 7 frn
const VSYM = (:rn, :xmen, :dn, :un, :tn, :ern, :frn)

struct Ctx
    ngrs::Int; ngre::Int; neqn::Int; mu::Int; lrad::Bool; k0::Int; nb::Int; mgr::Int
    ix::NTuple{7,Int}                 # equation/variable index of each code (0 = absent)
    code::Vector{Int}                 # code[j] for Newton variable j
    base::Matrix{Float64}             # base[k, code]
    nomu::Vector{Float64}; nomf::Vector{Float64}
    ro::Vector{Float64}; xmeo::Vector{Float64}; do_::Vector{Float64}; uo::Vector{Float64}
    ego::Vector{Float64}; ero::Vector{Float64}; fro::Vector{Float64}; egro::Vector{Float64}
    pgo::Vector{Float64}; aso::Vector{Float64}   # (no stored old limiter slopes: `oldslope!` recomputes them from the old cell values)
    chifo::Vector{Float64}; xkeo::Vector{Float64}; xkpo::Vector{Float64}; plfo::Vector{Float64}
    dvolo::Vector{Float64}; xmo::Vector{Float64}; fedd::Vector{Float64}; bri::Vector{Float64}
    lady::Vector{Int}; wt::Vector{Float64}; yscl::Vector{Float64}
    p::Par
    etotal::Bool                      # conservative total-energy row instead of gasnrg
    useint::Vector{Bool}              # dual-energy switch: zone k keeps the internal-energy row
    conv::Bool                        # MLT convective flux in the energy equation (Phase 7)
    alpha::Float64                    # mixing-length parameter
    cv::NTuple{4,Vector{Float64}}     # interface ∇_ad, c_P, δ, κ_R from the current iterate (held during the AD sweep)
end

const CTXBUF = Ref{Any}(nothing)
const BUFS = Dict{Symbol,Any}()     # reused context buffers (allocation-free context())
"Defaults for the energy formulation (Phase 4 item 2); ETOTAL[] = true selects the total-energy row."
const ETOTAL = Ref(false)
const ETASW = Ref(0.0)
"Added to e in the dual-energy test, for a material whose e has a negative zero point (Population III: −1e13 erg/g), where the ratio otherwise has no meaning and can be singular."
const ETASW_EOFF = Ref(0.0)
"""
ETASW_TH[] (`ETASWTH`, default 0 = off): a zone also keeps the internal-energy row where its thermal energy 3p/2ρ is less than this fraction of
3p/2ρ + E_r/ρ + K/m. There the total-energy row gives the gas energy as a small difference, and it is not the first law of the gas: total row minus
kinetic energy minus radiation-energy row leaves −u [∂(fE)/∂r + (3f − 1)E/r + χρF/c] on the gas, which vanishes only if the radiation momentum equation
holds in its static form. In transparent zones with a gradient of E_r it does not (the row is `FLUXTIE`'s there), and inflow against the gradient
cools the gas: 0.12 K/yr at 2.6 pc in the 146-zone run, to 17 K and then 11.7 K, where the run stopped twice (11 Nov 2026; docs/POPIII_RUNS.md).
"""
const ETASW_TH = Ref(0.0)
"""
ETASW_TIE[] (`ETASWTIE`, default 0 = off): a zone also keeps the internal-energy row where the FLUXTIE weight at its inner
interface is at least this value (0 < ETASWTIE <= 1). Where the radiation momentum row is the closure F = h c E instead of the
moment equation, total-energy row − kinetic energy − radiation-energy row leaves −u × (moment-equation residual) on the gas;
this selects the first law there by cause instead of by the size of the thermal energy (ETASWTH).
"""
const ETASW_TIE = Ref(0.0)
"""
GASROW[] (lrad = 1 only): energy row of the gas ALONE,
  Δ(m e_g)/dt + advection of e_g + p_g ∇·u + Q_visc + MLT flux = −(c κ_E E − 4π κ_P B) ρ dV,
instead of TITAN's row for gas + radiation, Δ(m(e_g + E/ρ))/dt + … (from which the gas energy is the difference
with the radiation row). The two are algebraically equal (the exchange term is the radiation row's, with
opposite sign, so the pair still conserves total energy), but differ in the discrete advection: TITAN advects
e_g + E/ρ and E/ρ with separate van Leer limiters. Where E/ρ ≫ e_g (thin infalling gas above a radiating star:
1e5–1e6) the difference of the two advected fluxes is dominated by limiter and round-off differences and the gas
temperature is garbage (zones cooling to a few hundred K inside a 6000 K radiation field). Default false (TITAN).
"""
const GASROW = Ref(false)
"""
COMPTON[] (lrad = 1, material with `material_gasne`): Compton energy exchange between gas and radiation,
  Λ_C = 4 σ_T n_e c E k_B (T − T_r)/(m_e c²),  T_r = (E/a)^{1/4}  (non-relativistic Kompaneets, grey radiation),
lost by the gas and gained by the radiation (in the radiation row and, with GASROW, in the gas row; TITAN's
combined row needs no term). θ-weighted between old and new state, like the other exchange terms. Default false (TITAN has none).
Matters in hot thin gas: behind a 1e8 K accretion shock in a 10 erg cm⁻³ radiation field t_C ≈ 1.5e6 s, far
shorter than the free-free cooling time.
"""
const COMPTON = Ref(false)
"""
INNER_SINK[] (lrad = 1): the inner wall becomes porous to mass. A mass flux `phil0` (< 0: out of the grid, onto the
central point mass, which grows through m(r) = m_tot − m_e) carries the gas enthalpy in `phil2` (set by the caller,
e.g. phil0·(e_g + p_g/ρ) of zone ngrs) and, with this switch, the radiation enthalpy (1 + f)E/ρ·phil0 in the radiation
row (and in TITAN's combined row). The wall itself does not move (u = 0 there). Default false.
"""
const INNER_SINK = Ref(false)
"Sign of the inner-boundary mass flux in the mass row (see the MASS rows of residual!). Default true (fixed)."
const FIX_PHIL0 = Ref(true)
"""
REGULAR_CENTRE[] (Population III, M1): interface ngrs is the centre, r(ngrs) = 0, u(ngrs) = 0 (inner wall, leibc = 1).
The radius row of ngrs becomes the fixed-value row R = x (the multiplicative correction r·eˣ has ∂R/∂x = 0 at r = 0),
and u/r at the centre takes its regular limit u(ngrs+1)/r(ngrs+1) (u ∝ r near the centre) in the viscosity and the
(1 − 3f) energy terms. Default false.
"""
const REGULAR_CENTRE = Ref(false)
"""
XMINT[] (Lagrangian grid, lgrid = 3): interior mass per steradian of each interface, fixed in time, used for gravity
instead of x_mtot − x_me. The exterior mass is known only to the Newton tolerance (δ ln x_me ≤ conv) times x_mtot,
far more than the central shells' masses; on a Lagrangian grid the interior mass of an interface never changes.
Empty (default) = TITAN's x_mtot − x_me.
"""
const XMINT = Ref(Float64[])
"""
PEXT[] (lgrid = 3, llobc = 1): free outer boundary. The outermost interface obeys the momentum equation of its half
zone with a constant external pressure PEXT[] outside (RHFC02 App. A2.3) instead of the piston u = uextr. NaN
(default) = TITAN's piston.
"""
const PEXT = Ref(NaN)
"""
GRIDVEL[] (lgrid = 2 only): `nothing` (default) = every interface fixed. Otherwise a vector of interface velocities w [cm/s] (driver
index, zero for an interface that stays): the grid row is r_new − r_old − w Δt with the Δt of the stage, so an interface follows the
same straight path r_old + w t through the TR, the BDF2 and the backward-Euler stage of a step. The driver sets w between steps
(GRIDEQ=0 in popiii_m2c.jl: the mesh relaxes to a prescribed target), so the row adds nothing to the Jacobian. All
else is TITAN's moving-grid form (u_rel = u − ṙ, fluxes through the moving faces, lab-frame gravity work).
"""
const GRIDVEL = Ref{Union{Nothing,Vector{Float64}}}(nothing)
"""
GRIDEQ[] (lgrid = 2 only; `nothing` = off): TITAN's implicit grid equation over the whole mesh (GRIDEQ in popiii_m2c.jl). A
NamedTuple (; k1, k2, w = (wρ, wP, wT, wE, wm), alph, gmax, tauv, sx): the interfaces k1 (the first above the regular centre) and
k2 (the outer edge) stay (with GEQ_LAGFIRST, ICMLOG, k1 is Lagrangian instead: it holds its enclosed mass, not its radius); every interface in between obeys n̂_j = n̂_{j−1} with n̂_j = x̃_j / sqrt(1 + G_j²), x̃_j = ν_j −
α(1+α) Δ²ν (ν_j = r̄_j/Δr_j, the logarithmic baseline; α = alph), time-smoothed per zone, x̃ → (1 + γ_j) x̃ − γ_j x̃_old with
γ_j = tauv[j]/Δt (tauv: a multiple of the zone's signal-crossing time, frozen for the step), and G_j² = ν_j² S_j capped smoothly at
gmax², S_j = Σ_q w_q (½ Δ ln q)² over q = ρ, P, T, E_r and the interior mass m (centred differences of the zone values; one-sided
at the ends) plus sx[j] (terms frozen for the step: the capped ln τ). All terms are logarithmic, so the rule has no scale.
"""
const GRIDEQ = Ref{Any}(nothing)
"""
GEQ_LAGFIRST[] (lgrid = 2 with GRIDEQ or GRIDVEL; default false): the first interface above the regular centre, k1 = ngrs + 1, is Lagrangian: its grid row holds
its enclosed mass fixed, (x_m,e,new − x_m,e,old) r_old / (zone mass of k1 − 1, old), instead of its radius (the other interface rows are unchanged). The
row is a function of the iterate only (the scale factors are old-level constants; r_old/zone mass makes it a length like the radius rows). ICMLOG sets it.
"""
const GEQ_LAGFIRST = Ref(false)
"A mesh with lgrid = 2 whose interfaces move (GRIDVEL or GRIDEQ)."
@inline gridmoves() = GRIDVEL[] !== nothing || GRIDEQ[] !== nothing
"""
AV_ZONE[] (G1): viscous length proportional to the local zone width. When > 0, zone k uses the length
l_k = AV_ZONE[] × (r[k+1] − r[k]) (θ-centred radii, ForwardDiff-differentiable) instead of ql0 + ql1 r̄; ql0, ql1 are then
ignored (cq1, cq2, cqvis still apply). On grids where Δr/r varies (geometric mass grids) l ∝ r does not follow the zone
width; l ∝ Δr keeps the shock spread over the same number of zones. Mapping to the von Neumann–Richtmyer form
q ≈ c_q ρ (Δv)² (Ripamonti et al. 2002 use c_q = 4): with qv = div u and l = c Δr the quadratic term is
q = cqvis·cq2·l²·ρ·qv·(du/dr) ≈ cqvis·cq2·c²·ρ (Δr div u)², i.e. c_q ≈ cqvis·cq2·c² (cqvis = 4/3 in our states, so
c_q = 4 needs c² ≈ 3/cq2). This mapping is derived here, not taken from the source. Default 0 = the old formula.
"""
const AV_ZONE = Ref(0.0)
"RAYQ[] (M3; empty = off): per zone, the energy the ray (IonRay) takes from the star's light [erg cm⁻³ s⁻¹]. The gas receives it in chemistry; here it leaves the gray field, so it is not counted twice."
const RAYQ = Ref(Float64[])
"""
GASQ[] (M3a; empty = off): per zone, energy moved from the gray field to the gas inside the hydro step [erg cm⁻³ s⁻¹], lagged over the step
(Lyman-Werner photons absorbed by H₂, PopIIIMaterial.STARLIGHT). Unlike RAYQ the gas receives it here, in the implicit energy equation, where
its radiative cooling is: given to the gas in the chemistry step (no cooling there in RAD runs) it heated a zone near the photosphere from
7300 to 18500 K over one step. Rad row +GASQ, gas-only row −GASQ, total rows nothing (an internal transfer)."""
const GASQ = Ref(Float64[])
"""
CQB[], YQB[] (empty = off): per zone, the three-level atom's chemistry heat cq and cycle heat yq [erg cm⁻³ s⁻¹] of the END state of the last chemistry step (PopIIIMaterial.chem_core!),
for the zones whose hydro rows do not compute them from the iterate. A zone that is coupled in the step (`cplz`, the flag of the whole step) takes cq, yq of the Newton iterate
(`cpQ`, `cpY` of `cpl_g`); every other zone takes these booked values, so each step applies the heat of each zone exactly once, whatever the zone's status in the neighbouring steps.
Radiation row +cq+yq, gas rows −yq (cq is inside the gas row's cooling `cpL` for coupled zones; for the other zones the cq of the booking leaves only the radiation), total/combined rows +cq.
"""
const CQB = Ref(Float64[]); const YQB = Ref(Float64[])
@inline cqheat(k, W) = cplz(k) ? W.cpQ[k] : (isempty(CQB[]) ? zero(eltype(W.cpQ)) : CQB[][k])
@inline yqheat(k, W) = cplz(k) ? W.cpY[k] : (isempty(YQB[]) ? zero(eltype(W.cpY)) : YQB[][k])
"""
CHEMB[] (S6, ENV BALSPLIT=1; empty = off): balanced operator splitting (Speth, Green, MacNamara & Strang 2013) of the gas energy between this
hydro step and the network chemistry. Per zone, c̄ [erg g⁻¹ s⁻¹] = the chemistry's net specific heating rate of the last step; the gas
energy row gets +ρ c̄ (here, all variants) and the chemistry integrates its heating − c̄ (PopIIIMaterial.chem!, network zones), so the sum
is the unsplit equations and in a steady state neither operator moves T. Without it, at 107.5 M_sun, the first network zone outside the
photosphere was heated 6227 → 8097 K in every hydro step (lagged three-level cycle heat GASQ against its radiative cooling) and cooled
back by its recombination losses in the chemistry: the 30 % excursion held Δt at the STOL limit (~2e3 s) for tens of thousands of steps.
c̄ ← c̄ + Δe_chem/Δt after each chemistry step in zones that stayed on the network (zero elsewhere). CHEMBON[k]: zone k was balanced.
"""
const CHEMB = Ref(Float64[]); const CHEMBON = Ref(Bool[])

"""
COUPLED (popiii_m2c COUPLED=1; default false): the slow species of the reduced network (ChemNet CHEMCK: x_H⁺, x_H₂ per H nucleus) are
unknowns 8 and 9 of the Newton system (neqn = 9, meqn ≥ 9), in logarithmic variables (x = base·e^δ, positive by construction), in the
network zones below the hot branch at the start of the step (`CPLZ`); other zones carry them frozen (row δ = 0). Per coupled zone
  species row  (ρ x ΔV)ⁿ⁺¹ − (ρ x ΔV)ⁿ)/Δt + ṁ-fluxes of the face values (advect!, the scheme and mass fluxes of the energy row)
               − [θ ρⁿ⁺¹ΔVⁿ⁺¹ S(ρ, T, x) + (1 − θ) ρⁿΔVⁿ Sⁿ] = 0,       S = dx/dt of the CK-expanded composition (ChemNet.ck_eval)
  energy row   + [θ Λ(ρ, T, x) + (1 − θ) Λⁿ] ΔV,                       Λ the zone's chemistry cooling (lines, HD, three-level heat, …)
and the zone's EOS is the network EOS of the CK expansion of the iterate, so the latent heat of ionisation and H₂ formation is in the
Jacobian. The material provides (p, e, S₁, S₂, Λ) of a zone with their derivatives with respect to (ρ, T, x₁, x₂) through `CPL_ZONE`
(a 4-wide local ForwardDiff Jacobian lifted onto the residual's duals). Old-state arrays: `XSO` (and its limiter slopes `XSSO`), the old
rates `XSRCO` and cooling `XLAMO` (θ < 1 stages only). Trbdf2 combines `XSO` on ρ x ΔV for its BDF2 stage.
"""
const COUPLED = Ref(false)
const XSN = Ref([Float64[], Float64[]]); const XSO = Ref([Float64[], Float64[]]); const XSSO = Ref([Float64[], Float64[]])
const XSRCO = Ref([Float64[], Float64[]]); const XLAMO = Ref(Float64[]); const CPLZ = Ref(Bool[])
const CPLW = Ref{Any}(nothing)                  # (xs1, xs2) of the residual's iterate while the material is evaluated
const CPL_ZONE = Ref{Any}(nothing)              # (k, ρ, T, x1, x2) -> (p, e, S1, S2, Λ), set by the material
"""
SLIMP (popiii_m2c `SLIMP=1`, with COUPLED): the star light of a coupled zone at the Newton iterate. `SL_SWEEP[](S, ...)` fills, from the
iterate's arrays, the state S = (q, t, fs, g, dp) of every coupled zone: q = W/W_start (dilution), t = T_c/T_c,start (colour temperature), fs
(the Lyman-Werner factor f_sh × its on-weight) and g (the angle-mean velocity gradient of Lyα escape); `CPL_ZONEQ[](k, ρ, T, x1, x2, q, t, fs, g)`
gives the zone's terms at that star light (PopIIIMaterial.cpl_zone_q) and `CPLSL[]` carries S (Duals) to the NLOIT opacity during the material sweep.
"""
const SL_SWEEP = Ref{Any}(nothing)
const SLIMP = Ref(false); const CPL_ZONEQ = Ref{Any}(nothing); const CPLSL = Ref{Any}(nothing)
const SLIMPDBG = Ref(parse(Int, get(ENV, "SLIMPDBG", "0")))   # (diagnostic) report the first n non-finite SLIMP values in the residual
const SLQBUF = Dict{Any,Any}()
const CPL_LOCAL = Ref{Any}(nothing)             # (k, dt) -> (x1, x2): the zone's species rows solved alone (initial guess), set by the material
"Weight floor of the species corrections in maxdel's convergence measure, and the clamp of a species correction per Newton iteration."
const CPL_XFLOOR = Ref(1e-8)
"""
Time weight of the chemistry source and cooling in the coupled rows (default 1: backward Euler in every stage, also in TR-BDF2's
trapezoidal stage; the residual's θ otherwise). A species relaxing in seconds within a step of 1e3–1e6 s has, under the trapezoidal
rule, x_γ ≈ −x_n: unreachable in the logarithmic species variables, so the TR stage failed and the step fell back to backward Euler.
"""
const CPL_THETA = Ref(parse(Float64, get(ENV, "CPLTHETA", "1")))
@inline cplz(k) = COUPLED[] && k <= length(CPLZ[]) && CPLZ[][k]
"Size the species arrays to mgr (keeps values where present)."
function cpl_alloc!(mgr)
    for R_ in (XSN, XSO, XSSO, XSRCO), i in 1:2
        length(R_[][i]) == mgr || (R_[][i] = zeros(mgr))
    end
    length(XLAMO[]) == mgr || (XLAMO[] = zeros(mgr)); length(CPLZ[]) == mgr || (CPLZ[] = falses(mgr) |> collect)
    nothing
end
"""
Newton's initial species of a coupled zone: its two species rows solved alone at the zone's (ρ, T), with the neighbours' species and
the mass fluxes of the last step frozen and first-order upwind faces (`CPL_LOCAL`, a 2×2 Newton in ln x). Without it the first Newton
corrections of ln x reached −25 … +440 in the pre-shock zones (the stiff chemistry's solution for the step lies far from the start
composition), and the ±1 clamp needed dozens of iterations.
"""
function cpl_local_all!(dt)
    CPL_LOCAL[] === nothing && return nothing
    ngrs = gi(:ngrs); ngre = gi(:ngre)
    # (PERF512) the zones are independent (every one sees the old neighbours; the results are applied together): on MATTHREADS tasks
    # (zone_sweep!), bit-identical; ~10 % of the hydro step of the 512-zone mesh when serial
    res = Vector{Union{Nothing,Tuple{Float64,Float64}}}(nothing, ngre); z0 = ZONE0[]
    try zone_sweep!(ngrs:ngre, cplz) do k
        y = CPL_LOCAL[](k, dt); y === nothing || (res[k] = (Float64(y[1]), Float64(y[2])))
    end
    finally ZONE0[] = z0 end                                             # (the serial loop did not touch ZONE[])
    for k in ngrs:ngre                                                   # (applied together: every zone saw the same neighbours)
        y = res[k]; y === nothing || (XSN[][1][k] = y[1]; XSN[][2][k] = y[2])
    end
    nothing
end
"Start of a step (TitanCore.peq!): old species := iterate, and their limiter slopes; then the initial guess (cpl_local_all!)."
function cpl_newstep!()
    ngrs = gi(:ngrs); ngre = gi(:ngre); cadv = gf(:cadv); eps = gf(:epsadv)
    for i in 1:2
        copyto!(XSO[][i], XSN[][i]); q = XSO[][i]; qs = XSSO[][i]; fill!(qs, 0.0)
        @inbounds for k in ngrs:ngre+1
            a = q[k] - q[k-1]; b = q[k+1] - q[k]; sdq = a + b; den = sdq + (sdq >= 0 ? eps : -eps)
            qs[k] = a * b >= 0 ? cadv * a * b / den : 0.0
        end
    end
    empty!(CPL_BAD)                                  # (only corrections of this step's Newton iterations count)
    cpl_local_all!(gf(:dtime))
    nothing
end
"""
Coupled zones whose raw species correction exceeded CPL_BADDELTA (100 in ln x) in a Newton iteration of a failed step: on the reset they are
handed back to the split chemistry for the rest of the step (CPLZ false; with ECHEMSPLIT treated as network, not table, zones). A zone next to
the shock whose x_H⁺ had fallen to 1e-9 asked for δ ln x = −1e76 in every retry, the step never converged and Δt collapsed (t_cpl4, step 123110).
"""
const CPL_BAD = Set{Int}(); const CPL_BADDELTA = Ref(parse(Float64, get(ENV, "CPLBADDELTA", "100"))); const CPL_NDECOUPLED = Ref(0)
function cpl_reset!()
    for i in 1:2; copyto!(XSN[][i], XSO[][i]); end
    if !isempty(CPL_BAD)
        a = gi(:ngrs)
        for k in CPL_BAD
            cplz(k) || continue
            CPLZ[][k] = false; CPL_NDECOUPLED[] += 1
            (CPL_NDECOUPLED[] <= 30 || CPL_NDECOUPLED[] % 100 == 0) && println("COUPLED: zone $(k - a + 1) decoupled for this step (species correction > $(CPL_BADDELTA[]) in ln x), $(CPL_NDECOUPLED[]) so far; chemistry/cycle heat now from the booking of the last chemistry: cq $(isempty(CQB[]) ? 0.0 : CQB[][k]), yq $(isempty(YQB[]) ? 0.0 : YQB[][k]) erg/cm3/s")
        end
        empty!(CPL_BAD)
    end
    cpl_local_all!(gf(:dtime)); nothing
end
"Apply the Newton correction (update!): x ← x·e^δ in the coupled zones."
function cpl_update!(rhs, a, b)
    @inbounds for k in a:b
        cplz(k) || continue
        XSN[][1][k] *= exp(rhs[8, k]); XSN[][2][k] *= exp(rhs[9, k])
    end
    nothing
end
"maxdel: the species' weighted maximum correction (convergence measure; with it TITAN's dtol damping of the whole correction limits a species step like every other)."
const CPLDM = Ref((0.0, 0, 0))     # (diagnostic, NWTRACE) the species' weighted convergence measure of the last cpl_dmax!, with its equation and zone
function cpl_dmax!(rhs, a, b)
    d = 0.0; f = CPL_XFLOOR[]; kx = TF.var(:kx); dx = TF.var(:dx); wm_ = 0.0; wi_ = 0; wk_ = 0
    for i in 8:9; kx[i] = 0; dx[i] = 0.0; end
    @inbounds for k in a:b
        cplz(k) || continue
        for (i, q) in ((8, XSN[][1]), (9, XSN[][2]))
            abs(rhs[i, k]) > abs(dx[i]) && (kx[i] = k - a + 1; dx[i] = rhs[i, k])     # (raw correction, for the divergence message)
            isfinite(rhs[i, k]) || (rhs[i, k] = 0.0; d = Inf)                        # a non-finite correction: no change, not converged
            abs(rhs[i, k]) > CPL_BADDELTA[] && push!(CPL_BAD, k)
            # weight by the larger of the iterate and the old abundance: a species the iterate drove to the floor (x_H⁺ 1e-40 next to the
            # shock, from 1e-3) is not converged; only species that are trace at both ends are left out
            xm_ = max(q[k], XSO[][i - 7][k]); d = max(d, abs(rhs[i, k]) * xm_ / (xm_ + f))
            w_ = abs(rhs[i, k]) * xm_ / (xm_ + f); w_ > wm_ && (wm_ = w_; wi_ = i; wk_ = k - a + 1)
        end
    end
    CPLDM[] = (wm_, wi_, wk_)
    d
end

"""
AV_TENSION[] (default true = TITAN): keep the viscous stress of a zone that shrinks in volume (div u < 0) while it is stretched
radially relative to the mean (dudr > 0, the traceless form). Then q < 0, a tension. With l = AV_ZONE·Δr and Δr ~ r (the stretched
zones outside an accretion shock) it can reach ~10× gravity and accelerate an infalling shell past the escape speed.
false: viscosity only under radial compression.
AV_TFRAC[], AV_TRAD[] (default 0, false = off; used with AV_TENSION[] = false): a reduced tension. AV_TFRAC = f in (0, 1): the stress of
the tension zones is kept, times f. AV_TRAD: it is kept only where the zone is compressed radially in the plain sense (∂u/∂r < 0: the
ramp of a shock whose zones are wide, Δr ~ r), not in smooth free fall (∂u/∂r > 0); times AV_TFRAC if that is set.
"""
const AV_TENSION = Ref(true)
const AV_TFRAC = Ref(0.0); const AV_TRAD = Ref(false)
"""
EXTMASS[] (Population III): static external mass distribution (e.g. the dark-matter halo of ABN02), a function
r -> M_ext(< r) per steradian, added to the gas interior mass in the gravity terms (evaluated at the θ-centred radius;
Float64 or Dual). `nothing` (default) = none.
"""
const EXTMASS = Ref{Any}(nothing)
"u/r at interface k, with the regular limit at the centre."
@inline uor(u, r, k, ngrs) = REGULAR_CENTRE[] && k == ngrs ? u[k+1] / r[k+1] : u[k] / r[k]
const MEC2 = 8.1871057769e-7            # m_e c² [erg]
"""
Dual-energy switch, evaluated on the OLD (start-of-step) state and frozen during the Newton
iteration: zone k keeps TITAN's internal-energy row when
η_k = m(e_g [+E_r/ρ]) / (m(e_g [+E_r/ρ]) + K_k) < ηsw  (K_k = m(u_k²+u_{k+1}²)/4).
"""
function dualswitch(ngrs, ngre, lrad, ηsw)
    mgr = TF.MGR; ui = get!(() -> fill(false, mgr), BUFS, :useint)::Vector{Bool}; fill!(ui, false)
    ηsw > 0 || return ui
    dO, dv, uo = v(:do), v(:dvolo), v(:uo); eo = lrad ? v(:egro) : v(:ego)
    for k in ngrs:ngre
        m = dO[k] * dv[k]; ein = m * max(eo[k] + ETASW_EOFF[], 0.0); K = 0.25m * (uo[k]^2 + uo[k+1]^2)
        ui[k] = ein / (ein + K) < ηsw
    end
    ηth = ETASW_TH[]
    if ηth > 0 && lrad
        pgo, ero = v(:pgo), v(:ero)
        for k in ngrs:ngre
            ui[k] && continue
            eth = 1.5 * pgo[k] / dO[k]
            ui[k] = eth < ηth * (eth + ero[k] / dO[k] + 0.25 * (uo[k]^2 + uo[k+1]^2))
        end
    end
    wtie = ETASW_TIE[]
    if wtie > 0 && lrad && TD.FLUXTIE[][2] > 0 && length(TD.TIEW) >= ngre
        W = TD.TIEW
        for k in ngrs:ngre
            ui[k] || (ui[k] = W[k] >= wtie)
        end
    end
    SDDSW[] === nothing || sd_dsw(ui, ngrs, dO, dv, uo, eo)
    ui
end
# stalldiag (SDDSW=lo:hi, default off): print the dual-energy row choice (ui = internal-energy row) of the zones lo:hi whenever it differs from the previous call
const SDTIEF = parse(Float64, get(ENV, "SDTIEF", "0"))   # (stalldiag, test) f of the SLABTIE correction replaced by this constant
const SDFIN = Ref(parse(Float64, get(ENV, "SDFIN", "1")))   # (stalldiag, test only) factor on the radiation-momentum inertia term (1/c^2) dF/dt of the F_r row
const SDFRZ = get(ENV, "SDFRZ", "0") == "1"; const SDFRZV = Float64[]
const SDW = let s = get(ENV, "SDW", ""); isempty(s) ? (0, 0, 0.0) : (f = split(s, ":"); (parse(Int, f[1]), parse(Int, f[2]), parse(Float64, f[3]))) end
const SDROW = get(ENV, "SDROW", "") != ""; const SDW0 = SDROW ? parse(Int, split(ENV["SDROW"], ":")[1]) : 0; const SDW1 = SDROW ? parse(Int, split(ENV["SDROW"], ":")[2]) : 0
const SDROWV = Dict{Int,NTuple{8,Float64}}()   # (SDROW=lo:hi) per interface: R_mom, B, w, dR_mom/dF (absorption part), dB/dF, dR_mom/dE_(k-1) (static f part), dB/dE_(k-1), dR_mom/dF (inertia part)
const SDHD = get(ENV, "SDHD", "0") == "1"; const SDHRAW = Dict{Int,NTuple{4,Float64}}()   # (SDHD=1) per interface: (raw h before max, h_FS, tau_i, a)
const SDDSW = Ref{Any}(let s = get(ENV, "SDDSW", ""); isempty(s) ? nothing : (parse(Int, split(s, ":")[1]):parse(Int, split(s, ":")[2])) end)
const SDDSWLAST = Ref("")
function sd_dsw(ui, ngrs, dO, dv, uo, eo)
    zs = SDDSW[]; bits = join((ui[ngrs+z-1] ? '1' : '0' for z in zs))
    if bits != SDDSWLAST[]
        etas = join((@sprintf("%d:%.3f/W%.2f", z, (m = dO[ngrs+z-1] * dv[ngrs+z-1]; ein = m * max(eo[ngrs+z-1], 0.0); ein / (ein + 0.25m * (uo[ngrs+z-1]^2 + uo[ngrs+z]^2))),
                      length(TD.TIEW) >= ngrs + z - 1 ? TD.TIEW[ngrs+z-1] : NaN) for z in zs), " ")
        println("SDDSW jstep ", Int(TF.getscalar(:jstep)), " thet ", TF.getscalar(:thet), " dtime ", TF.getscalar(:dtime), " ui[", first(zs), ":", last(zs), "] ", bits, "  eta/Wtie ", etas)
        SDDSWLAST[] = bits
    end
    nothing
end
v(s) = TF.var(s)
function context(; etotal::Bool = ETOTAL[], ηsw::Float64 = ETASW[], dt = nothing)
    ngrs = gi(:ngrs); ngre = gi(:ngre); neqn = gi(:neqn); lrad = gi(:lrad) == 1; mgr = TF.MGR
    lrad && compton_old!(ngrs, ngre)
    source_old!(ngrs, ngre)
    (gi(:ltur) == 0 && gi(:lhydr) == 1) || error("TitanResidFast: lhydr=1, ltur=0 only")
    MATERIAL[] === nothing && gi(:leos) != 2 && error("leos must be 2 (or install a material)")
    MATERIAL[] === nothing && lrad && !(gi(:lopac) in (3, 4)) && error("lopac must be 3 or 4 (or install a material)")
    lrad && gi(:ltran) != 1 && error("ltran must be 1")
    ix = Tuple(gi(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    code = get!(() -> zeros(Int, neqn), BUFS, Symbol(:code, neqn))::Vector{Int}; fill!(code, 0)
    for cde in 1:NV
        i = ix[cde]
        (1 <= i <= neqn) || continue
        (cde >= 6 && !lrad) && continue
        code[i] = cde
    end
    base = get!(() -> zeros(mgr, NV), BUFS, :base)::Matrix{Float64}
    for cde in 1:NV; copyto!(view(base, :, cde), v(VSYM[cde])); end
    P = Par(dt === nothing ? gf(:dtime) : Float64(dt), gf(:thet), gf(:cadv), gf(:epsadv), gf(:sigd), gf(:sige), gf(:phil0), gf(:phil2),
            gf(:phir0), gf(:phir2), gf(:cgrav), gf(:g), gf(:xmtot), gf(:alph), gf(:tau), gf(:xscale),
            gf(:ql0), gf(:ql1), gf(:cq1), gf(:cq2), gf(:cqvis), gf(:gam), gf(:gmu), gf(:ck), gf(:cm0),
            gf(:csige), gf(:chif0), gf(:ratio), gf(:cc), gf(:cpi), gf(:csigr), gf(:geddr), gf(:ximr),
            gf(:delml), gf(:delmr), gf(:delrl), gf(:delrr), gf(:dl), gf(:dr),
            gi(:ibet), gi(:ladx), gi(:lgrid), gi(:leibc), gi(:leobc), gi(:llibc), gi(:llobc), gi(:lribc),
            gi(:lrobc), gi(:lopac))
    Ctx(ngrs, ngre, neqn, gi(:mu), lrad, ngrs - 2, ngre - ngrs + 5, mgr, ix, code, base,
        v(:unom), v(:frnom),
        v(:ro), v(:xmeo), v(:do), v(:uo), v(:ego), v(:ero),
        v(:fro), v(:egro),
        v(:pgo), v(:aso), v(:chifo), v(:xkeo),
        v(:xkpo), v(:plfo), v(:dvolo), v(:xmo), v(:fedd), v(:bri),
        (la = v(:lady); lb = get!(() -> zeros(Int, length(la)), BUFS, :lady)::Vector{Int}; lb .= la; lb), v(:wt), v(:yscl), P, etotal, dualswitch(ngrs, ngre, lrad, ηsw),
        CONVECTION[], ALPHA_MLT[], convcoef(ngrs, ngre, base))
end

const OPAC_NOCACHE = Ref(false)
"Interface ∇_ad, c_P, δ, κ_R on the current iterate (base) for the MLT flux; zeros if convection is off."
function convcoef(ngrs, ngre, base)
    mgr = TF.MGR
    cv = Tuple(get!(() -> zeros(mgr), BUFS, nm)::Vector{Float64} for nm in (:cvnad, :cvcp, :cvdel, :cvkap))
    foreach(a -> fill!(a, 0.0), cv)
    CONVECTION[] || return cv
    m = MATERIAL[]; m === nothing && error("MLT convection needs a material with material_thermo (e.g. Z0Gas)")
    MLT_AD[] || @inbounds for k in ngrs+1:ngre       # lagged coefficients (not needed when MLT_AD)
        ρi = 0.5 * (base[k-1, 3] + base[k, 3]); Ti = 0.5 * (base[k-1, 5] + base[k, 5])
        (ρi > 0 && Ti > 0) || continue          # unphysical iterate: no MLT (the step will be rejected)
        ZONE[] = k                                   # composition of the zone above the interface
        nad, cp, δ = material_thermo(m, ρi, Ti)
        cv[1][k] = nad; cv[2][k] = cp; cv[3][k] = δ; cv[4][k] = material_opac(m, ρi, Ti)[1]
    end
    # MLT is a diffusion-limit theory: no convective flux above optical depth TAU_MLT (Rosseland, integrated from
    # the outer boundary on this iterate). Weight w(τ) per interface: a smooth step in ln τ from 0 at TAU_MLT/2 to 1
    # at 2·TAU_MLT (it was a sharp cut at TAU_MLT: a flux that jumps when τ crosses the cut makes Newton cycle).
    cvw = get!(() -> zeros(mgr), BUFS, :cvw)::Vector{Float64}; fill!(cvw, 0.0)
    # (PERF512) for a material whose Float64 opacity has no side effects (`opac_pure(m)`), the zones' χ first on MATTHREADS tasks, then the
    # serial sum. The serial loop calls material_opac without setting ZONE[] (a material that reads ZONE[] sees the zone left by the last
    # sweep, for every k), so the tasks use that same ZONE[]; a zone whose evaluation throws is evaluated again in the loop, which then throws
    # as before. The final zone visited by the downward loop (`klast`) is left to the serial loop so per-zone caches (NLOIT) end in the exact
    # same state. Bit-identical.
    χτ = get!(() -> zeros(mgr), BUFS, :cvchi)::Vector{Float64}; okτ = get!(() -> zeros(Bool, mgr), BUFS, :cvchiok)::Vector{Bool}; fill!(okτ, false)
    if opac_pure(m)
        z0 = ZONE0[]
        klast = 0
        @inbounds for k in ngrs+1:ngre
            if base[k, 3] > 0 && base[k, 5] > 0; klast = k; break; end
        end
        OPAC_NOCACHE[] = true
        try zone_sweep!(ngrs+1:ngre, k -> k != klast && @inbounds(base[k, 3] > 0 && base[k, 5] > 0)) do k
            ZONE[] = z0
            try
                y = @inbounds material_opac(m, base[k, 3], base[k, 5])[1]
                y isa Float64 && (@inbounds χτ[k] = y; @inbounds okτ[k] = true)
            catch
            end
        end
        finally OPAC_NOCACHE[] = false; ZONE0[] = z0 end
    end
    τ = 0.0
    @inbounds for k in ngre:-1:ngrs+1
        (base[k, 3] > 0 && base[k, 5] > 0) &&
            (τ += (okτ[k] ? χτ[k] : material_opac(m, base[k, 3], base[k, 5])[1]) * base[k, 3] * (base[k+1, 1] - base[k, 1]))
        x = clamp((log(max(τ, 1e-300)) - log(TAU_MLT[] / 2)) / log(4.0), 0.0, 1.0)
        cvw[k] = x * x * (3 - 2x)
        if MLT_NOINFALL[]          # no MLT flux in infalling gas: smooth step from inflow 10 km/s (MLT on) to 30 (off)
            y = clamp((-base[k, 4] - 1e6) / 2e6, 0.0, 1.0)
            cvw[k] *= 1 - y * y * (3 - 2y)
        end
    end
    cv
end

# ------------------------------------------------------------------ workspace
const WSNAMES = (:rn, :xmen, :dn, :un, :tn, :ern, :frn, :pgn, :egn, :asn, :chifn, :xken, :xkpn, :plfn,
    :egrn, :r, :u, :fr, :xm, :dmdt, :urel, :rmu, :rmup1, :d, :pg, :eg, :as, :er, :egr, :chif, :xke,
    :xkp, :plf, :dvol, :dvoln, :avchi, :qf, :dudr, :db, :dsn, :fl, :ub, :usn, :flu, :fle, :egb, :egsn,
    :dso, :uso, :egrso, :egso, :erso, :frso, :xso1, :xso2,       # old-level limiter slopes (oldslope!)
    :erb, :ersn, :frb, :frsn, :qe, :qeo, :qF, :qFo, :dfm, :dfe, :xnu, :xnc, :ss, :fcv, :egbg, :egsng, :xnen, :snuc, :srcw, :vcv, :lcv, :nabv, :nadv,
    :mnad, :mcp, :mdl, :mkp,                # MLT coefficients of the threaded pre-pass (MATTHREADS)
    :xs1, :xs2, :xb1, :xb2, :xsn1, :xsn2, :cpS1, :cpS2, :cpL, :cpQ, :cpY)     # COUPLED: species iterate, face values, slopes, rates, cooling
struct WS{T}
    a::NamedTuple{WSNAMES,NTuple{length(WSNAMES),Vector{T}}}
    xnuo::Vector{Float64}
end
const WSCACHE = Dict{DataType,Any}()
getws(::Type{T}, mgr) where {T} =
    get!(() -> WS{T}(NamedTuple{WSNAMES}(Tuple(zeros(T, mgr) for _ in WSNAMES)), zeros(mgr)), WSCACHE, T)::WS{T}

"""
Branch freezing (Phase 4 item 5 experiment): FREEZE[] = 0 normal, 1 record the branch
decisions (limiter a·b ≥ 0, sign of a+b, upwind direction, viscosity on/off) of this
evaluation, 2 replay the recorded decisions.  Default 0 = TITAN's scheme.
"""
const FREEZE = Ref(0)
const DEC = Dict{Tuple{Symbol,Int},Vector{Bool}}()
# stalldiag (SDDEC=1, default off): count, per (kind, id, zone), how often decide() returns a different value than at its previous call
# (successive residual evaluations: Newton iterates, stages, steps); SDLAST/SDFLIP indexed by the zone array index
const SDDEC = Ref(get(ENV, "SDDEC", "0") == "1"); const SDLOCK = ReentrantLock()
const SDLAST = Dict{Tuple{Symbol,Int},Vector{Int8}}(); const SDFLIP = Dict{Tuple{Symbol,Int},Vector{Int32}}()
function sd_dec(kind::Symbol, id::Int, k::Int, x::Bool)
    key = (kind, id); v = get(SDLAST, key, nothing)
    if v === nothing
        lock(SDLOCK) do
            haskey(SDLAST, key) || (SDLAST[key] = zeros(Int8, TF.MGR); SDFLIP[key] = zeros(Int32, TF.MGR))
        end
        v = SDLAST[key]
    end
    @inbounds begin
        o = v[k]; n = x ? Int8(1) : Int8(2)
        o != 0 && o != n && (SDFLIP[key][k] += 1)
        v[k] = n
    end
    nothing
end
const SDPREV = Dict{Tuple{Symbol,Int},Vector{Int32}}()
"flips of decide() in the zones ks since the last call, as text"
function sd_dec_text(ks)
    io = IOBuffer()
    for (key, f) in SDFLIP
        p = get!(() -> zeros(Int32, length(f)), SDPREV, key)
        d = [(k, f[k] - p[k]) for k in ks if f[k] != p[k]]
        isempty(d) || print(io, " ", key[1], ":", key[2], "{", join(("$(k - first(ks) + 1)x$n" for (k, n) in d), ","), "}")
        p .= f
    end
    String(take!(io))
end
@inline function decide(kind::Symbol, id::Int, k::Int, x::Bool)
    m = FREEZE[]
    m == 0 && (SDDEC[] && sd_dec(kind, id, k, x); return x)
    v = get!(() -> fill(false, TF.MGR), DEC, (kind, id))
    m == 1 && (@inbounds v[k] = x; return x)
    @inbounds return v[k]
end

"""
EXCAP[] > 0: factor on the gas-radiation exchange of the two energy rows, so that c κ_E ρ Δt never exceeds EXCAP[]. In an opaque zone that product can be
1e20 and the rows carry it as round-off; E = aT⁴ there holds to about 1/EXCAP anyway. The same factor in both rows keeps energy conserved. Differentiated with κ_E and ρ (it was a constant of the iteration: in the zones where it is below 1 the Jacobian then lacked the dependence of the exchange on κ_E; the product factor × κ_E E ρ = EXCAP E/(c Δt) does not depend on κ_E at all).
"""
const EXCAP = Ref(0.0)
"""
SUBZONE (popiii_m2c `SUBZONE=1[,τ_f]`, default off; off = bit-identical): optical depth of an interface across which an unresolved
front (photosphere / recombination front) lies, from a two-layer sub-zone picture instead of the harmonic mean of the two zones' χ.

TITAN's interface coefficient avchi = 2/(1/χ₋ + 1/χ₊) gives the interface optical depth Δτ_T = (M₋ + M₊) avchi (M± the half-zone
columns ½ ρ dV / r^μ). With a = χ₋M₋, b = χ₊M₊ the two half-zone optical depths (equal M: Δτ_T = 4ab/(a + b)), Δτ_T → 2b(1 + M₋/M₊)
for a ≫ b: across a front the opaque half-zone is given the transparent zone's opacity. The interface flux is then ∝ 1/χ of the
transparent zone, the zone that (outside a photosphere) is irradiated, steeply T dependent (χ ∝ T^q, q ~ 10–20) and lagged by the
operator-split chemistry, so the emergent flux of the hot surface follows the cold neighbour's opacity one to one.

Sub-zone picture (a front between the zone centres, thinner than any zone): the light leaving the opaque side originates within
an optical depth of order unity of its surface (Eddington–Barbier: τ = 2/3, e.g. Mihalas 1978, Stellar Atmospheres,
1978stat.book.....M), so the opaque half-zone contributes its own depth only up to τ_f, τ_f a/(a + τ_f), and the transparent
half-zone its own b. τ_f (default 2, `SUBZONE=1,τ_f`) is the depth of the opaque zone's centre below its unresolved surface, a
closure constant of order unity; it also sets the width of the front window s(x) below. On zdiagB's front τ_f = 2/3 does not
remove the orbit, τ_f = 1 and 2 do (docs/POPIII_SUBZONE_FRONT.md):
    D₋₊ = τ_f a/(a + τ_f) + b   (opaque side −),   D₊₋ = τ_f b/(b + τ_f) + a   (opaque side +).
Away from a front the two half-zones are in series (a + b), the form that keeps the thick diffusion limit (Larsen, Morel & Miller
1987, 1987JCoPh..69..283L); TITAN's harmonic mean is kept there unchanged. The front weight
    s(x) = x⁴/(x⁴ + τ_f⁴),  ω₋₊ = s(a)(1 − s(b)),  ω₊₋ = s(b)(1 − s(a))  (ω₋₊ + ω₊₋ ≤ 1)
is ≈ 1 only where one half-zone is opaque (≳ τ_f) and the other transparent (≲ τ_f), so
    Δτ = (1 − ω₋₊ − ω₊₋) Δτ_T + ω₋₊ D₋₊ + ω₊₋ D₊₋,   avchi_SZ = Δτ/(M₋ + M₊).
Δτ is a rational (C∞) function of χ±, ρ±, r: the Newton's Duals carry it, and the interface still couples only its two zones (the
block structure is unchanged). As the front moves into the transparent zone (b grows through τ_f), ω falls smoothly to 0 and the
next interface takes over. In thick–thick interiors the change is ≲ (τ_f/min(a,b))⁴, in thin–thin envelopes ≲ (max(a,b)/τ_f)⁴.
The same avchi enters the flux row and the radiation-force work of the energy rows (energy stays conserved).
"""
const SUBZONE = Ref(false); const SZTAUF = Ref(2.0)
"(SUBZONE diagnostic) per interface: front weight ω₋₊ + ω₊₋ and Δτ_SZ/Δτ_T (values) of the last residual evaluation."
const SZOMEGA = Float64[]; const SZRATIO = Float64[]
@inline function sz_tau(a, b, ΔτT, τf)
    τ4 = τf^4; a4 = a^4; b4 = b^4
    sa = a4 / (a4 + τ4); sb = b4 / (b4 + τ4)
    ωab = sa * (1 - sb); ωba = sb * (1 - sa)
    Dab = τf * a / (a + τf) + b; Dba = τf * b / (b + τf) + a
    return (1 - ωab - ωba) * ΔτT + ωab * Dab + ωba * Dba, ωab + ωba
end
"SUBZONE: replace avchi of the interior interfaces ngrs+1..ngre by the sub-zone front form (see SUBZONE); θ-level χ, ρ, dV, r^μ."
function sz_avchi!(avchi, chif, d, dvol, rmu, ngrs, ngre)
    τf = SZTAUF[]
    length(SZOMEGA) < ngre + 1 && (resize!(SZOMEGA, ngre + 1); resize!(SZRATIO, ngre + 1)); fill!(SZOMEGA, 0.0); fill!(SZRATIO, 1.0)
    @inbounds for k in ngrs+1:ngre
        Mm = 0.5 * d[k-1] * dvol[k-1] / rmu[k]; Mp = 0.5 * d[k] * dvol[k] / rmu[k]
        ΔτT = (Mm + Mp) * avchi[k]
        Δτ, ω = sz_tau(chif[k-1] * Mm, chif[k] * Mp, ΔτT, τf)
        avchi[k] = Δτ / (Mm + Mp)
        SZOMEGA[k] = val(ω); SZRATIO[k] = val(Δτ) / val(ΔτT)
    end
    nothing
end
"""
SLABFLUX (popiii_m2c `SLABFLUX=1[,n]`, default off; off = bit-identical): interface opacity from the exact steady solution of the
grey moment equations across the slab between the two zone centres (docs/POPIII_SLABFLUX.md), in place of TITAN's harmonic mean.

In the column-mass coordinate of the flux row (dm = ρ dr, M± = ½ρ±dV±/r^μ the half-zone columns, ΔM = M₋ + M₊, dτ = χ dm) the static
first moment d(q f E)/dτ = −q F/c with no source inside the slab (finite-volume: the exchange of each half-zone is in its zone's energy
rows) has F constant, so F = c Δ(q f E)/Δτ exactly, Δτ = ∫χ dm. The whole slab problem is then the optical depth between the centres
for an opacity profile through the two centre values. Profile family κ ∝ (τ + τ₀)^n (n ≥ 1; τ₀ fixed by χ± and ΔM), which for a grey
LTE layer at constant flux with χ ∝ T^q (T⁴ linear in τ) is exact with n = q/4; n → ∞ is the exponential in optical depth,
n = 1 the exponential in column (log-linear in m), n = 2 gives the geometric mean:
    Δτ = ΔM χ_n,   χ_n = χ₊ ℰ(x/n) / ℰ(−x(n−1)/n),   x = ln(χ₋/χ₊),   ℰ(y) = (e^y − 1)/y,
    n → ∞:  χ_∞ = χ₊ x/(1 − e^{−x}) = ln(χ₋/χ₊)/(1/χ₊ − 1/χ₋)     (the Scharfetter–Gummel / Bernoulli form, B(−x) = x/(1 − e^{−x})).
χ_n is symmetric in χ±, lies between the harmonic and the logarithmic mean (H ≤ χ_∞ ≤ G ≤ χ_1), equals χ when χ₋ = χ₊ (Δτ = a + b,
the series of the two half-zones: diffusion limit F = −(c/3χρ)∂E/∂r, Larsen, Morel & Miller 1987) and is ≤ max χ± (Δτ ≪ 1 stays
negligible against the r^μ Δ(q f E) term: the VEF/FLUXTIE free-streaming limit is untouched). It depends on χ± only (not on E, so the
flux stays diffusive in E: ∂F_k/∂E_k > 0 is never produced), is C∞ (ℰ by its series for |y| < 1e-2) and couples only the two zones.
Default n → ∞ (`SLABFLUX=1`), the smallest optical depth of the family (closest to TITAN's); `SLABFLUX=1,n` sets n.
The same avchi enters the flux row and the radiation-force work of the energy rows. With SUBZONE also set, SUBZONE's front form is
applied on top (its Δτ_T is then this Δτ).
"""
const SLABFLUX = Ref(false); const SFN = Ref(Inf)
"(SLABFLUX diagnostic) per interface: Δτ_SF/Δτ_T (values) of the last residual evaluation."
const SFRATIO = Float64[]
"ℰ(y) = (e^y − 1)/y, C∞ (series for |y| < 1e-2: truncation y⁷/8!, derivatives without the 1/y² cancellation)."
@inline function sf_E(y)
    abs(val(y)) < 1e-2 && return 1 + y * (1 / 2 + y * (1 / 6 + y * (1 / 24 + y * (1 / 120 + y * (1 / 720 + y / 5040)))))
    return expm1(y) / y
end
"SLABFLUX: the interface opacity χ_n of zones with χ₋ = χm, χ₊ = χp (n = Inf: exponential in optical depth)."
@inline function sf_mean(χm, χp, n)
    x = log(χm / χp)
    isinf(n) && return χp / sf_E(-x)
    return χp * sf_E(x / n) / sf_E(-x * (n - 1) / n)
end
"SLABFLUX: replace avchi of the interior interfaces ngrs+1..ngre by χ_n (see SLABFLUX); θ-level χ."
function sf_avchi!(avchi, chif, ngrs, ngre)
    n = SFN[]
    length(SFRATIO) < ngre + 1 && resize!(SFRATIO, ngre + 1); fill!(SFRATIO, 1.0)
    @inbounds for k in ngrs+1:ngre
        a0 = avchi[k]
        avchi[k] = sf_mean(chif[k-1], chif[k], n)
        SFRATIO[k] = val(avchi[k]) / val(a0)
    end
    nothing
end
"""
SLABTIE (popiii_m2c `SLABTIE=1`, needs SLABFLUX with n = ∞; default off; off = bit-identical): the outward donor factor of the FLUXTIE row
from the slab solution. FLUXTIE (F_k = h c E upwinded, where τ from the outer boundary ≤ 1) takes h = H/J of the formal solution at the
donor zone's *centre*. When the donor is the opaque side of an unresolved photosphere (zone Δτ 1e4, its centre τ ~ 1e4 below the
interface in the formal solution's piecewise-constant picture) h ≈ 1/(3τ) ~ 1e-5 and the interface flux is choked although, in the
exponential-in-τ slab profile of SLABFLUX, the donor centre lies only τ_i (~0.1–10) below the interface. SLABTIE uses the steady slab
solution from the donor centre to the interface (f E − f E_s = F τ_i/c) with a Marshak/Eddington free surface (F = c E_s/2):
    h_slab(τ) = f₋/(τ + 2 f₋),     h = max(h_FS + h_slab(τ_i) − h_slab(a), 0)
(h_FS the formal solution's donor h; a = χ₋M₋ the donor's half-zone depth in the formal solution's piecewise-constant picture, so the
correction is the change of the Marshak slab estimate when the donor centre is moved from depth a to the slab depth τ_i),
τ_i = Δτ_SF ln(1 + (R − 1) M₋/ΔM)/ln R (R = χ₋/χ₊; the column M₋ of the donor's half-zone in the χ ∝ e^{−βτ} profile; = a for R = 1).
No front (R = 1): h = h_FS exactly. Front (a ≫ τ_i): h ≈ h_FS + f/(τ_i + 2f): 1/(3τ_i + 2) for a thick slab (diffusion), ½ for τ_i → 0
(Marshak). The inward term is unchanged. Local (zones k−1, k), smooth in χ±, ρ±, r (Duals) except the max at h = 0.
"""
const SLABTIE = Ref(false)
"""
TIEFRAME (popiii_m2c `TIEFRAME=1`; provisional switch, default off = bit-identical; docs/POPIII_ENERGY.md §3): FLUXTIE and the outer
boundary row in the frame of the unknowns. E and F of the rows are comoving-frame moments (the energy row carries ∇·(uE) and P:∇u); the
formal solution that gives h = H/J, f = K/J and g is static, i.e. it describes the field of an observer at rest. To O(u/c)
(Mihalas & Mihalas 1984 §91; Krumholz, Klein & McKee 2007, ApJ 667, 626, eqs. 31–33; Mihalas & Auer 2001, JQSRT 71, 61)
    F_lab = F₀ + u (E₀ + P₀),   E_lab = E₀ + 2 u F₀/c²,
so the closure F_lab = h c E_lab is, for the unknowns, F₀ = [h c − (1 + f − 2h²) u] E₀ (h = f = 1: F₀ = c E₀, a beam stays a beam;
h = 0, f = 1/3: F₀ = −(4/3) u E₀, the flux a moving observer sees in an isotropic field). Without the u terms the row said F₀ = h c E₀:
in the CMB-dominated cloud (h ~ 1e-6…1e-3) the energy row's ∇·(uE) + (E/3)∇·u was then balanced by nothing and the infall compressed the
CMB like a trapped gas (E_r up to 6 E_CMB at 0.2 pc, −1e5 L☉ through the outer boundary, docs/POPIII_ENERGY.md §3).
Discrete form at interface k (all at the θ-level of the row, nothing lagged), chosen so that a uniform isotropic field is a steady state of
the energy row to round-off, also on a moving mesh (well-balanced), and the mesh-sweep rezoning term −V̇_k E_up is upwinded in the interface motion (M-matrix):
    r² F₀ = r_c² c (h E)_donor + ṁ_k (E/ρ)ᵇ_k − V̇_k E_up − r² u ½(f₋E₋ + f₊E₊) + 2 u r_c² (h² E)_donor,
ṁ_k (E/ρ)ᵇ_k the advective flux of the energy row itself (dmdt·erb, = −r² ρᵇ u_rel (E/ρ)ᵇ), V̇_k = (r_new³ − r_old³)/(3Δt) the volume the
interface sweeps per unit time (the same radii as ΔV_new − ΔV_old of the row: geometric conservation; E_up = E₊ for V̇_k ≥ 0, E₋ for V̇_k < 0); the two together are −r² u E.
Outer boundary: F_lab = g c E_s − π(2g + 1) I_in becomes F₀ = [g c E_s − π(2g + 1) I_in](1 + 2 g u/c) − (1 + f) u E.
"""
const TIEFRAME = Ref(false)
"""
ECOMPAT[] (default false; selected by ETOTAL=0 in popiii_m2c and phase6_run, radiation grey or lrad = 0): compatible energy rows in EVERY zone, no row switching and no E − K subtraction
(docs/POPIII_ENERGY.md §4). The gas row is the first law (GASROW's row) made the exact discrete adjoint of the momentum row, so that
Σ_k (gas row + radiation row)_k + Σ_j u_j (momentum row)_j is the conservative total-energy row sum identically:
  * viscous work q_k r̄_k^{3μ/2} (u_{k+1} r_{k+1}^{−μ/2} − u_k r_k^{−μ/2}), the adjoint of the momentum row's r^{−μ/2} Δ(r̄^{3μ/2} q);
  * the kinetic energy the momentum row's discretisation does not keep, as cancellation-free products of differences, returned as heat:
      time    −(2θ − 1) ¼ ½(m_n + m_o) Σ_faces (u_n − u_o)² / Δt + (θ² − θ + ½) ¼ (m_n − m_o) Σ_faces (u_n − u_o)² / Δt
              (BDF stage of TR-BDF2: + D_k/Δt, D = K(y*) − (A K_γ − B K_n), `KSTAR`: the reference state's kinetic energy is not the BDF2
              combination of the stage energies, and the step conserves only with the latter)
      upwind  ½(ṁ_k + ṁ_k+1) (u_k+1 − u_k) (½(u_k + u_k+1) − û_k) − (G_k+1 − G_k)    û the van Leer value of the momentum flux,
              G_j = ṁ_j (Δu)²_upwind/8 the flux of the sub-zone kinetic energy (zone kinetic energy ¼m(u_k² + u_k+1²) = ½m ū² + ⅛m(Δu)²)
              with the relative mass flux (`ksflux`). The first term alone is not sign-definite (third order in Δu, it cools where the
              velocity gradient steepens upstream); with the transport term the heating is ½|ṁ̄| [½(a − b)² + ¼(a − b)³/(a + b)] ≥ ⅛|ṁ̄|(a − b)²
              for monotone data (a, b the velocity differences of the zone and of its upwind neighbour) and ⅜|ṁ̄|a² + ⅛|ṁ̄|b² at an extremum
  * the old-level energy is the one the previous solve left, e_old = (e_g + E/ρ)_old − E_old/ρ_old, not a re-evaluation of the EOS
    (COUPLED zones: `cpl_prepare!` re-evaluates e_g(ρ, T, x) with the new step's expansion terms);
  * no gravity term: the work of gravity is what the momentum row's force does, Σ_j u_j g_j ½(m_j−1 + m_j).
The radiation row is the conservation law of the lab-frame energy in the comoving unknowns to O(u/c),
  Δ(E ΔV)/Δt + Δ(r² F₀) − Δ(ṁ (E/ρ)ᵇ) + Δ(r² u (fE)‾) + Σ_faces u_j χ̄_j F_j/c ½m_k = sources,
i.e. the comoving row plus u × (radiation momentum equation): the comoving work f E Δ(r²u) + (1 − 3f) E u/r ΔV is the adjoint of the
MOMENT equation's ∇·P, so the comoving row conserves only where that equation is the flux row; where FLUXTIE replaces it by a closure the
difference u (∇·P + χρF/c) was created (3.3e3 L☉ in the outer cloud of the 100-zone run). Use with TIEFRAME.
"""
const ECOMPAT = Ref(false)
"""
KSTAR[] (ECOMPAT; set by Trbdf2.bdf_old!, emptied by its restore!): per zone D_k = K(y*) − (A K_γ − B K_n) [erg sr⁻¹], the kinetic energy
¼m(u_k² + u_k+1²) of the BDF stage's reference state y* = A y_γ − B y_n (built on m and m u) minus the same combination of the stage
kinetic energies, = A B ⅛ m Σ_faces (u_γ − u_n)² for equal masses. A stage row that differences against y* conserves E(y*), the step conserves
A E_γ − B E_n; the gas row therefore takes D_k out of its old-level energy. Empty in every other stage.
"""
const KSTAR = Ref(Float64[])
"ECOMPAT: flux G_j = ṁ_j (u_up+1 − u_up)²/8 of the sub-zone kinetic energy through interior interface j (upwind zone: j for ṁ_j ≥ 0, else j − 1); 0 at the boundaries."
@inline function ksflux(j, u, dmdt, ngrs, ngre)
    ngrs < j <= ngre || return zero(eltype(u))
    up = decide(:up, 20, j, val(dmdt[j]) >= 0) ? j : j - 1          # (branch recorded/replayed with the other upwind decisions of J v)
    dmdt[j] * 0.125 * (u[up+1] - u[up])^2
end
"""
ERLAB[] (provisional, ENV ERLAB=1, default false; lrad = 1, grey): ECOMPAT's radiation half alone (docs/POPIII_ENERGY.md §18). The radiation
energy row is the lab-frame conservation law of ECOMPAT (comoving row + u × radiation momentum equation); the gas rows are the stack's
(dual-energy: total row where c.etotal && !useint, the GASROW first-law row elsewhere). Use with TIEFRAME (the driver switches it on).
"""
const ERLAB = Ref(false)
"""
GFV[] (provisional, ENV GFV=1, default false; μ = 2 self-gravity): finite-volume weight in the gas momentum row (docs/POPIII_ENERGY.md §19).
The row's gravity term is the node mass times the interface's acceleration, g_j ½(m_j−1 + m_j), exact where the interface moves with the gas.
GFV replaces it by ∫ρ g dV between the centres of the two zones for piecewise constant density and r² g = 4πG M_j across the node,
    4πG M_j ½[ρ_j−1 (r_j − r_j−1) + ρ_j (r_j+1 − r_j)],
the weight that the flux form of the inertia (ṁ Δū between the zone centres) balances in steady flow through the mesh. Zone k's part of
interface j's weight is `gwt(k, j, …)`; the uniform-field term (1 − μ/2) P.g keeps the node mass. Not applied to the PEXT free-surface row.
"""
const GFV = Ref(false)
"Zone k's part of the mass that the gravity term of interface j's momentum row multiplies by g_j: ½ m_k, or with GFV ½ r_j^μ ρ_k (r_k+1 − r_k)."
@inline gwt(k, j, d, dvol, r, rmu) = GFV[] ? 0.5 * rmu[j] * d[k] * (r[k+1] - r[k]) : 0.5 * d[k] * dvol[k]
"""
GVAR[] (provisional, ENV GVAR=1, default false; μ = 2 self-gravity): the gravity force of the gas momentum row as the gradient of the exact
binding energy of the piecewise-constant density (docs/POPIII_ENERGY.md §19),
    F_j = ∂W/∂r_j at fixed zone masses,   W = −4πG Σ_k ∫_zone k M(r) ρ_k r dr,   M(r) = M_j + ρ_k (r³ − r_j³)/3 on either side of interface j,
which is ∫ρ g ξ_j dV with the displacement field ξ_j of interface j moving at uniform zone densities (r² ξ_j linear in r³ in both zones):
    F_j = 4πG r_j² [ρ_j−1 ∫ M (r³ − r_j−1³)/(r² (r_j³ − r_j−1³)) dr + ρ_j ∫ M (r_j+1³ − r³)/(r² (r_j+1³ − r_j³)) dr]      (`gvar_in`, `gvar_out`).
A state-only force with a binding energy when the interfaces move with the gas (the finite-volume weight GFV has none), and within a few
per cent of the finite-volume weight where mass passes. The external mass (EXTMASS) and the uniform field P.g keep the node mass.
Not applied to the PEXT free-surface row.
"""
const GVAR = Ref(false)
"GVAR: r_b² ρ ∫_a^b M(r) (r³ − a³)/(r² (b³ − a³)) dr, M(r) = M − ρ (b³ − r³)/3: the zone below interface b (M the mass enclosed by b), per 4πG."
@inline function gvar_in(a, b, ρ, M)
    A = M - ρ * b^3 / 3
    b^2 * ρ * (b - a) * (A * (b + 2a) / (2b) + ρ / 30 * (2b^3 + 4b^2 * a + 6b * a^2 + 3a^3)) / (b^2 + a * b + a^2)
end
"GVAR: r_a² ρ ∫_a^b M(r) (b³ − r³)/(r² (b³ − a³)) dr, M(r) = M + ρ (r³ − a³)/3: the zone above interface a (M the mass enclosed by a), per 4πG."
@inline function gvar_out(a, b, ρ, M)
    A = M - ρ * a^3 / 3
    a^2 * ρ * (b - a) * (A * (2b + a) / (2a) + ρ / 30 * (3b^3 + 6b^2 * a + 4b * a^2 + 2a^3)) / (b^2 + a * b + a^2)
end
"""
GFV/GVAR: zone k's part of the self-gravity force of interface j's momentum row (j = k: the zone above the interface, j = k + 1: below it);
G4 = 4πG, xmj the enclosed mass of the row, xej its external part.
"""
@inline function gfz(k, j, d, dvol, r, rmu, G4, xmj, xej)
    GVAR[] || return G4 * xmj / rmu[j] * gwt(k, j, d, dvol, r, rmu)
    (j == k ? G4 * gvar_out(r[k], r[k+1], d[k], xmj - xej) : G4 * gvar_in(r[k], r[k+1], d[k], xmj - xej)) + G4 * xej / rmu[j] * 0.5 * d[k] * dvol[k]
end
@inline gxe(rj) = EXTMASS[] === nothing ? zero(rj) : EXTMASS[](rj)
"SLABTIE: optical depth from the donor centre (χm, column Mm) to the interface in SLABFLUX's exponential-in-τ profile."
@inline function st_tauint(χm, χp, Mm, Mp)
    ΔM = Mm + Mp; φ = Mm / ΔM; x = log(χm / χp)
    Δτ = ΔM * sf_mean(χm, χp, Inf)
    abs(val(x)) < 1e-6 && return Δτ * (φ + 0.5 * φ * (1 - φ) * x)
    return Δτ * log1p(expm1(x) * φ) / x
end
"""
TIESMOOTH (popiii_m2c `TIESMOOTH=ε`, with SLABTIE; default 0 = off, off = bit-identical; docs/POPIII_MULTIGROUP.md §5, Newton fix 1): the
max(h, 0) of SLABTIE's donor factor replaced by the C∞ ramp s(h) = ½(h + √(h² + ε²)) (≥ 0, → h for h ≫ ε, s(0) = ε/2), so that the
row's derivative is continuous where the slab correction drives h through 0.
"""
const TIESM = Ref(0.0)
@inline sramp(x, ε) = 0.5 * (x + sqrt(x * x + ε * ε))
"SLABTIE: the outward donor factor h (see SLABTIE)."
@inline function st_hdonor(hfs, fm, χm, χp, Mm, Mp)
    τi = st_tauint(χm, χp, Mm, Mp)
    TIESM[] > 0 && return sramp(hfs + fm / (τi + 2fm) - fm / (χm * Mm + 2fm), TIESM[])
    return max(hfs + fm / (τi + 2fm) - fm / (χm * Mm + 2fm), 0.0)
end
@inline function excap(cc, κ, ρ, dt)
    EXCAP[] > 0 || return 1.0
    return min(1.0, EXCAP[] / (cc * κ * ρ * dt))
end
include(joinpath(@__DIR__, "MGRad.jl"))                 # MGROUPS (default off)
"(MGROUPS with SLIMP) (q, t) of zone k at the iterate from the grey E_r and group MGSLG's E_g, F_g, F_g+1 (PopIIIMaterial.sl_qt_mg)."
const MGSLQ = Ref{Any}(nothing); const MGSLG = Ref(2)

"""
Limiter slope of the OLD level, from the old cell values with exactly the formula of the first loop of `advect!` (plain Float64: the old level is a constant
of the iteration, no branch decisions to record). Replaces the slopes TITAN stored in COMMON (dso, uso, egso, egrso, erso, frso, XSSO), which were the slopes of
the penultimate Newton iterate (store_side! in assemble!), zero for the gas energy under GASROW (zeroed by update!/peq), and 0 at the outer face for u.
"""
@inline function oldslope!(so, qo, cadv, eps, ngrs, ngre)
    @inbounds for k in ngrs:ngre+1
        a = qo[k] - qo[k-1]; b = qo[k+1] - qo[k]; sdq = a + b
        den = sdq + (sdq >= 0 ? eps : -eps)
        so[k] = a * b >= 0 ? cadv * a * b / den : 0.0
    end
    so
end

"van Leer interpolation (advectc: centred=true, advecti: centred=false)."
@inline function advect!(qb, qsn, qn, qo, qso, flow, θ, cadv, eps, ngrs, ngre, centred::Bool, id::Int = 0)
    T = eltype(qb)
    @inbounds for k in ngrs:ngre+1
        a = qn[k] - qn[k-1]; b = qn[k+1] - qn[k]; sdq = a + b
        den = sdq + (decide(:sgn, id, k, val(sdq) >= 0) ? eps : -eps)
        qsn[k] = decide(:lim, id, k, val(a * b) >= 0) ? cadv * a * b / den : zero(T)
    end
    if centred
        @inbounds for k in ngrs+1:ngre+1
            s = decide(:up, id, k, val(flow[k]) >= 0) ? 0.5 : -0.5
            qb[k] = (0.5 + s) * (θ * (qn[k-1] + 0.5qsn[k-1]) + (1 - θ) * (qo[k-1] + 0.5qso[k-1])) +
                    (0.5 - s) * (θ * (qn[k]   - 0.5qsn[k])   + (1 - θ) * (qo[k]   - 0.5qso[k]))
        end
    else
        @inbounds for k in ngrs:ngre
            s = decide(:up, id, k, val(flow[k]) >= 0) ? 0.5 : -0.5
            qb[k] = (0.5 + s) * (θ * (qn[k]   + 0.5qsn[k])   + (1 - θ) * (qo[k]   + 0.5qso[k])) +
                    (0.5 - s) * (θ * (qn[k+1] - 0.5qsn[k+1]) + (1 - θ) * (qo[k+1] - 0.5qso[k+1]))
        end
    end
    nothing
end

@inline function rescs(y, k, lady, yscl)
    lady == 1 && return (y[k+1] - y[k]) / yscl
    lady == 2 && return (y[k+1] - y[k]) / (y[k+1] + y[k])
    return y[k+1] / y[k] - y[k] / y[k+1]
end

# ------------------------------------------------------------------ pluggable materials
"""
Material used by the residual. `nothing` = TITAN's built-in ideal gas + lopac 3/4 (inline, as before).
Otherwise a material object for which `material_gas(m, ρ, T) -> (p_gas, e_gas, a_s)` and
`material_opac(m, ρ, T) -> (χ_flux, κ_E, κ_P)` have methods (see TitanMaterials); ρ, T may be Duals.
"""
const MATERIAL = Ref{Any}(nothing)
"Zone index of the current material call (set by `_material!` and by the fills), for materials with a composition per zone."
struct ZoneCell end
const ZONE = ZoneCell(); const ZONE0 = Ref(0); const ZTASK = Ref(false)
# `ZONE[]` is one global cell, except inside a threaded material sweep (ZTASK[]), where every task has its own
@inline Base.getindex(::ZoneCell) = ZTASK[] ? (task_local_storage(:titan_zone)::Int) : ZONE0[]
@inline Base.setindex!(::ZoneCell, k) = (ZTASK[] ? task_local_storage(:titan_zone, Int(k)) : (ZONE0[] = k); k)
"""
MATTHREADS = n > 1: the material sweep of the residual (`_material!`: EOS, nuclear source and opacity of every zone, in the Jacobian's
Dual numbers ~70 % of the implicit solve) on n tasks. The zones are independent and each writes only its own entries, so the result is
bit-identical to the serial sweep. 0 or 1 (the default): serial.
"""
const MAT_THREADS = Ref(parse(Int, get(ENV, "MATTHREADS", string(Threads.nthreads()))))
function material_gas end
"`material_gasne(m, ρ, T) -> (p_gas, e_gas, a_s, n_e)`; default: material_gas and n_e = 0 (no Compton exchange)."
material_gasne(m, ρ, T) = (material_gas(m, ρ, T)..., zero(ρ))
"Compton exchange rate per volume (erg cm⁻³ s⁻¹) of zone k at the new time level (gas → radiation)."
@inline function compton(k, W, P)
    COMPTON[] || return zero(W.tn[k])
    E = W.ern[k]; Tr = sqrt(sqrt(E * P.cc / (4P.csigr)))
    Λn = 4P.csige * W.xnen[k] * P.cc * E * P.ck * (W.tn[k] - Tr) / MEC2
    P.θ == 1 ? Λn : P.θ * Λn + (1 - P.θ) * LAMO[][k]
end
"Compton rate of the old (start-of-step) state, per zone, filled by `context()` (θ-weighting of the exchange)."
const LAMO = Ref(Float64[])
function compton_old!(ngrs, ngre)
    m = MATERIAL[]; buf = get!(() -> zeros(TF.MGR), BUFS, :lamo)::Vector{Float64}; fill!(buf, 0.0); LAMO[] = buf
    (COMPTON[] && m !== nothing) || return
    dO, tO, eO = v(:do), v(:to), v(:ero); cc, csigr, csige, ck = gf(:cc), gf(:csigr), gf(:csige), gf(:ck)
    @inbounds for k in ngrs:ngre
        (dO[k] > 0 && tO[k] > 0 && eO[k] > 0) || continue
        ne = material_gasne(m, dO[k], tO[k])[4]; Tr = sqrt(sqrt(eO[k] * cc / (4csigr)))
        buf[k] = 4csige * ne * cc * eO[k] * ck * (tO[k] - Tr) / MEC2
    end
end
function material_opac end
"(PERF512) true if `material_opac(m, ρ::Float64, T::Float64)` has no side effects (no per-zone caches), so it may run on several tasks with the same ZONE[]."
opac_pure(m) = false
"""
Nuclear (or other volumetric) energy source of a material, S4. `has_source(m)` (default false) switches it on;
`material_source(m, ρ, T) -> ε` [erg g⁻¹ s⁻¹, net of neutrino losses; ρ, T may be Duals; composition of zone `ZONE[]`,
held fixed during the Newton iteration] is added to the energy row of every formulation (combined gas+radiation, GASROW
gas row, total-energy row, pure-gas row) as −ε·ρ·ΔV, θ-weighted between the new level (implicit in T and ρ) and the old
level (`SRCO`, evaluated by `context()` on the start-of-step state), like the Compton exchange. Nothing is added, and no code path
changes, for materials without a source.
"""
has_source(m) = false
material_source(m, ρ, T) = zero(T)
"ε of the old (start-of-step) state per zone, filled by `context()` when the material has a source."
const SRCO = Ref(Float64[])
function source_old!(ngrs, ngre)
    m = MATERIAL[]; buf = get!(() -> zeros(TF.MGR), BUFS, :srco)::Vector{Float64}; fill!(buf, 0.0); SRCO[] = buf
    (m !== nothing && has_source(m)) || return
    dO, tO = v(:do), v(:to)
    # (PERF512) zones independent: on MATTHREADS tasks (zone_sweep!: the serial loop's ZONE[] and exception), bit-identical
    zone_sweep!(ngrs:ngre, k -> @inbounds(dO[k] > 0 && tO[k] > 0)) do k
        @inbounds buf[k] = material_source(m, dO[k], tO[k])
    end
end
"`material_thermo(m, ρ, T) -> (∇_ad, c_P, δ)` of gas + equilibrium radiation (Float64), for the MLT flux."
function material_thermo end
"""
MLT convection (Phase 7, Tom's choice: local, time-independent MLT). At interior interfaces
F_conv = ρ c_P T √(gδ) ℓ² ζ³ / (4√2 H_P^{3/2}),  ζ = √(U² + ∇ − ∇_ad) − U  (0 if ∇ ≤ ∇_ad),
ℓ = α H_P, U = 3acT³/(c_P ρ² κ ℓ²) √(8H_P/(gδ)); ∇ = Δln T/Δln P from the neighbouring zones (n level),
P = p_gas + f E_r. ∇_ad, c_P, δ, κ are evaluated on the current Newton iterate and held fixed during the
derivative sweep (lagged coefficients); the dependence on T, ρ, P and ∇ is differentiated exactly.
"""
const CONVECTION = Ref(false)
const ALPHA_MLT = Ref(0.8)
"Minimum Rosseland optical depth (from the outer boundary) for a nonzero MLT flux."
const TAU_MLT = Ref(1.0)
"""
MLT_AD[]: evaluate ∇_ad, c_P, δ and κ_R of the MLT flux on the current state inside the residual, so they enter the
Newton Jacobian (default true). false: take them from the Newton iterate and hold them constant in the Jacobian, as
before. Lagged coefficients converge to the same solution when Newton converges, but in the hydrogen ionisation zone
c_P changes by 10–30 % per 1 % in T and Newton then wanders instead of converging (finer adaptive run, 17 Oct).
"""
const MLT_AD = Ref(true)
"""
MLT_HS[]: in layers far from hydrostatic equilibrium (|Δln P| < ¾ of the hydrostatic −Δr/H_P, H_P = P/(ρg)) measure
∇ against the hydrostatic pressure scale, blending smoothly (default true); hydrostatic layers keep the actual Δln P.
Behind the accretion shock the pressure is set by the deceleration of the inflow, not by gravity, it is nearly flat,
and ∇ = Δln T/Δln P reached 110 and a convective flux of 2 L. (Using −Δr/H_P everywhere is not possible: deep
convection has ∇ − ∇_ad ~ 1e-6, below the 0.1–1 % discretisation difference between the two.)
"""
const MLT_HS = Ref(true)
"""
MLT_NOINFALL[]: switch the MLT flux off at interfaces where the gas falls inwards faster than 30 km/s (smoothly from
10 km/s). Mixing-length convection assumes a quasi-hydrostatic layer; in the supersonic free fall of a large inflow
domain (r_out = 30 R) the optical-depth cut alone left it on, and spurious convective fluxes of 10⁴–10⁶ L stalled the run.
Default off (the earlier runs, with r_out = 1.2 R, did not need it).
"""
const MLT_NOINFALL = Ref(false)
"""
MLT_RCAP[] (S5, default false = TITAN): mixing length ℓ = α (H_P⁻² + r⁻²)^(-1/2), i.e. α min(H_P, r) with a smooth minimum. At a regular
centre g = G m/r² → 0 and H_P = P/(ρ g) → ∞: ℓ = α H_P is then many stellar radii (24 R☉ at the first interface of a 4 R☉ star), an eddy larger
than the distance to the centre. F_conv ∝ ℓ² ζ³ H_P^(-3/2) √g does not diverge (ℓ² H_P^(-3/2) √g = α² √(P/ρ)), but the huge ℓ makes the flux a
function of ζ³ ∝ (∇ − ∇_ad)^(3/2) with ∇ − ∇_ad ~ 1e-8 in the whole inner core, stiff beyond the Newton tolerance (the S4 "L_conv > L" at
the first interface); with ℓ ≤ αr the superadiabaticity needed for the same flux is larger (∝ ℓ^(-4/3)) and the Newton problem well conditioned.
"""
const MLT_RCAP = Ref(false)
"""
MLT_HSDISC[] (S5, default false = TITAN): the hydrostatic Δln P of MLT_HS is taken from TITAN's own discrete momentum balance,
Δln P = −g M_k/(r_k² P̄), M_k = ½(ρ_{k−1}ΔV_{k−1} + ρ_k ΔV_k) (exactly the pressure difference the momentum row produces in equilibrium),
instead of −Δr/H_P with Δr = ½(r_{k+1} − r_{k−1}). The latter is the distance between the zone centres only for zones of similar size;
at the first interface of a regular centre (zone 1 is a sphere, r_{k−1} = 0) it is 2× too large, the ratio actual/hydrostatic falls into the
blending range (¼, ¾) of MLT_HS, and the blend weight (taken as a constant in the Jacobian) both spoiled ∇ at that interface (the S4 "L_conv > L") and
made Newton converge only linearly (factor ⅓ per iteration) in the central zone.
"""
const MLT_HSDISC = Ref(false)
"Function barrier: fills the n-level EOS/opacity workspace from the material (type-specialised loop)."
@inline function _matzone!(m, W, c, P, dn, tn, ern, k, hassrc)
    @inbounds begin
        if !(ForwardDiff.value(dn[k]) > 0 && ForwardDiff.value(tn[k]) > 0)   # unused ghost zones (e.g. ngrs-2)
            W.pgn[k] = zero(dn[k]); W.egn[k] = zero(dn[k]); W.asn[k] = zero(dn[k]); hassrc && (W.snuc[k] = zero(dn[k]))
            c.lrad && (W.chifn[k] = zero(dn[k]); W.xken[k] = zero(dn[k]); W.xkpn[k] = zero(dn[k]); W.plfn[k] = zero(dn[k]); W.egrn[k] = zero(dn[k]))
            return false
        end
        ZONE[] = k
        if COMPTON[]
            pg, eg, as_, ne = material_gasne(m, dn[k], tn[k]); W.xnen[k] = ne
        else
            pg, eg, as_ = material_gas(m, dn[k], tn[k])
        end
        W.pgn[k] = pg; W.egn[k] = eg; W.asn[k] = as_
        hassrc && (W.snuc[k] = material_source(m, dn[k], tn[k]))
        if c.lrad
            χ, κE, κP = material_opac(m, dn[k], tn[k])
            W.chifn[k] = χ; W.xken[k] = κE; W.xkpn[k] = κP
            t2 = tn[k] * tn[k]
            W.plfn[k] = P.csigr * (t2 * t2) / P.cpi
            W.egrn[k] = W.egn[k] + ern[k] / dn[k]
        end
    end
    true
end
"""
MATTHREADS pre-pass of the MLT loop: ∇_ad, c_P, δ and κ_R at the interface state (the means of the two zones, as in the
loop) of every interface with an MLT weight, on `nt` tasks. A zone whose evaluation throws is left unmarked and the serial
loop evaluates it itself, in order: the loop then behaves as without the pre-pass, exception included.
"""
function _mlt_pre!(m, W, cvw, mok, dn, tn, ngrs, ngre, nt)
    fill!(mok, false); next = Threads.Atomic{Int}(ngrs + 1)
    ZTASK[] = true
    try
        @sync for it in 1:nt
            Threads.@spawn begin
                while true
                    k = Threads.atomic_add!(next, 1)
                    k > ngre && break
                    @inbounds cvw[k] > 0 || continue
                    try
                        @inbounds begin
                            Ti = 0.5 * (tn[k-1] + tn[k]); ρi = 0.5 * (dn[k-1] + dn[k])
                            ZONE[] = k
                            nad, cpk, δk = material_thermo(m, ρi, Ti)
                            W.mnad[k] = nad; W.mcp[k] = cpk; W.mdl[k] = δk; W.mkp[k] = material_opac(m, ρi, Ti)[1]
                            mok[k] = true
                        end
                    catch
                    end
                end
            end
        end
    finally
        ZTASK[] = false
    end
    nothing
end
const SWEEP_WARMED = Set{DataType}()
const SWEEP_KERR = Int[]
const SWEEP_ERRS = Any[]

function _material!(m, W, c, P, dn, tn, ern, k0, ngre)
    hassrc = has_source(m)
    nzk = ngre + 3 - k0
    nt = min(MAT_THREADS[], Threads.nthreads(), nzk)
    (nt > 1 && nzk >= 2 && !ZTASK[]) && return _material_mt!(m, W, c, P, dn, tn, ern, k0, ngre, hassrc, nt)
    for k in k0:ngre+2
        _matzone!(m, W, c, P, dn, tn, ern, k, hassrc)
    end
    nothing
end
"""
The sweep of `_material!` on `nt` tasks, which take zones one at a time from a shared counter (the cost per zone is very uneven: stellar
EOS + network + opacity table in the hot core, a few table lookups outside). An exception in a zone is rethrown as in the serial sweep:
the one of the lowest zone, unwrapped, after all tasks have ended (the Newton iteration's `catch` then sees what it always saw).
"""
function _material_mt!(m, W, c, P, dn, tn, ern, k0, ngre, hassrc, nt)
    kend = ngre + 2; next = Threads.Atomic{Int}(k0)
    kerr = SWEEP_KERR; errs = SWEEP_ERRS
    length(kerr) < nt && (resize!(kerr, nt); resize!(errs, nt))
    @inbounds for it in 1:nt; kerr[it] = typemax(Int); errs[it] = nothing; end
    ZTASK[] = true
    try
        @sync for it in 1:nt
            Threads.@spawn begin
                while true
                    k = Threads.atomic_add!(next, 1)
                    k > kend && break
                    try
                        _matzone!(m, W, c, P, dn, tn, ern, k, hassrc)
                    catch e
                        k < kerr[it] && (kerr[it] = k; errs[it] = e)
                    end
                end
            end
        end
    finally
        ZTASK[] = false
    end
    kmin = typemax(Int); imin = 1
    @inbounds for it in 1:nt
        kerr[it] < kmin && (kmin = kerr[it]; imin = it)
    end
    if kmin < typemax(Int)
        ZONE0[] = kmin
        throw(errs[imin])
    end
    kl = kend; @inbounds while kl >= k0 && !(ForwardDiff.value(dn[kl]) > 0 && ForwardDiff.value(tn[kl]) > 0); kl -= 1; end
    kl >= k0 && (ZONE0[] = kl)                # what the serial sweep leaves behind
    nothing
end
"""
PERF512: `zone_sweep!(f, ks, ok)` calls `f(k)` for every zone `k in ks` with `ok(k)` true, on MATTHREADS tasks when MATTHREADS > 1
(zones taken one at a time from a shared counter), serially otherwise or when already inside a threaded sweep. For the Float64
EOS/opacity fills of the state arrays (`TitanMaterials.fill_eos!/fill_opac!` of a material whose zones are independent: each `f(k)` writes
only entries k). Bit-identical to the serial loop: the same values, ZONE[] left at the last zone of the serial order, and an exception
rethrown as the serial loop would have (the one of the lowest zone, after all tasks ended; the serial loop would not have evaluated the
zones above it, whose entries are then written here — the callers' state is reset after a failed evaluation).
"""
# (FIX-RACE) ZONESWEEP=0: every zone_sweep! serial (the PERF512 threaded loops off; MATTHREADS keeps the residual's material sweep
# threaded as before PERF512). Bit-identical either way.
const ZONESWEEP = Ref(get(ENV, "ZONESWEEP", "1") != "0")
function zone_sweep!(f::F, ks::UnitRange{Int}, ok::G) where {F,G}
    nt = ZONESWEEP[] ? min(MAT_THREADS[], Threads.nthreads(), length(ks)) : 1
    if !(nt > 1 && length(ks) >= 2 && !ZTASK[])
        for k in ks
            ok(k) || continue
            ZONE[] = k; f(k)
        end
        return nothing
    end
    # (FIX-RACE) the first zone serially once per closure type F, before any task starts: lazily built shared tables (photo-rate tables,
    # compiled networks, ...) that every zone needs are then built by one task, not raced for. Subsequent sweeps dispatch all zones to tasks.
    k1 = first(ks); while k1 <= last(ks) && !ok(k1); k1 += 1; end
    k1 > last(ks) && return nothing
    if !(F in SWEEP_WARMED)
        ZONE[] = k1; f(k1)
        push!(SWEEP_WARMED, F)
        kstart = k1 + 1
    else
        kstart = k1
    end
    next = Threads.Atomic{Int}(kstart); kend = last(ks)
    kerr = SWEEP_KERR; errs = SWEEP_ERRS
    length(kerr) < nt && (resize!(kerr, nt); resize!(errs, nt))
    @inbounds for it in 1:nt; kerr[it] = typemax(Int); errs[it] = nothing; end
    ZTASK[] = true
    try
        @sync for it in 1:nt
            Threads.@spawn begin
                while true
                    k = Threads.atomic_add!(next, 1)
                    k > kend && break
                    ok(k) || continue
                    try
                        ZONE[] = k; f(k)
                    catch e
                        k < kerr[it] && (kerr[it] = k; errs[it] = e)
                    end
                end
            end
        end
    finally
        ZTASK[] = false
    end
    kmin = typemax(Int); imin = 1
    @inbounds for it in 1:nt
        kerr[it] < kmin && (kmin = kerr[it]; imin = it)
    end
    if kmin < typemax(Int)
        ZONE0[] = kmin
        throw(errs[imin])
    end
    kl = kend; while kl >= first(ks) && !ok(kl); kl -= 1; end
    kl >= first(ks) && (ZONE0[] = kl)
    nothing
end

# ------------------------------------------------------------------ residual

"Grid-equation rows of GRIDEQ (see there) for the interfaces k1+1 … k2−1."
function grideq_rows!(R, ir, k0, rn, ro, dn, tn, pg, er, xmen, xmtot, dt, ge, chi = nothing; un = nothing, G4 = 0.0, cs = nothing, Gc = 0.0)
    k1, k2 = ge.k1, ge.k2; zl, zh = k1, k2 - 1                  # zones between the fixed interfaces
    α = ge.alph * (1 + ge.alph); wρ, wP, wT, wE, wm = ge.w
    zc(j) = clamp(j, zl, zh)
    ν(j) = (j = zc(j); (rn[j+1] + rn[j]) / (2 * (rn[j+1] - rn[j]))); νo(j) = (j = zc(j); (ro[j+1] + ro[j]) / (2 * (ro[j+1] - ro[j])))
    lm(j) = (j = zc(j); log(0.5 * ((xmtot - xmen[j]) + (xmtot - xmen[j+1]))))
    # Relaxation time of zone j from the Newton iterate (smooth, differentiated with the rows):
    #   tau = smax( smin(kappa t_cross, f_dyn t_dyn), t_floor ),  smin/smax(a,b) = (a^-p + b^-p)^(-1/p) / (a^p + b^p)^(1/p), p = 4.
    # t_cross = r_c / (|u|_c + c_s) (the driver's geq_tcross: r_c mid-zone radius, |u|_c mean of the two interface speeds, |u| -> sqrt(u^2 + eps^2),
    # eps = 1e-3 c_s; c_s = sqrt(gamma p/rho) = W.asn); t_dyn = 1/sqrt(G rho) (G = P.cgrav); t_floor: the speed floor of the driver (start-of-step).
    function tau(j)
        (cs === nothing || !hasproperty(ge, :tfl)) && return (ge.tauv[j], zero(ge.tauv[j]))
        c_ = cs[j]; ε = 1e-3 * c_
        um = 0.5 * (sqrt(un[j]^2 + ε^2) + sqrt(un[j+1]^2 + ε^2))
        tc_ = 0.5 * (rn[j] + rn[j+1]) / (um + c_)
        a_ = ge.kap * tc_; b_ = ge.fdyn / sqrt(Gc * dn[j])
        tm = (a_^-4 + b_^-4)^-0.25
        tf_ = log(rn[j+1] / rn[j]) * tc_ * ge.flc          # speed floor from the iterate: a zone moving its own width ln(r_j+1/r_j) in tau goes at most fv (|u|+c_s)  (ge.flc = 1/min(fmove-ish cap, fvr))
        (tm, tf_)
    end
    # Optical depth of the Newton iterate (GEQW entry 6, w_tau): tau_k = sum_{j>=k} chi_j rho_j (r_j+1 - r_j) from the outer zone zh inward (zones k1..zh),
    # ln tau confined to [ln tau_lo, ln tau_hi] by a SMOOTH window (smooth max/min, sqrt((x-L)^2 + e^2), e = 0.5), zone value = mean of the two interface values,
    # term w_tau h ((dlt_j+1 - dlt_j)^2 + (dlt_j - dlt_j-1)^2) in the same adjacent-interface stencil as the other monitor terms. Local to the iterate and differentiated.
    wt_ = hasproperty(ge, :wtau) && ge.wtau > 0 && chi !== nothing
    ltc = if wt_
        ta = Vector{typeof(rn[zl] * chi[zl] * dn[zl])}(undef, zh - zl + 2); ta[end] = zero(ta[end])
        for k in zh:-1:zl; ta[k-zl+1] = ta[k-zl+2] + chi[k] * dn[k] * (rn[k+1] - rn[k]); end
        llo, lhi = log(ge.tlo), log(ge.thi); e2 = 0.25
        smax_(x, L) = 0.5 * (x + L + sqrt((x - L)^2 + e2)); smin_(x, L) = 0.5 * (x + L - sqrt((x - L)^2 + e2))
        lw = [smin_(smax_(log(t + 1e-30), llo), lhi) for t in ta]
        [0.5 * (lw[j-zl+1] + lw[j-zl+2]) for j in zl:zh]
    else
        nothing
    end
    # monitors of zone j (physics M_p = sqrt(1+G2 capped), floor M_f), the smoothed node density xnt and its old value, all from the iterate
    function mons(j)
        xnt = ν(j) - α * (ν(j-1) - 2ν(j) + ν(j+1)); xnto = νo(j) - α * (νo(j-1) - 2νo(j) + νo(j+1))
        jm, jp = zc(j - 1), zc(j + 1); h = jp - jm == 2 ? 0.5 : 1.0
        d2(v) = h * (log(v[jp] / v[j])^2 + log(v[j] / v[jm])^2)
        S = wρ * d2(dn) + wP * d2(pg) + wT * d2(tn) +
            wE * d2(er) + wm * h * ((lm(jp) - lm(j))^2 + (lm(j) - lm(jm))^2) + ge.sx[j]
        wt_ && (S += ge.wtau * h * ((ltc[jp-zl+1] - ltc[j-zl+1])^2 + (ltc[j-zl+1] - ltc[jm-zl+1])^2))
        # (GEQW entry 8, wχ > 0: the extinction coefficient χρ of the Newton iterate, adjacent interface jumps: local, not lagged)
        (hasproperty(ge, :wχ) && ge.wχ > 0 && chi !== nothing) && (S += ge.wχ * h * (log((chi[jp] * dn[jp]) / (chi[j] * dn[j]))^2 + log((chi[j] * dn[j]) / (chi[jm] * dn[jm]))^2))
        # (GEQGP = w_g, L_ref; provisional, default off; docs/POPIII_ENERGY.md §14.6) the density jump weighted with the gravity power of the
        # zone, w_g (∂ln ρ/∂j)² ½(|u_j| g_j + |u_j+1| g_j+1) m_j / L_ref, of the Newton iterate (local, not lagged): where the advected interface
        # density differs from the mean density of the node AND gravity does work, i.e. where the two discrete gravity powers differ (G5)
        if hasproperty(ge, :gp) && ge.gp[1] > 0 && un !== nothing
            gq(i) = G4 * (xmtot - xmen[i] + (EXTMASS[] === nothing ? zero(rn[i]) : EXTMASS[](rn[i]))) / rn[i]^2
            S += ge.gp[1] * d2(dn) * 0.5 * (abs(un[j]) * gq(j) + abs(un[j+1]) * gq(j + 1)) * (xmen[j] - xmen[j+1]) / ge.gp[2]
        end
        G2 = ν(j)^2 * S
        G2 = G2 * ge.gmax^2 / (G2 + ge.gmax^2)
        # floor monitor (GEQDEC): M_f^2 = ω_r² + ω_m² ν² (Δ ln m)² (enclosed-mass term; adjacent-interface stencil, clamped to the shells from the first interface outward)
        Gf2 = ge.fr + ge.fm * ν(j)^2 * h * ((lm(jp) - lm(j))^2 + (lm(j) - lm(jm))^2)
        (sqrt(1 + G2), sqrt(Gf2), xnt, xnto)
    end
    # Mixture of separately normalised monitors (GEQDEC entry 3 = a): the row of zone j is x̂_j/M_j with
    #   M_j = a (I_p/I_f) M_f,j + (1-a) M_p,j,   I_x = Σ_j M_x,j Δξ_j,  Δξ_j = 1/ν_j (ν = (r_j+1 + r_j)/(2(r_j+1 - r_j)) ≈ 1/Δln r_j),
    # over the zones zl..zh; = I_p [a M_f/I_f + (1-a) M_p/I_p]; a = 0 gives the plain physics row. Equidistribution (the rows) makes M_j Δξ_j the same for all zones;
    # I_p, I_f are global sums of the iterate (Dual-safe, dense dependence only through J·v).
    ms = [mons(j) for j in zl:zh]
    Ipp = sum(m[1] / ν(j) for (j, m) in zip(zl:zh, ms)); Iff = sum(m[2] / ν(j) for (j, m) in zip(zl:zh, ms))
    amix = hasproperty(ge, :amix) ? ge.amix : 0.0
    # Inside-out ceiling on the relaxation time (smooth, differentiated, no lag): tau~_j = smin(tau_j, q tau~_{j-1}), tau~_zl = tau_zl, q = 1 + 5/N (5 % per zone).
    # A zone's own tau (t_dyn, t_cross, speed floor) is 1e3-1e8 dt once the step is set by the collapsing core. In the frozen limit the chain n̂_j = n̂_j-1 shares a change
    # of the span ln(r_N/r_zl) (r_zl collapses with the Lagrangian first interface) among the zones in proportion to M_j/tau_j, so the core zones alone absorb all of it:
    # dlnr_2..5 grow 0.3 -> 0.9, mass ratios 3 (n_c 1e16) -> 9. A nearly uniform tau (ceiling) shares it over all N zones (~0.03 e-fold each).
    τraw = [tau(j) for j in zl:zh]; τce = similar(first.(τraw)); τct = similar(τce); qτ = 1 + 5 / length(τraw)
    @inbounds for i in eachindex(τraw)
        τce[i] = i == 1 ? τraw[i][1] : (τraw[i][1]^-4 + (qτ * τce[i-1])^-4)^-0.25
        τct[i] = (τce[i]^4 + τraw[i][2]^4)^0.25        # speed floor AFTER the ceiling (smooth smax(tau~_j, t_floor,j)): no zone is slower than fv (|u|+c_s) allows; the ceiling chain itself is floor-free
    end
    nhat(j) = begin
        m = ms[j-zl+1]; Mj = amix > 0 ? amix * (Ipp / Iff) * m[2] + (1 - amix) * m[1] : m[1]
        γ = τct[j-zl+1] / dt
        ((1 + γ) * m[3] - γ * m[4]) / Mj
    end
    @inbounds for k in k1+1:k2-1
        R[ir, k - k0 + 1] = nhat(k) - nhat(k - 1)
        if GEQDBG[] && !isfinite(ForwardDiff.value(R[ir, k - k0 + 1]))     # (GEQDBG=1: report the first non-finite grid row)
            GEQDBG[] = false; j = k
            println("GEQDBG: grid row $k non-finite: r ", ForwardDiff.value.(rn[j-2:j+2]), " rho ", ForwardDiff.value.(dn[j-2:j+1]), " T ", ForwardDiff.value.(tn[j-2:j+1]),
                    " P ", ForwardDiff.value.(pg[j-2:j+1]), " E ", ForwardDiff.value.(er[j-2:j+1]), " chi ", chi === nothing ? nothing : ForwardDiff.value.(chi[j-2:j+1]),
                    " sx ", ge.sx[j-1:j], " tauv ", ge.tauv[j-1:j], " dt ", dt)
        end
    end
    nothing
end
const GEQDBG = Ref(get(ENV, "GEQDBG", "0") in ("1", "2"))
const GEQDBG2 = Ref(get(ENV, "GEQDBG", "0") == "2" ? 5 : 0)

const SLState{T} = NamedTuple{(:q, :t, :fs, :g, :ts, :wq, :dp), NTuple{7, Vector{T}}}
const AD_JV = Ref(false)
const EDDFAC_DUAL = Ref{Any}(nothing)

function residual!(R::AbstractMatrix{T}, ws::WS{T}, c::Ctx, x::AbstractVector{T}) where {T}
    ngrs, ngre, neqn, mu, k0, nb = c.ngrs, c.ngre, c.neqn, c.mu, c.k0, c.nb
    P = c.p; θ = P.θ; dt = P.dt; xmu = Float64(mu); xmup1 = xmu + 1; cc = P.cc
    ir, im, id, iu, it, ie, iF = c.ix
    W = ws.a
    rn, xmen, dn, un, tn, ern, frn = W.rn, W.xmen, W.dn, W.un, W.tn, W.ern, W.frn
    newv = (rn, xmen, dn, un, tn, ern, frn)
    @inbounds for cde in 1:NV, k in k0:ngre+2
        newv[cde][k] = c.base[k, cde]
    end
    @inbounds for j in 1:neqn
        cde = c.code[j]; cde == 0 && continue; a = newv[cde]
        if cde == 4
            for bb in 1:nb; k = k0 + bb - 1; a[k] = c.base[k, 4] + x[(bb-1)*neqn+j] * c.nomu[k]; end
        elseif cde == 7
            for bb in 1:nb; k = k0 + bb - 1; a[k] = c.base[k, 7] + x[(bb-1)*neqn+j] * c.nomf[k]; end
        else
            for bb in 1:nb; k = k0 + bb - 1; a[k] = c.base[k, cde] * exp(x[(bb-1)*neqn+j]); end
        end
    end
    if AD_JV[]
        mup1 = mu + 1; rxm1 = 1.0 / xmup1
        k = ngrs - 1
        if P.leibc == 1
            tn[k] = tn[k+1]; dn[k] = dn[k+1]
            rn[k] = rn[k+1] * 2.0 - rn[k+2]; xmen[k] = xmen[k+1] * 2.0 - xmen[k+2]; un[k] = -un[k+2]
        elseif P.leibc == 2
            tn[k] = T(gf(:tl)); dn[k] = T(P.dl); rn[k] = rn[k+1] - P.delrl
            xmen[k] = xmen[k+1] + P.dl * rxm1 * (rn[k+1]^mup1 - rn[k]^mup1); un[k] = T(gf(:ul))
        end
        if P.llibc == 1
            tn[k] = tn[k+1]; dn[k] = dn[k+1]
            rn[k] = (rn[k+1]^mup1 - xmup1 * P.delml / dn[k])^rxm1
            xmen[k] = xmen[k+1] + P.delml; un[k] = un[k+1]
        end
        if c.lrad
            ern[k] = ern[k+1]
            P.lribc == 1 && (frn[k] = +frn[k+1])
            P.lribc == 2 && (frn[k] = -frn[k+2])
            P.lribc == 3 && (frn[k] = frn[k+1] * (rn[k+1] / rn[k])^mu)
        end
        k = ngre + 1
        if P.leobc == 1
            tn[k] = tn[k-1]; tn[k+1] = tn[k]; dn[k] = dn[k-1]; dn[k+1] = dn[k]
            rn[k+1] = rn[k] * 2.0 - rn[k-1]; xmen[k+1] = xmen[k] * 2.0 - xmen[k-1]; un[k+1] = -un[k-1]
        elseif P.leobc == 2
            tr_ = T(gf(:tr)); dr_ = T(P.dr)
            tn[k] = tr_; tn[k+1] = tr_; dn[k] = dr_; dn[k+1] = dr_
            rn[k+1] = rn[k] + P.delrr; xmen[k+1] = xmen[k] - P.dr * rxm1 * (rn[k+1]^mup1 - rn[k]^mup1); un[k+1] = T(gf(:ur))
        elseif P.leobc == 3
            tn[k] = tn[k-1]; tn[k+1] = tn[k]; dn[k] = dn[k-1]; dn[k+1] = dn[k]
            rn[k+1] = rn[k] + P.delrr; xmen[k+1] = xmen[k] - dn[k] * rxm1 * (rn[k+1]^mup1 - rn[k]^mup1); un[k+1] = un[k]
        end
        if P.llobc > 0
            tn[k] = tn[k-1]; tn[k+1] = tn[k]; dn[k] = dn[k-1]; dn[k+1] = dn[k]
            rn[k+1] = (rn[k]^mup1 + xmup1 * P.delmr / dn[k])^rxm1; xmen[k+1] = xmen[k] - P.delmr; un[k+1] = un[k]
        end
        if c.lrad
            P.lrobc == 1 && (ern[k] = ern[k-1]; ern[k+1] = ern[k]; frn[k+1] = +frn[k])
            P.lrobc == 2 && (ern[k] = ern[k-1]; ern[k+1] = ern[k-2]; frn[k+1] = -frn[k-1])
            P.lrobc == 3 && (ern[k] = ern[k-1]; ern[k+1] = ern[k]; frn[k+1] = frn[k] * (rn[k] / rn[k+1])^mu)
        end
    end
    Mg = MG[] ? mgws(T, c.mgr) : nothing          # MGROUPS: the group unknowns of the iterate
    MG[] && mg_iterate!(Mg, x, neqn, nb, k0, ngre)
    # EOS / opacity (TITAN built-in physics unless a material is installed, see TitanMaterials)
    Rg = P.ck / (P.cm0 * P.gmu); γ = P.gam
    pgn, egn, asn, chifn, xken, xkpn, plfn, egrn = W.pgn, W.egn, W.asn, W.chifn, W.xken, W.xkpn, W.plfn, W.egrn
    if COUPLED[]
        xs1, xs2 = W.xs1, W.xs2; xn1, xn2 = XSN[]
        @inbounds for k in k0:ngre+2; xs1[k] = xn1[k]; xs2[k] = xn2[k]; end
        @inbounds for bb in 1:nb
            k = k0 + bb - 1; xs1[k] = xn1[k] * exp(x[(bb-1)*neqn+8]); xs2[k] = xn2[k] * exp(x[(bb-1)*neqn+9])
        end
        CPLW[] = (xs1, xs2)
        slst::Union{Nothing,SLState{T}} = if SLIMP[] # SLIMP: the star-light state of the iterate (new-level r, ρ, T, u, E_r, F_r, species)
            s_ = get!(() -> (; q = Vector{T}(undef, length(xs1)), t = Vector{T}(undef, length(xs1)), fs = Vector{T}(undef, length(xs1)),
                               g = Vector{T}(undef, length(xs1)), ts = Vector{T}(undef, length(xs1)), wq = Vector{T}(undef, length(xs1)), dp = Vector{T}(undef, length(xs1))), SLQBUF, (T, length(xs1)))::SLState{T}
            @inbounds for k in eachindex(s_.q); s_.q[k] = one(T); s_.t[k] = one(T); s_.fs[k] = zero(T); s_.g[k] = one(T); s_.ts[k] = zero(T); s_.wq[k] = zero(T); s_.dp[k] = zero(T); end
            SL_SWEEP[](s_, ngrs, ngre, rn, dn, tn, un, ern, frn, xs1, xs2, cplz;
                         mgqt = Mg === nothing ? nothing : (k -> (gq = MGSLG[]; MGSLQ[](k, ern[k], Mg.egn[gq][k], Mg.fgn[gq][k], Mg.fgn[gq][k+1]))))
            if SLIMPDBG[] > 0
                fin(v) = isfinite(ForwardDiff.value(v)) && (!(v isa ForwardDiff.Dual) || all(isfinite, ForwardDiff.partials(v)))
                for k in ngrs:ngre
                    (fin(s_.q[k]) && fin(s_.t[k]) && fin(s_.fs[k]) && fin(s_.g[k])) && continue
                    SLIMPDBG[] -= 1; println("SLIMPDBG star-light state non-finite zone $(k - ngrs + 1): q $(s_.q[k]) t $(s_.t[k]) fs $(s_.fs[k]) g $(s_.g[k]) E $(ern[k]) F $(frn[k]) $(frn[k+1])"); break
                end
            end
            CPLSL[] = s_
            s_
        else
            nothing
        end
        try
            MATERIAL[] === nothing || _material!(MATERIAL[], W, c, P, dn, tn, ern, k0, ngre)
        finally
            CPLW[] = nothing; CPLSL[] = nothing
        end
        # (PERF512) the coupled zones' rates in the Jacobian's Dual numbers on MATTHREADS tasks (zones independent; the material sweep above
        # left their values and local Jacobians in the zones' cache entries); bit-identical. ZONE[] is not used here and is kept.
        let z0 = ZONE0[]
            try zone_sweep!(ngrs:ngre, k -> true) do k
                @inbounds if cplz(k)
                    y = slst === nothing ? CPL_ZONE[](k, dn[k], tn[k], xs1[k], xs2[k]) :
                        CPL_ZONEQ[](k, dn[k], tn[k], xs1[k], xs2[k], slst.q[k], slst.t[k], slst.fs[k], slst.g[k], slst.dp[k])
                    W.cpS1[k] = y[3]; W.cpS2[k] = y[4]; W.cpL[k] = y[5]
                    if SLIMPDBG[] > 0 && !all(v -> isfinite(ForwardDiff.value(v)) && (!(v isa ForwardDiff.Dual) || all(isfinite, ForwardDiff.partials(v))), y)
                        SLIMPDBG[] -= 1; println("SLIMPDBG cpl zone $(k - ngrs + 1) non-finite: y ", y, " state ", slst === nothing ? 1 : (slst.q[k], slst.t[k], slst.fs[k], slst.g[k]))
                    end
                    W.cpQ[k] = y[6]; W.cpY[k] = y[7]
                else
                    W.cpS1[k] = zero(T); W.cpS2[k] = zero(T); W.cpL[k] = zero(T); W.cpQ[k] = zero(T); W.cpY[k] = zero(T)
                end
            end
            finally ZONE0[] = z0 end
        end
    else
        MATERIAL[] === nothing || _material!(MATERIAL[], W, c, P, dn, tn, ern, k0, ngre)
    end
    MATERIAL[] === nothing && @inbounds for k in k0:ngre+2
        pgn[k] = Rg * dn[k] * tn[k]; egn[k] = Rg * tn[k] / (γ - 1); asn[k] = sqrt(γ * pgn[k] / dn[k])
        if c.lrad
            chifn[k] = P.lopac == 3 ? P.chif0 / dn[k] : T(P.csige / P.cm0)
            xken[k] = P.ratio * chifn[k]; xkpn[k] = P.ratio * chifn[k]
            t2 = tn[k] * tn[k]
            plfn[k] = P.csigr * (t2 * t2) / P.cpi
            egrn[k] = egn[k] + ern[k] / dn[k]
        end
    end
    # time-centred
    r, u, fr, xm, dmdt, urel, rmu, rmup1 = W.r, W.u, W.fr, W.xm, W.dmdt, W.urel, W.rmu, W.rmup1
    @inbounds for k in ngrs-1:ngre+2
        r[k] = θ * rn[k] + (1 - θ) * c.ro[k]
        u[k] = θ * un[k] + (1 - θ) * c.uo[k]
        fr[k] = θ * frn[k] + (1 - θ) * c.fro[k]
        xm[k] = isempty(XMINT[]) ? θ * (P.xmtot - xmen[k]) + (1 - θ) * c.xmo[k] : T(XMINT[][k])
        EXTMASS[] === nothing || (xm[k] += EXTMASS[](r[k]))
        dmdt[k] = -(xmen[k] - c.xmeo[k]) / dt
        urel[k] = u[k] - (rn[k] - c.ro[k]) / dt
        rmu[k] = r[k]^mu; rmup1[k] = r[k]^(mu + 1)
    end
    d, pg, eg, as, er, egr, chif, xke, xkp, plf, dvol, dvoln =
        W.d, W.pg, W.eg, W.as, W.er, W.egr, W.chif, W.xke, W.xkp, W.plf, W.dvol, W.dvoln
    @inbounds for k in ngrs-1:ngre+1
        d[k] = θ * dn[k] + (1 - θ) * c.do_[k]
        pg[k] = θ * pgn[k] + (1 - θ) * c.pgo[k]
        eg[k] = θ * egn[k] + (1 - θ) * c.ego[k]
        as[k] = θ * asn[k] + (1 - θ) * c.aso[k]
        if c.lrad
            er[k] = θ * ern[k] + (1 - θ) * c.ero[k]
            egr[k] = θ * egrn[k] + (1 - θ) * c.egro[k]
            chif[k] = θ * chifn[k] + (1 - θ) * c.chifo[k]
            xke[k] = θ * xken[k] + (1 - θ) * c.xkeo[k]
            xkp[k] = θ * xkpn[k] + (1 - θ) * c.xkpo[k]
            plf[k] = θ * plfn[k] + (1 - θ) * c.plfo[k]
        end
        dvol[k] = (rmup1[k+1] - rmup1[k]) / xmup1
        dvoln[k] = (rn[k+1]^(mu + 1) - rn[k]^(mu + 1)) / xmup1
    end
    avchi = W.avchi
    c.lrad && @inbounds for k in ngrs:ngre+1; avchi[k] = 2 / (1 / chif[k-1] + 1 / chif[k]); end
    c.lrad && SLABFLUX[] && sf_avchi!(avchi, chif, ngrs, ngre)
    c.lrad && SUBZONE[] && sz_avchi!(avchi, chif, d, dvol, rmu, ngrs, ngre)
    # viscosity
    qf, dudr = W.qf, W.dudr; xmu4 = 0.25xmu
    @inbounds for k in ngrs:ngre
        dudr[k] = (u[k+1] - u[k]) / (r[k+1] - r[k]) - xmu4 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
        ql = AV_ZONE[] > 0 ? AV_ZONE[] * (r[k+1] - r[k]) : P.ql0 + P.ql1 * 0.5 * (r[k] + r[k+1])
        div = (rmu[k+1] * u[k+1] - rmu[k] * u[k]) / dvol[k]
        if AV_TFRAC[] > 0 || AV_TRAD[]            # reduced tension (see AV_TFRAC)
            rc = val(div) < 0 && val(dudr[k]) < 0
            ten = val(div) < 0 && !rc && (!AV_TRAD[] || val(u[k+1]) < val(u[k]))
            qv = decide(:vis, 0, k, rc) ? div : decide(:vis, 1, k, ten) ? (AV_TFRAC[] > 0 ? AV_TFRAC[] : 1.0) * div : zero(T)
        else
            qv = decide(:vis, 0, k, val(div) < 0 && (AV_TENSION[] || val(dudr[k]) < 0)) ? div : zero(T)
        end
        qm = ql * (P.cq1 * as[k] - P.cq2 * ql * qv)
        qf[k] = -P.cqvis * d[k] * qm * dudr[k]
    end
    # advection
    db, dsn, fl, ub, usn, flu, fle, egb, egsn = W.db, W.dsn, W.fl, W.ub, W.usn, W.flu, W.fle, W.egb, W.egsn
    @inbounds for k in ngrs+1:ngre+1; fl[k] = urel[k]; end
    advect!(db, dsn, dn, c.do_, oldslope!(W.dso, c.do_, P.cadv, P.eps, ngrs, ngre), fl, θ, P.cadv, P.eps, ngrs, ngre, true, 1)
    @inbounds for k in ngrs+1:ngre+1; fle[k] = -dmdt[k]; end
    @inbounds for k in ngrs:ngre; flu[k] = -(dmdt[k] + dmdt[k+1]); end
    advect!(ub, usn, un, c.uo, oldslope!(W.uso, c.uo, P.cadv, P.eps, ngrs, ngre), flu, θ, P.cadv, P.eps, ngrs, ngre, false, 2)
    if c.lrad
        advect!(egb, egsn, egrn, c.egro, oldslope!(W.egrso, c.egro, P.cadv, P.eps, ngrs, ngre), fle, θ, P.cadv, P.eps, ngrs, ngre, true, 3)
        GASROW[] && advect!(W.egbg, W.egsng, egn, c.ego, oldslope!(W.egso, c.ego, P.cadv, P.eps, ngrs, ngre), fle, θ, P.cadv, P.eps, ngrs, ngre, true, 6)
        qe, qeo, qF, qFo = W.qe, W.qeo, W.qF, W.qFo
        @inbounds for k in ngrs-1:ngre+2
            qe[k] = ern[k] / dn[k]; qeo[k] = c.ero[k] / c.do_[k]
            qF[k] = frn[k] * 2 / (dn[k] + dn[k-1]); qFo[k] = c.fro[k] * 2 / (c.do_[k] + c.do_[k-1])
        end
        advect!(W.erb, W.ersn, qe, qeo, oldslope!(W.erso, qeo, P.cadv, P.eps, ngrs, ngre), fle, θ, P.cadv, P.eps, ngrs, ngre, true, 4)
        advect!(W.frb, W.frsn, qF, qFo, oldslope!(W.frso, qFo, P.cadv, P.eps, ngrs, ngre), flu, θ, P.cadv, P.eps, ngrs, ngre, false, 5)
    else
        advect!(egb, egsn, egn, c.ego, oldslope!(W.egso, c.ego, P.cadv, P.eps, ngrs, ngre), fle, θ, P.cadv, P.eps, ngrs, ngre, true, 6)
    end
    erb, frb = W.erb, W.frb
    MG[] && mg_theta!(Mg, c, W, θ, ngrs, ngre, fle, flu)     # MGROUPS: θ-level group opacities, emission, χ̄_g, face values
    if COUPLED[]
        advect!(W.xb1, W.xsn1, W.xs1, XSO[][1], oldslope!(W.xso1, XSO[][1], P.cadv, P.eps, ngrs, ngre), fle, θ, P.cadv, P.eps, ngrs, ngre, true, 7)
        advect!(W.xb2, W.xsn2, W.xs2, XSO[][2], oldslope!(W.xso2, XSO[][2], P.cadv, P.eps, ngrs, ngre), fle, θ, P.cadv, P.eps, ngrs, ngre, true, 8)
    end
    # diffusion
    dfm, dfe = W.dfm, W.dfe
    @inbounds for k in ngrs:ngre+1
        dfm[k] = P.sigd == 0 ? zero(T) : 2P.sigd * rmu[k] * (d[k] - d[k-1]) / (r[k+1] - r[k-1])
        dfe[k] = P.sige == 0 ? zero(T) : P.sige * rmu[k] * (d[k] + d[k-1]) * (eg[k] - eg[k-1]) / (r[k+1] - r[k-1])
    end
    fedd = c.fedd; HEl = TD.HEDD; HXl = TD.HXEDD; FXl = TD.FXEDD; TWl = TD.TIEW; gedr = P.geddr
    if AD_JV[] && c.lrad && T <: ForwardDiff.Dual && EDDFAC_DUAL[] !== nothing
        fedd, HEl, HXl, FXl, TWl, gedr = EDDFAC_DUAL[](c, P, ngrs, ngre, chifn, xkpn, plfn, dn, tn, ern, frn, rn, W.xs1, W.xs2, T, DBUF, SLQBUF, MATERIAL[], SLIMP[], cplz)
    elseif TD.CLOSJ[] && c.lrad && TD.CSENS.valid && TD.CSENS.a == ngrs && T <: ForwardDiff.Dual
        fedd, HEl, HXl, FXl, TWl, gedr = closure_lift(c, P, ngrs, ngre, chifn, xkpn, plfn, dn, tn, ern, rn, W.xs1, W.xs2, T)
    end
    fill!(R, zero(T))
    # ---------------------------------------------------------------- GRID
    kb(k) = k - k0 + 1
    X(i, k) = AD_JV[] ? zero(T) : x[(kb(k) - 1) * neqn + i]
    if P.lgrid == 1
        α = P.alph * (1 + P.alph); gam = (P.tau / dt)^P.ibet
        xnu, xnc, ss, xnuo = W.xnu, W.xnc, W.ss, ws.xnuo
        @inbounds for k in ngrs:ngre
            if P.ladx == 1
                xnu[k] = P.xscale / (rn[k+1] - rn[k]); xnuo[k] = P.xscale / (c.ro[k+1] - c.ro[k])
            else
                xnu[k] = (rn[k+1] + rn[k]) / (rn[k+1] - rn[k]) / 2
                xnuo[k] = (c.ro[k+1] + c.ro[k]) / (c.ro[k+1] - c.ro[k]) / 2
            end
        end
        @inbounds for k in ngrs+1:ngre-1
            xnt = xnu[k] - α * (xnu[k-1] - 2xnu[k] + xnu[k+1])
            xnto = xnuo[k] - α * (xnuo[k-1] - 2xnuo[k] + xnuo[k+1])
            xnc[k] = (1 + gam) * xnt - gam * xnto
        end
        @inbounds for k in ngrs:ngre-1; ss[k] = zero(T); end
        for l in 1:11
            (c.wt[l] > 0 && c.lady[l] > 0) || continue
            l == 8 && error("lady(8) not ported")
            @inbounds for k in ngrs:ngre-1
                cs = if l == 1
                    xs = c.yscl[1]
                    c.lady[1] == 1 ? (xmen[k+1] - xmen[k]) / xs :
                    c.lady[1] == 2 ? (xmen[k+1] - xmen[k]) / (2xs - (xmen[k+1] + xmen[k])) :
                    (xs - xmen[k+1]) / (xs - xmen[k]) - (xs - xmen[k]) / (xs - xmen[k+1])
                elseif l == 11; rescs(xmen, k, c.lady[11], c.yscl[11])
                elseif l == 2;  rescs(dn, k, c.lady[2], c.yscl[2])
                elseif l == 3;  rescs(tn, k, c.lady[3], c.yscl[3])
                elseif l == 4;  rescs(ern, k, c.lady[4], c.yscl[4])
                elseif l == 5;  rescs(pgn, k, c.lady[5], c.yscl[5])
                elseif l == 6;  rescs(egn, k, c.lady[6], c.yscl[6])
                elseif l == 7;  rescs(chifn, k, c.lady[7], c.yscl[7])
                elseif l == 9;  rescs(frn, k, c.lady[9], c.yscl[9])
                else            T(rescs(c.bri, k, c.lady[10], c.yscl[10]))
                end
                ss[k] += c.wt[l] * cs^2
            end
        end
        R[ir, kb(ngrs+1)] = xnu[ngrs+1] - xnu[ngrs]
        rrm = sqrt(1 + xnu[ngrs+1]^2 * ss[ngrs+1])
        @inbounds for k in ngrs+2:ngre-1
            rrk = sqrt(1 + xnu[k]^2 * ss[k])
            R[ir, kb(k)] = xnc[k] / rrk - xnc[k-1] / rrm
            rrm = rrk
        end
        R[ir, kb(ngre)] = xnu[ngre] - xnu[ngre-1]
        P.leibc > 0 && (R[ir, kb(ngrs)] = rn[ngrs] - c.ro[ngrs])
        P.llibc == 1 && (R[ir, kb(ngrs)] = (rn[ngrs] - c.ro[ngrs]) / dt - u[ngrs])
        P.leobc > 0 && (R[ir, kb(ngre+1)] = rn[ngre+1] - c.ro[ngre+1])
        P.llobc > 0 && (R[ir, kb(ngre+1)] = (rn[ngre+1] - c.ro[ngre+1]) / dt - u[ngre+1])
    elseif P.lgrid == 2
        gv = GRIDVEL[]
        if gv === nothing
            @inbounds for k in ngrs:ngre+1; R[ir, kb(k)] = rn[k] - c.ro[k]; end
        else
            @inbounds for k in ngrs:ngre+1; R[ir, kb(k)] = rn[k] - c.ro[k] - gv[k] * dt; end
        end
        GRIDEQ[] === nothing || grideq_rows!(R, ir, k0, rn, c.ro, dn, tn, W.pgn, ern, xmen, P.xmtot, dt, GRIDEQ[], W.chifn; un, G4 = 2 * Float64(c.mu) * P.cpi * P.cgrav, cs = W.asn, Gc = Float64(P.cgrav))
        GEQ_LAGFIRST[] && gridmoves() && (R[ir, kb(ngrs+1)] = (xmen[ngrs+1] - c.xmeo[ngrs+1]) * (c.ro[ngrs+1] / (c.xmeo[ngrs] - c.xmeo[ngrs+1])))
    else
        @inbounds for k in ngrs:ngre+1; R[ir, kb(k)] = (rn[k] - c.ro[k]) / dt - u[k]; end
    end
    REGULAR_CENTRE[] && (R[ir, kb(ngrs)] = X(ir, ngrs))
    k = ngrs - 1
    if P.leibc == 1; R[ir, kb(k)] = rn[k] - 2rn[k+1] + rn[k+2]
    elseif P.leibc == 2; R[ir, kb(k)] = rn[k] - rn[k+1] + P.delrl
    elseif P.llibc == 1; R[ir, kb(k)] = -P.delml / dn[k] + (rn[k+1]^(mu+1) - rn[k]^(mu+1)) / xmup1
    end
    k = ngre + 2
    if P.leobc == 1; R[ir, kb(k)] = rn[k] - 2rn[k-1] + rn[k-2]
    elseif P.leobc in (2, 3); R[ir, kb(k)] = rn[k] - rn[k-1] - P.delrr
    elseif P.llobc > 0; R[ir, kb(k)] = -P.delmr / dn[k-1] + (rn[k]^(mu+1) - rn[k-1]^(mu+1)) / xmup1
    end
    R[ir, 1] = X(ir, k0)
    # ---------------------------------------------------------------- MASS
    @inbounds for k in ngrs+1:ngre+1; R[im, kb(k)] = xmen[k-1] - xmen[k] - dn[k-1] * dvoln[k-1]; end
    # TITAN 1.2 (mass.f) has `+ rmu·phil0·dt`: the opposite sign to its own continuity row (contin.f: inflow rmu·phil0
    # into zone ngrs) and to peq (dmdt(ngrs) = −r^mu·phil0). Never exercised (phil0 = 0 in every test problem); with
    # mass leaving through the inner boundary it moves the exterior mass the wrong way, so m(r) loses what the core
    # should gain and every advective flux derived from dm_e/dt is off by 2r²·phil0. FIX_PHIL0[] = false restores 1.2.
    R[im, kb(ngrs)] = xmen[ngrs] - c.xmeo[ngrs] + (FIX_PHIL0[] ? -1 : 1) * rmu[ngrs] * P.phil0 * dt
    R[im, 1] = X(im, k0)
    # Lagrangian with fixed interior masses: x_me is frozen (the continuity row with u_rel = 0 keeps every shell's
    # mass). TITAN's row x_me(k−1) − x_me(k) = ρΔV takes the shell mass as a difference of exterior masses ~x_mtot,
    # which loses log10(x_mtot/Δm) digits in the density of the central shells (1e-8 in a static sphere).
    !isempty(XMINT[]) && @inbounds for k in ngrs:ngre+1; R[im, kb(k)] = X(im, k); end
    k = ngrs - 1
    if P.leibc == 1; R[im, kb(k)] = xmen[k] - 2xmen[k+1] + xmen[k+2]
    elseif P.leibc == 2; R[im, kb(k)] = -xmen[k] + xmen[k+1] + P.dl * dvoln[k]
    elseif P.llibc > 0; R[im, kb(k)] = -xmen[k] + xmen[k+1] + P.delml
    else R[im, kb(k)] = X(im, k)
    end
    k = ngre + 2
    if P.leobc == 1; R[im, kb(k)] = xmen[k-2] - 2xmen[k-1] + xmen[k]
    elseif P.leobc == 2; R[im, kb(k)] = xmen[k-1] - xmen[k] - P.dr * dvoln[k-1]
    elseif P.leobc == 3; R[im, kb(k)] = xmen[k-1] - xmen[k] - dn[k-1] * dvoln[k-1]
    elseif P.llobc > 0; R[im, kb(k)] = xmen[k-1] - xmen[k] + P.delmr
    else R[im, kb(k)] = X(im, k)
    end
    # ---------------------------------------------------------------- CONTIN
    @inbounds for k in ngrs:ngre
        fin  = k == ngrs ? rmu[k] * P.phil0 : rmu[k] * urel[k] * db[k]
        fout = k == ngre ? (P.leobc == 3 ? rmu[k+1] * P.phir0 + rmu[k+1] * urel[k+1] * db[k+1] :
                            rmu[k+1] * P.phir0) : rmu[k+1] * urel[k+1] * db[k+1]
        R[id, kb(k)] = (dn[k] * dvoln[k] - c.do_[k] * c.dvolo[k]) / dt + fout - fin - dfm[k+1] + dfm[k]
    end
    R[id, 1] = X(id, k0)
    k = ngrs - 1
    R[id, kb(k)] = (P.leibc == 1 || P.llibc > 0) ? dn[k+1] - dn[k] : X(id, k)
    for k in ngre+1:ngre+2
        R[id, kb(k)] = P.leobc == 1 ? -dn[k-1] + dn[k] : X(id, k)
    end
    # ---------------------------------------------------------------- GAS MOMENTUM
    @inbounds for k in ngrs+1:ngre
        mk = 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1])
        r1 = r[k]^(-0.5xmu)
        r3k = ((r[k] + r[k+1]) / 2)^(1.5xmu); r3m = ((r[k-1] + r[k]) / 2)^(1.5xmu)
        fgk = if GFV[] || GVAR[]                # gravity force with the finite-volume weight / as the gradient of the binding energy
            xek = gxe(r[k]); G4_ = 2xmu * P.cpi * P.cgrav
            gfz(k - 1, k, d, dvol, r, rmu, G4_, xm[k], xek) + gfz(k, k, d, dvol, r, rmu, G4_, xm[k], xek) + (1 - 0.5xmu) * P.g * mk
        else
            (2xmu * P.cpi * P.cgrav * xm[k] / rmu[k] + (1 - 0.5xmu) * P.g) * mk
        end
        w = (un[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
             c.uo[k] * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2dt) -
            0.5 * ((dmdt[k] + dmdt[k+1]) * ub[k] - (dmdt[k-1] + dmdt[k]) * ub[k-1]) +
            rmu[k] * (pg[k] - pg[k-1]) +
            fgk +
            r1 * (r3k * qf[k] - r3m * qf[k-1])
        c.lrad && (w -= (Mg === nothing ? avchi[k] * fr[k] : Mg.frc[k]) / cc * mk)   # (MGROUPS: Σ_g χ̄_g F_g)
        R[iu, kb(k)] = w
    end
    R[iu, kb(ngrs)] = X(iu, ngrs); R[iu, kb(ngre+1)] = X(iu, ngre + 1); R[iu, 1] = X(iu, k0)
    if !isnan(PEXT[])
        # free surface: the interior momentum row of interface ngre+1 with an empty zone outside (no mass, no
        # viscous stress, pressure PEXT[]): half the mass of zone ngre, pressure jump PEXT − p(ngre)
        (P.lgrid == 3 && P.llobc == 1) || error("PEXT needs lgrid = 3, llobc = 1")
        k = ngre + 1
        mk = 0.5 * d[k-1] * dvol[k-1]
        r3m = ((r[k-1] + r[k]) / 2)^(1.5xmu)
        w = (un[k] * dn[k-1] * dvoln[k-1] - c.uo[k] * c.do_[k-1] * c.dvolo[k-1]) / (2dt) +
            0.5 * (dmdt[k-1] + dmdt[k]) * ub[k-1] +
            rmu[k] * (PEXT[] - pg[k-1]) +
            (2xmu * P.cpi * P.cgrav * xm[k] / rmu[k] + (1 - 0.5xmu) * P.g) * mk -
            r[k]^(-0.5xmu) * r3m * qf[k-1]
        c.lrad && (w -= avchi[k] * fr[k] / cc * mk)
        R[iu, kb(k)] = w
    end
    (P.leibc in (1, 2) || P.llibc == 1) || error("gasmom inner BC")
    (P.leobc < 3 || P.llobc == 1) || error("gasmom outer BC")
    k = ngrs - 1
    R[iu, kb(k)] = P.leibc == 1 ? un[k+2] + un[k] : P.llibc == 1 ? un[k] - un[k+1] : X(iu, k)
    k = ngre + 2
    R[iu, kb(k)] = P.leobc == 1 ? un[k-2] + un[k] : P.leobc == 3 ? -un[k-1] + un[k] : X(iu, k)
    # ---------------------------------------------------------------- MLT convective flux (Phase 7)
    fcv = W.fcv
    @inbounds for k in ngrs:ngre+1; fcv[k] = zero(T); W.vcv[k] = zero(T); W.lcv[k] = zero(T); W.nabv[k] = zero(T); W.nadv[k] = zero(T); end
    if c.conv
        nadv, cpv, δv, κv = c.cv; cvw = BUFS[:cvw]::Vector{Float64}
        arad = 4P.csigr / cc
        ntm = min(MAT_THREADS[], Threads.nthreads())
        mpre = MLT_AD[] && ntm > 1 && ngre - ngrs >= 4ntm && !ZTASK[]
        mok = get!(() -> fill(false, length(cvw)), BUFS, :mltok)::Vector{Bool}
        mpre && _mlt_pre!(MATERIAL[], W, cvw, mok, dn, tn, ngrs, ngre, ntm)
        @inbounds for k in ngrs+1:ngre
            Pm = pgn[k-1] + (AD_JV[] ? fedd[k-1] : c.fedd[k-1]) * ern[k-1]; Pp = pgn[k] + (AD_JV[] ? fedd[k] : c.fedd[k]) * ern[k]
            wτ = cvw[k]; wτ > 0 || continue          # optically thin: no MLT flux
            Ti = 0.5 * (tn[k-1] + tn[k]); ρi = 0.5 * (dn[k-1] + dn[k]); Pi = 0.5 * (Pm + Pp)
            gk = 2xmu * P.cpi * P.cgrav * xm[k] / rmu[k] + (1 - 0.5xmu) * P.g
            ForwardDiff.value(gk) > 0 || continue          # no gravity → no buoyancy
            HP = Pi / (ρi * gk); l = MLT_RCAP[] ? c.alpha / sqrt(1 / HP^2 + 1 / r[k]^2) : c.alpha * HP
            # ∇ = Δln T/Δln P with the actual pressure difference. MLT_HS: where the layer is far from hydrostatic
            # (|Δln P| below half of the hydrostatic −Δr/H_P, Δr between zone centres) blend smoothly to the
            # hydrostatic value: behind the accretion shock the pressure is nearly flat, Δln P → 0 and ∇ reached 110.
            # Hydrostatic layers (ratio ≈ 1) are untouched: deep convection has ∇ − ∇_ad ~ 1e-6, so ∇ must keep the
            # discretisation of the pressure difference that balances the momentum equation.
            lp = log(Pp / Pm)
            if MLT_HS[]
                lph = MLT_HSDISC[] ? ((GFV[] || GVAR[]) ? -(gfz(k - 1, k, d, dvol, r, rmu, 2xmu * P.cpi * P.cgrav, xm[k], gxe(r[k])) +
                                                           gfz(k, k, d, dvol, r, rmu, 2xmu * P.cpi * P.cgrav, xm[k], gxe(r[k])) +
                                                           (1 - 0.5xmu) * P.g * 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1])) / (rmu[k] * Pi) :
                                      -gk * 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) / (rmu[k] * Pi)) : -0.5 * (rn[k+1] - rn[k-1]) / HP
                q = clamp((ForwardDiff.value(lp / lph) - 0.25) / 0.5, 0.0, 1.0); sw = q * q * (3 - 2q)
                sw < 1 && (lp = sw * lp + (1 - sw) * lph)
            end
            abs(ForwardDiff.value(lp)) < 1e-12 && continue
            ∇ = log(tn[k] / tn[k-1]) / lp
            if MLT_AD[]                              # coefficients on the current state, differentiated
                ZONE[] = k
                if mpre && mok[k]                    # from the threaded pre-pass: the same calls on the same arguments
                    nad, cpk, δk, κk = W.mnad[k], W.mcp[k], W.mdl[k], W.mkp[k]
                else
                    th = material_thermo(MATERIAL[], ρi, Ti); nad, cpk, δk = th
                    κk = material_opac(MATERIAL[], ρi, Ti)[1]
                end
            else                                     # lagged: from the Newton iterate, constant in the Jacobian
                nad, cpk, δk, κk = nadv[k], cpv[k], δv[k], κv[k]
                cpk > 0 || continue
            end
            W.nabv[k] = ∇; W.nadv[k] = nad
            sup = ∇ - nad                            # superadiabaticity (not `x`: that is the Newton vector)
            ForwardDiff.value(sup) > 0 || continue
            U = 3arad * cc * Ti^3 / (cpk * ρi^2 * κk * l^2) * sqrt(8HP / (gk * δk))
            ζ = sqrt(U^2 + sup) - U
            fcv[k] = wτ * ρi * cpk * Ti * sqrt(gk * δk) * l^2 * ζ^3 / (4sqrt(2.0) * HP * sqrt(HP))
            W.vcv[k] = l * sqrt(gk * δk / (8HP)) * ζ; W.lcv[k] = l       # KW eq. 7.x: v̄ = ℓ (gδ/8H_P)^½ (∇−∇')^½, ζ² = ∇−∇'
        end
    end
    # ---------------------------------------------------------------- GAS ENERGY
    hassrc = MATERIAL[] !== nothing && has_source(MATERIAL[]); srco = SRCO[]
    @inbounds for k in ngrs:ngre
        aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
        # S4: nuclear source ε ρ ΔV (erg s⁻¹ per steradian), θ-weighted new/old; subtracted in every variant below
        src = hassrc ? (θ * W.snuc[k] + (1 - θ) * srco[k]) * d[k] * dvol[k] : zero(T)
        hassrc && (W.srcw[k] = src)
        if ECOMPAT[]                           # compatible first-law row in every zone (see ECOMPAT)
            egbg = c.lrad ? W.egbg : W.egb        # advected gas energy (lrad = 0: the only energy)
            eold = c.lrad ? c.egro[k] - c.ero[k] / c.do_[k] : c.ego[k]
            mn_ = dn[k] * dvoln[k]; mo_ = c.do_[k] * c.dvolo[k]
            r3k = ((r[k] + r[k+1]) / 2)^(1.5xmu)
            vlo = k == ngrs ? rmu[k] * u[k] : r3k * r[k]^(-0.5xmu) * u[k]
            vhi = (k == ngre && isnan(PEXT[])) ? rmu[k+1] * u[k+1] : r3k * r[k+1]^(-0.5xmu) * u[k+1]
            du2 = (un[k] - c.uo[k])^2 + (un[k+1] - c.uo[k+1])^2
            w = (mn_ * egn[k] - mo_ * eold) / dt -
                dfe[k+1] + dfe[k] + rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] +
                pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * (vhi - vlo) +
                (-(2θ - 1) * 0.125 * (mn_ + mo_) + (θ^2 - θ + 0.5) * 0.25 * (mn_ - mo_)) * du2 / dt +
                0.5 * (dmdt[k] + dmdt[k+1]) * (u[k+1] - u[k]) * (0.5 * (u[k] + u[k+1]) - W.ub[k]) -
                (ksflux(k + 1, u, dmdt, ngrs, ngre) - ksflux(k, u, dmdt, ngrs, ngre)) + (isempty(KSTAR[]) ? 0.0 : KSTAR[][k] / dt) -
                (c.lrad ? (Mg === nothing ? excap(cc, xke[k], d[k], dt) * (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) : Mg.xch[k]) * d[k] * dvol[k] -
                          compton(k, W, P) * dvol[k] : zero(T)) - src
            isempty(GASQ[]) || (w -= GASQ[][k] * dvol[k])
            w -= yqheat(k, W) * dvol[k]                         # cycle heat: the iterate's (coupled zone) or the booked one
            if k == ngrs
                w += -dmdt[k+1] * egbg[k+1] - P.phil2 * rmu[k]
            elseif k == ngre
                w += dmdt[k] * egbg[k]
                w += P.leobc != 3 ? rmu[k+1] * P.phir2 : -dmdt[k+1] * egbg[k+1]
            else
                w += -dmdt[k+1] * egbg[k+1] + dmdt[k] * egbg[k]
            end
        elseif c.etotal && !c.useint[k]
            w = etot_row(k, c, W, T, AD_JV[] ? fedd : c.fedd) + rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] - src
            isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])     # the gas gets the ray's energy in chemistry: it leaves the total here (as the lrad total row)
            w += cqheat(k, W) * dvol[k]                         # chemistry heat: the iterate's (coupled zone; inside cpL as heating of the gas) or the booked one
        elseif c.lrad && GASROW[]
            egbg = W.egbg
            w = (dn[k] * egn[k] * dvoln[k] - c.do_[k] * c.ego[k] * c.dvolo[k]) / dt -
                dfe[k+1] + dfe[k] + rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] +
                pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * dudr[k] * dvol[k] -
                (Mg === nothing ? excap(cc, xke[k], d[k], dt) * (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) : Mg.xch[k]) * d[k] * dvol[k] +
                compton(k, W, P) * dvol[k] - src                   # (MGROUPS: Σ_g of the group exchange)
            isempty(GASQ[]) || (w -= GASQ[][k] * dvol[k])
            w -= yqheat(k, W) * dvol[k]                         # cycle heat: the iterate's (coupled zone) or the booked one
            if k == ngrs
                w += -dmdt[k+1] * egbg[k+1] - P.phil2 * rmu[k]
            elseif k == ngre
                w += dmdt[k] * egbg[k]
                w += P.leobc != 3 ? rmu[k+1] * P.phir2 : -dmdt[k+1] * egbg[k+1]
            else
                w += -dmdt[k+1] * egbg[k+1] + dmdt[k] * egbg[k]
            end
        elseif c.lrad
            w = (dn[k] * egrn[k] * dvoln[k] - c.do_[k] * c.egro[k] * c.dvolo[k]) / dt +
                rmu[k+1] * fr[k+1] - rmu[k] * fr[k] - dfe[k+1] + dfe[k] +
                rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] +
                (pg[k] + er[k] * fedd[k]) * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur * dvol[k] + qf[k] * dudr[k] * dvol[k] - src
            isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])
            w += cqheat(k, W) * dvol[k]
            if k == ngrs
                w += -dmdt[k+1] * egb[k+1] - (P.phil2 + er[k] * urel[k]) * rmu[k]
                INNER_SINK[] && (w -= rmu[k] * P.phil0 * (1 + fedd[k]) * er[k] / d[k])
            elseif k == ngre
                w += dmdt[k] * egb[k]
                w += P.leobc != 3 ? rmu[k+1] * (P.phir2 + urel[k+1] * er[k]) : -dmdt[k+1] * egb[k+1]
            else
                w += -dmdt[k+1] * egb[k+1] + dmdt[k] * egb[k]
            end
        else
            w = (dn[k] * egn[k] * dvoln[k] - c.do_[k] * c.ego[k] * c.dvolo[k]) / dt -
                dfe[k+1] + dfe[k] + pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * dudr[k] * dvol[k] +
                rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] - src
            if k == ngrs
                w += -dmdt[k+1] * egb[k+1] - P.phil2 * rmu[k]
            elseif k == ngre
                w += dmdt[k] * egb[k]
                w += P.leobc != 3 ? rmu[k+1] * P.phir2 : -dmdt[k+1] * egb[k+1]
            else
                w += -dmdt[k+1] * egb[k+1] + dmdt[k] * egb[k]
            end
        end
        isempty(CHEMB[]) || (w -= d[k] * CHEMB[][k] * dvol[k])     # BALSPLIT: the chemistry's lagged mean heating (see CHEMB)
        cplz(k) && (θc = max(θ, CPL_THETA[]); w += (θc * W.cpL[k] + (1 - θc) * XLAMO[][k]) * dvol[k])   # COUPLED: the chemistry cooling, implicit
        R[it, kb(k)] = w
    end
    R[it, 1] = X(it, k0)
    k = ngrs - 1
    R[it, kb(k)] = (P.leibc == 1 || P.llibc == 1) ? egn[k+1] - egn[k] : X(it, k)
    for k in ngre+1:ngre+2
        R[it, kb(k)] = (P.leobc in (1, 3) || P.llobc > 0) ? -egn[k-1] + egn[k] : X(it, k)
    end
    if COUPLED[]
        # ------------------------------------------------------------ SPECIES (COUPLED)
        @inbounds for kb_ in 1:nb, i in 8:9; R[i, kb_] = AD_JV[] ? zero(T) : x[(kb_-1)*neqn+i]; end
        # no radiation (lrad = 0, the collapse stages): rows 6–7 (E_r, F_r) are not formed above, but neqn = 9 puts the species at 8–9;
        # make them trivial (δ = 0) instead of leaving the buffer's stale values (a singular block: all Newton corrections 0, step 1 fails)
        c.lrad || @inbounds for kb_ in 1:nb, i in 6:7; R[i, kb_] = AD_JV[] ? zero(T) : x[(kb_-1)*neqn+i]; end
        @inbounds for k in ngrs:ngre
            cplz(k) || continue
            for (i, q, qo, qb, S_, So) in ((8, W.xs1, XSO[][1], W.xb1, W.cpS1, XSRCO[][1]), (9, W.xs2, XSO[][2], W.xb2, W.cpS2, XSRCO[][2]))
                θc = max(θ, CPL_THETA[])
                w = (dn[k] * q[k] * dvoln[k] - c.do_[k] * qo[k] * c.dvolo[k]) / dt -
                    (θc * dn[k] * dvoln[k] * S_[k] + (1 - θc) * c.do_[k] * c.dvolo[k] * So[k])
                if k == ngrs
                    w += -dmdt[k+1] * qb[k+1]
                elseif k == ngre
                    w += dmdt[k] * qb[k]
                    P.leobc == 3 && (w += -dmdt[k+1] * qb[k+1])
                else
                    w += -dmdt[k+1] * qb[k+1] + dmdt[k] * qb[k]
                end
                R[i, kb(k)] = w
            end
        end
    end
    if c.lrad
        # ------------------------------------------------------------ RADIATION ENERGY
        @inbounds for k in ngrs:ngre
            aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
            if ECOMPAT[] || ERLAB[]            # lab-frame conservation law: flux of u P and the work on the gas (see ECOMPAT)
                fEl = k == ngrs ? fedd[k] * er[k] : 0.5 * (fedd[k-1] * er[k-1] + fedd[k] * er[k])
                fEu = k == ngre ? fedd[k] * er[k] : 0.5 * (fedd[k] * er[k] + fedd[k+1] * er[k+1])
                wk = rmu[k+1] * u[k+1] * fEu - rmu[k] * u[k] * fEl
                k > ngrs && (wk += u[k] * avchi[k] * fr[k] / cc * 0.5 * d[k] * dvol[k])
                k < ngre && (wk += u[k+1] * avchi[k+1] * fr[k+1] / cc * 0.5 * d[k] * dvol[k])
                w = (ern[k] * dvoln[k] - c.ero[k] * c.dvolo[k]) / dt + rmu[k+1] * fr[k+1] - rmu[k] * fr[k] + wk +
                    excap(cc, xke[k], d[k], dt) * (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) * d[k] * dvol[k] - compton(k, W, P) * dvol[k]
            else
            w = (ern[k] * dvoln[k] - c.ero[k] * c.dvolo[k]) / dt + rmu[k+1] * fr[k+1] - rmu[k] * fr[k] +
                fedd[k] * er[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur * dvol[k] +
                excap(cc, xke[k], d[k], dt) * (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) * d[k] * dvol[k] - compton(k, W, P) * dvol[k]
            end
            isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])
            isempty(GASQ[]) || (w += GASQ[][k] * dvol[k])
            w += (cqheat(k, W) + yqheat(k, W)) * dvol[k]          # chemistry and cycle heat
            if k == ngrs
                w += -dmdt[k+1] * erb[k+1] - urel[k] * rmu[k] * er[k]
                INNER_SINK[] && (w -= rmu[k] * P.phil0 * (1 + fedd[k]) * er[k] / d[k])
            elseif k == ngre
                w += dmdt[k] * erb[k]
                w += P.leobc != 3 ? urel[k+1] * rmu[k+1] * er[k] : -dmdt[k+1] * erb[k+1]
            else
                w += -dmdt[k+1] * erb[k+1] + dmdt[k] * erb[k]
            end
            R[ie, kb(k)] = w
        end
        R[ie, 1] = X(ie, k0)
        R[ie, kb(ngrs-1)] = ern[ngrs] - ern[ngrs-1]
        P.lrobc == 1 || error("lrobc must be 1")
        R[ie, kb(ngre+1)] = -ern[ngre] + ern[ngre+1]
        R[ie, kb(ngre+2)] = -ern[ngre+1] + ern[ngre+2]
        # ------------------------------------------------------------ RADIATION MOMENTUM
        qF = W.qF
        tie = TD.FLUXTIE[][2] > 0 && length(TD.TIEW) >= ngre + 1 && length(TD.HEDD) >= ngre + 1
        sph = TD.VEF_SPH[]                   # Auer's sphericality-factor form of the static terms (TitanDriver.VEF_SPH)
        @inbounds for k in ngrs+1:ngre
            qFo = c.fro[k] * 2 / (c.do_[k] + c.do_[k-1])
            if !sph
            R[iF, kb(k)] =
                SDFIN[] * (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                 qFo * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
                ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2) +
                rmu[k] * (fedd[k] * er[k] - fedd[k-1] * er[k-1]) +
                0.5 * rmu[k] * (fr[k] / cc^2) * (u[k+1] - u[k-1]) +
                (0.25xmu / r[k]) * ((3fedd[k-1] - 1) * er[k-1] * dvol[k-1] + (3fedd[k] - 1) * er[k] * dvol[k]) +
                0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] * fr[k] / cc
            else                                     # VEF_SPH: static terms in the sphericality-factor form (TitanDriver.VEF_SPH)
                qs = vefq(fedd, r, xmu, k)
                R[iF, kb(k)] =
                    SDFIN[] * (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                     qFo * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
                    ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2) +
                    rmu[k] * (qs * fedd[k] * er[k] - fedd[k-1] * er[k-1] / qs) +
                    0.5 * rmu[k] * (fr[k] / cc^2) * (u[k+1] - u[k-1]) +
                    0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] * fr[k] / cc
            end
            wt = tie ? TWl[k] : 0.0
            SDW[1] > 0 && SDW[1] <= k - ngrs + 1 <= SDW[2] && tie && (wt = SDW[3])      # stalldiag (SDW=lo:hi:w): tie weight forced at interfaces lo..hi (zone numbers)
            if wt > 0                                # FLUXTIE: F = h c E (TitanDriver.FLUXTIE)
                rcm = (0.5 * (r[k-1] + r[k]))^mu; rcp = (0.5 * (r[k] + r[k+1]))^mu
                # upwind (donor cell): outward light from the zone inside, inward from the zone outside, diluted as r^-2. The centred
                # average ½(hE_(k−1) + hE_k) cannot see an odd–even E mode, which nothing damps in thin gas (E → 0 at zone 214, step 87126)
                # TIEX (TitanDriver.TIEX): the closure is for the excess over the isotropic background, h_x (E − E_b)
                tx = TD.TIEX[] && TD.VEF_BG[] && length(TD.HXEDD) >= ngre + 1 && length(TD.FXEDD) >= ngre + 1
                HT = tx ? HXl : HEl; Ebg = tx ? 4P.cpi * P.ximr / cc : 0.0; em = er[k-1] - Ebg; ep = er[k] - Ebg
                if !SLABTIE[]
                B = rmu[k] * fr[k] / cc - (rcm * max(HT[k-1], 0.0) * em + rcp * min(HT[k], 0.0) * ep)
                else                                 # SLABTIE: outward donor h from the slab depth of the donor's centre below the interface
                    hd = st_hdonor(max(HT[k-1], 0.0), SDTIEF > 0 ? SDTIEF : (tx ? FXl[k-1] : fedd[k-1]), chif[k-1], chif[k], 0.5 * d[k-1] * dvol[k-1] / rmu[k], 0.5 * d[k] * dvol[k] / rmu[k])
                    if SDFRZ    # stalldiag (test only): the slab correction h_slab(tau_i) - h_slab(a) frozen at its first value (no feedback through chi, rho, dV)
                        _fm = tx ? FXl[k-1] : fedd[k-1]; _Mm = 0.5 * d[k-1] * dvol[k-1] / rmu[k]; _Mp = 0.5 * d[k] * dvol[k] / rmu[k]
                        _ti = st_tauint(chif[k-1], chif[k], _Mm, _Mp)
                        length(SDFRZV) < TF.MGR && (resize!(SDFRZV, TF.MGR); fill!(SDFRZV, NaN))
                        isnan(SDFRZV[k]) && (SDFRZV[k] = val(_fm / (_ti + 2_fm) - _fm / (chif[k-1] * _Mm + 2_fm)))
                        hd = max(max(HT[k-1], 0.0) + SDFRZV[k], 0.0)
                    end
                    B = rmu[k] * fr[k] / cc - (rcm * hd * em + rcp * min(HT[k], 0.0) * ep)
                    if SDHD    # stalldiag: the argument of SLABTIE's max(h, 0) (raw h), h_FS, h_slab terms at the last residual evaluation
                        _fm = tx ? FXl[k-1] : fedd[k-1]; _Mm = 0.5 * d[k-1] * dvol[k-1] / rmu[k]; _Mp = 0.5 * d[k] * dvol[k] / rmu[k]
                        _ti = st_tauint(chif[k-1], chif[k], _Mm, _Mp)
                        SDHRAW[k] = (val(max(HT[k-1], 0.0) + _fm / (_ti + 2_fm) - _fm / (chif[k-1] * _Mm + 2_fm)), val(HT[k-1]), val(_ti), val(chif[k-1] * _Mm))
                    end
                end
                if TIEFRAME[]                        # the unknown F is the comoving flux: F₀ = F_lab − (1 + f) u E + 2 h² u E (see TIEFRAME)
                    ho = SLABTIE[] ? hd : max(HT[k-1], 0.0); hi = min(HT[k], 0.0)
                    eup = decide(:up, 21, k, val(rn[k]) >= c.ro[k]) ? er[k] : er[k-1]
                    B += (-dmdt[k] * erb[k] + (rn[k]^(mu + 1) - c.ro[k]^(mu + 1)) / (xmup1 * dt) * eup +
                          rmu[k] * u[k] * 0.5 * (fedd[k-1] * er[k-1] + fedd[k] * er[k]) -
                          2u[k] * (rcm * ho^2 * em + rcp * hi^2 * ep)) / cc
                end
                if SDROW && SDW0 <= k - ngrs + 1 <= SDW1     # stalldiag: the two mixed rows and their coefficients on F and E (values)
                    cF_m = (qFo * 0 + 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] / cc)
                    SDROWV[k] = (val(R[iF, kb(k)]), val(B), val(wt), val(cF_m), val(rmu[k] / cc), val(rmu[k] * fedd[k-1]), val(rcm * (SLABTIE[] ? hd : max(HT[k-1], 0.0))), val(SDFIN[] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) / (2cc^2 * dt) * 2 / (dn[k] + dn[k-1])))
                end
                R[iF, kb(k)] = (1 - wt) * R[iF, kb(k)] + wt * B
            end
        end
        P.lribc in (2, 3) || error("lribc must be 2 or 3")
        R[iF, kb(ngrs)] = X(iF, ngrs)
        R[iF, kb(ngre+1)] = P.cpi * (2gedr + 1) * P.ximr + frn[ngre+1] - cc * gedr * (sph ? vef_esurf(W, (length(FXl) >= ngre + 1 ? FXl[ngre] : fedd[ngre]), ngre, xmu, cc, 4P.cpi * P.ximr / cc) : ern[ngre])
        if TIEFRAME[]                                # the boundary condition holds for the lab-frame flux (see TIEFRAME)
            ub_ = un[ngre+1]
            R[iF, kb(ngre+1)] = (R[iF, kb(ngre+1)] - frn[ngre+1]) * (1 + 2gedr * ub_ / cc) + frn[ngre+1] + (1 + fedd[ngre]) * ub_ * ern[ngre]
        end
        R[iF, 1] = X(iF, k0)
        k = ngrs - 1
        R[iF, kb(k)] = P.lribc == 2 ? frn[k+2] + frn[k] : rn[k+1]^mu * frn[k+1] - rn[k]^mu * frn[k]
        R[iF, kb(ngre+2)] = -frn[ngre+1] + frn[ngre+2]
        Mg === nothing || mg_rows!(R, Mg, c, W, x, kb, P, dt, θ, xmu, cc, sph, tie)   # MGROUPS: sum rows and the group rows
    end
    nothing
end

"""
VEF_SPH: e^{s/2} at interface k, s = (μ/4)(3 − 1/f_k−1 + 3 − 1/f_k) ln(r_c,k/r_c,k−1), r_c the zone midpoints (TitanDriver.VEF_SPH).
The momentum row's static terms are r_k² [e^{s/2} f_k E_k − e^{−s/2} f_k−1 E_k−1]: the ratio of Auer's sphericality factor at the
zone centres to its value at the interface, so that q f E is differenced across the interface. f ≥ 0.01 only guards 1/f.
"""
@inline function vefq(fedd, r, xmu, k)
    gq(f) = 3 - 1 / max(f, 0.01)
    rcm = 0.5 * (r[k-1] + r[k]); rcp = 0.5 * (r[k] + r[k+1])
    exp(0.125xmu * (gq(fedd[k-1]) + gq(fedd[k])) * log(rcp / rcm))
end
"""
VEF_SPH: energy density at the outer boundary radius r_N+1 from the outer zone's centre value (new time level, as TITAN's boundary
row). The incoming background (the CMB, intensity I_in = ximr, E_b = 4π I_in/c) is isotropic, K_b = E_b/3, and carried unchanged
to r_N+1; only the excess E_x = E_N − E_b is carried by the static momentum equation d(q f_x E_x)/dr = −q χρ F/c over the half zone
r_c,N → r_N+1 with its own Eddington factor f_x = K_x/J_x of the excess, taken from the formal solution of the excess (TD.FXEDD, a Dual through closure_lift;
fedd[N] without the background tie), held over the half zone. (It is not reconstructed as (f_N E_N − E_b/3)/E_x: where E_x ≪ E_b that is a
difference of two nearly equal numbers, ~2.2 at the quadrature error of f, and any clamp of it removes the E_r(N) column of the boundary row.)
E_s = E_b + e^{−s} E_x − e^{−s/2} χ_N ρ_N F_N+1 (r_N+1 − r_c,N)/(c f_x), s = (μ/2)(3 − 1/f_x) ln(r_N+1/r_c,N).
"""
@inline function vef_esurf(W, fx, N, xmu, cc, Eb)
    rn = W.rn; rc = 0.5 * (rn[N] + rn[N+1]); Ex = W.ern[N] - Eb
    s = 0.5xmu * (3 - 1 / fx) * log(rn[N+1] / rc)
    Eb + exp(-s) * Ex - exp(-0.5s) * W.chifn[N] * W.dn[N] * W.frn[N+1] * (rn[N+1] - rc) / (cc * fx)
end
# ------------------------------------------------------------------ total energy (item 2)
"Total pressure (gas + radiation fE) of zone k at the θ-level."
@inline ptot(k, c, W, fedd = c.fedd) = c.lrad ? W.pg[k] + fedd[k] * W.er[k] : W.pg[k]
"Viscous stress weight of zone k at interface k' (momentum equation: r1(k')·r3(k)·qf(k))."
@inline r3(k, W, xmu) = ((W.r[k] + W.r[k+1]) / 2)^(1.5xmu)
"""
Energy flux through interface k (positive outward, θ-centred), consistent with
TITAN's discretisation of each piece:
  advection  −ṁ_k (ê_k + u_k²/2)        ê = upwinded e_g (+E_r/ρ) exactly as in gasnrg
  work        r_k^μ u_k P̄_k              P̄ = (P_{k-1}+P_k)/2,  P = p_g (+ f E_r)
  viscous     u_k r_k^{-μ/2} (r3_{k-1} q_{k-1} + r3_k q_k)/2   (momentum-equation weights)
  radiation   r_k^μ F_k ;  diffusion  −dfe_k
Boundary interfaces use TITAN's own boundary fluxes (gasnrg/totnrg) plus the matching
kinetic-energy, work and viscous pieces of the adjacent zone.
"""
@inline function eflux(k, c, W, ::Type{T}, fedd = c.fedd) where {T}
    P = c.p; ngrs, ngre = c.ngrs, c.ngre; xmu = Float64(c.mu)
    u = W.u[k]; rmu = W.rmu[k]; dmdt = W.dmdt[k]
    Fl = zero(T)
    c.lrad && (Fl += rmu * W.fr[k])
    Fl -= W.dfe[k]
    if k == ngrs
        Fl += c.lrad ? (P.phil2 + W.er[k] * W.urel[k]) * rmu : P.phil2 * rmu
        Fl += -dmdt * 0.5 * u^2 + rmu * u * ptot(k, c, W, fedd) + rmu * u * W.qf[k]
    elseif k == ngre + 1
        if P.leobc != 3
            Fl += c.lrad ? rmu * (P.phir2 + W.urel[k] * W.er[k-1]) : rmu * P.phir2
        else
            Fl += -dmdt * W.egb[k]
        end
        Fl += -dmdt * 0.5 * u^2 + rmu * u * ptot(k - 1, c, W, fedd) + rmu * u * W.qf[k-1]
    else
        Fl += -dmdt * (W.egb[k] + 0.5 * u^2) + rmu * u * 0.5 * (ptot(k - 1, c, W, fedd) + ptot(k, c, W, fedd)) +
              u * W.r[k]^(-0.5xmu) * 0.5 * (r3(k - 1, W, xmu) * W.qf[k-1] + r3(k, W, xmu) * W.qf[k])
    end
    Fl
end
"ETOT_GRAV[]: the total-energy row includes the work of gravity (self-gravity or P.g). Default false, as in the tests and ledgers that book the potential energy themselves."
const ETOT_GRAV = Ref(false)
"Run f with ECOMPAT off: the conservative total-energy rows and the plain internal rows of the diagnostics (the independent check of the compatible rows)."
function nocompat(f)
    e = ECOMPAT[]; l = ERLAB[]; ECOMPAT[] = false; ERLAB[] = false
    try f() finally ECOMPAT[] = e; ERLAB[] = l end
end
"""
Gravity term of the total-energy row at interior face j, the share of the zone above (lower = false) or below (lower = true) the face. With the potential at the interfaces
φ_j = −G4 (xm_j − xe_j)/r_j + (1 − μ/2) g₀ r_j (G4 = 4πG for μ = 2; xe the external mass) the exact rate of change of the binding energy W = Σ_j m̄_j φ_j of the node masses
m̄_j = ½(m_{j-1} + m_j), with the masses changing by the mass fluxes, is  Σ_j [ g_j m̄_j ṙ_j + f_j Δ_j ]  with  ṙ_j = u_j − u_rel,j,  f_j = r_j^μ ρᵇ_j u_rel,j  (the continuity row's flux) and
Δ_j = ½(φ_{j+1} − φ_{j-1}) + G4 m̄_j/r_j + g_ext,j ½(r_{j+1} − r_{j-1}); the middle term is the change of the enclosed mass M_j when mass crosses the face (docs/POPIII_ENERGY.md §23.4).
Spherical geometry, or no self-gravity.
"""
@inline function gravface(j, c, W, lower::Bool = false)
    P = c.p; xmu = Float64(c.mu); r = W.r; xm = W.xm; G4 = 2xmu * P.cpi * P.cgrav
    xe(i) = EXTMASS[] === nothing ? zero(r[i]) : EXTMASS[](r[i])
    xmu == 2 || P.cgrav == 0 || error("gravity of the total-energy row: spherical geometry (or no self-gravity) only")
    φ(i) = (P.cgrav == 0 || val(r[i]) <= 0 ? zero(r[i]) : -G4 * (xm[i] - xe(i)) / r[i]) + (1 - 0.5xmu) * P.g * r[i]      # (regular centre: M/r → 0)
    mb = 0.5 * (W.d[j-1] * W.dvol[j-1] + W.d[j] * W.dvol[j])
    gj = G4 * xm[j] / W.rmu[j] + (1 - 0.5xmu) * P.g
    gext = P.cgrav == 0 ? zero(r[j]) : G4 * xe(j) / W.rmu[j]
    Δ = 0.5 * (φ(j + 1) - φ(j - 1)) + (P.cgrav == 0 ? zero(r[j]) : G4 * mb / r[j]) + gext * 0.5 * (r[j+1] - r[j-1])
    f = W.rmu[j] * W.db[j] * W.urel[j]
    wl = (r[j] - r[j-1]) / (r[j+1] - r[j-1])          # half-width weights: the zone below gets (r_j − r_j-1)/(r_j+1 − r_j-1)
    (gj * mb * (W.u[j] - W.urel[j]) + f * Δ) * (lower ? wl : 1 - wl)
end
"Conservative total-energy row of zone k: Δ[m(e_g [+E_r/ρ]) + m(u_k²+u_{k+1}²)/4]/dt + Φ_{k+1} − Φ_k."
@inline function etot_row(k, c, W, ::Type{T}, fedd = c.fedd) where {T}
    dn, un, dvoln = W.dn, W.un, W.dvoln
    en = c.lrad ? W.egrn[k] : W.egn[k]; eo = c.lrad ? c.egro[k] : c.ego[k]
    w = (dn[k] * dvoln[k] * (en + 0.25 * (un[k]^2 + un[k+1]^2)) -
         c.do_[k] * c.dvolo[k] * (eo + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))) / c.p.dt +
        eflux(k + 1, c, W, T, fedd) - eflux(k, c, W, T, fedd)
    if ETOT_GRAV[]
        # Gravity of the monitor and of the dual-energy stack's total rows: the exact rate of change of the discrete binding energy W = Σ_j m̄_j φ_j of the node masses
        # m̄_j = ½(m_{j-1} + m_j) at the interface potentials φ_j (docs/POPIII_ENERGY.md §23.4), whose gradient at fixed masses is the momentum row's force g_j m̄_j.
        # Face j contributes g_j m̄_j (u_j − u_rel,j) [the mesh moves through the gas] + f_j Δ_j [the mass flux f_j = r_j^μ ρᵇ_j u_rel,j carries the potential drop Δ_j], shared
        # between the two zones in proportion to their half widths.
        k > c.ngrs && (w += gravface(k, c, W))
        k < c.ngre && (w += gravface(k + 1, c, W, true))
        if k == c.ngre && !isnan(PEXT[])          # free surface (the momentum row of the last interface has the gas of zone ngre alone): the mesh carries half of its mass
            xmu = Float64(c.mu); P = c.p; j = k + 1
            w += (2xmu * P.cpi * P.cgrav * W.xm[j] / W.rmu[j] + (1 - 0.5xmu) * P.g) * 0.5 * W.d[k] * W.dvol[k] * (W.u[j] - W.urel[j])
        end
    end
    w
end

"""
    energy_imbalance(dt) -> (; dE, flux_in, flux_out, imb)

Call right after a converged step (old arrays = start of step; pass the dt actually
taken).  Evaluates the conservative total-energy rows at the converged state: their
sum × dt is exactly ΔE_total + dt(Φ_out − Φ_in), i.e. the step's discrete total-energy
non-conservation, whichever formulation produced the step.
"""
function energy_imbalance(dt)
    c = context(; etotal = true, ηsw = 0.0, dt)
    nb = c.nb; neqn = c.neqn
    R = zeros(neqn, nb); ws = getws(Float64, c.mgr)
    nocompat(() -> residual!(R, ws, c, zeros(neqn * nb)))
    it = c.ix[5]; W = ws.a
    imb = sum(R[it, k - c.k0 + 1] for k in c.ngrs:c.ngre) * dt
    E(k) = W.dn[k] * W.dvoln[k] * ((c.lrad ? W.egrn[k] : W.egn[k]) + 0.25 * (W.un[k]^2 + W.un[k+1]^2))
    Eo(k) = c.do_[k] * c.dvolo[k] * ((c.lrad ? c.egro[k] : c.ego[k]) + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))
    dE = sum(E(k) - Eo(k) for k in c.ngrs:c.ngre)
    (; dE, Enew = sum(E(k) for k in c.ngrs:c.ngre), fin = eflux(c.ngrs, c, W, Float64) * dt,
       fout = eflux(c.ngre + 1, c, W, Float64) * dt, imb)
end

"""
    etot_monitor(dt) -> (; imb, thr, G, gnet, zmax, vmax, zshare, vshare, dE, fin, fout)

(ETOTMON, diagnostic; no influence on the solution.) Right after a converged implicit solve (`TitanDriver.STAGE_HOOK`, old arrays = the solve's start,
`dt` = the solve's step): one Float64 residual evaluation of the CONSERVATIVE total-energy rows of every zone (ηsw = 0, ETOT_GRAV on, whatever
formulation produced the solve): row_k = Δ[m(e_g + E_r/ρ + K)]_k/dt + Φ_k+1 − Φ_k (`etot_row`: the fluxes Φ are `eflux`, including the boundary fluxes at the
first and last interface; the nuclear source, the star-light/ray heating RAYQ and the chemistry cooling are in the row as in the solve). The sum over the
zones × dt is therefore the global budget  ΔE_tot + dt (Φ_out − Φ_in) − ∫(sources),  `imb`  [erg sr⁻¹], the interior fluxes telescoping.
Scales: `thr` = Σ_k (|ΔE_k| + dt (|Φ_k| + |Φ_k+1|)), the solve's energy throughput, so  G = imb/thr; `gnet` = imb/(|ΔE_tot| + dt(|Φ_in| + |Φ_out|)).
Per zone: v_k = |row_k| dt/thr_k with thr_k = |ΔE_k| + dt (|Φ_k| + |Φ_k+1|) (`zmax`, `vmax`: the worst, relative to its own throughput) and
s_k = |row_k| dt/thr (`zshare`, `vshare`: the zone that contributes most to the global budget). Zone numbers count from the first zone (1).
"""
function etot_monitor(dt)
    eg0 = ETOT_GRAV[]; ETOT_GRAV[] = true
    # (context() fills reused buffers in place: the solver's row-switch flags and base state are restored, the monitor has no influence on the solve)
    us = haskey(BUFS, :useint) ? copy(BUFS[:useint]::Vector{Bool}) : nothing; bs = haskey(BUFS, :base) ? copy(BUFS[:base]::Matrix{Float64}) : nothing
    try
        c = context(; etotal = true, ηsw = 0.0, dt)
        nb, neqn = c.nb, c.neqn
        R = zeros(neqn, nb); ws = getws(Float64, c.mgr)
        nocompat(() -> residual!(R, ws, c, zeros(neqn * nb))); W = ws.a
        it = c.ix[5]; a, b = c.ngrs, c.ngre
        E(k) = W.dn[k] * W.dvoln[k] * ((c.lrad ? W.egrn[k] : W.egn[k]) + 0.25 * (W.un[k]^2 + W.un[k+1]^2))
        Eo(k) = c.do_[k] * c.dvolo[k] * ((c.lrad ? c.egro[k] : c.ego[k]) + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))
        Φ = [eflux(k, c, W, Float64) for k in a:b+1]
        res = [R[it, k - c.k0 + 1] * dt for k in a:b]
        thrk = [abs(E(k) - Eo(k)) + dt * (abs(Φ[k-a+1]) + abs(Φ[k-a+2])) for k in a:b]
        thr = sum(thrk); imb = sum(res); dE = sum(E(k) - Eo(k) for k in a:b)
        fin, fout = Φ[1] * dt, Φ[end] * dt
        v = abs.(res) ./ max.(thrk, 1e-300); zm = argmax(v)
        s = abs.(res) ./ max(thr, 1e-300); zs = argmax(s)
        (; imb, thr, G = imb / max(thr, 1e-300), gnet = imb / max(abs(dE) + abs(fin) + abs(fout), 1e-300), zmax = zm, vmax = v[zm], zshare = zs, vshare = s[zs], dE, fin, fout)
    finally
        ETOT_GRAV[] = eg0
        us === nothing || copyto!(BUFS[:useint]::Vector{Bool}, us); bs === nothing || copyto!(BUFS[:base]::Matrix{Float64}, bs)
    end
end

"""
    step_diag(dt) -> (; srcw, fcv, vcv, lcv)

Right after a converged step (old arrays = start of step): one Float64 residual evaluation at the converged state (as
`energy_imbalance`) returning per zone/interface (copies) the source power srcw = ε_θ ρ ΔV [erg s⁻¹ per steradian] that the
energy row injected, and the MLT flux F_conv (weight included), mean convective velocity v and mixing length ℓ of each interface.
"""
function step_diag(dt)
    c = context(; dt)
    R = zeros(c.neqn, c.nb); ws = getws(Float64, c.mgr)
    residual!(R, ws, c, zeros(c.neqn * c.nb)); W = ws.a
    (; srcw = copy(W.srcw), fcv = copy(W.fcv), vcv = copy(W.vcv), lcv = copy(W.lcv), nab = copy(W.nabv), nad = copy(W.nadv))
end

"""
    stage_diag(dt) -> (; srcw, fcv, vcv, lcv, dE, fin, fout, imb)

`step_diag` and `energy_imbalance` in one residual evaluation (conservative total-energy rows): the source power injected, the MLT
side outputs, and the stage's energy change and boundary fluxes (× dt, per steradian). Used by the S5 energy ledger at every
converged implicit solve (θ-step, TR stage, BDF stage), see `TitanDriver.STAGE_HOOK`.
"""
function stage_diag(dt)
    c = context(; etotal = true, ηsw = 0.0, dt)
    nb = c.nb; neqn = c.neqn
    R = zeros(neqn, nb); ws = getws(Float64, c.mgr)
    nocompat(() -> residual!(R, ws, c, zeros(neqn * nb))); W = ws.a
    it = c.ix[5]
    imb = sum(R[it, k - c.k0 + 1] for k in c.ngrs:c.ngre) * dt
    E(k) = W.dn[k] * W.dvoln[k] * ((c.lrad ? W.egrn[k] : W.egn[k]) + 0.25 * (W.un[k]^2 + W.un[k+1]^2))
    Eo(k) = c.do_[k] * c.dvolo[k] * ((c.lrad ? c.egro[k] : c.ego[k]) + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))
    dE = sum(E(k) - Eo(k) for k in c.ngrs:c.ngre)
    (; srcw = copy(W.srcw), fcv = copy(W.fcv), vcv = copy(W.vcv), lcv = copy(W.lcv), nab = copy(W.nabv), nad = copy(W.nadv), dE,
       fin = eflux(c.ngrs, c, W, Float64) * dt, fout = eflux(c.ngre + 1, c, W, Float64) * dt, imb)
end

"""
    zone_diag(dt) -> NamedTuple of per-zone / per-interface vectors

(Diagnostic, EDIAG=2; grey rows.) Right after a converged stage: the residual T_k of the conservative total-energy row of every zone (ηsw = 0)
at the converged state [erg s⁻¹ sr⁻¹] (its sum × dt is `stage_diag`'s imb), split exactly into the terms by which the total row differs from
the rows an internal-energy formulation solves (docs/POPIII_ENERGY.md §2). With G_k the gas internal-energy row (GASROW), E_k the radiation
energy row, M_j the gas momentum row of interface j split into the parts that belong to the zone below and above it (M_j = M_j^U(j−1) + M_j^L(j):
half the inertia and weight of each zone, the momentum flux Φ_k = ½(ṁ_k + ṁ_k+1) û_k of the zone against ṁ_j u_j of the face, the zone's pressure
and viscous stress against the face means that `eflux` uses), and u the θ-level velocity:

    T_k = G_k + E_k + ½(u_k M_k + u_k+1 M_k+1)                      rows of the solve (0 where the solve used them, to Newton tolerance)
        + (Ψ_k − Ψ_k+1) + (Ψr_k − Ψr_k+1)                           a different energy flux, not a source: Ψ_j = ½ u_j (M_j^L − M_j^U)
        + K2a + K2b + K3a (+ Kbnd)                                  kinetic energy: time discretisation, mass change, upwinded momentum flux
        + prs + V4 + G5 + adv                                       pressure work (0), viscous work, gravity power, advection limiters
        + R6a + R6b + R6c (+ Rbnd)                                  radiation: ½Σ u_j × (momentum row in moment form; VEF_SPH − TITAN static
                                                                    terms; −its inertia terms)
        + old                                                       the old arrays' e_gas+rad against e_gas + E/ρ (see below)
        + rem                                                       round-off

    K2a = −(2θ − 1) ¼ m̂_k Σ_j (u_j^new − u_j^old)²/dt              m̂ = ½(m_new + m_old), j = k, k+1
    K2b = (θ² − θ + ½) ¼ (m_new − m_old)_k Σ_j (u_j^new − u_j^old)²/dt
    K3a = ½(ṁ_k + ṁ_k+1) (u_k+1 − u_k) (½(u_k + u_k+1) − û_k)       û = `ub`, the van Leer face value of the momentum row
    V4  = q_k [r̄_k^{3μ/2} (u_k+1 r_k+1^{−μ/2} − u_k r_k^{−μ/2}) − (du/dr)_k ΔV_k]
    G5  = (ETOT_GRAV terms of the total row) − ½ m_k (u_k g_k + u_k+1 g_k+1)

old = −m_old (egro − ego − ero/ρ_old)/dt: the total row differences `egro`, the gas row `ego`. They differ in COUPLED zones, where `ego` is
re-evaluated between the steps with the species of the finished step (cpl_prepare!) and `egro` is the solve's own value (docs/POPIII_ENERGY.md §2.3).
Kbnd, Rbnd: what is left in the boundary zones, whose outer faces have no momentum row (0 in interior zones: the check of the formulas).
Also returned: `Kn`, `Ko` (zone kinetic energy ¼ m (u_k² + u_k+1²) of the new and of the old arrays; in a BDF stage the old arrays are y*, and
K(y*) − (A K_γ − B K_n) is the part of the step's kinetic-energy budget no stage row sees), which zones the solve closed with the internal-energy
row (`useint`), the FLUXTIE weights, and the radiation momentum row in moment form `rmom` (tie off).
"""
function zone_diag(dt)
    c0 = context(; dt); ui = copy(c0.useint)
    neqn, nb = c0.neqn, c0.nb; z = zeros(neqn * nb); ws = getws(Float64, c0.mgr)
    cg = context(; etotal = false, dt); Rg = zeros(neqn, nb); residual!(Rg, ws, cg, z)          # internal-energy rows in every zone
    Rc = Rg                                                                                     # (ECOMPAT: the compatible rows; the split below is of the plain ones)
    (ECOMPAT[] || ERLAB[]) && (Rg = zeros(neqn, nb); nocompat(() -> residual!(Rg, ws, cg, z)))
    c = context(; etotal = true, ηsw = 0.0, dt)
    R = zeros(neqn, nb); nocompat(() -> residual!(R, ws, c, z)); W = ws.a
    iu, it, ie, iF = c.ix[4], c.ix[5], c.ix[6], c.ix[7]; a, b = c.ngrs, c.ngre; kb(k) = k - c.k0 + 1
    imb = [R[it, kb(k)] for k in a:b]; Gr = [Rg[it, kb(k)] for k in a:b]; Er = c.lrad ? [R[ie, kb(k)] for k in a:b] : zeros(b - a + 1)
    Gc = [Rc[it, kb(k)] for k in a:b]; Ec = c.lrad ? [Rc[ie, kb(k)] for k in a:b] : zeros(b - a + 1)
    dm = W.dmdt; fa(j) = c.lrad && GASROW[] ? val(dm[j] * (W.egb[j] - W.egbg[j] - W.erb[j])) : 0.0
    adv = [(k > a ? fa(k) : 0.0) - (k < b ? fa(k + 1) : 0.0) for k in a:b]
    kin = [val(0.25 * (W.dn[k] * W.dvoln[k] * (W.un[k]^2 + W.un[k+1]^2) - c.do_[k] * c.dvolo[k] * (c.uo[k]^2 + c.uo[k+1]^2)) / dt) for k in a:b]
    tw = length(TD.TIEW) >= b + 1 ? copy(TD.TIEW[a:b+1]) : zeros(b - a + 2)
    ft = TD.FLUXTIE[]; TD.FLUXTIE[] = (ft[1], 0.0)
    rm = try
        R2 = zeros(neqn, nb); nocompat(() -> residual!(R2, ws, c, z)); [R2[iF, kb(k)] for k in a:b+1]
    finally
        TD.FLUXTIE[] = ft
    end
    # ---- the split
    P = c.p; θ = P.θ; xmu = Float64(c.mu); cc = P.cc; fedd = c.fedd; lrad = c.lrad
    r, u, rmu, d, dvol, pg, er, qf, dudr, ub, fr, avchi = W.r, W.u, W.rmu, W.d, W.dvol, W.pg, W.er, W.qf, W.dudr, W.ub, W.fr, W.avchi
    un, uo = W.un, c.uo
    mn(k) = W.dn[k] * W.dvoln[k]; mo(k) = c.do_[k] * c.dvolo[k]; m(k) = d[k] * dvol[k]
    int(j) = a < j <= b
    g_(j) = 2xmu * P.cpi * P.cgrav * W.xm[j] / rmu[j] + (1 - 0.5xmu) * P.g
    r1(j) = r[j]^(-0.5xmu); r3_(k) = r3(k, W, xmu); fE(k) = lrad ? fedd[k] * er[k] : 0.0
    Φ(k) = 0.5 * (dm[k] + dm[k+1]) * ub[k]
    kfl(j) = -dm[j] * 0.5 * u[j]^2
    pfg(j) = rmu[j] * u[j] * (j == a ? pg[a] : j == b + 1 ? pg[b] : 0.5 * (pg[j-1] + pg[j]))
    pfr(j) = rmu[j] * u[j] * (j == a ? fE(a) : j == b + 1 ? fE(b) : 0.5 * (fE(j - 1) + fE(j)))
    vfl(j) = j == a ? rmu[j] * u[j] * qf[a] : j == b + 1 ? rmu[j] * u[j] * qf[b] : u[j] * r1(j) * 0.5 * (r3_(j - 1) * qf[j-1] + r3_(j) * qf[j])
    # zone-k parts of the momentum row of interior face j (j = k: the zone above the face, j = k + 1: the zone below it)
    function mp(k, j)
        lo = j == k
        ine = 0.5 * (un[j] * mn(k) - uo[j] * mo(k)) / dt
        ad = lo ? -Φ(k) + dm[j] * u[j] : Φ(k) - dm[j] * u[j]
        pt = 0.5 * (pg[j-1] + pg[j]); pr_ = lo ? rmu[j] * (pg[k] - pt) : rmu[j] * (pt - pg[k])
        gr = (GFV[] || GVAR[]) ? gfz(k, j, d, dvol, r, rmu, 2xmu * P.cpi * P.cgrav, W.xm[j], gxe(r[j])) + (1 - 0.5xmu) * P.g * 0.5 * m(k) : g_(j) * 0.5 * m(k)
        sj =0.5 * (r3_(j - 1) * qf[j-1] + r3_(j) * qf[j]); vi = lo ? r1(j) * (r3_(k) * qf[k] - sj) : r1(j) * (sj - r3_(k) * qf[k])
        ra = lrad ? -avchi[j] * fr[j] / cc * 0.5 * m(k) : 0.0
        (ine, ad, pr_, gr, vi, ra)
    end
    nz = b - a + 1; nf = b - a + 2; o = 1 - a
    Mrow = zeros(nf); Ψ = zeros(nf); Ψr = zeros(nf); mchk = zeros(nf); sns = zeros(nf); ssp = zeros(nf); dyn = zeros(nf); rchk = zeros(nf)
    sph = TD.VEF_SPH[]
    for j in a+1:b
        Mrow[j+o] = R[iu, kb(j)]
        L = sum(mp(j, j)); U = sum(mp(j - 1, j)); Ψ[j+o] = 0.5 * u[j] * (L - U); mchk[j+o] = L + U - Mrow[j+o]
        if lrad
            ft_ = 0.5 * (fE(j - 1) + fE(j)); gm(k) = (0.25xmu / r[j]) * (3fedd[k] - 1) * er[k] * dvol[k]; cf(k) = 0.5 * m(k) * avchi[j] * fr[j] / cc
            RL = rmu[j] * (fE(j) - ft_) + gm(j) + cf(j); RU = rmu[j] * (ft_ - fE(j - 1)) + gm(j - 1) + cf(j - 1)
            Ψr[j+o] = 0.5 * u[j] * (RL - RU)
            sns[j+o] = rmu[j] * (fE(j) - fE(j - 1)) + gm(j) + gm(j - 1)
            ssp[j+o] = sph ? (qs = vefq(fedd, r, xmu, j); rmu[j] * (qs * fE(j) - fE(j - 1) / qs)) : sns[j+o]
            qFo = c.fro[j] * 2 / (c.do_[j] + c.do_[j-1])
            dyn[j+o] = (W.qF[j] * (mn(j - 1) + mn(j)) - qFo * (mo(j - 1) + mo(j))) / (2cc^2 * dt) -
                       ((dm[j] + dm[j+1]) * W.frb[j] - (dm[j-1] + dm[j]) * W.frb[j-1]) / (2cc^2) + 0.5 * rmu[j] * (fr[j] / cc^2) * (u[j+1] - u[j-1])
            rchk[j+o] = rm[j+o] - (dyn[j+o] + ssp[j+o] + cf(j) + cf(j - 1))
        end
    end
    v0() = zeros(nz)
    uM, dΨ, dΨr, K2a, K2b, K3a, Kbnd, prs, V4, G5, R6a, R6b, R6c, Rbnd, rem, Kn, Ko, old, cres = (v0() for _ in 1:19)
    dKs = v0(); kst = v0()
    eg0 = ETOT_GRAV[]
    for k in a:b
        i = k + o
        Kn[i] = 0.25 * mn(k) * (un[k]^2 + un[k+1]^2); Ko[i] = 0.25 * mo(k) * (uo[k]^2 + uo[k+1]^2)
        Kraw = (Kn[i] - Ko[i]) / dt + kfl(k + 1) - kfl(k)
        Δru = rmu[k+1] * u[k+1] - rmu[k] * u[k]
        p_ = pfg(k + 1) - pfg(k) - pg[k] * Δru
        vi = vfl(k + 1) - vfl(k) - qf[k] * dudr[k] * dvol[k]
        gr = 0.0
        if eg0
            wg = etot_row(k, c, W, Float64); ETOT_GRAV[] = false
            try gr = wg - etot_row(k, c, W, Float64) finally ETOT_GRAV[] = eg0 end
        end
        aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, a))
        ra = lrad ? pfr(k + 1) - pfr(k) - fE(k) * Δru - (1 - 3fedd[k]) * 0.5xmu * er[k] * aur * dvol[k] : 0.0
        S = 0.0
        for j in (k, k + 1)
            int(j) || continue
            q = mp(k, j); Kraw -= u[j] * (q[1] + q[2]); p_ -= u[j] * q[3]; gr -= u[j] * q[4]; vi -= u[j] * q[5]; ra -= u[j] * q[6]; S += u[j] * sum(q)
            uM[i] += 0.5 * u[j] * Mrow[j+o]
            R6a[i] += 0.5 * u[j] * rm[j+o]; R6b[i] += 0.5 * u[j] * (sns[j+o] - ssp[j+o]); R6c[i] -= 0.5 * u[j] * dyn[j+o]
        end
        dΨ[i] = Ψ[i] - Ψ[i+1]; dΨr[i] = Ψr[i] - Ψr[i+1]
        du2 = (un[k] - uo[k])^2 + (un[k+1] - uo[k+1])^2
        K2a[i] = -(2θ - 1) * 0.25 * 0.5 * (mn(k) + mo(k)) * du2 / dt
        K2b[i] = (θ^2 - θ + 0.5) * 0.25 * (mn(k) - mo(k)) * du2 / dt
        K3a[i] = 0.5 * (dm[k] + dm[k+1]) * (u[k+1] - u[k]) * (0.5 * (u[k] + u[k+1]) - ub[k])
        Kbnd[i] = Kraw - K2a[i] - K2b[i] - K3a[i]
        prs[i] = p_; V4[i] = vi; G5[i] = gr
        Rbnd[i] = ra - (R6a[i] + R6b[i] + R6c[i] + dΨr[i])
        lrad && (old[i] = (c.ero[k] * c.dvolo[k] - mo(k) * (c.egro[k] - c.ego[k]) + mn(k) * (W.egrn[k] - W.egn[k]) - W.ern[k] * W.dvoln[k]) / dt)
        rem[i] = imb[i] - Gr[i] - Er[i] - (Kraw + p_ + vi + gr + ra + adv[i]) - S - old[i]
        dKs[i] = ksflux(k + 1, u, dm, a, b) - ksflux(k, u, dm, a, b)           # sub-zone kinetic-energy flux of the compatible row (telescopes)
        kst[i] = isempty(KSTAR[]) ? 0.0 : KSTAR[][k] / dt                      # BDF stage: K(y*) − (A K_γ − B K_n), taken out by the compatible row
        # ECOMPAT: total row − (compatible rows + u × momentum rows), gravity and advection forms apart
        cres[i] = imb[i] - gr - adv[i] - p_ - S - Gc[i] - Ec[i] - (ECOMPAT[] ? dKs[i] - kst[i] : 0.0)
    end
    (; imb, adv, kin, useint = ui[a:b], tiew = tw, rmom = rm, u = copy(W.u[a:b+1]), r = copy(W.rn[a:b+1]), E = copy(W.ern[a:b]), T = copy(W.tn[a:b]),
       F = copy(W.frn[a:b+1]), ρ = copy(W.dn[a:b]), θ, Gr, Er, uM, dΨ, dΨr, K2a, K2b, K3a, Kbnd, prs, V4, G5, R6a, R6b, R6c, Rbnd, rem, Kn, Ko, old, Gc, Ec, cres, dKs, kst, ecg,
       mchk, rchk, Mrow, ub = copy(ub[a:b]), dmdt = copy(dm[a:b+1]), qf = copy(qf[a:b]), pg = copy(W.pgn[a:b]), eg = copy(W.egn[a:b]),
       un = copy(un[a:b+1]), uo = copy(uo[a:b+1]), mn = [mn(k) for k in a:b], mo = [mo(k) for k in a:b], fedd = copy(fedd[a:b]),
       hedd = length(TD.HEDD) >= b ? copy(TD.HEDD[a:b]) : zeros(nz), hxedd = length(TD.HXEDD) >= b ? copy(TD.HXEDD[a:b]) : zeros(nz), urel = copy(W.urel[a:b+1]), db = copy(W.db[a:b+1]), xm = copy(W.xm[a:b+1]))
end

# ------------------------------------------------------------------ assembly
struct JTag end
# JTag is a hand-made tag (not a ForwardDiff.Tag), so ForwardDiff cannot order it against tags created by
# derivatives taken *inside* the residual (e.g. the Newton polish in the Z0 Saha solve). The residual's
# sweep is always the outer differentiation: JTag ≺ any ForwardDiff.Tag.
ForwardDiff.:≺(::Type{JTag}, ::Type{<:ForwardDiff.Tag}) = true
ForwardDiff.:≺(::Type{<:ForwardDiff.Tag}, ::Type{JTag}) = false
const DBUF = Dict{Any,Any}()

"Residual + coloured AD Jacobian -> rhs, em2..ep2 in COMMON; side outputs -> COMMON."
function assemble!(c::Ctx; ncol::Int = 6)
    neqn = c.neqn; nb = c.nb; k0 = c.k0; N = neqn * nb
    ncol * neqn == 54 ? _assemble!(c, Val(54), ncol) :
    ncol * neqn == 63 ? _assemble!(c, Val(63), ncol) :
    ncol * neqn == 42 ? _assemble!(c, Val(42), ncol) :
    ncol * neqn == 35 ? _assemble!(c, Val(35), ncol) :
    ncol * neqn == 25 ? _assemble!(c, Val(25), ncol) :
    ncol * neqn == 49 ? _assemble!(c, Val(49), ncol) : _assemble!(c, Val(ncol * neqn), ncol)
end
function _assemble!(c::Ctx, ::Val{PP}, ncol) where {PP}
    neqn = c.neqn; nb = c.nb; k0 = c.k0; N = neqn * nb
    D = ForwardDiff.Dual{JTag,Float64,PP}
    xd = get!(() -> Vector{D}(undef, N), DBUF, (:x, D, N))::Vector{D}
    Rd = get!(() -> Matrix{D}(undef, neqn, nb), DBUF, (:R, D, neqn, nb))::Matrix{D}
    @inbounds for kb in 1:nb, j in 1:neqn
        slot = mod(kb - 1, ncol) * neqn + j
        xd[(kb-1)*neqn+j] = D(0.0, ForwardDiff.Partials(ntuple(i -> i == slot ? 1.0 : 0.0, Val(PP))))
    end
    ws = getws(D, c.mgr)
    residual!(Rd, ws, c, xd)
    if GEQDBG2[] > 0                          # (GEQDBG=2: the first non-finite residual value or Jacobian entry, row and block, once per few calls)
        fnd = false
        @inbounds for kb in 1:nb, i in 1:neqn
            fnd && break
            vv = ForwardDiff.value(Rd[i, kb]); pp = ForwardDiff.partials(Rd[i, kb])
            if !isfinite(vv) || !all(isfinite, pp)
                GEQDBG2[] -= 1
                println("GEQDBG2: row ", i, " block ", kb, " (k ", c.k0 + kb - 1, ") value ", vv, " non-finite partials ", count(!isfinite, pp), " of ", length(pp))
                fnd = true
            end
        end
    end
    E = (v(:em2), v(:em1), v(:e00), v(:ep1), v(:ep2)); rhs = v(:rhs)
    _scatter!(E, rhs, Rd, neqn, nb, k0, ncol)
    store_side!(c, ws)
    nothing
end
function _scatter!(E, rhs, Rd::Matrix{D}, neqn, nb, k0, ncol) where {D}
    em2, em1, e00, ep1, ep2 = E
    @inbounds for kb in 1:nb, i in 1:neqn
        k = k0 + kb - 1
        rhs[i, k] = ForwardDiff.value(Rd[i, kb])
        pr = ForwardDiff.partials(Rd[i, kb])
        for o in -2:2
            M = o == -2 ? em2 : o == -1 ? em1 : o == 0 ? e00 : o == 1 ? ep1 : ep2
            cb = kb + o
            if 1 <= cb <= nb
                base = mod(cb - 1, ncol) * neqn
                for j in 1:neqn; M[i, j, k] = pr[base+j]; end
            else
                for j in 1:neqn; M[i, j, k] = 0.0; end
            end
        end
    end
    nothing
end
function store_side!(c::Ctx, ws::WS{D}) where {D}
    a, b = c.ngrs, c.ngre; W = ws.a
    put!(dst, src, rng) = (@inbounds for k in rng; dst[k] = ForwardDiff.value(src[k]); end)
    put!(v(:dsn), W.dsn, a:b+1); put!(v(:usn), W.usn, a:b)
    if c.lrad
        put!(v(:egrsn), W.egsn, a:b+1); put!(v(:ersn), W.ersn, a:b+1); put!(v(:frsn), W.frsn, a:b+1)
        put!(v(:egrb), W.egb, a+1:b+1); put!(v(:erb), W.erb, a+1:b+1); put!(v(:frb), W.frb, a:b)
        GASROW[] && (put!(v(:egsn), W.egsng, a:b+1); put!(v(:egb), W.egbg, a+1:b+1))
    else
        put!(v(:egsn), W.egsn, a:b+1); put!(v(:egb), W.egb, a+1:b+1)
    end
    put!(v(:qf), W.qf, a:b); put!(v(:dudr), W.dudr, a:b); v(:qv)[a:b] .= 0.0
    put!(v(:db), W.db, a+1:b+1); put!(v(:ub), W.ub, a:b)
    c.p.lgrid == 1 && put!(v(:xnu), W.xnu, a:b)
    nothing
end


"""
CLOSJ: the closure arrays of the formal solution as Duals: value + sum over i = m-1..m+1 of d(out_m)/dS_i (S_i^Dual - S_i) + d(out_m)/d(dtau_i)(dtau_i^Dual - dtau_i)
(TitanDriver.CSENS; S_i, dtau_i built from the iterate exactly as TitanCore.eddfac! builds them). Returns (fedd, HEDD, HXEDD, FXEDD, TIEW, geddr).
"""
function closure_lift(c, P, a, b, chifn, xkpn, plfn, dn, tn, ern, rn, xs1, xs2, ::Type{T}) where {T}
    cs = TD.CSENS; xnen = v(:xnen); csige = P.csige; cc = P.cc; cpi = P.cpi
    Sd = get!(() -> zeros(T, b + 2), DBUF, (:clos_Sd, T, b + 2))::Vector{T}; fill!(Sd, zero(T))
    Td = get!(() -> zeros(T, b + 2), DBUF, (:clos_Td, T, b + 2))::Vector{T}; fill!(Td, zero(T))
    z0_ = ZONE0[]; w0_ = CPLW[]; sl0_ = CPLSL[]
    COUPLED[] && (CPLW[] = (xs1, xs2)); CPLSL[] = nothing
    try
        if MATERIAL[] !== nothing && SLIMP[]
            loc = get(SLQBUF, (:locop, T, length(xs1)), nothing)
            if loc !== nothing && length(loc[1]) >= b
                chiL, kpL = loc::Tuple{Vector{T},Vector{T}}
            else
                chiL = get!(() -> zeros(T, b + 2), DBUF, (:clos_chiL, T, b + 2))::Vector{T}; fill!(chiL, zero(T))
                kpL = get!(() -> zeros(T, b + 2), DBUF, (:clos_kpL, T, b + 2))::Vector{T}; fill!(kpL, zero(T))
                zone_sweep!(a:b, _ -> true) do i
                    ch_, _, kp_ = material_opac(MATERIAL[], dn[i], tn[i])
                    @inbounds chiL[i] = ch_; @inbounds kpL[i] = kp_
                end
            end
        end
        @inbounds for i in a:b
            chi, kp = (MATERIAL[] !== nothing && SLIMP[]) ? (chiL[i], kpL[i]) : (chifn[i], xkpn[i])
            pl = plfn[i]; ρ = dn[i]; E = ern[i]
            ne = xnen[i] * (ρ / val(ρ)) * (cplz(i) && val(xs1[i]) > 0 ? xs1[i] / val(xs1[i]) : one(T))
            S = if cs.msrc
                (cs.scat && chi * ρ > kp * ρ + csige * ne) ? (kp * pl + (chi - kp) * E * cc / (4cpi)) / chi :
                                                              (kp * ρ * pl + csige * ne * E * cc / (4cpi)) / (kp * ρ + csige * ne)
            else
                ((chi * ρ - csige * ne) * pl + csige * ne * E * cc / (4cpi)) / (chi * ρ)
            end
            τ = chi * ρ * (rn[i+1] - rn[i])
            Sd[i] = S - val(S); Td[i] = τ - val(τ)
        end
    finally
        ZONE0[] = z0_; CPLW[] = w0_; CPLSL[] = sl0_
    end
    function lift(name::Symbol, base::AbstractVector)
        out = get!(() -> Vector{T}(undef, length(base)), DBUF, (:clos_out, name, T, length(base)))::Vector{T}
        @inbounds for k in eachindex(base); out[k] = T(base[k]); end
        haskey(cs.S, name) || return out
        MS = cs.S[name]; MT = cs.T[name]
        @inbounds for m in a:min(b, length(base))
            acc = zero(T)
            for i in max(a, m - 1):min(b, m + 1); acc += MS[m, i - m + 2] * Sd[i] + MT[m, i - m + 2] * Td[i]; end
            out[m] += acc
        end
        out
    end
    g = T(P.geddr)
    @inbounds for i in max(a, b - 1):b; g += cs.gS[i - (b - 1) + 1] * Sd[i] + cs.gT[i - (b - 1) + 1] * Td[i]; end
    (lift(:f, c.fedd), lift(:h, TD.HEDD), lift(:hx, TD.HXEDD), lift(:fx, TD.FXEDD), lift(:w, TD.TIEW), g)
end

const RESID_ONLY = Ref(false)
function matgen!(; ncol = 6)
    if RESID_ONLY[]
        c = context(); neqn, nb, k0 = c.neqn, c.nb, c.k0
        ws = getws(Float64, c.mgr); R = get!(() -> zeros(neqn, nb), DBUF, (:Rf, neqn, nb))::Matrix{Float64}
        xz = get!(() -> zeros(neqn * nb), DBUF, (:xzf, neqn * nb))::Vector{Float64}
        residual!(R, ws, c, xz); copyto!(view(v(:rhs), 1:neqn, k0:k0+nb-1), R)
        return nothing
    end
    TD.CLOSJ[] && TD.CSENS_HOOK[] !== nothing && gi(:lrad) == 1 && (!TD.CSENS.valid || haskey(ENV, "NKCHECK")) && TD.CSENS_HOOK[]()
    assemble!(context(); ncol)
    nothing
end

"""
COUPLED test (popiii_m2c CPLJCHECK=n1,n2,…: zone numbers from the inner boundary, at the first Newton iteration of the first step): the
AD Jacobian blocks of the energy and species rows of zones k−2…k+2 with respect to T, x_H⁺, x_H₂ of zone k, against central finite
differences of the Float residual (step h in the logarithmic variables). Prints the worst relative difference per column.
"""
function cpl_jcheck(kz; h = parse(Float64, get(ENV, "CPLJCHH", "1e-6")))
    c = context(); assemble!(c)
    E = (v(:em2), v(:em1), v(:e00), v(:ep1), v(:ep2)); neqn = c.neqn; nb = c.nb; k0 = c.k0; it = c.ix[5]
    ws = getws(Float64, c.mgr); R1 = zeros(neqn, nb); R2 = zeros(neqn, nb); xz = zeros(neqn * nb)
    kbz = kz - k0 + 1
    for j in (it, 8, 9)
        idx = (kbz - 1) * neqn + j
        xz[idx] = h; residual!(R1, ws, c, xz); xz[idx] = -h; residual!(R2, ws, c, xz); xz[idx] = 0.0
        worst = 0.0; wl = ""
        for kr in kz-2:kz+2, i in (it, 8, 9)
            kbr = kr - k0 + 1; (1 <= kbr <= nb) || continue
            fd = (R1[i, kbr] - R2[i, kbr]) / 2h; ad = E[kz - kr + 3][i, j, kr]
            sc = max(abs(fd), abs(ad)); sc == 0 && continue
            rel = abs(fd - ad) / sc
            # relative to the row's largest Jacobian entry (a tiny entry may differ by round-off)
            rowmax = maximum(abs(E[o][i, jj, kr]) for o in 1:5 for jj in 1:neqn)
            rel = abs(fd - ad) / max(sc, 1e-6 * rowmax)
            rel > worst && (worst = rel; wl = "row $i zone $kr: AD $(round(ad, sigdigits = 6)) FD $(round(fd, sigdigits = 6))")
        end
        println("CPLJCHECK zone $kz column $j (cplz $(cplz(kz))): worst relative difference $(round(worst, sigdigits = 3))  ($wl)")
    end
    nothing
end

"""
SZJCHECK test (popiii_m2c SZJCHECK=n1,n2,…, at the first Newton iteration of the first step): the AD Jacobian of the gas-energy,
radiation-energy and radiation-flux rows of zones k−2…k+2 with respect to ρ, T, E_r, F_r (and x_H⁺, x_H₂ if COUPLED) of zone k,
against central finite differences of the value residual (step h in the scaled variables). Prints the worst relative difference
per column (relative to max(|AD|, |FD|, 1e-6 × the row's largest entry)). Written for SUBZONE, whose avchi enters the flux row.
"""
function rad_jcheck(kz; h = parse(Float64, get(ENV, "SZJCHH", "1e-6")))
    c = context(); assemble!(c)
    E = (v(:em2), v(:em1), v(:e00), v(:ep1), v(:ep2)); neqn = c.neqn; nb = c.nb; k0 = c.k0
    ir, im, id, iu, it, ie, iF = c.ix
    ws = getws(Float64, c.mgr); R1 = zeros(neqn, nb); R2 = zeros(neqn, nb); xz = zeros(neqn * nb)
    kbz = kz - k0 + 1
    cols = neqn >= 9 ? (id, it, ie, iF, 8, 9) : (id, it, ie, iF)
    for j in cols
        idx = (kbz - 1) * neqn + j
        xz[idx] = h; residual!(R1, ws, c, xz); xz[idx] = -h; residual!(R2, ws, c, xz); xz[idx] = 0.0
        worst = 0.0; wl = ""
        for kr in kz-2:kz+2, i in (neqn >= 9 ? (it, ie, iF, 8, 9) : (it, ie, iF))
            kbr = kr - k0 + 1; (1 <= kbr <= nb) || continue
            fd = (R1[i, kbr] - R2[i, kbr]) / 2h; ad = E[kz - kr + 3][i, j, kr]
            rowmax = maximum(abs(E[o][i, jj, kr]) for o in 1:5 for jj in 1:neqn)
            sc = max(abs(fd), abs(ad), 1e-6 * rowmax); sc == 0 && continue
            rel = abs(fd - ad) / sc
            rel > worst && (worst = rel; wl = "row $i zone $kr: AD $(round(ad, sigdigits = 6)) FD $(round(fd, sigdigits = 6))")
        end
        println("SZJCHECK zone $kz column $j (SUBZONE $(SUBZONE[]), SLIMP $(SLIMP[])): worst relative difference $(round(worst, sigdigits = 3))  ($wl)")
    end
    nothing
end

end # module
