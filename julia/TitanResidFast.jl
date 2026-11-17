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
    egso::Vector{Float64}; egrso::Vector{Float64}; dso::Vector{Float64}; uso::Vector{Float64}
    erso::Vector{Float64}; frso::Vector{Float64}; pgo::Vector{Float64}; aso::Vector{Float64}
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
k2 (the outer edge) stay; every interface in between obeys n̂_j = n̂_{j−1} with n̂_j = x̃_j / sqrt(1 + G_j²), x̃_j = ν_j −
α(1+α) Δ²ν (ν_j = r̄_j/Δr_j, the logarithmic baseline; α = alph), time-smoothed per zone, x̃ → (1 + γ_j) x̃ − γ_j x̃_old with
γ_j = tauv[j]/Δt (tauv: a multiple of the zone's signal-crossing time, frozen for the step), and G_j² = ν_j² S_j capped smoothly at
gmax², S_j = Σ_q w_q (½ Δ ln q)² over q = ρ, P, T, E_r and the interior mass m (centred differences of the zone values; one-sided
at the ends) plus sx[j] (terms frozen for the step: the capped ln τ). All terms are logarithmic, so the rule has no scale.
"""
const GRIDEQ = Ref{Any}(nothing)
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
SLIMP (popiii_m2c `SLIMP=1`, with COUPLED; default off, off = bit-identical): the star-light dilution of a coupled zone at the Newton
iterate. `SLQ_FN[](k, E_k, F_k, F_k+1)` gives (q_k, t_k) = (W_k, T_c,k at the iterate)/(their start-of-step values) (local: zone k's E_r and its
two interface fluxes), `CPL_ZONEQ[](k, ρ, T, x1, x2, q, t)` the zone's terms at that star light (PopIIIMaterial.cpl_zone_q), and `CPLQ[]`
carries (q, t) (Duals) to the NLOIT opacity during the material sweep. SLIMP=1: W only (t ≡ 1); SLIMP=2: also the local part of T_c.
"""
const SLON_FN = Ref{Any}(k -> true)   # (SLIMP) zone k's star light follows the iterate (else the SLIMP=0 path)
const SLIMP = Ref(false); const SLQ_FN = Ref{Any}(nothing); const CPL_ZONEQ = Ref{Any}(nothing); const CPLQ = Ref{Any}(nothing)
const SLIMPDBG = Ref(parse(Int, get(ENV, "SLIMPDBG", "0")))   # (diagnostic) report the first n non-finite SLIMP values in the residual
const CPLQT = Ref{Any}(nothing)   # (SLIMP) t = T_c/T_c,start of the last residual's material sweep, for the coupled-zone sweep that follows
const SLQBUF = Dict{Any,Any}()
const CPL_LOCAL = Ref{Any}(nothing)             # (k, dt) -> (x1, x2): the zone's species rows solved alone (initial guess), set by the material
"Weight floor of the species corrections in maxdel's convergence measure, and the clamp of a species correction per Newton iteration."
const CPL_XFLOOR = Ref(1e-8); const CPL_CLAMP = Ref(parse(Float64, get(ENV, "CPLCLAMP", "1")))
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
            (CPL_NDECOUPLED[] <= 30 || CPL_NDECOUPLED[] % 100 == 0) && println("COUPLED: zone $(k - a + 1) decoupled for this step (species correction > $(CPL_BADDELTA[]) in ln x), $(CPL_NDECOUPLED[]) so far")
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
"maxdel: clamp the species corrections to ±CPL_CLAMP and return their weighted maximum (convergence measure)."
function cpl_dmax!(rhs, a, b)
    d = 0.0; cl = CPL_CLAMP[]; f = CPL_XFLOOR[]; kx = TF.var(:kx); dx = TF.var(:dx)
    for i in 8:9; kx[i] = 0; dx[i] = 0.0; end
    @inbounds for k in a:b
        cplz(k) || continue
        for (i, q) in ((8, XSN[][1]), (9, XSN[][2]))
            abs(rhs[i, k]) > abs(dx[i]) && (kx[i] = k - a + 1; dx[i] = rhs[i, k])     # (raw correction, for the divergence message)
            isfinite(rhs[i, k]) || (rhs[i, k] = 0.0; d = Inf)                        # a non-finite correction: no change, not converged
            abs(rhs[i, k]) > CPL_BADDELTA[] && push!(CPL_BAD, k)
            rhs[i, k] = clamp(rhs[i, k], -cl, cl)
            # weight by the larger of the iterate and the old abundance: a species the iterate drove to the floor (x_H⁺ 1e-40 next to the
            # shock, from 1e-3) is not converged; only species that are trace at both ends are left out
            xm_ = max(q[k], XSO[][i - 7][k]); d = max(d, abs(rhs[i, k]) * xm_ / (xm_ + f))
        end
    end
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
    ui
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
        v(:fro), v(:egro), v(:egso), v(:egrso), v(:dso), v(:uso),
        v(:erso), v(:frso), v(:pgo), v(:aso), v(:chifo), v(:xkeo),
        v(:xkpo), v(:plfo), v(:dvolo), v(:xmo), v(:fedd), v(:bri),
        (la = v(:lady); lb = get!(() -> zeros(Int, length(la)), BUFS, :lady)::Vector{Int}; lb .= la; lb), v(:wt), v(:yscl), P, etotal, dualswitch(ngrs, ngre, lrad, ηsw),
        CONVECTION[], ALPHA_MLT[], convcoef(ngrs, ngre, base))
end

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
    # as before. Bit-identical.
    χτ = get!(() -> zeros(mgr), BUFS, :cvchi)::Vector{Float64}; okτ = get!(() -> zeros(Bool, mgr), BUFS, :cvchiok)::Vector{Bool}; fill!(okτ, false)
    if opac_pure(m)
        z0 = ZONE0[]
        try zone_sweep!(ngrs+1:ngre, k -> @inbounds(base[k, 3] > 0 && base[k, 5] > 0)) do k
            ZONE[] = z0
            try
                y = @inbounds material_opac(m, base[k, 3], base[k, 5])[1]
                y isa Float64 && (@inbounds χτ[k] = y; @inbounds okτ[k] = true)
            catch
            end
        end
        finally ZONE0[] = z0 end
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
@inline function decide(kind::Symbol, id::Int, k::Int, x::Bool)
    m = FREEZE[]
    m == 0 && return x
    v = get!(() -> fill(false, TF.MGR), DEC, (kind, id))
    m == 1 && (@inbounds v[k] = x; return x)
    @inbounds return v[k]
end

"""
EXCAP[] > 0: factor on the gas-radiation exchange of the two energy rows, so that c κ_E ρ Δt never exceeds EXCAP[]. In an opaque zone that product can be
1e20 and the rows carry it as round-off; E = aT⁴ there holds to about 1/EXCAP anyway. The same factor in both rows keeps energy conserved. A constant of the iteration (no derivative).
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
    return min(1.0, EXCAP[] / (cc * val(κ) * val(ρ) * val(dt)))
end
include(joinpath(@__DIR__, "MGRad.jl"))                 # MGROUPS (default off)
"""
(MGROUPS, stage 2c) The star light of a coupled zone is its group fields at the Newton iterate. `MGSLQ[](Mg, k, E_k)` fills Mg.qs[1:G, k] with
c u_g/c u_g(start of step) and returns (q of group MGSLG, T_c/T_c,start) (PopIIIMaterial.sl_qt_mg); `MGH2FN[](Mg, dn, tn, rn, un, x_H2, t, a, b)`
the H₂ Lyman–Werner closure of the iterate (Mg.kh2, Mg.kh2g, Mg.qs[G+1, :]; PopIIIMaterial.mg_h2iter!); `MGZONEQ[](k, ρ, T, x1, x2, qv, t)` the
zone's terms in those fields (PopIIIMaterial.cpl_zone_qg). `CPLQ[]` carries (q, t, Mg.qs) to the NLOIT opacity during the material sweep.
"""
const MGSLQ = Ref{Any}(nothing); const MGSLG = Ref(2); const MGH2FN = Ref{Any}(nothing); const MGZONEQ = Ref{Any}(nothing)

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
const MAT_THREADS = Ref(parse(Int, get(ENV, "MATTHREADS", "0")))
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
"MGROUPS: `material_opac`, and the zone's group record (MGRad.jl) into `o` in the number type of (ρ, T). Default: the grey record."
material_opacg!(o, m, ρ, T) = (mg_default!(o); material_opac(m, ρ, T))
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
            χ, κE, κP = MG[] ? material_opacg!(mg_gcol(eltype(W.chifn), k), m, dn[k], tn[k]) : material_opac(m, dn[k], tn[k])
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
function _material!(m, W, c, P, dn, tn, ern, k0, ngre)
    hassrc = has_source(m)
    nt = min(MAT_THREADS[], Threads.nthreads())
    (nt > 1 && ngre + 3 - k0 >= 4nt && !ZTASK[]) && return _material_mt!(m, W, c, P, dn, tn, ern, k0, ngre, hassrc, nt)
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
    kerr = fill(typemax(Int), nt); errs = Vector{Any}(nothing, nt)
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
    kmin, imin = findmin(kerr)
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
    nt = ZONESWEEP[] ? min(MAT_THREADS[], Threads.nthreads()) : 1
    if !(nt > 1 && length(ks) >= 4nt && !ZTASK[])
        for k in ks
            ok(k) || continue
            ZONE[] = k; f(k)
        end
        return nothing
    end
    # (FIX-RACE) the first zone serially, before any task starts: lazily built shared tables (photo-rate tables, compiled networks, ...)
    # that every zone needs are then built by one task, not raced for (the restart z100s50b crashed with UndefRefError in the first step,
    # where cpl_prepare!'s first photo-rate evaluations of the process ran on the threads). Same order of values: bit-identical.
    k1 = first(ks); while k1 <= last(ks) && !ok(k1); k1 += 1; end
    k1 > last(ks) && return nothing
    ZONE[] = k1; f(k1)
    next = Threads.Atomic{Int}(k1 + 1); kend = last(ks)
    kerr = fill(typemax(Int), nt); errs = Vector{Any}(nothing, nt)
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
    kmin, imin = findmin(kerr)
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
function grideq_rows!(R, ir, k0, rn, ro, dn, tn, pg, er, xmen, xmtot, dt, ge, chi = nothing)
    k1, k2 = ge.k1, ge.k2; zl, zh = k1, k2 - 1                  # zones between the fixed interfaces
    α = ge.alph * (1 + ge.alph); wρ, wP, wT, wE, wm = ge.w
    zc(j) = clamp(j, zl, zh)
    ν(j) = (j = zc(j); (rn[j+1] + rn[j]) / (2 * (rn[j+1] - rn[j]))); νo(j) = (j = zc(j); (ro[j+1] + ro[j]) / (2 * (ro[j+1] - ro[j])))
    lm(j) = (j = zc(j); log(0.5 * ((xmtot - xmen[j]) + (xmtot - xmen[j+1]))))
    function nhat(j)
        xnt = ν(j) - α * (ν(j-1) - 2ν(j) + ν(j+1)); xnto = νo(j) - α * (νo(j-1) - 2νo(j) + νo(j+1))
        jm, jp = zc(j - 1), zc(j + 1); h = jp - jm == 2 ? 0.5 : 1.0
        S = wρ * (h * log(dn[jp] / dn[jm]))^2 + wP * (h * log(pg[jp] / pg[jm]))^2 + wT * (h * log(tn[jp] / tn[jm]))^2 +
            wE * (h * log(er[jp] / er[jm]))^2 + wm * (h * (lm(jp) - lm(jm)))^2 + ge.sx[j]
        # (GEQW entry 8, wχ > 0: the extinction coefficient χρ of the Newton iterate, ln(χ_jp ρ_jp/χ_jm ρ_jm): the opacity cliff, local, not lagged)
        (hasproperty(ge, :wχ) && ge.wχ > 0 && chi !== nothing) && (S += ge.wχ * (h * log((chi[jp] * dn[jp]) / (chi[jm] * dn[jm])))^2)
        G2 = ν(j)^2 * S; G2 = G2 * ge.gmax^2 / (G2 + ge.gmax^2)
        γ = ge.tauv[j] / dt
        ((1 + γ) * xnt - γ * xnto) / sqrt(1 + G2)
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
        slq = nothing
        if SLIMP[]                                   # SLIMP: q_k of the iterate (new-level E_r, F_r)
            slq = get!(() -> Vector{T}(undef, length(xs1)), SLQBUF, (T, length(xs1), 1))::Vector{T}
            slt = get!(() -> Vector{T}(undef, length(xs1)), SLQBUF, (T, length(xs1), 2))::Vector{T}
            @inbounds for k in eachindex(slq); slq[k] = one(T); slt[k] = one(T); end
            if Mg === nothing
                @inbounds for k in ngrs:ngre; (cplz(k) && SLON_FN[](k)) && ((slq[k], slt[k]) = SLQ_FN[](k, ern[k], frn[k], frn[k+1])); end
            else                                     # MGROUPS: every group's field at the iterate against its start-of-step value, and the H₂ LW closure
                @inbounds for k in ngrs:ngre; (cplz(k) && SLON_FN[](k)) && ((slq[k], slt[k]) = MGSLQ[](Mg, k, ern[k])); end
                MGH2FN[] === nothing || MGH2FN[](Mg, dn, tn, rn, un, xs2, slt, ngrs, ngre)
            end
            if SLIMPDBG[] > 0
                fin(v) = isfinite(ForwardDiff.value(v)) && (!(v isa ForwardDiff.Dual) || all(isfinite, ForwardDiff.partials(v)))
                for k in ngrs:ngre
                    (fin(slq[k]) && fin(slt[k])) && continue
                    SLIMPDBG[] -= 1; println("SLIMPDBG q/t non-finite zone $(k - ngrs + 1): q $(slq[k]) t $(slt[k]) E $(ern[k]) F $(frn[k]) $(frn[k+1])"); break
                end
            end
            CPLQ[] = Mg === nothing ? (slq, slt) : (slq, slt, Mg.qs)
        end
        try
            MATERIAL[] === nothing || _material!(MATERIAL[], W, c, P, dn, tn, ern, k0, ngre)
        finally
            CPLW[] = nothing; slq === nothing || (CPLQT[] = CPLQ[][2]); CPLQ[] = nothing
        end
        # (PERF512) the coupled zones' rates in the Jacobian's Dual numbers on MATTHREADS tasks (zones independent; the material sweep above
        # left their values and local Jacobians in the zones' cache entries); bit-identical. ZONE[] is not used here and is kept.
        let z0 = ZONE0[]
            try zone_sweep!(ngrs:ngre, k -> true) do k
                @inbounds if cplz(k)
                    y = (slq === nothing || !SLON_FN[](k)) ? CPL_ZONE[](k, dn[k], tn[k], xs1[k], xs2[k]) :
                        Mg === nothing ? CPL_ZONEQ[](k, dn[k], tn[k], xs1[k], xs2[k], slq[k], CPLQT[][k]) :
                        MGZONEQ[](k, dn[k], tn[k], xs1[k], xs2[k], view(Mg.qs, :, k), CPLQT[][k])
                    W.cpS1[k] = y[3]; W.cpS2[k] = y[4]; W.cpL[k] = y[5]
                    if SLIMPDBG[] > 0 && !all(v -> isfinite(ForwardDiff.value(v)) && (!(v isa ForwardDiff.Dual) || all(isfinite, ForwardDiff.partials(v))), y)
                        SLIMPDBG[] -= 1; println("SLIMPDBG cpl zone $(k - ngrs + 1) non-finite: y ", y, " q ", slq === nothing ? 1 : slq[k], " t ", slq === nothing ? 1 : CPLQT[][k])
                    end
                    W.cpQ[k] = y[6]; W.cpY[k] = y[7]
                    Mg === nothing || (Mg.lcl[k] = y[8])
                else
                    W.cpS1[k] = zero(T); W.cpS2[k] = zero(T); W.cpL[k] = zero(T); W.cpQ[k] = zero(T); W.cpY[k] = zero(T)
                    Mg === nothing || (Mg.lcl[k] = zero(T))
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
    advect!(db, dsn, dn, c.do_, c.dso, fl, θ, P.cadv, P.eps, ngrs, ngre, true, 1)
    @inbounds for k in ngrs+1:ngre+1; fle[k] = -dmdt[k]; end
    @inbounds for k in ngrs:ngre; flu[k] = -(dmdt[k] + dmdt[k+1]); end
    advect!(ub, usn, un, c.uo, c.uso, flu, θ, P.cadv, P.eps, ngrs, ngre, false, 2)
    if c.lrad
        advect!(egb, egsn, egrn, c.egro, c.egrso, fle, θ, P.cadv, P.eps, ngrs, ngre, true, 3)
        GASROW[] && advect!(W.egbg, W.egsng, egn, c.ego, c.egso, fle, θ, P.cadv, P.eps, ngrs, ngre, true, 6)
        qe, qeo, qF, qFo = W.qe, W.qeo, W.qF, W.qFo
        @inbounds for k in ngrs-1:ngre+2
            qe[k] = ern[k] / dn[k]; qeo[k] = c.ero[k] / c.do_[k]
            qF[k] = frn[k] * 2 / (dn[k] + dn[k-1]); qFo[k] = c.fro[k] * 2 / (c.do_[k] + c.do_[k-1])
        end
        advect!(W.erb, W.ersn, qe, qeo, c.erso, fle, θ, P.cadv, P.eps, ngrs, ngre, true, 4)
        advect!(W.frb, W.frsn, qF, qFo, c.frso, flu, θ, P.cadv, P.eps, ngrs, ngre, false, 5)
    else
        advect!(egb, egsn, egn, c.ego, c.egso, fle, θ, P.cadv, P.eps, ngrs, ngre, true, 6)
    end
    erb, frb = W.erb, W.frb
    MG[] && mg_theta!(Mg, c, W, θ, ngrs, ngre, fle, flu)     # MGROUPS: θ-level group opacities, emission, χ̄_g, face values
    if COUPLED[]
        advect!(W.xb1, W.xsn1, W.xs1, XSO[][1], XSSO[][1], fle, θ, P.cadv, P.eps, ngrs, ngre, true, 7)
        advect!(W.xb2, W.xsn2, W.xs2, XSO[][2], XSSO[][2], fle, θ, P.cadv, P.eps, ngrs, ngre, true, 8)
    end
    # diffusion
    dfm, dfe = W.dfm, W.dfe
    @inbounds for k in ngrs:ngre+1
        dfm[k] = P.sigd == 0 ? zero(T) : 2P.sigd * rmu[k] * (d[k] - d[k-1]) / (r[k+1] - r[k-1])
        dfe[k] = P.sige == 0 ? zero(T) : P.sige * rmu[k] * (d[k] + d[k-1]) * (eg[k] - eg[k-1]) / (r[k+1] - r[k-1])
    end
    fedd = c.fedd
    fill!(R, zero(T))
    # ---------------------------------------------------------------- GRID
    kb(k) = k - k0 + 1
    X(i, k) = x[(kb(k) - 1) * neqn + i]
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
        GRIDEQ[] === nothing || grideq_rows!(R, ir, k0, rn, c.ro, dn, tn, W.pgn, ern, xmen, P.xmtot, dt, GRIDEQ[], W.chifn)
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
        w = (un[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
             c.uo[k] * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2dt) -
            0.5 * ((dmdt[k] + dmdt[k+1]) * ub[k] - (dmdt[k-1] + dmdt[k]) * ub[k-1]) +
            rmu[k] * (pg[k] - pg[k-1]) +
            (2xmu * P.cpi * P.cgrav * xm[k] / rmu[k] + (1 - 0.5xmu) * P.g) * mk +
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
            Pm = pgn[k-1] + c.fedd[k-1] * ern[k-1]; Pp = pgn[k] + c.fedd[k] * ern[k]
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
                lph = MLT_HSDISC[] ? -gk * 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) / (rmu[k] * Pi) : -0.5 * (rn[k+1] - rn[k-1]) / HP
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
        if c.etotal && !c.useint[k]
            w = etot_row(k, c, W, T) + rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] - src
            isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])     # the gas gets the ray's energy in chemistry: it leaves the total here (as the lrad total row)
        elseif c.lrad && GASROW[]
            egbg = W.egbg
            w = (dn[k] * egn[k] * dvoln[k] - c.do_[k] * c.ego[k] * c.dvolo[k]) / dt -
                dfe[k+1] + dfe[k] + rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] +
                pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * dudr[k] * dvol[k] -
                (Mg === nothing ? excap(cc, xke[k], d[k], dt) * (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) : Mg.xch[k]) * d[k] * dvol[k] +
                compton(k, W, P) * dvol[k] - src                   # (MGROUPS: Σ_g of the group exchange)
            isempty(GASQ[]) || (w -= GASQ[][k] * dvol[k])
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
        @inbounds for kb_ in 1:nb, i in 8:9; R[i, kb_] = x[(kb_-1)*neqn+i]; end
        # no radiation (lrad = 0, the collapse stages): rows 6–7 (E_r, F_r) are not formed above, but neqn = 9 puts the species at 8–9;
        # make them trivial (δ = 0) instead of leaving the buffer's stale values (a singular block: all Newton corrections 0, step 1 fails)
        c.lrad || @inbounds for kb_ in 1:nb, i in 6:7; R[i, kb_] = x[(kb_-1)*neqn+i]; end
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
            w = (ern[k] * dvoln[k] - c.ero[k] * c.dvolo[k]) / dt + rmu[k+1] * fr[k+1] - rmu[k] * fr[k] +
                fedd[k] * er[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur * dvol[k] +
                excap(cc, xke[k], d[k], dt) * (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) * d[k] * dvol[k] - compton(k, W, P) * dvol[k]
            isempty(RAYQ[]) || (w += RAYQ[][k] * dvol[k])
            isempty(GASQ[]) || (w += GASQ[][k] * dvol[k])
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
                (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                 qFo * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
                ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2) +
                rmu[k] * (fedd[k] * er[k] - fedd[k-1] * er[k-1]) +
                0.5 * rmu[k] * (fr[k] / cc^2) * (u[k+1] - u[k-1]) +
                (0.25xmu / r[k]) * ((3fedd[k-1] - 1) * er[k-1] * dvol[k-1] + (3fedd[k] - 1) * er[k] * dvol[k]) +
                0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] * fr[k] / cc
            else                                     # VEF_SPH: static terms in the sphericality-factor form (TitanDriver.VEF_SPH)
                qs = vefq(fedd, r, xmu, k)
                R[iF, kb(k)] =
                    (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                     qFo * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
                    ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2) +
                    rmu[k] * (qs * fedd[k] * er[k] - fedd[k-1] * er[k-1] / qs) +
                    0.5 * rmu[k] * (fr[k] / cc^2) * (u[k+1] - u[k-1]) +
                    0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] * fr[k] / cc
            end
            wt = tie ? TD.TIEW[k] : 0.0
            if wt > 0                                # FLUXTIE: F = h c E (TitanDriver.FLUXTIE)
                rcm = (0.5 * (r[k-1] + r[k]))^mu; rcp = (0.5 * (r[k] + r[k+1]))^mu
                # upwind (donor cell): outward light from the zone inside, inward from the zone outside, diluted as r^-2. The centred
                # average ½(hE_(k−1) + hE_k) cannot see an odd–even E mode, which nothing damps in thin gas (E → 0 at zone 214, step 87126)
                if !SLABTIE[]
                B = rmu[k] * fr[k] / cc - (rcm * max(TD.HEDD[k-1], 0.0) * er[k-1] + rcp * min(TD.HEDD[k], 0.0) * er[k])
                else                                 # SLABTIE: outward donor h from the slab depth of the donor's centre below the interface
                    hd = st_hdonor(max(TD.HEDD[k-1], 0.0), fedd[k-1], chif[k-1], chif[k], 0.5 * d[k-1] * dvol[k-1] / rmu[k], 0.5 * d[k] * dvol[k] / rmu[k])
                    B = rmu[k] * fr[k] / cc - (rcm * hd * er[k-1] + rcp * min(TD.HEDD[k], 0.0) * er[k])
                end
                R[iF, kb(k)] = (1 - wt) * R[iF, kb(k)] + wt * B
            end
        end
        P.lribc in (2, 3) || error("lribc must be 2 or 3")
        R[iF, kb(ngrs)] = X(iF, ngrs)
        R[iF, kb(ngre+1)] = P.cpi * (2P.geddr + 1) * P.ximr + frn[ngre+1] - cc * P.geddr * (sph ? vef_esurf(W, fedd, ngre, xmu, cc, 4P.cpi * P.ximr / cc) : ern[ngre])
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
r_c,N → r_N+1 with its own Eddington factor f_x = (f_N E_N − E_b/3)/E_x (K = K_x + K_b), held over the half zone and clamped to
[1/3, 1] (where E_x ≪ E_b, f_x is the difference of two nearly equal numbers; the clamp bounds the correction to the dilution of a
beam, e^{−2Δ}, so that the quadrature error of f, ~1e-3, is not multiplied by E_b/E_x ~ 1e5 in the boundary flux):
E_s = E_b + e^{−s} E_x − e^{−s/2} χ_N ρ_N F_N+1 (r_N+1 − r_c,N)/(c f_x), s = (μ/2)(3 − 1/f_x) ln(r_N+1/r_c,N).
"""
@inline function vef_esurf(W, fedd, N, xmu, cc, Eb)
    rn = W.rn; rc = 0.5 * (rn[N] + rn[N+1]); Ex = W.ern[N] - Eb
    fx = ForwardDiff.value(Ex) != 0 ? clamp(ForwardDiff.value((fedd[N] * W.ern[N] - Eb / 3) / Ex), 1 / 3, 1.0) : 1 / 3
    s = 0.5xmu * (3 - 1 / fx) * log(rn[N+1] / rc)
    Eb + exp(-s) * Ex - exp(-0.5s) * W.chifn[N] * W.dn[N] * W.frn[N+1] * (rn[N+1] - rc) / (cc * fx)
end
# ------------------------------------------------------------------ total energy (item 2)
"Total pressure (gas + radiation fE) of zone k at the θ-level."
@inline ptot(k, c, W) = c.lrad ? W.pg[k] + c.fedd[k] * W.er[k] : W.pg[k]
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
@inline function eflux(k, c, W, ::Type{T}) where {T}
    P = c.p; ngrs, ngre = c.ngrs, c.ngre; xmu = Float64(c.mu)
    u = W.u[k]; rmu = W.rmu[k]; dmdt = W.dmdt[k]
    Fl = zero(T)
    c.lrad && (Fl += rmu * W.fr[k])
    Fl -= W.dfe[k]
    if k == ngrs
        Fl += c.lrad ? (P.phil2 + W.er[k] * W.urel[k]) * rmu : P.phil2 * rmu
        Fl += -dmdt * 0.5 * u^2 + rmu * u * ptot(k, c, W) + rmu * u * W.qf[k]
    elseif k == ngre + 1
        if P.leobc != 3
            Fl += c.lrad ? rmu * (P.phir2 + W.urel[k] * W.er[k-1]) : rmu * P.phir2
        else
            Fl += -dmdt * W.egb[k]
        end
        Fl += -dmdt * 0.5 * u^2 + rmu * u * ptot(k - 1, c, W) + rmu * u * W.qf[k-1]
    else
        Fl += -dmdt * (W.egb[k] + 0.5 * u^2) + rmu * u * 0.5 * (ptot(k - 1, c, W) + ptot(k, c, W)) +
              u * W.r[k]^(-0.5xmu) * 0.5 * (r3(k - 1, W, xmu) * W.qf[k-1] + r3(k, W, xmu) * W.qf[k])
    end
    Fl
end
"ETOT_GRAV[]: the total-energy row includes the work of gravity (self-gravity or P.g). Default false, as in the tests and ledgers that book the potential energy themselves."
const ETOT_GRAV = Ref(false)
"ETOT_GRAV_GRID[] (ENV ETOTGRAVGRID=1; default off): with ETOT_GRAV[], also the work of gravity on mass carried by moving faces. Without it the total-energy row of a Lagrangian zone has no gravity work at all, and the kinetic energy gravity gives the gas is taken from its thermal energy (cloud-phase test of 2026-11-09)."
const ETOT_GRAV_GRID = Ref(false)
"Conservative total-energy row of zone k: Δ[m(e_g [+E_r/ρ]) + m(u_k²+u_{k+1}²)/4]/dt + Φ_{k+1} − Φ_k."
@inline function etot_row(k, c, W, ::Type{T}) where {T}
    dn, un, dvoln = W.dn, W.un, W.dvoln
    en = c.lrad ? W.egrn[k] : W.egn[k]; eo = c.lrad ? c.egro[k] : c.ego[k]
    w = (dn[k] * dvoln[k] * (en + 0.25 * (un[k]^2 + un[k+1]^2)) -
         c.do_[k] * c.dvolo[k] * (eo + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))) / c.p.dt +
        eflux(k + 1, c, W, T) - eflux(k, c, W, T)
    if ETOT_GRAV[]
        # Work of gravity on the mass that actually crosses each face (the continuity row's flux r² u ρ_b), over half a zone on each side: the rows sum to
        # Σ_faces (mass flux) g Δr, the change of Σ m Φ, so gas + kinetic + potential energy is conserved. The momentum row's force g m_k with m_k the mean of the
        # neighbouring zones differs from this where ρ is steep across a face (the accretion shock on a fixed mesh), and that difference must not create energy.
        xmu = Float64(c.mu); P = c.p
        gi_(j) = 2xmu * P.cpi * P.cgrav * W.xm[j] / W.rmu[j] + (1 - 0.5xmu) * P.g
        Δr = 0.5 * (W.r[k+1] - W.r[k])
        # (moving interfaces on lgrid = 2, GRIDVEL or GRIDEQ: the lab-frame mass flux r² u ρ_b; zones carried by the grid also change Σ m Φ)
        uf(j) = gridmoves() ? W.u[j] : W.urel[j]
        k > c.ngrs && (w += W.rmu[k] * uf(k) * W.db[k] * gi_(k) * Δr)
        k < c.ngre && (w += W.rmu[k+1] * uf(k + 1) * W.db[k+1] * gi_(k + 1) * Δr)
        if ETOT_GRAV_GRID[] && !gridmoves()
            # Mass carried by moving faces (Lagrangian zones: urel = 0, so the term above vanishes): the momentum row's force g_j (m_{j-1} + m_j)/2 does the work
            # u_grid,j g_j (m_{j-1} + m_j)/2 at face j, half of it in each neighbouring zone, matching the zone's kinetic energy m(u_k² + u_{k+1}²)/4. Zero on a fixed mesh.
            mh = 0.5 * W.d[k] * W.dvol[k]
            k > c.ngrs && (w += mh * (W.u[k] - W.urel[k]) * gi_(k))
            w += mh * (W.u[k+1] - W.urel[k+1]) * gi_(k + 1)
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
    residual!(R, ws, c, zeros(neqn * nb))
    it = c.ix[5]; W = ws.a
    imb = sum(R[it, k - c.k0 + 1] for k in c.ngrs:c.ngre) * dt
    E(k) = W.dn[k] * W.dvoln[k] * ((c.lrad ? W.egrn[k] : W.egn[k]) + 0.25 * (W.un[k]^2 + W.un[k+1]^2))
    Eo(k) = c.do_[k] * c.dvolo[k] * ((c.lrad ? c.egro[k] : c.ego[k]) + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))
    dE = sum(E(k) - Eo(k) for k in c.ngrs:c.ngre)
    (; dE, Enew = sum(E(k) for k in c.ngrs:c.ngre), fin = eflux(c.ngrs, c, W, Float64) * dt,
       fout = eflux(c.ngre + 1, c, W, Float64) * dt, imb)
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
    residual!(R, ws, c, zeros(neqn * nb)); W = ws.a
    it = c.ix[5]
    imb = sum(R[it, k - c.k0 + 1] for k in c.ngrs:c.ngre) * dt
    E(k) = W.dn[k] * W.dvoln[k] * ((c.lrad ? W.egrn[k] : W.egn[k]) + 0.25 * (W.un[k]^2 + W.un[k+1]^2))
    Eo(k) = c.do_[k] * c.dvolo[k] * ((c.lrad ? c.egro[k] : c.ego[k]) + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))
    dE = sum(E(k) - Eo(k) for k in c.ngrs:c.ngre)
    (; srcw = copy(W.srcw), fcv = copy(W.fcv), vcv = copy(W.vcv), lcv = copy(W.lcv), nab = copy(W.nabv), nad = copy(W.nadv), dE,
       fin = eflux(c.ngrs, c, W, Float64) * dt, fout = eflux(c.ngre + 1, c, W, Float64) * dt, imb)
end

"""
    zone_diag(dt) -> (; imb, useint, tiew, rmom, u, ...)

(Diagnostic, EDIAG=2.) Right after a converged stage: per zone the residual of the conservative total-energy row (all zones, ηsw = 0) at the
converged state [erg s⁻¹ sr⁻¹] (its sum × dt is `stage_diag`'s imb), which zones the solve actually closed with the internal-energy row
(`useint` of the run's own context), the FLUXTIE weight per interface, and the residual of the radiation-momentum row in its moment-equation
form (FLUXTIE off) at every interface [dyn sr⁻¹]: where the solve used the tie, total row − kinetic − radiation row leaves −u × that residual
on the gas (ETASW_TIE), so Σ u_k R_mom,k is the leak that bookkeeping predicts.
"""
function zone_diag(dt)
    c0 = context(; dt); ui = copy(c0.useint)
    c = context(; etotal = true, ηsw = 0.0, dt)
    nb = c.nb; neqn = c.neqn; R = zeros(neqn, nb); ws = getws(Float64, c.mgr)
    residual!(R, ws, c, zeros(neqn * nb)); W = ws.a
    it = c.ix[5]; iF = c.ix[7]; a, b = c.ngrs, c.ngre
    imb = [R[it, k - c.k0 + 1] for k in a:b]
    # the part of R_tot − R_gas − R_E that is the advection-limiter mismatch: the total row advects e_g + E/ρ with one van Leer limiter
    # (egb), the gas row e_g (egbg) and the radiation row E/ρ (erb) with their own; face term −ṁ (egb − egbg − erb)
    dm = W.dmdt; fa(j) = val(dm[j] * (W.egb[j] - W.egbg[j] - W.erb[j]))
    adv = [(k > a ? fa(k) : 0.0) - (k < b ? fa(k + 1) : 0.0) for k in a:b]
    kin = [val(0.25 * (W.dn[k] * W.dvoln[k] * (W.un[k]^2 + W.un[k+1]^2) - c.do_[k] * c.dvolo[k] * (c.uo[k]^2 + c.uo[k+1]^2)) / dt) for k in a:b]
    tw = length(TD.TIEW) >= b + 1 ? copy(TD.TIEW[a:b+1]) : zeros(b - a + 2)
    ft = TD.FLUXTIE[]; TD.FLUXTIE[] = (ft[1], 0.0)
    rm = try
        R2 = zeros(neqn, nb); residual!(R2, ws, c, zeros(neqn * nb)); [R2[iF, k - c.k0 + 1] for k in a:b+1]
    finally
        TD.FLUXTIE[] = ft
    end
    W2 = ws.a
    (; imb, adv, kin, useint = ui[a:b], tiew = tw, rmom = rm, u = copy(W2.u[a:b+1]), r = copy(W2.rn[a:b+1]), E = copy(W2.ern[a:b]), T = copy(W2.tn[a:b]),
       F = copy(W2.frn[a:b+1]), ρ = copy(W2.dn[a:b]))
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

matgen!(; ncol = 6) = (assemble!(context(); ncol); nothing)

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
