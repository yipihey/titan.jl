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
AV_ZONE[] (G1): viscous length proportional to the local zone width. When > 0, zone k uses the length
l_k = AV_ZONE[] × (r[k+1] − r[k]) (θ-centred radii, ForwardDiff-differentiable) instead of ql0 + ql1 r̄; ql0, ql1 are then
ignored (cq1, cq2, cqvis still apply). On grids where Δr/r varies (geometric mass grids) l ∝ r does not follow the zone
width; l ∝ Δr keeps the shock spread over the same number of zones. Mapping to the von Neumann–Richtmyer form
q ≈ c_q ρ (Δv)² (Ripamonti et al. 2002 use c_q = 4): with qv = div u and l = c Δr the quadratic term is
q = cqvis·cq2·l²·ρ·qv·(du/dr) ≈ cqvis·cq2·c²·ρ (Δr div u)², i.e. c_q ≈ cqvis·cq2·c² (cqvis = 4/3 in our states, so
c_q = 4 needs c² ≈ 3/cq2). This mapping is derived here, not taken from the source. Default 0 = the old formula.
"""
const AV_ZONE = Ref(0.0)
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
        m = dO[k] * dv[k]; ein = m * eo[k]; K = 0.25m * (uo[k]^2 + uo[k+1]^2)
        ui[k] = ein / (ein + K) < ηsw
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
    τ = 0.0
    @inbounds for k in ngre:-1:ngrs+1
        (base[k, 3] > 0 && base[k, 5] > 0) &&
            (τ += material_opac(m, base[k, 3], base[k, 5])[1] * base[k, 3] * (base[k+1, 1] - base[k, 1]))
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
    :erb, :ersn, :frb, :frsn, :qe, :qeo, :qF, :qFo, :dfm, :dfe, :xnu, :xnc, :ss, :fcv, :egbg, :egsng, :xnen, :snuc, :srcw, :vcv, :lcv, :nabv, :nadv)
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
const ZONE = Ref(0)
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
    @inbounds for k in ngrs:ngre
        (dO[k] > 0 && tO[k] > 0) || continue
        ZONE[] = k; buf[k] = material_source(m, dO[k], tO[k])
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
function _material!(m, W, c, P, dn, tn, ern, k0, ngre)
    pgn, egn, asn, chifn, xken, xkpn, plfn, egrn = W.pgn, W.egn, W.asn, W.chifn, W.xken, W.xkpn, W.plfn, W.egrn
    hassrc = has_source(m)
    @inbounds for k in k0:ngre+2
        if !(ForwardDiff.value(dn[k]) > 0 && ForwardDiff.value(tn[k]) > 0)   # unused ghost zones (e.g. ngrs-2)
            pgn[k] = zero(dn[k]); egn[k] = zero(dn[k]); asn[k] = zero(dn[k]); hassrc && (W.snuc[k] = zero(dn[k]))
            c.lrad && (chifn[k] = zero(dn[k]); xken[k] = zero(dn[k]); xkpn[k] = zero(dn[k]); plfn[k] = zero(dn[k]); egrn[k] = zero(dn[k]))
            continue
        end
        ZONE[] = k
        if COMPTON[]
            pg, eg, as_, ne = material_gasne(m, dn[k], tn[k]); W.xnen[k] = ne
        else
            pg, eg, as_ = material_gas(m, dn[k], tn[k])
        end
        pgn[k] = pg; egn[k] = eg; asn[k] = as_
        hassrc && (W.snuc[k] = material_source(m, dn[k], tn[k]))
        if c.lrad
            χ, κE, κP = material_opac(m, dn[k], tn[k])
            chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
            t2 = tn[k] * tn[k]
            plfn[k] = P.csigr * (t2 * t2) / P.cpi
            egrn[k] = egn[k] + ern[k] / dn[k]
        end
    end
    nothing
end

# ------------------------------------------------------------------ residual
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
        cde = c.code[j]; a = newv[cde]
        if cde == 4
            for bb in 1:nb; k = k0 + bb - 1; a[k] = c.base[k, 4] + x[(bb-1)*neqn+j] * c.nomu[k]; end
        elseif cde == 7
            for bb in 1:nb; k = k0 + bb - 1; a[k] = c.base[k, 7] + x[(bb-1)*neqn+j] * c.nomf[k]; end
        else
            for bb in 1:nb; k = k0 + bb - 1; a[k] = c.base[k, cde] * exp(x[(bb-1)*neqn+j]); end
        end
    end
    # EOS / opacity (TITAN built-in physics unless a material is installed, see TitanMaterials)
    Rg = P.ck / (P.cm0 * P.gmu); γ = P.gam
    pgn, egn, asn, chifn, xken, xkpn, plfn, egrn = W.pgn, W.egn, W.asn, W.chifn, W.xken, W.xkpn, W.plfn, W.egrn
    MATERIAL[] === nothing || _material!(MATERIAL[], W, c, P, dn, tn, ern, k0, ngre)
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
    # viscosity
    qf, dudr = W.qf, W.dudr; xmu4 = 0.25xmu
    @inbounds for k in ngrs:ngre
        dudr[k] = (u[k+1] - u[k]) / (r[k+1] - r[k]) - xmu4 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
        ql = AV_ZONE[] > 0 ? AV_ZONE[] * (r[k+1] - r[k]) : P.ql0 + P.ql1 * 0.5 * (r[k] + r[k+1])
        div = (rmu[k+1] * u[k+1] - rmu[k] * u[k]) / dvol[k]
        qv = decide(:vis, 0, k, val(div) < 0) ? div : zero(T)
        qm = ql * (P.cq1 * as[k] - P.cq2 * ql * qv)
        qf[k] = -P.cqvis * d[k] * qm * dudr[k]
    end
    # advection
    db, dsn, fl, ub, usn, flu, fle, egb, egsn = W.db, W.dsn, W.fl, W.ub, W.usn, W.flu, W.fle, W.egb, W.egsn
    @inbounds for k in ngrs+1:ngre+1; fl[k] = urel[k]; fle[k] = -dmdt[k]; end
    @inbounds for k in ngrs:ngre; flu[k] = -(dmdt[k] + dmdt[k+1]); end
    advect!(db, dsn, dn, c.do_, c.dso, fl, θ, P.cadv, P.eps, ngrs, ngre, true, 1)
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
        @inbounds for k in ngrs:ngre+1; R[ir, kb(k)] = rn[k] - c.ro[k]; end
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
        c.lrad && (w -= avchi[k] * fr[k] / cc * mk)
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
                th = material_thermo(MATERIAL[], ρi, Ti); nad, cpk, δk = th
                κk = material_opac(MATERIAL[], ρi, Ti)[1]
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
        elseif c.lrad && GASROW[]
            egbg = W.egbg
            w = (dn[k] * egn[k] * dvoln[k] - c.do_[k] * c.ego[k] * c.dvolo[k]) / dt -
                dfe[k+1] + dfe[k] + rmu[k+1] * fcv[k+1] - rmu[k] * fcv[k] +
                pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * dudr[k] * dvol[k] -
                (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) * d[k] * dvol[k] + compton(k, W, P) * dvol[k] - src
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
        R[it, kb(k)] = w
    end
    R[it, 1] = X(it, k0)
    k = ngrs - 1
    R[it, kb(k)] = (P.leibc == 1 || P.llibc == 1) ? egn[k+1] - egn[k] : X(it, k)
    for k in ngre+1:ngre+2
        R[it, kb(k)] = (P.leobc in (1, 3) || P.llobc > 0) ? -egn[k-1] + egn[k] : X(it, k)
    end
    if c.lrad
        # ------------------------------------------------------------ RADIATION ENERGY
        @inbounds for k in ngrs:ngre
            aur = 0.5 * (u[k+1] / r[k+1] + uor(u, r, k, ngrs))
            w = (ern[k] * dvoln[k] - c.ero[k] * c.dvolo[k]) / dt + rmu[k+1] * fr[k+1] - rmu[k] * fr[k] +
                fedd[k] * er[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur * dvol[k] +
                (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) * d[k] * dvol[k] - compton(k, W, P) * dvol[k]
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
        @inbounds for k in ngrs+1:ngre
            qFo = c.fro[k] * 2 / (c.do_[k] + c.do_[k-1])
            R[iF, kb(k)] =
                (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                 qFo * (c.do_[k-1] * c.dvolo[k-1] + c.do_[k] * c.dvolo[k])) / (2cc^2 * dt) -
                ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2) +
                rmu[k] * (fedd[k] * er[k] - fedd[k-1] * er[k-1]) +
                0.5 * rmu[k] * (fr[k] / cc^2) * (u[k+1] - u[k-1]) +
                (0.25xmu / r[k]) * ((3fedd[k-1] - 1) * er[k-1] * dvol[k-1] + (3fedd[k] - 1) * er[k] * dvol[k]) +
                0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] * fr[k] / cc
        end
        P.lribc in (2, 3) || error("lribc must be 2 or 3")
        R[iF, kb(ngrs)] = X(iF, ngrs)
        R[iF, kb(ngre+1)] = P.cpi * (2P.geddr + 1) * P.ximr + frn[ngre+1] - cc * P.geddr * ern[ngre]
        R[iF, 1] = X(iF, k0)
        k = ngrs - 1
        R[iF, kb(k)] = P.lribc == 2 ? frn[k+2] + frn[k] : rn[k+1]^mu * frn[k+1] - rn[k]^mu * frn[k]
        R[iF, kb(ngre+2)] = -frn[ngre+1] + frn[ngre+2]
    end
    nothing
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
"Conservative total-energy row of zone k: Δ[m(e_g [+E_r/ρ]) + m(u_k²+u_{k+1}²)/4]/dt + Φ_{k+1} − Φ_k."
@inline function etot_row(k, c, W, ::Type{T}) where {T}
    dn, un, dvoln = W.dn, W.un, W.dvoln
    en = c.lrad ? W.egrn[k] : W.egn[k]; eo = c.lrad ? c.egro[k] : c.ego[k]
    (dn[k] * dvoln[k] * (en + 0.25 * (un[k]^2 + un[k+1]^2)) -
     c.do_[k] * c.dvolo[k] * (eo + 0.25 * (c.uo[k]^2 + c.uo[k+1]^2))) / c.p.dt +
    eflux(k + 1, c, W, T) - eflux(k, c, W, T)
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

end # module
