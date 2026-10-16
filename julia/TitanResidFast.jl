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
end

const CTXBUF = Ref{Any}(nothing)
const BUFS = Dict{Symbol,Any}()     # reused context buffers (allocation-free context())
"Defaults for the energy formulation (Phase 4 item 2); ETOTAL[] = true selects the total-energy row."
const ETOTAL = Ref(false)
const ETASW = Ref(0.0)
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
    (gi(:ltur) == 0 && gi(:lhydr) == 1) || error("TitanResidFast: lhydr=1, ltur=0 only")
    gi(:leos) == 2 || error("leos must be 2")
    lrad && !(gi(:lopac) in (3, 4)) && error("lopac must be 3 or 4")
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
        (la = v(:lady); lb = get!(() -> zeros(Int, length(la)), BUFS, :lady)::Vector{Int}; lb .= la; lb), v(:wt), v(:yscl), P, etotal, dualswitch(ngrs, ngre, lrad, ηsw))
end

# ------------------------------------------------------------------ workspace
const WSNAMES = (:rn, :xmen, :dn, :un, :tn, :ern, :frn, :pgn, :egn, :asn, :chifn, :xken, :xkpn, :plfn,
    :egrn, :r, :u, :fr, :xm, :dmdt, :urel, :rmu, :rmup1, :d, :pg, :eg, :as, :er, :egr, :chif, :xke,
    :xkp, :plf, :dvol, :dvoln, :avchi, :qf, :dudr, :db, :dsn, :fl, :ub, :usn, :flu, :fle, :egb, :egsn,
    :erb, :ersn, :frb, :frsn, :qe, :qeo, :qF, :qFo, :dfm, :dfe, :xnu, :xnc, :ss)
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
function material_gas end
function material_opac end
"Function barrier: fills the n-level EOS/opacity workspace from the material (type-specialised loop)."
function _material!(m, W, c, P, dn, tn, ern, k0, ngre)
    pgn, egn, asn, chifn, xken, xkpn, plfn, egrn = W.pgn, W.egn, W.asn, W.chifn, W.xken, W.xkpn, W.plfn, W.egrn
    @inbounds for k in k0:ngre+2
        if !(ForwardDiff.value(dn[k]) > 0 && ForwardDiff.value(tn[k]) > 0)   # unused ghost zones (e.g. ngrs-2)
            pgn[k] = zero(dn[k]); egn[k] = zero(dn[k]); asn[k] = zero(dn[k])
            c.lrad && (chifn[k] = zero(dn[k]); xken[k] = zero(dn[k]); xkpn[k] = zero(dn[k]); plfn[k] = zero(dn[k]); egrn[k] = zero(dn[k]))
            continue
        end
        pg, eg, as_ = material_gas(m, dn[k], tn[k])
        pgn[k] = pg; egn[k] = eg; asn[k] = as_
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
        xm[k] = θ * (P.xmtot - xmen[k]) + (1 - θ) * c.xmo[k]
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
        dudr[k] = (u[k+1] - u[k]) / (r[k+1] - r[k]) - xmu4 * (u[k+1] / r[k+1] + u[k] / r[k])
        ql = P.ql0 + P.ql1 * 0.5 * (r[k] + r[k+1])
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
    R[im, kb(ngrs)] = xmen[ngrs] - c.xmeo[ngrs] + rmu[ngrs] * P.phil0 * dt
    R[im, 1] = X(im, k0)
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
    (P.leibc in (1, 2) || P.llibc == 1) || error("gasmom inner BC")
    (P.leobc < 3 || P.llobc == 1) || error("gasmom outer BC")
    k = ngrs - 1
    R[iu, kb(k)] = P.leibc == 1 ? un[k+2] + un[k] : P.llibc == 1 ? un[k] - un[k+1] : X(iu, k)
    k = ngre + 2
    R[iu, kb(k)] = P.leobc == 1 ? un[k-2] + un[k] : P.leobc == 3 ? -un[k-1] + un[k] : X(iu, k)
    # ---------------------------------------------------------------- GAS ENERGY
    @inbounds for k in ngrs:ngre
        aur = 0.5 * (u[k+1] / r[k+1] + u[k] / r[k])
        if c.etotal && !c.useint[k]
            w = etot_row(k, c, W, T)
        elseif c.lrad
            w = (dn[k] * egrn[k] * dvoln[k] - c.do_[k] * c.egro[k] * c.dvolo[k]) / dt +
                rmu[k+1] * fr[k+1] - rmu[k] * fr[k] - dfe[k+1] + dfe[k] +
                (pg[k] + er[k] * fedd[k]) * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur * dvol[k] + qf[k] * dudr[k] * dvol[k]
            if k == ngrs
                w += -dmdt[k+1] * egb[k+1] - (P.phil2 + er[k] * urel[k]) * rmu[k]
            elseif k == ngre
                w += dmdt[k] * egb[k]
                w += P.leobc != 3 ? rmu[k+1] * (P.phir2 + urel[k+1] * er[k]) : -dmdt[k+1] * egb[k+1]
            else
                w += -dmdt[k+1] * egb[k+1] + dmdt[k] * egb[k]
            end
        else
            w = (dn[k] * egn[k] * dvoln[k] - c.do_[k] * c.ego[k] * c.dvolo[k]) / dt -
                dfe[k+1] + dfe[k] + pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * dudr[k] * dvol[k]
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
            aur = 0.5 * (u[k+1] / r[k+1] + u[k] / r[k])
            w = (ern[k] * dvoln[k] - c.ero[k] * c.dvolo[k]) / dt + rmu[k+1] * fr[k+1] - rmu[k] * fr[k] +
                fedd[k] * er[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur * dvol[k] +
                (cc * xke[k] * er[k] - 4P.cpi * xkp[k] * plf[k]) * d[k] * dvol[k]
            if k == ngrs
                w += -dmdt[k+1] * erb[k+1] - urel[k] * rmu[k] * er[k]
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
