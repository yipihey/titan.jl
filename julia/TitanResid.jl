"""
    TitanResid

Julia port of TITAN's complete Newton residual (matgen: grid, mass, contin, viscous,
gasmom{h,rh}, gasnrg{h,rh}, radnrg rh, radmom rh, with advectc/advecti/diffuse) as a
single function of the Newton variables, so the whole block-pentadiagonal Jacobian
comes from one stencil-coloured ForwardDiff sweep.

Newton variables (TITAN's): for block k and equation/variable j
  log-type (r, m_ext, ρ, T, E):  value = base · exp(x)      (x = δ ln, about the current state)
  u, F:                          value = base + x · nom      (nom = unom / frnom)
x = 0 is the current state.  Rows that TITAN "prescribes" (rhs = 0, e00 = 1) are R = x.

Supported (everything the SOD/SUBC/SUPC benchmarks use, plus a few neighbours):
  lgrid 1 (ladx 1|2; lady 1|2|3 for xm, ρ, T, E, p, e, χ, F, bri; 11), 2, 3
  leibc 0|1, leobc 0|1, llibc 0|1, llobc 0, lribc 2|3, lrobc 1,
  leos 2, lopac 3|4, ltran 1 (Eddington factors frozen), ltur 0, lhydr 1.
Anything else -> error (the Fortran matgen remains available).

Note: `etotal = true` replaces the gas(+radiation) internal-energy equation by a
conservative TOTAL energy equation (see total_energy_rows below; Phase 4 item 2).
"""
module TitanResid

using ForwardDiff
import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

gi(s) = Int(TD.g(s)); gf(s) = Float64(TD.g(s))
# Branch on VALUES only, as TITAN does.  ForwardDiff ≥ 1.0 compares Duals including their
# partials when the values tie (e.g. a*b == 0 at a limiter switch), which would silently
# take the other branch than the Fortran at exact ties (found at SUBC start-up, Phase 4).
val(x) = ForwardDiff.value(x)
cp_(s) = copy(TF.var(s))

"Everything constant during one Newton iteration (read from COMMON)."
struct Ctx
    ngrs::Int; ngre::Int; neqn::Int; mu::Int; lrad::Bool; mgr::Int
    ix::NTuple{7,Int}                     # ir im id iu it ie if  (0 if absent)
    var::Vector{Symbol}                   # var[j] = :rn, :xmen, ...
    islog::Vector{Bool}
    base::Dict{Symbol,Vector{Float64}}    # current new-level arrays
    nom::Dict{Symbol,Vector{Float64}}
    o::Dict{Symbol,Vector{Float64}}       # old-level / frozen arrays
    p::Dict{Symbol,Float64}               # scalar parameters
    ip::Dict{Symbol,Int}                  # integer switches
    lady::Vector{Int}; wt::Vector{Float64}; yscl::Vector{Float64}
    etotal::Bool
end

const OLD = (:ro, :xmeo, :do, :uo, :to, :ego, :ero, :fro, :egro, :egso, :egrso, :dso, :uso,
             :erso, :frso, :pgo, :aso, :chifo, :xkeo, :xkpo, :plfo, :dvolo, :xmo, :fedd, :bri,
             :egrto)
const PSC = (:dtime, :thet, :cadv, :epsadv, :sigd, :sige, :phil0, :phil2, :phir0, :phir2, :cgrav, :g,
             :xmtot, :alph, :tau, :xscale, :ql0, :ql1, :cq1, :cq2, :cqvis, :gam, :gmu, :ck, :cm0, :csige,
             :chif0, :ratio, :cc, :cpi, :car, :csigr, :geddl, :geddr, :ximr, :xipl, :delml, :delmr,
             :delrl, :delrr, :dl, :dr)
const ISW = (:ibet, :ladx, :lgrid, :leibc, :leobc, :llibc, :llobc, :lribc, :lrobc, :leos, :lopac, :ltran)

function context(; etotal = false)
    ngrs = gi(:ngrs); ngre = gi(:ngre); neqn = gi(:neqn); lrad = gi(:lrad) == 1
    ix = Tuple(gi(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    gi(:ltur) == 0 && gi(:lhydr) == 1 || error("TitanResid: only lhydr=1, ltur=0")
    var = Vector{Symbol}(undef, neqn); islog = fill(true, neqn)
    for (s, i, lg) in ((:rn, ix[1], true), (:xmen, ix[2], true), (:dn, ix[3], true), (:un, ix[4], false),
                       (:tn, ix[5], true), (:ern, ix[6], true), (:frn, ix[7], false))
        (1 <= i <= neqn) || continue
        (s in (:ern, :frn) && !lrad) && continue
        var[i] = s; islog[i] = lg
    end
    base = Dict(s => cp_(s) for s in var)
    lrad || (base[:ern] = cp_(:ern); base[:frn] = cp_(:frn))
    nom = Dict(:un => cp_(:unom), :frn => cp_(:frnom))
    o = Dict(s => cp_(s) for s in OLD)
    p = Dict(s => gf(s) for s in PSC)
    ip = Dict(s => gi(s) for s in ISW)
    ip[:leos] == 2 || error("TitanResid: leos must be 2 (ideal gas)")
    lrad && !(ip[:lopac] in (3, 4)) && error("TitanResid: lopac must be 3 or 4")
    lrad && ip[:ltran] != 1 && error("TitanResid: ltran must be 1")
    Ctx(ngrs, ngre, neqn, gi(:mu), lrad, TF.MGR, ix, var, islog, base, nom, o, p, ip,
        Int.(TF.var(:lady)), copy(TF.var(:wt)), copy(TF.var(:yscl)), etotal)
end

# ------------------------------------------------------------------ helpers
"van Leer interface/centred interpolation, shared by advectc (centred=true) and advecti."
function advect!(qb, qsn, qn, qo, qso, flow, θ, cadv, eps, ngrs, ngre, centred::Bool)
    T = eltype(qb)
    dq(k) = qn[k+1] - qn[k]
    @inbounds for k in ngrs:ngre+1
        a = dq(k - 1); b = dq(k); sdq = a + b
        den = sdq + (val(sdq) >= 0 ? eps : -eps)
        qsn[k] = val(a * b) >= 0 ? cadv * a * b / den : zero(T)
    end
    if centred      # advectc: qb(k) at interface k from zones k-1 | k, k = ngrs+1..ngre+1
        @inbounds for k in ngrs+1:ngre+1
            s = val(flow[k]) >= 0 ? 0.5 : -0.5
            qb[k] = (0.5 + s) * (θ * (qn[k-1] + 0.5qsn[k-1]) + (1 - θ) * (qo[k-1] + 0.5qso[k-1])) +
                    (0.5 - s) * (θ * (qn[k]   - 0.5qsn[k])   + (1 - θ) * (qo[k]   - 0.5qso[k]))
        end
    else            # advecti: qb(k) at zone k from interfaces k | k+1, k = ngrs..ngre
        @inbounds for k in ngrs:ngre
            s = val(flow[k]) >= 0 ? 0.5 : -0.5
            qb[k] = (0.5 + s) * (θ * (qn[k]   + 0.5qsn[k])   + (1 - θ) * (qo[k]   + 0.5qso[k])) +
                    (0.5 - s) * (θ * (qn[k+1] - 0.5qsn[k+1]) + (1 - θ) * (qo[k+1] - 0.5qso[k+1]))
        end
    end
    nothing
end

"cs(k) for one resolution variable (grid.f), new-level values y."
function rescs(y, k, lady, yscl)
    lady == 1 && return (y[k+1] - y[k]) / yscl
    lady == 2 && return (y[k+1] - y[k]) / (y[k+1] + y[k])
    lady == 3 && return y[k+1] / y[k] - y[k] / y[k+1]
    error("lady=$lady")
end

# ------------------------------------------------------------------ residual
"""
    residual!(R, c, x)

R[i, kb] for equation i and block kb (k = ngrs-2+kb-1), x = Newton-variable vector
(layout (kb-1)*neqn + j).  Returns also the Float64 side quantities when T == Float64.
"""
function residual!(R::AbstractMatrix{T}, c::Ctx, x::AbstractVector{T}) where {T}
    ngrs, ngre, neqn, mu, mgr = c.ngrs, c.ngre, c.neqn, c.mu, c.mgr
    k0 = ngrs - 2; nb = ngre - ngrs + 5
    o, p, ip = c.o, c.p, c.ip
    ir, im, id, iu, it, ie, iF = c.ix
    θ = p[:thet]; dt = p[:dtime]; xmu = Float64(mu); xmup1 = xmu + 1
    # ---- new-level arrays from the Newton variables
    A = Dict{Symbol,Vector{T}}()
    for s in (:rn, :xmen, :dn, :un, :tn, :ern, :frn)
        A[s] = haskey(c.base, s) ? Vector{T}(c.base[s]) : zeros(T, mgr)
    end
    X(j, kb) = x[(kb - 1) * neqn + j]
    for j in 1:neqn
        s = c.var[j]; a = A[s]; b = c.base[s]
        if c.islog[j]
            @inbounds for kb in 1:nb; k = k0 + kb - 1; a[k] = b[k] * exp(X(j, kb)); end
        else
            nm = c.nom[s]
            @inbounds for kb in 1:nb; k = k0 + kb - 1; a[k] = b[k] + X(j, kb) * nm[k]; end
        end
    end
    rn, xmen, dn, un, tn, ern, frn = A[:rn], A[:xmen], A[:dn], A[:un], A[:tn], A[:ern], A[:frn]
    # ---- EOS (ideal gas) and opacity, new level, k = ngrs-2 .. ngre+2
    Rg = p[:ck] / (p[:cm0] * p[:gmu]); γ = p[:gam]
    pgn = zeros(T, mgr); egn = zeros(T, mgr); asn = zeros(T, mgr)
    chifn = zeros(T, mgr); xken = zeros(T, mgr); xkpn = zeros(T, mgr); plfn = zeros(T, mgr)
    egrn = zeros(T, mgr)
    @inbounds for k in k0:ngre+2
        pgn[k] = Rg * dn[k] * tn[k]; egn[k] = Rg * tn[k] / (γ - 1)
        asn[k] = sqrt(γ * pgn[k] / dn[k])
        if c.lrad
            if ip[:lopac] == 3
                chifn[k] = p[:chif0] / dn[k]
            else
                chifn[k] = T(p[:csige] / p[:cm0])
            end
            xken[k] = p[:ratio] * chifn[k]; xkpn[k] = p[:ratio] * chifn[k]
            plfn[k] = p[:csigr] * tn[k]^4 / p[:cpi]
            egrn[k] = egn[k] + ern[k] / dn[k]
        end
    end
    # ---- time-centred quantities (update*/peq*)
    r = zeros(T, mgr); u = zeros(T, mgr); fr = zeros(T, mgr); xm = zeros(T, mgr)
    dmdt = zeros(T, mgr); drdt = zeros(T, mgr); urel = zeros(T, mgr); rmu = zeros(T, mgr)
    rmup1 = zeros(T, mgr)
    @inbounds for k in ngrs-1:ngre+2
        r[k] = θ * rn[k] + (1 - θ) * o[:ro][k]
        u[k] = θ * un[k] + (1 - θ) * o[:uo][k]
        fr[k] = θ * frn[k] + (1 - θ) * o[:fro][k]
        xm[k] = θ * (p[:xmtot] - xmen[k]) + (1 - θ) * o[:xmo][k]
        dmdt[k] = -(xmen[k] - o[:xmeo][k]) / dt
        drdt[k] = (rn[k] - o[:ro][k]) / dt
        urel[k] = u[k] - drdt[k]
        rmu[k] = r[k]^mu; rmup1[k] = r[k]^(mu + 1)
    end
    d = zeros(T, mgr); tt = zeros(T, mgr); pg = zeros(T, mgr); eg = zeros(T, mgr); as = zeros(T, mgr)
    er = zeros(T, mgr); egr = zeros(T, mgr); chif = zeros(T, mgr); xke = zeros(T, mgr); xkp = zeros(T, mgr)
    plf = zeros(T, mgr); dvol = zeros(T, mgr); dvoln = zeros(T, mgr)
    @inbounds for k in ngrs-1:ngre+1
        d[k] = θ * dn[k] + (1 - θ) * o[:do][k]
        pg[k] = θ * pgn[k] + (1 - θ) * o[:pgo][k]
        eg[k] = θ * egn[k] + (1 - θ) * o[:ego][k]
        as[k] = θ * asn[k] + (1 - θ) * o[:aso][k]
        if c.lrad
            er[k] = θ * ern[k] + (1 - θ) * o[:ero][k]
            egr[k] = θ * egrn[k] + (1 - θ) * o[:egro][k]
            chif[k] = θ * chifn[k] + (1 - θ) * o[:chifo][k]
            xke[k] = θ * xken[k] + (1 - θ) * o[:xkeo][k]
            xkp[k] = θ * xkpn[k] + (1 - θ) * o[:xkpo][k]
            plf[k] = θ * plfn[k] + (1 - θ) * o[:plfo][k]
        end
        dvol[k] = (rmup1[k+1] - rmup1[k]) / xmup1
        dvoln[k] = (rn[k+1]^(mu + 1) - rn[k]^(mu + 1)) / xmup1
    end
    avchi = zeros(T, mgr)
    c.lrad && @inbounds for k in ngrs:ngre+1; avchi[k] = 2 / (1 / chif[k-1] + 1 / chif[k]); end
    # ---- artificial viscosity (viscous.f), zones ngrs..ngre
    qf = zeros(T, mgr); dudr = zeros(T, mgr)
    xmu4 = 0.25xmu
    f1(xk, xp, uk, up) = (up - uk) / (xp - xk) - xmu4 * (up / xp + uk / xk)
    @inbounds for k in ngrs:ngre
        dudr[k] = f1(r[k], r[k+1], u[k], u[k+1])
        ql = p[:ql0] + p[:ql1] * 0.5 * (r[k] + r[k+1])
        div = (rmu[k+1] * u[k+1] - rmu[k] * u[k]) / dvol[k]
        qv = val(div) < 0 ? div : zero(T)
        qm = ql * (p[:cq1] * as[k] - p[:cq2] * ql * qv)
        qf[k] = -p[:cqvis] * d[k] * qm * dudr[k]
    end
    # ---- advection (interface/centred interpolations)
    cadv, eps = p[:cadv], p[:epsadv]
    db = zeros(T, mgr); dsn = zeros(T, mgr); fl = zeros(T, mgr)
    @inbounds for k in ngrs+1:ngre+1; fl[k] = urel[k]; end
    advect!(db, dsn, dn, o[:do], o[:dso], fl, θ, cadv, eps, ngrs, ngre, true)
    ub = zeros(T, mgr); usn = zeros(T, mgr); flu = zeros(T, mgr)
    @inbounds for k in ngrs:ngre; flu[k] = -(dmdt[k] + dmdt[k+1]); end
    advect!(ub, usn, un, o[:uo], o[:uso], flu, θ, cadv, eps, ngrs, ngre, false)
    fle = zeros(T, mgr); @inbounds for k in ngrs+1:ngre+1; fle[k] = -dmdt[k]; end
    egb = zeros(T, mgr); egsn = zeros(T, mgr)          # (gas or gas+rad) specific energy
    if c.lrad
        advect!(egb, egsn, egrn, o[:egro], o[:egrso], fle, θ, cadv, eps, ngrs, ngre, true)
    else
        advect!(egb, egsn, egn, o[:ego], o[:egso], fle, θ, cadv, eps, ngrs, ngre, true)
    end
    erb = zeros(T, mgr); ersn = zeros(T, mgr); frb = zeros(T, mgr); frsn = zeros(T, mgr)
    if c.lrad
        qe = zeros(T, mgr); qeo = zeros(T, mgr)
        @inbounds for k in ngrs-1:ngre+2; qe[k] = ern[k] / dn[k]; qeo[k] = o[:ero][k] / o[:do][k]; end
        advect!(erb, ersn, qe, qeo, o[:erso], fle, θ, cadv, eps, ngrs, ngre, true)
        qf_ = zeros(T, mgr); qfo = zeros(T, mgr)
        @inbounds for k in ngrs-1:ngre+2
            qf_[k] = frn[k] * 2 / (dn[k] + dn[k-1]); qfo[k] = o[:fro][k] * 2 / (o[:do][k] + o[:do][k-1])
        end
        advect!(frb, frsn, qf_, qfo, o[:frso], flu, θ, cadv, eps, ngrs, ngre, false)
    end
    # ---- diffusion fluxes
    dfm = zeros(T, mgr); dfe = zeros(T, mgr)
    if p[:sigd] != 0
        @inbounds for k in ngrs:ngre+1
            dfm[k] = 2p[:sigd] * rmu[k] * (d[k] - d[k-1]) / (r[k+1] - r[k-1])
        end
    end
    if p[:sige] != 0
        @inbounds for k in ngrs:ngre+1
            dfe[k] = p[:sige] * rmu[k] * (d[k] + d[k-1]) * (eg[k] - eg[k-1]) / (r[k+1] - r[k-1])
        end
    end
    fedd = o[:fedd]; cc = p[:cc]
    fill!(R, zero(T))
    kb(k) = k - k0 + 1
    pres!(i, k) = (R[i, kb(k)] = X(i, kb(k)))           # TITAN's "rhs = 0, e00 = 1" rows
    # ================================================================ GRID
    lg = ip[:lgrid]
    if lg == 1
        α = p[:alph] * (1 + p[:alph]); gam = (p[:tau] / dt)^ip[:ibet]
        xnu = zeros(T, mgr); xnuo = zeros(Float64, mgr)
        ro = o[:ro]
        @inbounds for k in ngrs:ngre
            if ip[:ladx] == 1
                xnu[k] = p[:xscale] / (rn[k+1] - rn[k]); xnuo[k] = p[:xscale] / (ro[k+1] - ro[k])
            else
                xnu[k] = (rn[k+1] + rn[k]) / (rn[k+1] - rn[k]) / 2; xnuo[k] = (ro[k+1] + ro[k]) / (ro[k+1] - ro[k]) / 2
            end
        end
        xnc = zeros(T, mgr)
        @inbounds for k in ngrs+1:ngre-1
            xnt = xnu[k] - α * (xnu[k-1] - 2xnu[k] + xnu[k+1])
            xnto = xnuo[k] - α * (xnuo[k-1] - 2xnuo[k] + xnuo[k+1])
            xnc[k] = (1 + gam) * xnt - gam * xnto
        end
        ss = zeros(T, mgr)
        ys = Dict(2 => dn, 3 => tn, 4 => ern, 5 => pgn, 6 => egn, 7 => chifn, 9 => frn, 10 => o[:bri])
        for l in 1:min(length(c.lady), 11)
            (c.wt[l] > 0 && c.lady[l] > 0) || continue
            @inbounds for k in ngrs:ngre-1
                if l == 1
                    xs = c.yscl[1]
                    cs = c.lady[1] == 1 ? (xmen[k+1] - xmen[k]) / xs :
                         c.lady[1] == 2 ? (xmen[k+1] - xmen[k]) / (2xs - (xmen[k+1] + xmen[k])) :
                         (xs - xmen[k+1]) / (xs - xmen[k]) - (xs - xmen[k]) / (xs - xmen[k+1])
                elseif l == 11
                    cs = c.lady[11] == 1 ? (xmen[k+1] - xmen[k]) / c.yscl[11] : rescs(xmen, k, c.lady[11], c.yscl[11])
                elseif l == 8
                    error("TitanResid: lady(8) (viscosity resolution) not ported")
                else
                    cs = rescs(ys[l], k, c.lady[l], c.yscl[l])
                end
                ss[k] += c.wt[l] * cs^2
            end
        end
        rr(k) = sqrt(1 + xnu[k]^2 * ss[k])
        R[ir, kb(ngrs+1)] = xnu[ngrs+1] - xnu[ngrs]
        @inbounds for k in ngrs+2:ngre-1
            R[ir, kb(k)] = xnc[k] / rr(k) - xnc[k-1] / rr(k-1)
        end
        R[ir, kb(ngre)] = xnu[ngre] - xnu[ngre-1]
        k = ngrs
        ip[:leibc] > 0 && (R[ir, kb(k)] = rn[k] - ro[k])
        ip[:llibc] == 1 && (R[ir, kb(k)] = (rn[k] - ro[k]) / dt - u[k])
        k = ngre + 1
        ip[:leobc] > 0 && (R[ir, kb(k)] = rn[k] - ro[k])
        ip[:llobc] > 0 && (R[ir, kb(k)] = (rn[k] - ro[k]) / dt - u[k])
        if c.p[:dtime] > 0          # store grid concentration for writout (Float64 only)
            T === Float64 && (copyto!(TF.var(:xnu), xnu))
        end
    elseif lg == 2
        @inbounds for k in ngrs:ngre+1; R[ir, kb(k)] = rn[k] - o[:ro][k]; end
    else
        @inbounds for k in ngrs:ngre+1; R[ir, kb(k)] = (rn[k] - o[:ro][k]) / dt - u[k]; end
    end
    k = ngrs - 1                                           # general BC rows
    if ip[:leibc] == 1;  R[ir, kb(k)] = rn[k] - 2rn[k+1] + rn[k+2]
    elseif ip[:leibc] == 2; R[ir, kb(k)] = rn[k] - rn[k+1] + p[:delrl]
    elseif ip[:llibc] == 1; R[ir, kb(k)] = -p[:delml] / dn[k] + (rn[k+1]^(mu+1) - rn[k]^(mu+1)) / xmup1
    end
    k = ngre + 2
    if ip[:leobc] == 1;  R[ir, kb(k)] = rn[k] - 2rn[k-1] + rn[k-2]
    elseif ip[:leobc] in (2, 3); R[ir, kb(k)] = rn[k] - rn[k-1] - p[:delrr]
    elseif ip[:llobc] > 0; R[ir, kb(k)] = -p[:delmr] / dn[k-1] + (rn[k]^(mu+1) - rn[k-1]^(mu+1)) / xmup1
    end
    pres!(ir, ngrs - 2)
    # ================================================================ MASS
    @inbounds for k in ngrs+1:ngre+1
        R[im, kb(k)] = xmen[k-1] - xmen[k] - dn[k-1] * dvoln[k-1]
    end
    k = ngrs; R[im, kb(k)] = xmen[k] - o[:xmeo][k] + rmu[k] * p[:phil0] * dt
    pres!(im, ngrs - 2)
    k = ngrs - 1
    if ip[:leibc] == 1;  R[im, kb(k)] = xmen[k] - 2xmen[k+1] + xmen[k+2]
    elseif ip[:leibc] == 2; R[im, kb(k)] = -xmen[k] + xmen[k+1] + p[:dl] * dvoln[k]
    elseif ip[:llibc] > 0; R[im, kb(k)] = -xmen[k] + xmen[k+1] + p[:delml]
    else pres!(im, k)
    end
    k = ngre + 2
    if ip[:leobc] == 1;  R[im, kb(k)] = xmen[k-2] - 2xmen[k-1] + xmen[k]
    elseif ip[:leobc] == 2; R[im, kb(k)] = xmen[k-1] - xmen[k] - p[:dr] * dvoln[k-1]
    elseif ip[:leobc] == 3; R[im, kb(k)] = xmen[k-1] - xmen[k] - dn[k-1] * dvoln[k-1]
    elseif ip[:llobc] > 0; R[im, kb(k)] = xmen[k-1] - xmen[k] + p[:delmr]
    else pres!(im, k)
    end
    # ================================================================ CONTIN
    @inbounds for k in ngrs:ngre
        fin  = k == ngrs ? rmu[k] * p[:phil0] : rmu[k] * urel[k] * db[k]
        fout = k == ngre ? (ip[:leobc] == 3 ? rmu[k+1] * p[:phir0] + rmu[k+1] * urel[k+1] * db[k+1] :
                            rmu[k+1] * p[:phir0]) : rmu[k+1] * urel[k+1] * db[k+1]
        R[id, kb(k)] = (dn[k] * dvoln[k] - o[:do][k] * o[:dvolo][k]) / dt + fout - fin - dfm[k+1] + dfm[k]
    end
    pres!(id, ngrs - 2)
    k = ngrs - 1
    (ip[:leibc] == 1 || ip[:llibc] > 0) ? (R[id, kb(k)] = dn[k+1] - dn[k]) : pres!(id, k)
    for k in ngre+1:ngre+2
        ip[:leobc] == 1 ? (R[id, kb(k)] = -dn[k-1] + dn[k]) : pres!(id, k)
    end
    # ================================================================ GAS MOMENTUM
    r1(k) = r[k]^(-0.5xmu); r3(k) = ((r[k] + r[k+1]) / 2)^(1.5xmu)
    @inbounds for k in ngrs+1:ngre
        mk = 0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1])
        v = (un[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
             o[:uo][k] * (o[:do][k-1] * o[:dvolo][k-1] + o[:do][k] * o[:dvolo][k])) / (2dt) -
            0.5 * ((dmdt[k] + dmdt[k+1]) * ub[k] - (dmdt[k-1] + dmdt[k]) * ub[k-1]) +
            rmu[k] * (pg[k] - pg[k-1]) +
            (2xmu * p[:cpi] * p[:cgrav] * xm[k] / rmu[k] + (1 - 0.5xmu) * p[:g]) * mk +
            r1(k) * (r3(k) * qf[k] - r3(k-1) * qf[k-1])
        c.lrad && (v -= avchi[k] * fr[k] / cc * mk)
        R[iu, kb(k)] = v
    end
    (ip[:leibc] in (1, 2) || ip[:llibc] == 1) ? pres!(iu, ngrs) : error("gasmom inner BC")
    (ip[:leobc] < 3 || ip[:llobc] == 1) ? pres!(iu, ngre + 1) : error("gasmom outer BC")
    pres!(iu, ngrs - 2)
    k = ngrs - 1
    if ip[:leibc] == 1; R[iu, kb(k)] = un[k+2] + un[k]
    elseif ip[:llibc] == 1; R[iu, kb(k)] = un[k] - un[k+1]
    else pres!(iu, k)
    end
    k = ngre + 2
    if ip[:leobc] == 1; R[iu, kb(k)] = un[k-2] + un[k]
    elseif ip[:leobc] == 3; R[iu, kb(k)] = -un[k-1] + un[k]
    else pres!(iu, k)
    end
    # ================================================================ GAS ENERGY
    aur(k) = 0.5 * (u[k+1] / r[k+1] + u[k] / r[k])        # <u/r> (gasnrgrh f1, not viscous f1)
    if c.etotal
        total_energy_rows!(R, c, (; it, kb, dn, dvoln, egn, egrn, ern, un, rmu, u, fr, pg, er, qf, dmdt, egb,
                                   urel, d, eg, r, dfe, dvol, k0, θ, fedd), T)
    else
        @inbounds for k in ngrs:ngre
            if c.lrad
                v = (dn[k] * egrn[k] * dvoln[k] - o[:do][k] * o[:egro][k] * o[:dvolo][k]) / dt +
                    rmu[k+1] * fr[k+1] - rmu[k] * fr[k] - dfe[k+1] + dfe[k] +
                    (pg[k] + er[k] * fedd[k]) * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                    (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur(k) * dvol[k] + qf[k] * dudr[k] * dvol[k]
                if k == ngrs
                    v += -dmdt[k+1] * egb[k+1] - (p[:phil2] + er[k] * urel[k]) * rmu[k]
                elseif k == ngre
                    v += dmdt[k] * egb[k]
                    ip[:leobc] != 3 ? (v += rmu[k+1] * (p[:phir2] + urel[k+1] * er[k])) : (v -= dmdt[k+1] * egb[k+1])
                else
                    v += -dmdt[k+1] * egb[k+1] + dmdt[k] * egb[k]
                end
            else
                v = (dn[k] * egn[k] * dvoln[k] - o[:do][k] * o[:ego][k] * o[:dvolo][k]) / dt -
                    dfe[k+1] + dfe[k] + pg[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) + qf[k] * dudr[k] * dvol[k]
                if k == ngrs
                    v += -dmdt[k+1] * egb[k+1] - p[:phil2] * rmu[k]
                elseif k == ngre
                    v += dmdt[k] * egb[k]
                    ip[:leobc] != 3 ? (v += rmu[k+1] * p[:phir2]) : (v -= dmdt[k+1] * egb[k+1])
                else
                    v += -dmdt[k+1] * egb[k+1] + dmdt[k] * egb[k]
                end
            end
            R[it, kb(k)] = v
        end
    end
    pres!(it, ngrs - 2)
    k = ngrs - 1
    (ip[:leibc] == 1 || ip[:llibc] == 1) ? (R[it, kb(k)] = egn[k+1] - egn[k]) : pres!(it, k)
    for k in ngre+1:ngre+2
        (ip[:leobc] in (1, 3) || ip[:llobc] > 0) ? (R[it, kb(k)] = -egn[k-1] + egn[k]) : pres!(it, k)
    end
    if c.lrad
        # ============================================================ RADIATION ENERGY
        @inbounds for k in ngrs:ngre
            v = (ern[k] * dvoln[k] - o[:ero][k] * o[:dvolo][k]) / dt +
                rmu[k+1] * fr[k+1] - rmu[k] * fr[k] +
                fedd[k] * er[k] * (rmu[k+1] * u[k+1] - rmu[k] * u[k]) +
                (1 - 3fedd[k]) * 0.5 * xmu * er[k] * aur(k) * dvol[k] +
                (cc * xke[k] * er[k] - 4p[:cpi] * xkp[k] * plf[k]) * d[k] * dvol[k]
            if k == ngrs
                v += -dmdt[k+1] * erb[k+1] - urel[k] * rmu[k] * er[k]
            elseif k == ngre
                v += dmdt[k] * erb[k]
                ip[:leobc] != 3 ? (v += urel[k+1] * rmu[k+1] * er[k]) : (v -= dmdt[k+1] * erb[k+1])
            else
                v += -dmdt[k+1] * erb[k+1] + dmdt[k] * erb[k]
            end
            R[ie, kb(k)] = v
        end
        pres!(ie, ngrs - 2)
        k = ngrs - 1; R[ie, kb(k)] = ern[k+1] - ern[k]
        ip[:lrobc] == 1 || error("radnrg: lrobc must be 1")
        k = ngre + 1; R[ie, kb(k)] = -ern[k-1] + ern[k]
        R[ie, kb(k+1)] = -ern[k] + ern[k+1]
        # ============================================================ RADIATION MOMENTUM
        qF = zeros(T, mgr)
        @inbounds for k in ngrs-1:ngre+2; qF[k] = frn[k] * 2 / (dn[k] + dn[k-1]); end
        @inbounds for k in ngrs+1:ngre
            qFo = o[:fro][k] * 2 / (o[:do][k] + o[:do][k-1])
            R[iF, kb(k)] =
                (qF[k] * (dn[k-1] * dvoln[k-1] + dn[k] * dvoln[k]) -
                 qFo * (o[:do][k-1] * o[:dvolo][k-1] + o[:do][k] * o[:dvolo][k])) / (2cc^2 * dt) -
                ((dmdt[k] + dmdt[k+1]) * frb[k] - (dmdt[k-1] + dmdt[k]) * frb[k-1]) / (2cc^2) +
                rmu[k] * (fedd[k] * er[k] - fedd[k-1] * er[k-1]) +
                0.5 * rmu[k] * (fr[k] / cc^2) * (u[k+1] - u[k-1]) +
                (0.25xmu / r[k]) * ((3fedd[k-1] - 1) * er[k-1] * dvol[k-1] + (3fedd[k] - 1) * er[k] * dvol[k]) +
                0.5 * (d[k] * dvol[k] + d[k-1] * dvol[k-1]) * avchi[k] * fr[k] / cc
        end
        k = ngrs
        ip[:lribc] in (2, 3) ? pres!(iF, k) : error("radmom: lribc must be 2 or 3")
        k = ngre + 1
        R[iF, kb(k)] = p[:cpi] * (2p[:geddr] + 1) * p[:ximr] + frn[k] - cc * p[:geddr] * ern[k-1]
        pres!(iF, ngrs - 2)
        k = ngrs - 1
        if ip[:lribc] == 2; R[iF, kb(k)] = frn[k+2] + frn[k]
        else R[iF, kb(k)] = rn[k+1]^mu * frn[k+1] - rn[k]^mu * frn[k]
        end
        k = ngre + 2; R[iF, kb(k)] = -frn[k-1] + frn[k]
    end
    # side outputs (old slopes etc.) for the Float64 evaluation
    if T === Float64
        return (; dsn, usn, egsn, ersn, frsn, qf, dudr, db, ub, egb, erb, frb)
    end
    nothing
end

total_energy_rows!(args...) = error("total-energy rows: see TitanResidE.jl")

# ------------------------------------------------------------------ Jacobian
struct JTag end
"""
    assemble!(c) -> side

Evaluate residual and the full block-pentadiagonal Jacobian with one stencil-coloured
ForwardDiff sweep and write rhs, em2..ep2 into COMMON.  7 colours (not 5): the residual
reaches ±3 blocks in one place (radmom depends on ρ(k−3) through the F/ρ interface
interpolation, ~1e-8 of the row scale); 7 colours keep that out of the ±2 entries, and the
±3 entries are dropped, exactly as TITAN's pentadiagonal system does.
"""
function assemble!(c::Ctx; check_dense = false)
    neqn = c.neqn; nb = c.ngre - c.ngrs + 5; k0 = c.ngrs - 2; N = neqn * nb
    NC = 7; P = NC * neqn
    D = ForwardDiff.Dual{JTag,Float64,P}
    xd = Vector{D}(undef, N)
    @inbounds for kb in 1:nb, j in 1:neqn
        slot = mod(kb - 1, NC) * neqn + j
        xd[(kb-1)*neqn+j] = D(0.0, ForwardDiff.Partials(ntuple(i -> i == slot ? 1.0 : 0.0, P)))
    end
    Rd = Matrix{D}(undef, neqn, nb)
    residual!(Rd, c, xd)
    E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2)); rhs = TF.var(:rhs)
    @inbounds for kb in 1:nb, i in 1:neqn
        k = k0 + kb - 1
        rhs[i, k] = ForwardDiff.value(Rd[i, kb])
        pr = ForwardDiff.partials(Rd[i, kb])
        for o in -2:2
            cb = kb + o
            for j in 1:neqn
                E[o+3][i, j, k] = (1 <= cb <= nb) ? pr[mod(cb - 1, NC) * neqn + j] : 0.0
            end
        end
    end
    # Float64 pass for side outputs (old slopes for the next step, viscosity for totnrg/writout)
    Rf = zeros(neqn, nb)
    side = residual!(Rf, c, zeros(N))
    side
end

"Write the side outputs TITAN's matgen leaves in COMMON (slopes used as old values next step)."
function store_side!(c::Ctx, s)
    a, b = c.ngrs, c.ngre
    v(sym) = TF.var(sym)
    v(:dsn)[a:b+1] .= s.dsn[a:b+1]; v(:usn)[a:b] .= s.usn[a:b]
    if c.lrad
        v(:egrsn)[a:b+1] .= s.egsn[a:b+1]; v(:ersn)[a:b+1] .= s.ersn[a:b+1]; v(:frsn)[a:b+1] .= s.frsn[a:b+1]
        v(:egrb)[a+1:b+1] .= s.egb[a+1:b+1]; v(:erb)[a+1:b+1] .= s.erb[a+1:b+1]; v(:frb)[a:b] .= s.frb[a:b]
    else
        v(:egsn)[a:b+1] .= s.egsn[a:b+1]; v(:egb)[a+1:b+1] .= s.egb[a+1:b+1]
    end
    v(:qf)[a:b] .= s.qf[a:b]; v(:dudr)[a:b] .= s.dudr[a:b]; v(:qv)[a:b] .= 0.0
    v(:db)[a+1:b+1] .= s.db[a+1:b+1]; v(:ub)[a:b] .= s.ub[a:b]
    nothing
end

"Drop-in replacement for Fortran matgen."
function matgen!(; etotal = false)
    c = context(; etotal)
    s = assemble!(c)
    store_side!(c, s)
    nothing
end

end # module
