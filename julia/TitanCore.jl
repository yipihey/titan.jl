"""
    TitanCore

Phase 5: Julia translations of the remaining per-step Fortran control routines
(reset, update{h,rh}, maxdel, timstp, totnrg). They operate on the mapped COMMON
blocks exactly like the Fortran (same arrays, same index ranges), so they can be
swapped one at a time into `TitanDriver.ROUTINES` and checked in flight against the
Fortran with `PortHarness`. Printing to the .out file is not reproduced (diagnostics only);
`totnrg!` writes the .nrg record itself when an IO is given.

Not bitwise by design: integer powers use Julia `^`, and real powers use Julia `pow`.
"""
module TitanCore

import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

"Cached zero-copy views of COMMON arrays (the blocks never move, so the views stay valid)."
const VC = Dict{Symbol,Array}()
v(s::Symbol) = get!(() -> TF.var(s), VC, s)
push!(TF.ON_SWITCH, () -> empty!(VC))
# concretely typed accessors (type-stable loops; an untyped Array costs ~50x in update!)
vf(s::Symbol) = v(s)::Vector{Float64}
vm(s::Symbol) = v(s)::Matrix{Float64}
vi(s::Symbol) = v(s)::Vector{Int32}

"Real-exponent power. `PARITY[] = true` routes it to glibc libm (as the Fortran) for translation debugging."
const PARITY = Ref(false)
"""
FREE_OUTER[] (Population III, M1): the outermost interface moves under its own momentum equation (TitanResidFast.PEXT):
peq keeps its extrapolated velocity instead of imposing ur/uextr. Default false.
"""
const FREE_OUTER = Ref(false)
"""
DT_DYN[] (Population III, M1): cap the time step at DT_DYN[] × min over interfaces of the dynamical time
(r³/Gm)^½ (RHFC02 App. A2.2, there C_dyn = 1), m the interior mass. 0 (default) = TITAN's control only.
Needs gravity (cgrav > 0). TITAN's own control (stol) limits relative changes per step and lets dt grow 10× per
step; in a free-fall collapse that is far too coarse (M1a).
"""
const DT_DYN = Ref(0.0)
"Relative change, 0 where the old value is 0 (the regular centre, r = 0)."
@inline relchg(n, o) = o == 0.0 ? 0.0 : (n - o) / o
rpow(x::Float64, y::Float64) = PARITY[] ? ccall((:pow, "libm.so.6"), Float64, (Float64, Float64), x, y) : x^y
gi(s)::Int = Int(TD.g(s)); gf(s)::Float64 = Float64(TD.g(s))
s!(s, x) = TD.s!(s, x)

"BLAS isamax: 1-based index of the first element with maximal |x| among n elements of x[i0:inc:...]."
function isamax(n, x, i0, inc)
    n < 1 && return 0
    best = 1; m = abs(x[i0])
    @inbounds for i in 2:n
        a = abs(x[i0 + (i - 1) * inc])
        a > m && (best = i; m = a)
    end
    best
end

# ------------------------------------------------------------------ reset (reset.f)
const RESET_PAIRS1 = ((:tn, :to), (:dn, :do), (:rn, :ro), (:rmun, :rmuo), (:rmup1n, :rmup1o), (:rmum1n, :rmum1o),
                      (:dvoln, :dvolo), (:xmen, :xmeo), (:xmn, :xmo), (:un, :uo), (:egn, :ego), (:etn, :eto),
                      (:ern, :ero), (:frn, :fro), (:drdt, :drdto), (:dmdt, :dmdto), (:urel, :urelo), (:xnt, :xnto))
const RESET_PAIRS2 = ((:dsn, :dso), (:usn, :uso), (:pgn, :pgo), (:asn, :aso), (:egsn, :egso), (:etsn, :etso),
                      (:ersn, :erso), (:egrn, :egro), (:egrsn, :egrso), (:egtn, :egto), (:egtsn, :egtso),
                      (:egrtn, :egrto), (:egrtsn, :egrtso), (:frsn, :frso), (:chifn, :chifo), (:xken, :xkeo),
                      (:xkpn, :xkpo), (:xnen, :xneo), (:plfn, :plfo), (:cpn, :cpo), (:betan, :betao), (:delsn, :delso))
function reset!()
    a, b = gi(:ngrs), gi(:ngre)
    s!(:timen, gf(:timeo))
    for (n, o) in ((:tes, :teso), (:tew, :tewo), (:tel, :telo), (:teq, :teqo)); s!(n, gf(o)); end
    for (n, o) in RESET_PAIRS1; copyto!(view(vf(n), a-1:b+2), view(vf(o), a-1:b+2)); end
    for (n, o) in RESET_PAIRS2; copyto!(view(vf(n), a-1:b+1), view(vf(o), a-1:b+1)); end
    copyto!(view(vf(:avchin), a:b+1), view(vf(:avchio), a:b+1))
    nothing
end

# ------------------------------------------------------------------ maxdel (maxdel.f)
function maxdel!()
    a, b = gi(:ngrs), gi(:ngre); meqn = TF.PARAMS[:meqn]
    lhydr, lrad, ltur = gi(:lhydr), gi(:lrad), gi(:ltur)
    ni = b - a + 2; nc = b - a + 1
    rhs = vm(:rhs); kx = vi(:kx); dx = vf(:dx); cx = vf(:cx); rn = vf(:rn)
    fill!(kx, 0); fill!(dx, 0.0); cx[a:b] .= 0.0
    ir, im, id, iu, it, ie, iff = (gi(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    ltur == 0 || error("maxdel: ltur ≥ 1 not supported")
    flat = vec(rhs); base(i) = (a - 1) * meqn + i          # rhs(i, ngrs) in column-major
    pick!(i, n) = (kx[i] = isamax(n, flat, base(i), meqn); dx[i] = rhs[i, a + kx[i] - 1])
    if lhydr == 1 && lrad == 1
        pick!(ir, ni); pick!(im, ni); pick!(id, nc); pick!(iu, ni); pick!(it, nc); pick!(ie, nc); pick!(iff, ni)
        dmax = max(abs(dx[ir]), abs(dx[im]), abs(dx[id]), abs(dx[it]), abs(dx[ie])) + 1.0e-30
    elseif lhydr == 1 && lrad == 0
        pick!(ir, ni); pick!(im, ni); pick!(id, nc); pick!(iu, ni); pick!(it, nc)
        dmax = max(abs(dx[ir]), abs(dx[im]), abs(dx[id]), abs(dx[it])) + 1.0e-30
    elseif lhydr == 0 && lrad == 1
        pick!(ir, ni); pick!(it, nc); pick!(ie, nc); pick!(iff, ni)
        dmax = max(abs(dx[it]), abs(dx[ie])) + 1.0e-30
    else
        error("maxdel: unsupported lhydr/lrad")
    end
    s!(:dmax, dmax)
    cmax = 1.0e-30
    if gi(:lgrid) == 1
        @inbounds for k in a:b
            cx[k] = (rhs[ir, k+1] * rn[k+1] - rhs[ir, k] * rn[k]) / (rn[k+1] - rn[k])
        end
        kcx = isamax(nc, cx, a, 1)
        cmax = abs(cx[a + kcx - 1]) + 1.0e-30
    end
    s!(:cmax, cmax)
    lboos = gi(:lboos)
    lboos == 1 && return nothing
    lboos > 1 && error("maxdel01")
    sc = min(min(1.0, gf(:dtol) / dmax), min(1.0, gf(:ctol) / cmax))
    if sc < 1.0
        neqn = gi(:neqn)
        @inbounds for k in a:b+1, i in 1:neqn; rhs[i, k] *= sc; end
    end
    nothing
end

# ------------------------------------------------------------------ timstp (timstp.f)
function timstp!()
    a, b = gi(:ngrs), gi(:ngre); neqn = gi(:neqn)
    lhydr, lrad, ltur = gi(:lhydr), gi(:lrad), gi(:ltur)
    ltur == 0 || error("timstp: ltur ≥ 1 not supported")
    ni = b - a + 2; nc = b - a + 1
    sx = vm(:sx); kx = vi(:kx); dx = vf(:dx); mgr = size(sx, 1)
    kx[1:neqn] .= 0; dx[1:neqn] .= 0.0
    ir, im, id, iu, it, ie, iff = (gi(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    rn, ro, xmen, xmeo, dn, do_, un, uo, unom, tn, to = (vf(s) for s in (:rn, :ro, :xmen, :xmeo, :dn, :do, :un, :uo, :unom, :tn, :to))
    flat = vec(sx); base(i) = (i - 1) * mgr + a          # sx(ngrs, i)
    pick!(i, n) = (kx[i] = isamax(n, flat, base(i), 1); dx[i] = sx[a + kx[i] - 1, i])
    if lhydr == 1
        ern, ero, frn, fro, frnom = lrad == 1 ? (vf(:ern), vf(:ero), vf(:frn), vf(:fro), vf(:frnom)) : (nothing, nothing, nothing, nothing, nothing)
        @inbounds for k in a:b+1
            sx[k, ir] = relchg(rn[k], ro[k])
            sx[k, im] = (xmen[k] - xmeo[k]) / xmeo[k]
            sx[k, id] = (dn[k] - do_[k]) / do_[k]
            sx[k, iu] = (un[k] - uo[k]) / unom[k]
            sx[k, it] = (tn[k] - to[k]) / to[k]
            if lrad == 1
                sx[k, ie] = (ern[k] - ero[k]) / ero[k]
                sx[k, iff] = (frn[k] - fro[k]) / frnom[k]
            end
        end
        pick!(ir, ni); pick!(im, ni); pick!(id, nc); pick!(iu, ni); pick!(it, nc)
        if lrad == 1
            pick!(ie, nc); pick!(iff, ni)
            smax = max(abs(dx[ir]), abs(dx[im]), abs(dx[id]), abs(dx[it]), abs(dx[ie])) + 1.0e-30
        else
            smax = max(abs(dx[ir]), abs(dx[im]), abs(dx[id]), abs(dx[it])) + 1.0e-30
        end
    elseif lrad == 1
        ern, ero, frn, fro, frnom = vf(:ern), vf(:ero), vf(:frn), vf(:fro), vf(:frnom)
        @inbounds for k in a:b+1
            sx[k, it] = (tn[k] - to[k]) / to[k]
            sx[k, ie] = (ern[k] - ero[k]) / ero[k]
            sx[k, iff] = (frn[k] - fro[k]) / frnom[k]
        end
        pick!(it, nc); pick!(ie, nc); pick!(iff, ni)
        smax = max(abs(dx[it]), abs(dx[ie])) + 1.0e-30
    else
        error("timstp: unsupported lhydr/lrad")
    end
    s!(:smax, smax)
    x = gf(:stol) / smax
    phi1 = x > 1.0 ? sqrt(x) : x^2
    phi2 = max(0.3, min(3.0, phi1))
    tfac = max(0.1, min(10.0, sqrt(gf(:tfac) * phi2)))
    dt0 = gf(:dtime); dtn = tfac * dt0
    if DT_DYN[] > 0 && gf(:cgrav) > 0
        G4π = gf(:cgrav) * 4π; rn, xmn = vf(:rn), vf(:xmn); tdyn = Inf
        @inbounds for k in a:b+1
            (rn[k] > 0 && xmn[k] > 0) && (tdyn = min(tdyn, sqrt(rn[k]^3 / (G4π * xmn[k]))))
        end
        dtn = min(dtn, DT_DYN[] * tdyn); tfac = dtn / dt0
    end
    s!(:tfac, tfac)
    s!(:dtime, dtn)
    nothing
end

# ------------------------------------------------------------------ update (updateh.f / updaterh.f)
"""
update!(; eos!, opac!): apply the Newton correction and refresh all derived and
time-centred quantities. `eos!`/`opac!` default to the routines installed in the
driver (Fortran via ccall, or the Julia hooks).
"""
function update!(; eos! = () -> TD.EOS[](), opac! = () -> TD.OPAC[]())
    a, b = gi(:ngrs), gi(:ngre); lrad = gi(:lrad) == 1
    gi(:lhydr) == 1 && gi(:ltur) == 0 || error("update!: only lhydr=1, ltur=0 ported")
    rhs = vm(:rhs)
    ir, im, id, iu, it, ie, iff = (gi(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    rn, xmen, un, dn, tn, unom = (vf(s) for s in (:rn, :xmen, :un, :dn, :tn, :unom))
    ern, frn, frnom = vf(:ern), vf(:frn), vf(:frnom)
    mu = gi(:mu); mup1 = gi(:mup1); mum1 = gi(:mum1)
    xmup1, rxm1 = gf(:xmup1), gf(:rxm1)
    @inbounds for k in a:b+1
        rn[k] = (1.0 + rhs[ir, k]) * rn[k]
        xmen[k] = (1.0 + rhs[im, k]) * xmen[k]
        un[k] = un[k] + rhs[iu, k] * unom[k]
        lrad && (frn[k] = frn[k] + rhs[iff, k] * frnom[k])
    end
    @inbounds for k in a:b
        dn[k] = (1.0 + rhs[id, k]) * dn[k]
        tn[k] = (1.0 + rhs[it, k]) * tn[k]
        lrad && (ern[k] = (1.0 + rhs[ie, k]) * ern[k])
    end
    _ghost_and_derived!(:update, eos!, opac!)
end

"Eddington-factor routine called by peq (always) and update (ltran = 0). [Fortran until ported]"
const EDDFAC = Ref{Function}(() -> TD.fcall(:eddfac_))

"""
Shared tail of update{h,rh} and peq{h,rh}: ghost-zone boundary values, geometry, EOS/opacity,
Eddington factors and time-centred / zeroed arrays. `mode` ∈ (:update, :peq) selects the few
places where the Fortran routines differ (eddfac call, plfn operation order, extra zeroing).
"""
function _ghost_and_derived!(mode::Symbol, eos!, opac!)
    a, b = gi(:ngrs), gi(:ngre); lrad = gi(:lrad) == 1
    rn, xmen, un, dn, tn, unom = (vf(s) for s in (:rn, :xmen, :un, :dn, :tn, :unom))
    ern, frn, frnom = vf(:ern), vf(:frn), vf(:frnom)
    mu = gi(:mu); mup1 = gi(:mup1); mum1 = gi(:mum1)
    xmup1, rxm1 = gf(:xmup1), gf(:rxm1)
    leibc, llibc, lribc = gi(:leibc), gi(:llibc), gi(:lribc)
    leobc, llobc, lrobc = gi(:leobc), gi(:llobc), gi(:lrobc)
    # ---- inner boundary
    k = a - 1
    if leibc == 1
        tn[k] = tn[k+1]; dn[k] = dn[k+1]
        rn[k] = rn[k+1] * 2.0 - rn[k+2]; xmen[k] = xmen[k+1] * 2.0 - xmen[k+2]; un[k] = -un[k+2]
    end
    if leibc == 2
        tn[k] = gf(:tl); dn[k] = gf(:dl); rn[k] = rn[k+1] - gf(:delrl)
        xmen[k] = xmen[k+1] + gf(:dl) * rxm1 * (rn[k+1]^mup1 - rn[k]^mup1); un[k] = gf(:ul)
    end
    if llibc == 1
        tn[k] = tn[k+1]; dn[k] = dn[k+1]
        rn[k] = rpow(rn[k+1]^mup1 - xmup1 * gf(:delml) / dn[k], rxm1)
        xmen[k] = xmen[k+1] + gf(:delml); un[k] = un[k+1]
    end
    if lrad
        ern[k] = ern[k+1]
        lribc == 1 && (frn[k] = +frn[k+1])
        lribc == 2 && (frn[k] = -frn[k+2])
        lribc == 3 && (frn[k] = frn[k+1] * (rn[k+1] / rn[k])^mu)
    end
    # ---- outer boundary
    k = b + 1
    if leobc == 1
        tn[k] = tn[k-1]; tn[k+1] = tn[k]; dn[k] = dn[k-1]; dn[k+1] = dn[k]
        rn[k+1] = rn[k] * 2.0 - rn[k-1]; xmen[k+1] = xmen[k] * 2.0 - xmen[k-1]; un[k+1] = -un[k-1]
    end
    if leobc == 2
        tr, dr = gf(:tr), gf(:dr)
        tn[k] = tr; tn[k+1] = tr; dn[k] = dr; dn[k+1] = dr
        rn[k+1] = rn[k] + gf(:delrr); xmen[k+1] = xmen[k] - dr * rxm1 * (rn[k+1]^mup1 - rn[k]^mup1); un[k+1] = gf(:ur)
    end
    if leobc == 3
        tn[k] = tn[k-1]; tn[k+1] = tn[k]; dn[k] = dn[k-1]; dn[k+1] = dn[k]
        rn[k+1] = rn[k] + gf(:delrr); xmen[k+1] = xmen[k] - dn[k] * rxm1 * (rn[k+1]^mup1 - rn[k]^mup1); un[k+1] = un[k]
    end
    if llobc > 0
        tn[k] = tn[k-1]; tn[k+1] = tn[k]; dn[k] = dn[k-1]; dn[k+1] = dn[k]
        rn[k+1] = rpow(rn[k]^mup1 + xmup1 * gf(:delmr) / dn[k], rxm1); xmen[k+1] = xmen[k] - gf(:delmr); un[k+1] = un[k]
    end
    if lrad
        lrobc == 1 && (ern[k] = ern[k-1]; ern[k+1] = ern[k]; frn[k+1] = +frn[k])
        lrobc == 2 && (ern[k] = ern[k-1]; ern[k+1] = ern[k-2]; frn[k+1] = -frn[k-1])
        lrobc == 3 && (ern[k] = ern[k-1]; ern[k+1] = ern[k]; frn[k+1] = frn[k] * (rn[k] / rn[k+1])^mu)
    end
    # ---- geometry
    rmun, rmup1n, rmum1n, xmn, dvoln = (vf(s) for s in (:rmun, :rmup1n, :rmum1n, :xmn, :dvoln))
    xmtot = gf(:xmtot)
    @inbounds for k in a-1:b+2
        rmun[k] = rn[k]^mu; rmup1n[k] = rn[k]^mup1; rmum1n[k] = rn[k]^mum1; xmn[k] = xmtot - xmen[k]
    end
    @inbounds for k in a-1:b+1; dvoln[k] = (rmup1n[k+1] - rmup1n[k]) / xmup1; end
    eos!()
    egn = vf(:egn)
    if lrad
        opac!()
        egrn, plfn, dplf, chifn, avchin = (vf(s) for s in (:egrn, :plfn, :dplfdltn, :chifn, :avchin))
        csigr, cpi = gf(:csigr), gf(:cpi)
        @inbounds for k in a-1:b+2
            egrn[k] = egn[k] + (ern[k] / dn[k])
            plfn[k] = mode === :peq ? csigr / cpi * tn[k]^4 : csigr * tn[k]^4 / cpi
            dplf[k] = 4.0 * plfn[k]
        end
        @inbounds for k in a:b+1; avchin[k] = 2.0 / (1.0 / chifn[k-1] + 1.0 / chifn[k]); end
        (mode === :peq || gi(:ltran) == 0) && EDDFAC[]()
    end
    # ---- time-centred quantities
    θ = gf(:thet); dt = gf(:dtime)
    r, u, xm, xme, drdt, dmdt, urel, rmu, rmup1, rmum1 =
        (vf(s) for s in (:r, :u, :xm, :xme, :drdt, :dmdt, :urel, :rmu, :rmup1, :rmum1))
    ro, uo, xmo, xmeo = vf(:ro), vf(:uo), vf(:xmo), vf(:xmeo)
    fr, fro = vf(:fr), vf(:fro)
    @inbounds for k in a-1:b+2
        r[k] = θ * rn[k] + (1.0 - θ) * ro[k]
        u[k] = θ * un[k] + (1.0 - θ) * uo[k]
        xm[k] = θ * xmn[k] + (1.0 - θ) * xmo[k]
        xme[k] = θ * xmen[k] + (1.0 - θ) * xmeo[k]
        lrad && (fr[k] = θ * frn[k] + (1.0 - θ) * fro[k])
        dmdt[k] = -(xmen[k] - xmeo[k]) / dt
        drdt[k] = (rn[k] - ro[k]) / dt
        urel[k] = u[k] - drdt[k]
        rmu[k] = r[k]^mu; rmup1[k] = r[k]^mup1; rmum1[k] = r[k]^mum1
        if lrad
            frnom[k] = abs(frn[k]) + 1.0e-30           # (the Teff/xlum choice is overwritten in the Fortran too)
            unom[k] = abs(un[k]) + 1.0e-30
        else
            unom[k] = max(abs(un[k]), 1.0)
        end
    end
    d, t, eg, pg, as_, dvol = (vf(s) for s in (:d, :t, :eg, :pg, :as, :dvol))
    do_, to, ego, pgn, pgo, asn, aso = (vf(s) for s in (:do, :to, :ego, :pgn, :pgo, :asn, :aso))
    @inbounds for k in a-1:b+1
        d[k] = θ * dn[k] + (1.0 - θ) * do_[k]
        t[k] = θ * tn[k] + (1.0 - θ) * to[k]
        eg[k] = θ * egn[k] + (1.0 - θ) * ego[k]
        pg[k] = θ * pgn[k] + (1.0 - θ) * pgo[k]
        as_[k] = θ * asn[k] + (1.0 - θ) * aso[k]
        dvol[k] = (rmup1[k+1] - rmup1[k]) / xmup1
    end
    if lrad
        for (c, n, o) in ((:er, :ern, :ero), (:egr, :egrn, :egro), (:chif, :chifn, :chifo), (:xke, :xken, :xkeo),
                          (:xkp, :xkpn, :xkpo), (:plf, :plfn, :plfo))
            cc, nn, oo = vf(c), vf(n), vf(o)
            @inbounds for k in a-1:b+1; cc[k] = θ * nn[k] + (1.0 - θ) * oo[k]; end
        end
        chif, avchi = vf(:chif), vf(:avchi)
        @inbounds for k in a:b+1; avchi[k] = 2.0 / (1.0 / chif[k-1] + 1.0 / chif[k]); end
    end
    # ---- zeroed / aliased arrays (turbulence and, for pure hydro, radiation)
    rng = a-1:b+2
    zero!(s) = (view(vf(s), rng) .= 0.0)
    foreach(zero!, (:eto, :et, :etn, :etso, :etb, :etsn))
    mode === :peq && foreach(zero!, (:egto, :egt, :egtn, :egtso, :egtb, :egtsn))
    if lrad
        foreach(zero!, (:egso, :egb, :egsn))
    else
        foreach(zero!, (:ero, :er, :ern, :erso, :erb, :ersn))
        view(vf(:egro), rng) .= view(vf(:ego), rng); view(vf(:egr), rng) .= view(vf(:eg), rng)
        view(vf(:egrn), rng) .= view(vf(:egn), rng)
        foreach(zero!, (:egrso, :egrb, :egrsn))
    end
    view(vf(:egrto), rng) .= view(vf(:egro), rng); view(vf(:egrt), rng) .= view(vf(:egr), rng)
    view(vf(:egrtn), rng) .= view(vf(:egrn), rng)
    foreach(zero!, (:egrtso, :egrtb, :egrtsn))
    if !lrad
        foreach(zero!, (:fro, :fr, :frn, :frso, :frb, :frsn, :chifo, :chif, :chifn, :avchio, :avchi, :avchin,
                        :xkeo, :xke, :xken, :xkpo, :xkp, :xkpn, :plfo, :plf, :plfn, :fedd, :frnom))
        mode === :peq && foreach(zero!, (:xneo, :xne, :xnen))
    end
    view(vf(:bri), rng) .= view(unom, rng)
    nothing
end

# ------------------------------------------------------------------ peq (peq.f → peqh.f / peqrh.f)
const PEQ_OLD1 = ((:to, :tn), (:do, :dn), (:ro, :rn), (:rmuo, :rmun), (:rmup1o, :rmup1n), (:rmum1o, :rmum1n),
                  (:dvolo, :dvoln), (:xmeo, :xmen), (:xmo, :xmn), (:uo, :un), (:ego, :egn), (:drdto, :drdt),
                  (:dmdto, :dmdt), (:urelo, :urel), (:xnto, :xnt))
const PEQ_OLD2_H  = ((:dso, :dsn), (:uso, :usn), (:pgo, :pgn), (:aso, :asn), (:egso, :egsn))
const PEQ_OLD2_RH = ((:dso, :dsn), (:uso, :usn), (:pgo, :pgn), (:aso, :asn), (:erso, :ersn), (:egso, :egsn),
                     (:egro, :egrn), (:egrso, :egrsn), (:chifo, :chifn), (:frso, :frsn), (:xkeo, :xken),
                     (:xkpo, :xkpn), (:xneo, :xnen), (:plfo, :plfn))

"Prepare a time step (peq{h,rh}): new → old, extrapolate r and m, boundary values, derived quantities."
function peq!(; eos! = () -> TD.EOS[](), opac! = () -> TD.OPAC[]())
    a, b = gi(:ngrs), gi(:ngre); lrad = gi(:lrad) == 1
    gi(:lhydr) == 1 && gi(:ltur) == 0 || error("peq!: only lhydr=1, ltur=0 ported")
    # ---- replace old by new
    s!(:timeo, gf(:timen)); s!(:timen, gf(:timeo) + gf(:dtime))
    for (o, n) in ((:teso, :tes), (:tewo, :tew), (:telo, :tel), (:teqo, :teq)); s!(o, gf(n)); end
    for (o, n) in PEQ_OLD1; copyto!(view(vf(o), a-1:b+2), view(vf(n), a-1:b+2)); end
    lrad && for (o, n) in ((:ero, :ern), (:fro, :frn)); copyto!(view(vf(o), a-1:b+2), view(vf(n), a-1:b+2)); end
    for (o, n) in (lrad ? PEQ_OLD2_RH : PEQ_OLD2_H); copyto!(view(vf(o), a-1:b+1), view(vf(n), a-1:b+1)); end
    lrad && copyto!(view(vf(:avchio), a:b+1), view(vf(:avchin), a:b+1))
    # ---- extrapolate radii and masses to the new time level (with TITAN's halving retry)
    rn, xmen, ro, xmeo, un = vf(:rn), vf(:xmen), vf(:ro), vf(:xmeo), vf(:un)
    drdt, dmdt = vf(:drdt), vf(:dmdt)
    mu = gi(:mu)
    leibc, llibc, leobc, llobc = gi(:leibc), gi(:llibc), gi(:leobc), gi(:llobc)
    jext = 1; nnext = gi(:next); s!(:jext, jext)
    while true
        leibc > 0 && (drdt[a] = 0.0)
        leibc > 0 && (dmdt[a] = -rn[a]^mu * gf(:phil0))
        leobc > 0 && (drdt[b+1] = 0.0)
        leobc < 3 && (dmdt[b+1] = -rn[b+1]^mu * gf(:phir0))
        llibc > 0 && (drdt[a] = un[a]; dmdt[a] = 0.0)
        llobc > 0 && (drdt[b+1] = un[b+1]; dmdt[b+1] = 0.0)
        dt = gf(:dtime)
        @inbounds for k in a:b+1
            rn[k] = ro[k] + drdt[k] * dt
            xmen[k] = xmeo[k] - dmdt[k] * dt
        end
        mono = true
        @inbounds for k in a:b
            (rn[k+1] <= rn[k] || xmen[k+1] > xmen[k]) && (mono = false; break)
        end
        mono && break
        s!(:tfac, 0.5); s!(:dtime, 0.5 * gf(:dtime)); s!(:timen, gf(:timeo) + gf(:dtime))
        jext += 1; s!(:jext, jext)
        jext <= nnext || error("peq: grid extrapolation nonmonotonic after $nnext attempts")
    end
    # ---- boundary values of u and F
    leibc <= 2 && (un[a] = gf(:ul)); llibc == 1 && (un[a] = gf(:uextl))
    FREE_OUTER[] || (leobc <= 2 && (un[b+1] = gf(:ur)); llobc == 1 && (un[b+1] = gf(:uextr)))
    if lrad
        frn, ern = vf(:frn), vf(:ern)
        cpi, cc, csigr, xmu = gf(:cpi), gf(:cc), gf(:csigr), gf(:xmu)
        lribc, lrobc = gi(:lribc), gi(:lrobc)
        lribc == 1 && (frn[a] = +(2.0 * gf(:geddl) + 1.0) * cpi * gf(:xipl) - cc * gf(:geddl) * ern[a])
        lribc == 2 && (frn[a] = 0.0)
        lribc == 3 && (frn[a] = +csigr * gf(:teff)^4 * (1.0 - 0.5xmu) + 0.5 * gf(:xlum) * xmu / (4.0 * cpi * rn[a]^mu))
        lrobc == 1 && (frn[b+1] = -(2.0 * gf(:geddr) + 1.0) * cpi * gf(:ximr) + cc * gf(:geddr) * ern[b])
        lrobc == 2 && (frn[b+1] = 0.0)
        lrobc == 3 && (frn[b+1] = frn[b] * (rn[b] / rn[b+1])^mu)
    end
    _ghost_and_derived!(:peq, eos!, opac!)
end

# ------------------------------------------------------------------ eddfac (eddfac.f, eddf.f, eddg.f)
"""
`FIX_EDDFAC_PLANAR[] = true` fixes a TITAN 1.2 bug in the planar (μ = 0) reflecting-inner-boundary
branch (lribc = 2): the return sweep `do 9 kp = ngre, ngrs, -1` uses the stale outer index k = ngre+1
instead of kp. Off by default (faithful translation). Not used by SUBC/SUPC, which are spherical.
"""
const FIX_EDDFAC_PLANAR = Ref(false)
"""
`EDD_POSITIVE[]` (spherical formal solution): keep the linear-in-τ source reconstruction and the ray intensities
non-negative, and use f = 1/3 if the mean intensity vanishes (the moments are accumulated over decreasing μ, so
xmom0 < 0 is the normal sign). Without it a steep source-function jump made
intensities negative, f = K/J came out 0 in a run of zones above the accretion shock, and the radiation field there
broke up zone by zone (run K1 stopped at 9.5e7 s). On whenever a tabulated material is installed
(`TitanDriver.EDD_MATERIAL_SOURCE`, default :auto), or set true/false explicitly; off for TITAN's built-in physics (bit-reproducible benchmarks).
"""
const EDD_POSITIVE = Ref{Any}(:auto)        # :auto (follow EDD_MATERIAL_SOURCE), true or false
edpos() = EDD_POSITIVE[] === :auto ? TD.EDD_MATERIAL_SOURCE[] : EDD_POSITIVE[]::Bool

"eddf.f: accumulate the 0th and 2nd angular moments of xi over ang[1:n] into xmom0[k], xmom2[k]."
@inline function _eddf!(n, ang, xi, xmom0, xmom2, k)
    s0 = xmom0[k]; s2 = xmom2[k]
    @inbounds for i in 1:n-1
        a, b = ang[i], ang[i+1]; x, y = xi[i], xi[i+1]
        s0 += 0.5 * (x + y) * (b - a)
        s2 += (x * b - y * a) * (a^2 + a * b + b^2) / 3.0 + (y - x) * (a + b) * (a^2 + b^2) / 4.0
    end
    xmom0[k] = s0; xmom2[k] = s2
    nothing
end
"eddg.f: first moment / zeroth moment of xi over ang[1:n]."
@inline function _eddg(n, ang, xi)
    s0 = 0.0; s1 = 0.0
    @inbounds for i in 1:n-1
        a, b = ang[i], ang[i+1]; x, y = xi[i], xi[i+1]
        s0 += 0.5 * (x + y) * (b - a)
        s1 += (x * b - y * a) * (a + b) / 2.0 + (y - x) * (a^2 + a * b + b^2) / 3.0
    end
    s1 / s0
end
cvmgp(a, b, c) = c < 0.0 ? b : a
cvmgm(a, b, c) = c >= 0.0 ? b : a

"Eddington factors f = K/J and boundary factors g from a formal solution (eddfac.f)."
function eddfac!()
    a, b = gi(:ngrs), gi(:ngre); fedd = vf(:fedd)
    if gi(:ltran) > 1
        view(fedd, a-1:b+2) .= 1.0 / 3.0; s!(:geddl, 0.5); s!(:geddr, 0.5); return nothing
    end
    mu = gi(:mu)
    sf, dtau, dsf, bsf, asf = vf(:sf), vf(:dtau), vf(:dsf), vf(:bsf), vf(:asf)
    chifn, dn, xnen, plfn, ern, frn, rn = vf(:chifn), vf(:dn), vf(:xnen), vf(:plfn), vf(:ern), vf(:frn), vf(:rn)
    xim, xip, xi0, ang = vf(:xim), vf(:xip), vf(:xi0), vf(:ang)
    xmom0, xmom2 = vf(:xmom0), vf(:xmom2)
    csige, cc, cpi, eps = gf(:csige), gf(:cc), gf(:cpi), gf(:epsedf)
    nang = gi(:nang); lribc, lrobc = gi(:lribc), gi(:lrobc)
    xkpn = vf(:xkpn); matsrc = TD.EDD_MATERIAL_SOURCE[]
    source!(k) = (sf[k] = matsrc ?
                      (xkpn[k] * dn[k] * plfn[k] + csige * xnen[k] * ern[k] * cc / (4.0 * cpi)) / (xkpn[k] * dn[k] + csige * xnen[k]) :
                      ((chifn[k] * dn[k] - csige * xnen[k]) * plfn[k] + csige * xnen[k] * ern[k] * cc / (4.0 * cpi)) /
                          (chifn[k] * dn[k]);
                  dtau[k] = chifn[k] * dn[k] * (rn[k+1] - rn[k]))
    function slopes!(withasf)
        @inbounds for k in a:b-1; dsf[k] = sf[k] - sf[k+1]; end
        @inbounds for k in a+1:b-1
            dxm = 0.5 * (dtau[k-1] + dtau[k]); dxp = 0.5 * (dtau[k] + dtau[k+1])
            sm = dxm * dsf[k] + dxp * dsf[k-1]; prod = dsf[k-1] * dsf[k]
            denom = cvmgm(sm - eps, sm + eps, sm)
            bsf[k] = cvmgp(2.0 * prod / denom, 0.0, prod)
            withasf && (asf[k] = sf[k] - 0.5 * dtau[k] * bsf[k])
        end
        bsf[a] = bsf[a+1]; bsf[b] = bsf[b-1]
    end
    if mu == 0
        for k in a:b; source!(k); end
        slopes!(true)
        asf[a] = sf[a] - 0.5 * dtau[a] * bsf[a]; asf[b] = sf[b] - 0.5 * dtau[b] * bsf[b]
        if lribc == 2
            cjk, ejk, sjk = vf(:cjk), vf(:ejk), vf(:sjk)
            for j in 2:nang; cjk[j] = 0.0; ejk[j] = -1.0; sjk[j] = 0.0; end
            for k in a:b, j in 2:nang
                dt = dtau[k] / ang[j]; ex = exp(-dt)
                cjk[j] = -ex; ejk[j] = -cjk[j] * ejk[j]
                sjk[j] = -cjk[j] * sjk[j] + asf[k] * (1.0 - ex) + bsf[k] * ang[j] * (1.0 - (1.0 + dt) * ex)
            end
            for j in 2:nang; cjk[j] = -1.0; ejk[j] = -cjk[j] * ejk[j]; sjk[j] = -cjk[j] * sjk[j]; end
            kstale = b + 1                                   # Fortran DO index after `do 6 k = ngrs, ngre`
            for kp in b:-1:a, j in 2:nang
                k = FIX_EDDFAC_PLANAR[] ? kp : kstale
                dt = dtau[k] / ang[j]; ex = exp(-dt)
                cjk[j] = -ex; ejk[j] = -cjk[j] * ejk[j]
                sjk[j] = -cjk[j] * sjk[j] + asf[k] * (1.0 - ex) + bsf[k] * ang[j] * (dt - 1.0 + ex)
            end
            for j in 2:nang; xim[j] = sjk[j] / (1.0 + ejk[j]); end
            xim[1] = xim[2]
        end
        lribc == 1 && (xim[1:nang] .= gf(:xipl))
        if lribc == 3
            for j in 1:nang
                xim[j] = (cc * ern[a] + 1.5 * frn[a] * dtau[a]) / (4.0 * cpi) + 3.0 * ang[j] * frn[a] / (4.0 * cpi)
            end
        end
        for k in a:b
            for j in 2:nang
                dt = 0.5 * dtau[k] / ang[j]; ex = exp(-dt)
                xi0[j] = sf[k] * (1.0 - ex) + bsf[k] * ang[j] * (1.0 - (1.0 + dt) * ex) + xim[j] * ex
                dt = 2.0 * dt; ex = exp(-dt)
                xip[j] = asf[k] * (1.0 - ex) + bsf[k] * ang[j] * (1.0 - (1.0 + dt) * ex) + xim[j] * ex
            end
            xi0[1] = xi0[2]; xip[1] = xip[2]
            for j in 1:nang; xim[j] = xip[j]; end
            xmom0[k] = 0.0; xmom2[k] = 0.0
            _eddf!(nang, ang, xi0, xmom0, xmom2, k)
            k == b && s!(:geddr, _eddg(nang, ang, xip))
        end
        lrobc == 1 && (xip[1:nang] .= gf(:ximr))
        for k in b:-1:a
            for j in 2:nang
                dt = 0.5 * dtau[k] / ang[j]; ex = exp(-dt)
                xi0[j] = asf[k] * (1.0 - ex) + bsf[k] * ang[j] * (dt - 1.0 + ex) + xip[j] * ex
                dt = 2.0 * dt; ex = exp(-dt)
                xim[j] = asf[k] * (1.0 - ex) + bsf[k] * ang[j] * (dt - 1.0 + ex) + xip[j] * ex
            end
            xi0[1] = xi0[2]; xim[1] = xim[2]
            for j in 1:nang; xip[j] = xim[j]; end
            _eddf!(nang, ang, xi0, xmom0, xmom2, k)
            fedd[k] = max(0.0, min(xmom2[k] / xmom0[k], 1.0))
            k == a && s!(:geddl, _eddg(nang, ang, xim))
            s!(:geddl, max(0.5, min(gf(:geddl), 1.0)))
        end
        return nothing
    end
    mu == 2 || error("eddfac!: mu must be 0 or 2")
    rp, sl, sr, p, xang = vf(:rp), vf(:sl), vf(:sr), vf(:p), vf(:xang)
    ncor = gi(:ncor)
    for k in b:-1:a; source!(k); end
    slopes!(false)
    @inbounds for k in a:b
        l = 2 * (k - a) + 1
        rp[l] = rn[k]; rp[l+1] = 0.5 * (rn[k] + rn[k+1])
        sl[l] = sf[k] + 0.5 * dtau[k] * bsf[k]; sr[l] = sf[k]
        sl[l+1] = sf[k]; sr[l+1] = sf[k] - 0.5 * dtau[k] * bsf[k]
        if edpos()     # a source function reconstructed linearly in τ must not go negative
            sl[l] = max(sl[l], 0.0); sr[l+1] = max(sr[l+1], 0.0)
        end
    end
    lmax = 2 * (b - a + 1) + 1
    rp[lmax] = rn[b+1]
    for j in 1:ncor; p[j] = Float64(j - 1) * rp[1] / Float64(ncor); end
    nray = (ncor + 1) + 2 * (b + 1 - a)
    for j in ncor+1:nray; p[j] = rp[j - ncor]; end
    lrobc == 1 && (xim[1:nray] .= gf(:ximr))
    # inward sweep along tangent rays
    for l in lmax-1:-1:1
        k = a + (l - 1) ÷ 2; nj = ncor + l
        @inbounds for j in 1:nj
            xjl = sqrt(rp[l]^2 - p[j]^2); xjlp1 = sqrt(rp[l+1]^2 - p[j]^2)
            dtjl = dn[k] * chifn[k] * (xjlp1 - xjl); ex = exp(-dtjl)
            xip[j] = (1.0 - ex) * sl[l] + (1.0 - (1.0 + dtjl) * ex) * (sr[l] - sl[l]) / dtjl + ex * xim[j]
            edpos() && (xip[j] = max(xip[j], 0.0))
            xang[j] = sqrt(1.0 - p[j]^2 / rp[l]^2)
        end
        for j in 1:nj; xim[j] = xip[j]; end
        if iseven(l)
            k = a + (l - 1) ÷ 2
            xmom0[k] = 0.0; xmom2[k] = 0.0
            _eddf!(nj, xang, xip, xmom0, xmom2, k)
        end
        l == 1 && s!(:geddl, _eddg(nj, xang, xip))
    end
    if lribc == 1
        for j in 1:ncor+1; xang[j] = sqrt(1.0 - p[j]^2 / rp[1]^2); xim[j] = gf(:xipl); end
    end
    if lribc == 3                                    # (uses the planar quadrature ang(j), as in TITAN 1.2)
        for j in 1:ncor+1
            xang[j] = sqrt(1.0 - p[j]^2 / rp[1]^2)
            xim[j] = (cc * ern[a] + 1.5 * frn[a] * dtau[a]) / (4.0 * cpi) + 3.0 * ang[j] * frn[a] / (4.0 * cpi)
        end
    end
    # outward sweep
    for l in 2:lmax
        k = a + (l - 2) ÷ 2; nj = ncor + l - 1
        @inbounds for j in 1:nj
            xjl = sqrt(rp[l]^2 - p[j]^2); xjlm1 = sqrt(rp[l-1]^2 - p[j]^2)
            dtjl = dn[k] * chifn[k] * (xjl - xjlm1); ex = exp(-dtjl)
            xip[j] = (1.0 - ex) * sr[l-1] + (1.0 - (1.0 + dtjl) * ex) * (sl[l-1] - sr[l-1]) / dtjl + ex * xim[j]
            edpos() && (xip[j] = max(xip[j], 0.0))
            xang[j] = sqrt(1.0 - p[j]^2 / rp[l]^2)
        end
        for j in 1:nj; xim[j] = xip[j]; end
        if iseven(l)
            k = a + (l - 1) ÷ 2
            xang[nj+1] = 0.0
            _eddf!(nj + 1, xang, xim, xmom0, xmom2, k)
            fedd[k] = edpos() && !(abs(xmom0[k]) > 0 && isfinite(xmom2[k] / xmom0[k])) ? 1.0 / 3.0 :
                      max(0.0, min(xmom2[k] / xmom0[k], 1.0))
        end
        l == lmax && s!(:geddr, _eddg(nj, xang, xip))
        s!(:geddr, max(0.5, min(gf(:geddr), 1.0)))
    end
    nothing
end

# ------------------------------------------------------------------ totnrg (totnrg.f)
"Energy bookkeeping. If `io` is given, appends the .nrg record (same columns as TITAN)."
function totnrg!(io = nothing)
    a, b = gi(:ngrs), gi(:ngre)
    mu, mup1, mum1 = gi(:mu), gi(:mup1), gi(:mum1); xmu, xmup1 = gf(:xmu), gf(:xmup1)
    xmtot, g, cgrav, cpi, dt = gf(:xmtot), gf(:g), gf(:cgrav), gf(:cpi), gf(:dtime)
    rn, ro, xmen, xmeo, un, uo, frn, fro = (vf(s) for s in (:rn, :ro, :xmen, :xmeo, :un, :uo, :frn, :fro))
    r, xm, u, fr, fc, ft, rmu, rmup1, rmum1 = (vf(s) for s in (:r, :xm, :u, :fr, :fc, :ft, :rmu, :rmup1, :rmum1))
    @inbounds for k in a:b+1
        r[k] = (rn[k] + ro[k]) / 2.0; xm[k] = xmtot - (xmen[k] + xmeo[k]) / 2.0
        u[k] = (un[k] + uo[k]) / 2.0; fr[k] = (frn[k] + fro[k]) / 2.0; fc[k] = 0.0; ft[k] = 0.0
        rmu[k] = r[k]^mu; rmup1[k] = r[k]^mup1; rmum1[k] = r[k]^mum1
    end
    d, eg, pg, er, et, egrt, dvol = (vf(s) for s in (:d, :eg, :pg, :er, :et, :egrt, :dvol))
    for (c, n, o) in ((d, :dn, :do), (eg, :egn, :ego), (pg, :pgn, :pgo), (er, :ern, :ero), (et, :etn, :eto), (egrt, :egrtn, :egrto))
        nn, oo = vf(n), vf(o)
        @inbounds for k in a:b; c[k] = 0.5 * (nn[k] + oo[k]); end
    end
    @inbounds for k in a:b; dvol[k] = (rmup1[k+1] - rmup1[k]) / xmup1; end
    dn, dvoln, egrtn, xmn, tescr, tescr0 = (vf(s) for s in (:dn, :dvoln, :egrtn, :xmn, :tescr, :tescr0))
    tee = 0.0
    @inbounds for k in a:b
        tescr[k] = dn[k] * dvoln[k] * (egrtn[k] + 0.25 * (un[k]^2 + un[k+1]^2) +
                   (0.5 - 0.25xmu) * g * (rn[k] + rn[k+1]) -
                   xmu * cpi * cgrav * (mor(xmn[k], rn[k]) + mor(xmn[k+1], rn[k+1]))) - tescr0[k]
        tee += tescr[k]                                     # (TITAN: BLAS ssum, same order)
    end
    s!(:tee, tee)
    u0, r0, xm0, fedd, qf, qv = vf(:u0), vf(:r0), vf(:xm0), vf(:fedd), vf(:qf), vf(:qv)
    urel_ne = u[b+1] - (rn[b+1] - ro[b+1]) / dt
    xmdot = gf(:phir0)
    gi(:leobc) == 3 && (xmdot = d[b] * urel_ne)
    phil0 = gf(:phil0)
    bracket(k, kc) = egrt[kc] + 0.5 * (u[k]^2 - u0[k]^2) + g * (0.5 - 0.25xmu) * (r[k] - r0[k]) -
                     cgrav * cpi * 2.0 * xmu * (mor(xm[k], r[k]) - mor(xm0[k], r0[k]))
    tes = gf(:teso) + dt * rmu[b+1] * xmdot * bracket(b + 1, b)
    tes = tes - dt * rmu[a] * phil0 * bracket(a, a)
    tew = gf(:tewo) + dt * (rmu[b+1] * u[b+1] * (pg[b] + fedd[b] * er[b] + (2.0 / 3.0) * d[b] * et[b]) -
                            rmu[a] * u[a] * (pg[a] + fedd[a] * er[a] + (2.0 / 3.0) * d[a] * et[a]))
    tel = gf(:telo) + dt * (rmu[b+1] * (fr[b+1] + fc[b+1] + ft[b+1]) - rmu[a] * (fr[a] + fc[a] + ft[a]))
    teq = gf(:teqo) + dt * (rmu[b+1] * u[b+1] * (qf[b+1] + qv[b+1]) - rmu[a] * u[a] * (qf[a] + qv[a]))
    tetot = tee + tes + tew + tel - teq
    s!(:tes, tes); s!(:tew, tew); s!(:tel, tel); s!(:teq, teq); s!(:tetot, tetot)
    if io !== nothing
        println(io, rpad(gi(:jstep), 5), join((fmt(x, 12, 4) for x in (gf(:timen), dt)), ""),
                fmt(tetot, 16, 8), join((fmt(x, 12, 4) for x in (tee, tes, tew, tel, teq)), ""))
    end
    nothing
end
"m/r, 0 at the regular centre (r = 0, m = 0)."
mor(m, r) = r == 0.0 ? 0.0 : m / r
fmt(x, w, d) = lpad(replace(string(round(x; sigdigits = d + 1)), "e" => "E"), w)

end # module
