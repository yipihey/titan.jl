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
# concretely typed accessors (type-stable loops; an untyped Array costs ~50x in update!)
vf(s::Symbol) = v(s)::Vector{Float64}
vm(s::Symbol) = v(s)::Matrix{Float64}
vi(s::Symbol) = v(s)::Vector{Int32}

"Real-exponent power. `PARITY[] = true` routes it to glibc libm (as the Fortran) for translation debugging."
const PARITY = Ref(false)
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
            sx[k, ir] = (rn[k] - ro[k]) / ro[k]
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
    s!(:tfac, tfac)
    s!(:dtime, tfac * gf(:dtime))
    nothing
end

# ------------------------------------------------------------------ update (updateh.f / updaterh.f)
"""
update!(; eos!, opac!): apply the Newton correction and refresh all derived and
time-centred quantities. `eos!`/`opac!` default to the routines installed in the
driver (Fortran via ccall, or the Julia hooks).
"""
function update!(; eos! = () -> TD.fcall(:eos_), opac! = () -> TD.fcall(:opac_))
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
            egrn[k] = egn[k] + (ern[k] / dn[k]); plfn[k] = csigr * tn[k]^4 / cpi; dplf[k] = 4.0 * plfn[k]
        end
        @inbounds for k in a:b+1; avchin[k] = 2.0 / (1.0 / chifn[k-1] + 1.0 / chifn[k]); end
        gi(:ltran) == 0 && TD.fcall(:eddfac_)              # [Fortran] formal solution, not yet ported
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
    end
    view(vf(:bri), rng) .= view(unom, rng)
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
                   xmu * cpi * cgrav * (xmn[k] / rn[k] + xmn[k+1] / rn[k+1])) - tescr0[k]
        tee += tescr[k]                                     # (TITAN: BLAS ssum, same order)
    end
    s!(:tee, tee)
    u0, r0, xm0, fedd, qf, qv = vf(:u0), vf(:r0), vf(:xm0), vf(:fedd), vf(:qf), vf(:qv)
    urel_ne = u[b+1] - (rn[b+1] - ro[b+1]) / dt
    xmdot = gf(:phir0)
    gi(:leobc) == 3 && (xmdot = d[b] * urel_ne)
    phil0 = gf(:phil0)
    bracket(k, kc) = egrt[kc] + 0.5 * (u[k]^2 - u0[k]^2) + g * (0.5 - 0.25xmu) * (r[k] - r0[k]) -
                     cgrav * cpi * 2.0 * xmu * (xm[k] / r[k] - xm0[k] / r0[k])
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
fmt(x, w, d) = lpad(replace(string(round(x; sigdigits = d + 1)), "e" => "E"), w)

end # module
