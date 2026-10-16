"""
    QuasiStarSetup

Phase 7d: a TITAN initial state for the outer shell of a quasi-star, built from a steady
QuasiStarInfall envelope (`scripts/export_profile.jl` there writes the CSV).

Domain: r_in (a fraction of the photospheric radius R; the mass inside r_in, BH + inner envelope,
acts as a point mass) out through the photosphere into a grey hydrostatic atmosphere of `nH`
density e-foldings. Eulerian grid (lgrid = 2) for the hydrostatic test: zero-flux walls at both
ends (leibc = leobc = 1), imposed luminosity at the inner boundary (lribc = 3, xlum = L, all of it
radiative at the wall), transmitting radiation outer boundary (lrobc = 1, no incident intensity).
Material: Z0Gas (Saha H/He, tabulated Rosseland/Planck means), local MLT flux on.

TITAN conventions used here: masses are per steradian (gravity = 4π G xm / r² with xm = xmtot − xmen),
F is the physical flux (erg cm⁻² s⁻¹) on interfaces, ρ, T, E_r are zone averages.
TITAN's physical constants are replaced by those of the Z0 physics (and the steady model), because
the outer envelope has Γ ≈ 1: the net force on the gas is ~1 % of gravity, so a 2·10⁻⁴ mismatch in G
or a_rad would already be a few per cent of it.
"""
module QuasiStarSetup

using Printf
import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..TitanMaterials, ..Z0Physics
const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
const TM = TitanMaterials; const Z0 = Z0Physics

# ------------------------------------------------------------------ steady profile
"Read a profile CSV from QuasiStarInfall/scripts/export_profile.jl (outer → inner)."
function read_profile(path)
    lines = readlines(path)
    hdr = Dict{Symbol,Float64}()
    for kv in split(lstrip(lines[1], ['#', ' ']))
        k, v = split(kv, "="); hdr[Symbol(k)] = parse(Float64, v)
    end
    names = Symbol.(split(lines[2], ","))
    data = [parse.(Float64, split(l, ",")) for l in lines[3:end] if !isempty(l)]
    cols = Dict(n => [d[i] for d in data] for (i, n) in enumerate(names))
    return (hdr = hdr, cols = cols)
end

"Gas density from gas pressure and temperature (Saha), by bisection in log ρ."
function rho_from_pg(pg, T)
    lo, hi = log(1e-30), log(1e2)
    for _ in 1:200
        mid = 0.5(lo + hi)
        Z0.gas_p(exp(mid), T) > pg ? (hi = mid) : (lo = mid)
        hi - lo < 1e-13 && break
    end
    exp(0.5(lo + hi))
end

"""
Grey Eddington atmosphere above the photosphere r = R (τ = 2/3, T = Tb), hydrostatic with the
radiation force of the full flux on the Rosseland opacity:
  dp_g/dr = −ρ (g − κ F/c),  dτ/dr = −κρ,  T⁴ = ¾ T_b⁴ (τ + ⅔).
Integrated outward until ρ has fallen by e^nH. Returns vectors (r, T, ρ, p_g).
"""
function atmosphere(R, Tb, ρb, M, L; nH = 12.0, tab = Z0.table())
    r = R; pg = Z0.gas_p(ρb, Tb); τ = 2 / 3
    rs = [r]; Ts = [Tb]; ρs = [ρb]; ps = [pg]
    ρstop = ρb * exp(-nH)
    T = Tb; ρ = ρb
    while ρ > ρstop
        g = Z0.G * M / r^2; F = L / (4π * r^2)
        κ = Z0.kap(tab, ρ, T)
        geff = g - κ * F / Z0.c
        geff > 0 || error("atmosphere: Γ ≥ 1 at r/R = $(r/R), T = $T (no static atmosphere)")
        H = pg / (ρ * geff)
        dr = 0.02H
        # Heun step in (ln p_g, τ)
        f(r, lp, τ) = begin
            Tl = (0.75 * Tb^4 * (max(τ, 0.0) + 2 / 3))^0.25
            ρl = rho_from_pg(exp(lp), Tl); κl = Z0.kap(tab, ρl, Tl)
            gl = Z0.G * M / r^2 - κl * L / (4π * r^2 * Z0.c)
            (-ρl * gl / exp(lp), -κl * ρl)
        end
        lp = log(pg)
        k1 = f(r, lp, τ); k2 = f(r + dr, lp + dr * k1[1], τ + dr * k1[2])
        lp += 0.5dr * (k1[1] + k2[1]); τ += 0.5dr * (k1[2] + k2[2]); r += dr
        pg = exp(lp); T = (0.75 * Tb^4 * (max(τ, 0.0) + 2 / 3))^0.25; ρ = rho_from_pg(pg, T)
        push!(rs, r); push!(Ts, T); push!(ρs, ρ); push!(ps, pg)
    end
    return (r = rs, T = Ts, ρ = ρs, pg = ps, τtop = τ)
end

"""
Free-fall inflow above the photosphere R: ρ = Ṁ/(4πr²v_ff), T from Lucy's dilute radiation field with the
inflow's own (LTE, Rosseland) optical depth, T⁴ = T_L⁴ [W(r) + ¾τ], T_L⁴ = L/(4πσR²) (stellar L only:
the accretion luminosity is produced by the code at the shock). Same construction as the steady model's
boundary (QuasiStarInfall/src/boundary.jl).
"""
function inflow_region(R, M, L, Mdot; rout = 2.0, n = 4000, tab = Z0.table())
    TL4 = L / (4π * Z0.σSB * R^2)
    rs = exp.(range(log(rout * R), log(R), length = n))
    τ = 0.0; Ts = Float64[]; ρs = Float64[]
    for (i, r) in enumerate(rs)
        ρ = Mdot / (4π * r^2 * sqrt(2Z0.G * M / r))
        W = 0.5 * (1 - sqrt(max(1 - (R / r)^2, 0.0)))
        T = (TL4 * (W + 0.75τ))^0.25
        push!(Ts, T); push!(ρs, ρ)
        i < n && (τ += Z0.kap(tab, ρ, T) * ρ * (r - rs[i+1]))
    end
    return (r = reverse(rs), T = reverse(Ts), ρ = reverse(ρs), τtop = τ)
end

rout_actual(atm, infall, rout, R) = infall ? rout * R : atm.r[end]

"""
n+1 interfaces on [a, b] with equal arc length s = ∫ √(1/xscale² + Σ (d y/dr)²) dr of the tabulated columns `ys`
(here ln ρ, ln T), then limited to a neighbour spacing ratio ≤ `qmax` (TITAN's grid smoothing does the same).
"""
function equidistribute(rt, ys, a, b, n; xscale, qmax = 1.3, nfine = 40000)
    rf = collect(range(a, b, length = nfine))
    yf = [[linterp(rt, y, r) for r in rf] for y in ys]
    s = zeros(nfine)
    for i in 2:nfine
        dr = rf[i] - rf[i-1]
        s[i] = s[i-1] + sqrt((dr / xscale)^2 + sum((y[i] - y[i-1])^2 for y in yf))
    end
    rI = [linterp(s, rf, x) for x in range(0, s[end], length = n + 1)]
    rI[1] = a; rI[end] = b
    for _ in 1:2000                      # limit the spacing ratio by relaxing towards neighbour means
        d = diff(rI); bad = false
        for k in 2:n
            if d[k] > qmax * d[k-1] || d[k-1] > qmax * d[k]
                rI[k] = 0.5(rI[k-1] + rI[k+1]); bad = true; d = diff(rI)
            end
        end
        bad || break
    end
    rI
end

"n+1 interfaces from a to b, first spacing ≈ d0, growing geometrically."
function geometric(a, b, d0, n)
    f(q) = d0 * (q^n - 1) / (q - 1) - (b - a)
    lo, hi = 1.0 + 1e-12, 2.0
    for _ in 1:200; mid = 0.5(lo + hi); f(mid) > 0 ? (hi = mid) : (lo = mid); end
    q = 0.5(lo + hi)
    r = [a]; d = d0
    for _ in 1:n; push!(r, r[end] + d); d *= q; end
    r[end] = b; r
end

"Piecewise-linear interpolation of y(x) for increasing x (clamped)."
function linterp(x::Vector{Float64}, y::Vector{Float64}, xq)
    xq <= x[1] && return y[1]; xq >= x[end] && return y[end]
    i = searchsortedlast(x, xq)
    w = (xq - x[i]) / (x[i+1] - x[i]); (1 - w) * y[i] + w * y[i+1]
end

# ------------------------------------------------------------------ the initial state
"""
    setup!(profile_csv; template, rin = 0.4, nenv = 140, natm = 60, nH = 12, alpha = 0.8)

Load `template` (a TITANST1 state providing all the parameters we don't touch: SUBC), then overwrite
geometry, boundary conditions, constants, grid and the initial structure. Returns a NamedTuple with
the profile and the derived numbers. The state is on the :julia backend; `TF.save_state` stores it.
"""
function setup!(profile_csv; template = joinpath(TF.ROOT, "reference_runs", "init", "SUBC.tstate.gz"),
                rin = 0.4, nenv = 140, natm = 60, nH = 12.0, alpha = 0.8, xmext_fac = 10.0, grey = false,
                Mdot = 0.0, rout = 2.0, dr0 = 2e-4, tmax_inflow = 3e9, Mdot_max = Mdot, lgrid = 2, lady = (2, 2, 2),
                wcontact = 0.0, equi = false)
    ncell = nenv + natm
    iseven(ncell) || error("ncell must be even")
    P = read_profile(profile_csv); h = P.hdr; C = P.cols
    M, L, R = h[:M], h[:L], h[:R]
    tab = Z0.table()
    TF.load_state!(template)
    s!(n, x) = haskey(TF.JSTATE, n) ? (a = TF.var(n); a[] = eltype(a)(x)) : nothing   # (ncell is not in COMMON)
    v(n) = TF.var(n)::Vector{Float64}

    # ---- constants: those of the Z0 physics / steady model
    s!(:cgrav, Z0.G); s!(:cc, Z0.c); s!(:csigr, Z0.σSB); s!(:car, Z0.arad)
    # ---- grid size
    ngrs = Int(TF.getscalar(:ngrs)); ngre = ngrs + ncell - 1
    ngre + 2 <= length(v(:rn)) || error("ncell too large for mgr")
    s!(:ncell, ncell); s!(:ngre, ngre)
    # ---- physics switches
    s!(:lgeom, 2); s!(:g, 0.0); s!(:lrad, 1); s!(:ltran, 1); s!(:lam, 0)
    s!(:lribc, 3); s!(:lrobc, 1); s!(:xipl, 0.0); s!(:ximl, 0.0); s!(:xipr, 0.0); s!(:ximr, 0.0)
    s!(:teff, 0.0); s!(:xlum, L)
    s!(:lgrid, lgrid)                                 # 2 Eulerian, 1 adaptive (Dorfi–Drury)
    infall = Mdot > 0
    s!(:leibc, 1); s!(:leobc, infall ? 2 : 1); s!(:llibc, 0); s!(:llobc, 0)
    for n in (:phil0, :phil1, :phil2, :phir0, :phir1, :phir2, :ul, :ur, :uextl, :uextr, :pextr); s!(n, 0.0); end
    s!(:leos, 2); s!(:lopac, 3)                       # placeholders; the material overrides both

    # ---- composite structure table vs r (inner envelope … photosphere … atmosphere)
    r_e = reverse(C[:r]); T_e = reverse(C[:T]); ρ_e = reverse(C[:rho]); P_e = reverse(C[:P])
    frad_e = reverse(C[:nabla] ./ C[:nabla_rad])     # F_rad / F_tot
    m_e = reverse(C[:m])
    atm = infall ? inflow_region(R, M, L, Mdot; rout, n = 4000) : atmosphere(R, T_e[end], ρ_e[end], M, L; nH, tab)
    r_in = rin * R
    r_in > r_e[1] || error("rin below the first profile point")
    # interfaces: envelope uniform in ln P_tot, atmosphere uniform in r
    lP = log.(P_e)                                    # decreasing with r
    lPq = range(linterp(r_e, lP, r_in), lP[end], length = nenv + 1)
    rlP = reverse(r_e); lPr = reverse(lP)             # ln P increasing → r decreasing
    r_env = [linterp(lPr, rlP, x) for x in lPq]; r_env[1] = r_in; r_env[end] = R
    r_atm = infall ? geometric(R, rout * R, dr0 * R, natm) : collect(range(R, atm.r[end], length = natm + 1))
    rI = vcat(r_env, r_atm[2:end])                    # ncell + 1 interfaces
    # combined (r, ln T, ln ρ) table; with wcontact > 0 the envelope–inflow contact is spread over [R, R(1+wcontact)]
    # (log-linear in ρ and T, linear in u): a smooth "settling layer" instead of a jump, so an adaptive grid can resolve it
    keep = infall && wcontact > 0 ? findall(>=(R * (1 + wcontact)), atm.r) : collect(2:length(atm.r))
    rt = vcat(r_e, atm.r[keep]); lTt = log.(vcat(T_e, atm.T[keep])); lρt = log.(vcat(ρ_e, atm.ρ[keep]))
    fradt = vcat(frad_e, ones(length(keep)))
    if equi            # equidistribute: equal arc length in (r/xscale, ln ρ, ln T), TITAN's resolution measure for lady = log
        rI = equidistribute(rt, (lρt, lTt), r_in, rout_actual(atm, infall, rout, R), ncell; xscale = R)
    end

    if lgrid == 1      # TITAN's grid equation sets the first two and the last two zones equal (no time smoothing there)
        rI[2] = 0.5(rI[1] + rI[3]); rI[end-1] = 0.5(rI[end] + rI[end-2])
    end
    rn, dn, tn, un = v(:rn), v(:dn), v(:tn), v(:un)
    ern, frn = v(:ern), v(:frn)
    mu = 2; xmup1 = 3.0
    for (j, k) in enumerate(ngrs:ngre+1); rn[k] = rI[j]; end
    for k in ngrs:ngre
        rc = ((rn[k]^3 + rn[k+1]^3) / 2)^(1 / 3)       # volume midpoint
        tn[k] = exp(linterp(rt, lTt, rc)); dn[k] = exp(linterp(rt, lρt, rc))
        ern[k] = Z0.arad * tn[k]^4
    end
    for k in ngrs:ngre+1
        un[k] = infall && rn[k] > R ? -sqrt(2Z0.G * M / rn[k]) * (wcontact > 0 ? clamp((rn[k] / R - 1) / wcontact, 0.0, 1.0) : 1.0) : 0.0
        frn[k] = L / (4π * rn[k]^2) * linterp(rt, fradt, rn[k])
    end
    frn[ngrs] = L / (4π * rn[ngrs]^2)                 # lribc = 3: all of L enters as radiation
    # ---- masses (per steradian). Interior point mass = m(r_in) from the steady model.
    m_in = linterp(r_e, m_e, r_in)
    dvol(k) = (rn[k+1]^3 - rn[k]^3) / xmup1
    xmn = v(:xmn)
    xmn[ngrs] = m_in / 4π
    for k in ngrs+1:ngre+1; xmn[k] = xmn[k-1] + dn[k-1] * dvol(k - 1); end
    shell = xmn[ngre+1] - xmn[ngrs]
    # ghost zones (leibc = leobc = 1: mirror)
    k = ngrs - 1
    rn[k] = 2rn[k+1] - rn[k+2]; dn[k] = dn[k+1]; tn[k] = tn[k+1]; un[k] = 0.0; ern[k] = ern[k+1]
    frn[k] = frn[k+1] * (rn[k+1] / rn[k])^mu
    xmn[k] = xmn[k+1] - dn[k] * dvol(k)
    k = ngre + 1
    if infall      # leobc = 2: fixed external inflow state (ρ, T, u) in the ghost zones
        rext = rn[k]; vff = sqrt(2Z0.G * M / rext)
        ρx = Mdot / (4π * rext^2 * vff); Tx = exp(linterp(rt, lTt, rext))
        dn[k] = ρx; dn[k+1] = ρx; tn[k] = Tx; tn[k+1] = Tx
        rn[k+1] = 2rn[k] - rn[k-1]; un[k] = -vff; un[k+1] = -vff
        s!(:ur, -vff); s!(:dr, ρx); s!(:tr, Tx)
        s!(:phir0, -ρx * vff)                                   # mass flux through r_out (outward positive)
        s!(:phir2, -ρx * vff * Z0.gas_e(ρx, Tx))               # gas internal-energy flux
    else
        dn[k] = dn[k-1]; dn[k+1] = dn[k]; tn[k] = tn[k-1]; tn[k+1] = tn[k]
        rn[k+1] = 2rn[k] - rn[k-1]; un[k+1] = 0.0
    end
    ern[k] = ern[k-1]; ern[k+1] = ern[k]; frn[k+1] = frn[k]
    xmn[k+1] = xmn[k] + dn[k] * dvol(k)
    top = dn[ngre] * dvol(ngre)
    # exterior mass: must stay positive while the inflow removes Ṁ/4π per second from it
    xmext = infall ? 10 * max(Mdot, Mdot_max) / 4π * tmax_inflow : xmext_fac * top
    xmass = xmn[ngre+2]; xmtot = xmass + xmext
    s!(:xmass, xmass); s!(:xmext, xmext); s!(:xmtot, xmtot)
    if lgrid == 1   # resolution: linear in r (scale R) plus logarithmic in m (exterior), ρ, T (lady(1:3))
        s!(:ladx, 1); s!(:xscale, R); s!(:alph, 1.5); s!(:ibet, 1); s!(:tau, parse(Float64, get(ENV, "QS_GRIDTAU", "1e3")))
        la = TF.var(:lady); wt = TF.var(:wt); ys = TF.var(:yscl)
        la .= 0; wt .= 0.0
        for (i, l) in enumerate(lady); la[i] = l; wt[i] = l > 0 ? 1.0 : 0.0; end
        ys[1] = xmtot
    end
    xmen = v(:xmen)
    for k in ngrs-1:ngre+2; xmen[k] = xmtot - xmn[k]; end
    s!(:delml, dn[ngrs] * dvol(ngrs)); s!(:delmr, dn[ngre] * dvol(ngre))
    s!(:delrl, rn[ngrs+1] - rn[ngrs]); s!(:delrr, rn[ngre+2] - rn[ngre+1])
    infall || (s!(:dr, dn[ngre]); s!(:tr, tn[ngre]))
    for k in ngrs-1:ngre+2
        v(:rmun)[k] = rn[k]^2; v(:rmup1n)[k] = rn[k]^3; v(:rmum1n)[k] = rn[k]
    end
    for k in ngrs-1:ngre+1; v(:dvoln)[k] = dvol(k); end

    # ---- zero the bookkeeping arrays, set old = new, fill EOS/opacity/Eddington factors
    rng = ngrs-1:ngre+2
    for n in (:dsn, :egsn, :etn, :etsn, :ersn, :egtn, :egtsn, :egrsn, :egrtn, :egrtsn, :usn, :urel, :frsn,
              :dmdt, :drdt, :qx, :xnto, :xnt, :tescr, :tescr0)
        view(v(n), rng) .= 0.0
    end
    for n in (:tee, :tes, :tew, :tel, :teq, :tetot, :teso, :tewo, :telo, :teqo); s!(n, 0.0); end
    s!(:geddl, 0.0); s!(:geddr, 0.0)
    for (o, n) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(n), rng)); end
    copyto!(view(v(:ero), rng), view(ern, rng)); copyto!(view(v(:fro), rng), view(frn, rng))
    for k in rng; v(:frnom)[k] = abs(frn[k]) + 1e-30; v(:unom)[k] = abs(un[k]) + 1e-30; end

    TX.CONVECTION[] = true; TX.ALPHA_MLT[] = alpha
    TM.use!(TM.Z0Gas(tab; grey)); TC.EDDFAC[] = TC.eddfac!
    TD.EOS[](); TD.OPAC[]()
    for (o, n) in TC.PEQ_OLD2_RH; copyto!(view(v(o), ngrs-1:ngre+1), view(v(n), ngrs-1:ngre+1)); end
    s!(:dtime, 1.0)
    TC._ghost_and_derived!(:peq, () -> TD.EOS[](), () -> TD.OPAC[]())
    for (o, n) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(n), rng)); end
    for (o, n) in TC.PEQ_OLD2_RH; copyto!(view(v(o), ngrs-1:ngre+1), view(v(n), ngrs-1:ngre+1)); end
    copyto!(view(v(:avchio), ngrs:ngre+1), view(v(:avchin), ngrs:ngre+1))
    for k in ngrs:ngre+1
        v(:xm0)[k] = xmn[k]; v(:r0)[k] = rn[k]; v(:u0)[k] = un[k]
    end
    s!(:timen, 0.0); s!(:timeo, 0.0); s!(:tfac, 1.0); s!(:dtime, 1e3)
    s!(:cmax, 0.0); s!(:dmax, 0.0); s!(:smax, 0.0)
    return (; M, L, R, r_in, m_in, shell = 4π * shell, atm, ngrs, ngre, rI, prof = P)
end

"Hydrostatic/thermal diagnostics of the current state (zone-centred and interface quantities)."
function diagnostics(; L = TF.getscalar(:xlum))
    ngrs, ngre = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    v(n) = TF.var(n)::Vector{Float64}
    rn, dn, tn, un, frn, ern = v(:rn), v(:dn), v(:tn), v(:un), v(:frn), v(:ern)
    (; r = rn[ngrs:ngre+1], d = dn[ngrs:ngre], t = tn[ngrs:ngre], u = un[ngrs:ngre+1],
       Lrad = [4π * rn[k]^2 * frn[k] for k in ngrs:ngre+1])
end

end # module
