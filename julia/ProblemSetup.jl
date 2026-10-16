"""
    ProblemSetup

Generic construction of a TITAN state from profiles, for test problems written in Julia (no Fortran `start*.f`).
A captured benchmark state supplies every parameter that is not overwritten: SOD for pure hydrodynamics (neqn = 5),
SUBC for radiation hydrodynamics (neqn = 7).

    build!(; rI, rho, T, u, geom = 2, rad = false, Er = nothing, Fr = nothing,
             inner = :wall, outer = :wall, inflow = nothing, material = TM.Builtin(), gam = 5/3, gmu = 0.5,
             Mc = 0.0, lgrid = 2, gridtau = 1e3, lady = (0, 2, 2), xscale = rI[end], av = (;), xmext = nothing)

- `rI` interfaces (n+1), `rho`, `T` zone values (n), `u` interface velocities (n+1); `Er` zones, `Fr` interfaces.
- `geom` 0 planar, 2 spherical. `Mc` central point mass (g; spherical gravity on with the Z0 value of G).
- `inner`/`outer` ∈ (:wall, :inflow): Eulerian zero-flux wall (le?bc = 1) or a fixed external state (le?bc = 2);
  `inflow = (rho, T, u)` gives that state for the outer boundary (u < 0 for inflow). Mass and gas-energy fluxes through
  the boundary follow from it: phir0 = ρu, phir2 = ρ u e_g(ρ, T).
- `av` overrides artificial-viscosity parameters, e.g. `(ql0 = …, ql1 = …, cq1 = …, cq2 = …)`.
Radiation boundaries (rad = true): inner reflecting (lribc = 2) unless `lum` is given (lribc = 3, imposed L),
outer transmitting (lrobc = 1).
"""
module ProblemSetup

import ..TitanDriver, ..TitanCore, ..TitanResidFast, ..TitanMaterials, ..Z0Physics
const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
const TM = TitanMaterials; const Z0 = Z0Physics

const INIT = joinpath(@__DIR__, "..", "reference_runs", "init")

function build!(; rI, rho, T, u, geom = 2, rad = false, Er = nothing, Fr = nothing, lum = nothing,
                inner = :wall, outer = :wall, inflow = nothing, material = TM.Builtin(), gam = 5 / 3, gmu = 0.5,
                Mc = 0.0, lgrid = 2, gridtau = 1e3, lady = (0, 2, 2), xscale = rI[end] - rI[1], av = (;), xmext = nothing,
                convection = false, opac = nothing)
    n = length(rho); length(rI) == n + 1 && length(T) == n && length(u) == n + 1 || error("profile lengths")
    iseven(n) || error("ncell must be even")
    TF.load_state!(joinpath(INIT, rad ? "SUBC.tstate.gz" : "SOD.tstate.gz"))
    s!(k, x) = haskey(TF.JSTATE, k) ? (a = TF.var(k); a[] = eltype(a)(x)) : nothing
    v(k) = TF.var(k)::Vector{Float64}
    ngrs = Int(TF.getscalar(:ngrs)); ngre = ngrs + n - 1
    ngre + 2 <= length(v(:rn)) || error("too many zones for mgr")
    s!(:ngre, ngre)
    mu = geom; s!(:lgeom, geom); s!(:mu, mu); s!(:mup1, mu + 1); s!(:mum1, mu - 1)
    s!(:xmu, mu); s!(:xmup1, mu + 1); s!(:xmum1, mu - 1); s!(:rxm1, 1 / (mu + 1))
    xmup1 = mu + 1.0
    s!(:cgrav, Mc > 0 ? Z0.G : 0.0); s!(:g, 0.0)
    rad && (s!(:cc, Z0.c); s!(:csigr, Z0.σSB); s!(:car, Z0.arad))
    s!(:leos, 2); s!(:gam, gam); s!(:gmu, gmu)
    opac === nothing || for (k, x) in pairs(opac); s!(k, x); end      # e.g. (lopac = 3, chif0 = …)
    s!(:lgrid, lgrid)
    s!(:leibc, inner === :wall ? 1 : 2); s!(:leobc, outer === :wall ? 1 : 2); s!(:llibc, 0); s!(:llobc, 0)
    for k in (:phil0, :phil1, :phil2, :phir0, :phir1, :phir2, :ul, :ur, :uextl, :uextr, :pextr); s!(k, 0.0); end
    for (k, x) in pairs(av); s!(k, x); end
    if rad
        s!(:lrad, 1); s!(:ltran, 1); s!(:lrobc, 1); s!(:ximr, 0.0); s!(:xipl, 0.0); s!(:teff, 0.0)
        lum === nothing ? (s!(:lribc, 2); s!(:xlum, 0.0)) : (s!(:lribc, 3); s!(:xlum, lum))
    end

    rn, dn, tn, un = v(:rn), v(:dn), v(:tn), v(:un)
    for (j, k) in enumerate(ngrs:ngre+1); rn[k] = rI[j]; un[k] = u[j]; end
    for (j, k) in enumerate(ngrs:ngre); dn[k] = rho[j]; tn[k] = T[j]; end
    if rad
        ern, frn = v(:ern), v(:frn)
        for (j, k) in enumerate(ngrs:ngre); ern[k] = Er === nothing ? Z0.arad * T[j]^4 : Er[j]; end
        for (j, k) in enumerate(ngrs:ngre+1); frn[k] = Fr === nothing ? 0.0 : Fr[j]; end
        lum === nothing || (frn[ngrs] = lum / (4π * rn[ngrs]^2))
    end
    dvol(k) = (rn[k+1]^(mu + 1) - rn[k]^(mu + 1)) / xmup1
    # ghost zones
    k = ngrs - 1
    rn[k] = 2rn[k+1] - rn[k+2]; dn[k] = dn[k+1]; tn[k] = tn[k+1]; un[k] = inner === :wall ? 0.0 : un[k+1]
    k = ngre + 1
    rn[k+1] = 2rn[k] - rn[k-1]
    if outer === :inflow
        inflow === nothing && error("outer = :inflow needs inflow = (rho, T, u)")
        ρx, Tx, ux = inflow
        dn[k] = ρx; dn[k+1] = ρx; tn[k] = Tx; tn[k+1] = Tx; un[k] = ux; un[k+1] = ux
        # TITAN's ideal gas (JuliaPhysics.eos!): e = (k/(m0 μ)) T/(γ−1); Z0 gas: translational + ionization
        ex = material isa TM.Builtin ? TF.getscalar(:ck) / (TF.getscalar(:cm0) * gmu) * Tx / (gam - 1) : Z0.gas_e(ρx, Tx)
        s!(:ur, ux); s!(:dr, ρx); s!(:tr, Tx); s!(:phir0, ρx * ux); s!(:phir2, ρx * ux * ex)
    else
        dn[k] = dn[k-1]; dn[k+1] = dn[k]; tn[k] = tn[k-1]; tn[k+1] = tn[k]; un[k+1] = 0.0
        s!(:dr, dn[ngre]); s!(:tr, tn[ngre])
    end
    if rad
        ern, frn = v(:ern), v(:frn)
        ern[ngrs-1] = ern[ngrs]; frn[ngrs-1] = -frn[ngrs+1]
        ern[ngre+1] = ern[ngre]; ern[ngre+2] = ern[ngre]; frn[ngre+2] = frn[ngre+1]
    end
    # masses (per steradian in spherical geometry): xmn = Mc/4π + mass inside
    xmn = v(:xmn)
    xmn[ngrs] = mu == 2 ? Mc / 4π : 0.0
    for k in ngrs+1:ngre+1; xmn[k] = xmn[k-1] + dn[k-1] * dvol(k - 1); end
    xmn[ngrs-1] = xmn[ngrs] - dn[ngrs-1] * dvol(ngrs - 1); xmn[ngre+2] = xmn[ngre+1] + dn[ngre+1] * dvol(ngre + 1)
    top = dn[ngre] * dvol(ngre)
    xme = xmext === nothing ? 1e3 * (xmn[ngre+2] - xmn[ngrs]) : xmext
    xmtot = xmn[ngre+2] + xme
    s!(:xmass, xmn[ngre+2]); s!(:xmext, xme); s!(:xmtot, xmtot)
    xmen = v(:xmen); for k in ngrs-1:ngre+2; xmen[k] = xmtot - xmn[k]; end
    s!(:delml, dn[ngrs] * dvol(ngrs)); s!(:delmr, dn[ngre] * dvol(ngre))
    s!(:delrl, rn[ngrs+1] - rn[ngrs]); s!(:delrr, rn[ngre+2] - rn[ngre+1])
    for k in ngrs-1:ngre+2; v(:rmun)[k] = rn[k]^mu; v(:rmup1n)[k] = rn[k]^(mu + 1); v(:rmum1n)[k] = rn[k]^(mu - 1); end
    for k in ngrs-1:ngre+1; v(:dvoln)[k] = dvol(k); end
    if lgrid == 1
        s!(:ladx, 1); s!(:xscale, xscale); s!(:alph, 1.5); s!(:ibet, 1); s!(:tau, gridtau)
        la = TF.var(:lady); wt = TF.var(:wt); ys = TF.var(:yscl); la .= 0; wt .= 0.0
        for (i, l) in enumerate(lady); la[i] = l; wt[i] = l > 0 ? 1.0 : 0.0; end
        ys[1] = xmtot
    end
    # bookkeeping, old = new, derived quantities
    rng = ngrs-1:ngre+2
    for k in (:dsn, :egsn, :etn, :etsn, :ersn, :egtn, :egtsn, :egrsn, :egrtn, :egrtsn, :usn, :urel, :frsn,
              :dmdt, :drdt, :qx, :xnto, :xnt, :tescr, :tescr0)
        haskey(TF.JSTATE, k) && (view(v(k), rng) .= 0.0)
    end
    for k in (:tee, :tes, :tew, :tel, :teq, :tetot, :teso, :tewo, :telo, :teqo); s!(k, 0.0); end
    s!(:geddl, 0.0); s!(:geddr, 0.0)
    TX.CONVECTION[] = convection
    TM.use!(material); TC.EDDFAC[] = TC.eddfac!
    for (o, nn) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(nn), rng)); end
    if rad
        copyto!(view(v(:ero), rng), view(v(:ern), rng)); copyto!(view(v(:fro), rng), view(v(:frn), rng))
        for k in rng; v(:frnom)[k] = abs(v(:frn)[k]) + 1e-30; end
    end
    for k in rng; v(:unom)[k] = rad ? abs(un[k]) + 1e-30 : max(abs(un[k]), 1.0); end
    TD.EOS[](); rad && TD.OPAC[]()
    old2 = rad ? TC.PEQ_OLD2_RH : TC.PEQ_OLD2_H
    for (o, nn) in old2; copyto!(view(v(o), ngrs-1:ngre+1), view(v(nn), ngrs-1:ngre+1)); end
    s!(:dtime, 1.0)
    TC._ghost_and_derived!(:peq, () -> TD.EOS[](), () -> TD.OPAC[]())
    for (o, nn) in TC.PEQ_OLD1; copyto!(view(v(o), rng), view(v(nn), rng)); end
    for (o, nn) in old2; copyto!(view(v(o), ngrs-1:ngre+1), view(v(nn), ngrs-1:ngre+1)); end
    rad && copyto!(view(v(:avchio), ngrs:ngre+1), view(v(:avchin), ngrs:ngre+1))
    for k in ngrs:ngre+1; v(:xm0)[k] = xmn[k]; v(:r0)[k] = rn[k]; v(:u0)[k] = un[k]; end
    s!(:timen, 0.0); s!(:timeo, 0.0); s!(:tfac, 1.0); s!(:cmax, 0.0); s!(:dmax, 0.0); s!(:smax, 0.0)
    s!(:jstep, 1); s!(:jstepe, 10^7)
    nothing
end

"Install every routine of the time step in Julia (as phase6_run.jl)."
function install_julia!(BS)
    TD.ROUTINES[:matgen] = () -> TX.matgen!()
    TD.ROUTINES[:matslv] = BS.Penta.matslv!
    for (s, f) in ((:peq, TC.peq!), (:reset, TC.reset!), (:update, TC.update!), (:maxdel, TC.maxdel!),
                   (:timstp, TC.timstp!), (:totnrg, TC.totnrg!))
        TD.ROUTINES[s] = f
    end
    TC.EDDFAC[] = TC.eddfac!
    nothing
end

"Snapshot of the profile: interfaces r, u; zones ρ, T, p."
function snap()
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    v(k) = copy(TF.var(k)::Vector{Float64})
    (; t = TF.getscalar(:timen), dt = TF.getscalar(:dtime), r = v(:rn)[a:b+1], u = v(:un)[a:b+1],
       rho = v(:dn)[a:b], T = v(:tn)[a:b], p = v(:pgn)[a:b])
end

end # module
