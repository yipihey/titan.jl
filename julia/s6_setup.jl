# S6a (approach A): TITAN state of an accreting protostar: a convective (isentropic) star with a regular centre, the accretion shock at its surface and a
# free-fall inflow region out to rout·R*, on TITAN's adaptive grid (lgrid = 1), inflow outer boundary (leobc = 2), StarGas (burning, mixing).
# Included by s6_run.jl / s6_probe*.jl (after s4_common.jl and s6_star.jl).
const YR = 3.15576e7

"Abundance vector (mol/g) of the accreted gas: hydrogen mass fraction X (p + D), D/H and ³He/H by number, the rest ⁴He."
function accreted_Y(X = 0.75, DH = 2.5e-5, He3H = 1e-5)
    AA = NN.AA; Y = zeros(NN.NS)
    x = X / (AA[SM.IH] + AA[SM.ID] * DH)
    Y[SM.IH] = x; Y[SM.ID] = DH * x; Y[NN.IDX[:he3]] = He3H * x
    Y[NN.IDX[:he4]] = (1 - sum(AA[i] * Y[i] for i in 1:NN.NS)) / AA[NN.IDX[:he4]]
    Y
end

"Free-fall inflow above the shock at radius R for enclosed mass M and Ṁ: ρ(r), u(r) = −v_ff and the Lucy radiation temperature T(r) of the luminosity Ltot (diluted + optically thick part)."
function inflow_profile(r, R, M, Mdot, Ltot, mat; X = 0.75)
    n = length(r)
    ρ = [Mdot / (4π * ri^2 * sqrt(2GRAV * M / ri)) for ri in r]
    TL4 = Ltot / (4π * Z0Physics.σSB * R^2)
    T = zeros(n); τ = 0.0
    for i in n:-1:1      # r increasing; integrate τ from the outside inwards
        W = 0.5 * (1 - sqrt(max(1 - (R / r[i])^2, 0.0)))
        T[i] = max((TL4 * (W + 0.75τ))^0.25, 2.5e3)
        i > 1 && (τ += StellarOpacity.kappa_R(ρ[i], T[i], X) * ρ[i] * (r[i] - r[i-1]))
    end
    (; ρ, T, u = -[sqrt(2GRAV * M / ri) for ri in r])
end

"Geometric grid from a to b with n zones, first width d0."
function geom_grid(a, b, d0, n)
    f(q) = d0 * (q^n - 1) / (q - 1) - (b - a)
    lo, hi = 1.0 + 1e-12, 3.0
    for _ in 1:200; mid = 0.5(lo + hi); f(mid) > 0 ? (hi = mid) : (lo = mid); end
    q = 0.5(lo + hi); r = [a]; d = d0
    for _ in 1:n; push!(r, r[end] + d); d *= q; end
    r[end] = b; r
end

"""
    build_accreting(; rhoc, Mdot, Tc, ...) -> (; mat, ...)

Star: isentrope of central temperature Tc and mass M, ending where P = P_ram (contact); nin zones geometric in mass (ratio qc) + nenv zones log-spaced in the
exterior mass fraction (qsplit → mfloor); inflow region of ninf zones from R to rout·R, first width dr0·R. Initial luminosity Lstar everywhere in the star, Lacc = G M Ṁ/R
added in the inflow. TITAN: lgrid = 1, ladx = 1, xscale = R, regular centre, leibc = 1, leobc = 2, lrobc = 1 (no incident radiation).
"""
function build_accreting(; rhoc, Mdot, Tc, X = 0.75, DH = 2.5e-5, He3H = 1e-5, nin = 60, nenv = 40, ninf = 70, rout = 4.0, dr0 = 1e-3, qc = 1.15, qsplit = 0.05, mfloor = 1e-9,
                         Lstar = 3 * Lsun, alpha = 1.5, nuc = true, gridtau = 1e7, lady = (2, 2, 2), gridw = (1.0, 1.0, 1.0), xmext_fac = 100.0,
                         efrac = 0.5, dXmax = 4e-3, burnT = :path, xscale_fac = 1.0, mixing = true, fcoupling = true, kcap = Inf)
    n = nin + nenv + ninf; iseven(n) || error("zone number must be even")
    # ---- star
    s = isentropic_star(0.5 * Msun, Tc, rhoc; X, Mdot, dlm = 0.02)
    (s.stopped && !s.bad) || error("isentrope did not reach P_ram (Mstop = $(s.Mstop / Msun) Msun)")
    R = s.R; Mst = s.Mstop
    @printf("isentropic star: M = %.4f Msun, R = %.3f Rsun, rho_c = %.4g, T_s = %.3e, rho_s = %.3e\n", Mst / Msun, R / Rsun, s.ρ[1], s.T[end], s.ρ[end])
    w = [qc^(j - 1) for j in 1:nin]; mint = (1 - qsplit) .* [0; cumsum(w)] ./ sum(w)
    qe = 10 .^ range(log10(qsplit), log10(mfloor); length = nenv + 1)[2:end]
    mfrac = vcat(mint, 1 .- qe); mfrac[end] = 1.0
    mI = Mst .* mfrac                                    # interface masses of the star (nin + nenv + 1)
    lin(x, y, xq) = (xq <= x[1] ? y[1] : xq >= x[end] ? y[end] : (i = searchsortedlast(x, xq); w_ = (xq - x[i]) / (x[i+1] - x[i]); (1 - w_) * y[i] + w_ * y[i+1]))
    rI = [cbrt(lin(s.m, s.r .^ 3, mm)) for mm in mI]; rI[1] = 0.0; rI[end] = R       # (r³ is linear in m near the centre)
    # TITAN's grid equation makes the first two zones equally wide: r₂ = r₃/2, and the mass of zone 1 follows from the isentrope (r³ linear in m there)
    rI[2] = 0.5 * rI[3]; mI[2] = lin(s.r .^ 3, s.m, rI[2]^3)
    nst = nin + nenv
    dvol(a, b) = (b^3 - a^3) / 3
    ρz = [(mI[j+1] - mI[j]) / (4π * dvol(rI[j], rI[j+1])) for j in 1:nst]
    Tz = [lin(s.m, s.T, 0.5 * (mI[j] + mI[j+1])) for j in 1:nst]
    # ---- inflow
    rinf = geom_grid(R, rout * R, dr0 * R, ninf)
    rc = [((rinf[j]^3 + rinf[j+1]^3) / 2)^(1 / 3) for j in 1:ninf]
    Lacc = GRAV * Mst * Mdot / R
    inf = inflow_profile(rc, R, Mst, Mdot, Lacc + Lstar, nothing; X)
    rI = vcat(rI, rinf[2:end])
    # TITAN's grid equation forces the last two zones to equal width as well
    rI[end-1] = 0.5(rI[end] + rI[end-2])
    ρ = vcat(ρz, inf.ρ); T = vcat(Tz, inf.T)
    # recompute the two zones whose widths changed so that the mass of every zone stays (the isentrope's) one
    for j in (1, 2, n - 1, n)
        j <= nst && (ρ[j] = (mI[j+1] - mI[j]) / (4π * dvol(rI[j], rI[j+1])))
    end
    u = zeros(n + 1)
    for j in nst+2:n+1; u[j] = -sqrt(2GRAV * Mst / rI[j]); end
    u[nst+1] = 0.0
    Er = Z0Physics.arad .* T .^ 4
    Fr = [j == 1 ? 0.0 : (j <= nst + 1 ? Lstar : Lacc + Lstar) / (4π * rI[j]^2) for j in 1:n+1]
    Fr[nst+1] = (Lstar) / (4π * rI[nst+1]^2)
    mat = SM.StarGas(TF.MGR; X, Y = 1 - X, Z = 0.0, nuc, opac = :stellar, efrac, dXmax, burnT, screen = :chugunov, mixing, kcap, fcoupling)
    Y0 = accreted_Y(X, DH, He3H)
    for k in 1:size(mat.Y, 2); SM.setY!(mat, k, Y0); end
    ρx, Tx, ux = inf.ρ[end], inf.T[end], inf.u[end]
    xmext = xmext_fac * Mst / (4π)
    PS.build!(; rI, rho = ρ, T, u, geom = 2, rad = true, Er, Fr, lgrid = 1, inner = :wall, outer = :inflow, inflow = (ρx, Tx, ux), material = mat, gam = 5 / 3, gmu = 1.0,
              gridtau, lady, xscale = xscale_fac * R, xmext, convection = true, selfgrav = true)
    # weights of the grid terms (PS.build! sets 1 for every active term)
    wt = TF.var(:wt); la = TF.var(:lady); for (i, l) in enumerate(lady); wt[i] = l > 0 ? gridw[i] : 0.0; end
    TX.REGULAR_CENTRE[] = true
    TX.ALPHA_MLT[] = alpha
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    # Z0 constants are set by PS.build! (cgrav = Z0.G, ...). Re-fill derived state with the StarGas EOS
    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    (; mat, M = Mst, R, rI, ρ, T, nst, n, Lacc, Y0, mI)
end
