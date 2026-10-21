# S6a, approach B (production): Lagrangian accreting protostar on a self-similar prescribed mass grid (TitanResidFast.ALE_MDOT) with SPS86/OP03-type surface
# boundary conditions. The star is TITAN's lgrid = 3 problem with a free outer boundary (as S5); mass is added at the surface (phir0 = −Ṁ/4πR²) and the mass grid is
# m_k ∝ M(t): the matter flows through the zones, composition is advected (StarMaterial.advect!).
# Surface BC: pressure P_ext = ram pressure Ṁ v_ff(R)/(4πR²) of the post-shock gas; incident radiation = radiation field E_s of the optically thick free-fall precursor
# (Lucy: T_s⁴ = T_eff⁴(½ + ¾τ_prec), L_tot = L* + G M Ṁ/R); the accreted gas enters with the internal energy of the top zone ("cold accretion").
# Included by s6b_run.jl (after s4_common.jl and s6_star.jl, s6_setup.jl for accreted_Y).

"""
    precursor(R, M, Mdot, Ltot; X, rmax = 50R, n = 300) -> (; Ts, tau, Tph, rph, Es)

Optically thick free-fall precursor of the accretion shock at radius R: ρ(r) = Ṁ/(4πr²v_ff), Lucy's radiation temperature T⁴ = T_L⁴ (W(r) + ¾τ(r)),
T_L⁴ = L_tot/(4πσR²), W = ½(1 − √(1 − R²/r²)), τ(r) the Rosseland optical depth from infinity (evaluated with the gas density and T(r), floor 2.5e3 K).
Ts = T(R) (the radiation temperature at the shock), tau = τ(R), the photosphere r_ph where τ = 2/3 with its T_eff,ph = (L_tot/4πr_ph²σ)^¼ (dilution included).
"""
function precursor(R, M, Mdot, Ltot; X = 0.75, rmax = 50.0, n = 300)
    r = exp.(range(log(rmax * R), log(R), length = n))        # decreasing
    TL4 = Ltot / (4π * Z0Physics.σSB * R^2)
    τ = 0.0; T = 0.0; rph = NaN; Tph = NaN
    for i in 1:n
        ρ = Mdot / (4π * r[i]^2 * sqrt(2GRAV * M / r[i]))
        W = 0.5 * (1 - sqrt(max(1 - (R / r[i])^2, 0.0)))
        T = max((TL4 * (W + 0.75τ))^0.25, 2.5e3)
        if i < n
            dτ = StellarOpacity.kappa_R(ρ, T, X) * ρ * (r[i] - r[i+1])
            if isnan(rph) && τ + dτ >= 2 / 3
                w = (2 / 3 - τ) / dτ; rph = r[i] - w * (r[i] - r[i+1])
                Tph = (Ltot / (4π * rph^2 * Z0Physics.σSB))^0.25
            end
            τ += dτ
        end
    end
    isnan(rph) && (rph = R; Tph = (Ltot / (4π * R^2 * Z0Physics.σSB))^0.25)       # optically thin: the photosphere is the shock
    (; Ts = T, tau = τ, Tph, rph, Es = Z0Physics.arad * T^4)
end

"""
    build_ale_star(; rhoc, Mdot, Tc, ...) -> (; mat, M, R, ...)

Lagrangian star (isentrope of central temperature Tc and density rhoc, ending where P = P_ram), nin interior zones geometric in mass + nenv envelope zones log-spaced
in the exterior mass fraction (qsplit → mfloor), regular centre, free outer boundary at pressure P_ram, StarGas (accreted composition), MLT.
"""
function build_ale_star(; rhoc, Mdot, Tc, X = 0.75, DH = 2.5e-5, He3H = 1e-5, nin = 70, nenv = 50, qc = 1.12, qsplit = 0.05, mfloor = 3e-7, Lstar = 3 * Lsun,
                        alpha = 1.5, nuc = true, efrac = 0.5, dXmax = 4e-3, burnT = :path, mixing = true, fcoupling = true, kcap = Inf, xmext_fac = 1e-3, opac = :stellar)
    n = nin + nenv; iseven(n) || error("zone number must be even")
    s = isentropic_star(0.6 * Msun, Tc, rhoc; X, Mdot, dlm = 0.02)
    (s.stopped && !s.bad) || error("isentrope did not reach P_ram (Mstop = $(s.Mstop / Msun) Msun)")
    R = s.R; M = s.Mstop
    @printf("isentropic star: M = %.4f Msun, R = %.3f Rsun, rho_c = %.4g, T_s = %.3e, rho_s = %.3e\n", M / Msun, R / Rsun, s.ρ[1], s.T[end], s.ρ[end])
    w = [qc^(j - 1) for j in 1:nin]; mint = (1 - qsplit) .* [0; cumsum(w)] ./ sum(w)
    qe = 10 .^ range(log10(qsplit), log10(mfloor); length = nenv + 1)[2:end]
    mfrac = vcat(mint, 1 .- qe); mfrac[end] = 1.0
    mI = M .* mfrac
    lin(x, y, xq) = (xq <= x[1] ? y[1] : xq >= x[end] ? y[end] : (i = searchsortedlast(x, xq); w_ = (xq - x[i]) / (x[i+1] - x[i]); (1 - w_) * y[i] + w_ * y[i+1]))
    rI = [cbrt(lin(s.m, s.r .^ 3, mm)) for mm in mI]; rI[1] = 0.0; rI[end] = R
    dvol(a, b) = (b^3 - a^3) / 3
    ρz = [(mI[j+1] - mI[j]) / (4π * dvol(rI[j], rI[j+1])) for j in 1:n]
    Tz = [lin(s.m, s.T, 0.5 * (mI[j] + mI[j+1])) for j in 1:n]
    Er = Z0Physics.arad .* Tz .^ 4
    Fr = [j == 1 ? 0.0 : Lstar / (4π * rI[j]^2) for j in 1:n+1]
    mat = SM.StarGas(TF.MGR; X, Y = 1 - X, Z = 0.0, nuc, opac, efrac, dXmax, burnT, screen = :chugunov, mixing, kcap, fcoupling)
    Y0 = accreted_Y(X, DH, He3H)
    for k in 1:size(mat.Y, 2); SM.setY!(mat, k, Y0); end
    Pram = Mdot * sqrt(2GRAV * M / R) / (4π * R^2)
    PS.build!(; rI, rho = ρz, T = Tz, u = zeros(n + 1), geom = 2, rad = true, Er, Fr, lgrid = 3, inner = :wall, outer = :free, pext = Pram, selfgrav = true, material = mat,
              convection = true, gam = 5 / 3, gmu = 1.0, xmext = xmext_fac * M / (4π))
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = Float64[]
    TX.ALPHA_MLT[] = alpha
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    (; mat, M, R, rI, ρz, Tz, n, Y0, mI, Lstar)
end
