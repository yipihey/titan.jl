# S4: Eddington standard model (n = 3 polytrope) of a radiation-gas star as a TITAN rad-hydro initial state
# (spherical, regular centre, Lagrangian mass grid, free surface, StarGas, MLT). Included by s4_a.jl, s4_c.jl, s4_d.jl.
include(joinpath(@__DIR__, "s4_common.jl"))

"Lane–Emden n = 3: arrays ξ, θ, μ(ξ) = −ξ²θ' (RK4 to the first zero ξ1)."
function lane_emden3(; h = 2e-4)
    ξ = [1e-8]; θ = [1.0]; μ = [0.0]
    f(x, th, ph) = (ph, -th^3 - 2ph / x)       # θ'' = −θ^n − 2θ'/ξ; state (θ, θ')
    y = (1.0 - 1e-16, 0.0); x = 1e-8
    while y[1] > 0
        k1 = f(x, y...); k2 = f(x + h / 2, y[1] + h / 2 * k1[1], y[2] + h / 2 * k1[2])
        k3 = f(x + h / 2, y[1] + h / 2 * k2[1], y[2] + h / 2 * k2[2]); k4 = f(x + h, y[1] + h * k3[1], y[2] + h * k3[2])
        yn = (y[1] + h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]), y[2] + h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]))
        x += h
        if yn[1] <= 0                       # linear interpolation to the zero
            xz = x - h + h * y[1] / (y[1] - yn[1]); push!(ξ, xz); push!(θ, 0.0); push!(μ, -xz^2 * (y[2] + (yn[2] - y[2]) * y[1] / (y[1] - yn[1]))); break
        end
        y = yn; push!(ξ, x); push!(θ, y[1]); push!(μ, -x^2 * y[2])
    end
    ξ, θ, μ
end
"Solve the Eddington quartic (1−β)/β⁴ = 3 (k/(μ m_u))⁴ /(a (πG)^{... }) in the form given M: returns β."
function eddington_beta(M, μmol; mu1 = 2.0182359519173048)
    # M = 4π (K/(πG))^{3/2} μ1 with K³ = (k/(μ m_u))⁴ 3(1−β)/(a β⁴)  ⇒ (1−β)/β⁴ = a (πG)³ (M/(4π μ1))² (μ m_u/k)⁴/3
    G = Z0.G; k = Z0.kB; a = Z0.arad;
    rhs = a * (π * G)^3 * (M / (4π * mu1))^2 * (μmol * SE.amu / k)^4 / 3
    lo, hi = 1e-6, 1 - 1e-12
    for _ in 1:200; mid = 0.5(lo + hi); (1 - mid) / mid^4 > rhs ? (lo = mid) : (hi = mid); end
    0.5(lo + hi)
end

"""
    build_star(; M, R, X = 0.75, n = 100, nenv = 24, qc = 1.12, mfloor = 1e-10, nuc, ...) -> (; mat, ...)

Builds the TITAN state of the Eddington standard model of mass M [g] and radius R [cm] (fully ionised H/He of mass fraction X,
gas+radiation, β from the Eddington quartic). Mass grid: `n - nenv` zones geometric (ratio qc) in the interior, then `nenv` zones
log-spaced in the exterior mass fraction down to `mfloor` (surface zone mass = mfloor·M). T is floored at Tfloor
(the polytrope's surface is T → 0, not an atmosphere); F_r = L_rad(m)/4πr² with L_rad from hydrostatic radiative equilibrium.
"""
function build_star(; M = 100Msun, R = 5Rsun, X = 0.75, n = 100, nenv = 24, qc = 1.12, mfloor = 1e-10, Tfloor = 6e4, nuc = true,
                    alpha = 1.5, opac = :stellar, efrac = 0.05, dXmax = 1e-3, burnT = :end, theta = nothing, mixing = true,
                    Yinit = nothing, screen = :chugunov, kind = :star, kcap = Inf, wmax = 40.0, qsplit = 0.3)
    ξ, θp, μp = lane_emden3(); ξ1 = ξ[end]; μ1 = μp[end]
    μmol = 1 / (2X + 0.75(1 - X)); β = eddington_beta(M, μmol)
    G = Z0.G; k = Z0.kB; a = Z0.arad
    K = ((k / (μmol * SE.amu))^4 * 3 * (1 - β) / (a * β^4))^(1 / 3)
    al = R / ξ1; ρc = (K / (π * G * al^2))^1.5; Mpoly = 4π * al^3 * ρc * μ1
    Pc = K * ρc^(4 / 3); Tc = (3 * (1 - β) * Pc / a)^0.25
    # interior mass fractions of interfaces: n+1 values, centre 0 ... surface M
    # mass grid: interior zones geometric (ratio qc, capped at wmax·Δm_1) up to exterior fraction qsplit; envelope zones log-spaced in
    # the exterior mass fraction from qsplit to mfloor (neighbour mass ratio (qsplit/mfloor)^(1/nenv))
    nin = n - nenv
    w = [min(qc^(j - 1), wmax) for j in 1:nin]; mint = (1 - qsplit) .* [0; cumsum(w)] ./ sum(w)
    qe = 10 .^ range(log10(qsplit), log10(mfloor); length = nenv + 1)[2:end]
    mfrac = vcat(mint, 1 .- qe); mfrac[end] = 1.0
    m = M .* mfrac
    ξof(mm) = begin                                    # invert μ(ξ) = m/(4π al³ ρc)
        t = mm / (4π * al^3 * ρc); t >= μ1 && return ξ1
        i = searchsortedlast(μp, t); i >= length(μp) && return ξ1
        ξ[i] + (t - μp[i]) / (μp[i+1] - μp[i]) * (ξ[i+1] - ξ[i])
    end
    θof(x) = (i = searchsortedlast(ξ, x); i >= length(ξ) ? 0.0 : θp[i] + (x - ξ[i]) / (ξ[i+1] - ξ[i]) * (θp[i+1] - θp[i]))
    rI = [al * ξof(mm) for mm in m]; rI[1] = 0.0; rI[end] = R
    dvol(j) = (rI[j+1]^3 - rI[j]^3) / 3
    ρz = [(m[j+1] - m[j]) / (4π * dvol(j)) for j in 1:n]
    mid(j) = 0.5 * (m[j] + m[j+1])
    ρpoly = [ρc * max(θof(ξof(mid(j))), 0.0)^3 for j in 1:n]
    Tz = Float64[]
    for j in 1:n
        P = K * ρpoly[j]^(4 / 3); push!(Tz, max((3 * (1 - β) * P / a)^0.25, Tfloor))
    end
    Er = a .* Tz .^ 4
    get(ENV,"S4DBG","0")=="1" && (println("rho z: ", ρz[[1,2,3,n-nenv,n-2,n-1,n]], "\nT z: ", Tz[[1,2,3,n-nenv,n-2,n-1,n]], "\nr: ", rI[[1,2,3,n-1,n,n+1]]); flush(stdout))
    mat = kind === :star ? SM.StarGas(TF.MGR; X = X, Y = 1 - X, Z = 0.0, nuc, opac, efrac, dXmax, burnT, screen, mixing, kcap) :
          StellarMaterial.StellarGas(TF.MGR; X = X, Y = 1 - X, Z = 0.0)
    # radiative-equilibrium luminosity: L(m) = 4π c G m (P_rad/P_tot)/κ
    Lz = zeros(n + 1)
    for j in 2:n+1
        j <= n || continue
        ρi = ρz[min(j, n)]; Ti = Tz[min(j, n)]
        TX.ZONE[] = 1; κ = SO_kap(mat, ρi, Ti, X)
        Pr = a * Ti^4 / 3; Pg = SE.gas_pe(ρi, Ti; X, Y = 1 - X)[1]
        Lz[j] = 4π * Z0.c * G * m[j] * Pr / (Pr + Pg) / κ
    end
    Fz = [j == 1 ? 0.0 : Lz[j] / (4π * rI[j]^2) for j in 1:n+1]; Fz[end] = Lz[n] / (4π * R^2)
    PS.build!(; rI, rho = ρz, T = Tz, u = zeros(n + 1), geom = 2, rad = true, Er, Fr = Fz, lgrid = 3, inner = :wall, outer = :free,
              pext = 0.0, selfgrav = true, material = mat, convection = true, gam = 5 / 3, gmu = 1.0)
    TX.REGULAR_CENTRE[] = true
    TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TX.ALPHA_MLT[] = alpha
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    theta === nothing || TF.setscalar!(:thet, theta)
    if Yinit !== nothing; for kk in 1:size(mat.Y, 2); SM.setY!(mat, kk, Yinit); end; end
    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    (; mat, M, R, β, K, ρc, Tc, Mpoly, μmol, ξ1, al, rI, m, ρz, Tz, n)
end
SO_kap(mat, ρ, T, X) = StellarOpacity.kappa_R(ρ, T, X)

"Stellar global quantities from the current TITAN state (L_surface = 4π r² F at the outer interface)."
function star_diag(st)
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    rn, dn, tn, frn, ern, egn, pgn = (TF.var(s) for s in (:rn, :dn, :tn, :frn, :ern, :egn, :pgn))
    Ls = 4π * rn[b+1]^2 * frn[b+1]
    R = rn[b+1]; Teff = (Ls / (4π * R^2 * Z0.σSB))^0.25
    (; t = TF.getscalar(:timen), R, Ls, Teff, Tc = tn[a], ρc = dn[a], Pc = pgn[a] + ern[a] / 3, Lnuc = SM.nuclear_power(st.mat))
end
