# Population III project, M2c: 1D collapse of the (provisional) RHFC02 P100 cloud with the RHFC02 chemistry and local
# cooling (RA04's simplified transfer: thick H₂ lines × min(1, (n/8e9)^−0.45) + CIE with escape factor) — the
# configuration of RA04's own 1D test (RA04 §2.3.3, Figs. 4–5), here in TITAN.
#   TITAN_NOLIB=1 julia --project=env popiii_m2c.jl [cooling] [nshell]     cooling ∈ ra04 | ra04lines | thin
# Initial state (provisional, see POPIII_NOTES M2b): n = 1.5 polytrope shape, ρ_c/m_H = 2.6e6, cut at 1e-3 ρ_c, scaled
# to 100 M☉; T = 270 K (ρ/ρ_c)^{2/3}; at rest; f_H2 = 5e-4, x_e = 1e-10. Geometric shell masses with the central shell
# 3e-7 M☉ (RHFC02 §3.1). Regular centre, Lagrangian, self-gravity, free surface at the edge pressure.
# Log: lib/popiii/m2c_<cooling>_n<nshell>.log (central track every 0.25 dex of n_c, ON98 Fig. 3a for comparison).
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl"))
include(joinpath(@__DIR__, "SCVH.jl")); include(joinpath(@__DIR__, "PopIIIEqEOS.jl"))
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl"))
include(joinpath(@__DIR__, "PopIIIMaterial.jl")); include(joinpath(@__DIR__, "H2Lines.jl"))
include(joinpath(@__DIR__, "StellarOpacity.jl"))
# S6-C1: HOT=1: stellar EOS, nuclear network and convective mixing in the hot core (T > HOTT1), see PopIIIMaterial.HOT
const HOTON = get(ENV, "HOT", "0") == "1"
if HOTON
    for f in ("StellarEOS", "NucNet", "StellarMaterial", "StarMaterial"); include(joinpath(@__DIR__, f * ".jl")); end
end
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const PC = PrimordialChem; const PM = PopIIIMaterial; const Z0 = Z0Physics
using ForwardDiff
# INTEGRATOR=trbdf2 (default) | theta: second-order L-stable TR-BDF2 (Trbdf2.jl, G2) or TITAN's theta scheme
include(joinpath(@__DIR__, "Trbdf2.jl"))
# default since G2b: θ = 0.55 lets the z = 19 core oscillate at 1e11–1e13 (three-body heating phase) where TR-BDF2 is smooth
get(ENV, "INTEGRATOR", "trbdf2") == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
const Msun = 1.98847e33

cool = isempty(ARGS) ? "ra04" : ARGS[1]; nsh = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100
cf = cool == "thin" ? PC.cool_H2_thin : cool == "ra04lines" ? PC.cool_H2_ra04 :
     cool == "lines" ? PC.cool_CIE_ra04 :                     # H₂ lines by transfer (H2Lines), continuum: RA04 CIE
     (ρ, T, y) -> PC.cool_H2_ra04(ρ, T, y) + PC.cool_CIE_ra04(ρ, T, y)
haskey(ENV, "TCMB") && (cf = PC.cool_floor(cf, parse(Float64, ENV["TCMB"])))
# S6-B: RAD=1 (with RADFROM=checkpoint) continues a run with TITAN's radiation (E_r, F_r, Eddington factors) and the
# StellarOpacity continuum; the RA04 CIE escape-factor term is then not used (the continuum is in the transfer), the
# H₂ lines stay as Λext
const RAD = get(ENV, "RAD", "0") == "1"
RAD && (cf = (ρ, T, y) -> zero(T))
io = open(joinpath(TF.ROOT, "lib", "popiii", "m2c_$(get(ENV, "IC", "p100") == "p100" ? "" : get(ENV, "IC", "") * "_")$(cool)_n$(nsh)$(get(ENV, "TAG", "")).log"), haskey(ENV, "RESTART") ? "a" : "w")
pr(f, a...) = (s = Printf.format(Printf.Format(f), a...); print(s); print(io, s); flush(io))

# ---- n = 1.5 Lane–Emden θ(ξ), integrated to θ^1.5 = 1e-3
function lane_emden(n; cut = 1e-3)
    f(ξ, y) = (y[2], -y[1]^n - 2y[2] / ξ)
    ξ = 1e-4; y = (1 - ξ^2 / 6, -ξ / 3); h = 1e-4; out = [(0.0, 1.0, 0.0)]
    while max(y[1], 0.0)^n > cut
        k1 = f(ξ, y); k2 = f(ξ + h / 2, y .+ h / 2 .* k1); k3 = f(ξ + h / 2, y .+ h / 2 .* k2); k4 = f(ξ + h, y .+ h .* k3)
        y = y .+ h / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); ξ += h
        push!(out, (ξ, y[1], -ξ^2 * y[2]))          # (ξ, θ, ξ²|θ'| ∝ m)
    end
    out
end
ic = get(ENV, "IC", "p100")
dm1 = 3e-7Msun
q_of(M, n) = let lo = 1.0001, hi = 2.0
    f(q) = dm1 * (q^n - 1) / (q - 1) - M
    for _ in 1:200; m_ = (lo + hi) / 2; f(m_) > 0 ? (hi = m_) : (lo = m_); end
    (lo + hi) / 2
end
if ic in ("abn02", "abn02hs", "abn02z19")
    # ABN02 (Abel, Bryan & Norman 2002, Science 295, 93; arXiv astro-ph/0112088) Fig. 2, profile "9 Myr later"
    # (red dotted with circles), read off the 400-dpi rendering [read-src, approximate ±0.05 dex]: particle number
    # density n(r) (panel A), T(M), H₂ mass fraction (panel C, lower curve), v_r(M) (panel E). ρ = n μ m_H, μ = 1.22.
    z19 = startswith(ic, "abn02z19")
    if z19
        # z = 19 profile (blue solid lines of ABN02 Fig. 2), read off the 400-dpi rendering [read-src, approximate]
        lr_tab = [-3.0, -1.24, -0.59, -0.27, 0.06, 0.38, 0.71, 0.95, 1.20, 1.44, 1.68, 1.97]   # log r [pc]
        ln_tab = [3.05, 3.05, 3.05, 2.90, 2.61, 2.01, 1.27, 0.60, 0.30, -0.07, -0.59, -1.48]   # log n [cm⁻³]
        lM_tab = [-1.0, 1.19, 2.28, 2.64, 2.86, 3.07, 3.29, 3.51, 3.72, 3.91, 4.16, 4.34]       # log M [M☉]
        T_tab  = [178.0, 180.0, 178.0, 178.0, 194.0, 247.0, 416.0, 650.0, 855.0, 1000.0, 890.0, 725.0]
        lf_tab = [-3.18, -3.14, -3.18, -3.20, -3.25, -3.33, -3.45, -3.60, -3.76, -3.90, -4.0, -4.1]
        v_tab  = [0.0, 0.08, 0.22, 0.14, 0.0, -0.38, -1.0, -0.54, -0.49, -0.27, -1.51, -2.6]    # km/s
        lrd = [-2.0, -1.24, -0.76, -0.43, 0.63, 1.20, 1.52, 2.0]                                # DM, log r [pc]
        lρd = [5.0, 4.39, 4.17, 3.65, 2.23, 1.71, 1.19, -0.59]                                  # log ρ_DM [GeV cm⁻³]
    else
        lr_tab = [-3.0, -1.57, -1.28, -1.08, -0.84, -0.59, -0.35, 0.06, 0.87]          # log r [pc]
        ln_tab = [5.83, 5.81, 5.40, 5.25, 4.73, 3.99, 3.47, 2.87, 0.97]                # log n [cm⁻³]
        lM_tab = [-2.0, 1.5, 2.0, 2.5, 3.0, 3.4, 3.7, 4.0]                             # log M [M☉]
        T_tab  = [252.0, 252.0, 229.0, 205.0, 164.0, 355.0, 1120.0, 1000.0]
        lf_tab = [-2.90, -2.90, -2.98, -3.11, -3.18, -3.37, -3.65, -4.0]
        v_tab  = [-0.03, -0.32, -0.49, -0.19, -0.81, -0.70, -1.73, -0.95]               # km/s
        lrd = [-2.0, -0.35, 0.06, 0.47, 0.87, 1.28, 1.60, 2.0]; lρd = [5.07, 3.42, 2.90, 2.68, 2.16, 1.49, 0.67, -0.59]
    end
    # monotone cubic Hermite (Fritsch–Carlson) through the read-off points: no kinks (a piecewise-linear log n(r) put
    # spurious pressure jumps of up to 40× gravity at the knots, M2e)
    function pchip(xs, ys)
        n = length(xs); d = diff(ys) ./ diff(xs); m = zeros(n)
        m[1] = d[1]; m[n] = d[n-1]
        for i in 2:n-1; m[i] = d[i-1] * d[i] <= 0 ? 0.0 : (d[i-1] + d[i]) / 2; end
        for i in 1:n-1
            d[i] == 0 && (m[i] = 0.0; m[i+1] = 0.0; continue)
            a, b = m[i] / d[i], m[i+1] / d[i]; h = a^2 + b^2
            h > 9 && (t = 3 / sqrt(h); m[i] = t * a * d[i]; m[i+1] = t * b * d[i])
        end
        x -> begin
            x <= xs[1] && return ys[1] + m[1] * (x - xs[1]) * 0   # constant outside
            x >= xs[n] && return ys[n] + m[n] * (x - xs[n])        # linear (log-log power law) beyond the last point
            j = searchsortedlast(xs, x); h_ = xs[j+1] - xs[j]; t = (x - xs[j]) / h_
            (2t^3 - 3t^2 + 1) * ys[j] + (t^3 - 2t^2 + t) * h_ * m[j] + (-2t^3 + 3t^2) * ys[j+1] + (t^3 - t^2) * h_ * m[j+1]
        end
    end
    interp(x, xs, ys) = pchip(xs, ys)(x)
    pc = 3.0857e18; μ_abn = 1.22
    nfit = pchip(lr_tab, ln_tab); nr(r) = 10^nfit(log10(r / pc))
    # dark matter (ABN02 Fig. 2A, the upper line between 1e-2 and 200 pc: DM density in GeV cm⁻³, read off) [read-src,
    # approximate]; inside 0.45 pc continued with slope −1 (not legible under the gas curves)
    dmfit = pchip(lrd, lρd); GeV = 1.78266e-24
    if get(ENV, "DM", "1") == "1"
        rgd = [pc * 10^x for x in range(-4.0, 2.3, length = 20001)]; Mgd = zeros(length(rgd))
        ρdm(r) = 10^dmfit(log10(r / pc)) * GeV
        for j in 2:length(rgd)
            rm = 0.5 * (rgd[j] + rgd[j-1]); Mgd[j] = Mgd[j-1] + 4π * rm^2 * ρdm(rm) * (rgd[j] - rgd[j-1])
        end
        Mgd .+= 4π / 3 * rgd[1]^3 * ρdm(rgd[1])
        lrg = log10.(rgd); lMg = log10.(Mgd ./ 4π)                         # per steradian
        TX.EXTMASS[] = r -> (x = log10(max(ForwardDiff.value(r), rgd[1])); j = clamp(searchsortedlast(lrg, x), 1, length(lrg) - 1);
                             10^(lMg[j] + (x - lrg[j]) / (lrg[j+1] - lrg[j]) * (lMg[j+1] - lMg[j])) * one(r))
    end
    # M(r) by integrating ρ(r) on a fine log grid
    rg = [pc * 10^x for x in range(-4.0, 1.5, length = 20001)]; Mg = zeros(length(rg))
    for j in 2:length(rg)
        rm = 0.5 * (rg[j] + rg[j-1]); Mg[j] = Mg[j-1] + 4π * rm^2 * nr(rm) * μ_abn * PC.mH * (rg[j] - rg[j-1])
    end
    Mg .+= 4π / 3 * rg[1]^3 * nr(rg[1]) * μ_abn * PC.mH
    M = parse(Float64, get(ENV, "MCLOUD", z19 ? "1e4" : "1000")) * Msun
    if ic in ("abn02hs", "abn02z19")
        # hydrostatic version: ABN02's central density, T(M) and dark matter; dP/dr = −ρ G (M_gas + M_DM)/r² with
        # P = ρ c²(M), c² = k T(M)/(μ m_H): dρ/dr = −[ρ G M_tot/r² + ρ² 4π r² dc²/dM]/c². RK4 in ln r.
        Tfit = pchip(lM_tab, T_tab); c2(m) = PC.kB * Tfit(log10(max(m / Msun, 1e-2))) / (μ_abn * PC.mH)
        dc2dM(m) = (c2(m * 1.001) - c2(m / 1.001)) / (m * (1.001 - 1 / 1.001))
        Mdm(r) = TX.EXTMASS[] === nothing ? 0.0 : 4π * TX.EXTMASS[](r)
        function rhs_hs(r, ρ, m)
            dm = 4π * r^2 * ρ
            dρ = -(ρ * PC.G * (m + Mdm(r)) / r^2 + ρ * dm * dc2dM(m)) / c2(m)
            (dρ, dm)
        end
        function integrate_hs(Mmax)
            r = 1e-4 * pc; ρ = 10^ln_tab[1] * μ_abn * PC.mH; m = 4π / 3 * r^3 * ρ
            rs = [0.0, r]; Ms = [0.0, m]; ρs = [ρ, ρ]
            while m < Mmax
                h = 0.002r
                k1 = rhs_hs(r, ρ, m); k2 = rhs_hs(r + h / 2, ρ + h / 2 * k1[1], m + h / 2 * k1[2])
                k3 = rhs_hs(r + h / 2, ρ + h / 2 * k2[1], m + h / 2 * k2[2]); k4 = rhs_hs(r + h, ρ + h * k3[1], m + h * k3[2])
                ρ += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]); m += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]); r += h
                push!(rs, r); push!(Ms, m); push!(ρs, ρ)
            end
            rs, Ms, ρs
        end
        global rg, Mg, ρg = integrate_hs(z19 ? 1e3Msun : 1.2M)
        if z19
            # outside r_join the gas is not static in ABN02 (virial-shock heating, infall −0.4…−2.6 km/s): splice the
            # read-off n(r) on continuously at r_join (default 10^0.87 = 7.4 pc, where the two agree within 0.12 dex) and
            # integrate the mass outward; velocities there from v(M) (zero inside)
            # G8: r_join automatic (LRJOIN=auto, default): the outermost radius inside which the hydrostatic interior stays
            # within DLNJOIN (default 0.15 dex, the read-off accuracy of ABN02 Fig. 2) of ABN02's read-off n(r) everywhere ("best agreement" alone is degenerate: in
            # the flat core both agree trivially); the slope mismatch there is printed
            lj = get(ENV, "LRJOIN", "auto")
            if lj == "auto"
                lnh(x) = (j = clamp(searchsortedfirst(rg, pc * 10^x), 2, length(rg)); log10(ρg[j] / (μ_abn * PC.mH)))
                xs = range(log10(rg[3] / pc) + 0.1, log10(rg[end] / pc) - 0.1, length = 2000)
                δj = parse(Float64, get(ENV, "DLNJOIN", "0.15"))
                kbad = findfirst(x -> abs(lnh(x) - nfit(x)) > δj, xs)
                xbest = kbad === nothing ? xs[end] : xs[max(kbad - 1, 1)]
                dsl = (lnh(xbest + 0.05) - lnh(xbest - 0.05)) / 0.1 - (nfit(xbest + 0.05) - nfit(xbest - 0.05)) / 0.1
                println("G8 automatic splice: log r_join = ", round(xbest, digits = 3), " (|Δ log n| ≤ $δj inside; slope mismatch ", round(dsl, digits = 2), ")")
                rj = pc * 10^xbest
            else
                rj = pc * 10^parse(Float64, lj)
            end
            k = searchsortedlast(rg, rj); rg = rg[1:k]; Mg = Mg[1:k]; ρg = ρg[1:k]
            fac = ρg[end] / nr(rg[end])
            while Mg[end] < 1.2M && rg[end] < 3e3pc
                r0 = rg[end]; r1 = r0 * 1.002; rm = 0.5 * (r0 + r1); ρm = fac * nr(rm)
                push!(rg, r1); push!(ρg, fac * nr(r1)); push!(Mg, Mg[end] + 4π * rm^2 * ρm * (r1 - r0))
            end
            global Mjoin = Mg[k]
            println("splice at r = ", round(rj / pc, digits = 2), " pc, M_join = ", round(Mjoin / Msun, digits = 1), " Msun, density factor ", round(fac / (μ_abn * PC.mH), digits = 3))
        end
        for x in (z19 ? (-0.84, -0.35, 0.06, 0.87) : (-1.57, -1.08, -0.59, 0.06, 0.87))   # compare with the read-off
            j = searchsortedfirst(rg, pc * 10^x); j <= length(rg) &&
                println("hydrostatic vs ABN02: log r = $x  log n = ", round(log10(ρg[j] / (μ_abn * PC.mH)), digits = 2), "  (ABN02 ", round(nfit(x), digits = 2), ")  log M = ", round(log10(Mg[j] / Msun), digits = 2))
        end
        println("hydrostatic: M = 1000 Msun at r = ", rg[searchsortedfirst(Mg, M)] / pc, " pc")
    end
    q = q_of(M, nsh)
    mI = [0.0; cumsum([dm1 * q^(j - 1) for j in 1:nsh])]; mI[end] = M
    r_of_m(m) = m <= 0 ? 0.0 : (j = searchsortedfirst(Mg, m); j = clamp(j, 2, length(Mg));
        rg[j-1] + (m - Mg[j-1]) / (Mg[j] - Mg[j-1]) * (rg[j] - rg[j-1]))
    rI = [r_of_m(m) for m in mI]; rI[1] = 0.0
    rho = [(mI[j+1] - mI[j]) / (4π / 3 * (rI[j+1]^3 - rI[j]^3)) for j in 1:nsh]
    mz = [0.5 * (mI[j] + mI[j+1]) / Msun for j in 1:nsh]
    Tz = [interp(log10(m), lM_tab, T_tab) for m in mz]
    fz = [10^interp(log10(m), lM_tab, lf_tab) for m in mz]
    uI = ic == "abn02z19" ? [m < Mjoin ? 0.0 : 1e5 * interp(log10(m / Msun), lM_tab, v_tab) * min(1.0, (m - Mjoin) / Mjoin) for m in mI] :
         ic == "abn02hs" ? zeros(nsh + 1) : [1e5 * interp(log10(max(m / Msun, 1e-2)), lM_tab, v_tab) for m in mI]; uI[1] = 0.0
    ρc = rho[1]; Tc0 = Tz[1]
    xe0 = parse(Float64, get(ENV, "XE0", z19 ? "1e-4" : "1e-7"))   # x_e [assumed; ABN02 does not show it]
    yz = [(xe0, 0.0, f / 2, 0.0) for f in fz]
else
le = lane_emden(1.5); ξs, _, μs = le[end]
le = lane_emden(1.5); ξs, _, μs = le[end]
M = 100Msun; ρc = 2.6e6 * PC.mH; Tc0 = 270.0
a_ = cbrt(M / (4π * ρc * μs))                      # radius unit: M = 4π ρc a³ ξ²|θ'| at the cut
q = q_of(M, nsh)
mI = [0.0; cumsum([dm1 * q^(j - 1) for j in 1:nsh])]; mI[end] = M
# radius of each interface from m(ξ) (linear interpolation in the tabulated profile)
mξ = [4π * ρc * a_^3 * e[3] for e in le]; ξt = [e[1] for e in le]
function ξ_of_m(m)
    m <= 0 && return 0.0
    j = searchsortedfirst(mξ, m); j = clamp(j, 2, length(mξ))
    ξt[j-1] + (m - mξ[j-1]) / (mξ[j] - mξ[j-1]) * (ξt[j] - ξt[j-1])
end
rI = [a_ * ξ_of_m(m) for m in mI]; rI[1] = 0.0
rho = [(mI[j+1] - mI[j]) / (4π / 3 * (rI[j+1]^3 - rI[j]^3)) for j in 1:nsh]
Tz = [Tc0 * (ρ / ρc)^(2 / 3) for ρ in rho]
uI = zeros(nsh + 1); yz = nothing
end
y0 = (1e-10, 0.0, 2.5e-4, 0.0)
k5s = get(ENV, "K5", "pss")
k5f = k5s == "db" ? (T -> PC.k5_db(T)) : k5s == "ra04" ? PC.k5_ra04 : PC.k5_pss
ystore = fill(y0, 2000)
if yz !== nothing
    a0 = 5                                         # ngrs of the captured SOD state used by build!
    for j in 1:nsh; ystore[a0+j-1] = yz[j]; end
    for k in 1:a0-1; ystore[k] = yz[1]; end; for k in a0+nsh:length(ystore); ystore[k] = yz[end]; end
end
const IMPCOOL = get(ENV, "IMPCOOL", "1") in ("1", "2")               # G7: transfer cooling implicit in T (linearised in ln T)
const DTCOOL = parse(Float64, get(ENV, "DTCOOL", "1.0"))
const IMPALWAYS = get(ENV, "IMPCOOL", "1") == "2"                  # IMPCOOL=2: slope every step
const NIMP = Ref(0)                                                  # steps with implicit line cooling
const CHEM = get(ENV, "CHEM", "rhfc02")     # rhfc02: PrimordialGas (RHFC02 network, ρ_eq switch); net / netmodern / netmodernli: NetGas (G5)
if CHEM == "rhfc02"
    global mat = PM.PrimordialGas(ystore, PC.rates_rhfc02(; k5 = k5f), cf, parse(Float64, get(ENV, "RHOEQ", "5e-11")), get(ENV, "EQEOS", "scvh") == "scvh")
else
    # G5: ChemNet network; equilibrium table where the ideal-gas pressure ≥ PSWITCH (EOS validity, RHFC02's 1e4 dyn cm⁻²)
    netw, phot = CHEM == "net" ? (ChemRates.rhfc02(; k5 = k5f), ChemNet.Photo[]) : ChemRates.modern(; li = CHEM == "netmodernli")
    tcmb = parse(Float64, get(ENV, "TCMB", "0"))
    Bν_(ν, T) = 2ChemNet.hP * ν^3 / ChemNet.c^2 / expm1(ChemNet.hP * ν / (ChemNet.kB * T))
    envn = (Trad = tcmb, J = tcmb > 0 ? (ν -> Bν_(ν, tcmb)) : (ν -> 0.0), static = true)   # same J in every zone: tabulated photo-rates
    ytup(d) = (d[:Hp], d[:Hm], d[:H2], d[:H2p])
    coolnet = (ρ, T, d) -> cf(ρ, T, d === nothing ? PopIIIEqEOS.composition(ρ, T) : ytup(d))
    xinit(y) = [s === :Hp ? y[1] : s === :Hm ? y[2] : s === :H2 ? y[3] : s === :H2p ? y[4] : s === :Lip ? 0.0 : 0.0 for s in netw.free]
    global mat = PM.NetGas([xinit(y) for y in ystore], netw, phot, coolnet; env = envn, pswitch = parse(Float64, get(ENV, "PSWITCH", "1e4")))
end
"Composition of zone k as the PrimordialChem tuple (x_H⁺, x_H⁻, x_H₂, x_H₂⁺), for either material."
ycomp(k) = mat isa PM.NetGas ? (d = ChemNet.full(mat.x[k], mat.net); (d[:Hp], d[:Hm], d[:H2], d[:H2p])) : mat.y[k]
"Equilibrium-regime flag of zone k."
eqflag(k) = mat isa PM.NetGas ? mat.eqzone[k] : TF.var(:dn)[k] >= mat.ρeq
# artificial viscosity: viscous length l = ql1 r ≈ 2 Δr (the SOD state's ql0 = 1e-4 cm is zero at these scales, M2g);
# Δr/r ≈ 0.04 (core) … 0.12 (envelope) on the geometric mass grid; quadratic term only (cq2), as PORTING_NOTES P7
av = (ql0 = 0.0, ql1 = parse(Float64, get(ENV, "AVQL1", "0.15")), cq1 = 0.0, cq2 = parse(Float64, get(ENV, "AVCQ2", "1.0")))
PS.build!(; rI, rho, T = Tz, u = uI, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = 0.0,
          selfgrav = true, material = mat, av)
TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))    # G1: viscous length = AVZONE × zone width (0 = ql1 r̄ above)
TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
TX.PEXT[] = TF.var(:pgn)[b]
if get(ENV, "FORCECHECK", "0") == "1"           # initial force balance: −r²Δp / (G m m̄/r²) per interface
    pg = TF.var(:pgn); dn = TF.var(:dn); rn = TF.var(:rn); dv = TF.var(:dvoln); xm = TX.XMINT[]
    for k in a+1:4:b
        mext = TX.EXTMASS[] === nothing ? 0.0 : TX.EXTMASS[](rn[k])
        mb = 0.5 * (dn[k-1] * dv[k-1] + dn[k] * dv[k]); g = 4π * Z0.G * (xm[k] + mext) / rn[k]^2 * mb
        pr("force k=%3d m=%.3e Msun r=%.3e pc  pressure/gravity = %.3f  T=%.1f n=%.3e\n", k, 4π * xm[k] / Msun, rn[k] / 3.0857e18,
           -rn[k]^2 * (pg[k] - pg[k-1]) / g, TF.var(:tn)[k], dn[k] / (1.22 * PC.mH))
    end
end
TF.setscalar!(:thet, parse(Float64, get(ENV, "THETA", "0.55")))   # 0.5 lets the centre ring (M2c)
TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.01"))
tff = sqrt(3π / (32Z0.G * ρc))
TF.setscalar!(:dtime, 1e-3tff); TF.setscalar!(:tmax, 1e3tff)
pr("M2c collapse [IC %s]: cooling %s, %d shells (q = %.4f, dm1 = 3e-7 Msun), R = %.3e cm, rho_c/rho_edge = %.1f, T_edge = %.1f K, p_ext = %.3e\n",
   ic, cool, nsh, q, rI[end], rho[1] / rho[end], Tz[end], TX.PEXT[])
pr("theta = %.2f, C_dyn = %g, k5 = %s, chemistry %s, %s, AV ql1 = %g cq2 = %g, central t_ff0 = %.3e yr\n", TF.getscalar(:thet), TC.DT_DYN[], k5s, CHEM, mat isa PM.NetGas ? @sprintf("equilibrium table where p ≥ %g", mat.pswitch) : @sprintf("rho_eq = %g, eq EOS %s", mat.ρeq, mat.eqtable ? "Saha+SCVH95" : "Saha"), av.ql1, av.cq2, tff / 3.156e7)
# ON98 Fig. 3a (run A), read off: (log n_c, log T_c)
on98 = [(6.3, 2.40), (8.0, 2.60), (10.0, 2.75), (12.0, 3.05), (15.0, 3.30), (16.0, 3.35), (20.0, 3.80)]
on98T(ln) = (j = findfirst(p -> p[1] >= ln, on98); j === nothing ? NaN : j == 1 ? 10^on98[1][2] :
             10^(on98[j-1][2] + (ln - on98[j-1][1]) / (on98[j][1] - on98[j-1][1]) * (on98[j][2] - on98[j-1][2])))
pr("%7s %8s %10s %9s %9s %10s %10s %9s %9s\n", "step", "log n_c", "t [yr]", "T_c", "ON98 T_c", "f_H2,c", "u_c [km/s]", "max dT/T", "dt [s]")
nextln = Ref(log10(rho[1] / PC.mH) + 0.25)
using Serialization
const SNAPFILE = joinpath(TF.ROOT, "lib", "popiii", "snaps_$(get(ENV, "IC", "p100"))_$(cool)_n$(nsh)$(get(ENV, "TAG", "")).jls")
const SNAPS = get(ENV, "SNAPS", "0") == "1" ? Any[] : nothing
const SNAPNEXT = Ref(log10(rho[1] / PC.mH)); const SNAPEVERY = parse(Int, get(ENV, "SNAPEVERY", "0"))
# NLTE (line source functions): 1 = GP98/HM79 factor on S (old default), 0 = LTE, ep = statistical equilibrium with escape
# probabilities (H2Levels), ali = statistical equilibrium coupled to the transfer (ALI, G6c), hyb = ALI every NALI steps
# and EP × (ALI/EP per zone, from the last ALI) in between, the interval adaptive (see linecool)
const ALIWMIN = parse(Float64, get(ENV, "ALIWMIN", "0")); const ALIGROW = parse(Float64, get(ENV, "ALIGROW", "0.01")); const ALISHRINK = parse(Float64, get(ENV, "ALISHRINK", "0.05")); const NLTE = get(ENV, "NLTE", "1"); const NALI = parse(Int, get(ENV, "NALI", "20")); const ALIFAC = Ref{Vector{Float64}}(Float64[]); const ALILAST = Ref(-10^9); const ALIINT = Ref(1); const NALICALL = Ref(0)
NLTE in ("ep", "ali", "hyb") && include(joinpath(@__DIR__, "H2Levels.jl"))
const HLV = NLTE in ("ep", "ali", "hyb") ? H2Levels.load() : nothing
"""
S6: H₂ line transfer outside the continuum photosphere only (RAD runs, LTAUCUT > 0 [10]): the zones deeper than continuum
optical depth LTAUCUT from the surface are merged into one transparent central zone without H₂ for the line transfer (their
line cooling is 0: the continuum carries the energy there and H₂ is dissociated); the cost ∝ zones² × bins drops (the hot core
also set the velocity span of the frequency grid). Returns (Λ, Lr) on the full grid.
"""
function linecool(rr, uu, dd, tt, nh2, k)
    jc = 0
    if RAD && LTAUCUT > 0
        ch = TF.var(:chifn); τ = 0.0
        for j in nsh:-1:1
            τ += ch[a+j-1] * dd[j] * (rr[j+1] - rr[j])
            τ > LTAUCUT && (jc = j; break)
        end
        jc = jc >= nsh - 4 ? 0 : jc                  # (keep at least a few zones)
    end
    jc < 2 && return linecool_full(rr, uu, dd, tt, nh2, k, 0)
    r2 = [0.0; rr[jc+1:end]]; u2 = [0.0; uu[jc+1:end]]
    m_c = sum(dd[j] * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:jc); d2 = [m_c / (4π / 3 * rr[jc+1]^3); dd[jc+1:end]]
    t2 = [tt[jc]; tt[jc+1:end]]; h2 = [0.0; nh2[jc+1:end]]
    Λ2, L2, x2 = linecool_full(r2, u2, d2, t2, h2, k, jc)
    Λ = zeros(nsh); Λ[jc+1:end] .= Λ2[2:end]; Lr = zeros(nsh + 1); Lr[jc+1:end] .= L2[2:end]
    Λ, Lr, x2
end
const LTAUCUT = parse(Float64, get(ENV, "LTAUCUT", "10"))
"Line cooling per zone (and Lr) for zone temperatures tt; step k (for NLTE=hyb); jc = number of core zones merged away."
function linecool_full(rr, uu, dd, tt, nh2, k, jc)
    nz_ = length(dd); zk(j) = j == 1 && jc > 0 ? a + jc - 1 : a + j - 1 + max(jc - 1, 0)   # TITAN zone of reduced zone j
    if NLTE in ("ep", "ali", "hyb")
        yc = [ycomp(zk(j)) for j in 1:nz_]; nHn = dd ./ PC.μH
        nl = (LV = HLV, nH = [PC.xH(yc[j]) * nHn[j] for j in 1:nz_], nHe = PC.yHe .* nHn, nHp = [yc[j][1] * nHn[j] for j in 1:nz_])
        NLTE != "hyb" && return H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = NLTE == "ali" ? :ali : :ep)
        tq = time(); Λe, Lr, x = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ep)
        TIMING && k >= 0 && pr("        EP call %.1f s (nep %d)\n", time() - tq, H2Lines.LAST[].nep)
        # ALI refresh interval adaptive (≤ NALI): halved when the ALI/EP factor moved by > 5 % since the last refresh,
        # doubled when < 1 % — a factor frozen while the lines turn thick jumps at the refresh and kicks the core
        if k >= 0 && (k - ALILAST[] >= ALIINT[] || length(ALIFAC[]) != nz_)
            tq = time(); Λa, Lra, xa = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ali)
            TIMING && pr("        ALI call %.1f s (nep %d, ALI iterations %d, converged %s)\n", time() - tq, H2Lines.LAST[].nep, H2Lines.LAST[].iters, H2Lines.LAST[].converged)
            fnew = [Λe[j] != 0 && Λa[j] / Λe[j] > 0 ? Λa[j] / Λe[j] : 1.0 for j in 1:nz_]
            if length(ALIFAC[]) == nz_
                # (S6-C2, ALIWMIN > 0, e.g. 1e-4) only over zones carrying ≥ ALIWMIN of the line luminosity: tiny zones (rezoned shock) with
                # negligible cooling otherwise keep the interval at 1 and the ALI (≈ 50× an EP step) runs every step
                wv = [abs(Λe[j]) * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:nz_]; wmin = ALIWMIN * sum(wv)
                δ = maximum((abs(log(fnew[j] / ALIFAC[][j])) for j in 1:nz_ if wv[j] >= wmin); init = 0.0)
                ALIINT[] = δ > ALISHRINK ? max(1, ALIINT[] ÷ 2) : δ < ALIGROW ? min(NALI, 2ALIINT[]) : ALIINT[]
            end
            ALIFAC[] = fnew; ALILAST[] = k; NALICALL[] += 1
            return Λa, Lra, xa
        end
        return Λe .* (length(ALIFAC[]) == nz_ ? ALIFAC[] : 1.0), Lr, x     # (perturbed call right after the reduced grid changed: EP)
    end
    # NLTE=1: source function × f = 1/(1 + Λ_LTE/(n_H Λ_0)) (the GP98/HM79 combination of PrimordialChem.cool_H2_thin)
    fnl = [1 / (1 + PC.Λ_H2_LTE_HM(tt[j]) / max(PC.xH(ycomp(zk(j))) * dd[j] / PC.μH * PC.Λ_H2_0_GP(tt[j]), 1e-300)) for j in 1:nz_]
    H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, fnlte = NLTE == "1" ? fnl : nothing)
end
const HLL = cool == "lines" ? H2Lines.load() : nothing; const NBIN = parse(Int, get(ENV, "NBIN", "30")); const LLINE = Ref(0.0); const JCMIN = Ref(1)
Tstop = parse(Float64, get(ENV, "TSTOP", "6000")); nstop = parse(Float64, get(ENV, "NSTOP", "1e18"))
# Checkpoint / restart (S6-A). SAVE=name: checkpoint at the end and every CKEVERY [500] steps (alternating
# <name>_auto0/_auto1); RESTART=name: the set-up above runs as usual (same ENV), then the TITAN state, the material's
# per-zone state, the driver's NLTE/cooling bookkeeping, the TR-BDF2 counters and the step number are replaced by the
# checkpoint's, so the continuation is the uninterrupted run (checked bit for bit). Files: CKDIR [lib/popiii/ck].
const CKDIR = get(ENV, "CKDIR", joinpath(TF.ROOT, "lib", "popiii", "ck"))
const KOFF = Ref(0); const SG = Ref{Any}(nothing)
matstate() = mat isa PM.NetGas ? (; x = deepcopy(mat.x), eqzone = copy(mat.eqzone), Λext = copy(mat.Λext), Tref = copy(mat.Tref), slope = copy(mat.slope), Y = HOTON ? copy(SG[].Y) : nothing) :
             error("checkpoint: only the NetGas material (CHEM=net*) is supported")
function checkpoint(name, k)
    mkpath(CKDIR)
    TF.save_state(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    serialize(joinpath(CKDIR, "m2c_$(name).aux"), (; k, mat = matstate(), ALIFAC = copy(ALIFAC[]), ALILAST = ALILAST[], ALIINT = ALIINT[],
        NALICALL = NALICALL[], NIMP = NIMP[], LLINE = LLINE[], nextln = nextln[], SNAPNEXT = SNAPNEXT[], SNAPS,
        trb = (; STATS = copy(Trbdf2.STATS), BE_LEFT = Trbdf2.BE_LEFT[]), PEXT = TX.PEXT[], env = Dict(k_ => v_ for (k_, v_) in ENV if k_ in CKENV)))
end
# settings that must agree between a run and its restart (the set-up and the physics); the rest may change
const CKENV = ("IC", "TCMB", "AVZONE", "AVQL1", "AVCQ2", "CHEM", "NLTE", "NALI", "NBIN", "K5", "MCLOUD", "LRJOIN", "DLNJOIN", "DM", "PSWITCH",
               "IMPCOOL", "DTCOOL", "INTEGRATOR", "THETA", "XE0", "RAD")   # (not CDYN: the accretion phase relaxes the dynamical Δt limit at a restart)
function restart!(name; allow = ())
    aux = deserialize(joinpath(CKDIR, "m2c_$(name).aux"))
    for (k_, v_) in aux.env
        k_ in allow && continue
        get(ENV, k_, nothing) == v_ || error("restart: ENV $k_ = $(repr(get(ENV, k_, nothing))) differs from the checkpoint's $(repr(v_))")
    end
    for k_ in CKENV; k_ in allow && continue; haskey(ENV, k_) && !haskey(aux.env, k_) && error("restart: ENV $k_ is set but was not in the checkpoint's run"); end
    TF.load_state!(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1   # (zone count may differ after rezoning)
    s = aux.mat; mat.x .= deepcopy(s.x); mat.eqzone .= s.eqzone; mat.Λext .= s.Λext; mat.Tref .= s.Tref; mat.slope .= s.slope
    HOTON && hasproperty(s, :Y) && s.Y !== nothing && (SG[].Y .= s.Y)      # (a checkpoint from before HOT=1 keeps the initial abundances)
    ALIFAC[] = copy(aux.ALIFAC); ALILAST[] = aux.ALILAST; ALIINT[] = aux.ALIINT; NALICALL[] = aux.NALICALL; NIMP[] = aux.NIMP
    LLINE[] = aux.LLINE; nextln[] = aux.nextln; SNAPNEXT[] = aux.SNAPNEXT
    SNAPS === nothing || aux.SNAPS === nothing || (empty!(SNAPS); append!(SNAPS, aux.SNAPS))
    empty!(Trbdf2.STATS); merge!(Trbdf2.STATS, aux.trb.STATS); Trbdf2.BE_LEFT[] = aux.trb.BE_LEFT
    TX.PEXT[] = aux.PEXT; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TF.setscalar!(:tmax, 1e3tff)                   # a run stopped by TSTOP/NSTOP has tmax = its last time
    KOFF[] = aux.k
    pr("restarted from checkpoint %s at step %d, t = %.6e yr, log n_c = %.3f, T_c = %.1f K\n", name, aux.k, TF.getscalar(:timen) / 3.156e7,
       log10(TF.var(:dn)[a] / PC.mH), TF.var(:tn)[a])
end
if HOTON
    const SM = StarMaterial; const NN = NucNet; const SE = StellarEOS
    const HX = 1 / (1 + 4PC.yHe)
    const DH = parse(Float64, get(ENV, "DH", "2.5e-5"))
    SG[] = SM.StarGas(TF.MGR + 2; X = HX, Y = 1 - HX, Z = 0.0, opac = :stellar)
    let Y0 = zeros(NN.NS)                             # primordial abundances (mol/g), as S5: D/H, ³He/H; ⁴He takes the rest of Y
        Yp = HX / (1 + 2DH); Y0[NN.IDX[:p]] = Yp; Y0[NN.IDX[:d]] = DH * Yp
        Y3 = 1e-5 * Yp; Y0[NN.IDX[:he3]] = Y3; Y0[NN.IDX[:he4]] = (1 - HX - 3Y3) / 4
        for k in axes(SG[].Y, 2); SG[].Y[:, k] .= Y0; end
    end
    let g = (ρ, T, k) -> (st = SE.gas_state(ρ, T, SM.zonecomp(SG[], k)); (st.p, st.e)),
        th = (ρ, T, k) -> (st = SE.eos(ρ, T; comp = SM.zonecomp(SG[], k), radiation = true); (st.∇ad, st.cP, st.δ))
        ρr, Tr = 1e-4, 1e6                            # zero point of the two EOS in the fully ionised dilute limit
        Δ0 = g(ρr, Tr, 1)[2] - PopIIIEqEOS.eos(ρr, Tr)[1]
        PM.HOT[] = (; T1 = parse(Float64, get(ENV, "HOTT1", "3e5")), T2 = parse(Float64, get(ENV, "HOTT2", "9e5")), Δ0, gas = g, thermo = th)
    end
    pr("HOT: stellar branch above T1 = %.3g K (full at %.3g K), Delta0 = %.4e erg/g, D/H = %.2e\n", PM.HOT[].T1, PM.HOT[].T2, PM.HOT[].Δ0, DH)
    "T with the installed material's e_gas(ρ, T) = e in zone k (safeguarded Newton in ln T)."
    function T_hot(ρ, e, T0, k)
        TX.ZONE[] = k; lT = log(T0); lo = lT - 1.5; hi = lT + 1.5
        for _ in 1:60
            f = TX.material_gas(mat, ρ, exp(lT))[2] - e
            abs(f) <= 1e-12 * abs(e) && break
            f > 0 ? (hi = lT) : (lo = lT)
            d = (TX.material_gas(mat, ρ, exp(lT + 1e-6))[2] - TX.material_gas(mat, ρ, exp(lT))[2]) / 1e-6
            lN = lT - f / d; lT = (d > 0 && lo < lN < hi) ? lN : 0.5 * (lo + hi)
        end
        exp(lT)
    end
    const EBURN = Ref(0.0); const HOTEFRAC = parse(Float64, get(ENV, "HOTEFRAC", "0.02"))
    """
    S6-C1: after a converged step and chem!: in zones above T1 the network over Δt at (ρ, T) (operator split), the release into e at
    fixed ρ, T from the EOS; then convective mixing of the abundances (D from the MLT of `step_diag`, StarMaterial.mix_species!; e fixed).
    Δt ≤ HOTEFRAC e/|ε| of the zones. Total energy (gas + cumulative release EBURN) conserved to the T inversion. Good for D burning; the
    stiff CNO ignition needs the implicit source coupling of S5 (solve_step!) still.
    """
    function hotburn!(dt)
        h = PM.HOT[]; dn, tn, egn, rn, dv = TF.var(:dn), TF.var(:tn), TF.var(:egn), TF.var(:rn), TF.var(:dvoln)
        any(tn[k] > h.T1 for k in a:b) || return 0.0
        sg = SG[]; frac = 0.0; ΣE = 0.0
        for k in a:b
            tn[k] > h.T1 || continue
            res = NN.burn(sg.Y[:, k], dn[k], tn[k], dt; screen = sg.screen)
            sg.Y[:, k] .= res.Y; de = res.ε * dt; ΣE += de * 4π * dn[k] * dv[k]
            frac = max(frac, abs(de) / abs(egn[k]))
            egn[k] += de; tn[k] = T_hot(dn[k], egn[k], tn[k], k)
        end
        EBURN[] += ΣE
        sd = TX.step_diag(dt); n = b - a + 1; Di = zeros(n + 1)
        for j in 2:n; k = a + j - 1; (sd.fcv[k] > 0 && sd.vcv[k] > 0 && tn[k] > h.T1 && tn[k-1] > h.T1) && (Di[j] = sd.vcv[k] * sd.lcv[k] / 3); end
        if any(>(0), Di)
            dm = [4π * dn[a+j-1] * dv[a+j-1] for j in 1:n]; ar = zeros(n + 1)
            for j in 2:n; k = a + j - 1; ar[j] = 4π * rn[k]^2 * 0.5 * (dn[k-1] + dn[k]); end
            Yz = sg.Y[:, a:b]; SM.mix_species!(Yz, Di, dm, ar, dt); sg.Y[:, a:b] .= Yz
            for j in 1:n; k = a + j - 1; tn[k] > h.T1 && (tn[k] = T_hot(dn[k], egn[k], tn[k], k)); end     # (composition changed at fixed e)
        end
        SM.ghosts!(sg)
        frac > 0 && TF.setscalar!(:dtime, min(TF.getscalar(:dtime), dt * HOTEFRAC / frac))
        ΣE / dt
    end
end
"""
S6-B: settings of a radiating run (after a RAD rebuild and after restarting a RAD checkpoint): StellarOpacity continuum for
NetGas, CMB intensity at the outer boundary, Eddington-factor mode, thin-envelope E_r out of the Δt rule, Δt cap.
"""
const DTRHOOK = Ref(false)
function radsettings!()
    X = 1 / (1 + 4PC.yHe)
    PM.OPAC[] = (ρ, T) -> (StellarOpacity.kappa_R(ρ, T, X), StellarOpacity.kappa_P(ρ, T, X))
    Tb = parse(Float64, get(ENV, "TCMB", "0"))
    TF.setscalar!(:ximr, Z0.σSB * Tb^4 / π)
    TF.setscalar!(:ltran, parse(Int, get(ENV, "LTRAN", "1")))      # 1: formal solution for f, g; 2: Eddington approximation f = 1/3, g = 1/2
    TC.ER_TAUMIN[] = parse(Float64, get(ENV, "ERTAUMIN", "1.0"))   # E_r of the optically thin envelope does not limit Δt
    # Δt cap: steps much longer than the light-crossing time of the nearly free-streaming layers outside the photosphere
    # (r ~ 1e14 cm: r/c ~ 3e3 s) fail in Newton (Eddington factors lagged over the step); default 1500 s
    if HOTON                                          # convection (hot branch only: PM.material_thermo), as the S5 stars
        TX.CONVECTION[] = true; TX.ALPHA_MLT[] = parse(Float64, get(ENV, "ALPHAMLT", "1.5")); TX.MLT_HSDISC[] = true; TX.MLT_RCAP[] = true
    end
    dtr = parse(Float64, get(ENV, "DTRADMAX", "1500"))
    dtr > 0 && !DTRHOOK[] && (push!(TC.DT_HOOKS, () -> dtr); DTRHOOK[] = true)
    TD.MSGIO[] = get(ENV, "MSG", "0") == "1" ? stdout : devnull
end
"""
Rebuild the TITAN state from a zone profile (interfaces rI, u, Fr; zones ρ, T, Er; `matz` = per-zone material state, NamedTuple
of vectors x, eqzone, Λext, Tref, slope) with PS.build! (rad = true), keeping time, Δt, step number, the solver controls and
PEXT. Used by the RAD rebuild (S6-B) and by the rezoning (S6-C2).
"""
function rebuild!(rI, ρz, Tz, uI, Er, Fr, matz)
    t0, dt0, js = TF.getscalar(:timen), TF.getscalar(:dtime), TF.getscalar(:jstep)
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    pext = TX.PEXT[]
    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lgrid = 3, inner = :wall, outer = :free, pext = 0.0,
              selfgrav = true, material = mat, av)
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    for j in 1:nsh
        k = a + j - 1
        mat.x[k] = copy(matz.x[j]); mat.eqzone[k] = matz.eqzone[j]; mat.Λext[k] = matz.Λext[j]; mat.Tref[k] = matz.Tref[j]; mat.slope[k] = matz.slope[j]
    end
    if HOTON
        for j in 1:nsh; SG[].Y[:, a+j-1] .= matz.Y[j]; end
        SM.ghosts!(SG[])
    end
    for k in 1:a-1; mat.x[k] = copy(mat.x[a]); mat.eqzone[k] = mat.eqzone[a]; end
    for k in b+1:length(mat.x); mat.x[k] = copy(mat.x[b]); mat.eqzone[k] = mat.eqzone[b]; mat.Λext[k] = 0.0; end
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    PS.install_julia!(BlockSolver)
    TX.PEXT[] = pext
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, t0); TF.setscalar!(:timeo, t0); TF.setscalar!(:dtime, dt0); TF.setscalar!(:jstep, js); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    empty!(ALIFAC[])                                  # per-zone NLTE factors: recomputed at the next line-transfer call
end
"""
S6-B: continue a hydrodynamics-only checkpoint with TITAN's radiation. The checkpoint is restored (`restart!`), its profile
(r, ρ, T, u) rebuilt as a radiating state (`rebuild!`: the SUBC template, neqn = 7). Initial radiation field:
E_r = a [T_CMB⁴ + (T⁴ − T_CMB⁴)(1 − e^{−τ_R})] (τ_R from the surface), F_r the diffusion flux −c ΔE_r/(3χρ Δr)
flux-limited (|F| ≤ c E_r), 0 at the centre; outside the photosphere the steady free-streaming field of the photospheric
luminosity on top of the CMB.
"""
function radbuild!(name)
    restart!(name; allow = ("RAD",))
    radsettings!()
    rI = TF.var(:rn)[a:b+1]; ρz = TF.var(:dn)[a:b]; Tz = TF.var(:tn)[a:b]; uI = TF.var(:un)[a:b+1]
    old = matstate(); aold = a
    n = length(ρz); Tb = parse(Float64, get(ENV, "TCMB", "0")); ar = Z0.arad
    χ = [PM.OPAC[](ρz[j], Tz[j])[1] * ρz[j] for j in 1:n]
    τ = zeros(n); for j in n-1:-1:1; τ[j] = τ[j+1] + 0.5 * (χ[j] * (rI[j+1] - rI[j]) + χ[j+1] * (rI[j+2] - rI[j+1])); end
    Er = [ar * (Tb^4 + (Tz[j]^4 - Tb^4) * (-expm1(-τ[j]))) for j in 1:n]
    rc = [0.5 * (rI[j] + rI[j+1]) for j in 1:n]
    Fr = zeros(n + 1)
    for j in 2:n
        gE = (Er[j] - Er[j-1]) / (rc[j] - rc[j-1]); Em = 0.5 * (Er[j] + Er[j-1])
        Fr[j] = -Z0.c * gE / (3 * 0.5 * (χ[j] + χ[j-1]) + abs(gE) / Em) * (-expm1(-τ[j]))   # flux-limited diffusion, switched off where thin (isotropic CMB: F = 0)
    end
    Fr[n+1] = Fr[n]
    # outside the photosphere (τ_R < 1) the steady free-streaming field of the photospheric luminosity on top of the CMB
    # (starting from the CMB alone, the core's light front sweeps the thin envelope and Newton fails for Δt ≫ r/c)
    jph = findfirst(<(1.0), τ); Lph = 0.0
    if jph !== nothing && jph > 1
        Lph = 4π * rI[jph]^2 * max(Fr[jph], 0.0)
        for j in jph:n; Er[j] = max(Er[j], ar * Tb^4 + Lph / (4π * rc[j]^2 * Z0.c)); end
        for j in jph:n+1; Fr[j] = Lph / (4π * rI[j]^2); end
    end
    zr = aold:aold+n-1
    rebuild!(rI, ρz, Tz, uI, Er, Fr, (; x = old.x[zr], eqzone = old.eqzone[zr], Λext = old.Λext[zr], Tref = old.Tref[zr], slope = old.slope[zr],
                                         Y = HOTON ? [old.Y[:, k] for k in zr] : nothing))
    pr("RAD: rebuilt with radiation from %s: ngrs %d -> %d, tau_R(centre) = %.3e, E_r/aT^4 at centre %.3f, photosphere at zone %s (r = %.3e cm), L_ph = %.3e erg/s, ximr = %.3e\n",
       name, aold, a, τ[1], Er[1] / (ar * Tz[1]^4), string(jph), jph === nothing ? NaN : rI[jph], Lph, TF.getscalar(:ximr))
end
include(joinpath(@__DIR__, "Rezone.jl"))
"Gas energy of zone k's material at (ρ, T) with composition (x, eq)."
gasu(ρ, T, x, eq) = (PM.HOT[] !== nothing && T > PM.HOT[].T1) ? TX.material_gas(mat, ρ, T)[2] : eq ? EQ_eos(ρ, T)[1] : ChemNet.eos(ρ, T, x, mat.net)[1]
EQ_eos(ρ, T) = PopIIIEqEOS.eos(ρ, T)
"T with gasu(ρ, T) = e (bracketed Newton in ln T)."
function T_of_e(ρ, e, T0, x, eq)
    lo, hi = log(T0) - 3, log(T0) + 3; lT = log(T0)
    for _ in 1:100
        f = gasu(ρ, exp(lT), x, eq) - e
        abs(f) <= 1e-12 * abs(e) && break
        f > 0 ? (hi = lT) : (lo = lT)
        h = 1e-6; d = (gasu(ρ, exp(lT + h), x, eq) - gasu(ρ, exp(lT), x, eq)) / h
        lT1 = lT - f / d
        lT = (d > 0 && lo < lT1 < hi) ? lT1 : 0.5 * (lo + hi)
    end
    exp(lT)
end
"Current zone profile for Rezone (E_r, F_r from the radiating state)."
function profile_now()
    rr = TF.var(:rn)[a:b+1]; rr[1] = 0.0
    Rezone.Profile(rr, TF.var(:un)[a:b+1], TF.var(:frn)[a:b+1], TF.var(:dn)[a:b], TF.var(:egn)[a:b], TF.var(:tn)[a:b], TF.var(:ern)[a:b], TF.var(:asn)[a:b], TF.var(:pgn)[a:b])
end
const REZ = (on = get(ENV, "REZONE", "0") == "1", every = parse(Int, get(ENV, "REZEVERY", "20")),
             jsplit = parse(Float64, get(ENV, "REZJSPLIT", "0.3")), jmerge = parse(Float64, get(ENV, "REZJMERGE", "0.05")),
             dmmin = parse(Float64, get(ENV, "REZDMMIN", "1e-7")) * Msun, qmax = parse(Float64, get(ENV, "REZQMAX", "2.5")),
             nmax = parse(Int, get(ENV, "REZNMAX", "260")), nhard = TF.MGR - 12, fmin = parse(Float64, get(ENV, "REZFRAC", "0.01")), keepcentre = parse(Int, get(ENV, "REZKEEPC", "5")))
PM.PSWHYS[] = parse(Float64, get(ENV, "PSWHYS", "1"))
const REZW = haskey(ENV, "REZWIN") ? parse.(Float64, split(ENV["REZWIN"], ",")) : nothing; const NREZ = Ref(0); const REZNULL = get(ENV, "REZNULL", "0") == "1"
"""
S6-C2: rezone (Rezone.plan/apply) and rebuild. The material state of a split zone is copied; of a merged zone, mass-weighted
(x), equilibrium flag if either was, Λext volume-weighted (it is per volume), Tref/slope of the heavier zone. Returns the
numbers of splits and merges.
"""
function rezone!(k; force = nothing)
    p = profile_now()
    # moving window (REZWIN=flo,fhi,frac): zones with m_enc in [flo, fhi] × m★ (m★: inside the first interface falling faster than
    # 2 km/s) and heavier than frac × m_enc are split while they are still smooth, ahead of the star
    fl = nothing
    if REZW !== nothing
        dm_ = Rezone.dmass(p); me_ = cumsum(dm_); jc_ = findfirst(<(-2e5), p.u); mst = jc_ === nothing ? 0.0 : me_[max(jc_ - 1, 1)]
        fl = [REZW[1] * mst <= me_[j] <= REZW[2] * mst && dm_[j] > REZW[3] * me_[j] for j in eachindex(dm_)]
    end
    sp, mg = REZNULL ? (Int[], Int[]) : force === nothing ? Rezone.plan(p; REZ.jsplit, REZ.jmerge, REZ.dmmin, REZ.qmax, REZ.nmax, REZ.nhard, REZ.keepcentre, REZ.fmin, flag = fl) : force
    (isempty(sp) && isempty(mg)) && !REZNULL && return (0, 0)    # REZNULL=1: rebuild from the unchanged profile (test of the rebuild alone)
    if get(ENV, "REZDBG", "0") == "1"
        J = Rezone.jumps(p); dm = Rezone.dmass(p); me = cumsum(dm)
        pr("    splits (zone, m_enc, dm/m_enc, J): %s\n    merges: %s\n", join(("($i, $(round(me[i]/Msun, sigdigits=3)), $(round(dm[i]/me[i], sigdigits=2)), $(round(J[i], sigdigits=2)))" for i in sp), " "),
           join(("($i, $(round(me[i]/Msun, sigdigits=3)))" for i in mg), " "))
    end
    old = matstate(); zx(i) = a + i - 1
    Tofe(ρ, e, T0, srcj) = (TX.ZONE[] = zx(srcj[1][1]); eq = any(old.eqzone[zx(i)] for (i, _) in srcj);
                           x = sum(w .* old.x[zx(i)] for (i, w) in srcj); T_of_e(ρ, e, T0, x, eq))
    function eosT(ρ, T, srcj)                         # (e, P) of the material with the composition of srcj
        TX.ZONE[] = zx(srcj[1][1]); eq = any(old.eqzone[zx(i)] for (i, _) in srcj); x = sum(w .* old.x[zx(i)] for (i, w) in srcj)
        if PM.HOT[] !== nothing && T > PM.HOT[].T1; p_, u_ = TX.material_gas(mat, ρ, T); return (u_, p_); end
        u_, p_ = eq ? PopIIIEqEOS.eos(ρ, T) : ChemNet.eos(ρ, T, x, mat.net); (u_, p_)
    end
    q, src = Rezone.apply(p, sp, mg; Tofe, eosT = get(ENV, "REZPSPLIT", "1") == "1" ? eosT : nothing)
    # conservation over the range of zones that changed (whole-grid totals are dominated by the 1e4 M☉ envelope)
    ch = [j for (j, s_) in enumerate(src) if length(s_) > 1 || s_[1][1] in sp]
    isempty(ch) && (ch = [1, length(q.ρ)])
    lo_n, hi_n = max(first(ch) - 1, 1), min(last(ch) + 1, length(q.ρ))
    lo_o, hi_o = src[lo_n][1][1], src[hi_n][end][1]
    t0 = Rezone.totals(p, lo_o, hi_o); t1 = Rezone.totals(q, lo_n, hi_n)
    V(rI, j) = 4π / 3 * (rI[j+1]^3 - rI[j]^3)
    matz = (; x = [sum(w .* old.x[zx(i)] for (i, w) in s) for s in src],
              eqzone = [any(old.eqzone[zx(i)] for (i, _) in s) for s in src],
              Λext = [sum(old.Λext[zx(i)] * V(p.rI, i) for (i, _) in s) / V(q.rI, j) for (j, s) in enumerate(src)],
              Tref = [old.Tref[zx(s[argmax([w for (_, w) in s])][1])] for s in src],
              slope = [old.slope[zx(s[argmax([w for (_, w) in s])][1])] for s in src],
              Y = HOTON ? [sum(w .* old.Y[:, zx(i)] for (i, w) in s) for s in src] : nothing)
    nold = nsh
    rebuild!(q.rI, q.ρ, q.T, q.u, q.Er, q.Fr, matz)
    # the energies after the rebuild (e from the EOS at the new T): the T inversion's error is part of the check
    pb = profile_now(); t2 = Rezone.totals(pb, lo_n, hi_n)
    NREZ[] += 1
    pr("  rezone at step %d: %d splits, %d merges, zones %d -> %d; zones %d-%d: dM/M %.1e, dE/(Eint+Ekin+Erad) %.1e (after rebuild %.1e), dP/P %.1e\n",
       k, length(sp), length(mg), nold, nsh, lo_o, hi_o, (t1.M - t0.M) / t0.M,
       (t1.Eint + t1.Ekin + t1.Erad - t0.Eint - t0.Ekin - t0.Erad) / (abs(t0.Eint) + t0.Ekin + t0.Erad),
       (t2.Eint + t2.Ekin + t2.Erad - t0.Eint - t0.Ekin - t0.Erad) / (abs(t0.Eint) + t0.Ekin + t0.Erad), (t2.P - t0.P) / abs(t0.P))
    (length(sp), length(mg))
end
"""
S6-C2: static pre-refinement while the gas is still smooth (REFINE=mlo,mhi,frac [M☉, M☉, Δm/m_enc]): zones with enclosed mass in
[mlo, mhi] heavier than frac × m_enc are split, repeatedly, with the REZQMAX grading at the edges. Refining an already formed core
makes it ring (its discrete hydrostatic balance belongs to the coarse grid); a core formed on the fine grid does not need it.
"""
function prerefine!(k, mlo, mhi, frac)
    for pass in 1:8
        p = profile_now(); dm = Rezone.dmass(p); me = cumsum(dm)
        fl = [mlo * Msun <= me[j] <= mhi * Msun && dm[j] > frac * me[j] for j in eachindex(dm)]
        any(fl) || break
        sp, mg = Rezone.plan(p; jsplit = Inf, REZ.qmax, nhard = REZ.nhard, flag = fl)
        rezone!(k; force = (sp, mg))
    end
    pr("REFINE: zones %d after pre-refinement of %g-%g Msun to dm/m <= %g\n", nsh, mlo, mhi, frac)
end
if haskey(ENV, "RADFROM")
    RAD || error("RADFROM needs RAD=1")
    radbuild!(ENV["RADFROM"])
    haskey(ENV, "REFINE") && prerefine!(KOFF[], parse.(Float64, split(ENV["REFINE"], ","))...)
elseif haskey(ENV, "RESTART")
    restart!(ENV["RESTART"])
    RAD && Int(TF.getscalar(:lrad)) == 1 && (radsettings!(); TD.EOS[](); TD.OPAC[]())
end
const CHEMDIAG = get(ENV, "CHEMDIAG", "0") == "1"; const LNUC = Ref(0.0); const TIMING = get(ENV, "TIMING", "0") == "1"; const TSTEP = Ref(0.0); const TCHEM = Ref(0.0); const TLINE = Ref(0.0); const TOTH = Ref(0.0)
const TMARK = Ref(time()); const TOUT = Ref(0.0)
const PRINTEVERY = parse(Int, get(ENV, "PRINTEVERY", "500")); const CKEVERY = parse(Int, get(ENV, "CKEVERY", "500")); const KLAST = Ref(0)
try
    TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = function (k0)
        k = k0 + KOFF[]; KLAST[] = k
        tnow = time(); TSTEP[] += tnow - TMARK[] - TOUT[]; TOUT[] = 0.0
        tc0 = time(); Tb4 = CHEMDIAG ? copy(TF.var(:tn)) : nothing; ub4 = CHEMDIAG ? copy(TF.var(:egn)) : nothing; dT = PM.chem!(mat)
        if CHEMDIAG && dT > 0.01                      # the zones chem! moved most
            tt_ = TF.var(:tn); jj = sort(collect(a:b); by = kk -> -abs(tt_[kk] / Tb4[kk] - 1))[1:3]
            for kk in jj
                pr("        chem! zone %d: T %.1f -> %.1f, m %.4e Msun, n %.3e, p %.3e, eq %s, e_before %.5e e_eos_after %.5e, dt %.2e\n", kk - a + 1, Tb4[kk], tt_[kk], 4π * TX.XMINT[][kk+1] / Msun,
                   TF.var(:dn)[kk] / (1.22PC.mH), TF.var(:pgn)[kk], string(mat.eqzone[kk]), ub4[kk], TF.var(:egn)[kk], TF.getscalar(:timen) - TF.getscalar(:timeo))
            end
        end
        HOTON && (LNUC[] = hotburn!(TF.getscalar(:timen) - TF.getscalar(:timeo)))
        TCHEM[] += time() - tc0
        RAD && TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())   # T changed in chem!: Planck function, opacity, e_g + E_r/ρ
        if SNAPS !== nothing                         # profiles for post-processing: every 0.25 dex of n_c, and every SNAPEVERY steps
            lnc = log10(TF.var(:dn)[a] / PC.mH)
            if lnc >= SNAPNEXT[] || (SNAPEVERY > 0 && k % SNAPEVERY == 0)
                push!(SNAPS, (; step = k, t = TF.getscalar(:timen), r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1],
                    ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b], y = [ycomp(a + j - 1) for j in 1:nsh],
                    m = 4π .* TX.XMINT[][a:b+1], Λext = copy(mat.Λext[a:b]), eq = [eqflag(a + j - 1) for j in 1:nsh]))
                lnc >= SNAPNEXT[] && (SNAPNEXT[] = lnc + 0.25)
                serialize(SNAPFILE, SNAPS)
            end
        end
        if cool == "lines"
            # H₂ line transfer on the new state (RHFC02 A2.4); Λ held fixed over the next step; Δt ≤ 0.1 u ρ/|Λ| (A2.2)
            rr = TF.var(:rn)[a:b+1]; uu = TF.var(:un)[a:b+1]; dd = TF.var(:dn)[a:b]; tt = TF.var(:tn)[a:b]
            nh2 = [dd[j] / PC.μH * ycomp(a + j - 1)[3] for j in 1:nsh]
            tl0 = time(); Λl, Lr, _ = linecool(rr, uu, dd, tt, nh2, k); TLINE[] += time() - tl0
            tcmin = Inf
            for j in 1:nsh
                mat.Λext[a+j-1] = Λl[j]; mat.Tref[a+j-1] = tt[j]; mat.slope[a+j-1] = 0.0
                ue = TF.var(:egn)[a+j-1] + PC.UOFF_TH                 # thermal-ish energy scale (positive)
                Λl[j] != 0 && abs(ue * dd[j] / Λl[j]) < tcmin && (tcmin = abs(ue * dd[j] / Λl[j]); JCMIN[] = j)
            end
            # G7 (IMPCOOL=1, default): where the explicit cap Δt ≤ 0.1 t_cool would bind, the slope d ln Λ/d ln T from a
            # second transfer call with all T × (1 + ε) makes Λ_j (T/T_j)^s implicit in the chemistry step, and the cap
            # relaxes to Δt ≤ DTCOOL t_cool (default 1); otherwise Λ is held fixed (explicit, accurate there)
            if IMPCOOL && (IMPALWAYS || 0.1tcmin < TF.getscalar(:dtime))
                ε = 0.02; tl0 = time(); Λp, _, _ = linecool(rr, uu, dd, tt .* (1 + ε), nh2, -1); TLINE[] += time() - tl0
                for j in 1:nsh
                    Λl[j] != 0 && Λp[j] / Λl[j] > 0 && (mat.slope[a+j-1] = clamp(log(Λp[j] / Λl[j]) / log(1 + ε), -10.0, 30.0))
                end
                NIMP[] += 1
            end
            LLINE[] = Lr[end]
            if get(ENV, "DTDIAG", "0") == "1"
                j_ = JCMIN[]; pr("        line cap: t_cool %.3e s at zone %d (m %.4e Msun, T %.1f, n %.3e, e+U %.3e, Lambda %.3e); dtime before %.3e, smax %.3e\n", tcmin, j_, 4π * TX.XMINT[][a+j_] / Msun, tt[j_], dd[j_] / (1.22PC.mH),
                                TF.var(:egn)[a+j_-1] + PC.UOFF_TH, Λl[j_], TF.getscalar(:dtime), TF.getscalar(:smax))
            end
            TF.setscalar!(:dtime, min(TF.getscalar(:dtime), (IMPCOOL ? DTCOOL : 0.1) * tcmin))
        end
        if RAD && haskey(ENV, "ERDIAG") && k in parse.(Int, split(ENV["ERDIAG"], ","))     # radiation field across the grid
            ern = TF.var(:ern); frn = TF.var(:frn); fe = TF.var(:fedd); ch = TF.var(:chifn); dd = TF.var(:dn); tt = TF.var(:tn); rr = TF.var(:rn)
            Tb = parse(Float64, get(ENV, "TCMB", "0")); τ = 0.0; τs = zeros(b)
            for kk in b:-1:a; τ += ch[kk] * dd[kk] * (rr[kk+1] - rr[kk]); τs[kk] = τ; end
            pr("  radiation at step %d (geddr %.3f):\n", k, TF.getscalar(:geddr))
            for kk in a:6:b
                pr("    z %3d m %.3e r %.3e tau %.3e dtau %.2e E/aT4 %.3e E/aTcmb4 %.3e F/cE %+.3e f %.3f T %.1f\n", kk - a + 1, 4π * TX.XMINT[][kk+1] / Msun, rr[kk+1],
                   τs[kk], ch[kk] * dd[kk] * (rr[kk+1] - rr[kk]), ern[kk] / (Z0.arad * tt[kk]^4), ern[kk] / (Z0.arad * Tb^4), frn[kk+1] / (Z0.c * ern[kk]), fe[kk], tt[kk])
            end
        end
        if haskey(ENV, "PROFAT") && k in parse.(Int, split(ENV["PROFAT"], ","))
            rr = TF.var(:rn); uu = TF.var(:un); dd = TF.var(:dn); tt = TF.var(:tn)
            pr("  profile at step %d, t = %.4e yr:\n", k, TF.getscalar(:timen) / 3.156e7)
            for j in 1:8:nsh
                kk = a + j - 1
                pr("    j %3d m %.3e r %.3e pc  n %.3e  T %7.1f  u %7.3f km/s  f_H2 %.2e\n", j, 4π * TX.XMINT[][kk+1] / Msun, rr[kk+1] / 3.0857e18,
                   dd[kk] / (1.22 * PC.mH), tt[kk], uu[kk+1] / 1e5, 2ycomp(kk)[3])
            end
        end
        ρc_ = TF.var(:dn)[a]; Tc_ = TF.var(:tn)[a]; ln = log10(ρc_ / PC.mH)
        if ln >= nextln[] || k % PRINTEVERY == 0
            umin = minimum(TF.var(:un)[a:b+1]) / 1e5
            pr("%7d %8.3f %10.4e %9.1f %9.1f %10.3e %10.3f %9.2e %9.2e%s\n", k, ln, TF.getscalar(:timen) / 3.156e7, Tc_, on98T(ln),
               2ycomp(a)[3], umin, dT, TF.getscalar(:dtime), cool == "lines" ? @sprintf("  L_lines %.3e  bins %d, thick lines %d", LLINE[], H2Lines.LAST[].nbins, H2Lines.LAST[].nthick) : "")
            ln >= nextln[] && (nextln[] += 0.25)
            if get(ENV, "UDIAG", "0") == "1"            # where the fastest infall is
                uu = TF.var(:un)[a:b+1]; ku = argmin(uu); rr = TF.var(:rn)[a:b+1]
                pr("        u_min %.3f km/s at interface %d of %d, r = %.3e pc, m = %.3e Msun, T(zone) = %.1f, n = %.3e\n", uu[ku] / 1e5, ku, nsh + 1,
                   rr[ku] / 3.0857e18, 4π * TX.XMINT[][a+ku-1] / Msun, TF.var(:tn)[a+min(ku, nsh)-1], TF.var(:dn)[a+min(ku, nsh)-1] / (1.22 * PC.mH))
            end
            if get(ENV, "DTDIAG", "0") == "1"           # which variable/zone limits Δt (timstp: kx, dx per equation)
                kx = TF.var(:kx); dx = TF.var(:dx); names = ("r", "m", "rho", "u", "T", "E_r", "F_r"); ne_ = RAD ? 6 : 5   # (F_r is not in smax)
                j = argmax(abs.(dx[[1, 2, 3, 5, 6][1:ne_-1]])); j = (1, 2, 3, 5, 6)[j]; xs = TF.getscalar(:stol) / TF.getscalar(:smax)
                pr("        dt limit: %s in zone %d (rel. change %.2e, stol/smax = %.2f)%s\n", names[j], kx[j], dx[j], xs,
                   RAD ? @sprintf("; E_r: zone %d %.2e, rho: zone %d %.2e", kx[6], dx[6], kx[3], dx[3]) : "")
            end
        end
        (Tc_ > Tstop || ρc_ / PC.mH > nstop) && TF.setscalar!(:tmax, TF.getscalar(:timen))
        REZ.on && RAD && k % REZ.every == 0 && rezone!(k)
        if TIMING && k % PRINTEVERY == 0
            pr("        timing since last report: step (hydro+radiation) %.1f s, chemistry %.1f s, H2 lines %.1f s, other output %.1f s\n",
               TSTEP[], TCHEM[], TLINE[], time() - tnow - 0.0 + TOTH[] - TCHEM[] - TLINE[])
            TSTEP[] = TCHEM[] = TLINE[] = TOTH[] = 0.0
        else
            TOTH[] += time() - tnow
        end
        TMARK[] = time()
        haskey(ENV, "SAVE") && CKEVERY > 0 && k % CKEVERY == 0 && checkpoint(ENV["SAVE"] * "_auto" * string((k ÷ CKEVERY) % 2), k)
    end)
catch e
    pr("stopped: %s\n", sprint(showerror, e)[1:min(end, 400)])
    get(ENV, "DEBUG", "0") == "1" && Base.display_error(stderr, e, catch_backtrace())
end
haskey(ENV, "SAVE") && (checkpoint(ENV["SAVE"], KLAST[]); pr("checkpoint %s at step %d\n", ENV["SAVE"], KLAST[]))
IMPCOOL && pr("implicit line cooling on %d steps\n", NIMP[])
NLTE == "hyb" && pr("hybrid NLTE: %d ALI calls\n", NALICALL[])
TD.STEPFN[] === Trbdf2.step! && pr("integrator TR-BDF2, statistics %s\n", string(Trbdf2.STATS))
close(io)
