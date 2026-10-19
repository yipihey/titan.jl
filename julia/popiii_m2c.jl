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
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const PC = PrimordialChem; const PM = PopIIIMaterial; const Z0 = Z0Physics
using ForwardDiff
# INTEGRATOR=trbdf2: second-order L-stable TR-BDF2 (Trbdf2.jl, G2) instead of TITAN's theta scheme
include(joinpath(@__DIR__, "Trbdf2.jl"))
get(ENV, "INTEGRATOR", "theta") == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
const Msun = 1.98847e33

cool = isempty(ARGS) ? "ra04" : ARGS[1]; nsh = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100
cf = cool == "thin" ? PC.cool_H2_thin : cool == "ra04lines" ? PC.cool_H2_ra04 :
     cool == "lines" ? PC.cool_CIE_ra04 :                     # H₂ lines by transfer (H2Lines), continuum: RA04 CIE
     (ρ, T, y) -> PC.cool_H2_ra04(ρ, T, y) + PC.cool_CIE_ra04(ρ, T, y)
haskey(ENV, "TCMB") && (cf = PC.cool_floor(cf, parse(Float64, ENV["TCMB"])))
io = open(joinpath(TF.ROOT, "lib", "popiii", "m2c_$(get(ENV, "IC", "p100") == "p100" ? "" : get(ENV, "IC", "") * "_")$(cool)_n$(nsh)$(get(ENV, "TAG", "")).log"), "w")
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
            rj = pc * 10^parse(Float64, get(ENV, "LRJOIN", "0.87"))
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
const IMPCOOL = get(ENV, "IMPCOOL", "1") == "1"                 # G7: transfer cooling implicit in T (linearised in ln T)
const DTCOOL = parse(Float64, get(ENV, "DTCOOL", "1.0"))
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
# and EP × (ALI/EP per zone, from the last ALI) in between
const NLTE = get(ENV, "NLTE", "1"); const NALI = parse(Int, get(ENV, "NALI", "20")); const ALIFAC = Ref{Vector{Float64}}(Float64[])
NLTE in ("ep", "ali", "hyb") && include(joinpath(@__DIR__, "H2Levels.jl"))
const HLV = NLTE in ("ep", "ali", "hyb") ? H2Levels.load() : nothing
"Line cooling per zone (and Lr) for zone temperatures tt; step k (for NLTE=hyb)."
function linecool(rr, uu, dd, tt, nh2, k)
    if NLTE in ("ep", "ali", "hyb")
        yc = [ycomp(a + j - 1) for j in 1:nsh]; nHn = dd ./ PC.μH
        nl = (LV = HLV, nH = [PC.xH(yc[j]) * nHn[j] for j in 1:nsh], nHe = PC.yHe .* nHn, nHp = [yc[j][1] * nHn[j] for j in 1:nsh])
        NLTE != "hyb" && return H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = NLTE == "ali" ? :ali : :ep)
        Λe, Lr, x = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ep)
        if k % NALI == 0 || length(ALIFAC[]) != nsh
            Λa, Lra, xa = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ali)
            ALIFAC[] = [Λe[j] != 0 && Λa[j] / Λe[j] > 0 ? Λa[j] / Λe[j] : 1.0 for j in 1:nsh]
            return Λa, Lra, xa
        end
        return Λe .* ALIFAC[], Lr, x
    end
    # NLTE=1: source function × f = 1/(1 + Λ_LTE/(n_H Λ_0)) (the GP98/HM79 combination of PrimordialChem.cool_H2_thin)
    fnl = [1 / (1 + PC.Λ_H2_LTE_HM(tt[j]) / max(PC.xH(ycomp(a + j - 1)) * dd[j] / PC.μH * PC.Λ_H2_0_GP(tt[j]), 1e-300)) for j in 1:nsh]
    H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, fnlte = NLTE == "1" ? fnl : nothing)
end
const HLL = cool == "lines" ? H2Lines.load() : nothing; const NBIN = parse(Int, get(ENV, "NBIN", "30")); const LLINE = Ref(0.0)
Tstop = parse(Float64, get(ENV, "TSTOP", "6000")); nstop = parse(Float64, get(ENV, "NSTOP", "1e18"))
try
    TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = function (k)
        dT = PM.chem!(mat)
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
            Λl, Lr, _ = linecool(rr, uu, dd, tt, nh2, k)
            tcmin = Inf
            for j in 1:nsh
                mat.Λext[a+j-1] = Λl[j]; mat.Tref[a+j-1] = tt[j]; mat.slope[a+j-1] = 0.0
                ue = TF.var(:egn)[a+j-1] + PC.UOFF_TH                 # thermal-ish energy scale (positive)
                Λl[j] != 0 && (tcmin = min(tcmin, abs(ue * dd[j] / Λl[j])))
            end
            # G7 (IMPCOOL=1, default): where the explicit cap Δt ≤ 0.1 t_cool would bind, the slope d ln Λ/d ln T from a
            # second transfer call with all T × (1 + ε) makes Λ_j (T/T_j)^s implicit in the chemistry step, and the cap
            # relaxes to Δt ≤ DTCOOL t_cool (default 1); otherwise Λ is held fixed (explicit, accurate there)
            if IMPCOOL && 0.1tcmin < TF.getscalar(:dtime)
                ε = 0.02; Λp, _, _ = linecool(rr, uu, dd, tt .* (1 + ε), nh2, -1)
                for j in 1:nsh
                    Λl[j] != 0 && Λp[j] / Λl[j] > 0 && (mat.slope[a+j-1] = clamp(log(Λp[j] / Λl[j]) / log(1 + ε), -10.0, 30.0))
                end
                NIMP[] += 1
            end
            LLINE[] = Lr[end]
            TF.setscalar!(:dtime, min(TF.getscalar(:dtime), (IMPCOOL ? DTCOOL : 0.1) * tcmin))
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
        if ln >= nextln[] || k % 500 == 0
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
                kx = TF.var(:kx); dx = TF.var(:dx); names = ("r", "m", "rho", "u", "T")
                j = argmax(abs.(dx[1:5])); xs = TF.getscalar(:stol) / TF.getscalar(:smax)
                pr("        dt limit: %s in zone %d (rel. change %.2e, stol/smax = %.2f)\n", names[j], kx[j], dx[j], xs)
            end
        end
        (Tc_ > Tstop || ρc_ / PC.mH > nstop) && TF.setscalar!(:tmax, TF.getscalar(:timen))
    end)
catch e
    pr("stopped: %s\n", sprint(showerror, e)[1:min(end, 400)])
    get(ENV, "DEBUG", "0") == "1" && Base.display_error(stderr, e, catch_backtrace())
end
IMPCOOL && pr("implicit line cooling on %d steps\n", NIMP[])
TD.STEPFN[] === Trbdf2.step! && pr("integrator TR-BDF2, statistics %s\n", string(Trbdf2.STATS))
close(io)
