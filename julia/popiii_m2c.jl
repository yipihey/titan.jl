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
# LINEGPU=cuda|metal|cpu: formal solution of the thick H2 lines through KernelAbstractions (H2LinesKA.jl); needs an
# environment with KernelAbstractions and CUDA or Metal. Off unless set.
haskey(ENV, "LINEGPU") && include(joinpath(@__DIR__, "h2lines_gpu.jl"))
include(joinpath(@__DIR__, "StellarOpacity.jl")); include(joinpath(@__DIR__, "IonRay.jl"))
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
# POPIII_OUT: directory for the run log and the snapshot file (default lib/popiii); checkpoints follow CKDIR
const OUTDIR = get(ENV, "POPIII_OUT", joinpath(TF.ROOT, "lib", "popiii")); mkpath(OUTDIR)
io = open(joinpath(OUTDIR, "m2c_$(get(ENV, "IC", "p100") == "p100" ? "" : get(ENV, "IC", "") * "_")$(cool)_n$(nsh)$(get(ENV, "TAG", "")).log"), haskey(ENV, "RESTART") || haskey(ENV, "RADFROM") ? "a" : "w")   # (a continued run appends)
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
const IMPALWAYS = get(ENV, "IMPCOOL", "1") == "2" || get(ENV, "IMPALWAYS", "0") == "1"   # IMPCOOL=2 (or IMPALWAYS=1, which a restart may change): slope every step
# TCOOLTH=1: the line-cooling time of the Δt rules below from the thermal energy (3/2 p/ρ) instead of e + UOFF_TH. The latter
# is the energy above cold H₂, which in atomic gas is ~12× the thermal energy (the H₂ binding energy, 2.24 eV per H), so
# while the gas turns molecular the explicit Λ was used up to Δt ≈ 1.2 thermal cooling times: beyond ~0.4 of it (2/(d ln Λ/
# d ln T)) T overshoots from step to step and the core rings. Default off (the accretion production runs are unchanged).
const TCOOLTH = get(ENV, "TCOOLTH", "0") == "1"
const NIMP = Ref(0)                                                  # steps with implicit line cooling
const CHEM = get(ENV, "CHEM", "rhfc02")     # rhfc02: PrimordialGas (RHFC02 network, ρ_eq switch); net / netmodern / netmodernli: NetGas (G5)
if CHEM == "rhfc02"
    global mat = PM.PrimordialGas(ystore, PC.rates_rhfc02(; k5 = k5f), cf, parse(Float64, get(ENV, "RHOEQ", "5e-11")), get(ENV, "EQEOS", "scvh") == "scvh")
else
    # G5: ChemNet network; equilibrium table where the ideal-gas pressure ≥ PSWITCH (EOS validity, RHFC02's 1e4 dyn cm⁻²)
    netw, phot = CHEM == "net" ? (ChemRates.rhfc02(; k5 = k5f), ChemNet.Photo[]) : ChemRates.modern(; li = CHEM == "netmodernli")
    # SL3=1 (with STARLIGHT=1): three-level H atom in the lit network zones (PM.tl3_reactions; rates are zero elsewhere)
    get(ENV, "SL3", "0") == "1" && (netw = PM.tl3_network(netw))
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
# SHOCKRAD=1: sub-grid radiative shock (TitanResidFast.SHOCK_RAD); SHOCKDU="du1,du2" in km/s
TX.SHOCK_RAD[] = get(ENV, "SHOCKRAD", "0") == "1"
TX.AV_TENSION[] = get(ENV, "AVTENSION", "1") == "1"       # 1 (TITAN) | 0 | f in (0,1) | radial | radial,f  (TX.AV_TFRAC)
let p = split(get(ENV, "AVTENSION", "1"), ","), f = tryparse(Float64, p[end])
    TX.AV_TRAD[] = p[1] == "radial"; TX.AV_TFRAC[] = f !== nothing && 0 < f < 1 ? f : 0.0
end
TX.AV_CAP[] = parse(Float64, get(ENV, "AVCAP", "0"))
haskey(ENV, "SHOCKDU") && (TX.SHOCK_DU[] = Tuple(1e5 .* parse.(Float64, split(ENV["SHOCKDU"], ","))))
TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))    # G1: viscous length = AVZONE × zone width (0 = ql1 r̄ above)
TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
PS.install_julia!(BlockSolver); TD.MSGIO[] = get(ENV, "MSG", "0") == "1" ? stdout : devnull   # (Newton messages also before RAD switches on)
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
const SNAPFILE = joinpath(OUTDIR, "snaps_$(get(ENV, "IC", "p100"))_$(cool)_n$(nsh)$(get(ENV, "TAG", "")).jls")
# RUNREC=<dir>: record the run for movies, interactive panels and cost accounting (off by default). Each process segment
# (start step k0, pid) writes snaps_<k0>_<pid>.jls, one serialized snapshot every SNAPEVERY steps, appended and never
# rewritten (a restart opens a new file; readers keep the surviving branch by scanning back in time), and meta_<k0>_<pid>.txt
# (the settings the code reads, git state, threads). perf.csv gets one row per PRINTEVERY steps (needs TIMING=1): the wall-time
# split, Newton/TR-BDF2 counters, ALI calls, chemistry fallbacks, GC time and allocations, the Δt limiter and the core state.
const RUNREC = get(ENV, "RUNREC", "")
isempty(RUNREC) || mkpath(RUNREC)
const REC = (on = !isempty(RUNREC), io = Ref{Any}(nothing), wall0 = time(), klast = Ref(-1), gct = Ref(Base.gc_num().total_time), gcb = Ref(Base.gc_bytes()),
             dtlim = Ref(""), mstar = Ref(NaN), rstar = Ref(NaN), spn = Ref{Vector{Symbol}}(Symbol[]))
"ENV keys any of the julia/ sources reads (only these go into the meta file)."
function recenvkeys()
    ks = Set{String}()
    for f in readdir(@__DIR__; join = true)
        endswith(f, ".jl") || continue
        for m in eachmatch(r"ENV(?:, |\[)\"([A-Z][A-Z0-9_]*)\"", read(f, String)); push!(ks, m[1]); end
    end
    ks
end
function recopen!(k)
    tag = @sprintf("%09d_%d", k, getpid())
    REC.io[] = open(joinpath(RUNREC, "snaps_$(tag).jls"), "w")
    open(joinpath(RUNREC, "meta_$(tag).txt"), "w") do io
        gitc(c) = try strip(read(Cmd(vcat(["git", "-C", TF.ROOT], c)), String)) catch; "?" end
        println(io, "start_step ", k, "\nstart_t_s ", TF.getscalar(:timen), "\nhost ", gethostname(), "\nthreads ", Threads.nthreads(), "\njulia ", VERSION)
        println(io, "git ", gitc(["rev-parse", "HEAD"]), " ", gitc(["rev-parse", "--abbrev-ref", "HEAD"]))
        println(io, "git_dirty_julia ", isempty(gitc(["status", "--short", "julia"])) ? "no" : "yes")
        ks = recenvkeys()
        for key in sort(collect(keys(ENV))); key in ks && println(io, "env ", key, "=", ENV[key]); end
    end
    pf = joinpath(RUNREC, "perf.csv")
    isfile(pf) || open(pf, "w") do io
        println(io, "step,t_yr,dt_s,wall_s,pid,nsteps,t_step,t_chem,t_lines,t_other,gc_s,alloc_GB,ndiv,nrepeat,tr_steps,fail_tr,fail_bdf,be_fallback,nali,nchemfb,nthreads,nz,log_nc,Tc,mstar,rstar,umin_kms,dtlim")
    end
end
"One snapshot into the RUNREC stream: everything needed to redraw any profile panel offline."
function recsnap!(k)
    REC.io[] === nothing && recopen!(k)
    d = [ChemNet.full(mat.x[j], mat.net) for j in a:b]
    isempty(REC.spn[]) && (REC.spn[] = sort(collect(keys(d[1]))))
    sp = REC.spn[]
    sd = PM.SLDIAG[]
    snap = (; step = k, t = TF.getscalar(:timen), dt = TF.getscalar(:dtime), a, b,
            r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1], m = 4π .* (ENVEUL === nothing ? TX.XMINT[][a:b+1] : envmenc()),
            ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b], pg = TF.var(:pgn)[a:b],
            er = TF.var(:ern)[a:b], fr = TF.var(:frn)[a:b+1], chi = TF.var(:chifn)[a:b], kp = TF.var(:xkpn)[a:b],
            fedd = TF.var(:fedd)[a:b], qf = TF.var(:qf)[a:b], eq = [eqflag(j) for j in a:b], Λext = copy(mat.Λext[a:b]),
            spn = sp, Y = Float64[Float64(ForwardDiff.value(d[j][s])) for s in sp, j in eachindex(d)],
            gasq = copy(TX.GASQ[]), hotY = HOTON ? SG[].Y[:, a:b] : nothing, Lnuc = LNUC[], Llines = LLINE[],
            nchemfb = PM.NCHEMFB[], kph = sd === nothing ? 0 : sd.kph, Teff = sd === nothing ? NaN : Float64(sd.Te))
    mat isa PM.NetGas && PM.HDCHEM[] && (snap = merge(snap, (; hd = PM.hd_snap(mat, a, b))))     # HDCHEM=1: x_D+, x_HD, HD cooling
    serialize(REC.io[], snap); flush(REC.io[])
end
"One perf.csv row (called at the TIMING report, before the accumulators are reset)."
function recperf!(k, tnow)
    REC.io[] === nothing && recopen!(k)
    gct = Base.gc_num().total_time; gcb = Base.gc_bytes()
    st = Trbdf2.STATS; un_ = TF.var(:un)[a:b+1]
    oth = time() - tnow + TOTH[] - TCHEM[] - TLINE[]
    ns = REC.klast[] < 0 ? 0 : k - REC.klast[]
    open(joinpath(RUNREC, "perf.csv"), "a") do io
        @printf(io, "%d,%.8e,%.4e,%.1f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%.6e,%.6g,%.6g,%.3f,%s\n",
                k, TF.getscalar(:timen) / 3.156e7, TF.getscalar(:dtime), time() - REC.wall0, getpid(), ns, TSTEP[], TCHEM[], TLINE[], oth,
                (gct - REC.gct[]) / 1e9, (gcb - REC.gcb[]) / 2^30, TD.NDIV[], TD.NREP[], get(st, :steps, 0), get(st, :fail_tr, 0), get(st, :fail_bdf, 0),
                get(st, :be_fallback, 0), NALICALL[], PM.NCHEMFB[], Threads.nthreads(), b - a + 1, log10(TF.var(:dn)[a] / PC.mH), TF.var(:tn)[a],
                REC.mstar[], REC.rstar[], minimum(un_) / 1e5, REC.dtlim[])
    end
    REC.klast[] = k; REC.gct[] = gct; REC.gcb[] = gcb
end
# M3: the star's ionizing light (IonRay). RAY=1: spectrum from the photosphere (τ_R = 2/3), a blackbody at T_eff from L and r there; no source below RAYTMIN.
# RAYTEST="Q,Teff" (a stress test): that many photons/s above 13.6 eV from that blackbody, whatever the star is doing.
if get(ENV, "RAY", "0") == "1" && mat isa PM.NetGas
    const RG = IonRay.groups(; species = [:H, :H2])
    const RAYTMIN = parse(Float64, get(ENV, "RAYTMIN", "1.5e4"))
    const RAYTEST = haskey(ENV, "RAYTEST") ? parse.(Float64, split(ENV["RAYTEST"], ",")) : nothing
    const RAYOUT = Ref((0, 0.0, 0.0, 0.0))                     # photosphere zone, L, T_eff, photons/s above 13.6 eV (last step)
    function raystart(a, b)
        ch, dd, rr, fr = TF.var(:chifn), TF.var(:dn), TF.var(:rn), TF.var(:frn)
        τ = 0.0; k = b
        while k > a
            τ += ch[k] * dd[k] * (rr[k+1] - rr[k]); τ >= 2 / 3 && break; k -= 1
        end
        L = 4π * rr[k+1]^2 * fr[k+1]; Te = L > 0 ? (L / (4π * rr[k+1]^2 * IonRay.σSB))^0.25 : 0.0
        if RAYTEST !== nothing
            N = IonRay.blackbody(RG, 1.0, RAYTEST[2]); N .*= RAYTEST[1] / sum(N); Te = RAYTEST[2]
        else
            Te < RAYTMIN && (RAYOUT[] = (k, L, Te, 0.0); return nothing)
            N = IonRay.blackbody(RG, L, Te)
        end
        RAYOUT[] = (k, L, Te, sum(N))
        k, N
    end
    TX.RAYQ[] = zeros(length(mat.x))
    PM.RAY[] = (; start = raystart, source = (ρ, r1, r2, Nin) -> IonRay.source(ρ, r1, r2, Nin, RG, mat.net), env = merge(mat.env, (; caseB = true)),
                Ec = RG.Ec, Q = TX.RAYQ[], xmin = parse(Float64, get(ENV, "RAYXMIN", "1e-6")))
end
# M3a: the star's light below the Lyman edge in the network zones (PM.STARLIGHT): H⁻ photodetachment, H(n ≥ 2) photoionization, Lyman-Werner
# H₂ dissociation with shielding. STARLIGHT=1; SLSHIELD = static (default) | sobolev | none (no shielding) | off (no LW); SLH=1 adds H(n ≥ 2) (energy-unsafe, see PM.STARLIGHT).
if get(ENV, "STARLIGHT", "0") == "1" && mat isa PM.NetGas
    let ph = mat.photos, pH = ph[findfirst(p -> p.R === :H, ph)]
        PM.STARLIGHT[] = (; pHm = get(ENV, "SLHM", "1") == "1" ? ph[findfirst(p -> p.R === :Hm, ph)] : nothing,
                          pH = get(ENV, "SLH", "0") == "1" ? ChemNet.Photo(pH.name * " (n ≥ 2)", :H, pH.P, pH.levels[2:end], pH.Z, false) : nothing,
                          shield = Symbol(get(ENV, "SLSHIELD", "static")), tcmb = parse(Float64, get(ENV, "TCMB", "0")),
                          tl3 = get(ENV, "SL3", "0") == "1")
    end
    (PM.STARLIGHT[].shield !== :off || PM.STARLIGHT[].tl3) && (TX.GASQ[] = zeros(length(mat.x)))
    PM.STARLIGHT[].tl3 && isempty(TX.RAYQ[]) && (TX.RAYQ[] = zeros(length(mat.x)))
    if get(ENV, "BALSPLIT", "0") == "1"
        TX.CHEMB[] = zeros(length(mat.x)); TX.CHEMBON[] = fill(false, length(mat.x))
        pr("BALSPLIT: balanced splitting of the gas energy between hydro and network chemistry (lagged mean chemistry heating in both)\n")
    end
    PM.STARLIGHT[].tl3 && PM.STARLIGHT[].pH !== nothing && error("SL3=1 replaces SLH=1 (the Boltzmann n ≥ 2 photoionization); set one")
    PM.STARLIGHT[].tl3 && pr("STARLIGHT: three-level H (Peebles C factor, Sobolev Lya escape, e collisions) in the lit network zones\n")
    PM.STARLIGHT[].tl3 && PM.LYART[] && pr("LYART: Lya transfer (Voigt, CRD core / coherent wings, ALI with the three-level atom) from tau >= %g outward\n", PM.LYATAU[])
    PM.LYART[] && !PM.STARLIGHT[].tl3 && error("LYART=1 needs STARLIGHT=1 SL3=1")
    PM.STARLIGHT[].tl3 && PM.LYABB[] && pr("LYABB: three-level H bound-bound (Lya pumping / collisional de-excitation, collisional excitation / Lya + 2-photon escape) heat in the chemistry\n")
    # NLOPAC=1 (with SL3): network zones take the non-LTE opacity of their own composition (PM.NLOPAC: Saha/Boltzmann table × lagged ratio, frozen-composition T slope)
    PM.NLOPAC[] = PM.STARLIGHT[].tl3 && get(ENV, "NLOPAC", "0") == "1"
    PM.PSWAGREE[] = PM.STARLIGHT[].tl3 ? parse(Float64, get(ENV, "PSWAGREE", "0")) : 0.0
    PM.SLBLEND[] = get(ENV, "SLBLEND", "0") == "1"
    PM.NLOPACJ[] = PM.NLOPAC[] && get(ENV, "NLOPACJ", "0") == "1"
    PM.NLOPACJ[] && pr("NLOPACJ: kappa_E of network zones weighted with the star light's B_nu(T_c)\n")
    PM.SLEB[] = PM.SLBLEND[] && get(ENV, "SLEB", "0") == "1"
    PM.SLEB[] && pr("SLEB: star-light colour temperature T_r at tau + 2/3 (Eddington-Barbier), not the photosphere's T_eff\n")
    PM.SLBLEND[] && pr("SLBLEND: star-light colour temperature T_eff (tau <= 0.3) to T_r (tau >= 3), smooth in log tau\n")
    PM.PSWSTAY[] = !(PM.STARLIGHT[].tl3 && get(ENV, "PSWSTAY", "1") == "0")
    PM.PSWSTAY[] || pr("PSWSTAY=0: below the pressure switch a table zone stays on the table only by the tau rule\n")
    PM.PSWAGREE[] > 0 && pr("PSWAGREE: below the pressure switch, table only where it agrees with the three-level network to %.3f in T\n", PM.PSWAGREE[])
    PM.NLOPAC[] && pr("NLOPAC: non-LTE opacity in the three-level network zones (no H bound-free in kappa_P; n=2 of the three-level atom in kappa_R)\n")
    pr("STARLIGHT: H- detachment%s, LW (shield %s) from the star's light in network zones\n", PM.STARLIGHT[].pH === nothing ? "" : ", H(n>=2) ionization", String(PM.STARLIGHT[].shield))
end
const SNAPS = get(ENV, "SNAPS", "0") == "1" ? Any[] : nothing
const SNAPNEXT = Ref(log10(rho[1] / PC.mH)); const SNAPEVERY = parse(Int, get(ENV, "SNAPEVERY", "0"))
# NLTE (line source functions): 1 = GP98/HM79 factor on S (old default), 0 = LTE, ep = statistical equilibrium with escape
# probabilities (H2Levels), ali = statistical equilibrium coupled to the transfer (ALI, G6c), hyb = ALI every NALI steps
# and EP × (ALI/EP per zone, from the last ALI) in between, the interval adaptive (see linecool)
const ALIWMIN = parse(Float64, get(ENV, "ALIWMIN", "0")); const ALIGROW = parse(Float64, get(ENV, "ALIGROW", "0.01")); const ALISHRINK = parse(Float64, get(ENV, "ALISHRINK", "0.05")); const NLTE = get(ENV, "NLTE", "1"); const NALI = parse(Int, get(ENV, "NALI", "20")); const ALICAP = (c_ = parse(Int, get(ENV, "ALICAP", "0")); c_ > 0 ? c_ : NALI); const ALIFAC = Ref{Vector{Float64}}(Float64[]); const ALILAST = Ref(-10^9); const ALIINT = Ref(1); const NALICALL = Ref(0)
# LINECOLL=1 (S6): the gas term of the line transfer is the net collisional exchange with the H₂ levels at the converged ALI populations
# (H2Levels.collcool; = their net radiative rate in statistical equilibrium), zone by zone, not −∇·F of the formal solution: the latter
# is the difference of the line flux across a zone, and where a large flux crosses dilute zones its discretisation error (pumped,
# fluorescent and maser lines included) was booked as gas heating or cooling up to 1e5 × the real exchange (91001: 1e-19 against
# 1e-24 erg cm⁻³ s⁻¹ at 3e7 R_sun; ΣΛV 4.5e35 against 1.7e35). Between ALI calls Λ = Λ_EP × (Λcoll/Λ_EP) as before, except where the
# two differ in sign (heated zones: escape probability can only cool): there Λcoll of the last ALI call is held (ALIFIX).
const LINECOLL = get(ENV, "LINECOLL", "0") == "1"; const ALIFIX = Ref{Vector{Float64}}(Float64[])
# LINEINC=1: the H₂ levels also see the radiation that falls on the gas from outside the lines (`H2Lines.incident`): the CMB (TCMB)
# and the continuum of the photosphere, B_ν(T_eff) from the sphere of the photosphere's radius (STARLIGHT's, of the last chemistry
# step), diluted. LINEINC=cmb: the CMB alone. Default 0: the lines see only their own radiation, so in the cold outer envelope the
# levels are excited by collisions alone and the line cooling can take the gas below T_CMB. With the incident field the levels are
# pumped as well and the gas term (LINECOLL: the net collisional exchange) becomes heating where the radiation is warmer than the
# gas. The escape-probability rate between ALI calls is then net of the absorbed background and can change sign; a zone where the
# absorbed background exceeds half its emission holds the collisional rate of the last ALI call (ALIFIX), as the heated zones do.
const LINEINC = get(ENV, "LINEINC", "0"); const LINEINCDIAG = get(ENV, "LINEINCDIAG", "0") == "1"
function lineinc()
    LINEINC in ("1", "cmb") || return nothing
    Tb = parse(Float64, get(ENV, "TCMB", "0")); sd = PM.SLDIAG[]
    (LINEINC == "cmb" || sd === nothing || !(sd.Te > 0)) && return (Tbg = Tb, Tstar = 0.0, Rph = 0.0)
    (Tbg = Tb, Tstar = Float64(sd.Te), Rph = Float64(TF.var(:rn)[sd.kph+1]))
end
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
    # the merged core is NCORE transparent zones with interfaces uniform in p² (LINECORE, S6): impact parameters through the core, so the
    # inner nodes have rays over all μ (with one core zone the inner edge saw a single ray through the centre); u linear in r inside
    nc = LINECORE; rin = rr[jc+1]
    rcz = [rin * sqrt((i - 1) / nc) for i in 1:nc]
    r2 = [rcz; rr[jc+1:end]]; u2 = [uu[jc+1] .* rcz ./ rin; uu[jc+1:end]]
    m_c = sum(dd[j] * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:jc); d2 = [fill(m_c / (4π / 3 * rin^3), nc); dd[jc+1:end]]
    t2 = [fill(tt[jc], nc); tt[jc+1:end]]; h2 = [zeros(nc); nh2[jc+1:end]]
    Λ2, L2, x2 = linecool_full(r2, u2, d2, t2, h2, k, jc, nc)
    Λ = zeros(nsh); Λ[jc+1:end] .= Λ2[nc+1:end]; Lr = zeros(nsh + 1); Lr[jc+1:end] .= L2[nc+1:end]
    Λ, Lr, x2
end
const LTAUCUT = parse(Float64, get(ENV, "LTAUCUT", "10")); const LINECORE = parse(Int, get(ENV, "LINECORE", "1"))
"Line cooling per zone (and Lr) for zone temperatures tt; step k (for NLTE=hyb); jc = number of core zones merged away."
function linecool_full(rr, uu, dd, tt, nh2, k, jc, ncz = jc > 0 ? 1 : 0)
    nz_ = length(dd); zk(j) = j <= ncz ? a + jc - 1 : a + (j - ncz) + jc - 1   # TITAN zone of reduced zone j (ncz merged-core zones first)
    if NLTE in ("ep", "ali", "hyb")
        yc = [ycomp(zk(j)) for j in 1:nz_]; nHn = dd ./ PC.μH
        nl = (LV = HLV, nH = [PC.xH(yc[j]) * nHn[j] for j in 1:nz_], nHe = PC.yHe .* nHn, nHp = [yc[j][1] * nHn[j] for j in 1:nz_])
        inc = lineinc()
        NLTE != "hyb" && return H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = NLTE == "ali" ? :ali : :ep, inc)
        tq = time(); Λe, Lr, x = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ep, inc)
        Λg = inc === nothing ? nothing : copy(H2Lines.LAST[].Λem)       # the EP emission alone (Λe is net of the absorbed background)
        TIMING && k >= 0 && pr("        EP call %.1f s (nep %d)\n", time() - tq, H2Lines.LAST[].nep)
        # (ALICAP=n, Nov 2026: the largest interval, default NALI; a separate variable because NALI is part of the restart check)
        # ALI refresh interval adaptive (≤ NALI): halved when the ALI/EP factor moved by > 5 % since the last refresh,
        # doubled when < 1 % — a factor frozen while the lines turn thick jumps at the refresh and kicks the core
        # (fixed mesh: the photosphere only moves the inner end of the reduced grid; the factors stay with their zones, aligned from the outside)
        if ENVEUL !== nothing && k >= 0 && !isempty(ALIFAC[]) && length(ALIFAC[]) != nz_
            o_ = ALIFAC[]; ALIFAC[] = [j + length(o_) - nz_ >= 1 ? o_[j+length(o_)-nz_] : 1.0 for j in 1:nz_]
            if length(ALIFIX[]) == length(o_)
                q_ = ALIFIX[]; ALIFIX[] = [j + length(q_) - nz_ >= 1 ? q_[j+length(q_)-nz_] : NaN for j in 1:nz_]
            end
        end
        if k >= 0 && (k - ALILAST[] >= ALIINT[] || length(ALIFAC[]) != nz_)
            tq = time(); Λa, Lra, xa = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ali, inc)
            TIMING && pr("        ALI call %.1f s (nep %d, ALI iterations %d, converged %s)\n", time() - tq, H2Lines.LAST[].nep, H2Lines.LAST[].iters, H2Lines.LAST[].converged)
            Λfs = Λa                                     # −∇·F of the formal solution (diagnostics)
            if LINECOLL
                P_ = H2Lines.LAST[].pops
                Λa = [H2Levels.collcool(HLV, tt[j], view(P_, :, j), nl.nH[j], nh2[j], nl.nHe[j], nl.nHp[j]) for j in 1:nz_]
                Lra = zeros(nz_ + 1); for j in 1:nz_; Lra[j+1] = Lra[j] + Λa[j] * 4π / 3 * (rr[j+1]^3 - rr[j]^3); end
                ALIFIX[] = [Λe[j] != 0 && Λa[j] / Λe[j] > 0 ? NaN : Λa[j] for j in 1:nz_]
                if Λg !== nothing
                    for j in 1:nz_; abs(Λe[j] - Λg[j]) > 0.5Λg[j] && (ALIFIX[][j] = Λa[j]); end
                    LINEINCDIAG && pr("        LINEINC k %d: T_bg %.1f T_star %.0f R_ph %.3e R_sun node %d of %d; held zones %d; sum V*(Λcoll %.3e, heating part %.3e)\n", k,
                                      inc.Tbg, inc.Tstar, inc.Rph / 6.957e10, H2Lines.LAST[].kq, nz_ + 1, count(!isnan, ALIFIX[]), Lra[end],
                                      sum(min(Λa[j], 0.0) * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:nz_))
                end
            end
            fnew = [Λe[j] != 0 && Λa[j] / Λe[j] > 0 && !(LINECOLL && Λg !== nothing && !isnan(ALIFIX[][j])) ? Λa[j] / Λe[j] : 1.0 for j in 1:nz_]
            if length(ALIFAC[]) == nz_
                # (S6-C2, ALIWMIN > 0, e.g. 1e-4) only over zones carrying ≥ ALIWMIN of the line luminosity: tiny zones (rezoned shock) with
                # negligible cooling otherwise keep the interval at 1 and the ALI (≈ 50× an EP step) runs every step
                wv = [abs(Λe[j]) * 4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:nz_]; wmin = ALIWMIN * sum(wv)
                δ = maximum((abs(log(fnew[j] / ALIFAC[][j])) for j in 1:nz_ if wv[j] >= wmin); init = 0.0)
                ALIINT[] = δ > ALISHRINK ? max(1, ALIINT[] ÷ 2) : δ < ALIGROW ? min(ALICAP, 2ALIINT[]) : ALIINT[]
                if TIMING
                    jd = argmax([wv[j] >= wmin ? abs(log(fnew[j] / ALIFAC[][j])) : 0.0 for j in 1:nz_])
                    pr("        ALI interval %d: δ %.3f at reduced zone %d of %d (r %.1f R_sun, T %.0f, n_H2 %.2e, share of |Λ|V %.2e, f %.3g -> %.3g), %d zones with δ > ALISHRINK\n", ALIINT[], δ, jd, nz_, rr[jd] / 6.957e10, tt[jd], nh2[jd], wv[jd] / sum(wv), ALIFAC[][jd], fnew[jd],
                       count(j -> wv[j] >= wmin && abs(log(fnew[j] / ALIFAC[][j])) > ALISHRINK, 1:nz_))
                end
            end
            if get(ENV, "LLDIAG", "0") == "1"           # how the ALI correction factor f = Λa/Λe closes the line luminosity
                vol = [4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:nz_]; bad = [Λe[j] == 0 || Λa[j] / Λe[j] <= 0 for j in 1:nz_]
                pr("        LLDIAG k %d: sum(Le) %.3e sum(La) %.3e sum(f Le) %.3e Lr_e %.3e Lr_a %.3e; sign-flip zones %d: sum(Le) %.3e sum(La) %.3e; +La %.3e -La %.3e; f median %.3f\n", k,
                   sum(Λe .* vol), sum(Λa .* vol), sum(Λe .* fnew .* vol), Lr[end], Lra[end], count(bad), sum(Λe[bad] .* vol[bad]; init = 0.0), sum(Λa[bad] .* vol[bad]; init = 0.0),
                   sum(max.(Λa, 0) .* vol), sum(min.(Λa, 0) .* vol), sort(fnew)[nz_÷2])
                ib = sortperm(abs.((Λe .* fnew .- Λa) .* vol); rev = true)[1:min(5, nz_)]
                for j in ib; pr("          zone %d (TITAN %d) r %.3e T %.1f n %.3e nH2 %.3e: Le %.3e La %.3e flag %s\n", j, zk(j), rr[j+1], tt[j], dd[j] / PC.mH, nh2[j], Λe[j] * vol[j], Λa[j] * vol[j], bad[j]); end
            end
            if get(ENV, "LINECONS", "0") == "1"         # per-zone consistency: −∇·F (Λa), net radiative rate (Λrate), net collisional (gas) exchange
                P_ = H2Lines.LAST[].pops; Λr_ = H2Lines.LAST[].Λrate; vol = [4π / 3 * (rr[j+1]^3 - rr[j]^3) for j in 1:nz_]
                Λc_ = [H2Levels.collcool(HLV, tt[j], view(P_, :, j), nl.nH[j], nh2[j], nl.nHe[j], nl.nHp[j]) for j in 1:nz_]
                pr("        LINECONS k %d: sum V*(Λa %.3e, Λrate %.3e, Λcoll %.3e) Lr %.3e; sum|Λa-Λcoll|V %.3e; maser lines %d, ALI it %d conv %s\n", k,
                   sum(Λfs .* vol), sum(Λr_ .* vol), sum(Λc_ .* vol), Lra[end], sum(abs.(Λfs .- Λc_) .* vol), H2Lines.LAST[].nmaser, H2Lines.LAST[].iters, H2Lines.LAST[].converged)
                for j in sortperm(abs.((Λfs .- Λc_) .* vol); rev = true)[1:min(6, nz_)]
                    pr("          zone %d (TITAN %d) r %.3e R_sun T %.1f nH2 %.2e: Λa %+.3e Λrate %+.3e Λcoll %+.3e\n", j, zk(j), rr[j] / 6.957e10, tt[j], nh2[j], Λfs[j], Λr_[j], Λc_[j])
                end
            end
            ALIFAC[] = fnew; ALILAST[] = k; NALICALL[] += 1
            return Λa, Lra, xa
        end
        length(ALIFAC[]) == nz_ || return Λe, Lr, x                        # (perturbed call right after the reduced grid changed: EP)
        Λs = Λe .* ALIFAC[]; Lrs = zeros(nz_ + 1)
        if LINECOLL && length(ALIFIX[]) == nz_
            for j in 1:nz_; isnan(ALIFIX[][j]) || (Λs[j] = ALIFIX[][j]); end
        end                          # Lr of the Λ actually applied (the EP Lr overstates L_lines by ~2)
        for j in 1:nz_; Lrs[j+1] = Lrs[j] + Λs[j] * 4π / 3 * (rr[j+1]^3 - rr[j]^3); end
        return Λs, Lrs, x
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
matstate() = mat isa PM.NetGas ? (; x = deepcopy(mat.x), eqzone = copy(mat.eqzone), Λext = copy(mat.Λext), Tref = copy(mat.Tref), slope = copy(mat.slope), Y = HOTON ? copy(SG[].Y) : nothing,
                                     star = HOTON ? (; fscale = copy(SG[].fscale), fslope = copy(SG[].fslope), fcur = copy(SG[].fcur), Dmix = copy(SG[].Dmix), dtcap = SG[].dtcap, capE = SG[].capE, ledger = SG[].ledger, Yprev = copy(SG[].Yprev)) : nothing,
                                     # lagged over one step (set in chemistry, used by the next hydro step): without them a restart is not exact
                                     lag = (; NLOPC = deepcopy(PM.NLOPC[]), GASQ = copy(TX.GASQ[]), RAYQ = copy(TX.RAYQ[]), SLDIAG = deepcopy(PM.SLDIAG[]), CHEMB = copy(TX.CHEMB[]))) :
             error("checkpoint: only the NetGas material (CHEM=net*) is supported")
function checkpoint(name, k)
    mkpath(CKDIR)
    TF.save_state(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    serialize(joinpath(CKDIR, "m2c_$(name).aux"), (; k, mat = matstate(), ALIFAC = copy(ALIFAC[]), ALIFIX = copy(ALIFIX[]), ALILAST = ALILAST[], ALIINT = ALIINT[],
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
    if HOTON && hasproperty(s, :star) && s.star !== nothing
        t_ = s.star; g_ = SG[]; g_.fscale .= t_.fscale; g_.fslope .= t_.fslope; g_.fcur .= t_.fcur; g_.Dmix = t_.Dmix; g_.dtcap = t_.dtcap; g_.capE = t_.capE; g_.ledger = t_.ledger; g_.Yprev .= t_.Yprev
    end
    if hasproperty(s, :lag)                        # (older checkpoints: these start empty, so the first step after a restart differs)
        l_ = s.lag; PM.NLOPC[] = deepcopy(l_.NLOPC); PM.SLDIAG[] = deepcopy(l_.SLDIAG)
        length(l_.GASQ) == length(TX.GASQ[]) && (TX.GASQ[] .= l_.GASQ); length(l_.RAYQ) == length(TX.RAYQ[]) && (TX.RAYQ[] .= l_.RAYQ)
        hasproperty(l_, :CHEMB) && length(l_.CHEMB) == length(TX.CHEMB[]) && (TX.CHEMB[] .= l_.CHEMB)
    end
    ALIFAC[] = copy(aux.ALIFAC); ALIFIX[] = hasproperty(aux, :ALIFIX) ? copy(aux.ALIFIX) : Float64[]; ALILAST[] = aux.ALILAST; ALIINT[] = aux.ALIINT; NALICALL[] = aux.NALICALL; NIMP[] = aux.NIMP
    LLINE[] = aux.LLINE; nextln[] = aux.nextln; SNAPNEXT[] = aux.SNAPNEXT
    SNAPS === nothing || aux.SNAPS === nothing || (empty!(SNAPS); append!(SNAPS, aux.SNAPS))
    empty!(Trbdf2.STATS); merge!(Trbdf2.STATS, aux.trb.STATS); Trbdf2.BE_LEFT[] = aux.trb.BE_LEFT
    TX.PEXT[] = aux.PEXT; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TF.setscalar!(:tmax, 1e3tff)                   # a run stopped by TSTOP/NSTOP has tmax = its last time
    KOFF[] = aux.k
    pr("restarted from checkpoint %s at step %d, t = %.6e yr, log n_c = %.3f, T_c = %.1f K\n", name, aux.k, TF.getscalar(:timen) / 3.156e7,
       log10(TF.var(:dn)[a] / PC.mH), TF.var(:tn)[a])
end
"""
THCHAIN=1 (default): the stellar EOS's p, e and MLT coefficients (∇_ad, c_P, δ) inside the Jacobian sweep are evaluated at the zone's values with a
2-partial derivative in (ρ, T), and passed into the outer Dual by the chain rule. Same values and same first derivatives as evaluating
StellarEOS.eos on the outer Dual (whose own nested derivative then carries all of the sweep's partials through every operation), at a
fraction of the cost (it was ~40 % of the run time). THCHAIN=0: the direct evaluation.
"""
const THCHAIN = get(ENV, "THCHAIN", "1") == "1"
struct ThTag end
function thchain(f, ρ, T, k)
    (ρ isa ForwardDiff.Dual || T isa ForwardDiff.Dual) || return f(ρ, T, k)
    Dout = typeof(ρ + T)
    (ForwardDiff.valtype(Dout) === Float64 && (!(ρ isa ForwardDiff.Dual) || typeof(ρ) === Dout) && (!(T isa ForwardDiff.Dual) || typeof(T) === Dout)) || return f(ρ, T, k)   # nested / mixed tags: direct
    ρv, Tv = ForwardDiff.value(ρ), ForwardDiff.value(T)
    D2 = ForwardDiff.Dual{typeof(ForwardDiff.Tag(ThTag(), Float64)),Float64,2}
    o = f(D2(ρv, ForwardDiff.Partials((1.0, 0.0))), D2(Tv, ForwardDiff.Partials((0.0, 1.0))), k)
    pρ = ρ isa ForwardDiff.Dual ? ForwardDiff.partials(ρ) : zero(ForwardDiff.partials(T))
    pT = T isa ForwardDiff.Dual ? ForwardDiff.partials(T) : zero(ForwardDiff.partials(ρ))
    map(x -> Dout(ForwardDiff.value(x), ForwardDiff.partials(x, 1) * pρ + ForwardDiff.partials(x, 2) * pT), o)
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
    let g0 = (ρ, T, k) -> (st = SE.gas_state(ρ, T, SM.zonecomp(SG[], k)); (st.p, st.e)),
        g = THCHAIN ? ((ρ, T, k) -> thchain(g0, ρ, T, k)) : g0,
        th0 = (ρ, T, k) -> (st = SE.eos(ρ, T; comp = SM.zonecomp(SG[], k), radiation = true); (st.∇ad, st.cP, st.δ)),
        th = THCHAIN ? ((ρ, T, k) -> thchain(th0, ρ, T, k)) : th0
        ρr, Tr = 1e-4, 1e6                            # zero point of the two EOS in the fully ionised dilute limit
        Δ0 = g(ρr, Tr, 1)[2] - PopIIIEqEOS.eos(ρr, Tr)[1]
        ee = (ρ, T, k) -> SG[].fcur[k] * SM.eps_nuc(ρ, T, view(SG[].Y, :, k); screen = SG[].screen)
        PM.HOT[] = (; T1 = parse(Float64, get(ENV, "HOTT1", "3e5")), T2 = parse(Float64, get(ENV, "HOTT2", "9e5")), Δ0, gas = g, thermo = th, eps = ee)
    end
    pr("HOT: stellar branch above T1 = %.3g K (full at %.3g K), Delta0 = %.4e erg/g, D/H = %.2e\n", PM.HOT[].T1, PM.HOT[].T2, PM.HOT[].Δ0, DH)
    "e_gas of the blended material for an explicit stellar composition (the zone's, or a trial one in mix!)."
    function e_blend(ρ, T, comp)
        h = PM.HOT[]; T <= h.T1 && return PopIIIEqEOS.eos(ρ, T)[1]
        es = SE.gas_state(ρ, T, comp).e - h.Δ0; T >= h.T2 && return es
        u = PopIIIEqEOS.eos(ρ, T)[1]; u + PM.hotw(T, h.T1, h.T2) * (es - u)
    end
    "T with e_blend(ρ, T, comp) = e (safeguarded Newton in ln T), and the relative error of e."
    function T_blend(ρ, e, T0, comp)
        lT = log(T0); lo = lT - 1.5; hi = lT + 1.5; f = Inf
        for _ in 1:60
            e0 = e_blend(ρ, exp(lT), comp); f = e0 - e
            abs(f) <= 1e-12 * abs(e) && break
            f > 0 ? (hi = lT) : (lo = lT)
            d = (e_blend(ρ, exp(lT + 1e-6), comp) - e0) / 1e-6
            lN = lT - f / d; lT = (d > 0 && lo < lN < hi) ? lN : 0.5 * (lo + hi)
        end
        exp(lT), abs(f) / abs(e)
    end
    SM.TSOLVE[] = T_blend; SM.ZONEOK[] = k -> TF.var(:tn)[k] > PM.HOT[].T1
    # S5's coupling: every step through solve_step! (level factor f_k on ε iterated against the network), the stage ledger, Δt caps on composition
    # and energy; after the step burn! (network along the T, ρ history of the step, the difference to what the residual injected into e) and mix!
    SG[].nuc = true; SG[].dXmax = parse(Float64, get(ENV, "DXMAX", "4e-3")); SG[].efrac = parse(Float64, get(ENV, "EFRAC", "5"))
    SM.REGION_CAP[] = true; SM.install_stage_hook!(SG[])
    push!(TC.DT_HOOKS, () -> min(SG[].dtcap, SM.predict_cap(SG[])))
    const SPLIT = (tol = parse(Float64, get(ENV, "SPLITTOL", "3e-3")), maxit = parse(Int, get(ENV, "SPLITIT", "4")), ztol = parse(Float64, get(ENV, "ZSPLITTOL", "3e-3")))
    const LASTSPLIT = Ref((NaN, 0, NaN)); const STEP0 = TD.STEPFN[]
    TD.STEPFN[] = () -> (LASTSPLIT[] = SM.solve_step!(SG[], STEP0; SPLIT...); nothing)
    if get(ENV, "S5NUM", "0") == "1"                  # S5's numerics for quasi-static stars: T, E_r predictor, limited grid extrapolation, gentle Δt cut after a failure
        SM.install_state_predictor!(); TC.EXTRAP_FRAC[] = 0.3; TC.PRED_LAGRHO[] = true; TD.FAIL_TFAC[] = 0.5; Trbdf2.TFAC[] = 0.5
        TC.DT_RULE[] = function (dt0)                 # S5's rule: Δt on max |Δ ln T, ρ, r| per step, growth ≤ 1.6
            tn, to, dn, do_, rn, ro = (TF.var(s_) for s_ in (:tn, :to, :dn, :do, :rn, :ro)); d = 0.0
            for k in a:b; d = max(d, abs(log(tn[k] / to[k])), abs(log(dn[k] / do_[k]))); end
            for k in a+1:b+1; d = max(d, abs(log(rn[k] / ro[k]))); end
            f = clamp(0.9 * parse(Float64, get(ENV, "S5TOL", "0.03")) / max(d, 1e-12), 0.5, 1.6)
            TF.getscalar(:iter) > 20 && (f = min(f, 0.8))
            f * dt0
        end
    end
    const HOTR = Ref{Any}((; Ebur = 0.0, corr = 0.0)); const HOTM = Ref{Any}((; nconv = 0)); const EBURN = Ref(0.0)
    "S6-C1: after a converged step and chem!: burn! and mix! of the hot zones (S5)."
    function hotburn!(dt)
        sg = SG[]
        if !any(TF.var(:tn)[k] > PM.HOT[].T1 for k in a:b); sg.acc.nstage = 0; return 0.0; end
        sd = TX.step_diag(dt)
        HOTR[] = SM.burn!(sg; sd); EBURN[] += HOTR[].Ebur
        HOTM[] = SM.mix!(sg; sd, capmix = false, bridge = 4, tmax = 0.02)
        dt > 0 ? HOTR[].Ebur / dt : 0.0      # a step of zero duration (Δt below the Float64 spacing of t) burns nothing
    end
end
"""
S6-B: settings of a radiating run (after a RAD rebuild and after restarting a RAD checkpoint): StellarOpacity continuum for
NetGas, CMB intensity at the outer boundary, Eddington-factor mode, thin-envelope E_r out of the Δt rule, Δt cap.
"""
const DTRHOOK = Ref(false); const XMEXT = Ref{Any}(nothing)
function radsettings!()
    X = 1 / (1 + 4PC.yHe)
    # OPACSENS="fcond,fext" (sensitivity test, default 1,1): κ_cond × fcond, and κ_rad × fext beyond OPAL's log R = 1 (over one dex in log R)
    os_ = parse.(Float64, split(get(ENV, "OPACSENS", "1,1"), ","))
    PM.OPAC[] = all(==(1.0), os_) ? ((ρ, T) -> (StellarOpacity.kappa_R(ρ, T, X), StellarOpacity.kappa_P(ρ, T, X))) : ((ρ, T) -> begin
        κ = StellarOpacity.kappa_R(ρ, T, X); κr = StellarOpacity.kappa_R_rad(ρ, T, X); ic = max(1 / κ - 1 / κr, zero(κ))     # 1/κ_cond,eff
        lR = log10(ρ) - 3 * (log10(T) - 6); w = clamp(lR - 1, 0.0, 1.0); w = w * w * (3 - 2w)
        (1 / (1 / (κr * os_[2]^w) + ic / os_[1]), StellarOpacity.kappa_P(ρ, T, X))
    end)
    Tb = parse(Float64, get(ENV, "TCMB", "0"))
    TF.setscalar!(:ximr, Z0.σSB * Tb^4 / π)
    TF.setscalar!(:ltran, parse(Int, get(ENV, "LTRAN", "1")))      # 1: formal solution for f, g; 2: Eddington approximation f = 1/3, g = 1/2
    TC.ER_TAUCONV[] = parse(Float64, get(ENV, "ERTAUCONV", "1.0"))   # causality floor of E_r outside the photosphere
    TC.ER_TAUMIN[] = parse(Float64, get(ENV, "ERTAUMIN", "1.0"))   # E_r of the optically thin envelope does not limit Δt
    # Δt cap: steps much longer than the light-crossing time of the nearly free-streaming layers outside the photosphere
    # (r ~ 1e14 cm: r/c ~ 3e3 s) fail in Newton (Eddington factors lagged over the step); default 1500 s
    if HOTON                                          # convection (hot branch only: PM.material_thermo), as the S5 stars
        TX.CONVECTION[] = true; TX.ALPHA_MLT[] = parse(Float64, get(ENV, "ALPHAMLT", "1.5")); TX.MLT_HSDISC[] = true; TX.MLT_RCAP[] = true
        PM.COLDMLT[] = parse(Float64, get(ENV, "COLDMLT", "0")); PM.COLDMLT[] > 0 && (TX.MLT_NOINFALL[] = true)
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
              selfgrav = true, material = mat, av, xmext = XMEXT[])
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    for j in 1:nsh
        k = a + j - 1
        mat.x[k] = copy(matz.x[j]); mat.eqzone[k] = matz.eqzone[j]; mat.Λext[k] = matz.Λext[j]; mat.Tref[k] = matz.Tref[j]; mat.slope[k] = matz.slope[j]
    end
    if HOTON
        g_ = SG[]
        for j in 1:nsh; g_.Y[:, a+j-1] .= matz.Y[j]; end
        SM.ghosts!(g_); g_.Yprev .= g_.Y; g_.acc.nstage = 0
        fill!(g_.fscale, 1.0); fill!(g_.fslope, 0.0); fill!(g_.fcur, 1.0); g_.Dmix = zeros(length(TF.var(:dn)))
        if hasproperty(matz, :fslope)                 # (rezoning: level-factor slopes and the mixed regions follow the zones)
            for j in 1:nsh; g_.fscale[a+j-1] = matz.fscale[j]; g_.fslope[a+j-1] = matz.fslope[j]; g_.Dmix[a+j-1] = matz.Dmix[j]; end
        end
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
PM.PSWTMIN[] = parse(Float64, get(ENV, "PSWTMIN", "0"))
PM.PSWTEX[] = parse(Float64, get(ENV, "PSWTEX", "0")); PM.PSWTAU[] = parse(Float64, get(ENV, "PSWTAU", "0"))
const REZW = haskey(ENV, "REZWIN") ? parse.(Float64, split(ENV["REZWIN"], ",")) : nothing; const NREZ = Ref(0); const REZNULL = get(ENV, "REZNULL", "0") == "1"
"""
S6-C2: rezone (Rezone.plan/apply) and rebuild. The material state of a split zone is copied; of a merged zone, mass-weighted
(x), equilibrium flag if either was, Λext volume-weighted (it is per volume), Tref/slope of the heavier zone. Returns the
numbers of splits and merges.
"""
const REZDR = parse(Float64, get(ENV, "REZDR", "0"))
Rezone.SLOPE_LNR[] = get(ENV, "REZSLOPE", "") == "lnr"      # slopes of a split per unit of ln r (Rezone.SLOPE_LNR)
function rezone!(k; force = nothing)
    p = profile_now()
    # moving window (REZWIN=flo,fhi,frac): zones with m_enc in [flo, fhi] × m★ (m★: inside the first interface falling faster than
    # 2 km/s) and heavier than frac × m_enc are split while they are still smooth, ahead of the star
    fl = nothing
    if REZW !== nothing
        dm_ = Rezone.dmass(p); me_ = cumsum(dm_); jc_ = findfirst(<(-2e5), p.u); mst = jc_ === nothing ? 0.0 : me_[max(jc_ - 1, 1)]
        fl = [REZW[1] * mst <= me_[j] <= REZW[2] * mst && dm_[j] > REZW[3] * me_[j] for j in eachindex(dm_)]
        # REZDR=x: inside the same outer mass limit, zones from the surface outwards are also split while Δr/r̄ > x, so that the
        # infall is resolved in radius (a falling Lagrangian shell stretches, Δr/r ∝ dm r^(-3/2) in steady free fall)
        if REZDR > 0 && jc_ !== nothing
            for j in max(jc_ - 1, 1):length(dm_)
                me_[j] <= REZW[2] * mst && 2 * (p.rI[j+1] - p.rI[j]) / (p.rI[j+1] + p.rI[j]) > REZDR && (fl[j] = true)
            end
        end
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
              Y = HOTON ? [sum(w .* old.Y[:, zx(i)] for (i, w) in s) for s in src] : nothing,
              fscale = HOTON ? [old.star.fscale[zx(s[1][1])] for s in src] : nothing, fslope = HOTON ? [old.star.fslope[zx(s[1][1])] for s in src] : nothing,
              Dmix = HOTON ? [old.star.Dmix[zx(s[1][1])] for s in src] : nothing)
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
"""
S6 validation: continue an S5 star (s5_run.jl checkpoint `STARFROM`, directory S5DIR) on the NetGas hot path: same profile, abundances, solver
controls and time; then T_c, L, X_c(t) can be compared with S5's own continuation.
"""
function starfrom!(name)
    dir = get(ENV, "S5DIR", "/tmp/s5"); TF.load_state!(joinpath(dir, "s5_$(name).tstate.gz"))
    aux = deserialize(joinpath(dir, "s5_$(name).aux"))
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre))
    rI = TF.var(:rn)[a:b+1]; ρz = TF.var(:dn)[a:b]; Tz = TF.var(:tn)[a:b]; uI = TF.var(:un)[a:b+1]; Er = TF.var(:ern)[a:b]; Fr = TF.var(:frn)[a:b+1]
    n = b - a + 1; XMEXT[] = TF.getscalar(:xmext); TX.PEXT[] = 0.0; aold = a
    matz = (; x = [PM.x_from_table(mat, ρz[j], min(Tz[j], 1e6)) for j in 1:n], eqzone = trues(n), Λext = zeros(n), Tref = ones(n), slope = zeros(n),
              Y = [aux.Y[:, aold+j-1] for j in 1:n], fscale = ones(n), fslope = [aux.fslope[aold+j-1] for j in 1:n], Dmix = [aux.Dmix[aold+j-1] for j in 1:n])
    rebuild!(rI, ρz, Tz, uI, Er, Fr, matz)
    SG[].dtcap = aux.dtcap; SG[].capE = aux.capE; TF.setscalar!(:tmax, 1e30)
    pr("STARFROM %s: %d zones, t = %.6e yr, T_c = %.5e K, rho_c = %.3f, X_c = %.5f\n", name, nsh, TF.getscalar(:timen) / 3.156e7, TF.var(:tn)[a], TF.var(:dn)[a], SG[].Y[SM.IH, a])
end
# ---------------------------------------------------------------------------------------- fixed logarithmic mesh for the infall
"""
ENVEUL=rin,dlnr_in,dlnr_out,rb (radii in R☉): remap the restarted state onto a fixed mesh (lgrid = 2) from r_in, inside the star and
so inside the shock, to the cloud edge: Δln r = dlnr_in below rb, growing 2 % per zone to dlnr_out above. Two optional further entries,
ENVEUL=rin,dlnr_in,dlnr_out,rb,gout,gin, set the factor per zone by which Δln r grows above rb (default 1.02) and shrinks from the core
spacing to dlnr_in above the core boundary of ENVCORE (default 0.95), for a mesh with fewer zones in the envelope and the core. Mass inside r_in becomes the
central point mass; the inner wall passes mass onto it with TITAN's INNER_SINK, holding the mass below the shock at its initial value
(relaxation time ENVSINKTAU s). Gas and radiation are advected by TITAN's implicit rows (GASROW, as the quasi-star runs). Species are not
TITAN unknowns: `advect_species!` moves them with the converged mass fluxes.
"""
const ENVEUL = haskey(ENV, "ENVEUL") ? parse.(Float64, split(ENV["ENVEUL"], ",")) : nothing
const Rsun = 6.957e10
const ENVCENTRE = ENVEUL !== nothing && ENVEUL[1] == 0     # r_in = 0: the mesh runs to a regular centre (ENVCORE=rmin,dlnr,rcb: first interface [R☉], Δln r below rcb); no point mass, wall or sink
const ENVCORE = parse.(Float64, split(get(ENV, "ENVCORE", "0.1,0.08,10"), ","))
const ENVCONS = get(ENV, "ENVCONS", "0") == "1"       # the ENVEUL remap takes the zone densities from the source's enclosed mass (mass-conserving); default: log-interpolated at the zone centres
const ENVSINKTAU = parse(Float64, get(ENV, "ENVSINKTAU", "1e6"))
include(joinpath(@__DIR__, "GridMap.jl"))
"""
GRIDMAP=N[,hr,h1,h2,eps,fc,dmin,g,tau,fmove,nrelax,ctol,ftar] (with ENVEUL; default off): one mesh of N zones for the whole range of scales that
follows the star. Every zone covers the same arc length in ln(r + r_c)/hr, ln(m + fc M★)/h1 and asinh((m − M★)/(eps M★))/h2 (GridMap.jl:
radius where the enclosed mass does not change, mass where it does, counted from the centre and from the surface of the star; M★ from
the velocity profile), no zone narrower than dmin in ln r, neighbours within the factor g. At the hand-over from the Lagrangian mesh
this is the starting mesh (ENVCORE's first interface and the outer edge are kept); afterwards the target is computed again after every
step and each interface moves towards it at the velocity (r_target − r)/max(tau, nrelax Δt), at most fmove of its smaller neighbouring
zone per step (TitanResidFast.GRIDVEL: prescribed velocities, TITAN's moving-grid residual; abundances move with the mass that crosses
each interface, advect_species!). Δt is limited so that the target moves at most ftar zone widths per step. Defaults hr 0.3, h1 0.3,
h2 0.4, eps 1e-5, fc 1e-3, dmin 0.02, g 1.2, tau 0, fmove 0.4, nrelax 3, ctol 0.5 (maxdel: relative zone-width change per Newton
iteration), ftar 0.2. The first interface above the centre and the outer edge
never move. A checkpoint of such a run restarts with the same setting (any lgrid = 2 mesh relaxes to the target).
"""
const GRIDMAPP = haskey(ENV, "GRIDMAP") ? parse.(Float64, split(ENV["GRIDMAP"], ",")) : nothing
gmp(i, d) = GRIDMAPP !== nothing && length(GRIDMAPP) >= i ? GRIDMAPP[i] : d
gridmap_kw() = (; hr = gmp(2, 0.3), h1 = gmp(3, 0.3), h2 = gmp(4, 0.4), eps = gmp(5, 1e-5), fc = gmp(6, 1e-3), dmin = gmp(7, 0.02), g = gmp(8, 1.2))
const GMPREV = Ref(Float64[])      # the target of the previous step (for its speed)
const GMLAST = Ref((0.0, 0, 0.0, NaN, 0, 0.0, 0.0, Inf))    # (move planned for the next step / zone width, interface, M★ per steradian, scale λ, zones at the floor, largest distance from the target / zone width)
xmi(k) = ENVEUL === nothing ? TX.XMINT[][k] : TF.getscalar(:xmtot) - TF.var(:xmen)[k]
const XMEPREV = Ref(Float64[]); const ENVMS0 = Ref(0.0); const ENVMDOT = Ref(0.0)
envmenc() = TF.getscalar(:xmtot) .- TF.var(:xmen)[a:b+1]              # interior mass per steradian at the interfaces
function envsurface()
    j = findfirst(<(-2e5), TF.var(:un)[a:b+1]); j === nothing ? nsh + 1 : j
end
function envsink!()
    ENVCENTRE && return
    xme = TF.var(:xmen); j = envsurface(); Ms = xme[a] - xme[a+j-1]
    Md = max(0.0, Ms - ENVMS0[]) / ENVSINKTAU; ENVMDOT[] = Md                  # per steradian
    φ0 = -Md / TF.var(:rn)[a]^2
    TF.setscalar!(:phil0, φ0); TF.setscalar!(:phil2, φ0 * (TF.var(:egn)[a] + TF.var(:pgn)[a] / TF.var(:dn)[a]))
end
# XADV2=1 (default with SL3=1): the abundances cross each interface with the van Leer interface value TITAN uses for e (`advect!`, centred,
# cadv; start-of-step values), explicitly, instead of the donor value. With x_e ~ 0.1–0.5 the first-order scheme carried a different amount of
# ionization energy across each interface than e did: up to ~1 eV per H a step, enough to leave e below u(10 K, x). Falls back to the
# Applied as an explicit correction to the implicit upwind step, at interfaces where less than the donor zone's mass crosses, and per species
# only if it keeps that species ≥ 0 everywhere; `NXADV1` counts the steps where some species lost the correction.
const XADV2 = get(ENV, "XADV2", get(ENV, "SL3", "0")) == "1"; const NXADV1 = Ref(0)
"Second-order correction to the implicit upwind result `xnew` (vectors per zone, from the start-of-step `X`): see the XADV2 note."
function xadv2_correct!(xnew, X, F, M, n)
    cadv = TF.getscalar(:cadv); nf = length(X[1])
    S = [zeros(nf) for _ in 1:n]
    for j in 2:n-1, i in 1:nf
        da = X[j][i] - X[j-1][i]; db = X[j+1][i] - X[j][i]
        S[j][i] = da * db > 0 ? cadv * da * db / (da + db) : 0.0
    end
    # antidiffusive flux through interior interface j (between zones j-1, j): F (x_face − x_donor), only where |F| < the donor's mass
    G = [zeros(nf) for _ in 1:n+1]
    for j in 2:n
        if F[j] >= 0
            F[j] < M[j-1] && (G[j] .= F[j] .* 0.5 .* S[j-1])
        else
            -F[j] < M[j] && (G[j] .= F[j] .* (-0.5) .* S[j])
        end
    end
    ok = trues(nf); xc = [xnew[j] .+ (G[j] .- G[j+1]) ./ M[j] for j in 1:n]
    for j in 1:n, i in 1:nf; xc[j][i] < 0 && (ok[i] = false); end
    for j in 1:n
        for i in 1:nf; ok[i] || (xc[j][i] = xnew[j][i]); end
        d = ChemNet.full(xc[j], mat.net)
        (d[:H] >= 0 && d[:e] >= 0) || (xc[j] = copy(xnew[j]))
    end
    all(ok) || (NXADV1[] += 1)
    for j in 1:n; xnew[j] = xc[j]; end
    nothing
end
"""
Implicit upwind transport of mat.x with the mass that crossed each interface in the converged step (change of the exterior mass m_e), so
that it is the discrete continuity equation applied to the mass fractions: a uniform composition stays uniform and elements are conserved.
"""
function advect_species!()
    mat isa PM.NetGas && PM.HDCHEM[] && PM.hd_pad!(mat)                       # (first step after a restart from a run without HD)
    xme = TF.var(:xmen); old = XMEPREV[]; n = nsh
    F = [xme[a+j-1] - old[a+j-1] for j in 1:n+1]                              # > 0: mass moved outwards through interface j
    M = [xme[a+j-1] - xme[a+j] for j in 1:n]
    lo = zeros(n); up = zeros(n); di = zeros(n)
    for j in 1:n
        di[j] = M[j] + max(F[j+1], 0.0) + max(-F[j], 0.0)
        j > 1 ? (lo[j] = -max(F[j], 0.0)) : (di[j] -= max(F[j], 0.0))         # (no inflow through a boundary: own composition)
        j < n ? (up[j] = -max(-F[j+1], 0.0)) : (di[j] -= max(-F[j+1], 0.0))
    end
    M0 = [old[a+j-1] - old[a+j] for j in 1:n]
    cp_ = zeros(n); den = zeros(n)
    den[1] = di[1]; cp_[1] = up[1] / den[1]
    for j in 2:n; den[j] = di[j] - lo[j] * cp_[j-1]; cp_[j] = up[j] / den[j]; end
    xnew = [M0[j] .* mat.x[a+j-1] for j in 1:n]
    xnew[1] = xnew[1] ./ den[1]
    for j in 2:n; xnew[j] = (xnew[j] .- lo[j] .* xnew[j-1]) ./ den[j]; end
    for j in n-1:-1:1; xnew[j] = xnew[j] .- cp_[j] .* xnew[j+1]; end
    XADV2 && mat isa PM.NetGas && xadv2_correct!(xnew, [copy(mat.x[a+j-1]) for j in 1:n], F, M, n)
    for j in 1:n; mat.x[a+j-1] = xnew[j]; end
    if HOTON                                                                 # the nuclear abundances move with the same fluxes
        Y = SG[].Y; ynew = M0[1] .* Y[:, a]; ys = Vector{Vector{Float64}}(undef, n); ys[1] = ynew ./ den[1]
        for j in 2:n; ys[j] = (M0[j] .* Y[:, a+j-1] .- lo[j] .* ys[j-1]) ./ den[j]; end
        for j in n-1:-1:1; ys[j] = ys[j] .- cp_[j] .* ys[j+1]; end
        for j in 1:n; Y[:, a+j-1] .= ys[j]; end
        SM.ghosts!(SG[])
    end
    copyto!(old, xme)
    nothing
end
function envflags!()
    TX.REGULAR_CENTRE[] = ENVCENTRE; empty!(TX.XMINT[]); TX.GASROW[] = true; TX.INNER_SINK[] = !ENVCENTRE
    PM.CHEMFALLBACK[] = true
    XMEPREV[] = copy(TF.var(:xmen)::Vector{Float64})
end
function envbuild!()
    rin, dl1, dl2, rb = ENVEUL[1] * Rsun, ENVEUL[2], ENVEUL[3], ENVEUL[4] * Rsun
    r0 = TF.var(:rn)[a:b+1]; u0 = TF.var(:un)[a:b+1]; d0 = TF.var(:dn)[a:b]; T0 = TF.var(:tn)[a:b]; E0 = TF.var(:ern)[a:b]; F0 = TF.var(:frn)[a:b+1]
    m0 = 4π .* TX.XMINT[][a:b+1]; x0 = [copy(mat.x[a+j-1]) for j in 1:nsh]; eq0 = [mat.eqzone[a+j-1] for j in 1:nsh]
    t0 = TF.getscalar(:timen); js = TF.getscalar(:jstep)
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    rcb = ENVCENTRE ? ENVCORE[3] * Rsun : 0.0
    lr = [log(ENVCENTRE ? ENVCORE[1] * Rsun : rin)]; dl = ENVCENTRE ? ENVCORE[2] : dl1
    gout_ = length(ENVEUL) >= 5 ? ENVEUL[5] : 1.02; gin_ = length(ENVEUL) >= 6 ? ENVEUL[6] : 0.95
    while lr[end] < log(r0[end])
        r_ = exp(lr[end])
        r_ > rb ? (dl = min(dl2, gout_ * dl)) : r_ > rcb && (dl = max(dl1, gin_ * dl))
        push!(lr, lr[end] + dl)
    end
    (length(lr) + (ENVCENTRE ? 1 : 0)) % 2 == 0 && pop!(lr)                    # (an even zone count: an odd number of interfaces)
    lr = lr[1] .+ (lr .- lr[1]) .* (log(r0[end]) - lr[1]) / (lr[end] - lr[1])
    if GRIDMAPP !== nothing      # the planned mesh (GridMap) from the source's enclosed mass in place of the ENVCORE/ENVEUL spacing; same first and last interface
        i1 = findfirst(j -> r0[j] > 0 && m0[j] > 0, eachindex(r0)); Nz = Int(GRIDMAPP[1]); (Nz + (ENVCENTRE ? 1 : 0)) % 2 == 0 || (Nz -= 1)
        Mst = GridMap.anchor(m0[i1:end], u0[i1:end])
        P_ = GridMap.target(r0[i1:end], m0[i1:end], Mst; N = Nz, rlo = exp(lr[1]), gridmap_kw()...)
        lr = log.(P_.r); GMLAST[] = (0.0, 0, Mst / 4π, P_.lam, P_.nfloor, 0.0, 0.0, Inf)
        pr("GRIDMAP: starting mesh of %d zones from the enclosed mass of the source, M* %.5f Msun, scale %.3f, %d zones at the floor\n", Nz, Mst / Msun, P_.lam, P_.nfloor)
    end
    rI = exp.(lr); rI[end] = r0[end]; ENVCENTRE && pushfirst!(rI, 0.0); n = length(rI) - 1
    rc = [0.5 * (r0[j] + r0[j+1]) for j in 1:nsh]; lc = log.(rc); lf = log.(max.(r0, 1.0))
    function interp(xs, ys, x)                                                 # linear in the first argument, flat beyond the ends
        x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
        i = searchsortedlast(xs, x); w = (x - xs[i]) / (xs[i+1] - xs[i])
        (1 - w) * ys[i] + w * ys[i+1]
    end
    rz = [0.5 * (rI[j] + rI[j+1]) for j in 1:n]; lz = log.(rz)
    ρz = [exp(interp(lc, log.(d0), l)) for l in lz]; Tz = [exp(interp(lc, log.(T0), l)) for l in lz]
    Er = [exp(interp(lc, log.(E0), l)) for l in lz]
    lfi = log.(rI); mI = [interp(lf, m0, l) for l in lfi]
    if ENVCONS      # zone means of a monotone cubic of ln m against ln r through the source's interfaces: the mass inside every source interface is kept
        i1 = findfirst(j -> r0[j] > 0 && m0[j] > 0, eachindex(r0)); xs = log.(r0[i1:end]); ys = log.(m0[i1:end]); N = length(xs)
        h = diff(xs); δ = diff(ys) ./ h; dd = zeros(N)
        for i in 2:N-1; dd[i] = δ[i-1] > 0 && δ[i] > 0 ? 2 / (1 / δ[i-1] + 1 / δ[i]) : 0.0; end        # harmonic mean of the two slopes: monotone
        dd[1] = min(3.0, 3δ[1]); dd[N] = δ[N-1]
        function lnm(x)
            x <= xs[1] && return ys[1] + 3 * (x - xs[1])                       # uniform density inside the first interface that encloses mass
            x >= xs[N] && return ys[N]
            i = searchsortedlast(xs, x); t = (x - xs[i]) / h[i]
            (2t^3 - 3t^2 + 1) * ys[i] + (t^3 - 2t^2 + t) * h[i] * dd[i] + (3t^2 - 2t^3) * ys[i+1] + (t^3 - t^2) * h[i] * dd[i+1]
        end
        mI = [rI[j] > 0 ? exp(lnm(log(rI[j]))) : 0.0 for j in eachindex(rI)]; mI[end] = m0[end]
        ρint = ρz; ρz = [(mI[j+1] - mI[j]) / (4π / 3 * (rI[j+1]^3 - rI[j]^3)) for j in 1:n]
        all(>(0), ρz) || error("ENVCONS: a zone without mass")
        jm = argmax(abs.(log.(ρz ./ ρint)))
        pr("ENVCONS: densities from the enclosed mass; largest change against the interpolated value: factor %.3g in zone %d (r %.4g .. %.4g R_sun)\n",
           ρz[jm] / ρint[jm], jm, rI[jm] / Rsun, rI[jm+1] / Rsun)
    end
    uI = [interp(lf, u0, l) for l in lfi]
    ncap = 0
    for j in eachindex(uI)                                                     # the state keeps shells past the escape speed: cap them
        ve = sqrt(2Z0.G * mI[j] / rI[j]); uI[j] < -ve && (uI[j] = -ve; ncap += 1)
    end
    uI[1] = 0.0; uI[end] = 0.0
    Ln = 4π .* r0 .^ 2 .* F0; Fr = [rI[j] > 0 ? interp(lf, Ln, lfi[j]) / (4π * rI[j]^2) : 0.0 for j in eachindex(rI)]
    Lin = abs(4π * rI[1]^2 * Fr[1])
    idx(l) = clamp(searchsortedlast(lc, l), 1, nsh - 1)
    xz = [(i = idx(l); w = clamp((l - lc[i]) / (lc[i+1] - lc[i]), 0.0, 1.0); (1 - w) .* x0[i] .+ w .* x0[i+1]) for l in lz]
    eqz = [eq0[clamp(searchsortedlast(lc, l), 1, nsh)] for l in lz]
    Mc = ENVCENTRE ? 0.0 : ENVCONS ? mI[1] : interp(lf, m0, lfi[1]); Mnew = sum(4π / 3 * ρz[j] * (rI[j+1]^3 - rI[j]^3) for j in 1:n)
    TX.GASROW[] = true
    # per-zone material on the new indices before build! fills the EOS (as in the ENVMERGE rebuild: else a new zone far out inherits the
    # table-EOS flag of the old zone with its index: "rho below the table" when the source has more zones in the star than the target
    # has in all, e.g. a refined stage A3 onto 64 zones); set again below on the rebuilt grid
    for k in eachindex(mat.x); j = clamp(k - a + 1, 1, n); mat.x[k] = copy(xz[j]); mat.eqzone[k] = eqz[j]; end
    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lum = ENVCENTRE ? nothing : Lin, lgrid = 2, inner = :wall, outer = :wall, Mc,
              selfgrav = true, material = mat, av, xmext = (Mc + Mnew) / 4π)
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    for j in 1:nsh
        k = a + j - 1; mat.x[k] = xz[j]; mat.eqzone[k] = eqz[j]; mat.Λext[k] = 0.0; mat.Tref[k] = Tz[j]; mat.slope[k] = 0.0
    end
    for k in 1:a-1; mat.x[k] = copy(mat.x[a]); mat.eqzone[k] = mat.eqzone[a]; end
    for k in b+1:length(mat.x); mat.x[k] = copy(mat.x[b]); mat.eqzone[k] = mat.eqzone[b]; mat.Λext[k] = 0.0; end
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    PS.install_julia!(BlockSolver)
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, t0); TF.setscalar!(:timeo, t0); TF.setscalar!(:dtime, parse(Float64, get(ENV, "ENVDT0", "100"))); TF.setscalar!(:jstep, js); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    envflags!()
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    empty!(ALIFAC[]); SNAPS === nothing || empty!(SNAPS)
    xme = TF.var(:xmen); ENVMS0[] = xme[a] - xme[a+envsurface()-1]
    pr("ENVEUL: %d zones, r %.3g R_sun .. %.3g pc, Mc %.4f Msun, grid mass %.1f (was %.1f) Msun, L_in %.3e erg/s, %d velocities capped at v_esc, surface zone %d, M_shell0 %.3e Msun\n",
       n, rI[1] / Rsun, rI[end] / 3.0857e18, Mc / Msun, Mnew / Msun, (m0[end] - Mc) / Msun, Lin, ncap, envsurface(), 4π * ENVMS0[] / Msun)
end
if haskey(ENV, "STARFROM")
    starfrom!(ENV["STARFROM"])
elseif haskey(ENV, "RADFROM")
    RAD || error("RADFROM needs RAD=1")
    radbuild!(ENV["RADFROM"])
    haskey(ENV, "REFINE") && prerefine!(KOFF[], parse.(Float64, split(ENV["REFINE"], ","))...)
elseif haskey(ENV, "RESTART")
    restart!(ENV["RESTART"]; allow = ENVEUL === nothing ? () : ("PSWITCH",))
    RAD && Int(TF.getscalar(:lrad)) != 1 && error("RAD=1, but checkpoint $(ENV["RESTART"]) has no radiation field: build it with RAD=1 RADFROM=<checkpoint>")
    RAD && (radsettings!(); TD.EOS[](); TD.OPAC[]())
    if ENVEUL !== nothing
        Int(TF.getscalar(:lgrid)) == 3 ? envbuild!() : (envflags!(); ENVMS0[] = parse(Float64, ENV["ENVMS0"]) * Msun / 4π)
    end
elseif RAD
    # RAD=1 switches the CIE continuum cooling off; without a rebuild there is no radiation field to take its place
    error("RAD=1 from the initial condition: run with RAD=0 first, then RAD=1 RADFROM=<checkpoint> (scripts/popiii_run/fullrun.sh)")
end
"""
ENVSPLIT=r1,r2,k (radii in R☉; with RESTART of a fixed-mesh run): split every zone whose centre lies in [r1, r2] into k zones of equal Δln r,
and the zones beyond either end into k−1, k−2, … 2 (graded), keeping the zone count even. The rest of the mesh is unchanged. Conservative
piecewise-constant remap: each sub-zone takes its parent's ρ, T, E_r, composition, nuclear abundances and per-zone material and driver state
(mass, internal and radiation energy unchanged); u and the luminosity are linear in r inside the parent. Restart the split run's own
checkpoints WITHOUT ENVSPLIT (it would split again).
"""
function envsplit!(r1, r2, kmax)
    kmax = Int(kmax); r0 = TF.var(:rn)[a:b+1]; n0 = nsh; a0 = a; b0 = b
    rc = [0.5 * (r0[j] + r0[j+1]) for j in 1:n0]
    inw = [r1 * Rsun <= rc[j] <= r2 * Rsun for j in 1:n0]; any(inw) || error("ENVSPLIT: no zone centre in [$r1, $r2] R_sun")
    j1 = findfirst(inw); j2 = findlast(inw)
    f = [max(1, kmax - (j < j1 ? j1 - j : j > j2 ? j - j2 : 0)) for j in 1:n0]
    if isodd(sum(f) - n0)                                                       # (an even zone count, as envbuild!)
        jx = findfirst(j -> j > j2 && f[j] == 1, 1:n0); jx === nothing && error("ENVSPLIT: window too close to the outer edge")
        f[jx] = 2
    end
    u0 = TF.var(:un)[a:b+1]; L0 = 4π .* r0 .^ 2 .* TF.var(:frn)[a:b+1]
    d0 = TF.var(:dn)[a:b]; T0 = TF.var(:tn)[a:b]; E0 = TF.var(:ern)[a:b]
    rI = [r0[1]]; uI = [u0[1]]; LI = [L0[1]]; P = Int[]
    for j in 1:n0
        for i in 1:f[j]
            w = i / f[j]; r_ = f[j] == 1 ? r0[j+1] : exp((1 - w) * log(r0[j]) + w * log(r0[j+1]))
            i == f[j] && (r_ = r0[j+1]); wl = (r_ - r0[j]) / (r0[j+1] - r0[j])
            push!(rI, r_); push!(uI, (1 - wl) * u0[j] + wl * u0[j+1]); push!(LI, (1 - wl) * L0[j] + wl * L0[j+1]); push!(P, j)
        end
    end
    n = length(P); n + 2 <= TF.MGR || error("ENVSPLIT: $n zones exceed mgr = $(TF.MGR)")
    Fr = [rI[j] > 0 ? LI[j] / (4π * rI[j]^2) : 0.0 for j in eachindex(rI)]
    # per-zone state of the old mesh (TITAN index k = a + j − 1; ghosts beyond b keep their offset from b)
    kmap(knew) = knew < a0 ? knew : knew <= a0 + n - 1 ? a0 + P[knew-a0+1] - 1 : b0 + (knew - (a0 + n - 1))
    remap(v) = (o = copy(v); [o[clamp(kmap(k), 1, length(o))] for k in eachindex(o)])
    mx = deepcopy(mat.x); meq = copy(mat.eqzone); mΛ = copy(mat.Λext); mTr = copy(mat.Tref); msl = copy(mat.slope)
    lagF = Dict(:GASQ => TX.GASQ[], :RAYQ => TX.RAYQ[], :CHEMB => TX.CHEMB[], :CHEMBON => TX.CHEMBON[])
    lagold = Dict(s_ => copy(v_) for (s_, v_) in lagF)
    nlo = copy(PM.NLOPC[])
    xmext = TF.getscalar(:xmext); xmt0 = TF.getscalar(:xmtot)
    t0 = TF.getscalar(:timen); dt0 = TF.getscalar(:dtime); js = TF.getscalar(:jstep)
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    ρz = [d0[P[j]] for j in 1:n]; Tz = [T0[P[j]] for j in 1:n]; Er = [E0[P[j]] for j in 1:n]
    Mold = sum(4π / 3 * d0[j] * (r0[j+1]^3 - r0[j]^3) for j in 1:n0); Mnew = sum(4π / 3 * ρz[j] * (rI[j+1]^3 - rI[j]^3) for j in 1:n)
    Eold = sum(4π / 3 * E0[j] * (r0[j+1]^3 - r0[j]^3) for j in 1:n0); Enew = sum(4π / 3 * Er[j] * (rI[j+1]^3 - rI[j]^3) for j in 1:n)
    TX.GASROW[] = true
    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lum = ENVCENTRE ? nothing : abs(LI[1]), lgrid = 2, inner = :wall, outer = :wall,
              Mc = ENVCENTRE ? 0.0 : 4π * (xmt0 - TF.var(:xmen)[a0]) , selfgrav = true, material = mat, av, xmext)
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    a == a0 || error("ENVSPLIT: ngrs moved")
    for k in eachindex(mat.x)
        ko = clamp(kmap(k), 1, length(mx))
        mat.x[k] = copy(mx[ko]); mat.eqzone[k] = meq[ko]; mat.Λext[k] = mΛ[ko]; mat.Tref[k] = mTr[ko]; mat.slope[k] = msl[ko]
    end
    for (s_, v_) in lagF                                                        # (in place: PM.RAY etc. hold these arrays)
        o = lagold[s_]; isempty(o) && continue
        for k in eachindex(v_); v_[k] = o[clamp(kmap(k), 1, length(o))]; end
    end
    length(nlo) == length(TF.var(:dn)) && (PM.NLOPC[] = Any[nlo[clamp(kmap(k), 1, length(nlo))] for k in eachindex(nlo)])
    if HOTON
        g_ = SG[]; Yo = copy(g_.Y); Ypo = copy(g_.Yprev)
        for k in axes(g_.Y, 2); ko = clamp(kmap(k), 1, size(Yo, 2)); g_.Y[:, k] .= Yo[:, ko]; g_.Yprev[:, k] .= Ypo[:, ko]; end
        for v_ in (g_.fscale, g_.fslope, g_.fcur); o = copy(v_); for k in eachindex(v_); v_[k] = o[clamp(kmap(k), 1, length(o))]; end; end
        g_.Dmix = zeros(length(TF.var(:dn))); g_.acc.nstage = 0
        SM.ghosts!(g_)
    end
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    PS.install_julia!(BlockSolver)
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, t0); TF.setscalar!(:timeo, t0); TF.setscalar!(:dtime, min(dt0, parse(Float64, get(ENV, "ENVDT0", "1e4")))); TF.setscalar!(:jstep, js); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    envflags!()
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    empty!(ALIFAC[]); empty!(ALIFIX[])
    pr("ENVSPLIT: zones %d -> %d (old zones %d..%d with centres in %.4g..%.4g R_sun into %d; graded %s), M %.6e -> %.6e Msun, E_rad %.6e -> %.6e erg, xmtot %.6e -> %.6e\n",
       n0, n, j1, j2, rc[j1] / Rsun, rc[j2] / Rsun, kmax, string(f[max(1, j1 - kmax):min(n0, j2 + kmax)]), Mold / Msun, Mnew / Msun, Eold, Enew, xmt0, TF.getscalar(:xmtot))
end
"""
ENVMERGE=r1,r2,k (radii in R☉; with RESTART of a fixed-mesh run): the inverse of ENVSPLIT. Consecutive zones whose centres lie in [r1, r2]
are merged in groups of k (the last group may be smaller; one group gives a zone back if the zone count would be odd). The rest of the mesh is
unchanged, so the zone width jumps by k at the two ends of the window. Conservative remap: a merged zone has the mass, the radiation energy
and the gas internal energy of its parents (ρ = M/V, E_r volume-weighted, composition and nuclear abundances mass-weighted, T from the
equation of state such that M e(ρ, T, x) = Σ m e); the per-zone material and driver state comes from the heaviest parent; u and the luminosity
are kept at the interfaces that remain (kinetic energy is not conserved). Restart the merged run's own checkpoints WITHOUT ENVMERGE.
"""
function envmerge!(r1, r2, kmax)
    kmax = Int(kmax); r0 = TF.var(:rn)[a:b+1]; n0 = nsh; a0 = a; b0 = b
    rc = [0.5 * (r0[j] + r0[j+1]) for j in 1:n0]
    inw = [r1 * Rsun <= rc[j] <= r2 * Rsun for j in 1:n0]; any(inw) || error("ENVMERGE: no zone centre in [$r1, $r2] R_sun")
    j1 = findfirst(inw); j2 = findlast(inw)
    G = UnitRange{Int}[]; j = 1
    while j <= n0
        e_ = j1 <= j <= j2 ? min(j + kmax - 1, j2) : j
        push!(G, j:e_); j = e_ + 1
    end
    if isodd(length(G))                                                         # (an even zone count, as envbuild!)
        i = findlast(g -> length(g) > 1, G); i === nothing && error("ENVMERGE: nothing to merge")
        g = G[i]; splice!(G, i, [first(g):last(g)-1, last(g):last(g)])
    end
    n = length(G)
    u0 = TF.var(:un)[a:b+1]; L0 = 4π .* r0 .^ 2 .* TF.var(:frn)[a:b+1]
    d0 = TF.var(:dn)[a:b]; T0 = TF.var(:tn)[a:b]; E0 = TF.var(:ern)[a:b]; e0 = TF.var(:egn)[a:b]
    V0 = [(r0[j+1]^3 - r0[j]^3) / 3 for j in 1:n0]; m0 = d0 .* V0
    rI = vcat([r0[first(g)] for g in G], r0[end]); uI = vcat([u0[first(g)] for g in G], u0[end]); LI = vcat([L0[first(g)] for g in G], L0[end])
    Mg = [sum(m0[g]) for g in G]; Vg = [sum(V0[g]) for g in G]
    ρz = Mg ./ Vg; Er = [sum(E0[g] .* V0[g]) for g in G] ./ Vg
    Tz = [sum(T0[g] .* m0[g]) for g in G] ./ Mg; et = [sum(e0[g] .* m0[g]) for g in G] ./ Mg
    P = [g[argmax(m0[g])] for g in G]                                          # heaviest parent: per-zone material and driver state
    Fr = [rI[j] > 0 ? LI[j] / (4π * rI[j]^2) : 0.0 for j in eachindex(rI)]
    kmap(knew) = knew < a0 ? knew : knew <= a0 + n - 1 ? a0 + P[knew-a0+1] - 1 : b0 + (knew - (a0 + n - 1))
    mx = deepcopy(mat.x); meq = copy(mat.eqzone); mΛ = copy(mat.Λext); mTr = copy(mat.Tref); msl = copy(mat.slope)
    lagF = Dict(:GASQ => TX.GASQ[], :RAYQ => TX.RAYQ[], :CHEMB => TX.CHEMB[], :CHEMBON => TX.CHEMBON[])
    lagold = Dict(s_ => copy(v_) for (s_, v_) in lagF)
    nlo = copy(PM.NLOPC[])
    xmext = TF.getscalar(:xmext); xmt0 = TF.getscalar(:xmtot)
    t0 = TF.getscalar(:timen); dt0 = TF.getscalar(:dtime); js = TF.getscalar(:jstep)
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    Mold = 4π * sum(m0); Eold = 4π * sum(E0 .* V0); Egold = 4π * sum(e0 .* m0)
    Mc = ENVCENTRE ? 0.0 : 4π * (xmt0 - TF.var(:xmen)[a0])
    TX.GASROW[] = true
    # per-zone material on the new indices before build! fills the EOS (else a zone far out can inherit the table-EOS flag of the
    # old zone with its index: "rho below the table" with a 3-zone merge); set again below, with the merged compositions
    for k in eachindex(mat.x); ko = clamp(kmap(k), 1, length(mx)); mat.x[k] = copy(mx[ko]); mat.eqzone[k] = meq[ko]; end
    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lum = ENVCENTRE ? nothing : abs(LI[1]), lgrid = 2, inner = :wall, outer = :wall,
              Mc, selfgrav = true, material = mat, av, xmext)
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    a == a0 || error("ENVMERGE: ngrs moved")
    for k in eachindex(mat.x)
        ko = clamp(kmap(k), 1, length(mx))
        mat.x[k] = copy(mx[ko]); mat.eqzone[k] = meq[ko]; mat.Λext[k] = mΛ[ko]; mat.Tref[k] = mTr[ko]; mat.slope[k] = msl[ko]
    end
    for (j, g) in enumerate(G)                                                  # mass-weighted composition of the merged zones
        length(g) > 1 || continue
        mat.x[a+j-1] = sum(m0[i] .* mx[a0+i-1] for i in g) ./ Mg[j]
    end
    for (s_, v_) in lagF                                                        # (in place: PM.RAY etc. hold these arrays)
        o = lagold[s_]; isempty(o) && continue
        for k in eachindex(v_); v_[k] = o[clamp(kmap(k), 1, length(o))]; end
    end
    length(nlo) == length(TF.var(:dn)) && (PM.NLOPC[] = Any[nlo[clamp(kmap(k), 1, length(nlo))] for k in eachindex(nlo)])
    if HOTON
        g_ = SG[]; Yo = copy(g_.Y); Ypo = copy(g_.Yprev)
        for k in axes(g_.Y, 2); ko = clamp(kmap(k), 1, size(Yo, 2)); g_.Y[:, k] .= Yo[:, ko]; g_.Yprev[:, k] .= Ypo[:, ko]; end
        for (j, g) in enumerate(G)
            length(g) > 1 || continue
            g_.Y[:, a+j-1] .= sum(m0[i] .* Yo[:, a0+i-1] for i in g) ./ Mg[j]; g_.Yprev[:, a+j-1] .= sum(m0[i] .* Ypo[:, a0+i-1] for i in g) ./ Mg[j]
        end
        for v_ in (g_.fscale, g_.fslope, g_.fcur); o = copy(v_); for k in eachindex(v_); v_[k] = o[clamp(kmap(k), 1, length(o))]; end; end
        g_.Dmix = zeros(length(TF.var(:dn))); g_.acc.nstage = 0
        SM.ghosts!(g_)
    end
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    PS.install_julia!(BlockSolver)
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, t0); TF.setscalar!(:timeo, t0); TF.setscalar!(:dtime, min(dt0, parse(Float64, get(ENV, "ENVDT0", "1e4")))); TF.setscalar!(:jstep, js); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    envflags!()
    # T of the merged zones from the gas internal energy (secant in T; e can be negative, its zero is the neutral atomic gas at T = 0)
    tn = TF.var(:tn); eg = TF.var(:egn); jm = [j for j in 1:n if length(G[j]) > 1]
    TD.EOS[]()
    Tp = [tn[a+j-1] for j in jm]; ep = [eg[a+j-1] for j in jm]; dTmax = 0.0; emax = 0.0
    for (i, j) in enumerate(jm); tn[a+j-1] = Tp[i] * 1.001; end
    for it in 1:40
        TD.EOS[](); emax = 0.0
        for (i, j) in enumerate(jm)
            k = a + j - 1; T_ = tn[k]; e_ = eg[k]
            T_ == Tp[i] && continue
            cv = (e_ - ep[i]) / (T_ - Tp[i]); cv > 0 || (cv = abs(e_) / T_ + 1e7)
            dT = (et[j] - e_) / cv; emax = max(emax, abs(dT) / T_)
            abs(dT) < 1e-13 * T_ && continue
            Tp[i] = T_; ep[i] = e_; tn[k] = clamp(T_ + dT, 0.7 * T_, 1.4 * T_)
        end
        emax < 1e-11 && break
    end
    for j in jm; dTmax = max(dTmax, abs(tn[a+j-1] / Tz[j] - 1)); end
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    empty!(ALIFAC[]); empty!(ALIFIX[])
    dn = TF.var(:dn); rn = TF.var(:rn); ern = TF.var(:ern)
    Vn = [(rn[a+j]^3 - rn[a+j-1]^3) / 3 for j in 1:n]
    Mnew = 4π * sum(dn[a+j-1] * Vn[j] for j in 1:n); Enew = 4π * sum(ern[a+j-1] * Vn[j] for j in 1:n); Egnew = 4π * sum(eg[a+j-1] * dn[a+j-1] * Vn[j] for j in 1:n)
    pr("ENVMERGE: zones %d -> %d (old zones %d..%d with centres in %.4g..%.4g R_sun in groups of %d), M %.9e -> %.9e Msun, E_rad %.9e -> %.9e erg, E_gas %.9e -> %.9e erg (last T correction %.1e, max |T/T_massweighted - 1| %.2e), xmtot %.9e -> %.9e\n",
       n0, n, j1, j2, rc[j1] / Rsun, rc[j2] / Rsun, kmax, Mold / Msun, Mnew / Msun, Eold, Enew, Egold, Egnew, emax, dTmax, xmt0, TF.getscalar(:xmtot))
end
"""
DTENV=r1,r2,C (R☉): cap Δt at C × min over the interfaces with r1 ≤ r ≤ r2 of the dynamical time (r³/Gm)^½ (m interior mass). TITAN's
DT_DYN (CDYN) takes the minimum over the whole mesh, so the stellar core sets it and the accretion runs relax it (CDYN = 1e5); this one
follows only the envelope (its edge collapses and re-expands on its dynamical time). Reported on the DTDIAG line.
"""
const DTENV = haskey(ENV, "DTENV") ? parse.(Float64, split(ENV["DTENV"], ",")) : nothing; const DTENVLAST = Ref((Inf, 0))
if DTENV !== nothing
    push!(TC.DT_HOOKS, function ()
        rn = TF.var(:rn); xmn = TF.var(:xmn); G4π = TF.getscalar(:cgrav) * 4π; td = Inf; kd = 0
        for k in a:b+1
            (DTENV[1] * Rsun <= rn[k] <= DTENV[2] * Rsun && xmn[k] > 0) || continue
            t_ = sqrt(rn[k]^3 / (G4π * xmn[k])); t_ < td && (td = t_; kd = k)
        end
        DTENVLAST[] = (DTENV[3] * td, kd); DTENV[3] * td
    end)
    pr("DTENV: dt <= %g x dynamical time over %g..%g R_sun\n", DTENV[3], DTENV[1], DTENV[2])
end
haskey(ENV, "ENVSPLIT") && ENVEUL !== nothing && haskey(ENV, "RESTART") && Int(TF.getscalar(:lgrid)) == 2 && envsplit!(parse.(Float64, split(ENV["ENVSPLIT"], ","))...)
haskey(ENV, "ENVMERGE") && ENVEUL !== nothing && haskey(ENV, "RESTART") && Int(TF.getscalar(:lgrid)) == 2 && for w_ in split(ENV["ENVMERGE"], ";"); envmerge!(parse.(Float64, split(w_, ","))...); end
"""
TRESET=k:T[,k:T...] (with RESTART; k = zone number counted from the first envelope zone, as in the recorded snapshots; T in K): set the gas
temperature of the listed zones once, right after the checkpoint is read. Density, velocity, radiation field and composition are untouched;
the gas energy and pressure follow from the equation of state, and the zone's held line-cooling rate is cleared (the transfer refills it
after the first step). This is a repair of zones whose temperature was changed by a discretization artefact (the total-energy row in the
transparent outer envelope, see docs/POPIII_RUNS.md): it is not conservative, and the gas energy it adds is printed. Use it on one
restart only; the checkpoints written afterwards are restarted WITHOUT it.
"""
function treset!(spec)
    tn = TF.var(:tn); eg = TF.var(:egn); dn = TF.var(:dn); rn = TF.var(:rn)
    TD.EOS[]()
    L = Tuple{Int,Float64,Float64,Float64}[]
    for it in split(spec, ",")
        ks, Ts = split(it, ":"); k = a + parse(Int, ks) - 1; Tnew = parse(Float64, Ts)
        a <= k <= b || error("TRESET: zone $ks is not on the mesh (1..$nsh)")
        Tnew > 0 || error("TRESET: T must be positive")
        push!(L, (k, tn[k], Tnew, eg[k])); tn[k] = Tnew
        mat.Λext[k] = 0.0; mat.Tref[k] = Tnew; mat.slope[k] = 0.0
    end
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    dE = 0.0
    for (k, Told, Tnew, eold) in L
        m_ = 4π * dn[k] * (rn[k+1]^3 - rn[k]^3) / 3; dE += m_ * (eg[k] - eold)
        pr("TRESET: zone %d (r = %.4g..%.4g pc, %.4g Msun): T %.2f -> %.2f K, gas energy %+.3e erg\n", k - a + 1, rn[k] / 3.0857e18, rn[k+1] / 3.0857e18,
           m_ / Msun, Told, Tnew, m_ * (eg[k] - eold))
    end
    pr("TRESET: %d zones, gas energy added %+.4e erg\n", length(L), dE)
end
haskey(ENV, "TRESET") && haskey(ENV, "RESTART") && treset!(ENV["TRESET"])
"""
GRIDWIN=r1,r2[,tau[,wd,wt[,alph[,ctol[,fmove[,gmax]]]]]] (R☉; fixed mesh only): the interfaces nearest r1 and r2 stay fixed, the ones in between
follow TITAN's adaptive grid equation (TitanResidFast.GRIDWIN: equal arc length in (ln r, √wd ln ρ, √wt ln T), spatial smoothing alph,
time smoothing tau [s]). Defaults tau 1e7, wd 0.3, wt 1, alph 1.5, ctol 0.3 (maxdel: relative zone-width change per Newton iteration),
fmove 0.3 (Δt such that no window interface moves more than fmove × its smaller neighbouring zone per step), gmax Inf
(cap on the resolved gradient: zones at most √(1+gmax²) finer than in smooth gas). Chemistry and nuclear
abundances already move with the mass crossing each interface (advect_species!); the total-energy row's gravity work uses the lab-frame
mass flux. Restart a window checkpoint with the same r1, r2 (the end interfaces never move, so they are found again).
"""
const GRIDWINP = haskey(ENV, "GRIDWIN") ? parse.(Float64, split(ENV["GRIDWIN"], ",")) : nothing
const GWLAST = Ref((0.0, 0, Inf))     # (max displacement / zone width over the last step, interface, Δt cap)
if GRIDWINP !== nothing
    Int(TF.getscalar(:lgrid)) == 2 || error("GRIDWIN needs the fixed mesh (lgrid = 2)")
    let p = GRIDWINP, rn = TF.var(:rn), gwp(i, d) = length(p) >= i ? p[i] : d
        k1 = a - 1 + argmin(abs.(rn[a:b+1] .- p[1] * Rsun)); k2 = a - 1 + argmin(abs.(rn[a:b+1] .- p[2] * Rsun))
        (k1 - 1 > a && k2 + 1 < b + 1 && k2 - k1 >= 4) || error("GRIDWIN: window $k1..$k2 too small or too close to the mesh ends")
        TX.GRIDWIN[] = (; k1, k2, tau = gwp(3, 1e7), wd = gwp(4, 0.3), wt = gwp(5, 1.0), alph = gwp(6, 1.5), gmax = gwp(9, Inf))
        TC.GRIDMOVE[] = true; TF.setscalar!(:ctol, gwp(7, 0.3)); fmove = gwp(8, 0.3)
        push!(TC.DT_HOOKS, function ()
            rn_ = TF.var(:rn); ro_ = TF.var(:ro); f = 0.0; kf = 0
            for k in k1+1:k2-1
                w = min(ro_[k] - ro_[k-1], ro_[k+1] - ro_[k]); q = abs(rn_[k] - ro_[k]) / w; q > f && (f = q; kf = k)
            end
            cap = f > 0 ? fmove / f * TF.getscalar(:dtime) : Inf
            GWLAST[] = (f, kf, cap); cap
        end)
        pr("GRIDWIN: adaptive interfaces %d..%d (r %.4g..%.4g R_sun, %d zones), tau %.3g s, weights ln rho %.3g ln T %.3g, alph %.3g, ctol %.3g, fmove %.3g, gmax %.3g\n",
           k1 + 1, k2 - 1, rn[k1] / Rsun, rn[k2] / Rsun, k2 - k1, TX.GRIDWIN[].tau, TX.GRIDWIN[].wd, TX.GRIDWIN[].wt, TX.GRIDWIN[].alph, TF.getscalar(:ctol), fmove, TX.GRIDWIN[].gmax)
    end
end
if GRIDMAPP !== nothing
    (ENVEUL !== nothing && Int(TF.getscalar(:lgrid)) == 2) || error("GRIDMAP needs the fixed-mesh envelope (ENVEUL, lgrid = 2)")
    GRIDWINP === nothing || error("GRIDMAP and GRIDWIN exclude each other")
    TX.GRIDVEL[] = zeros(length(TF.var(:rn))); TC.GRIDMOVE[] = true; TF.setscalar!(:ctol, gmp(12, 0.5))
    pr("GRIDMAP: %d zones follow the target mesh (hr %.3g, h1 %.3g, h2 %.3g, eps %.3g, fc %.3g, dmin %.3g, g %.3g), tau %.3g s, fmove %.3g, nrelax %.3g, ctol %.3g, ftar %.3g\n",
       nsh, gridmap_kw()..., gmp(9, 0.0), gmp(10, 0.4), gmp(11, 3.0), TF.getscalar(:ctol), gmp(13, 0.2))
end
"""
GRIDMAP: interface velocities for the next step from the target mesh of the present state (after the chemistry step, when Δt of the
next step is known). The mesh closes 1/nrelax of its distance from the target per step (at most fmove zone widths; not faster than
in the time tau), and Δt is limited so that the target itself moves at most ftar zone widths per step, judged from its move in the
step just taken: with ftar < fmove the mesh stays within about ftar × nrelax zone widths of its target.
"""
function gridmap_step!()
    rn = TF.var(:rn); i0 = a + (ENVCENTRE ? 1 : 0); N = b + 1 - i0
    r = rn[i0:b+1]; m = envmenc()[i0-a+1:end]; u = TF.var(:un)[i0:b+1]
    Mst = GridMap.anchor(m, u); P_ = GridMap.target(r, m, Mst; N, lam0 = GMLAST[][4], gridmap_kw()...)
    dt = TF.getscalar(:dtime); dtl = TC.TSTP_IN[][2]; qT = 0.0; cap = Inf; Pp = GMPREV[]
    if length(Pp) == N + 1 && dtl > 0
        for j in 2:N; qT = max(qT, abs(P_.r[j] - Pp[j]) / min(r[j] - r[j-1], r[j+1] - r[j])); end
        qT > 0 && (cap = max(dtl * gmp(13, 0.2) / qT, 0.3 * dtl))
        cap < dt && (dt = cap; TF.setscalar!(:dtime, dt); TF.setscalar!(:tfac, dt / dtl))
    end
    GMPREV[] = P_.r; trel = max(gmp(9, 0.0), gmp(11, 3.0) * dt); fmove = gmp(10, 0.4)
    w = TX.GRIDVEL[]; fill!(w, 0.0); f = 0.0; kf = 0; dmax = 0.0
    for j in 2:N
        wd = min(r[j] - r[j-1], r[j+1] - r[j]); dmax = max(dmax, abs(P_.r[j] - r[j]) / wd)
        wk = (P_.r[j] - r[j]) / trel; q = abs(wk) * dt / wd
        q > fmove && (wk *= fmove / q; q = fmove)
        w[i0+j-1] = wk; q > f && (f = q; kf = i0 + j - 1)
    end
    GMLAST[] = (f, kf, Mst, P_.lam, P_.nfloor, dmax, qT, cap)
    nothing
end
"Mesh state for the progress print."
function gridmap_line()
    rn = TF.var(:rn); i0 = a + (ENVCENTRE ? 1 : 0); dl = [log(rn[j+1] / rn[j]) for j in i0:b]; jm = argmin(dl); f, kf, Mst, lam, nfl, dmax, qT, cap = GMLAST[]
    @sprintf("        grid map: M* %.5f Msun, scale %.3f, %d zones at the floor; narrowest zone dln r %.2e at %.1f R_sun; mesh at most %.2f zone widths from its target, next step moves %.3f (interface %d); target moved %.3f in the last step, dt cap %.2e s\n",
             4π * Mst / Msun, lam, nfl, dl[jm], rn[i0+jm-1] / Rsun, dmax, f, kf, qT, cap)
end
"Window state for the progress print: narrowest zone (Δ ln r, r), widest, last-step displacement."
function gridwin_line()
    gw = TX.GRIDWIN[]; rn = TF.var(:rn)
    dl = [log(rn[j+1] / rn[j]) for j in gw.k1:gw.k2-1]; jm = argmin(dl); jx = argmax(dl)
    f, kf, cap = GWLAST[]
    @sprintf("        grid window: dln r min %.2e at %.0f R_sun (zone %d), max %.2e at %.0f; last step moved %.3f zone widths (interface %d, dt cap %.2e)\n",
             dl[jm], rn[gw.k1+jm-1] / Rsun, gw.k1 + jm - 1, dl[jx], rn[gw.k1+jx-1] / Rsun, f, kf, cap)
end
const CHEMDIAG = get(ENV, "CHEMDIAG", "0") == "1"; const LNUC = Ref(0.0); const TIMING = get(ENV, "TIMING", "0") == "1"; const TSTEP = Ref(0.0); const TCHEM = Ref(0.0); const TLINE = Ref(0.0); const TOTH = Ref(0.0)
const TMARK = Ref(time()); const TOUT = Ref(0.0)
# MSTOP [M☉]: stop when the mass inside the first interface falling faster than 2 km/s (the star, as REZWIN and the performance
# record define it) exceeds this: the hand-over from the Lagrangian core-growth stage to the accretion mesh (scripts/popiii_run)
const MSTOP = parse(Float64, get(ENV, "MSTOP", "0"))
const DTSTOP = parse(Float64, get(ENV, "DTSTOP", "0"))   # stop an unattended run whose Δt has collapsed, before the checkpoints are overwritten
const PRINTEVERY = parse(Int, get(ENV, "PRINTEVERY", "500")); const CKEVERY = parse(Int, get(ENV, "CKEVERY", "500")); const KLAST = Ref(0)
# Newton controls, applied after any restart: NITER (iterations per attempt), FAILTFAC (Δt factor after a failed attempt)
haskey(ENV, "NITER") && TF.setscalar!(:niter, parse(Float64, ENV["NITER"]))
haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))
# STOLCHEM=1: the step-size rule (stol on the relative change of r, m, rho, T, E_r, and the caps in DT_HOOKS) is evaluated again after
# the chemistry step, so that it acts on the change over the whole step (transport + chemistry) and not on the state between the two
# halves of the operator split, where the transport step has moved T at the old composition. Default off.
const STOLCHEM = get(ENV, "STOLCHEM", "0") == "1"
haskey(ENV, "NTRY") && TF.setscalar!(:ntry, parse(Float64, ENV["NTRY"]))
if haskey(ENV, "FAILTFAC"); TD.FAIL_TFAC[] = parse(Float64, ENV["FAILTFAC"]); Trbdf2.TFAC[] = TD.FAIL_TFAC[]; end
# EDDREAL=1: the moments of any radiation field obey K/J >= (H/J)^2, so f >= (F/cE)^2. Applied after the formal solution (f is lagged over the step anyway), with F, E of the current state: in thin
# zones f from the formal solution can be too low for the E the F row integrates from it (E too low by 3, F/cE = 3, with Eddington factor 0.6). No effect in opaque zones, where (F/cE)^2 << 1/3.
if get(ENV, "EDDREAL", "0") == "1"
    let ef0 = TC.EDDFAC[], EDDTAU = parse(Float64, get(ENV, "EDDTAU", "0.3")), EDDRELAX = parse(Float64, get(ENV, "EDDRELAX", "1")), prev = zeros(length(TF.var(:fedd)))
        TC.EDDFAC[] = function ()
            ef0()
            fe = TF.var(:fedd); fr_ = TF.var(:frn); er_ = TF.var(:ern); cl = TF.getscalar(:cc)
            chi_ = TF.var(:chifn); d_ = TF.var(:dn); r_ = TF.var(:rn); tau_ = 0.0     # only above tau = EDDTAU: at the unresolved photosphere F/cE across one zone is a grid effect, and f = 1 there would destroy diffusion
            for k in Int(TF.getscalar(:ngre)):-1:Int(TF.getscalar(:ngrs))
                tau_ += chi_[k] * d_[k] * (r_[k+1] - r_[k])
                tau_ > EDDTAU && break
                x = max(fr_[k], fr_[k+1], 0.0) / (cl * er_[k])
                raw = fe[k]; lim = max(raw, min(x * x, 1.0))
                prev[k] > 0 && (lim = prev[k] + EDDRELAX * (lim - prev[k]))      # under-relaxed over steps
                fe[k] = max(raw, lim); prev[k] = fe[k]
            end
        end
    end
end
haskey(ENV, "ERTHICK") && (TC.ER_THICK[] = parse(Float64, ENV["ERTHICK"]))
# ETOTAL=1: TITAN's conservative total-energy row (P4.2); the gas row stays where e/(e + KE) < ETASW, in the cold hypersonic infall, where T is a small difference
haskey(ENV, "ETASWTH") && (TX.ETASW_TH[] = parse(Float64, ENV["ETASWTH"]))     # internal-energy row also where 3p/2ρ < ETASWTH × (3p/2ρ + E_r/ρ + K/m) (off: as before)
haskey(ENV, "ETASWTIE") && (TX.ETASW_TIE[] = parse(Float64, ENV["ETASWTIE"]))   # internal-energy row also where the FLUXTIE weight at the zone's inner interface >= ETASWTIE (off: as before)
get(ENV, "ETOTGRAVGRID", "0") == "1" && (TX.ETOT_GRAV_GRID[] = true)     # gravity work on Lagrangian zones in the total-energy row (off: as before)
get(ENV, "ETOTAL", "0") == "1" && (TX.ETOTAL[] = true; TX.ETOT_GRAV[] = true; TX.ETASW_EOFF[] = PC.UOFF_TH; TX.ETASW[] = parse(Float64, get(ENV, "ETASW", "0.1")))
haskey(ENV, "EXCAP") && (TX.EXCAP[] = parse(Float64, ENV["EXCAP"]))
const LINEFRAC = parse(Float64, get(ENV, "LINEFRAC", "0")); const LINELAST = [-Inf, 0.0]
# NEWTKX=1: the largest Newton correction of each equation at every iteration from the third on, with the zone it is in (debugging)
if get(ENV, "NEWTKX", "0") == "1"
    let md0 = TD.ROUTINES[:maxdel]
        TD.ROUTINES[:maxdel] = function ()
            md0()
            it_ = Int(TF.getscalar(:iter))
            it_ >= 3 && pr("KX dt %.3e iter %d dmax %.2e: %s\n", TF.getscalar(:dtime), it_, TF.getscalar(:dmax),
                           join(("$(TF.var(:kx)[i]):$(round(TF.var(:dx)[i], sigdigits = 2))" for i in 1:Int(TF.getscalar(:neqn))), " "))
        end
    end
end
const JACSAVE = Ref{Any}(nothing)
# JACCHK="k,z,iter": at step k, Newton iteration iter, finite differences of the residual against the assembled Jacobian for the variables of absolute zone z (debugging)
if haskey(ENV, "JACCHK")
    let w = parse.(Int, split(ENV["JACCHK"], ",")), mg0 = TD.ROUTINES[:matgen], done = Ref(false)
        ms0 = TD.ROUTINES[:matslv]
        TD.ROUTINES[:matslv] = function ()                 # residual of the linear solve just done, per row, for the saved system
            sv = JACSAVE[]
            ms0()
            sv === nothing && return
            JACSAVE[] = nothing
            x = TF.var(:rhs); worst = zeros(sv.ne); wz = zeros(Int, sv.ne)
            for kb in 1:sv.nb, i in 1:sv.ne
                k = sv.k0 + kb - 1; sm = 0.0; sc = abs(sv.rhs[i, k])
                for o in -2:2
                    1 <= kb + o <= sv.nb || continue
                    for j in 1:sv.ne; t = sv.E[o + 3][i, j, k] * x[j, k + o]; sm += t; sc = max(sc, abs(t)); end
                end
                r = min(abs(sm - sv.rhs[i, k]), abs(sm + sv.rhs[i, k])) / (sc + 1e-300)
                r > worst[i] && (worst[i] = r; wz[i] = k)
            end
            pr("SOLVE worst |A dx - b| / largest term per row: %s\n", join(("$i: $(round(worst[i], sigdigits = 2)) (zone $(wz[i]))" for i in 1:sv.ne), ", "))
        end
        TD.ROUTINES[:matgen] = function ()
            mg0()
            (done[] || KLAST[] + 1 < w[1] || Int(TF.getscalar(:iter)) != w[3]) && return
            done[] = true
            TX = TitanResidFast; c = TX.context(); ne = c.neqn; nb = c.nb; x0 = zeros(ne * nb); ws = TX.getws(Float64, c.mgr)
            Rp = zeros(ne, nb); Rm = zeros(ne, nb); R0 = zeros(ne, nb); TX.residual!(R0, ws, c, x0)
            E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2))
            z = w[2]; kbz = z - c.k0 + 1
            pr("JACCHK step %d iter %d zone %d dt %.3e (rows: 1 r 2 m 3 rho 4 u 5 T 6 E 7 F; code of each variable %s)\n", KLAST[] + 1, w[3], z, TF.getscalar(:dtime), string(c.code))
            for j in 1:ne
                h = 1e-5; x = copy(x0); ix = (kbz - 1) * ne + j
                x[ix] = h; TX.residual!(Rp, ws, c, x); x[ix] = -h; TX.residual!(Rm, ws, c, x)
                for kb in max(1, kbz - 3):min(nb, kbz + 3), i in 1:ne
                    fd = (Rp[i, kb] - Rm[i, kb]) / 2h; o = kbz - kb
                    ad = abs(o) <= 2 ? E[o + 3][i, j, c.k0 + kb - 1] : 0.0
                    scale = maximum(abs, (E[3][i, jj, c.k0 + kb - 1] for jj in 1:ne)) + 1e-300
                    abs(fd - ad) > 1e-3 * scale && pr("   row %d zone %d var %d (offset %d): FD %.5e AD %.5e, largest diagonal entry of the row %.3e\n", i, c.k0 + kb - 1, j, o, fd, ad, scale)
                end
            end
            pr("JACCHK end\n")
            JACSAVE[] = (; E = map(copy, E), rhs = copy(TF.var(:rhs)), ne, nb, k0 = c.k0)
        end
    end
end
# JACDUMP="k,iter,file": at step k, Newton iteration iter (first stage tried), the assembled block Jacobian, right-hand side and state to a file (debugging;
# every later matgen of that step is appended to the same list, so the stages and retries of the step can be compared)
if haskey(ENV, "JACDUMP")
    let w = split(ENV["JACDUMP"], ","), mg0 = TD.ROUTINES[:matgen], acc = Any[]
        ks = parse(Int, w[1]); its = parse(Int, w[2])
        TD.ROUTINES[:matgen] = function ()
            mg0()
            (KLAST[] + 1 == ks && Int(TF.getscalar(:iter)) == its) || return
            c = TitanResidFast.context()
            push!(acc, (; E = map(copy, (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2))), rhs = copy(TF.var(:rhs)), ne = c.neqn, nb = c.nb, k0 = c.k0,
                         code = copy(c.code), dtime = TF.getscalar(:dtime), thet = TF.getscalar(:thet), jtry = Int(TF.getscalar(:jtry)),
                         r = copy(TF.var(:rn)), d = copy(TF.var(:dn)), T = copy(TF.var(:tn)), u = copy(TF.var(:un)), p = copy(TF.var(:pgn)), m = copy(TF.var(:xmen)),
                         ro = copy(TF.var(:ro)), d_o = copy(TF.var(:do)), To = copy(TF.var(:to)), uo = copy(TF.var(:uo)), unom = copy(TF.var(:unom))))
            Serialization.serialize(String(w[3]), acc)
        end
    end
end
TC.EDD_SCATTER[] = get(ENV, "EDDSCAT", "0") == "1"
if haskey(ENV, "FLUXTIE")
    let v = parse.(Float64, split(ENV["FLUXTIE"], ","))
        TD.FLUXTIE[] = (v[1], v[2])
        pr("FLUXTIE: radiation momentum row -> F = h c E (h = H/J of the formal solution) where tau from the outer boundary <= %g, blended to the moment equation by tau = %g\n", v[1], v[2])
    end
end
TC.EDDG_MU0[] = get(ENV, "EDDGMU0", "0") == "1"
# EDDCHI=1: extinction of the formal solution at most kappa_P + electron scattering. It takes S from kappa_P and scattering; with chi >> kappa_P (thin network zones: chi uses LTE tables, 50-1000 x kappa_P) it would emit
# chi B where the energy rows emit kappa_P B, and that isotropic light lowers f (0.8 at 1000-25000 R_sun), so E comes out below F/c when the F row integrates it inward.
if get(ENV, "EDDCHI", "0") == "1"
    let ef0 = TC.EDDFAC[], sv = Float64[]
        TC.EDDFAC[] = function ()
            chi, kp, ne, d_ = TF.var(:chifn), TF.var(:xkpn), TF.var(:xnen), TF.var(:dn); cs = TF.getscalar(:csige)
            resize!(sv, length(chi)); copyto!(sv, chi)
            for k in Int(TF.getscalar(:ngrs)):Int(TF.getscalar(:ngre)); chi[k] = min(chi[k], kp[k] + cs * ne[k] / d_[k]); end
            try ef0() finally copyto!(chi, sv) end
        end
    end
end
# EDDDBG=k: at step k, the formal solution's J (as 4 pi J / c) and f beside the moment system's E, F, with the opacities and the two emissivities, zone by zone (debugging)
if haskey(ENV, "EDDDBG")
    let ks = parse(Int, ENV["EDDDBG"]), ef0 = TC.EDDFAC[], done = Ref(false)
        TC.EDDFAC[] = function ()
            ef0()
            (done[] || KLAST[] + 1 < ks) && return
            done[] = true
            xm0, xm2 = TF.var(:xmom0), TF.var(:xmom2); fe, E, Fr = TF.var(:fedd), TF.var(:ern), TF.var(:frn)
            chi, ke, kp, dn_, tn_, rn_, pl, ne = TF.var(:chifn), TF.var(:xken), TF.var(:xkpn), TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:plfn), TF.var(:xnen)
            cl = TF.getscalar(:cc); cpi = TF.getscalar(:cpi); cs = TF.getscalar(:csige); tau = 0.0
            for k in Int(TF.getscalar(:ngre)):-1:Int(TF.getscalar(:ngrs))
                tau += chi[k] * dn_[k] * (rn_[k+1] - rn_[k])
                tau > 3 && break
                pr("EDD %d r %.3e T %.0f tau %.2e chi %.3e kE %.3e kP %.3e sige/rho %.2e | 4piB/c %.3e 4piJ/c %.3e E %.3e F/c %.3e | f %.3f K/J %.3f\n", k - 4, rn_[k] / 6.957e10, tn_[k], tau, chi[k], ke[k], kp[k], cs * ne[k] / dn_[k],
                   4cpi * pl[k] / cl, 4cpi * xm0[k] / cl, E[k], max(Fr[k], Fr[k+1]) / cl, fe[k], xm2[k] / xm0[k])
            end
        end
    end
end
# NEWTSEC="k1,itmin,nmax": from step k1, for up to nmax iterations >= itmin: R(x+dx) - R(x) against J dx row by row (what the Jacobian misses, including offsets beyond +-2),
# and R(x+dx; context of x) against R(x+dx; recomputed context) (what is lagged in the context). Debugging.
if haskey(ENV, "NEWTSEC")
    let w = parse.(Int, split(ENV["NEWTSEC"], ",")), mg0 = TD.ROUTINES[:matgen], sv = Ref{Any}(nothing), nout = Ref(0)
        TD.ROUTINES[:matgen] = function ()
            mg0()
            TX = TitanResidFast
            it_ = Int(TF.getscalar(:iter))
            (KLAST[] + 1 < w[1] || nout[] >= w[3]) && (sv[] = nothing; return)
            E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2))
            c = TX.context(); ne = c.neqn; nb = c.nb; ws = TX.getws(Float64, c.mgr)
            if sv[] !== nothing && it_ == sv[].it + 1
                s = sv[]; Rnew = zeros(ne, nb); TX.residual!(Rnew, ws, c, zeros(ne * nb))
                basenew = copy(c.base); cvwnew = copy(TX.BUFS[:cvw])
                xa = zeros(ne * nb); nbad = Ref(0)
                for kb in 1:nb, j in 1:ne
                    k = c.k0 + kb - 1; cde = c.code[j]; cde == 0 && continue
                    xa[(kb-1)*ne+j] = cde == 4 ? (basenew[k,4] - s.base[k,4]) / c.nomu[k] : cde == 7 ? (basenew[k,7] - s.base[k,7]) / c.nomf[k] : log(basenew[k,cde] / s.base[k,cde])
                    isfinite(xa[(kb-1)*ne+j]) || (xa[(kb-1)*ne+j] = 0.0; nbad[] += 1)
                end
                c.base .= s.base; TX.BUFS[:cvw] .= s.cvw
                Rlag = zeros(ne, nb); TX.residual!(Rlag, ws, c, xa)
                c.base .= basenew; TX.BUFS[:cvw] .= cvwnew
                jerr = zeros(ne, nb); lerr = zeros(ne, nb); sc = zeros(ne, nb); Jd = zeros(ne, nb); rs = zeros(ne, nb)
                for kb in 1:nb, i in 1:ne
                    k = c.k0 + kb - 1; Jdx = 0.0; m = abs(s.R[i,kb])
                    for o in -2:2, j in 1:ne
                        1 <= kb + o <= nb || continue
                        t = s.E[o+3][i,j,k] * xa[(kb+o-1)*ne+j]; Jdx += t; m = max(m, abs(t))
                    end
                    sc[i,kb] = m
                    Jd[i,kb] = Jdx; rs[i,kb] = maximum(abs, (s.E[3][i,jj,k] for jj in 1:ne))
                    jerr[i,kb] = Rlag[i,kb] - s.R[i,kb] - Jdx        # nonlinearity + missing Jacobian
                    lerr[i,kb] = Rnew[i,kb] - Rlag[i,kb]             # lagged context
                end
                pr("SEC k %d iter %d dt %.3e |x|max %.2e (nonfinite %d) cvw change %.1e floors so far %d\n", KLAST[]+1, it_, TF.getscalar(:dtime), maximum(abs, xa), nbad[], maximum(abs, cvwnew .- s.cvw), TC.ESRC_FLOORS[])
                for kb in 1:nb
                    k = c.k0 + kb - 1; kx = k - 4
                    112 <= kx <= 126 || continue
                    pr("   %d x:%s\n", kx, join((@sprintf("%9.1e", xa[(kb-1)*ne+j]) for j in 1:ne), ""))
                    pr("      miss/Jdx:%s\n", join((@sprintf("%9.1e", jerr[i,kb] / (Jd[i,kb] == 0 ? 1.0 : abs(Jd[i,kb]))) for i in 1:ne), ""))
                    pr("      R/row   :%s\n", join((@sprintf("%9.1e", Rnew[i,kb] / (rs[i,kb] + 1e-300)) for i in 1:ne), ""))
                end
                nout[] += 1
            end
            if it_ >= w[2] - 1
                R0 = zeros(ne, nb); TX.residual!(R0, ws, c, zeros(ne * nb))
                sv[] = (; it = it_, E = map(copy, E), R = R0, base = copy(c.base), cvw = copy(TX.BUFS[:cvw]))
            else
                sv[] = nothing
            end
        end
    end
end
const NEWTFLIP_OLD = Dict{Tuple{Symbol,Int},Vector{Bool}}()
# NEWTDBG="k1,k2,z1,z2": per-iteration state of zones z1:z2 during steps k1:k2 (debugging)
if haskey(ENV, "NEWTDBG")
    let w = parse.(Int, split(ENV["NEWTDBG"], ",")), up0 = TD.ROUTINES[:update]
        TD.ROUTINES[:update] = function ()
            up0()
            k = KLAST[] + 1
            if w[1] <= k <= w[2]
                pr("NEWT k %d iter %d jtry %d dt %.4e dmax %.3e boost %d\n", k, Int(TF.getscalar(:iter)), Int(TF.getscalar(:jtry)), TF.getscalar(:dtime), TF.getscalar(:dmax), Int(TF.getscalar(:lboos)))
                if get(ENV, "NEWTFLIP", "0") == "1"       # which branch decisions of the residual changed since the last iteration (record mode changes nothing)
                    TitanResidFast.FREEZE[] = 1
                    for (key, v) in TitanResidFast.DEC
                        old = get!(() -> copy(v), NEWTFLIP_OLD, key)
                        fl = findall(old .!= v)
                        isempty(fl) || pr("   FLIP %s/%d zones %s\n", key[1], key[2], join(fl, ","))
                        old .= v
                    end
                end
                if Int(TF.getscalar(:iter)) == 1 && Int(TF.getscalar(:jtry)) == 1 && mat isa PM.NetGas
                    for z in w[3]:w[4]
                        ρz = TF.var(:dn)[z]; Tz = TF.var(:tn)[z]; ez = TF.var(:egn)[z]
                        _, Tt = PopIIIEqEOS.step(ρz, ez, Tz, 0.0; cooling = nothing)
                        xn = ChemNet.full(mat.x[z], mat.net)
                        utab = PopIIIEqEOS.eos(ρz, Tz)[1]; unet = ChemNet.eos(ρz, Tz, mat.x[z], mat.net)[1]
                        pr("   CMP z %d p %.3e T %.1f e %.5e | T_tab(e) %.1f u_tab(T) %.5e u_net(T) %.5e | xH+ net %.4e tab %.4e | xH2 %.3e\n", z, TF.var(:pgn)[z], Tz, ez, Tt, utab, unet, xn[:Hp], PopIIIEqEOS.composition(ρz, Tz)[1], xn[:H2])
                    end
                end
                if Int(TF.getscalar(:iter)) == 1 || Int(TF.getscalar(:lboos)) == 1
                    dualu = TitanResidFast.dualswitch(Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)), true, TitanResidFast.ETASW[])
                    pr("   U: %s\n", join(("$(z): $(round(TF.var(:un)[z]/1e5, sigdigits=5)) gasrow=$(dualu[z]) pg=$(round(TF.var(:pgn)[z], sigdigits=4)) pr=$(round(TF.var(:ern)[z]/3, sigdigits=4))" for z in w[3]:w[4]), " | "))
                end
                for z in w[3]:w[4]
                    pr("   z %d T %.1f To %.1f rho %.4e E %.4e Eo %.4e F %.4e chi %.3e fedd %.3f eq %d\n", z, TF.var(:tn)[z], TF.var(:to)[z], TF.var(:dn)[z], TF.var(:ern)[z], TF.var(:ero)[z], TF.var(:frn)[z], TF.var(:chifn)[z], TF.var(:fedd)[z], eqflag(z))
                end
            end
        end
    end
end
const RUNEND = Ref{Any}(nothing)
try
    RUNEND[] = TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = function (k0)
        k = k0 + KOFF[]; KLAST[] = k
        tnow = time(); TSTEP[] += tnow - TMARK[] - TOUT[]; TOUT[] = 0.0
        if haskey(ENV, "ZDIAG")                     # S6: energy terms of the hydro step per gas mass in zones z1..z2 (before chem!)
            zz = parse.(Int, split(ENV["ZDIAG"], ",")); dt_ = TF.getscalar(:timen) - TF.getscalar(:timeo)
            for z in zz[1]:zz[2]
                ρz = TF.var(:dn)[z]; Tz = TF.var(:tn)[z]; kp = TF.var(:xkpn)[z]; Ez = TF.var(:ern)[z]
                ke = TF.var(:xken)[z]
                xch = (ke * ChemNet.c * Ez - kp * 4 * 5.670374e-5 * Tz^4) * dt_            # (κ_E c E − κ_P 4σT⁴) Δt, the hydro exchange [erg/g]
                gq = isempty(TX.GASQ[]) ? 0.0 : TX.GASQ[][z] * dt_ / ρz; rq = isempty(TX.RAYQ[]) ? 0.0 : TX.RAYQ[][z] * dt_ / ρz
                pr("        ZDIAG %d: T %.1f (To %.1f) e %.5e (eo %.5e) rho %.3e (do %.3e) | per g over dt %.2e: kE cE-kP 4sT4 %.3e, GASQ %.3e, RAYQ %.3e, Lambda_ext %.3e, -p d(1/rho) %.3e; kP %.3e kE %.3e chi %.3e eq %s, cbar %.3e\n",
                   z, Tz, TF.var(:to)[z], TF.var(:egn)[z], TF.var(:ego)[z], ρz, TF.var(:do)[z], dt_, xch, gq, rq, mat.Λext[z] * dt_ / ρz,
                   -TF.var(:pgn)[z] * (1 / ρz - 1 / TF.var(:do)[z]), kp, ke, TF.var(:chifn)[z], string(mat.eqzone[z]), isempty(TX.CHEMB[]) ? 0.0 : TX.CHEMB[][z] * dt_)
                if ENVEUL !== nothing && !isempty(XMEPREV[])     # Eulerian budget: Δ(M e) per g vs exchange + advection (F e_b through both faces) + p dV; rest = MLT, diffusion, viscosity, …
                    xn = TF.var(:xmen); xo = XMEPREV[]; Mn = xn[z] - xn[z+1]; Mo = xo[z] - xo[z+1]; eb = TF.var(:egb)
                    Fi = xn[z] - xo[z]; Fo = xn[z+1] - xo[z+1]
                    adv = (Fi * eb[z] - Fo * eb[z+1]) / Mn; dME = (Mn * TF.var(:egn)[z] - Mo * TF.var(:ego)[z]) / Mn
                    rr_ = TF.var(:rn); uu_ = TF.var(:un); pdv = -TF.var(:pgn)[z] * (rr_[z+1]^2 * uu_[z+1] - rr_[z]^2 * uu_[z]) * dt_ / Mn   # (fixed mesh: p ∇·u, new u)
                    pr("        ZBUD %d: d(Me)/M %.3e = exch %.3e + adv %.3e (e_b in %.4e out %.4e, e %.4e, F in/M %.2e out/M %.2e) + GASQ %.3e + pdV %.3e + rest %.3e | x_H2 %.3e x_p %.3e\n",
                       z, dME, xch, adv, eb[z], eb[z+1], TF.var(:egn)[z], Fi / Mn, Fo / Mn, gq, pdv, dME - xch - adv - gq - pdv, ChemNet.full(mat.x[z], mat.net)[:H2], ChemNet.full(mat.x[z], mat.net)[:Hp])
                end
            end
        end
        tc0 = time(); Tb4 = CHEMDIAG ? copy(TF.var(:tn)) : nothing; ub4 = CHEMDIAG ? copy(TF.var(:egn)) : nothing; ENVEUL !== nothing && advect_species!()
        EB = isempty(TX.CHEMB[]) ? nothing : copy(TF.var(:egn))
        SLZ_ = haskey(ENV, "SLZ") ? parse.(Int, split(ENV["SLZ"], ",")) : nothing
        if SLZ_ !== nothing; xb_ = [copy(mat.x[z]) for z in SLZ_[1]:SLZ_[2]]; Tb_ = TF.var(:tn)[SLZ_[1]:SLZ_[2]]; end
        dT = PM.chem!(mat)
        if STOLCHEM
            TF.setscalar!(:tfac, TC.TSTP_IN[][1]); TF.setscalar!(:dtime, TC.TSTP_IN[][2]); TD.r(:timstp)
        end
        GRIDMAPP === nothing || gridmap_step!()
        if SLZ_ !== nothing && PM.SLDIAG[] !== nothing      # S6: starlight and chemistry of zones z1..z2 over this step's chem! (molecular-edge onset)
            sd = PM.SLDIAG[]; dtq = TF.getscalar(:timen) - TF.getscalar(:timeo)
            pr("        SLZ k %d dt %.3e: photosphere zone %d T_eff %.0f\n", k, dtq, sd.kph, sd.Te)
            ld_ = PM.LYADIAG[]
            ld_ === nothing || pr("        SLZ Lya: zones %d..%d (%d layers; photosphere rebuilt in zones %s, top T %.0f), %d bins, %d ALI iterations (dS %.1e), %.2f s\n", ld_.kin, ld_.kin + ld_.n - 1, ld_.nlayer, string(ld_.ksub), ld_.Tsurf, ld_.nb, ld_.it, ld_.δ, ld_.time)
            for (i, z) in enumerate(SLZ_[1]:SLZ_[2])
                d0_ = ChemNet.full(xb_[i], mat.net); d1_ = ChemNet.full(mat.x[z], mat.net)
                pr("        SLZ %d: n %.3e T %.1f -> %.1f  xH2 %.3e -> %.3e  xp %.3e -> %.3e  xHm %.2e | kLW %.3e fsh %.3e kHm %.3e eq %s E_r %.3e\n", z, TF.var(:dn)[z] / PC.μH,
                   Tb_[i], TF.var(:tn)[z], d0_[:H2], d1_[:H2], d0_[:Hp], d1_[:Hp], d1_[:Hm], sd.kLW[z], sd.fsh[z], sd.kHm[z], string(mat.eqzone[z]), TF.var(:ern)[z])
                t3_ = sd.tl3[z]
                (t3_ !== nothing && hasproperty(t3_, :qbb)) && pr("        SLZ %d tl3: T_c %.0f P12 %.2e RL %.2e beta %.2e ne %.2e | chemq %.3e cycq %.3e Lya-bb %.3e erg/cm3/s, Lya-bb dt/(rho e) %.2e\n", z,
                    t3_.Tc, t3_.P12, t3_.RL, t3_.β, t3_.ne, t3_.chemq, t3_.cycq, t3_.qbb, t3_.qbb * dtq / (TF.var(:dn)[z] * TF.var(:egn)[z]))
                if ld_ !== nothing && t3_ !== nothing && ld_.kin <= z < ld_.kin + ld_.n && !isnan(ld_.nJ[z-ld_.kin+1])
                    i_ = z - ld_.kin + 1; Tx_(nb) = nb > 0 ? 10.19881 * 1.602176634e-12 / (1.380649e-16 * log1p(1 / nb)) : 0.0
                    pr("        SLZ %d Lya: nbar(Jbar) %.3e (T_exc %.0f) nbar(S_L) %.3e (T_exc %.0f) nbar(J_ext) %.3e | R_Lya %.3e (Sobolev %.3e) P12 %.3e (Sobolev %.3e)\n", z,
                       ld_.nJ[i_], Tx_(ld_.nJ[i_]), ld_.nS[i_], Tx_(ld_.nS[i_]), ld_.Jext[i_], 3 * 6.2648e8 * ld_.oml[i_], ld_.Rsob[i_], t3_.P12, ld_.Psob[i_])
                end
            end
        end
        if EB !== nothing                            # BALSPLIT: c̄ ← c̄ + Δe_chem/Δt where the zone stayed on the network
            eg_ = TF.var(:egn); dtc_ = TF.getscalar(:timen) - TF.getscalar(:timeo); cb_ = TX.CHEMB[]; on_ = TX.CHEMBON[]
            for kk in eachindex(cb_)
                cb_[kk] = (on_[kk] && dtc_ > 0) ? cb_[kk] + (eg_[kk] - EB[kk]) / dtc_ : 0.0
            end
        end
        if CHEMDIAG && dT > 0.01                      # the zones chem! moved most
            tt_ = TF.var(:tn); jj = sort(collect(a:b); by = kk -> -abs(tt_[kk] / Tb4[kk] - 1))[1:3]
            for kk in jj
                pr("        chem! zone %d: T %.1f -> %.1f, m %.4e Msun, n %.3e, p %.3e, eq %s, e_before %.5e e_eos_after %.5e, dt %.2e\n", kk - a + 1, Tb4[kk], tt_[kk], 4π * xmi(kk+1) / Msun,
                   TF.var(:dn)[kk] / (1.22PC.mH), TF.var(:pgn)[kk], string(mat.eqzone[kk]), ub4[kk], TF.var(:egn)[kk], TF.getscalar(:timen) - TF.getscalar(:timeo))
            end
        end
        HOTON && (LNUC[] = hotburn!(TF.getscalar(:timen) - TF.getscalar(:timeo)))
        TCHEM[] += time() - tc0
        RAD && TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())   # T changed in chem!: Planck function, opacity, e_g + E_r/ρ
        REC.on && SNAPEVERY > 0 && k % SNAPEVERY == 0 && recsnap!(k)
        if SNAPS !== nothing                         # profiles for post-processing: every 0.25 dex of n_c, and every SNAPEVERY steps
            lnc = log10(TF.var(:dn)[a] / PC.mH)
            if lnc >= SNAPNEXT[] || (SNAPEVERY > 0 && k % SNAPEVERY == 0)
                push!(SNAPS, (; step = k, t = TF.getscalar(:timen), r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1],
                    ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b], y = [ycomp(a + j - 1) for j in 1:nsh],
                    m = 4π .* (ENVEUL === nothing ? TX.XMINT[][a:b+1] : envmenc()), Λext = copy(mat.Λext[a:b]), eq = [eqflag(a + j - 1) for j in 1:nsh],
                    er = TF.var(:ern)[a:b], fr = TF.var(:frn)[a:b+1], pg = TF.var(:pgn)[a:b], chi = TF.var(:chifn)[a:b],
                    fedd = TF.var(:fedd)[a:b], qf = TF.var(:qf)[a:b], xp = [ycomp(a + j - 1)[1] for j in 1:nsh], xh2 = [ycomp(a + j - 1)[3] for j in 1:nsh], L_lines = LLINE[]))
                lnc >= SNAPNEXT[] && (SNAPNEXT[] = lnc + 0.25)
                serialize(SNAPFILE, SNAPS)
            end
        end
        if cool == "lines"
            # H₂ line transfer on the new state (RHFC02 A2.4); Λ held fixed over the next step; Δt ≤ 0.1 u ρ/|Λ| (A2.2)
            rr = TF.var(:rn)[a:b+1]; uu = TF.var(:un)[a:b+1]; dd = TF.var(:dn)[a:b]; tt = TF.var(:tn)[a:b]
            nh2 = [dd[j] / PC.μH * ycomp(a + j - 1)[3] for j in 1:nsh]
            # LINEFRAC > 0: keep the last Λ until LINEFRAC × (shortest t_cool) has elapsed; transfer is most of a step's cost
            if LINEFRAC > 0 && TF.getscalar(:timen) - LINELAST[1] < LINEFRAC * LINELAST[2]
                TF.setscalar!(:dtime, min(TF.getscalar(:dtime), (IMPCOOL ? DTCOOL : 0.1) * LINELAST[2]))
            else
            tl0 = time(); Λl, Lr, _ = linecool(rr, uu, dd, tt, nh2, k); TLINE[] += time() - tl0
            tcmin = Inf
            for j in 1:nsh
                mat.Λext[a+j-1] = Λl[j]; mat.Tref[a+j-1] = tt[j]; mat.slope[a+j-1] = 0.0
                ue = TCOOLTH ? 1.5 * TF.var(:pgn)[a+j-1] / dd[j] : TF.var(:egn)[a+j-1] + PC.UOFF_TH     # energy scale of the cooling time
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
                j_ = JCMIN[]; pr("        line cap: t_cool %.3e s at zone %d (m %.4e Msun, T %.1f, n %.3e, e+U %.3e, Lambda %.3e); dtime before %.3e, smax %.3e\n", tcmin, j_, 4π * xmi(a+j_) / Msun, tt[j_], dd[j_] / (1.22PC.mH),
                                TF.var(:egn)[a+j_-1] + PC.UOFF_TH, Λl[j_], TF.getscalar(:dtime), TF.getscalar(:smax))
            end
            TF.setscalar!(:dtime, min(TF.getscalar(:dtime), (IMPCOOL ? DTCOOL : 0.1) * tcmin))
            LINELAST[1] = TF.getscalar(:timen); LINELAST[2] = tcmin
            end
        end
        if get(ENV, "STEPKX", "0") == "1"             # which equation and zone set the step's relative change (timstp!)
            pr("STEP %d dt %.3e smax %.3e: %s\n", k, TF.getscalar(:dtime), TF.getscalar(:smax),
               join(("$(TF.var(:kx)[i]):$(round(TF.var(:dx)[i], sigdigits = 2))" for i in 1:Int(TF.getscalar(:neqn))), " "))
        end
        if RAD && haskey(ENV, "ERDIAG") && k in parse.(Int, split(ENV["ERDIAG"], ","))     # radiation field across the grid
            ern = TF.var(:ern); frn = TF.var(:frn); fe = TF.var(:fedd); ch = TF.var(:chifn); dd = TF.var(:dn); tt = TF.var(:tn); rr = TF.var(:rn)
            Tb = parse(Float64, get(ENV, "TCMB", "0")); τ = 0.0; τs = zeros(b)
            for kk in b:-1:a; τ += ch[kk] * dd[kk] * (rr[kk+1] - rr[kk]); τs[kk] = τ; end
            pr("  radiation at step %d (geddr %.3f):\n", k, TF.getscalar(:geddr))
            for kk in (haskey(ENV, "ERDIAGZ") ? (zr = parse.(Int, split(ENV["ERDIAGZ"], ",")); (a - 1 + zr[1]):zr[2]:(a - 1 + zr[3])) : (a:6:b))     # ERDIAGZ=lo,step,hi (zones)
                pr("    z %3d m %.3e r %.3e tau %.3e dtau %.2e E/aT4 %.3e E/aTcmb4 %.3e F/cE %+.3e f %.3f T %.1f\n", kk - a + 1, 4π * xmi(kk+1) / Msun, rr[kk+1],
                   τs[kk], ch[kk] * dd[kk] * (rr[kk+1] - rr[kk]), ern[kk] / (Z0.arad * tt[kk]^4), ern[kk] / (Z0.arad * Tb^4), frn[kk+1] / (Z0.c * ern[kk]), fe[kk], tt[kk])
            end
        end
        if haskey(ENV, "PROFAT") && k in parse.(Int, split(ENV["PROFAT"], ","))
            rr = TF.var(:rn); uu = TF.var(:un); dd = TF.var(:dn); tt = TF.var(:tn)
            pr("  profile at step %d, t = %.4e yr:\n", k, TF.getscalar(:timen) / 3.156e7)
            for j in 1:8:nsh
                kk = a + j - 1
                pr("    j %3d m %.3e r %.3e pc  n %.3e  T %7.1f  u %7.3f km/s  f_H2 %.2e\n", j, 4π * xmi(kk+1) / Msun, rr[kk+1] / 3.0857e18,
                   dd[kk] / (1.22 * PC.mH), tt[kk], uu[kk+1] / 1e5, 2ycomp(kk)[3])
            end
        end
        ρc_ = TF.var(:dn)[a]; Tc_ = TF.var(:tn)[a]; ln = log10(ρc_ / PC.mH)
        if ln >= nextln[] || k % PRINTEVERY == 0
            umin = minimum(TF.var(:un)[a:b+1]) / 1e5
            pr("%7d %8.3f %10.4e %9.1f %9.1f %10.3e %10.3f %9.2e %9.2e%s\n", k, ln, TF.getscalar(:timen) / 3.156e7, Tc_, on98T(ln),
               2ycomp(a)[3], umin, dT, TF.getscalar(:dtime), cool == "lines" ? @sprintf("  L_lines %.3e  bins %d, thick lines %d", LLINE[], H2Lines.LAST[].nbins, H2Lines.LAST[].nthick) : "")
            ln >= nextln[] && (nextln[] += 0.25)
            if HOTON
                g_ = SG[]; tn_ = TF.var(:tn); un_ = TF.var(:un); jc_ = findfirst(<(-2e5), un_[a:b+1]); ms_ = jc_ === nothing ? NaN : 4π * xmi(a+jc_-1) / Msun; rs_ = jc_ === nothing ? NaN : TF.var(:rn)[a+jc_-1] / 6.957e10
                REC.mstar[] = ms_; REC.rstar[] = rs_
                pr("        hot: t %.6e yr, m* %.5g Msun, R* %.4g Rsun, L_nuc %.3e erg/s, X_c %.5f, D_c/D_0 %.3e, CNO_c %.2e, hot zones %d, conv %d, split %.1e/%d/%.1e, corr %.1e\n", TF.getscalar(:timen) / 3.156e7, ms_, rs_, LNUC[],
                   g_.Y[SM.IH, a], g_.Y[NN.IDX[:d], a] / (DH * HX), SM.xyz(view(g_.Y, :, a))[3], count(k -> tn_[k] > PM.HOT[].T1, a:b), HOTM[].nconv, LASTSPLIT[]..., HOTR[].corr)
            end
            if get(ENV, "UDIAG", "0") == "1"            # where the fastest infall is
                uu = TF.var(:un)[a:b+1]; ku = argmin(uu); rr = TF.var(:rn)[a:b+1]
                pr("        u_min %.3f km/s at interface %d of %d, r = %.3e pc, m = %.3e Msun, T(zone) = %.1f, n = %.3e\n", uu[ku] / 1e5, ku, nsh + 1,
                   rr[ku] / 3.0857e18, 4π * xmi(a+ku-1) / Msun, TF.var(:tn)[a+min(ku, nsh)-1], TF.var(:dn)[a+min(ku, nsh)-1] / (1.22 * PC.mH))
            end
            if get(ENV, "DTDIAG", "0") == "1"           # which variable/zone limits Δt (timstp: kx, dx per equation)
                kx = TF.var(:kx); dx = TF.var(:dx); names = ("r", "m", "rho", "u", "T", "E_r", "F_r"); ne_ = RAD ? 6 : 5   # (F_r is not in smax)
                j = argmax(abs.(dx[[1, 2, 3, 5, 6][1:ne_-1]])); j = (1, 2, 3, 5, 6)[j]; xs = TF.getscalar(:stol) / TF.getscalar(:smax)
                REC.dtlim[] = string(names[j], "@", kx[j])
                pr("        dt limit: %s in zone %d (rel. change %.2e, stol/smax = %.2f)%s\n", names[j], kx[j], dx[j], xs,
                   (RAD ? @sprintf("; E_r: zone %d %.2e, rho: zone %d %.2e", kx[6], dx[6], kx[3], dx[3]) : "") *
                   (DTENV === nothing ? "" : @sprintf("; DTENV cap %.2e s at interface %d", DTENVLAST[]...)))
            end
        end
        if ENVEUL !== nothing
            envsink!()
            if k % PRINTEVERY == 0
                j_ = envsurface(); uu = TF.var(:un)[a:b+1]; ku = argmin(uu); rr = TF.var(:rn)[a:b+1]; mm = 4π .* envmenc() ./ Msun
                pr("        env: chem fallbacks %d, t %.6e yr, Mc+grid(<surf) %.4f Msun, sink %.2e Msun/yr, surface r %.1f R_sun, u_min %.1f km/s at r %.1f (v_esc %.1f), dt %.2e s\n", PM.NCHEMFB[], TF.getscalar(:timen) / 3.156e7, mm[j_],
                   4π * ENVMDOT[] / Msun * 3.156e7, rr[j_] / Rsun, uu[ku] / 1e5, rr[ku] / Rsun, sqrt(2Z0.G * mm[ku] * Msun / rr[ku]) / 1e5, TF.getscalar(:dtime))
                # radiation force on the infall: Γ = χF/c / (G m/r²) per zone outside 2 × the surface; Ṁ = −4π r² ρ u at the surface and at 1e5 R_sun
                let fr_ = TF.var(:frn)[a:b+1], ch_ = TF.var(:chifn)[a:b], dn_ = TF.var(:dn)[a:b], Γ = zeros(nsh)
                    for j in 1:nsh
                        rr[j] > 2rr[j_] || continue
                        rc = 0.5(rr[j] + rr[j+1]); Γ[j] = ch_[j] * 0.5(fr_[j] + fr_[j+1]) / ChemNet.c / (Z0.G * 0.5(mm[j] + mm[j+1]) * Msun / rc^2)
                    end
                    jg = argmax(Γ); w5 = findfirst(>=(1e5 * Rsun), rr); w4 = findfirst(>=(1e4 * Rsun), rr)
                    mdot(j) = j === nothing || j > nsh ? NaN : -4π * rr[j]^2 * 0.5(dn_[max(j-1, 1)] + dn_[min(j, nsh)]) * uu[j] / Msun * 3.156e7
                    pr("        feedback: Gamma_max %.3f at r %.0f R_sun; Gamma(1e4) %.3f, Gamma(1e5) %.3f; Mdot surface %.3e, 1e4 %.3e, 1e5 %.3e Msun/yr; xadv1 %d; to table %d, from table %d\n",
                       Γ[jg], 0.5(rr[jg] + rr[jg+1]) / Rsun, w4 === nothing || w4 > nsh ? NaN : Γ[w4], w5 === nothing || w5 > nsh ? NaN : Γ[w5],
                       mdot(j_), mdot(w4), mdot(w5), NXADV1[], PM.NSWIN[], PM.NSWOUT[])
                end
                TX.GRIDWIN[] === nothing || pr("%s", gridwin_line())
                GRIDMAPP === nothing || pr("%s", gridmap_line())
                # causality of the radiation field outside the surface: F/cE at the interfaces (FLUXTIE or not)
                let fz = TF.var(:frn), ez = TF.var(:ern), rz = TF.var(:rn), q = 0.0, kq = 0, k1 = 0
                    for k in a+1:b
                        rz[k] > 2rr[j_] || continue
                        v = fz[k] / (ChemNet.c * 0.5(ez[k-1] + ez[k])); v > q && (q = v; kq = k)
                    end
                    if TD.FLUXTIE[][2] > 0 && length(TD.TIEW) >= b
                        k1 = something(findfirst(k -> TD.TIEW[k] > 0, a+1:b), 0); k1 > 0 && (k1 += a)
                    end
                    pr("        causality: max F/cE %.3f at zone %d (r %.0f R_sun)%s\n", q, kq, kq > 0 ? rz[kq] / Rsun : NaN,
                       k1 > 0 ? @sprintf("; flux tie from zone %d (r %.0f R_sun, h %.3f)", k1, rz[k1] / Rsun, TD.HEDD[k1]) : "")
                    if haskey(ENV, "FLUXTIEDIAG") && length(TD.HEDD) >= b
                        zr = parse.(Int, split(ENV["FLUXTIEDIAG"], ","))     # first, step, last zone
                        pr("        h vs F/cE (zone centre): %s\n", join((@sprintf("%d:%.3f/%.3f", k, TD.HEDD[k], 0.5(fz[k] + fz[k+1]) / (ChemNet.c * ez[k])) for k in zr[1]:zr[2]:min(b, zr[3])), " "))
                    end
                end
                if PM.SLDIAG[] !== nothing
                    sd_ = PM.SLDIAG[]; rz_ = TF.var(:rn); xs_ = String[]
                    for lr in (3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
                        kk = findfirst(k -> rz_[k] / Rsun >= 10^lr, sd_.kph+1:length(sd_.kLW))
                        kk === nothing && continue; kk += sd_.kph
                        push!(xs_, @sprintf("r 1e%.1f: kHm %.1e kH %.1e kLW %.1e fsh %.1e", lr, sd_.kHm[kk], sd_.kH[kk], sd_.kLW[kk], sd_.fsh[kk]))
                    end
                    pr("        light: photosphere zone %d r %.0f R_sun, T_eff %.0f K; %s\n", sd_.kph, TF.var(:rn)[sd_.kph+1] / Rsun, sd_.Te, join(xs_, "; "))
                end
            end
        end
        (Tc_ > Tstop || ρc_ / PC.mH > nstop) && (pr("stopped: T_c %.4e K (TSTOP %.2e), n_c %.4e (NSTOP %.2e) at step %d\n", Tc_, Tstop, ρc_ / PC.mH, nstop, k); TF.setscalar!(:tmax, TF.getscalar(:timen)))
        TF.getscalar(:dtime) < DTSTOP && (pr("stopped: dt = %.2e s below DTSTOP at step %d\n", TF.getscalar(:dtime), k); TF.setscalar!(:tmax, TF.getscalar(:timen)))
        if MSTOP > 0
            jm_ = findfirst(<(-2e5), TF.var(:un)[a:b+1]); mm_ = jm_ === nothing ? 0.0 : 4π * xmi(a + jm_ - 1) / Msun
            mm_ > MSTOP && (pr("stopped: M_star %.4f Msun (MSTOP %.3f) at step %d\n", mm_, MSTOP, k); TF.setscalar!(:tmax, TF.getscalar(:timen)))
        end
        REZ.on && RAD && k % REZ.every == 0 && rezone!(k)
        if TIMING && k % PRINTEVERY == 0
            pr("        timing since last report: step (hydro+radiation) %.1f s, chemistry %.1f s, H2 lines %.1f s, other output %.1f s\n",
               TSTEP[], TCHEM[], TLINE[], time() - tnow - 0.0 + TOTH[] - TCHEM[] - TLINE[])
            REC.on && recperf!(k, tnow)
            TSTEP[] = TCHEM[] = TLINE[] = TOTH[] = 0.0
        else
            TOTH[] += time() - tnow
        end
        TMARK[] = time()
        haskey(ENV, "SAVE") && CKEVERY > 0 && k % CKEVERY == 0 && checkpoint(ENV["SAVE"] * "_auto" * string((k ÷ CKEVERY) % 2), k)
    end)
catch e
    pr("stopped: %s\n", first(sprint(showerror, e), 400))
    get(ENV, "DEBUG", "0") == "1" && Base.display_error(stderr, e, catch_backtrace())
end
pr("run ended: %s; jstep %d jstepe %g timen %.6e tmax %.6e\n", string(RUNEND[]), Int(TF.getscalar(:jstep)), Float64(TF.getscalar(:jstepe)), TF.getscalar(:timen), TF.getscalar(:tmax))
haskey(ENV, "SAVE") && (checkpoint(ENV["SAVE"], KLAST[]); pr("checkpoint %s at step %d\n", ENV["SAVE"], KLAST[]))
IMPCOOL && pr("implicit line cooling on %d steps\n", NIMP[])
Trbdf2.DAMP[] && pr("damped Newton: %s\n", string(Trbdf2.DSTATS))
NLTE == "hyb" && pr("hybrid NLTE: %d ALI calls\n", NALICALL[])
TD.STEPFN[] === Trbdf2.step! && pr("integrator TR-BDF2, statistics %s\n", string(Trbdf2.STATS))
close(io)
