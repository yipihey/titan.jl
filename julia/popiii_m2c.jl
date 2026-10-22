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
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl")); include(joinpath(@__DIR__, "AtomicHHe.jl"))
include(joinpath(@__DIR__, "PopIIIMaterial.jl")); include(joinpath(@__DIR__, "H2Lines.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const PC = PrimordialChem; const PM = PopIIIMaterial; const Z0 = Z0Physics
# E2E-A: performance recording (TitanDriver.Perf, Perf.jl): always accumulated (overhead < 1 %), written to lib/popiii/perf/<PERF>_* when PERF=<name> is set
const Perf = TD.Perf; const PERFNAME = get(ENV, "PERF", "")
Perf.reset!(isempty(PERFNAME) ? "m2c" : PERFNAME; stage = get(ENV, "PERFSTAGE", "collapse"), driver = "popiii_m2c.jl", args = join(ARGS, " "), ic = get(ENV, "IC", "p100"),
            chem = get(ENV, "CHEM", "rhfc02"), nlte = get(ENV, "NLTE", "1"), sink = get(ENV, "SINK", "0"), integrator = get(ENV, "INTEGRATOR", "trbdf2"), tag = get(ENV, "TAG", ""))
using ForwardDiff
# INTEGRATOR=trbdf2 (default) | theta: second-order L-stable TR-BDF2 (Trbdf2.jl, G2) or TITAN's theta scheme
include(joinpath(@__DIR__, "Trbdf2.jl"))
# default since G2b: θ = 0.55 lets the z = 19 core oscillate at 1e11–1e13 (three-body heating phase) where TR-BDF2 is smooth
get(ENV, "INTEGRATOR", "trbdf2") == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
# S7b (opt-in): coupled protostellar feedback on the envelope (RadTransfer fields, radiation force, exchange with the protostar process); see s7b_env.jl
const S7B = get(ENV, "S7B", "0") == "1"
if S7B
    include(joinpath(@__DIR__, "RadTransfer.jl")); include(joinpath(@__DIR__, "s7b_exchange.jl")); include(joinpath(@__DIR__, "s7b_feedback.jl")); include(joinpath(@__DIR__, "s7b_env.jl"))
end
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
# S6b: RESTART=<snapshot file>, RESIDX=<index, default last>: start from a stored profile (same shell number as the IC above; the IC
# only supplies the dark-matter halo). TMAXYR = run time [yr] after the restart (default 1e3 t_ff,c as before). Default off.
const RESTART = get(ENV, "RESTART", "")
const T_OFF = Ref(0.0)                                   # time of the snapshot [s]
if RESTART != ""
    using Serialization
    sn = deserialize(RESTART); ix = parse(Int, get(ENV, "RESIDX", "0")); s0 = sn[ix <= 0 ? length(sn) + ix : ix]
    length(s0.ρ) <= nsh || error("restart: snapshot has $(length(s0.ρ)) shells, run has $nsh")
    global rI = copy(s0.r); global uI = copy(s0.u); global rho = copy(s0.ρ); global Tz = copy(s0.T); global yz = copy(s0.y); T_OFF[] = s0.t
    global s0m = copy(s0.m)                              # interface masses [g] (include the point mass of a sink snapshot)
    # SPLIT = s > 1: every zone is divided into s zones of equal mass (r³ linear in m inside the zone: the density of the zone is kept, u linear in m; T and composition copied)
    # — finer zones near the sink boundary (the zone masses of the collapse grid are 20 % of the enclosed mass, i.e. a zone spans a factor 3–4 in radius there)
    SPLIT = parse(Int, get(ENV, "SPLIT", "1"))
    if SPLIT > 1
        n0 = length(rho); mI_ = Float64[]; rI_ = Float64[]; uI_ = Float64[]; rho_ = Float64[]; Tz_ = Float64[]; yz_ = eltype(yz)[]
        for j in 1:n0, i in 0:SPLIT-1
            f = i / SPLIT
            push!(mI_, s0m[j] + f * (s0m[j+1] - s0m[j])); push!(rI_, cbrt(rI[j]^3 + f * (rI[j+1]^3 - rI[j]^3))); push!(uI_, uI[j] + f * (uI[j+1] - uI[j]))
            push!(rho_, rho[j]); push!(Tz_, Tz[j]); push!(yz_, yz[j])
        end
        push!(mI_, s0m[end]); push!(rI_, rI[end]); push!(uI_, uI[end])
        global s0m = mI_; global rI = rI_; global uI = uI_; global rho = rho_; global Tz = Tz_; global yz = yz_
        println("SPLIT = $SPLIT: $n0 -> $(length(rho)) zones")
    end
    global nsh = length(rho)
    ρc = rho[1]; Tc0 = Tz[1]
    println("restart from ", RESTART, " step ", s0.step, " t = ", s0.t / 3.156e7, " yr, n_c = ", rho[1] / PC.mH)
end
# S6b mass sink (SINK=1, needs RESTART or not): the grid is cut at RSINK [AU, default 1]: zones whose outer interface lies inside it are
# removed and their mass added to a central point mass (gravity only); the innermost interface is Lagrangian with a prescribed
# velocity (free-fall scaling u(a+1) sqrt(r(a+1)/r(a)): TITAN's llibc = 1), i.e. the grid is open to supersonic inflow.
const SINK = get(ENV, "SINK", "0") == "1"
const AU = 1.495978707e13
# sink radius r_s(M) = RSINKAU (M_pt / RSMREF)^RSINKP AU (M_pt ≥ RSMREF); RSINKP = 0: fixed. With RSINKP = 1 the free-fall time from r_s,
# t_ff(r_s, M) ∝ M, stays a fixed fraction of the accretion time M/Ṁ (Ṁ ≈ const), and the time step with it.
const RSINK0 = parse(Float64, get(ENV, "RSINKAU", "1")) * AU; const RSINKP = parse(Float64, get(ENV, "RSINKP", "0")); const RSMREF = parse(Float64, get(ENV, "RSMREF", "1")) * Msun
const RSMFLOOR = parse(Float64, get(ENV, "RSMFLOOR", string(RSMREF / Msun))) * Msun      # r_s stays at its value for M_pt = RSMFLOOR below that
rsink(M) = RSINK0 * (max(M, RSMFLOOR) / RSMREF)^RSINKP
const MPT = Ref(0.0)                                    # central point mass [g]
const MSKIP = Ref(0)                                    # zones removed from the restart profile
if SINK
    global j0 = 1
    if RESTART != ""
        while true                                       # r_s depends on the point mass, which grows with j0
            jn = findfirst(>=(rsink(s0m[j0])), rI); jn === nothing && (jn = length(rI) - 8)
            jn <= j0 && break
            global j0 = jn
        end
    end
    (nsh - j0 + 1) % 2 == 1 && (global j0 += 1)                # build! needs an even zone number
    MPT[] = RESTART != "" ? s0m[j0] : 0.0
    MSKIP[] = j0 - 1
    global rI = rI[j0:end]; global uI = uI[j0:end]; global rho = rho[j0:end]; global Tz = Tz[j0:end]; global yz = yz[j0:end]
    global nsh = length(rho)
    println("sink at ", rsink(MPT[]) / AU, " AU: ", j0 - 1, " zones removed, point mass ", MPT[] / Msun, " Msun, ", nsh, " zones left; r_inner = ", rI[1] / 1.495978707e13, " AU")
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
const atomic_dt = (get(ENV, "CHEM", "") |> c -> startswith(c, "netatomic")) ? parse(Float64, get(ENV, "ATOMDT", "0")) : 0.0
const ATOMTC = Ref(Inf)
const NIMP = Ref(0)                                                  # steps with implicit line cooling
const CHEM = get(ENV, "CHEM", "rhfc02")     # rhfc02: PrimordialGas (RHFC02 network, ρ_eq switch); net / netmodern / netmodernli: NetGas (G5)
if CHEM == "rhfc02"
    global mat = PM.PrimordialGas(ystore, PC.rates_rhfc02(; k5 = k5f), cf, parse(Float64, get(ENV, "RHOEQ", "5e-11")), get(ENV, "EQEOS", "scvh") == "scvh")
else
    # G5: ChemNet network; equilibrium table where the ideal-gas pressure ≥ PSWITCH (EOS validity, RHFC02's 1e4 dyn cm⁻²)
    # S6c: CHEM=netatomic[li] adds He, He⁺, He²⁺ chemistry and the atomic H/He cooling (AtomicHHe) to the net of netmodern[li]
    atomic = startswith(CHEM, "netatomic")
    netw, phot = CHEM == "net" ? (ChemRates.rhfc02(; k5 = k5f), ChemNet.Photo[]) :
                 atomic ? AtomicHHe.atomic(; li = CHEM == "netatomicli") : ChemRates.modern(; li = CHEM == "netmodernli")
    S7B && (atomic || error("S7B needs CHEM=netatomic[li]"); (netw, phot) = S7bEnv.prepare_net(netw, phot))
    tcmb = parse(Float64, get(ENV, "TCMB", "0"))
    Bν_(ν, T) = 2ChemNet.hP * ν^3 / ChemNet.c^2 / expm1(ChemNet.hP * ν / (ChemNet.kB * T))
    envn = (Trad = tcmb, J = tcmb > 0 ? (ν -> Bν_(ν, tcmb)) : (ν -> 0.0), static = true)   # same J in every zone: tabulated photo-rates
    S7B && (envn = S7bEnv.make_env(phot, tcmb))        # S7b: per-zone radiation environment (RadTransfer.ZoneEnv), the CMB until the coupling starts
    ytup(d) = (d[:Hp], d[:Hm], d[:H2], d[:H2p])
    eqd(ρ, T) = (y = PopIIIEqEOS.composition(ρ, T); Dict(:H => 1 - y[1] - 2y[3], :Hp => y[1], :H2 => y[3], :He => ChemNet.yHe, :Hep => 0.0, :Hepp => 0.0, :e => y[1], :Hm => 0.0, :H2p => 0.0))
    coolnet = (ρ, T, d) -> cf(ρ, T, d === nothing ? PopIIIEqEOS.composition(ρ, T) : ytup(d))
    if atomic
        coolat = AtomicHHe.cooling(phot, envn; eqcomp = eqd)
        coolnet = let c0 = coolnet; (ρ, T, d) -> c0(ρ, T, d) + coolat(ρ, T, d); end
        atomic_off = get(ENV, "ATOMICOFF", "0") == "1"            # sensitivity: atomic cooling off (chemistry only)
        atomic_off && (coolnet = let c0 = (ρ, T, d) -> cf(ρ, T, d === nothing ? PopIIIEqEOS.composition(ρ, T) : ytup(d)); c0; end)
    end
    xinit(y) = [s === :Hp ? y[1] : s === :Hm ? y[2] : s === :H2 ? y[3] : s === :H2p ? y[4] : s === :Lip ? 0.0 : 0.0 for s in netw.free]
    global mat = PM.NetGas([xinit(y) for y in ystore], netw, phot, coolnet; env = envn, pswitch = parse(Float64, get(ENV, "PSWITCH", "1e4")), Teqmax = atomic ? parse(Float64, get(ENV, "TEQMAX", "1.2e4")) : Inf)
    atomic && (AtomicHHe.THICKCONT[] = parse(Int, get(ENV, "THICKCONT", "0")))      # E2E-A: RA04-type thick-gas factor on the photo-pair emission (see AtomicHHe.THICKCONT)
    if atomic && get(ENV, "ESCAPE", "1") == "1"
        # escape probabilities of lines and recombination continua from the zone's size and velocity gradient (Sobolev, capped by the column to r)
        mat.zonehook = function (k, d)
            rn = TF.var(:rn); un = TF.var(:un); dn = TF.var(:dn)
            r1, r2 = rn[k], rn[k+1]; rm = 0.5 * (r1 + r2)
            grad = max(abs((un[k+1] - un[k]) / (r2 - r1)), abs(0.5 * (un[k] + un[k+1]) / rm), 1e-30)
            dd_ = d === nothing ? eqd(dn[k], TF.var(:tn)[k]) : d
            AtomicHHe.set_zone!(rm, 1 / grad, dn[k] / ChemNet.μH(), dd_)
            if AtomicHHe.THICKCONT[] == 2           # E2E-A: Planck-mean continuum optical depth of the zone over r (only where it can matter: n_H > 1e12)
                nHk = dn[k] / ChemNet.μH(); Tk = TF.var(:tn)[k]
                AtomicHHe.ZS.thick = nHk < 1e12 ? 1.0 : (κP = ContinuumOpacity.means(dn[k], Tk; nH = nHk * get(dd_, :H, 0.0), nH2 = nHk * get(dd_, :H2, 0.0), nHe = nHk * get(dd_, :He, ChemNet.yHe), nHp = nHk * get(dd_, :Hp, 0.0),
                                                   ne = nHk * get(dd_, :e, 0.0), nHm = nHk * get(dd_, :Hm, 0.0), nHep = nHk * get(dd_, :Hep, 0.0), nHepp = nHk * get(dd_, :Hepp, 0.0), νs = ContinuumOpacity.grid(; n = 250, T = Tk)).κP;
                                                   AtomicHHe.βesc(κP * dn[k] * rm))
            end
        end
    end
end
S7B && S7bEnv.install!(mat)
"Composition of zone k as the PrimordialChem tuple (x_H⁺, x_H⁻, x_H₂, x_H₂⁺), for either material."
ycomp(k) = mat isa PM.NetGas ? (d = ChemNet.full(mat.x[k], mat.net); (d[:Hp], d[:Hm], d[:H2], d[:H2p])) : mat.y[k]
"Equilibrium-regime flag of zone k."
eqflag(k) = mat isa PM.NetGas ? mat.eqzone[k] : TF.var(:dn)[k] >= mat.ρeq
# artificial viscosity: viscous length l = ql1 r ≈ 2 Δr (the SOD state's ql0 = 1e-4 cm is zero at these scales, M2g);
# Δr/r ≈ 0.04 (core) … 0.12 (envelope) on the geometric mass grid; quadratic term only (cq2), as PORTING_NOTES P7
av = (ql0 = 0.0, ql1 = parse(Float64, get(ENV, "AVQL1", "0.15")), cq1 = 0.0, cq2 = parse(Float64, get(ENV, "AVCQ2", "1.0")))
SINK && (TD.INNER_GHOST_GEOM[] = true)
PS.build!(; rI, rho, T = Tz, u = uI, geom = 2, lgrid = 3, inner = SINK ? :lagr : :wall, outer = :free, pext = 0.0,
          selfgrav = true, material = mat, av, Mc = MPT[])
TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))    # G1: viscous length = AVZONE × zone width (0 = ql1 r̄ above)
TX.REGULAR_CENTRE[] = !SINK; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
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
PM.CHEM_FALLBACK[] = get(ENV, "CHEMFALLBACK", "0") == "1"; ChemNet.XFLOOR[] = parse(Float64, get(ENV, "XFLOOR", "1e-20"))
CHEM != "rhfc02" && startswith(CHEM, "netatomic") && (ChemNet.NESTED[] = get(ENV, "NESTED", "1") == "1")      # S6c: bisection-in-T fallback of the network step (stiff atomic cooling)
TF.setscalar!(:thet, parse(Float64, get(ENV, "THETA", "0.55")))   # 0.5 lets the centre ring (M2c)
TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.01"))
tff = sqrt(3π / (32Z0.G * ρc))
TF.setscalar!(:dtime, parse(Float64, get(ENV, "DT0", "1")) * 1e-3tff); TF.setscalar!(:tmax, haskey(ENV, "TMAXYR") ? T_OFF[] + parse(Float64, ENV["TMAXYR"]) * 3.156e7 : 1e3tff)
RESTART != "" && (TF.setscalar!(:timen, T_OFF[]); TF.setscalar!(:timeo, T_OFF[]))
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
const SNAPFILE = joinpath(get(ENV, "SNAPDIR", joinpath(TF.ROOT, "lib", "popiii")), "snaps_$(get(ENV, "IC", "p100"))_$(cool)_n$(nsh)$(get(ENV, "TAG", "")).jls")
const SNAPS = get(ENV, "SNAPS", "0") == "1" ? Any[] : nothing
const SNAPNEXT = Ref(log10(rho[1] / PC.mH)); const SNAPEVERY = parse(Int, get(ENV, "SNAPEVERY", "0"))
const HI21ON = get(ENV, "HI21", "0") == "1"; const HI21N = parse(Int, get(ENV, "HI21EVERY", "200"))   # S8b (opt-in): 21 cm diagnostic (HI21.jl) every HI21EVERY steps -> lib/popiii/hi21_<run>.csv
HI21ON && include(joinpath(@__DIR__, "HI21.jl"))
HI21ON && HI21.open_inline!(joinpath(TF.ROOT, "lib", "popiii"), get(ENV, "HI21RUN", "m2c_$(get(ENV, "IC", "p100"))_$(cool)_n$(nsh)$(get(ENV, "TAG", ""))"); Tcmb = parse(Float64, get(ENV, "TCMB", "54.5")))
# NLTE (line source functions): 1 = GP98/HM79 factor on S (old default), 0 = LTE, ep = statistical equilibrium with escape
# probabilities (H2Levels), ali = statistical equilibrium coupled to the transfer (ALI, G6c), hyb = ALI every NALI steps
# and EP × (ALI/EP per zone, from the last ALI) in between, the interval adaptive (see linecool)
const NLTE = get(ENV, "NLTE", "1"); const NALI = parse(Int, get(ENV, "NALI", "20")); const ALIFAC = Ref{Vector{Float64}}(Float64[]); const ALILAST = Ref(-10^9); const ALIINT = Ref(1); const NALICALL = Ref(0)
NLTE in ("ep", "ali", "hyb") && include(joinpath(@__DIR__, "H2Levels.jl"))
const HLV = NLTE in ("ep", "ali", "hyb") ? H2Levels.load() : nothing
"Line cooling per zone (and Lr) for zone temperatures tt; step k (for NLTE=hyb)."
function linecool(rr, uu, dd, tt, nh2, k)
    if NLTE in ("ep", "ali", "hyb")
        yc = [ycomp(a + j - 1) for j in 1:nsh]; nHn = dd ./ PC.μH
        nl = (LV = HLV, nH = [PC.xH(yc[j]) * nHn[j] for j in 1:nsh], nHe = PC.yHe .* nHn, nHp = [yc[j][1] * nHn[j] for j in 1:nsh])
        NLTE != "hyb" && return H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = NLTE == "ali" ? :ali : :ep)
        Λe, Lr, x = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ep)
        # ALI refresh interval adaptive (≤ NALI): halved when the ALI/EP factor moved by > 5 % since the last refresh,
        # doubled when < 1 % — a factor frozen while the lines turn thick jumps at the refresh and kicks the core
        if k >= 0 && (k - ALILAST[] >= ALIINT[] || length(ALIFAC[]) != nsh)
            Λa, Lra, xa = Perf.@section(:linecool_ali, H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ali))
            Perf.count!(:ali_calls)
            fnew = [Λe[j] != 0 && Λa[j] / Λe[j] > 0 ? Λa[j] / Λe[j] : 1.0 for j in 1:nsh]
            if length(ALIFAC[]) == nsh
                δ = maximum(abs.(log.(fnew ./ ALIFAC[])))
                ALIINT[] = δ > 0.05 ? max(1, ALIINT[] ÷ 2) : δ < 0.01 ? min(NALI, 2ALIINT[]) : ALIINT[]
            end
            ALIFAC[] = fnew; ALILAST[] = k; NALICALL[] += 1
            return Λa, Lra, xa
        end
        return Λe .* ALIFAC[], Lr, x
    end
    # NLTE=1: source function × f = 1/(1 + Λ_LTE/(n_H Λ_0)) (the GP98/HM79 combination of PrimordialChem.cool_H2_thin)
    fnl = [1 / (1 + PC.Λ_H2_LTE_HM(tt[j]) / max(PC.xH(ycomp(a + j - 1)) * dd[j] / PC.μH * PC.Λ_H2_0_GP(tt[j]), 1e-300)) for j in 1:nsh]
    H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, fnlte = NLTE == "1" ? fnl : nothing)
end
const HLL = cool == "lines" ? H2Lines.load() : nothing; const NBIN = parse(Int, get(ENV, "NBIN", "30")); const LLINE = Ref(0.0)
# S6b sink bookkeeping (SINK=1): zone removal, prescribed inner velocity, the mass flux through r = RSINK and the infall state there
const SINKIO = SINK ? open(joinpath(TF.ROOT, "lib", "popiii", "s6c_$(get(ENV, "TAG", "")).csv"), "w") : nothing
const INNERBC = get(ENV, "INNERBC", "ballistic")
const LINEEVERY = parse(Int, get(ENV, "LINEEVERY", "1")); const LASTLL = Ref{Any}(([], [])); const TCMIN = Ref(Inf)
const PERFLAST = Ref(time_ns()); const PERFFAIL = Ref(0)
const SINKREM = Ref(0); const TLINE = Ref(0.0); const TCHEM = Ref(0.0); const T0WALL = time()             # wall time in linecool [s]
const RPREV = Ref(Float64[]); const TPREV = Ref(0.0)
const EVIO = SINK ? open(joinpath(TF.ROOT, "lib", "popiii", "s6c_$(get(ENV, "TAG", ""))_events.csv"), "w") : nothing
"Interior mass [Msun], radius [AU] and Eulerian mass flux 4π r² ρ̄ |u| [Msun/yr] through interface i (ρ̄ = mean of the two adjacent zones): the steady-state Mdot at that mass coordinate (OHY10 eq. 13)."
function fluxcols(i, b_)
    i = min(i, b_ - 1); rn = TF.var(:rn); dn = TF.var(:dn); un = TF.var(:un)
    string(4π * TX.XMINT[][i] / Msun, ",", rn[i] / AU, ",", 4π * rn[i]^2 * 0.5 * (dn[i-1] + dn[i]) * abs(un[i]) / Msun * 3.15576e7)
end
function sink_update!(k)
    a_ = Int(TF.getscalar(:ngrs)); b_ = Int(TF.getscalar(:ngre))
    rn = TF.var(:rn); xm = TX.XMINT[]; dn = TF.var(:dn); un = TF.var(:un); dv = TF.var(:dvoln); t = TF.getscalar(:timen)
    k == 1 && println(EVIO, "step,t_cross_yr,m_interface,r_s_AU,u_kms,rho,T,fH2,xe,e,p,cs_kms,rho_up,T_up,u_up_kms")
    nrem = 0; rs = rsink(4π * xm[a_])
    while rn[a_+nrem] < rs && b_ - a_ - nrem > 8
        i = a_ + nrem                                    # interface i has crossed r_s: the mass inside it, m_i, is inside the sink
        rp = length(RPREV[]) >= i ? RPREV[][i] : rn[i]
        fr = rp > rn[i] ? clamp((rp - rs) / (rp - rn[i]), 0.0, 1.0) : 1.0
        tc = TPREV[] + fr * (t - TPREV[])
        j = min(i, b_); yc = ycomp(j)
        println(EVIO, k, ",", tc / 3.15576e7, ",", 4π * xm[i] / Msun, ",", rs / AU, ",", un[i] / 1e5, ",", dn[j], ",", TF.var(:tn)[j], ",", 2yc[3], ",", yc[1], ",",
                TF.var(:egn)[j], ",", TF.var(:pgn)[j], ",", TF.var(:asn)[j] / 1e5, ",", dn[min(j + 1, b_)], ",", TF.var(:tn)[min(j + 1, b_)], ",", un[min(i + 1, b_)] / 1e5)
        nrem += 1; rs = rsink(4π * xm[a_ + nrem])
    end
    flush(EVIO)
    if nrem > 0
        a_ += nrem; SINKREM[] += nrem
        TF.setscalar!(:ngrs, a_); global a = a_; global nsh = b_ - a_ + 1
        for _ in 1:nrem; length(ALIFAC[]) > 0 && popfirst!(ALIFAC[]); end
        TF.setscalar!(:delml, dn[a_] * dv[a_]); TF.setscalar!(:delrl, rn[a_+1] - rn[a_])
    end
    RPREV[] = copy(rn); TPREV[] = t
    # inner interface: ballistic (no pressure force: the grid is open to supersonic inflow), mean velocity over the next step u − ½ g dt. (Setting it from the
    # neighbour, u(a+1) sqrt(r(a+1)/r(a)), is a positive feedback through the pressure of the first zone and runs away.)
    TF.setscalar!(:uextl, INNERBC == "ballistic" ? un[a_] - 0.5 * Z0.G * 4π * xm[a_] / rn[a_]^2 * TF.getscalar(:dtime) : un[a_+1] * sqrt(rn[a_+1] / rn[a_]))
    rs = rsink(4π * xm[a_])
    Mpt = 4π * xm[a_]; ρs = dn[a_]; us = un[a_]
    Min = Mpt
    if MSTOP > 0 && Mpt >= MSTOP
        TF.setscalar!(:tmax, TF.getscalar(:timen)); pr("sink stop: M_pt = %.3f Msun >= MSTOPMSUN at step %d, t = %.5f yr\n", Mpt / Msun, k, (t - T_OFF[]) / 3.15576e7)
    end
    if MDOTSTOP > 0
        fl = parse(Float64, split(fluxcols(a_ + 3, b_), ",")[3]); MDOTLOW[] = fl < MDOTSTOP ? MDOTLOW[] + 1 : 0
        MDOTLOW[] >= MDOTSTOPN && (TF.setscalar!(:tmax, TF.getscalar(:timen)); pr("sink stop: mass flux %.3e Msun/yr < MDOTSTOP for %d steps (M_pt = %.3f Msun, t = %.5f yr)\n", fl, MDOTLOW[], Mpt / Msun, (t - T_OFF[]) / 3.15576e7))
    end
    if haskey(ENV, "SINKPROF") && k % parse(Int, ENV["SINKPROF"]) == 0
        pr("  sink profile step %d t = %.5f yr: first zones (k, r/AU, u km/s, rho, T, p, fH2, M_int/Msun)\n", k, (t - T_OFF[]) / 3.15576e7)
        for kk in a_:min(a_ + 11, b_)
            pr("    %3d %9.5f %8.2f %10.3e %9.0f %10.3e %6.3f %9.5f\n", kk, rn[kk] / AU, un[kk] / 1e5, dn[kk], TF.var(:tn)[kk], TF.var(:pgn)[kk], 2ycomp(kk)[3], 4π * xm[kk] / Msun)
        end
    end
    if k % parse(Int, get(ENV, "SINKLOG", "1")) == 0
        i = a_; j = min(a_ + 1, b_)
        k == 1 && println(SINKIO, "step,t_yr,dt_s,r_in_AU,r_s_AU,M_pt,M_in,rho,T,u_kms,cs_kms,fH2,xe,e,p,rho2,T2,u2_kms,fH2_2,ngrs,nrem,m3,r3_AU,Mdot3,m6,r6_AU,Mdot6")
        yc1 = ycomp(i); yc2 = ycomp(j)
        println(SINKIO, k, ",", t / 3.15576e7, ",", TF.getscalar(:dtime), ",", rn[a_] / AU, ",", rs / AU, ",", Mpt / Msun, ",", Min / Msun, ",",
                ρs, ",", TF.var(:tn)[i], ",", us / 1e5, ",", TF.var(:asn)[i] / 1e5, ",", 2yc1[3], ",", yc1[1], ",", TF.var(:egn)[i], ",", TF.var(:pgn)[i], ",",
                dn[j], ",", TF.var(:tn)[j], ",", un[j+1] / 1e5, ",", 2yc2[3], ",", a_, ",", SINKREM[], ",", fluxcols(a_ + 3, b_), ",", fluxcols(a_ + 6, b_))
        flush(SINKIO)
    end
end
"Profile snapshot of the current state (restart file of RESTART=, analysis)."
push_snap!(k) = push!(SNAPS, (; step = k, t = TF.getscalar(:timen), r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1],
    ρ = TF.var(:dn)[a:b], T = TF.var(:tn)[a:b], e = TF.var(:egn)[a:b], y = [ycomp(a + j - 1) for j in 1:nsh],
    m = 4π .* TX.XMINT[][a:b+1], Λext = copy(mat.Λext[a:b]), eq = [eqflag(a + j - 1) for j in 1:nsh]))
Tstop = parse(Float64, get(ENV, "TSTOP", "6000")); nstop = parse(Float64, get(ENV, "NSTOP", "1e18"))
# E2E-A stops of the sink run (default off): MSTOPMSUN = point mass [Msun] at which accretion is declared finished (a documented cap: without feedback the cloud
# keeps feeding the star); MDOTSTOP [Msun/yr] with MDOTSTOPN [200]: stop when the Eulerian mass flux at the third interface has stayed below it for that many steps
const MSTOP = parse(Float64, get(ENV, "MSTOPMSUN", "0")) * Msun
const MDOTSTOP = parse(Float64, get(ENV, "MDOTSTOP", "0")); const MDOTSTOPN = parse(Int, get(ENV, "MDOTSTOPN", "200")); const MDOTLOW = Ref(0)
PERFLAST[] = time_ns()
try
    TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = function (k)
        _tout = time_ns(); _tc = _tout; dT = PM.chem!(mat); TCHEM[] += (time_ns() - _tc) * 1e-9
        Perf.add!(:chem, time_ns() - _tc, 0); Perf.add!(:step, _tout - PERFLAST[], 0); Perf.count!(:steps)     # (:step = integrator time since the end of the last output call)
        Perf.count!(:chem_fail, PM.CHEM_FAILS[] - PERFFAIL[]); PERFFAIL[] = PM.CHEM_FAILS[]
        SINK && Perf.@section(:sink, sink_update!(k))
        S7B && S7bEnv.after_step(k)
        if atomic_dt > 0           # S6c: Δt ≤ ATOMDT × (thermal energy / atomic cooling rate) of the fastest-cooling zone (the cooling inside the chemistry step is implicit; this keeps the operator splitting accurate)
            tcmin_at = Inf
            for j in a:b
                (TF.var(:pgn)[j] >= mat.pswitch) && continue
                d_ = ChemNet.full(mat.x[j], mat.net); T_ = TF.var(:tn)[j]; ρ_ = TF.var(:dn)[j]
                Λ_ = abs(coolat(ρ_, T_, d_)); Λ_ > 0 || continue
                uth = 1.5 * ChemNet.kB * T_ * sum(values(d_)) / ChemNet.μH()
                tcmin_at = min(tcmin_at, uth * ρ_ / Λ_)
            end
            ATOMTC[] = tcmin_at; TF.setscalar!(:dtime, min(TF.getscalar(:dtime), atomic_dt * tcmin_at))
        end
        if SNAPS !== nothing                         # profiles for post-processing: every 0.25 dex of n_c, and every SNAPEVERY steps
            lnc = log10(TF.var(:dn)[a] / PC.mH)
            if lnc >= SNAPNEXT[] || (SNAPEVERY > 0 && k % SNAPEVERY == 0)
                push_snap!(k)
                lnc >= SNAPNEXT[] && (SNAPNEXT[] = lnc + 0.25)
                Perf.@section :io serialize(SNAPFILE, SNAPS)
            end
        end
        HI21ON && k % HI21N == 0 && HI21.inline!(k, TF.getscalar(:timen), TF.var(:rn)[a:b+1], TF.var(:un)[a:b+1], TF.var(:dn)[a:b], TF.var(:tn)[a:b], [ycomp(a + j - 1) for j in 1:nsh], 4π .* TX.XMINT[][a:b+1])
        if cool == "lines"
            # H₂ line transfer on the new state (RHFC02 A2.4); Λ held fixed over the next step; Δt ≤ 0.1 u ρ/|Λ| (A2.2)
            rr = TF.var(:rn)[a:b+1]; uu = TF.var(:un)[a:b+1]; dd = TF.var(:dn)[a:b]; tt = TF.var(:tn)[a:b]
            nh2 = [dd[j] / PC.μH * ycomp(a + j - 1)[3] for j in 1:nsh]
            _t0 = time_ns()
            if LINEEVERY == 1 || k % LINEEVERY == 1 || length(LASTLL[][1]) != nsh
                Λl, Lr, _ = Perf.@section(:linecool, linecool(rr, uu, dd, tt, nh2, k)); LASTLL[] = (Λl, Lr)
                TLINE[] += (time_ns() - _t0) * 1e-9
                tcmin = Inf
                for j in 1:nsh
                    mat.Λext[a+j-1] = Λl[j]; mat.Tref[a+j-1] = tt[j]; mat.slope[a+j-1] = 0.0
                    ue = TF.var(:egn)[a+j-1] + PC.UOFF_TH                 # thermal-ish energy scale (positive)
                    Λl[j] != 0 && (tcmin = min(tcmin, abs(ue * dd[j] / Λl[j])))
                end
                # G7 (IMPCOOL=1, default): where the explicit cap Δt ≤ 0.1 t_cool would bind, the slope d ln Λ/d ln T from a
                # second transfer call with all T × (1 + ε) makes Λ_j (T/T_j)^s implicit in the chemistry step, and the cap
                # relaxes to Δt ≤ DTCOOL t_cool (default 1); otherwise Λ is held fixed (explicit, accurate there)
                if IMPCOOL && (IMPALWAYS || 0.1tcmin < TF.getscalar(:dtime))
                    ε = 0.02; Λp, _, _ = Perf.@section(:linecool, linecool(rr, uu, dd, tt .* (1 + ε), nh2, -1))
                    for j in 1:nsh
                        Λl[j] != 0 && Λp[j] / Λl[j] > 0 && (mat.slope[a+j-1] = clamp(log(Λp[j] / Λl[j]) / log(1 + ε), -10.0, 30.0))
                    end
                    NIMP[] += 1
                end
                LLINE[] = Lr[end]; TCMIN[] = tcmin
            end                        # S6b: with LINEEVERY > 1 the line cooling (Λext, Tref, slope) is held from the last transfer call
            TF.setscalar!(:dtime, min(TF.getscalar(:dtime), (IMPCOOL ? DTCOOL : 0.1) * TCMIN[]))
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
        k % 50 == 0 && Perf.sample!(step = k, t_yr = (TF.getscalar(:timen) - T_OFF[]) / 3.156e7, dt_s = TF.getscalar(:dtime), zones = nsh, newton = Int(TF.getscalar(:iter)), log_nc = ln, Tc = Tc_,
                                    wall_s = Perf.elapsed())
        PERFLAST[] = time_ns(); Perf.add!(:output, PERFLAST[] - _tout, 0)     # (whole output callback, inclusive of :chem, :linecool, :sink, :io)
    end)
catch e
    pr("stopped: %s\n", sprint(showerror, e)[1:min(end, 400)])
    get(ENV, "DEBUG", "0") == "1" && Base.display_error(stderr, e, catch_backtrace())
end
if SNAPS !== nothing && get(ENV, "SNAPEND", "0") == "1"          # E2E-A: the final state as the last snapshot (restart of the next stage)
    push_snap!(-1); serialize(SNAPFILE, SNAPS)
end
S7B && S7bEnv.finish()
pr("wall time: total %.1f s, line transfer %.1f s, chemistry %.1f s (%d failed zone steps kept)\n", time() - T0WALL, TLINE[], TCHEM[], PM.CHEM_FAILS[])
IMPCOOL && pr("implicit line cooling on %d steps\n", NIMP[])
NLTE == "hyb" && pr("hybrid NLTE: %d ALI calls\n", NALICALL[])
TD.STEPFN[] === Trbdf2.step! && pr("integrator TR-BDF2, statistics %s\n", string(Trbdf2.STATS))
if !isempty(PERFNAME)
    Perf.write(); Perf.summary()
end
close(io)
