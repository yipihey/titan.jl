# Population III protostar formation and accretion from the ABN02 cloud: 1D radiation hydrodynamics (TITAN, pure Julia) with the
# coupled implicit primordial chemistry, H₂ line transfer, star light, the G2 network material and the GRIDEQ moving mesh — the one
# supported configuration is julia/popiii_stack.jl (docs/POPIII_STACK.md); launch with scripts/popiii_run/iceulrun.sh.
#   julia --project=<env> popiii_m2c.jl lines <IC shells>      (start from the IC; RESTART=<checkpoint> continues a run)
using Printf
include(joinpath(@__DIR__, "popiii_stack.jl"))     # the one supported configuration: defaults for every unset stack key
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "PrimordialChem.jl"))
include(joinpath(@__DIR__, "SCVH.jl")); include(joinpath(@__DIR__, "PopIIIEqEOS.jl"))
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "G2Opacity.jl")); include(joinpath(@__DIR__, "ChemNet.jl")); include(joinpath(@__DIR__, "ChemRates.jl"))
include(joinpath(@__DIR__, "PopIIIMaterial.jl")); include(joinpath(@__DIR__, "H2Lines.jl"))
# LINEGPU=cuda|metal|cpu: formal solution of the thick H2 lines through KernelAbstractions (H2LinesKA.jl); needs an
# environment with KernelAbstractions and CUDA or Metal. Off unless set.
haskey(ENV, "LINEGPU") && include(joinpath(@__DIR__, "h2lines_gpu.jl"))
include(joinpath(@__DIR__, "StellarOpacity.jl")); include(joinpath(@__DIR__, "IonRay.jl"))
# MGROUPS (docs/POPIII_MULTIGROUP.md; default off): the pivoted banded solver for the wider block system
haskey(ENV, "MGROUPS") && ENV["MGROUPS"] != "0" && include(joinpath(@__DIR__, "JuliaSolver.jl"))
# S6-C1: HOT=1: stellar EOS, nuclear network and convective mixing in the hot core (T > HOTT1), see PopIIIMaterial.HOT
const HOTON = get(ENV, "HOT", "0") == "1"
if HOTON
    for f in ("StellarEOS", "NucNet", "StellarMaterial", "StarMaterial"); include(joinpath(@__DIR__, f * ".jl")); end
end
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const PC = PrimordialChem; const PM = PopIIIMaterial; const Z0 = Z0Physics
using ForwardDiff
# INTEGRATOR=trbdf2 (default) | theta: second-order L-stable TR-BDF2 (Trbdf2.jl, G2) or TITAN's theta scheme
include(joinpath(@__DIR__, "Trbdf2.jl")); include(joinpath(@__DIR__, "NewtonKrylov.jl")); include(joinpath(@__DIR__, "StageCons.jl"))
# default since G2b: θ = 0.55 lets the z = 19 core oscillate at 1e11–1e13 (three-body heating phase) where TR-BDF2 is smooth
get(ENV, "INTEGRATOR", "trbdf2") == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
const Msun = 1.98847e33

cool = isempty(ARGS) ? "ra04" : ARGS[1]; nsh = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100
cf = cool == "thin" ? PC.cool_H2_thin : cool == "ra04lines" ? PC.cool_H2_ra04 :
     cool == "lines" ? PC.cool_CIE_ra04 :                     # H₂ lines by transfer (H2Lines), continuum: RA04 CIE
     (ρ, T, y) -> PC.cool_H2_ra04(ρ, T, y) + PC.cool_CIE_ra04(ρ, T, y)
haskey(ENV, "TCMB") && (cf = PC.cool_floor(cf, parse(Float64, ENV["TCMB"])))
# RAD=1 (the stack; from the IC): TITAN's radiation (E_r, F_r, Eddington factors) and the StellarOpacity continuum; the RA04
# CIE escape-factor term is then not used (the continuum is in the transfer), the H₂ lines stay as Λext
const RAD = get(ENV, "RAD", "0") == "1"
RAD && (cf = (ρ, T, y) -> zero(T))
# POPIII_OUT: directory for the run log and the snapshot file (default lib/popiii); checkpoints follow CKDIR
const OUTDIR = get(ENV, "POPIII_OUT", joinpath(TF.ROOT, "lib", "popiii")); mkpath(OUTDIR)
io = open(joinpath(OUTDIR, "m2c_$(get(ENV, "IC", "p100") == "p100" ? "" : get(ENV, "IC", "") * "_")$(cool)_n$(nsh)$(get(ENV, "TAG", "")).log"), haskey(ENV, "RESTART") ? "a" : "w")   # (a continued run appends)
pr(f, a...) = (s = Printf.format(Printf.Format(f), a...); print(s); print(io, s); flush(io))

ic = get(ENV, "IC", "abn02red")
dm1 = 3e-7Msun
q_of(M, n) = let lo = 1.0001, hi = 2.0
    f(q) = dm1 * (q^n - 1) / (q - 1) - M
    for _ in 1:200; m_ = (lo + hi) / 2; f(m_) > 0 ? (hi = m_) : (lo = m_); end
    (lo + hi) / 2
end
ic in ("abn02red", "abn02fin", "abn02z19") || error("IC=$ic: the supported initial conditions are abn02red (main), abn02fin, abn02z19")
begin
    # ABN02 (Abel, Bryan & Norman 2002, Science 295, 93; arXiv astro-ph/0112088) Fig. 2, profile "9 Myr later"
    # (red dotted with circles), read off the 400-dpi rendering [read-src, approximate ±0.05 dex]: particle number
    # density n(r) (panel A), T(M), H₂ mass fraction (panel C, lower curve), v_r(M) (panel E). ρ = n μ m_H, μ = 1.22.
    z19 = startswith(ic, "abn02z19")
    # IC=abn02fin: ABN02's FINAL output (z = 18.181164, "200 years later", black dotted with circles), taken exactly from the
    # vector paths of the arXiv PDF (lib/popiii/abn02fin: CSV, README, overlay; scripts/abn02_fig2). Enclosed mass M(r) from
    # panel B directly (no μ assumption; inside its first point, r = 4.8e-6 pc, uniform at the density that gives its mass);
    # T(M) panel D; v_r(M) panel E as the initial velocity (supersonic infall, NOT hydrostatic; inside the first tabulated
    # mass, 0.01 M☉, homologous v ∝ r to the regular centre); final-time DM line of panel A for TX.EXTMASS. Composition from
    # panel C: its two curves are the H₂ and atomic-H MASS fractions (f_H2 + f_HI = 0.756 ± 0.005 = X at every mass, so the
    # second curve is neutral H, as T. Abel (ABN02's first author) confirmed, Nov 2026), x_H₂ = f_H2 / (2 (f_H2 + f_HI)) per
    # H nucleus. x_e is not in the figure: x_H⁺ = x_e = XE0 = 2e-4 in the outer cloud (T. Abel, pers. comm., Nov 2026),
    # recombined for one local free-fall time inside (xefin below; XEFLAT=1 for the constant). He neutral (He⁺ = He²⁺ = 0).
    # IC=abn02red: the SECOND output ("nine Myrs later", RED dotted with circles; n_c = 6.5e5 cm⁻³), same pipeline
    # (lib/popiii/abn02red), T. Abel, 15 Nov 2026: the outer profile has stopped changing (it coincides with the final output
    # outside ~0.1 pc) and the start is free of the final output's strong initial transient. Differences from abn02fin:
    # (i) CONSTANT density, temperature and composition from the end of the central plateau to r = 0: the plateau ends at the
    # last vertex of the red n(r) path that equals its central value (log r = −2.034 pc; M_B = 0.090 M☉ there); the core
    # density is panel B's mean density inside that radius (mass-continuous with panel B, n = ρ/(μ m_H) 0.15 dex above panel A's
    # plateau 10^5.81 at the panel-C μ = 1.23, see the README), T, f_H2, f_HI and x_e those at the plateau edge (T(M), f(M) are
    # flat there anyway), v_r homologous (∝ r) from the plateau edge to v = 0 at the centre; (ii) dark matter: the final-time
    # DM line (the figure shows DM only at z = 19 and at the final time; the final one is the closer in time, 9 Myr vs 200 yr
    # later, and the DM inside ~1 pc barely changes over the last ~0.3 Myr); x_e: the abn02fin recipe (xefin, XEFLAT).
    fin = ic in ("abn02fin", "abn02red"); red = ic == "abn02red"
    if fin
        function readfig2(csv)
            tab = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
            for l in eachline(csv)
                (startswith(l, "#") || startswith(l, "curve")) && continue
                f_ = split(l, ","); k_ = String(f_[1])
                haskey(tab, k_) || (tab[k_] = (Float64[], Float64[]))
                push!(tab[k_][1], parse(Float64, f_[2])); push!(tab[k_][2], parse(Float64, f_[3]))
            end
            tab
        end
        fintab = readfig2(get(ENV, "ABNFINCSV", joinpath(TF.ROOT, "lib", "popiii", "abn02fin", "abn02_final_fig2.csv")))
        redtab = red ? readfig2(get(ENV, "ABNREDCSV", joinpath(TF.ROOT, "lib", "popiii", "abn02red", "abn02_red_fig2.csv"))) : fintab
        fincurve(k_, t = redtab) = (xs = t[k_][1]; all(>(0), diff(xs)) || error("$ic: $k_ not monotone in x"); (xs, t[k_][2]))
        lr_tab, ln_tab = fincurve("n_r"); lrd, lρd = fincurve("dm_r", fintab); lrB_fin, lMB_fin = fincurve("M_r")
        # abn02red, DMLINE=z19 (default) | final: which DM line of panel A. The final-time one has contracted with the collapsed
        # baryons: inside r_p it is 20–35× the red core's gas density (3.2 M☉ of DM inside 9e-3 pc against 0.09 M☉ of gas), and
        # the core then collapses on the DM free-fall time (8e3 yr; run abnred6). The z = 19 line (blue solid, from 0.087 pc,
        # 10^4.38 GeV cm⁻³ there, ~2 % of the core gas density) is used instead, held at that density inside 0.087 pc (the
        # pchip's constant continuation: a uniform DM core, no cusp; with DMCORE no point mass inside 1e-4 pc either).
        if red
            dmline = get(ENV, "DMLINE", "z19"); dmline in ("z19", "final") || error("DMLINE=$dmline: z19 or final")
            dmline == "z19" && ((lrd, lρd) = fincurve("dm_z19_r"))
            println("abn02red: DM line ", dmline, " (log rho_DM = ", lρd[1], " GeV cm^-3 at log r = ", lrd[1], " pc and inside)")
        end
        lMT_fin, lT_fin = fincurve("T_M"); lMH_fin, lfH2_fin = fincurve("fH2_M"); lMI_fin, lfHI_fin = fincurve("fHI_M")
        lMv_fin, v_fin = fincurve("vr_M")
    elseif z19
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
        # DMCORE=1 (default for abn02red only): inside the table's first radius (1e-4 pc) the DM mass falls as r³ (uniform
        # density, as Mgd[1] assumes) instead of staying at M_DM(1e-4 pc) — that constant is a central POINT MASS (2e-8 M☉
        # for abn02z19, 4e-6 M☉ for abn02fin/red) which, on a mesh reaching 0.01 R☉, dominates the gas gravity at the centre
        # (abnred5: zone 1 to −6.6 km/s and 2400 K in two steps)
        if get(ENV, "DMCORE", ic == "abn02red" ? "1" : "0") == "1"
            let f0 = TX.EXTMASS[], r1 = rgd[1], m1 = 10^lMg[1]
                TX.EXTMASS[] = r -> (rv = ForwardDiff.value(r); rv < r1 ? m1 * (rv / r1)^3 * one(r) : f0(r))
            end
        end
    end
    # M(r) by integrating ρ(r) on a fine log grid
    rg = [pc * 10^x for x in range(-4.0, 1.5, length = 20001)]; Mg = zeros(length(rg))
    for j in 2:length(rg)
        rm = 0.5 * (rg[j] + rg[j-1]); Mg[j] = Mg[j-1] + 4π * rm^2 * nr(rm) * μ_abn * PC.mH * (rg[j] - rg[j-1])
    end
    Mg .+= 4π / 3 * rg[1]^3 * nr(rg[1]) * μ_abn * PC.mH
    if fin
        # panel-B M(r) (pchip in log M – log r), uniform density inside its first point; the n(r)·1.22 m_H integral above
        # is printed against it only (μ = 1.22 is atomic gas: the core is molecular, μ ≈ 2.26 there; see the README)
        MBfin = pchip(lrB_fin, lMB_fin); rB0 = pc * 10^lrB_fin[1]; ρin_fin = 10^lMB_fin[1] * Msun / (4π / 3 * rB0^3)
        if red                                                               # constant-density core out to the plateau edge
            kpl = findlast(==(ln_tab[1]), ln_tab); rB0 = pc * 10^lr_tab[kpl]; ρin_fin = 10^MBfin(lr_tab[kpl]) * Msun / (4π / 3 * rB0^3)
            println("abn02red: central plateau of n(r) to log r = ", lr_tab[kpl], " pc (vertex ", kpl, "): n_A = 10^", ln_tab[1],
                    ", M_B = ", round(10^MBfin(lr_tab[kpl]), sigdigits = 4), " Msun; uniform core rho = ", round(ρin_fin, sigdigits = 4), " g cm^-3")
        end
        rgA, MgA = rg, Mg
        global rg = pc .* 10 .^ collect(range(-10.0, lrB_fin[end], length = 60001))
        global Mg = [r < rB0 ? 4π / 3 * r^3 * ρin_fin : 10^MBfin(log10(r / pc)) * Msun for r in rg]
        all(>(0), diff(Mg)) || error("abn02fin: panel-B M(r) not increasing")
        for x in (-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0)
            jA = searchsortedfirst(rgA, pc * 10^x); jB = searchsortedfirst(rg, pc * 10^x)
            println("$ic: log r = $x pc  log M panel B ", round(log10(Mg[jB] / Msun), digits = 3), "  n(r) 1.22 m_H integral ",
                    1 < jA <= length(rgA) ? round(log10(MgA[jA] / Msun), digits = 3) : NaN)
        end
        println("$ic: inner uniform core n = ", round(ρin_fin / PC.mH, sigdigits = 4), " m_H cm^-3 (panel A n_c = 10^", ln_tab[1],
                ", i.e. mu_c = ", round(ρin_fin / (PC.mH * 10^ln_tab[1]), digits = 3), ")")
        Tfin_ = pchip(lMT_fin, lT_fin); fH2fin_ = pchip(lMH_fin, lfH2_fin); fHIfin_ = pchip(lMI_fin, lfHI_fin); vfin_ = pchip(lMv_fin, v_fin)
        mcore_fin = red ? 4π / 3 * rB0^3 * ρin_fin / Msun : 0.0             # (abn02red: T, f, x_e constant inside the plateau mass)
        Tfin(m) = 10^Tfin_(log10(max(m, mcore_fin)))                         # m in M☉ (constant below 0.01 M☉)
        # x_e(M) (T. Abel, 15 Nov 2026): the outer cloud's 2e-4 (XE0) recombined for one local free-fall time at the IC's
        # density and temperature, so the dense regions do not start with an artificially high x_e (and H₂ formation through
        # H⁻): dx/dt = −k2(T) n_H x² ⇒ x = x0 / (1 + k2 n_H x0 t_ff), t_ff = sqrt(3π/(32 G ρ)) (as the run's t_ff), k2 = AAZN97
        # case A (the network's H⁺ + e rate), n_H = ρ/μ_H, ρ from panel B's dM/dr. Made non-decreasing outward
        # (x(m) = min over m' ≥ m) so wiggles of the differentiated panel-B mass cannot invert it (only the innermost ~0.01 M☉
        # is affected). XEFLAT=1: constant XE0 everywhere (the 15 Nov morning IC, runs abnfin*, z100fin).
        xefin_flat = get(ENV, "XEFLAT", "0") == "1"; xefin_0 = parse(Float64, get(ENV, "XE0", "2e-4"))
        lmxe_fin, lxe_fin = let jj = 2:length(Mg)-1
            xs = [begin ρ_ = (Mg[j+1] - Mg[j-1]) / (4π * rg[j]^2 * (rg[j+1] - rg[j-1])); m_ = Mg[j] / Msun
                      xefin_0 / (1 + PC.k2_aazn(Tfin(m_)) * (ρ_ / ChemNet.μH()) * xefin_0 * sqrt(3π / (32PC.G * ρ_))) end for j in jj]
            for j in length(xs)-1:-1:1; xs[j] = min(xs[j], xs[j+1]); end
            (log10.(Mg[jj] ./ Msun), log10.(xs))
        end
        xefin(m) = xefin_flat ? xefin_0 : 10^(let x = log10(m); x <= lmxe_fin[1] ? lxe_fin[1] : x >= lmxe_fin[end] ? lxe_fin[end] :
            (j = searchsortedlast(lmxe_fin, x); lxe_fin[j] + (x - lmxe_fin[j]) / (lmxe_fin[j+1] - lmxe_fin[j]) * (lxe_fin[j+1] - lxe_fin[j])) end)
        for lm_ in (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0)
            println("$ic: x_e(M = 10^", lm_, " Msun) = ", round(xefin(10^lm_), sigdigits = 3), xefin_flat ? " (XEFLAT)" : " (recombined for one t_ff from XE0)")
        end
        function finy(m)                                                     # (x_H⁺, x_H⁻, x_H₂, x_H₂⁺) per H nucleus at m [M☉]
            m = max(m, mcore_fin); f2 = 10^fH2fin_(log10(m)); f1 = 10^fHIfin_(log10(m))
            (xefin(m), 0.0, 0.5 * f2 / (f2 + f1), 0.0)
        end
        mv1_fin = max(10^lMv_fin[1], mcore_fin); vv1_fin = vfin_(log10(mv1_fin))   # (abn02fin: the first table point, v_fin[1])
        rM1fin = let m1 = mv1_fin * Msun, j = clamp(searchsortedfirst(Mg, m1), 2, length(Mg))
            rg[j-1] + (m1 - Mg[j-1]) / (Mg[j] - Mg[j-1]) * (rg[j] - rg[j-1])
        end
        vfin(m, r) = m >= mv1_fin ? 1e5 * vfin_(log10(m)) : 1e5 * vv1_fin * r / rM1fin   # [cm/s]; m [M☉], r [cm]
    end
    M = parse(Float64, get(ENV, "MCLOUD", z19 ? "1e4" : fin ? "2e4" : "1000")) * Msun
    if z19
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
    Tz = fin ? [Tfin(m) for m in mz] : [interp(log10(m), lM_tab, T_tab) for m in mz]
    fz = fin ? [10^fH2fin_(log10(m)) for m in mz] : [10^interp(log10(m), lM_tab, lf_tab) for m in mz]   # (fin: mass fraction, unused below)
    uI = fin ? [vfin(mI[j] / Msun, rI[j]) for j in 1:nsh+1] : ic == "abn02z19" ? [m < Mjoin ? 0.0 : 1e5 * interp(log10(m / Msun), lM_tab, v_tab) * min(1.0, (m - Mjoin) / Mjoin) for m in mI] :
         [1e5 * interp(log10(max(m / Msun, 1e-2)), lM_tab, v_tab) for m in mI]; uI[1] = 0.0
    ρc = rho[1]; Tc0 = Tz[1]
    xe0 = parse(Float64, get(ENV, "XE0", z19 ? "1e-4" : fin ? "2e-4" : "1e-7"))   # x_e [assumed; ABN02 does not show it; fin: T. Abel, outer cloud]
    yz = fin ? [finy(m) for m in mz] : [(xe0, 0.0, f / 2, 0.0) for f in fz]
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
const NIMP = Ref(0)                                                  # (kept in the checkpoint; the IMPCOOL slope it counted is gone)
const CHEM = get(ENV, "CHEM", "netmodernli")     # net / netmodern / netmodernli: NetGas (G5); the RHFC02 PrimordialGas material is gone
CHEM in ("net", "netmodern", "netmodernli") || error("CHEM=$CHEM: the network material net, netmodern or netmodernli")
begin
    # G5: ChemNet network; equilibrium table where the ideal-gas pressure ≥ PSWITCH (EOS validity, RHFC02's 1e4 dyn cm⁻²)
    netw, phot = CHEM == "net" ? (ChemRates.rhfc02(; k5 = k5f), ChemNet.Photo[]) : ChemRates.modern(; li = CHEM == "netmodernli")
    # SL3=1 (with STARLIGHT=1): three-level H atom in the lit network zones (PM.tl3_reactions; rates are zero elsewhere)
    get(ENV, "SL3", "0") == "1" && (netw = PM.tl3_network(netw))
    ChemNet.HM88 && (netw = ChemNet.hm88_network(netw))      # HM88=1: H recombination × Z_w/Z_el (ChemNet.HM88)
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
TX.AV_TENSION[] = get(ENV, "AVTENSION", "1") == "1"       # 1 (TITAN) | 0 | f in (0,1) | radial | radial,f  (TX.AV_TFRAC)
let p = split(get(ENV, "AVTENSION", "1"), ","), f = tryparse(Float64, p[end])
    TX.AV_TRAD[] = p[1] == "radial"; TX.AV_TFRAC[] = f !== nothing && 0 < f < 1 ? f : 0.0
end
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
        println(io, "step,t_yr,dt_s,wall_s,pid,nsteps,t_step,t_chem,t_lines,t_other,gc_s,alloc_GB,ndiv,nrepeat,tr_steps,fail_tr,fail_bdf,be_fallback,nali,nchemfb,nthreads,nz,log_nc,Tc,mstar,rstar,umin_kms,dtlim,etot_Gmax,etot_vmax,etot_zone")
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
# ETOTMON=1 (diagnostic, default off; no influence on the solution): the total-energy monitor of every converged implicit solve, TX.etot_monitor, via TD.STAGE_HOOK.
# `buf`: the stages since the last trace line (DTERRTRACE), `pG/pv/pz`: the largest |G|, the largest per-zone value and its zone since the last perf.csv row.
mutable struct EtotMon; on::Bool; buf::Vector{Any}; pG::Float64; pv::Float64; pz::Int; end
const EM = EtotMon(false, Any[], 0.0, 0.0, 0)
"The ETOT text of the DTERR trace line: the stages since the last line (G = global total-energy residual / energy throughput of a solve, see TX.etot_monitor)."
function etot_text()
    isempty(EM.buf) && return ""
    b_ = EM.buf; EM.buf = Any[]
    Gl = b_[end][4].G; Σi = sum(x[4].imb for x in b_); Σt = sum(x[4].thr for x in b_)
    iG = argmax([abs(x[4].G) for x in b_]); iv = argmax([x[4].vmax for x in b_]); is = argmax([x[4].vshare for x in b_])
    @sprintf(" | ETOT %d solves: G(last %s) %+.2e, sum %+.2e, max|G| %.2e (step %d %s); worst zone/own throughput z%d %.2e (step %d %s), largest share z%d %.2e (step %d %s)",
             length(b_), b_[end][2], Gl, Σi / max(Σt, 1e-300), abs(b_[iG][4].G), b_[iG][1], b_[iG][2], b_[iv][4].zmax, b_[iv][4].vmax, b_[iv][1], b_[iv][2],
             b_[is][4].zshare, b_[is][4].vshare, b_[is][1], b_[is][2])
end
"One perf.csv row (called at the TIMING report, before the accumulators are reset)."
function recperf!(k, tnow)
    REC.io[] === nothing && recopen!(k)
    gct = Base.gc_num().total_time; gcb = Base.gc_bytes()
    st = Trbdf2.STATS; un_ = TF.var(:un)[a:b+1]
    oth = time() - tnow + TOTH[] - TCHEM[] - TLINE[]
    ns = REC.klast[] < 0 ? 0 : k - REC.klast[]
    open(joinpath(RUNREC, "perf.csv"), "a") do io
        @printf(io, "%d,%.8e,%.4e,%.1f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%.6e,%.6g,%.6g,%.3f,%s,%.3e,%.3e,%d\n",
                k, TF.getscalar(:timen) / 3.156e7, TF.getscalar(:dtime), time() - REC.wall0, getpid(), ns, TSTEP[], TCHEM[], TLINE[], oth,
                (gct - REC.gct[]) / 1e9, (gcb - REC.gcb[]) / 2^30, TD.NDIV[], TD.NREP[], get(st, :steps, 0), get(st, :fail_tr, 0), get(st, :fail_bdf, 0),
                get(st, :be_fallback, 0), NALICALL[], PM.NCHEMFB[], Threads.nthreads(), b - a + 1, log10(TF.var(:dn)[a] / PC.mH), TF.var(:tn)[a],
                REC.mstar[], REC.rstar[], minimum(un_) / 1e5, REC.dtlim[], EM.pG, EM.pv, EM.pz)
    end
    EM.pG = 0.0; EM.pv = 0.0; EM.pz = 0
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
    PM.STARLIGHT[].tl3 && (TX.CQB[] = zeros(length(mat.x)); TX.YQB[] = zeros(length(mat.x)))
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
    PM.SLBLEND[] = get(ENV, "SLBLEND", "0") == "1"
    PM.NLOPACJ[] = PM.NLOPAC[] && get(ENV, "NLOPACJ", "0") == "1"
    PM.NLOPACJ[] && pr("NLOPACJ: kappa_E of network zones weighted with the star light's B_nu(T_c)\n")
    PM.SLEB[] = PM.SLBLEND[] && get(ENV, "SLEB", "0") == "1"
    PM.SLEB[] && pr("SLEB: star-light colour temperature T_r at tau + 2/3 (Eddington-Barbier), not the photosphere's T_eff\n")
    PM.SLBLEND[] && pr("SLBLEND: star-light colour temperature T_eff (tau <= 0.3) to T_r (tau >= 3), smooth in log tau\n")
    PM.PSWSTAY[] = !(PM.STARLIGHT[].tl3 && get(ENV, "PSWSTAY", "1") == "0")
    PM.PSWSTAY[] || pr("PSWSTAY=0: below the pressure switch a table zone stays on the table only by the tau rule\n")
    PM.NLOPAC[] && pr("NLOPAC: non-LTE opacity in the three-level network zones (no H bound-free in kappa_P; n=2 of the three-level atom in kappa_R)\n")
    PM.NLOPAC[] && PM.NLOIT && pr("NLOIT: network-zone opacities from the Newton iterate's own populations (table x continuum-code ratio at the iterate; NLOIMPT, NLODAMP, CPLNLOX, NLOFIX, NLOTMAX inactive there)\n")
    !PM.NLOPAC[] && PM.NLOIT && error("NLOIT=1 needs NLOPAC=1 (STARLIGHT=1 SL3=1)")
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
# LINEIMP (always on since the clean-up; the switch is gone): the line cooling of each zone is a function of the zone's state in the chemistry and in
# the hydro Newton (PopIIIMaterial.LINEIMP: populations re-solved at the iterate's T with the radiative rates β, J̄ of the last transfer
# call), not Λ held from the last call. Retires as devices: the IMPCOOL slope (second transfer call) and its [−10, 30] clamp, the
# Δt ≤ 0.1 (or DTCOOL) t_cool cap, TCOOLTH (only the cap used it). LINEFRAC still skips transfer calls (the model is then kept).
# Hybrid NLTE: on an ALI step the zone takes the ALI rates (β = 1 − Λ*, J̄ = J̄_fs − Λ* S); in between the EP rates of the step × the
# last ALI/EP factor (the EP rates themselves follow the iterate's populations, velocity gradient and incident field: PopIIIMaterial.li_ep), bounded to [1/LIFMAX, LIFMAX] (default 10); where LINECOLL holds the ALI collisional rate (ALIFIX: EP and ALI differ
# in sign), the zone holds the last ALI call's RATES instead, so its cooling still follows its T. Needs NLTE = ep | ali | hyb.
const LINEIMP = true; const LIFMAX = parse(Float64, get(ENV, "LIFMAX", "10"))
LINEIMP && !(NLTE in ("ep", "ali", "hyb")) && error("LINEIMP=1 needs NLTE=ep, ali or hyb")
LINEIMP && (PM.LINEIMP[] = true; pr("LINEIMP: H2 line cooling of each zone from its own state (populations at the iterate's T, radiative rates of the last transfer call); no t_cool cap, no IMPCOOL slope\n"))
const LIALI = Ref{Vector{Any}}(Any[]); const LISTAT = zeros(Int, 4)    # (zones with EP, ALI, held-ALI rates; f clamped)
const LIJCMB = Ref{Vector{Float64}}(Float64[])
"(diagnostic, NANSCAN=1) report the first zones whose state is not finite."
const NANSCAN = get(ENV, "NANSCAN", "0") == "1"; const FLUSHOUT = get(ENV, "FLUSHOUT", "0") == "1"
function nanscan(where, k)
    n = 0
    for nm in (:tn, :dn, :egn, :ern, :un, :frn, :rn, :pgn, :xkpn, :chifn)
        v = TF.var(nm)
        for kk in a:b+1
            isfinite(v[kk]) && continue
            n += 1; n <= 12 && pr("NANSCAN step %d %s: %s[%d] (zone %d) = %s; T %.4g rho %.4g eq %s x %s\n", k, where, String(nm), kk, kk - a + 1, string(v[kk]),
                                  TF.var(:tn)[min(kk, b)], TF.var(:dn)[min(kk, b)], string(mat.eqzone[min(kk, b)]), string(round.(mat.x[min(kk, b)], sigdigits = 3)))
        end
    end
    for kk in a:b
        all(isfinite, mat.x[kk]) && continue
        n += 1; n <= 12 && pr("NANSCAN step %d %s: composition of zone %d not finite: %s (T %.4g)\n", k, where, kk - a + 1, string(mat.x[kk]), TF.var(:tn)[kk])
    end
    nothing
end
"(LINEIMP) zones remapped: drop the line models (the held Λext of the remap is used until the next transfer call)."
li_inval!() = (PM.li_invalidate!(); fill!(LIALI[], nothing); nothing)
"(LINEIMP) per TITAN zone, the line model from the transfer call `Lc` (H2Lines.LAST[]) on the reduced grid (see LINEIMP above)."
function li_store!(Lc, fac, dd, nh2, nl, zk, ncz, nz_; ali = false, held = nothing, ep = nothing)
    hasproperty(Lc, :β) || return nothing
    nk = length(mat.x); length(PM.LIMODEL[]) == nk || (PM.LIMODEL[] = Any[nothing for _ in 1:nk]); length(LIALI[]) == nk || (LIALI[] = Any[nothing for _ in 1:nk])
    for kk in a:b; PM.LIMODEL[][kk] = nothing; end
    if isempty(LIJCMB[])
        Tb = parse(Float64, get(ENV, "TCMB", "0"))
        LIJCMB[] = [Tb > 0 ? H2Lines.Bν(ν_, Tb) : 0.0 for ν_ in HLL.ν]
    end
    jc_ = LIJCMB[]; fill!(LISTAT, 0)
    for j in ncz+1:nz_
        kk = zk(j)
        if held !== nothing && held(j) && LIALI[][kk] !== nothing
            PM.LIMODEL[][kk] = LIALI[][kk]; LISTAT[3] += 1; continue
        end
        β_ = Lc.β[:, j]; J_ = max.(Lc.Jeff[:, j], β_ .* jc_)
        f_ = fac(j); fc = clamp(f_, 1 / LIFMAX, LIFMAX); fc != f_ && (LISTAT[4] += 1)
        PM.li_set!(kk, HLV, β_, J_, nl.nH[j], nh2[j], nl.nHe[j], nl.nHp[j], dd[j], fc;
                   ep = ep === nothing ? nothing : (; jc = jc_, P = Lc.pops[:, j], T = ep.T[j], g0 = PM.sl_grad(ep.r[j], ep.r[j+1], ep.u[j], ep.u[j+1])))
        LISTAT[ali ? 2 : 1] += 1
        ali && (LIALI[][kk] = PM.LIMODEL[][kk])
    end
    nothing
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
        if NLTE != "hyb"
            out_ = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = NLTE == "ali" ? :ali : :ep, inc)
            LINEIMP && k >= 0 && li_store!(H2Lines.LAST[], j -> 1.0, dd, nh2, nl, zk, ncz, nz_; ep = NLTE == "ep" ? (; T = tt, r = rr, u = uu) : nothing)
            return out_
        end
        tq = time(); Λe, Lr, x = H2Lines.cooling(HLL, rr, uu, dd, tt, nh2; nbin = NBIN, nlte = nl, nlte_mode = :ep, inc)
        LEP = (LINEIMP && k >= 0) ? H2Lines.LAST[] : nothing
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
            LINEIMP && li_store!(H2Lines.LAST[], j -> 1.0, dd, nh2, nl, zk, ncz, nz_; ali = true)
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
            if Trbdf2.SDON       # stalldiag: line-cooling ALI factor change per TITAN zone in SDZONES at this call
                for j in 1:nz_
                    (zk(j) - Int(TF.getscalar(:ngrs)) + 1) in Trbdf2.SDZ || continue
                    pr("        SDALI jstep %d TITAN zone %d reduced %d: f %.4g -> %.4g, Lambda_ALI %.4e Lambda_EP %.4e\n", Int(TF.getscalar(:jstep)), zk(j) - Int(TF.getscalar(:ngrs)) + 1, j,
                       length(ALIFAC[]) == nz_ ? ALIFAC[][j] : NaN, fnew[j], Λa[j], Λe[j])
                end
            end
            ALIFAC[] = fnew; ALILAST[] = k; NALICALL[] += 1
            return Λa, Lra, xa
        end
        if LEP !== nothing                  # LINEIMP: EP rates × the ALI/EP factor; zones that hold the ALI collisional rate (ALIFIX) hold the ALI rates instead
            hf = length(ALIFAC[]) == nz_
            li_store!(LEP, j -> hf ? ALIFAC[][j] : 1.0, dd, nh2, nl, zk, ncz, nz_; ep = (; T = tt, r = rr, u = uu),
                      held = (LINECOLL && hf && length(ALIFIX[]) == nz_) ? (j -> !isnan(ALIFIX[][j])) : nothing)
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
# LINEIMP line models as plain data (the level data LV is the constant HLV, the LiCache is rebuilt empty): a checkpoint stays readable without the modules
lim_plain(m) = m === nothing ? nothing : (; β = copy(m.β), J = copy(m.J), nH = m.nH, nH2 = m.nH2, nHe = m.nHe, nHp = m.nHp, ρ = m.ρ, f = m.f)
lim_make(p) = p === nothing ? nothing : (; LV = HLV, β = copy(p.β), J = copy(p.J), nH = p.nH, nH2 = p.nH2, nHe = p.nHe, nHp = p.nHp, ρ = p.ρ, f = p.f,
                                         c = PM.LiCache(NaN, 0.0, 0.0, false, H2Levels.PopWS(HLV)))
matstate() = mat isa PM.NetGas ? (; x = deepcopy(mat.x), eqzone = copy(mat.eqzone), Λext = copy(mat.Λext), Tref = copy(mat.Tref), slope = copy(mat.slope), Y = HOTON ? copy(SG[].Y) : nothing,
                                     star = HOTON ? (; fscale = copy(SG[].fscale), fslope = copy(SG[].fslope), fcur = copy(SG[].fcur), Dmix = copy(SG[].Dmix), dtcap = SG[].dtcap, capE = SG[].capE, ledger = SG[].ledger, Yprev = copy(SG[].Yprev)) : nothing,
                                     # lagged over one step (set in chemistry, used by the next hydro step): without them a restart is not exact
                                     lag = (; NLOPC = deepcopy(PM.NLOPC[]), GASQ = copy(TX.GASQ[]), RAYQ = copy(TX.RAYQ[]), CQB = copy(TX.CQB[]), YQB = copy(TX.YQB[]), SLDIAG = deepcopy(PM.SLDIAG[]), CHEMB = copy(TX.CHEMB[]),
                                                 LIMODEL = [lim_plain(m_) for m_ in PM.LIMODEL[]], LIALI = [lim_plain(m_) for m_ in LIALI[]], NLOST = deepcopy(PM.NLOST[]))) :   # (lag.LIMODEL/LIALI: LINEIMP's H2 line models / held ALI rates, set by the last transfer call, used by the next hydro step)
             error("checkpoint: only the NetGas material (CHEM=net*) is supported")
# jinc: 1 for a checkpoint written inside the step loop's output callback, where TITAN's step counter jstep has not yet been advanced for the
# step just done (run_julia! increments it after output()); the restart then continues with jstep + 1, as the uninterrupted run would
# HELDDUMP=1 (diagnostic, default off): with every checkpoint also write <name>.held, the numeric content of every global (arrays, Refs, tuples, dicts) of the
# solver modules and of this driver, to find the Julia-side state that a restart does not reproduce (tools/ckheld.jl compares two such files)
function held_flat!(out, name, v, depth = 0)
    depth > 6 && return
    if v isa Number; v isa Bool || v isa Integer || v isa AbstractFloat ? push!(out, name => Float64[v]) : nothing
    elseif v isa AbstractArray{<:Number}; length(v) <= 20_000_000 && push!(out, name => Float64.(vec(collect(v))))
    elseif v isa Ref; isassigned(v) && held_flat!(out, name, v[], depth + 1)
    elseif v isa AbstractArray; length(v) <= 200_000 && for (i, e) in enumerate(v); held_flat!(out, "$name[$i]", e, depth + 1); end
    elseif v isa Union{Tuple,NamedTuple}; for (i, e) in enumerate(v); held_flat!(out, "$name.$i", e, depth + 1); end
    elseif v isa AbstractDict; length(v) <= 10_000 && for kk in sort(collect(keys(v)); by = string); held_flat!(out, "$name{$kk}", v[kk], depth + 1); end
    end
    nothing
end
function held_dump(path)
    out = Pair{String,Vector{Float64}}[]
    mods = Pair{String,Module}["Main" => Main]      # Main itself and every module that was include()d into it (not the registered packages)
    for nm in names(Main; all = true)
        isdefined(Main, nm) || continue; v = getfield(Main, nm)
        v isa Module && v !== Main && v !== Base && v !== Core && pathof(v) === nothing && push!(mods, string(nm) => v)
    end
    for (mn, m) in mods
        for nm in names(m; all = true)
            s_ = string(nm); (startswith(s_, "#") || nm in (:ans, :ARGS, :Core, :Base, :Main, :ENV, :TF)) && continue
            isdefined(m, nm) && !isconst(m, nm) && m === Main && continue
            v = try getfield(m, nm) catch; continue end
            (v isa Module || v isa Function || v isa Type) && continue
            try held_flat!(out, "$mn.$s_", v) catch; end
        end
    end
    serialize(path, Dict(out))
end
function checkpoint(name, k; jinc = 0)
    get(ENV, "HELDDUMP", "0") == "1" && held_dump(joinpath(CKDIR, "m2c_$(name).held"))
    mkpath(CKDIR)
    TF.save_state(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    ax = (; k, jnext = Int(TF.getscalar(:jstep)) + jinc, mat = matstate(), ALIFAC = copy(ALIFAC[]), ALIFIX = copy(ALIFIX[]), ALILAST = ALILAST[], ALIINT = ALIINT[],
        NALICALL = NALICALL[], NIMP = NIMP[], LLINE = LLINE[], nextln = nextln[], SNAPNEXT = SNAPNEXT[], SNAPS,
        trb = (; STATS = copy(Trbdf2.STATS), BE_LEFT = Trbdf2.BE_LEFT[]),
        # solver / controller / held-input state that is not in the TITAN COMMON (without it the first steps after a restart differ: no FSAL error estimate,
        # no PI history, no mesh feature tracking, a restarting table-to-network ramp)
        ctl = (; HNEXT = Trbdf2.HNEXT[], ERRPREV = Trbdf2.ERRPREV[], HEUR = Trbdf2.HEUR[], LASTACC = Trbdf2.LASTACC[], DSAL = deepcopy(Trbdf2.DSAL[]),
               ESTATS = copy(Trbdf2.ESTATS), LIMHIST = copy(Trbdf2.LIMHIST)),
        edd = (; HEDD = copy(TD.HEDD), HXEDD = copy(TD.HXEDD), FXEDD = copy(TD.FXEDD), TIEW = copy(TD.TIEW), HPLUS = copy(TC.HPLUS), GXEDD = TD.GXEDD[]),   # Eddington closure arrays held between steps
        geq = (; J0 = GEQJ0[], PREV = deepcopy(GEQPREV[]), HELD = deepcopy(TX.GRIDEQ[])), pswramp = (; J0 = PM.RAMPJ0[], K0 = PM.RAMPK0[]), PEXT = TX.PEXT[], env = Dict(k_ => v_ for (k_, v_) in ENV if k_ in CKENV))
    TX.MG[] && (ax = merge(ax, (; mg = (; edges = copy(TX.MGEDGE[]), E = deepcopy(TX.MGEN[]), F = deepcopy(TX.MGFN[])))))   # MGROUPS: the group fields
    serialize(joinpath(CKDIR, "m2c_$(name).aux"), ax)
end
const GEQHELD = Ref{Any}(nothing)
const JNEXT = Ref(-1); const CKAUX = Ref{Any}(nothing)
function edd_restore!(e)
    for (dst, src) in ((TD.HEDD, e.HEDD), (TD.HXEDD, e.HXEDD), (TD.FXEDD, e.FXEDD), (TD.TIEW, e.TIEW), (TC.HPLUS, e.HPLUS))
        resize!(dst, length(src)); copyto!(dst, src)
    end
    TD.GXEDD[] = e.GXEDD
end
const MGAUX = Ref{Any}(nothing)       # (MGROUPS) the group fields of the checkpoint a run restarted from
# settings that must agree between a run and its restart (the set-up and the physics); the rest may change
const CKENV = ("IC", "TCMB", "AVZONE", "AVQL1", "AVCQ2", "CHEM", "NLTE", "NALI", "NBIN", "K5", "MCLOUD", "LRJOIN", "DLNJOIN", "DM", "PSWITCH",
               "INTEGRATOR", "THETA", "XE0", "RAD")   # (not CDYN: the accretion phase relaxes the dynamical Δt limit at a restart)
function restart!(name; allow = ())
    aux = deserialize(joinpath(CKDIR, "m2c_$(name).aux"))
    for (k_, v_) in aux.env
        (k_ in allow || !(k_ in CKENV)) && continue        # (a key whose switch was removed in the clean-up: its value is the method now)
        get(ENV, k_, nothing) == v_ || error("restart: ENV $k_ = $(repr(get(ENV, k_, nothing))) differs from the checkpoint's $(repr(v_))")
    end
    for k_ in CKENV; k_ in allow && continue; haskey(ENV, k_) && !haskey(aux.env, k_) && error("restart: ENV $k_ is set but was not in the checkpoint's run"); end
    TF.load_state!(joinpath(CKDIR, "m2c_$(name).tstate.gz"))
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1   # (zone count may differ after rezoning)
    s = aux.mat; mat.x .= deepcopy(s.x); mat.eqzone .= s.eqzone; mat.Λext .= s.Λext; mat.Tref .= s.Tref; mat.slope .= s.slope
    ChemNet.HECHEM && mat isa PM.NetGas && PM.he_pad!(mat)                    # (restart from a checkpoint written without HECHEM)
    HOTON && hasproperty(s, :Y) && s.Y !== nothing && (SG[].Y .= s.Y)      # (a checkpoint from before HOT=1 keeps the initial abundances)
    if HOTON && hasproperty(s, :star) && s.star !== nothing
        t_ = s.star; g_ = SG[]; g_.fscale .= t_.fscale; g_.fslope .= t_.fslope; g_.fcur .= t_.fcur; g_.Dmix = t_.Dmix; g_.dtcap = t_.dtcap; g_.capE = t_.capE; g_.ledger = t_.ledger; g_.Yprev .= t_.Yprev
    end
    if hasproperty(s, :lag)                        # (older checkpoints: these start empty, so the first step after a restart differs)
        l_ = s.lag; PM.NLOPC[] = deepcopy(l_.NLOPC); PM.SLDIAG[] = deepcopy(l_.SLDIAG)
        PM.NLOIT && PM.nlost_restart!(mat)       # (older checkpoints: only the zones with non-LTE factors (NLOPC) get their composition and three-level parameters back)
        # NLOST (composition and three-level parameters at the end of the last chemistry = the populations of the next step's iterate-based opacity): the
        # rebuild above misses every zone when NLOPC is empty (NLOIT runs), so the first step after a restart had the plain table opacity in all network zones
        hasproperty(l_, :NLOST) && (PM.NLOST[] = Any[deepcopy(e_) for e_ in l_.NLOST])
        length(l_.GASQ) == length(TX.GASQ[]) && (TX.GASQ[] .= l_.GASQ); length(l_.RAYQ) == length(TX.RAYQ[]) && (TX.RAYQ[] .= l_.RAYQ)
        hasproperty(l_, :CQB) && length(l_.CQB) == length(TX.CQB[]) && (TX.CQB[] .= l_.CQB; TX.YQB[] .= l_.YQB)
        if !hasproperty(l_, :CQB) && !isempty(TX.CQB[]) && PM.SLDIAG[] !== nothing      # checkpoint from before CQB/YQB: cq, yq were booked inside RAYQ/GASQ; separate them (SLDIAG keeps them per zone)
            nconv = 0
            for (k_, t_) in enumerate(PM.SLDIAG[].tl3)
                t_ === nothing && continue
                TX.CQB[][k_] = t_.chemq; TX.YQB[][k_] = t_.cycq; TX.RAYQ[][k_] -= t_.chemq; TX.GASQ[][k_] -= t_.cycq; nconv += 1
            end
            pr("restart: checkpoint without CQB/YQB; chemistry and cycle heat of %d zones taken out of RAYQ/GASQ and booked as CQB/YQB\n", nconv)
        end
        if hasproperty(l_, :LIMODEL)
            PM.LIMODEL[] = Any[lim_make(m_) for m_ in l_.LIMODEL]; LIALI[] = Any[lim_make(m_) for m_ in l_.LIALI]
        else
            pr("WARNING: checkpoint %s has no H2 line models (older code): the first step after the restart has no line cooling until the first transfer call\n", name)
        end
        hasproperty(l_, :CHEMB) && length(l_.CHEMB) == length(TX.CHEMB[]) && (TX.CHEMB[] .= l_.CHEMB)
    end
    ALIFAC[] = copy(aux.ALIFAC); ALIFIX[] = hasproperty(aux, :ALIFIX) ? copy(aux.ALIFIX) : Float64[]; ALILAST[] = aux.ALILAST; ALIINT[] = aux.ALIINT; NALICALL[] = aux.NALICALL; NIMP[] = aux.NIMP
    LLINE[] = aux.LLINE; nextln[] = aux.nextln; SNAPNEXT[] = aux.SNAPNEXT
    SNAPS === nothing || aux.SNAPS === nothing || (empty!(SNAPS); append!(SNAPS, aux.SNAPS))
    empty!(Trbdf2.STATS); merge!(Trbdf2.STATS, aux.trb.STATS); Trbdf2.BE_LEFT[] = aux.trb.BE_LEFT
    if hasproperty(aux, :ctl)
        c_ = aux.ctl; Trbdf2.HNEXT[] = c_.HNEXT; Trbdf2.ERRPREV[] = c_.ERRPREV; Trbdf2.HEUR[] = c_.HEUR; Trbdf2.LASTACC[] = c_.LASTACC; Trbdf2.DSAL[] = deepcopy(c_.DSAL)
        empty!(Trbdf2.ESTATS); merge!(Trbdf2.ESTATS, c_.ESTATS); empty!(Trbdf2.LIMHIST); merge!(Trbdf2.LIMHIST, c_.LIMHIST)
    else
        pr("WARNING: checkpoint %s has no solver-controller state (written by older code): the first step after the restart has no TR-BDF2 error estimate and the controller history is empty\n", name)
    end
    if hasproperty(aux, :geq)
        GEQJ0[] = aux.geq.J0; GEQPREV[] = deepcopy(aux.geq.PREV)
        # the grid equation's held monitor (tauv, sx) was computed in the checkpointed step before chem!'s EOS/opacity refresh: not reproducible from the final state
        GEQHELD[] = aux.geq.HELD
         PM.RAMPJ0[] = aux.pswramp.J0; PM.RAMPK0[] = aux.pswramp.K0
    else
        pr("WARNING: checkpoint %s has no grid-equation / PSWRAMP state (older code): the mesh feature tracking restarts empty and the ramps restart\n", name)
    end
    TX.PEXT[] = aux.PEXT; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    hasproperty(aux, :mg) && (MGAUX[] = aux.mg)
    TF.setscalar!(:tmax, 1e3tff)                   # a run stopped by TSTOP/NSTOP has tmax = its last time
    CKAUX[] = aux
    JNEXT[] = hasproperty(aux, :jnext) ? aux.jnext : -1     # (set after the grideq_step! below; older checkpoints: jstep stays that of the last step, one behind the uninterrupted run)
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
const XMEXT = Ref{Any}(nothing)
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
    TC.ER_TAUMIN[] = parse(Float64, get(ENV, "ERTAUMIN", "1.0"))   # E_r of the optically thin envelope does not limit Δt
    if HOTON                                          # convection (hot branch only: PM.material_thermo), as the S5 stars
        TX.CONVECTION[] = true; TX.ALPHA_MLT[] = parse(Float64, get(ENV, "ALPHAMLT", "1.5")); TX.MLT_HSDISC[] = true; TX.MLT_RCAP[] = true
        PM.COLDMLT[] = parse(Float64, get(ENV, "COLDMLT", "0")); PM.COLDMLT[] > 0 && (TX.MLT_NOINFALL[] = true)
    end
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
    for k in b+1:length(mat.x); mat.x[k] = copy(mat.x[b]); mat.eqzone[k] = mat.eqzone[b]; mat.Λext[k] = 0.0; end; LINEIMP && li_inval!()
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
PM.PSWHYS[] = parse(Float64, get(ENV, "PSWHYS", "1"))
PM.PSWTMIN[] = parse(Float64, get(ENV, "PSWTMIN", "0"))
PM.PSWTEX[] = parse(Float64, get(ENV, "PSWTEX", "0")); PM.PSWTAU[] = parse(Float64, get(ENV, "PSWTAU", "0"))
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
const ENVLIN = parse(Int, get(ENV, "ENVLIN", "1"))   # with a regular centre (ENVEUL=0,...): ENVLIN zones of equal Δr from r = 0 to ENVCORE's first interface (default 1: one central zone)
const ENVCONS = get(ENV, "ENVCONS", "0") == "1"       # the ENVEUL remap takes the zone densities from the source's enclosed mass (mass-conserving); default: log-interpolated at the zone centres
const ENVSINKTAU = parse(Float64, get(ENV, "ENVSINKTAU", "1e6"))
"Grey optical depth from the outer edge at the interfaces i0..b+1 (the τ term of the GRIDEQ monitor)."
function gm_tau(i0)
    ch, dd, rr = TF.var(:chifn), TF.var(:dn), TF.var(:rn); t = zeros(b + 2 - i0); τ = 0.0
    for k in b:-1:i0; τ += ch[k] * dd[k] * (rr[k+1] - rr[k]); t[k-i0+1] = τ; end
    t
end
xmi(k) = ENVEUL === nothing ? TX.XMINT[][k] : TF.getscalar(:xmtot) - TF.var(:xmen)[k]
const XMEPREV = Ref(Float64[]); const ENVMS0 = Ref(0.0); const ENVMDOT = Ref(0.0)
envmenc() = TF.getscalar(:xmtot) .- TF.var(:xmen)[a:b+1]              # interior mass per steradian at the interfaces
include(joinpath(@__DIR__, "grideq_driver.jl"))                         # GRIDEQ (default off)
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
COUPLED (the method; the switch is gone): x_H⁺ and x_H₂ of every network zone are unknowns of the hydro Newton (TitanResidFast.COUPLED,
PopIIIMaterial cpl_*; implies CHEMCK). The chemistry after the step keeps only what the hydro solve does not carry (HD, the three-level
bookkeeping, the table switch); table and hot zones keep the split step.
"""
const COUPLED = true
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
"ln r of the fixed mesh's interfaces from the first interface (r_in, or ENVCORE's first interface for a regular centre) to rout (ENVEUL/ENVCORE spacing)."
function envlr(rout)
    rin, dl1, dl2, rb = ENVEUL[1] * Rsun, ENVEUL[2], ENVEUL[3], ENVEUL[4] * Rsun
    rcb = ENVCENTRE ? ENVCORE[3] * Rsun : 0.0
    lr = [log(ENVCENTRE ? ENVCORE[1] * Rsun : rin)]; dl = ENVCENTRE ? ENVCORE[2] : dl1
    gout_ = length(ENVEUL) >= 5 ? ENVEUL[5] : 1.02; gin_ = length(ENVEUL) >= 6 ? ENVEUL[6] : 0.95
    while lr[end] < log(rout)
        r_ = exp(lr[end])
        r_ > rb ? (dl = min(dl2, gout_ * dl)) : r_ > rcb && (dl = max(dl1, gin_ * dl))
        push!(lr, lr[end] + dl)
    end
    (length(lr) + (ENVCENTRE ? ENVLIN : 0)) % 2 == 0 && pop!(lr)                    # (an even zone count: an odd number of interfaces)
    lr[1] .+ (lr .- lr[1]) .* (log(rout) - lr[1]) / (lr[end] - lr[1])
end
# ICMLOG=m1 (default unset = unchanged): the IC mesh is logarithmic in enclosed mass instead of envlr's / the rule's (GRIDEQ) mesh:
# N = GEQN (default nsh) zones: the N interfaces above the centre at m_k = m1 (M_cloud/m1)^((k-1)/(N-1)), k = 1..N (m1 in Msun; the centre zone holds m1),
# r(m) the exact inverse of the IC's M(r) interpolation. Only for the start from the IC. The first interface is Lagrangian under GRIDEQ
# (TitanResidFast.GEQ_LAGFIRST, set from the presence of ICMLOG in the environment: restarts of the run follow).
# The start (every run without RESTART/STARFROM; formerly ICEUL=1): the initial condition directly on the fixed log mesh of
# ENVEUL/ENVCORE (lgrid = 2, regular centre), at t = 0. The mesh is `envlr`'s, from ENVCORE's first interface to the outer radius of the IC. Mass-conserving: the interior mass at
# every interface is that of the IC's tabulated M(r) (piecewise-uniform density between its nodes, so the zone masses add up to the
# cloud mass exactly); T and the H₂ fraction are the mass-weighted zone means of the IC's T(M), f_H2(M); u at the interfaces from the
# IC's v(M); x_e as the IC (XE0). Radiation: E_r = a [T_CMB⁴ + (T⁴ − T_CMB⁴)(1 − e^{−τ_R})], F_r flux-limited diffusion.
# All fixed-mesh settings (GASROW, REGULAR_CENTRE, no INNER_SINK, xmen/xmtot interior mass, CHEMFALLBACK, empty
# ALI factors). Only the abn02* initial conditions (they carry the tabulated M(r)).
"""
Zero point of the mass coordinate at the centre (fixed mesh with a regular centre; used by ICEUL, ICXMZ=0 turns it off). TITAN's
mass unknown is the exterior mass x_m,e = xmtot − m (per steradian), with xmtot ≥ the grid mass (PS.build!: grid + xmext), and the
zone masses are its differences. At the IC the central zones of the log mesh hold ~1e-17 of the cloud: zone 39 of the 100-zone mesh
(r ≈ 1.5e3 R☉) has 4.5e20 g sr⁻¹ against an ulp of 5.9e20 of x_m,e ≈ 3.2e36, so its mass row is round-off and Newton diverged at
step 1 whatever Δt (T, E_r corrections −9.4 there). Here xmtot := ε (half the central zone's mass) and x_m,e = ε − m: negative outside
the first interface, but the residual only uses differences of x_m,e, xmtot − x_m,e and relative (sign-preserving) corrections,
so the interior mass is carried at full relative precision at every radius. Recomputed from the interior masses m (summed from the
centre by PS.build!), not by shifting x_m,e.
"""
function xmzero!()
    xmen = TF.var(:xmen); xmeo = TF.var(:xmeo); dn = TF.var(:dn); dv = TF.var(:dvoln); xmt0 = TF.getscalar(:xmtot)
    TF.var(:rn)[a] == 0 || error("xmzero!: needs a regular centre (r = 0 at the first interface)")
    m = zeros(length(xmen)); for k in a+1:b+2; m[k] = m[k-1] + dn[k-1] * dv[k-1]; end      # interior mass, summed from the centre
    ε = 0.5 * m[a+1]
    for k in a:b+2; xmen[k] = ε - m[k]; end
    xmen[a-1] = 2xmen[a] - xmen[a+1]                                                   # (the regular centre's ghost, as TitanCore's bndry)
    for k in a-1:b+2; xmeo[k] = xmen[k]; end
    TF.setscalar!(:xmtot, ε); TF.setscalar!(:xmext, ε - m[b+2])
    pr("XMZERO: mass zero point at the centre, xmtot %.4e (was %.4e) g/sr; x_m,e at the interfaces %d..%d: %s; grid mass %.6f Msun\n", ε, xmt0, a - 1, a + 2,
       join((@sprintf("%.4e", xmen[k]) for k in a-1:a+2), " "), 4π * m[b+1] / Msun)
    nothing
end
function iceulbuild!(; rIset = nothing, icmlog = nothing)     # rIset: the interfaces (centre 0 … outer edge) in place of envlr's mesh (GRIDEQ's starting mesh); icmlog = (m1 [Msun], N): ICMLOG's mass mesh
    (RAD && ENVCENTRE) || error("the start needs RAD=1 and ENVEUL=0,... (regular centre, ENVCORE)")
    radsettings!()
    rout = TF.var(:rn)[b+1]; Mtot = M                                          # the IC's outer radius and cloud mass
    ctrl = Dict(s_ => TF.getscalar(s_) for s_ in (:niter, :ntry, :omega, :sigd, :sige, :stol, :conv, :thet))
    rg3 = rg .^ 3
    Mof(r) = (r <= 0 && return 0.0; r >= rg[end] && return Mg[end]; j = clamp(searchsortedlast(rg, r), 1, length(rg) - 1);
              Mg[j] + (r^3 - rg3[j]) / (rg3[j+1] - rg3[j]) * (Mg[j+1] - Mg[j]))
    # exact inverse of Mof (r^3 linear in m between the table nodes)
    rofM(m) = (j = clamp(searchsortedlast(Mg, m), 1, length(Mg) - 1); cbrt(rg3[j] + (m - Mg[j]) / (Mg[j+1] - Mg[j]) * (rg3[j+1] - rg3[j])))
    if icmlog !== nothing       # ICMLOG: the N interfaces above the centre at log-uniform enclosed mass from m1 to Mtot (N - 1 zones above the centre zone, which holds m1)
        m1_, N_ = icmlog[1] * Msun, icmlog[2]
        mk_ = [m1_ * (Mtot / m1_)^((k - 1) / (N_ - 1)) for k in 1:N_]
        rk_ = [rofM(m) for m in mk_]; rk_[end] = rout
        rk_[1] > 0 || error("ICMLOG: r(m1) <= 0")
        all(>(0), diff(rk_)) || error("ICMLOG: the interfaces are not monotone")
        for k in 1:N_-1; abs(Mof(rk_[k]) - mk_[k]) <= 1e-10 * mk_[k] || error("ICMLOG: Mof(r(m)) does not reproduce m at interface $(k+1)"); end
        rIset = vcat(0.0, rk_)
    end
    rI_ = rIset === nothing ? (r_ = exp.(envlr(rout)); r_[end] = rout; vcat([r_[1] * i / ENVLIN for i in 0:ENVLIN-1], r_)) : copy(rIset); n = length(rI_) - 1
    mI_ = [Mof(r) for r in rI_]; mI_[end] = Mtot
    all(>(0), diff(mI_)) || error("ICEUL: a zone without mass")
    if icmlog !== nothing
        dl_ = diff(log.(rI_[2:end])); i2_ = searchsortedfirst(mI_, 1e-2Msun) - 1       # (zone j spans interfaces j, j+1 of rI_; i2_ = the zone containing 1e-2 Msun)
        pr("ICMLOG: IC starting mesh, %d zones: first interface (Lagrangian) at m %.3g Msun, r %.4g R_sun, then log in enclosed mass to %.4g Msun (r %.4g R_sun), dln r %.3f .. %.3f, dln m %.4f\n",
           n, mI_[2] / Msun, rI_[2] / Rsun, Mtot / Msun, rI_[end] / Rsun, minimum(dl_), maximum(dl_), log(Mtot / mI_[2]) / (n - 1))
        pr("ICMLOG: 1e-2 Msun at r %.4g R_sun, in zone %d (of %d); interfaces 1,2,3,4 r [R_sun] %s\n", rofM(1e-2Msun) / Rsun, i2_, n,
           join((@sprintf("%.4g", rI_[k] / Rsun) for k in 2:5), " "))
    end
    ρz = [(mI_[j+1] - mI_[j]) / (4π / 3 * (rI_[j+1]^3 - rI_[j]^3)) for j in 1:n]
    # mass-weighted zone means of T(M), f_H2(M) (64 midpoints in m per zone; the IC tables are in log M [M☉], constant below their first node)
    zmean(f, j) = (ms = range(mI_[j], mI_[j+1], length = 65); sum(f(0.5 * (ms[i] + ms[i+1]) / Msun) for i in 1:64) / 64)
    fin_ = ic in ("abn02fin", "abn02red")
    Tz = fin_ ? [zmean(Tfin, j) for j in 1:n] : [zmean(m -> interp(log10(m), lM_tab, T_tab), j) for j in 1:n]
    fz = fin_ ? [zmean(m -> finy(m)[3], j) for j in 1:n] : [zmean(m -> 10^interp(log10(m), lM_tab, lf_tab), j) for j in 1:n]   # (fin: x_H2 per H)
    uI_ = fin_ ? [vfin(mI_[j] / Msun, rI_[j]) for j in 1:n+1] : ic == "abn02z19" ? [m < Mjoin ? 0.0 : 1e5 * interp(log10(m / Msun), lM_tab, v_tab) * min(1.0, (m - Mjoin) / Mjoin) for m in mI_] :
          [1e5 * interp(log10(max(m / Msun, 1e-2)), lM_tab, v_tab) for m in mI_]
    uI_[1] = 0.0; uI_[end] = 0.0                                               # regular centre; fixed outer wall (as the remap)
    xz = fin_ ? [xinit((zmean(m -> finy(m)[1], j), 0.0, fz[j], 0.0)) for j in 1:n] : [xinit((xe0, 0.0, f / 2, 0.0)) for f in fz]; eqz = falses(n)
    # initial radiation field: diffusion inside, the CMB and the free-streaming photospheric light outside
    Tb = parse(Float64, get(ENV, "TCMB", "0")); ar = Z0.arad
    χ = [PM.OPAC[](ρz[j], Tz[j])[1] * ρz[j] for j in 1:n]
    τ = zeros(n); for j in n-1:-1:1; τ[j] = τ[j+1] + 0.5 * (χ[j] * (rI_[j+1] - rI_[j]) + χ[j+1] * (rI_[j+2] - rI_[j+1])); end
    Er = [ar * (Tb^4 + (Tz[j]^4 - Tb^4) * (-expm1(-τ[j]))) for j in 1:n]
    rc = [0.5 * (rI_[j] + rI_[j+1]) for j in 1:n]; Fr = zeros(n + 1)
    for j in 2:n
        gE = (Er[j] - Er[j-1]) / (rc[j] - rc[j-1]); Em = 0.5 * (Er[j] + Er[j-1])
        Fr[j] = -Z0.c * gE / (3 * 0.5 * (χ[j] + χ[j-1]) + abs(gE) / Em) * (-expm1(-τ[j]))
    end
    Fr[n+1] = Fr[n]
    Mnew = sum(4π / 3 * ρz[j] * (rI_[j+1]^3 - rI_[j]^3) for j in 1:n)
    TX.GASROW[] = true
    for k in eachindex(mat.x); j = clamp(k - a + 1, 1, n); mat.x[k] = copy(xz[j]); mat.eqzone[k] = eqz[j]; end
    PS.build!(; rI = rI_, rho = ρz, T = Tz, u = uI_, geom = 2, rad = true, Er, Fr, lum = nothing, lgrid = 2, inner = :wall, outer = :wall, Mc = 0.0,
              selfgrav = true, material = mat, av, xmext = Mnew / 4π)
    global a = Int(TF.getscalar(:ngrs)); global b = Int(TF.getscalar(:ngre)); global nsh = b - a + 1
    get(ENV, "ICXMZ", "1") == "1" && xmzero!()
    for j in 1:nsh
        k = a + j - 1; mat.x[k] = xz[j]; mat.eqzone[k] = eqz[j]; mat.Λext[k] = 0.0; mat.Tref[k] = Tz[j]; mat.slope[k] = 0.0
    end
    for k in 1:a-1; mat.x[k] = copy(mat.x[a]); mat.eqzone[k] = mat.eqzone[a]; end
    for k in b+1:length(mat.x); mat.x[k] = copy(mat.x[b]); mat.eqzone[k] = mat.eqzone[b]; mat.Λext[k] = 0.0; end
    TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
    PS.install_julia!(BlockSolver)
    for (s_, v_) in ctrl; TF.setscalar!(s_, v_); end
    TF.setscalar!(:timen, 0.0); TF.setscalar!(:timeo, 0.0); TF.setscalar!(:jstep, 0)
    TF.setscalar!(:dtime, haskey(ENV, "ENVDT0") ? parse(Float64, ENV["ENVDT0"]) : 1e-3tff); TF.setscalar!(:tmax, 1e3tff)
    radsettings!()
    envflags!()
    TD.EOS[](); TD.OPAC[](); TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    empty!(ALIFAC[]); empty!(ALIFIX[]); SNAPS === nothing || empty!(SNAPS)
    xme = TF.var(:xmen); ENVMS0[] = xme[a] - xme[a+envsurface()-1]
    pr("ICEUL: initial condition on the fixed mesh: %d zones, r 0, %.3g R_sun .. %.4g pc (dln r %.4f), grid mass %.2f Msun (cloud %.2f), rho_c %.4e T_c %.1f, tau_R(centre) %.3e, %d interfaces falling, dt0 %.3e s, M_shell0 %.3e Msun\n",
       n, rI_[2] / Rsun, rI_[end] / 3.0857e18, log(rI_[3] / rI_[2]), Mnew / Msun, Mtot / Msun, ρz[1], Tz[1], τ[1], count(<(0), uI_), TF.getscalar(:dtime), 4π * ENVMS0[] / Msun)
end
if haskey(ENV, "STARFROM")
    starfrom!(ENV["STARFROM"])
elseif haskey(ENV, "RESTART")
    restart!(ENV["RESTART"]; allow = ENVEUL === nothing ? () : ("PSWITCH",))
    Int(TF.getscalar(:lrad)) == 1 || error("checkpoint $(ENV["RESTART"]) has no radiation field (a stage-chain checkpoint): not supported")
    if RAD
        # The checkpoint holds the converged step's complete COMMON; the EOS/opacity routines rebuild their internal tables but, run here on the
        # restored composition, give different opacities (xkpn/xken/chifn differ by up to 100 % in some zones: the in-step values come from the Newton
        # iterate's populations) and e (1e-5): a new-format checkpoint keeps its own arrays. Older checkpoints: the refresh as before.
        snap_ = JNEXT[] >= 0 ? TD.snapshot() : nothing
        radsettings!(); TD.EOS[](); TD.OPAC[]()
        snap_ === nothing || TD.restore!(snap_)
        hasproperty(CKAUX[], :edd) && edd_restore!(CKAUX[].edd)
    end
    if ENVEUL !== nothing
        Int(TF.getscalar(:lgrid)) == 3 && error("checkpoint $(ENV["RESTART"]) is on the Lagrangian mesh (stage chain): not supported")
        (envflags!(); ENVMS0[] = (ENVCENTRE && !haskey(ENV, "ENVMS0") ? 0.0 : parse(Float64, ENV["ENVMS0"])) * Msun / 4π)   # (regular centre: no sink, ENVMS0 unused)
    end
elseif haskey(ENV, "ICMLOG")      # ICMLOG=m1 [Msun]: the IC on a mesh logarithmic in enclosed mass (replaces the rule's starting mesh)
    iceulbuild!(; icmlog = (parse(Float64, ENV["ICMLOG"]), parse(Int, get(ENV, "GEQN", string(nsh)))))
else
    iceulbuild!()
    if GRIDEQP !== nothing           # GRIDEQ: build the IC again on the rule's mesh of the IC (planned on envlr's mesh), same first interface and edge
        Nq_ = parse(Int, get(ENV, "GEQN", string(nsh))); iseven(Nq_) || error("GEQN must be even")
        rq_ = geq_target(geq_state(); N = Nq_ - 1)
        pr("GRIDEQ: IC starting mesh, %d zones from %.4g R_sun, dln r %.3f .. %.3f (envlr's %.3f)\n", nsh, rq_[1] / Rsun,
           minimum(diff(log.(rq_))), maximum(diff(log.(rq_))), log(TF.var(:rn)[a+2] / TF.var(:rn)[a+1]))
        iceulbuild!(; rIset = vcat(0.0, rq_))
    end
end
if GRIDEQP !== nothing
    (ENVEUL !== nothing && ENVCENTRE && Int(TF.getscalar(:lgrid)) == 2) || error("GRIDEQ needs the fixed-mesh envelope with a regular centre (ENVEUL=0,..., lgrid = 2)")
    GRIDEQP in (0, 1) || error("GRIDEQ=0 (prescribed target) or 1 (implicit grid equation)")
    TC.GRIDMOVE[] = true; TF.setscalar!(:ctol, GEQC[8])
    TX.GEQ_LAGFIRST[] = haskey(ENV, "ICMLOG")       # ICMLOG: the first interface holds its enclosed mass (Lagrangian) instead of its radius; restarts read the same run.env
    if GRIDEQP == 0
        TX.GRIDVEL[] = zeros(length(TF.var(:rn)))
    else
        GEQHELD[] === nothing ? grideq_step!() : (TX.GRIDEQ[] = GEQHELD[])       # (a restart from a new-format checkpoint continues with the monitor the run had)
        let fmv = GEQC[5]
            fmcap = function ()        # no interface moves more than fmove × its smaller neighbouring zone per step
                rn_ = TF.var(:rn); ro_ = TF.var(:ro); f = 0.0
                for k in (TX.GEQ_LAGFIRST[] ? a+1 : a+2):b          # (ICMLOG: the Lagrangian first interface moves too)
                    wd = min(ro_[k] - ro_[k-1], ro_[k+1] - ro_[k]); f = max(f, abs(rn_[k] - ro_[k]) / wd)
                end
                f > 0 ? fmv / f * TF.getscalar(:dtime) : Inf
            end
            push!(TC.DT_HOOKS, fmcap)      # (TITAN's rule: diagnostic under DTERR)
            push!(Trbdf2.DT_BOUNDS, fmcap) # DTERR respects the mesh-motion bound
        end
    end
    pr("GRIDEQ: mode %d (%s), %d zones, weights rho,P,T,E,m,tau %s, gmax %.3g, alph %.3g, kappa %.3g, fv %.3g, fmove %.3g, tau %.3g..%.3g, ctol %.3g, kappa_f %.3g\n",
       GRIDEQP, GRIDEQP == 1 ? "implicit grid equation" : "prescribed target", nsh, join(GEQW, ","), GEQC[1:9]...)
    GEQGP[1] > 0 && pr("GEQGP (provisional): monitor term w_g %g x (dln rho/dj)^2 x gravity power of the zone / %g L_sun (implicit grid equation only)\n", GEQGP...)
end
JNEXT[] >= 0 && TF.setscalar!(:jstep, JNEXT[])      # restart: the step counter of the uninterrupted run (after the repeated grideq_step!, which ran at the checkpointed step's jstep)
const CHEMDIAG = get(ENV, "CHEMDIAG", "0") == "1"; const LNUC = Ref(0.0); const TIMING = get(ENV, "TIMING", "0") == "1"; const TSTEP = Ref(0.0); const TCHEM = Ref(0.0); const TLINE = Ref(0.0); const TOTH = Ref(0.0)
const TMARK = Ref(time()); const TOUT = Ref(0.0)
const DTSTOP = parse(Float64, get(ENV, "DTSTOP", "0"))   # stop an unattended run whose Δt has collapsed, before the checkpoints are overwritten
const PRINTEVERY = parse(Int, get(ENV, "PRINTEVERY", "500")); const CKEVERY = parse(Int, get(ENV, "CKEVERY", "500")); const KLAST = Ref(0)
# Newton controls, applied after any restart: NITER (iterations per attempt), FAILTFAC (Δt factor after a failed attempt)
haskey(ENV, "NITER") && TF.setscalar!(:niter, parse(Float64, ENV["NITER"]))
# NKRYLOV=1[,<eta>[,<mmax>[,<eps>[,<central>]]]] (provisional, default off = bit-identical): each Newton correction solved by GMRES on the true Jacobian (finite-difference products of the residual after a
# full update!), preconditioned by the band Jacobian; eta = linear tolerance relative to the band correction (1e-3), mmax = Krylov dimension (30), eps = size of the weighted perturbation (1e-5),
# central = 1 central / 0 one-sided difference. NewtonKrylov.jl, docs/POPIII_NEWTON_KRYLOV.md.
# With NWPROBE the Krylov solve is switched on at the probe state only (the steps before it run plain Newton, so the probe starts from the reference state).
# NKCHECK=<iter> (diagnostic, with NWPROBE): finite-difference J v against the band Jacobian at Newton iteration <iter> of the first attempt, then exit.
if get(ENV, "NKRYLOV", "0") != "0"
    NewtonKrylov.configure!(ENV["NKRYLOV"], activate = !haskey(ENV, "NWPROBE")); haskey(ENV, "NKCHECK") && (NewtonKrylov.OPT.check = parse(Int, ENV["NKCHECK"]))
    haskey(ENV, "NKSCAN") && (f_ = parse.(Int, split(ENV["NKSCAN"], ",")); NewtonKrylov.SCAN[] = (f_[1], f_[2:end]))      # NKSCAN=<att>,<it>[,<it>...]: diagnostic scan along the Krylov step (needs NWTRACE)
    atexit(() -> println(NewtonKrylov.summary_text()))
    pr("NKRYLOV: Newton corrections by GMRES on the true Jacobian (eta %g, Krylov dimension %d, eps %g), band Jacobian as preconditioner\n", NewtonKrylov.OPT.eta, NewtonKrylov.OPT.mmax, NewtonKrylov.OPT.eps)
end
# STAGECONS=n (diagnostic, default off = bit-identical): every n-th accepted step print `STAGECONS` lines testing stage consistency (start state with new := old: residual independent of θ, dt?; against the FSAL
# derivative; leftover residual of each converged stage; old-level vs new-level arrays), all as normalized rates [1/s] comparable with DTERRDIAG's (z2 - z1)/(y_n Δ). StageCons.jl. Run with DTERRDIAG=n DTERRTRACE=n.
haskey(ENV, "STAGECONS") && ENV["STAGECONS"] != "0" && (StageCons.configure!(parse(Int, ENV["STAGECONS"])); pr("STAGECONS: stage-consistency diagnostic every %d steps (StageCons.jl)\n", StageCons.EVERY[]))
haskey(ENV, "NWFREEZE") && (TC.NWFREEZE[] = String.(split(ENV["NWFREEZE"], ",")))     # (diagnostic, default off) see TitanCore.NWFREEZE
haskey(ENV, "NWTRACE") && (Trbdf2.NWTRACE[] = open(ENV["NWTRACE"], "a"))      # (diagnostic, default off) Newton trace, see Trbdf2.NWTRACE
haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))
haskey(ENV, "NTRY") && TF.setscalar!(:ntry, parse(Float64, ENV["NTRY"]))
if haskey(ENV, "FAILTFAC"); TD.FAIL_TFAC[] = parse(Float64, ENV["FAILTFAC"]); Trbdf2.TFAC[] = TD.FAIL_TFAC[]; end
haskey(ENV, "ERTHICK") && (TC.ER_THICK[] = parse(Float64, ENV["ERTHICK"]))
# ETOTAL=1: TITAN's conservative total-energy row (P4.2); the gas row stays where e/(e + KE) < ETASW, in the cold hypersonic infall, where T is a small difference
haskey(ENV, "ETASWTH") && (TX.ETASW_TH[] = parse(Float64, ENV["ETASWTH"]))     # internal-energy row also where 3p/2ρ < ETASWTH × (3p/2ρ + E_r/ρ + K/m) (off: as before)
haskey(ENV, "ETASWTIE") && (TX.ETASW_TIE[] = parse(Float64, ENV["ETASWTIE"]))   # internal-energy row also where the FLUXTIE weight at the zone's inner interface >= ETASWTIE (off: as before)
get(ENV, "ETOTAL", "0") == "1" && (TX.ETOTAL[] = true; TX.ETOT_GRAV[] = true; TX.ETASW_EOFF[] = PC.UOFF_TH; TX.ETASW[] = parse(Float64, get(ENV, "ETASW", "0.1")))
haskey(ENV, "EXCAP") && (TX.EXCAP[] = parse(Float64, ENV["EXCAP"]))
# SUBZONE=1[,τ_f] (default off): interface optical depth across an unresolved front from the two-layer sub-zone form (TX.SUBZONE)
if get(ENV, "SUBZONE", "0") != "0"
    let w_ = parse.(Float64, split(ENV["SUBZONE"], ","))
        TX.SUBZONE[] = w_[1] != 0; length(w_) > 1 && (TX.SZTAUF[] = w_[2])
        TX.SUBZONE[] && pr("SUBZONE: interface optical depth across unresolved fronts from the sub-zone form (opaque half-zone to tau_f = %.3g, transparent half-zone in series; front weight s(a)(1-s(b)), s = x^4/(x^4+tau_f^4))\n", TX.SZTAUF[])
    end
end
# SLABFLUX=1[,n] (default off): interface opacity from the exact steady slab solution between the zone centres, opacity profile
# κ ∝ (τ + τ₀)^n (n omitted or 0: n → ∞, exponential in optical depth) instead of the harmonic mean (TX.SLABFLUX, docs/POPIII_SLABFLUX.md)
if get(ENV, "SLABFLUX", "0") != "0"
    let w_ = parse.(Float64, split(ENV["SLABFLUX"], ","))
        TX.SLABFLUX[] = w_[1] != 0; length(w_) > 1 && w_[2] != 0 && (TX.SFN[] = w_[2])
        TX.SFN[] >= 1 || error("SLABFLUX=1,n needs n >= 1 (or 0 for n = Inf)")
        TX.SLABFLUX[] && pr("SLABFLUX: interface opacity from the steady slab solution, profile kappa ~ (tau + tau0)^n, n = %s%s\n", isinf(TX.SFN[]) ? "Inf (exponential in tau)" : @sprintf("%.3g", TX.SFN[]), TX.SUBZONE[] ? "; SUBZONE applied on top" : "")
    end
end
# SLABTIE=1 (default off; needs SLABFLUX=1 with n = Inf): FLUXTIE's outward donor h from the slab depth of the donor centre below the interface (TX.SLABTIE)
if get(ENV, "SLABTIE", "0") == "1"
    (TX.SLABFLUX[] && isinf(TX.SFN[])) || error("SLABTIE=1 needs SLABFLUX=1 (n = Inf)")
    TX.SLABTIE[] = true
    pr("SLABTIE: FLUXTIE outward donor h = h_FS + f/(tau_i + 2f) - f/(a + 2f), tau_i the slab depth of the donor centre below the interface, a its half-zone depth\n")
    if haskey(ENV, "TIESMOOTH") && ENV["TIESMOOTH"] != "0"      # (default off) Newton fix 1 of docs/POPIII_MULTIGROUP.md §5
        TX.TIESM[] = parse(Float64, ENV["TIESMOOTH"])
        pr("TIESMOOTH: SLABTIE's max(h, 0) replaced by the smooth ramp (h + sqrt(h^2 + eps^2))/2, eps = %g\n", TX.TIESM[])
    end
end
# TIEFRAME=1 (provisional, default off): FLUXTIE and the outer boundary row in the frame of the unknowns (comoving), TX.TIEFRAME, docs/POPIII_ENERGY.md §3
if get(ENV, "TIEFRAME", "0") == "1"
    TX.TIEFRAME[] = true
    pr("TIEFRAME: FLUXTIE closes the lab-frame flux, F_0 = [h c - (1 + f - 2 h^2) u] E_0 for the comoving unknowns; outer boundary row likewise\n")
end
# ETOTAL=0 (or ECOMPAT=1): compatible internal-energy rows in every zone, no row switching and no total-energy row (TX.ECOMPAT, docs/POPIII_ENERGY.md §4, §23); implies TIEFRAME
if get(ENV, "ECOMPAT", "0") == "1" || get(ENV, "ETOTAL", "") == "0"
    TX.GASROW[] || error("ETOTAL=0 / ECOMPAT=1 needs the gas-alone row (GASROW)")
    TX.ECOMPAT[] = true
    TX.TIEFRAME[] || (TX.TIEFRAME[] = true; pr("ECOMPAT: TIEFRAME switched on\n"))
    pr("ECOMPAT: first-law gas row adjoint to the momentum row in every zone (ETOTAL/ETASW* row switching unused), lab-frame radiation energy row\n")
end
# ERLAB=1 (provisional, default off): ECOMPAT's lab-frame radiation energy row with the stack's gas rows (TX.ERLAB, docs/POPIII_ENERGY.md §18); implies TIEFRAME
if get(ENV, "ERLAB", "0") == "1"
    TX.GASROW[] || error("ERLAB=1 needs the gas-alone row (GASROW)")
    TX.ECOMPAT[] && error("ERLAB=1: ECOMPAT=1 already has the lab-frame radiation row")
    TX.ERLAB[] = true
    TX.TIEFRAME[] || (TX.TIEFRAME[] = true; pr("ERLAB: TIEFRAME switched on\n"))
    pr("ERLAB: lab-frame radiation energy row; gas rows of the stack (dual-energy)\n")
end
# GFV=1 (provisional, default off): finite-volume gravity weight in the gas momentum row (TX.GFV, docs/POPIII_ENERGY.md §19)
if get(ENV, "GFV", "0") == "1"
    TX.GFV[] = true
    pr("GFV: gas momentum row with the finite-volume gravity weight 4 pi G M_j (rho_j-1 dr_j-1 + rho_j dr_j)/2\n")
end
# GVAR=1 (provisional, default off): gravity force of the gas momentum row = gradient of the exact binding energy of the piecewise-constant density (TX.GVAR, docs/POPIII_ENERGY.md §19)
if get(ENV, "GVAR", "0") == "1"
    TX.GFV[] && error("GVAR=1 and GFV=1 are alternatives")
    TX.GVAR[] = true
    pr("GVAR: gas momentum row with the gravity force dW/dr_j of the piecewise-constant density\n")
end

# NWFD (diagnostic): for the unknowns j of every zone z in zlo:zhi, the central difference (h = 1e-5) of (a) the residual with the Newton context frozen (what
# JACCHK does: everything the iterate's update! refreshes outside residual! -- EOS/opacity fills, Eddington factors f, g, h, the flux-tie weights -- held),
# and (b) the full nonlinear map: the state is perturbed with update! (EOS, opacity, formal solution refreshed), a fresh context is built and the residual
# evaluated at x = 0; the state is restored by the inverse update. Both against the assembled Jacobian of the iterate. Prints every (row, zone, var) with
# |FD - AD| > 1e-3 x (largest diagonal entry of the row) for rows in zones z-3:z+3, for (a) and (b) separately, and the (b)-(a) difference.
function nwfd!(zlo, zhi, iter)
    a_ = Int(TF.getscalar(:ngrs)); c = TX.context(); ne = c.neqn; nb = c.nb; k0 = c.k0; ws = TX.getws(Float64, c.mgr)
    E = (TF.var(:em2), TF.var(:em1), TF.var(:e00), TF.var(:ep1), TF.var(:ep2)); AD = [copy(e) for e in E]
    rhs = TF.var(:rhs); h = parse(Float64, get(ENV, "NWFDH", "1e-5"))
    pr("NWFD: step dt %.4e theta %.2f iteration %d, neqn %d, zones %d..%d (relative numbering; absolute = +%d), code %s\n", TF.getscalar(:dtime), TF.getscalar(:thet), iter, ne, zlo, zhi, a_ - 1, string(c.code))
    base_ = (; f = copy(TF.var(:fedd)), g = TF.getscalar(:geddr), h = copy(TD.HEDD), hx = copy(TD.HXEDD), fx = copy(TD.FXEDD), gx = TD.GXEDD[], tie = copy(TD.TIEW))
    FZ = Ref(String[])
    function resid_now()
        for n_ in FZ[]
            n_ == "f" && copyto!(TF.var(:fedd), base_.f); n_ == "g" && TF.setscalar!(:geddr, base_.g)
            n_ == "h" && copyto!(TD.HEDD, base_.h); n_ == "hx" && copyto!(TD.HXEDD, base_.hx); n_ == "fx" && copyto!(TD.FXEDD, base_.fx)
            n_ == "gx" && (TD.GXEDD[] = base_.gx); n_ == "tie" && copyto!(TD.TIEW, base_.tie)
        end
        cc = TX.context(); Rr = zeros(ne, nb); TX.residual!(Rr, ws, cc, zeros(ne * nb)); Rr
    end
    function bump!(za, j, δ)
        fill!(rhs, 0.0); rhs[j, za] = δ; TC.update!()
    end
    revert(δ, j) = (j <= 3 || j == 5 || j == 6) ? -δ / (1 + δ) : -δ
    Rf = zeros(ne, nb)
    R0 = resid_now()
    sets_ = [s_ == "none" ? String[] : String.(split(s_, ",")) for s_ in split(get(ENV, "NWFDSETS", "none"), ";")]
    for (iset, fz_) in enumerate(sets_)
    FZ[] = fz_; ACC = Dict{Tuple{Int,Int},Tuple{Float64,Int,Float64}}()
    pr("NWFD set %d: closure quantities held at the pre-perturbation values in the full FD: %s\n", iset, isempty(fz_) ? "none (full refresh)" : join(fz_, ","))
    for z in zlo:zhi
        za = a_ + z - 1; kbz = za - k0 + 1
        for j in 1:ne
            xv = zeros(ne * nb); ix = (kbz - 1) * ne + j
            Rp = zeros(ne, nb); Rm = zeros(ne, nb)
            xv[ix] = h; TX.residual!(Rp, ws, c, xv); xv[ix] = -h; TX.residual!(Rm, ws, c, xv)
            FDa = (Rp .- Rm) ./ 2h
            bump!(za, j, h); Rp2 = resid_now(); bump!(za, j, revert(h, j))
            bump!(za, j, -h); Rm2 = resid_now(); bump!(za, j, revert(-h, j))
            FDb = (Rp2 .- Rm2) ./ 2h
            for kb in max(1, kbz - 3):min(nb, kbz + 3), i in 1:ne
                o = kbz - kb; ad = abs(o) <= 2 ? AD[o + 3][i, j, k0 + kb - 1] : 0.0
                scale = maximum(abs, (AD[3][i, jj, k0 + kb - 1] for jj in 1:ne)) + 1e-300
                ea = abs(FDa[i, kb] - ad) > 1e-3 * scale; eb = abs(FDb[i, kb] - ad) > 1e-3 * scale; ec = abs(FDb[i, kb] - FDa[i, kb]) > 1e-3 * scale
                let v_ = abs(FDb[i, kb] - ad) / scale; o_ = get(ACC, (i, j), (0.0, 0, 0.0)); ACC[(i, j)] = (max(o_[1], v_), o_[2] + (v_ > 1e-3 ? 1 : 0), max(o_[3], abs(FDa[i, kb] - ad) / scale)); end
                (iset == 1 && (ea || eb || ec)) && pr("   row %d zone %d | col var %d zone %d (offset %d): FDfrozen %.5e FDfull %.5e AD %.5e diagscale %.3e | %s%s%s\n", i, kb + k0 - a_, j, z, o, FDa[i, kb], FDb[i, kb], ad, scale,
                                       ea ? "frozen-vs-AD " : "", eb ? "full-vs-AD " : "", ec ? "full-vs-frozen" : "")
            end
        end
    end
    for ((i_, j_), (mx_, n_, mf_)) in sort(collect(ACC)); (mx_ > 1e-3 || mf_ > 1e-3) && pr("   SET%d row %d col var %d: max |FDfull-AD|/diag %.4g (%d entries > 1e-3), max |FDfrozen-AD|/diag %.4g\n", iset, i_, j_, mx_, n_, mf_); end
    end
    pr("NWFD: done; state drift check, max |R(after) - R(before)| / (max |R|) = %.3e\n", maximum(abs, resid_now() .- R0) / (maximum(abs, R0) + 1e-300))
    nothing
end
const LINECAPV = Ref(NaN)
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
TC.EDD_SCATTER[] = get(ENV, "EDDSCAT", "0") == "1"
if haskey(ENV, "FLUXTIE")
    let v = parse.(Float64, split(ENV["FLUXTIE"], ","))
        TD.FLUXTIE[] = (v[1], v[2])
        pr("FLUXTIE: radiation momentum row -> F = h c E (h = H/J of the formal solution) where tau from the outer boundary <= %g, blended to the moment equation by tau = %g\n", v[1], v[2])
    end
end
TC.EDDG_MU0[] = get(ENV, "EDDGMU0", "0") == "1"
# VEFIMP=1 (G4, default off): consistent closure of the gray moments. The formal solution's f, g and h (h with FLUXTIE) at every Newton
# iterate (TitanDriver.VEF_ITER), the momentum row and outer boundary in Auer's sphericality-factor form (VEF_SPH), the CMB only as the
# isotropic background, added to the moments analytically (VEF_BG); FLUXTIE, if set, then ties F = h c E with the h of the iterate.
# VEFIMP=notie: the same with FLUXTIE (and ETASWTIE) off, the moment equations everywhere; nobg: notie without VEF_BG; it / sph: one part (tests).
if get(ENV, "VEFIMP", "0") != "0"
    let v = ENV["VEFIMP"]
        v in ("1", "notie", "nobg", "it", "sph") || error("VEFIMP must be 1, notie, nobg, it or sph")
        TD.VEF_ITER[] = v != "sph"; TD.VEF_SPH[] = v != "it"; TD.VEF_BG[] = v in ("1", "notie")
        if v != "1" && (TD.FLUXTIE[][2] > 0 || TX.ETASW_TIE[] > 0)
            pr("VEFIMP: FLUXTIE (%g,%g) and ETASWTIE %g switched off\n", TD.FLUXTIE[]..., TX.ETASW_TIE[])
            TD.FLUXTIE[] = (0.0, 0.0); TX.ETASW_TIE[] = 0.0
        end
        pr("VEFIMP=%s: Eddington factors f, g%s %s; radiation momentum row and outer boundary %s%s\n", v, TD.FLUXTIE[][2] > 0 ? " and the flux tie's h" : "",
           TD.VEF_ITER[] ? "of the formal solution at every Newton iterate (not lagged)" : "lagged (formal solution of the predictor)",
           TD.VEF_SPH[] ? "in Auer's sphericality-factor form (exact for chi = 0)" : "as TITAN",
           TD.VEF_BG[] ? "; the CMB only as the isotropic background, added to the formal solution's moments analytically" : "")
    end
end
# VEFBGCHK=n (diagnostic, default 0): at the first n formal solutions compare f, g (and h) of the background split (VEF_BG) with
# those of the formal solution that carries the CMB along the rays; the two agree to the quadrature error where nothing is clamped
if parse(Int, get(ENV, "VEFBGCHK", "0")) > 0 && TD.VEF_BG[]
    let ef0 = TC.EDDFAC[], nchk = Ref(parse(Int, ENV["VEFBGCHK"]))
        TC.EDDFAC[] = function ()
            nchk[] <= 0 && return ef0()
            nchk[] -= 1
            a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); fe = TF.var(:fedd)
            TD.VEF_BG[] = false; ef0(); f0 = fe[a_:b_]; g0 = TF.getscalar(:geddr); h0 = length(TD.HEDD) >= b_ ? TD.HEDD[a_:b_] : Float64[]
            TD.VEF_BG[] = true; ef0(); f1 = fe[a_:b_]; g1 = TF.getscalar(:geddr); h1 = length(TD.HEDD) >= b_ ? TD.HEDD[a_:b_] : Float64[]
            d = abs.(f1 .- f0); j = argmax(d)
            pr("VEFBGCHK: g %.8f / %.8f; max |df| %.3e at zone %d (f %.6f / %.6f, f_x %.4f); |df| > 1e-3 in zones %s%s\n", g0, g1, d[j], j, f0[j], f1[j],
               TD.FXEDD[a_+j-1], join(findall(>(1e-3), d), ","), isempty(h0) ? "" : @sprintf("; max |dh| %.3e", maximum(abs.(h1 .- h0))))
        end
    end
end
# TIEX=1 (provisional, default off): FLUXTIE on the excess over the isotropic background (TD.TIEX, docs/POPIII_ENERGY.md §15)
if get(ENV, "TIEX", "0") == "1"
    (TD.VEF_BG[] && TD.FLUXTIE[][2] > 0) || error("TIEX=1 needs VEFIMP=1 (VEF_BG) and FLUXTIE")
    TD.TIEX[] = true; pr("TIEX: FLUXTIE closes F = h_x c (E - E_b) with the flux factor h_x = H_x/J_x of the excess over the isotropic background\n")
end
# ERESTOL=<L_sun> (provisional, default 0 = off): residual-based stopping test of the stage Newton iterations in addition to the correction
# test (Trbdf2.ERESTOL, docs/POPIII_ENERGY.md §17)
if (v = parse(Float64, get(ENV, "ERESTOL", "0"))) > 0
    get(ENV, "INTEGRATOR", "trbdf2") == "trbdf2" || error("ERESTOL needs INTEGRATOR=trbdf2")
    Trbdf2.ERESTOL[] = v; pr("ERESTOL: stage Newton iterations continue until |sum of the energy rows + u x momentum rows| <= %g L_sun (or it stalls)\n", v)
end
const RUNEND = Ref{Any}(nothing)
"""
TABAUDIT=1 (diagnostic): for every table zone below the hot branch (and, with TABAUDIT=2, also every network zone) compare at the zone's (ρ, T)
the equilibrium table (PopIIIEqEOS composition, u, p; StellarOpacity κ_R, κ_P) with the network's own LTE state (ChemNet.steady with J = B(T),
which by detailed balance is also the thick-limit three-level state): x_e, x_H2, x_H-, He stages, u, p, c_v and Γ1 (derivatives at
equilibrium composition, finite differences of the steady state), and the continuum means of PopIIIMaterial.nlte_means for the network
composition and for the table composition. Prints one line per zone, then exits.
"""
function tabaudit()
    B(T) = ν -> 2 * ChemNet.hP * ν^3 / ChemNet.c^2 / expm1(ChemNet.hP * ν / (ChemNet.kB * T))
    net = mat.net; ph = filter(p -> true, mat.photos); dn_ = TF.var(:dn); tn_ = TF.var(:tn); rn_ = TF.var(:rn); h = PM.HOT[]
    lte(ρ, T, x0) = ChemNet.steady(x0, ρ, T, net; env = (Trad = T, J = B(T), static = true), photos = ph)
    E21 = PM.E21; kB = ChemNet.kB
    function netstate(ρ, T, x0)
        x = lte(ρ, T, x0); u, p = ChemNet.eos(ρ, T, x, net); x, u, p
    end
    function means(ρ, T, d)
        nH = ρ / ChemNet.μH(); n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
        nHe = nH * max(get(d, :He, PC.yHe), 0.0); nHmS = n1 * ne / ChemNet.Keq([:Hm], [:H, :e], T); nHm = nH * max(get(d, :Hm, 0.0), 0.0)
        PM.nlte_means(ρ, T; n1, n2 = 4n1 * exp(-E21 / (kB * T)), nH2 = nH * max(d[:H2], 0.0), nHe, nHp = np, ne, nHm, nHmS, hbfP = true)
    end
    @printf("TABAUDIT: zone r[Rsun] T[K] nH | x_e tab/net | x_H2 tab/net | x_H- tab/net | x_He+ x_He++ tab/net | u tab/net (rel diff) | p rel diff | c_v rel diff | Γ1 tab/net | κR_tab κP_tab | κ(E,P,R) nlte net comp | κ(E,P,R) nlte tab comp\n")
    if get(ENV, "TABAUDIT", "1") == "4"          # opacity audit of every zone below the hot branch (table and network zones)
        X = 1 / (1 + 4 * PC.yHe); dlr(ρ, T) = log10(ρ) - 3 * (log10(T) - 6)
        @printf("TABAUD4: zone r[Rsun] T nH eq | logT logR flag(1 OPAL 2 ext 3 blend 4 own) | OPAC κR κP | run κR κE κP | nlte(net comp) κE κP κR | nlte(LTE comp) κE κP κR | x_e x_H2 x_H- net\n")
        for k in a:b
            T = tn_[k]; (h !== nothing && T > h.T1) && continue
            ρ = dn_[k]; nH = ρ / ChemNet.μH(); lT = log10(T); lR = dlr(ρ, T)
            try
                fl = StellarOpacity.coverage(lT, lR, X); κR, κP = PM.OPAC[](ρ, T)
                TX.ZONE[] = k; rr = TX.material_opac(mat, ρ, T)
                dn2 = ChemNet.full(mat.x[k], net); dt_ = ChemNet.full(PM.x_from_table(mat, ρ, T), net)
                mn = means(ρ, T, dn2); mt = means(ρ, T, dt_); g(d, s_) = Float64(get(d, s_, 0.0))
                @printf("TABAUD4 %3d %9.4g %8.0f %.2e %d | %.3f %.2f %d | %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e | %.3e %.2e %.2e\n",
                        k - a + 1, rn_[k] / 6.957e10, T, nH, mat.eqzone[k], lT, lR, fl, κR, κP, rr..., mn..., mt..., g(dn2, :e), g(dn2, :H2), g(dn2, :Hm))
            catch e
                println("TABAUD4 $(k - a + 1) failed: ", first(sprint(showerror, e), 120))
            end
        end
        return nothing
    end
    # TABAUDIT=3: a (n_H, T) grid instead of the zones (where the ideal network departs from the table's non-ideal SCVH)
    pts = get(ENV, "TABAUDIT", "1") == "3" ? [(0, nH_ * ChemNet.μH(), T_) for nH_ in 10.0 .^ (parse(Float64, get(ENV, "TABAUDNH", "17")):0.5:24) for T_ in (6e3, 8e3, 1e4, 1.5e4, 2e4, 3e4, 5e4, 1e5, 2e5, 2.9e5)] :
          [(k, dn_[k], tn_[k]) for k in a:b if (mat.eqzone[k] || get(ENV, "TABAUDIT", "1") == "2") && !(h !== nothing && tn_[k] > h.T1)]
    for (k, ρ, T) in pts
        nH = ρ / ChemNet.μH()
        try
            xt = PM.x_from_table(mat, ρ, T); dt_ = ChemNet.full(xt, net)
            ut, pt = PopIIIEqEOS.eos(ρ, T)
            fT(t) = PopIIIEqEOS.eos(ρ, t); fR(r) = PopIIIEqEOS.eos(r, T)
            uT = ForwardDiff.derivative(t -> fT(t)[1], T); pT = ForwardDiff.derivative(t -> fT(t)[2], T)
            uR = ForwardDiff.derivative(r -> fR(r)[1], ρ); pR = ForwardDiff.derivative(r -> fR(r)[2], ρ)
            Γt = ρ / pt * (pR + pT * (pt / ρ^2 - uR) / uT)
            xn, un, pn = netstate(ρ, T, xt); dn2 = ChemNet.full(xn, net)
            δ = 1e-3
            _, u1, p1 = netstate(ρ, T * (1 + δ), xn); _, u0, p0 = netstate(ρ, T * (1 - δ), xn)
            _, u3, p3 = netstate(ρ * (1 + δ), T, xn); _, u2, p2 = netstate(ρ * (1 - δ), T, xn)
            nuT = (u1 - u0) / (2δ * T); npT = (p1 - p0) / (2δ * T); nuR = (u3 - u2) / (2δ * ρ); npR = (p3 - p2) / (2δ * ρ)
            Γn = ρ / pn * (npR + npT * (pn / ρ^2 - nuR) / nuT)
            κR, κP = PM.OPAC[](ρ, T); mn = means(ρ, T, dn2); mt = means(ρ, T, dt_)
            g(d, s) = Float64(get(d, s, 0.0))
            @printf("TABAUDIT %3d %8.3g %8.0f %.2e | %.3e %.3e | %.2e %.2e | %.2e %.2e | %.2e %.2e  %.2e %.2e | %.4e %.4e (%+.2e) | %+.2e | %+.2e | %.4f %.4f | %.3e %.3e | %.3e %.3e %.3e | %.3e %.3e %.3e\n",
                    k == 0 ? 0 : k - a + 1, k == 0 ? ρ : rn_[k] / 6.957e10, T, nH, g(dt_, :e), g(dn2, :e), g(dt_, :H2), g(dn2, :H2), g(dt_, :Hm), g(dn2, :Hm),
                    g(dt_, :Hep), g(dt_, :Hepp), g(dn2, :Hep), g(dn2, :Hepp), ut, un, (un - ut) / abs(ut), (pn - pt) / pt, (nuT - uT) / uT,
                    Γt, Γn, κR, κP, mn..., mt...)
        catch e
            println("TABAUDIT $(k == 0 ? 0 : k - a + 1) T $T nH $nH: failed: ", first(sprint(showerror, e), 120))
        end
    end
end
"""
CPLPROF=1 (diagnostic, default off; also without COUPLED): wall time of matgen (residual + AD Jacobian), matslv, maxdel, update and the
Newton iteration count, and (COUPLED) the local 4-wide Jacobians of PopIIIMaterial.cpl_zone (computed / reused, time); printed every 100 steps.
"""
const CPLPROF = get(ENV, "CPLPROF", "0") == "1"
const PROFT = Dict{Symbol,Float64}(); const PROFN = Dict{Symbol,Int}(); const PROFSTEP = Ref(-1)
const PROFB = Dict{Symbol,Float64}()            # (PERF) bytes allocated inside each routine (all threads, Base.gc_bytes)
const CPLPROFEVERY = parse(Int, get(ENV, "CPLPROFEVERY", "100"))   # (PERF) CPLPROF print interval in steps
if CPLPROF
    for s_ in (:matgen, :matslv, :maxdel, :update, :peq, :timstp)
        let f0 = TD.ROUTINES[s_], s_ = s_
            TD.ROUTINES[s_] = () -> (t0 = time(); b0 = Base.gc_bytes(); try f0() finally PROFT[s_] = get(PROFT, s_, 0.0) + time() - t0; PROFN[s_] = get(PROFN, s_, 0) + 1
                                     PROFB[s_] = get(PROFB, s_, 0.0) + (Base.gc_bytes() - b0) end)
        end
    end
end
"""
MGROUPS (docs/POPIII_MULTIGROUP.md): G radiation groups as Newton unknowns 10 … 9 + 2G (TitanResidFast MGRad.jl), group opacity ratios
before every hydro step (PopIIIMaterial.mg_ratios!), the star-light rates from the group fields (starlight_setup), the pivoted banded
solver. The group fields come from the checkpoint (MGAUX) or, the first time, from the grey field (CMB part at T_cmb, the excess at T_r
inside the photosphere and at the star's T_eff outside it).
"""
function mg_install!()
    COUPLED || error("MGROUPS needs COUPLED=1")
    Int(TF.getscalar(:lrad)) == 1 || error("MGROUPS needs the radiation field (RAD=1)")
    (TX.COMPTON[] || TX.INNER_SINK[] || !isnan(TX.PEXT[])) && error("MGROUPS: COMPTON, INNER_SINK, PEXT are not supported")
    edges = PM.mg_edges(ENV["MGROUPS"]); G = length(edges) + 1
    TF.PARAMS[:meqn] >= 9 + 2G || error("MGROUPS: meqn = $(TF.PARAMS[:meqn]) < $(9 + 2G)")
    tcmb = parse(Float64, get(ENV, "TCMB", "0"))
    TX.mg_setup!(edges, PM.EVERG / ChemNet.kB, tcmb, 9, length(TF.var(:dn)))
    TF.setscalar!(:neqn, 9 + 2G)
    PM.mg_roles!(); r_ = PM.MGROLE[]; TX.MGQG[] = r_.Lya; TX.MGSLG[] = r_.Hm
    PM.SLIMP && (TX.MGSLQ[] = PM.sl_qt_mg)
    get(ENV, "MGSOLVER", "band") == "band" && (TD.ROUTINES[:matslv] = JuliaSolver.matslv!)
    if get(ENV, "MGVEF", "1") == "1"               # stage 2a: Eddington factors per group (mg_eddfac!)
        TX.MGVEF[] = true
        let f0 = TC.EDDFAC[]; TC.EDDFAC[] = () -> mg_eddfac!(f0) end
    end
    h0 = TC.CPLH[]
    TC.CPLH[] = (; newstep = () -> (h0.newstep(); TX.mg_newstep!()), reset = () -> (h0.reset(); TX.mg_reset!()),
                 update = (rhs, a_, b_) -> (h0.update(rhs, a_, b_); TX.mg_update!(rhs, a_, b_)),
                 dmax = (rhs, a_, b_) -> max(h0.dmax(rhs, a_, b_), TX.mg_dmax!(rhs, a_, b_)))
    let f0 = TD.STEPFN[]
        TD.STEPFN[] = () -> (PM.mg_ratios!(mat); f0(); nothing)
    end
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    ax = MGAUX[]
    if ax !== nothing && ax.edges == TX.MGEDGE[]
        for g in 1:G; copyto!(TX.MGEN[][g], ax.E[g]); copyto!(TX.MGFN[][g], ax.F[g]); copyto!(TX.MGEO[][g], ax.E[g]); copyto!(TX.MGFO[][g], ax.F[g]); end
        TX.mg_newstep!(); src = "the checkpoint"
    else
        dn, rn, chi, er, fr = TF.var(:dn), TF.var(:rn), TF.var(:chifn), TF.var(:ern), TF.var(:frn)
        arad = 7.565723e-15; τ = zeros(length(dn)); t_ = 0.0
        for k in b_:-1:a_; t_ += chi[k] * dn[k] * (rn[k+1] - rn[k]); τ[k] = t_; end
        kph = a_; for k in b_:-1:a_; τ[k] >= 2 / 3 && (kph = k; break); end
        L = 4π * rn[kph+1]^2 * fr[kph+1]; Te = L > 0 ? (L / (4π * rn[kph+1]^2 * 5.670374419e-5))^0.25 : 0.0
        # the excess over the CMB has the star's colour only where it is the star's light (at least half the free-streaming estimate
        # L/(4π r² c)); farther out it is the cloud's own emission: T_s = the gas temperature there. (The T_eff shape everywhere put
        # 0.1–0.3 E_CMB of optical light into the 1–143 pc cloud, which per-group free streaming then drained at up to 5e8 L☉, t_mg28.)
        # the excess over the CMB has the star's colour only where it is the star's light (at least half the free-streaming estimate
        # L/(4π r² c)); farther out it is the cloud's own emission: T_s = the gas temperature there. (The T_eff shape everywhere put
        # 0.1–0.3 E_CMB of optical light into the 1–143 pc cloud, which per-group free streaming then drained at up to 5e8 L☉, t_mg28.)
        Ec_ = arad * parse(Float64, get(ENV, "TCMB", "0"))^4; tg = TF.var(:tn); rr = TF.var(:rn)
        lit(k) = Te > 0 && k > kph && k < length(rr) && (er[k] - Ec_) >= 0.5 * L / (4π * (0.5 * (rr[k] + rr[k+1]))^2 * 2.99792458e10)
        kend = kph + 1; while kend <= b_ && lit(kend); kend += 1; end          # one star-lit region from the photosphere out (no zone-to-zone switching)
        Ts = [(τ[k] >= 2 / 3 || Te <= 0) ? (max(er[k], 1e-300) / arad)^0.25 : (k > kph && k < kend) ? max(Te, (max(er[k], 1e-300) / arad)^0.25) : max(tg[k], 1.0)
              for k in eachindex(dn)]
        TX.mg_init_from_grey!(Ts); src = @sprintf("the grey field (T_eff %.0f K outside zone %d)", Te, kph - a_ + 1)
        # (Tried: the excess beyond the star-lit region at T_CMB, and a continuous starlight/CMB split capped at L/(4π r² c), e8d0b3b: both fail
        #  the first step (t_mg40–t_mg43); this partition leaves a draining optical reservoir in the far cloud, docs §14.)
    end
    PM.mg_ratios!(mat)
    pr("MGROUPS: %d groups, edges %s eV, unknowns 10..%d (neqn %d, meqn %d), q from group %d, chemistry energy from group %d, solver %s; fields from %s\n",
       G, join((@sprintf("%.4g", e) for e in edges), ","), 9 + 2G, 9 + 2G, TF.PARAMS[:meqn], TX.MGSLG[], TX.MGQG[], get(ENV, "MGSOLVER", "band"), src)
    nothing
end
"""
(MGROUPS, MGVEF, stage 2a) Eddington factors per group: after the grey formal solution `f0`, the same formal solution (TitanCore.eddfac!,
tangent rays, linear-in-τ source, VEF_BG background, EDDSCAT) once per group on the group's arrays, swapped into the state arrays it reads:
χ → χ_g = max(χ r_R,g, κ_E,g), absorption κ_P → κ_E,g = κ_E r_E,g + κ_1s,g, B → κ_P,g B_g/κ_E,g (so that the source is
(κ_P,g B_g + (χ_g − κ_E,g) J_g)/χ_g: emission κ_P,g B_g, absorption κ_E,g, the rest of χ_g scattering), E_r → E_g, F_r → F_g, I_cmb → I_cmb b_g(T_cmb).
Out: f_g, h_g, the tie weights of the group's own depth, g_g. Then the grey state is restored and fedd := Σ f_g E_g / Σ E_g.
"""
function mg_eddfac!(f0)
    f0()
    TX.MGVEF[] || return nothing
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); G = TX.MGN[]
    ins = (:chifn, :xkpn, :plfn, :ern, :frn, :dtau); sv = Dict(s_ => copy(TF.var(s_)) for s_ in (ins..., :fedd))
    hsv = copy(TD.HEDD); wsv = copy(TD.TIEW); fxsv = copy(TD.FXEDD); gxsv = TD.GXEDD[]
    gr0, gl0, xi0 = TF.getscalar(:geddr), TF.getscalar(:geddl), TF.getscalar(:ximr)
    chi, kp, pl, er, fr = (TF.var(s_) for s_ in (:chifn, :xkpn, :plfn, :ern, :frn)); xke = TF.var(:xken); tn = TF.var(:tn)
    try
        for g in 1:G
            RE, RP, RR, K1 = TX.MGRE[][g], TX.MGRP[][g], TX.MGRR[][g], TX.MGK1[][g]
            for k in a_-1:b_+1
                kE = max(xke[k] * RE[k] + K1[k], 1e-300)
                kP = sv[:xkpn][k] * RP[k] + TX.MGK1E[][g][k]; B = sv[:plfn][k] * TX.mg_bfrac(g, tn[k])
                chi[k] = max(sv[:chifn][k] * RR[k], kE); kp[k] = kE; pl[k] = kP * B / kE
            end
            copyto!(er, TX.MGEN[][g]); copyto!(fr, TX.MGFN[][g]); TF.setscalar!(:ximr, xi0 * TX.MGBCMB_[g])
            TC.HPON[] = true
            try f0() finally TC.HPON[] = false end
            hp = TX.MGHPLUS[][g]; @inbounds for k in eachindex(hp); hp[k] = k <= length(TC.HPLUS) ? TC.HPLUS[k] : 0.0; end
            copyto!(TX.MGFEDD[][g], TF.var(:fedd))
            length(TD.HEDD) >= b_ + 1 && (h = TX.MGHEDD[][g]; @inbounds for k in eachindex(h); h[k] = k <= length(TD.HEDD) ? TD.HEDD[k] : 0.0; end)
            length(TD.TIEW) >= b_ + 1 && (w = TX.MGTIEW[][g]; @inbounds for k in eachindex(w); w[k] = k <= length(TD.TIEW) ? TD.TIEW[k] : 0.0; end)
            TX.MGGEDDR[][g] = TF.getscalar(:geddr)
        end
    finally
        for s_ in (ins..., :fedd); copyto!(TF.var(s_), sv[s_]); end
        copyto!(TD.HEDD, hsv); copyto!(TD.TIEW, wsv); copyto!(TD.FXEDD, fxsv); TD.GXEDD[] = gxsv
        TF.setscalar!(:geddr, gr0); TF.setscalar!(:geddl, gl0); TF.setscalar!(:ximr, xi0)
    end
    fe = TF.var(:fedd)
    for k in a_-1:b_+2
        sE = 0.0; sF = 0.0
        for g in 1:G; e_ = TX.MGEN[][g][k]; sE += e_; sF += TX.MGFEDD[][g][k] * e_; end
        sE > 0 && (fe[k] = sF / sE)
    end
    TX.MGVEFOK[] = true
    nothing
end
"(MGROUPS) a record per snapshot in POPIII_OUT/mg_rec.jls: groups, grey field, state, ratios, star-light rates."
function mg_record(k)
    a_, b_ = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)); rg = a_-1:b_+2; sd = PM.SLDIAG[]
    rec = (; step = k, t = TF.getscalar(:timen), dt = TF.getscalar(:timen) - TF.getscalar(:timeo), r = TF.var(:rn)[rg], ρ = TF.var(:dn)[rg], T = TF.var(:tn)[rg],
           E = TF.var(:ern)[rg], F = TF.var(:frn)[rg], chi = TF.var(:chifn)[rg], kE = TF.var(:xken)[rg], kP = TF.var(:xkpn)[rg], fedd = TF.var(:fedd)[rg],
           Eg = [TX.MGEN[][g][rg] for g in 1:TX.MGN[]], Fg = [TX.MGFN[][g][rg] for g in 1:TX.MGN[]],
           rE = [TX.MGRE[][g][rg] for g in 1:TX.MGN[]], rP = [TX.MGRP[][g][rg] for g in 1:TX.MGN[]], rR = [TX.MGRR[][g][rg] for g in 1:TX.MGN[]],
           k1 = [TX.MGK1[][g][rg] for g in 1:TX.MGN[]], edges = copy(TX.MGEDGE[]),
           fg = TX.MGVEFOK[] ? [TX.MGFEDD[][g][rg] for g in 1:TX.MGN[]] : nothing, gr = copy(TX.MGGEDDR[]),
           hplus = TX.MGVEFOK[] ? [TX.MGHPLUS[][g][rg] for g in 1:TX.MGN[]] : nothing,
           n2 = length(PM.MGN2[]) >= b_ + 2 ? PM.MGN2[][rg] : nothing,
           sl = sd === nothing ? nothing : (; kph = sd.kph, Te = sd.Te, L = sd.L, kHm = sd.kHm[rg], kH = sd.kH[rg], kLW = sd.kLW[rg]),
           xHp = [ChemNet.full(mat.x[z], mat.net)[:Hp] for z in a_:b_], xH2 = [ChemNet.full(mat.x[z], mat.net)[:H2] for z in a_:b_],
           wrec = [z <= length(PM.MGWREC[]) ? PM.MGWREC[][z] : nothing for z in rg])
    open(io -> serialize(io, rec), joinpath(get(ENV, "POPIII_OUT", "."), "mg_rec.jls"), "a")
    nothing
end
const PROFSAMPLE = haskey(ENV, "CPLPROFILE") ? parse(Int, ENV["CPLPROFILE"]) : -1     # step at which a 200-step Profile sample starts
const PROFSAMPLEN = parse(Int, get(ENV, "CPLPROFILEN", "200"))                       # (PERF512) its length in steps
PROFSAMPLE >= 0 && @eval using Profile
# (PERF) CPLPROFALLOC=rate: an allocation profile (Profile.Allocs, sampling rate) over the same steps, summarised by the innermost
# repository frame and by the chain of the three innermost repository frames (bytes of the sampled allocations / rate)
const PROFALLOC = parse(Float64, get(ENV, "CPLPROFALLOC", "0"))
function alloc_report(io, rate)
    r = Profile.Allocs.fetch(); byf = Dict{String,Float64}(); byc = Dict{String,Float64}(); tot = 0.0
    for a in r.allocs
        fr = [f for f in a.stacktrace if occursin(".jl", string(f.file)) && !startswith(string(f.file), "@") && !occursin("/share/julia/", string(f.file)) &&
              !occursin("/.julia/", string(f.file))]
        key1 = isempty(fr) ? "(no repository frame)" : string(basename(string(fr[1].file)), ":", fr[1].line, " ", fr[1].func)
        keyc = isempty(fr) ? key1 : join((string(basename(string(f.file)), ":", f.line, " ", f.func) for f in fr[1:min(3, end)]), " < ")
        byf[key1] = get(byf, key1, 0.0) + a.size / rate; byc[keyc] = get(byc, keyc, 0.0) + a.size / rate; tot += a.size / rate
    end
    @printf(io, "allocation profile: %.1f GB estimated (%d samples, rate %g)\n\nby innermost repository frame:\n", tot / 1e9, length(r.allocs), rate)
    for (k, v) in sort(collect(byf); by = last, rev = true)[1:min(60, end)]; @printf(io, "%10.3f GB  %5.1f %%  %s\n", v / 1e9, 100v / tot, k); end
    println(io, "\nby the three innermost repository frames:")
    for (k, v) in sort(collect(byc); by = last, rev = true)[1:min(80, end)]; @printf(io, "%10.3f GB  %5.1f %%  %s\n", v / 1e9, 100v / tot, k); end
end
function prof_sample(js)
    if js == PROFSAMPLE
        Base.invokelatest(() -> (Profile.init(n = 10^7, delay = 0.01); Profile.clear(); PROFALLOC > 0 ? (Profile.Allocs.clear(); Profile.Allocs.start(sample_rate = PROFALLOC)) : Profile.start_timer()))
    elseif js == PROFSAMPLE + PROFSAMPLEN
        PROFALLOC > 0 && return Base.invokelatest() do
            Profile.Allocs.stop()
            open(io -> alloc_report(io, PROFALLOC), joinpath(get(ENV, "POPIII_OUT", "."), "cplallocs_$(js).txt"), "w")
        end
        Base.invokelatest() do
            Profile.stop_timer()
            open(joinpath(get(ENV, "POPIII_OUT", "."), "cplprofile_$(js).txt"), "w") do io
                Profile.print(IOContext(io, :displaysize => (100000, 250)); format = :flat, sortedby = :count, mincount = 50)
            end
            # (PERF512) also the call tree, all threads merged, with the module-qualified frames of the repository code
            open(joinpath(get(ENV, "POPIII_OUT", "."), "cplprofile_tree_$(js).txt"), "w") do io
                Profile.print(IOContext(io, :displaysize => (100000, 400)); format = :tree, mincount = 100, maxdepth = 60, noisefloor = 2.0, C = false, combine = true)
            end
        end
    end
end
function prof_report()
    js = Int(TF.getscalar(:jstep)); PROFSAMPLE >= 0 && prof_sample(js)
    (CPLPROF && js % CPLPROFEVERY == 0 && js != PROFSTEP[]) || return nothing
    PROFSTEP[] = js
    @printf("CPLPROF step %d:", js); for s_ in (:matgen, :matslv, :maxdel, :update, :peq, :timstp); @printf(" %s %.1f s/%d %.1f GB", s_, get(PROFT, s_, 0.0), get(PROFN, s_, 0), get(PROFB, s_, 0.0) / 1e9); end
    empty!(PROFB)
    @printf("  cpl_zone jacobians %d (%.1f s), reused %d\n", PM.CPLPN[], PM.CPLPT[] * 1e-9, PM.CPLPH[]); flush(stdout)
    empty!(PROFT); empty!(PROFN); PM.CPLPN[] = 0; PM.CPLPH[] = 0; PM.CPLPT[] = 0
    nothing
end
try
    if haskey(ENV, "TABAUDIT") && mat isa PM.NetGas  # (diagnostic, then exit) the table zones below the hot branch: table vs the network's LTE state
        tabaudit(); exit(0)
    end
    if haskey(ENV, "DTMAXTEST")                    # (test) Δt ≤ DTMAXTEST s, also for the first step: the Δt → 0 limit of a scheme
        dtmx = parse(Float64, ENV["DTMAXTEST"]); push!(TC.DT_HOOKS, () -> dtmx)
        TF.setscalar!(:dtime, min(TF.getscalar(:dtime), dtmx)); @printf("DTMAXTEST: Δt ≤ %.3g s\n", dtmx)
    end
    if COUPLED
        mat isa PM.NetGas || error("COUPLED=1 needs the network gas (CHEM=net...)")
        TF.PARAMS[:meqn] >= 9 || error("COUPLED=1: meqn = $(TF.PARAMS[:meqn]) < 9")
        ChemNet.CK[] = true; TX.COUPLED[] = true; TF.setscalar!(:neqn, 9); TX.cpl_alloc!(length(TF.var(:dn)))
        TX.CPL_ZONE[] = (k, ρ, T, a, b) -> PM.cpl_zone(mat, k, ρ, T, a, b)
        if PM.SLIMP                                # SLIMP>0: the star light of the coupled zones at the Newton iterate (PM.sl_sweep!, cpl_zone_q)
            PM.SLBLEND[] || error("SLIMP needs SLBLEND=1 (the colour temperature and dilution of the iterate are the blend of PM.sl_sweep!)")
            TX.SLIMP[] = true; TX.SL_SWEEP[] = PM.sl_sweep!; TX.CPL_ZONEQ[] = (k, ρ, T, a, b, q, t, fs, g, dp) -> PM.cpl_zone_q(mat, k, ρ, T, a, b, q, t, fs, g, dp)
            pr("SLIMP: star light of the coupled zones at the Newton iterate: dilution W, colour temperature, Lyman-Werner shielding, Lyα escape (local variables of the zone's chemistry and NLOIT opacity)\n")
        end
        get(ENV, "CPLLOCAL", "1") == "1" && (TX.CPL_LOCAL[] = (k, dt) -> PM.cpl_local!(mat, k, dt))
        TC.CPLH[] = (; newstep = TX.cpl_newstep!, reset = TX.cpl_reset!, update = TX.cpl_update!, dmax = TX.cpl_dmax!)
        if haskey(ENV, "CPLJCHECK")                # Jacobian test of the coupled rows (TitanResidFast.cpl_jcheck), once
            zs = parse.(Int, split(ENV["CPLJCHECK"], ",")); once = Ref(true); g0 = TD.ROUTINES[:matgen]
            TD.ROUTINES[:matgen] = () -> (once[] && (once[] = false; a_ = Int(TF.getscalar(:ngrs)); for z in zs; TX.cpl_jcheck(a_ + z - 1); TX.cplz(a_ + z - 1) && PM.cpl_gcheck(mat, a_ + z - 1); end); g0())
        end
        if haskey(ENV, "SZJCHECK")                 # Jacobian test of the radiation and energy rows (TitanResidFast.rad_jcheck), once
            zs = parse.(Int, split(ENV["SZJCHECK"], ",")); once = Ref(true); g0 = TD.ROUTINES[:matgen]
            TD.ROUTINES[:matgen] = () -> (once[] && (once[] = false; a_ = Int(TF.getscalar(:ngrs)); for z in zs; TX.rad_jcheck(a_ + z - 1); end); g0())
        end
        let f0 = TD.STEPFN[]
            TD.STEPFN[] = () -> (PM.cpl_prepare!(mat); f0(); nothing)
        end
        @printf("COUPLED: x_H+ and x_H2 of the network zones are unknowns 8 and 9 of the hydro Newton (neqn 9, meqn %d)\n", TF.PARAMS[:meqn])
    elseif Int(TF.getscalar(:neqn)) == 9
        TF.setscalar!(:neqn, 7)                    # (restart from a COUPLED checkpoint)
    end
    COUPLED && !PM.MGON && Int(TF.getscalar(:neqn)) > 9 && TF.setscalar!(:neqn, 9)   # (restart from an MGROUPS checkpoint without groups)
    PM.MGON && mg_install!()
    if get(ENV, "EDIAG", "0") == "2"               # (diagnostic, default off) EDIAG=1 plus the per-zone breakdown (TX.zone_diag) to ediagz.jls
        let f_ = joinpath(get(ENV, "POPIII_OUT", "."), "ediag.csv"), fz = joinpath(get(ENV, "POPIII_OUT", "."), "ediagz.jls")
            isfile(f_) || open(io -> println(io, "jstep,kind,dt,imb,dE,fin,fout"), f_, "w")
            TD.STAGE_HOOK[] = function (kind, dt)
                d_ = TX.stage_diag(dt)
                open(io -> @printf(io, "%d,%s,%.6e,%.6e,%.6e,%.6e,%.6e\n", Int(TF.getscalar(:jstep)), kind, dt, d_.imb, d_.dE, d_.fin, d_.fout), f_, "a")
                z_ = TX.zone_diag(dt)
                open(io -> serialize(io, merge((; jstep = Int(TF.getscalar(:jstep)), kind, dt), z_)), fz, "a")
                nothing
            end
        end
        pr("EDIAG=2: total-energy rows of every converged stage -> ediag.csv, per zone -> ediagz.jls\n")
    end
    if get(ENV, "EDIAG", "0") == "1"               # (diagnostic, default off) total-energy non-conservation of every converged implicit stage
        # TitanResidFast.stage_diag: the conservative total-energy rows of all zones (ηsw = 0) at the converged state of the stage, summed × the stage's
        # dt = ΔE_total + dt (Φ_out − Φ_in) − (what the rows of the solve actually conserved): one line per stage in POPIII_OUT/ediag.csv
        let f_ = joinpath(get(ENV, "POPIII_OUT", "."), "ediag.csv")
            isfile(f_) || open(io -> println(io, "jstep,kind,dt,imb,dE,fin,fout"), f_, "w")
            TD.STAGE_HOOK[] = function (kind, dt)
                d_ = TX.stage_diag(dt)
                open(io -> @printf(io, "%d,%s,%.6e,%.6e,%.6e,%.6e,%.6e\n", Int(TF.getscalar(:jstep)), kind, dt, d_.imb, d_.dE, d_.fin, d_.fout), f_, "a")
                nothing
            end
        end
        pr("EDIAG: total-energy rows of every converged stage -> ediag.csv\n")
    end
    if get(ENV, "ETOTMON", "0") != "0"              # (diagnostic, default off) total-energy monitor of every converged implicit solve, on the DTERR trace line and in perf.csv
        EM.on = true; prev_ = TD.STAGE_HOOK[]
        TD.STAGE_HOOK[] = function (kind, dt)
            prev_ === nothing || prev_(kind, dt)
            m_ = TX.etot_monitor(dt)
            push!(EM.buf, (Int(TF.getscalar(:jstep)), kind, dt, m_))
            abs(m_.G) > EM.pG && (EM.pG = abs(m_.G)); m_.vmax > EM.pv && (EM.pv = m_.vmax; EM.pz = m_.zmax)
            nothing
        end
        pr("ETOTMON: total-energy rows of every converged solve (TX.etot_monitor) -> ETOT text on the DTERR trace lines, etot_* columns of perf.csv\n")
    end
    Trbdf2.ERRCTL[] || error("DTERR must be set (the step control of the stack; TITAN's heuristic is diagnostic only)")
    Trbdf2.SDON && (Trbdf2.SD_FLAGS[] = k -> @sprintf("eq %d cpl %d", mat.eqzone[k], TX.cplz(k)))      # (stalldiag: SDZONES)
    if Trbdf2.ERRCTL[]                             # DTERR: coverage of the COUPLED zones and the line cap the heuristic would apply, on the trace lines
        Trbdf2.TRACE_EXTRA[] = function ()
            nnet = count(k -> !mat.eqzone[k], a:b); ncpl = count(k -> TX.cplz(k), a:b)
            @sprintf(" | network zones %d, coupled %d | line cap %.3e", nnet, ncpl, LINECAPV[]) * (TD.KRYLOV[] === nothing ? "" : NewtonKrylov.trace_text(Trbdf2.ESTATS[:steps])) * (EM.on ? etot_text() : "")
        end
        @printf("DTERR: error-controlled TR-BDF2 (rtol %g, x floor %g, F floor %g cE, E floor %g aT^4, %s norm, PI %s, safety/facmin/facmax %g/%g/%g, Newton-failure factor %g, BE retry %s); TITAN's heuristic and caps are diagnostic only\n",
                Trbdf2.EOPT.rtol, Trbdf2.EOPT.xfl, Trbdf2.EOPT.ffl, Trbdf2.EOPT.efl, Trbdf2.EOPT.maxnorm ? "max" : "rms", Trbdf2.EOPT.pi, Trbdf2.EOPT.safety, Trbdf2.EOPT.facmin, Trbdf2.EOPT.facmax, Trbdf2.EOPT.fail, Trbdf2.EOPT.be)
    end
    # NWPROBE=<dt1>,<dt2>,... (diagnostic, default off): at step NWPROBEAT+1 [41] of this process, instead of the step, one trapezoidal-stage attempt at each
    # of these step sizes from the live state (as a retry after a Newton failure starts, Trbdf2.step_err!), traced if NWTRACE is set; then the process ends.
    # NWFD=<zlo>,<zhi>,<iter>: with it, at Newton iteration <iter> of the first attempt, finite differences against the assembled Jacobian for the unknowns
    # of zones zlo..zhi (zone numbers as in the divergence message), see nwfd!.
    if haskey(ENV, "NWPROBE") && haskey(ENV, "NWDUMP")      # NWDUMP=zlo,zhi: table of the photospheric zones at the probe state
        Trbdf2.NWDUMP[] = function ()
            zl, zh = parse.(Int, split(ENV["NWDUMP"], ",")); a_ = Int(TF.getscalar(:ngrs)); b_ = Int(TF.getscalar(:ngre))
            rn = TF.var(:rn); dn = TF.var(:dn); tn = TF.var(:tn); ch = TF.var(:chifn); fe = TF.var(:fedd); er = TF.var(:ern); fr = TF.var(:frn)
            ro = TF.var(:ro); dmdt = TF.var(:dmdt); un = TF.var(:un)
            tau = zeros(length(rn)); for k in b_:-1:a_; tau[k] = (k < b_ ? tau[k+1] : 0.0) + ch[k] * dn[k] * (rn[k+1] - rn[k]); end
            pr("NWDUMP step %d t %.9e dt(last) %.4e; geddr %.5f GXEDD %.5f; zones %d..%d (relative numbering)\n", Int(TF.getscalar(:jstep)), TF.getscalar(:timeo), TF.getscalar(:dtime), TF.getscalar(:geddr), TD.GXEDD[], zl, zh)
            pr("  z   r/Rsun     T        rho       tau(z..out) chi*rho*dr  x_H+     x_H2    f        h(HEDD)  hx       fx       tiew(z) E_r/aT4  F/cE    u_km/s  urel_if_km/s dmdt_if  hmean/avchi eq\n")
            for z in zl:zh
                k = a_ + z - 1; y = ycomp(k); aT4 = 7.5657e-15 * tn[k]^4
                hm = 2 / (1 / ch[k-1] + 1 / ch[k]); sfr = (isdefined(TX, :SFRATIO) && k <= length(TX.SFRATIO)) ? TX.SFRATIO[k] : NaN
                pr("%3d %8.3f %8.1f %.3e %.3e %.3e %.3e %.3e %.5f %+.4f %+.4f %+.4f %.3f %.3e %.3f %+.2f %+.2f %+.3e %.4f %s\n", z, rn[k] / Rsun, tn[k], dn[k], tau[k], ch[k] * dn[k] * (rn[k+1] - rn[k]), y[1], y[3], fe[k],
                   k <= length(TD.HEDD) ? TD.HEDD[k] : NaN, k <= length(TD.HXEDD) ? TD.HXEDD[k] : NaN, k <= length(TD.FXEDD) ? TD.FXEDD[k] : NaN, k <= length(TD.TIEW) ? TD.TIEW[k] : NaN,
                   er[k] / aT4, fr[k] / (2.998e10 * er[k]), un[k] / 1e5, (un[k] - (rn[k] - ro[k]) / TF.getscalar(:dtime)) / 1e5, dmdt[k], sfr, mat.eqzone[k] ? "table" : "net")
            end
        end
    end
    if haskey(ENV, "NWPROBE")
        if haskey(ENV, "NWFD")
            let w = parse.(Int, split(ENV["NWFD"], ",")), mg0 = TD.ROUTINES[:matgen], done = Ref(false)
                TD.ROUTINES[:matgen] = function ()
                    mg0()
                    (!isempty(Trbdf2.NWPROBE[]) && !done[] && Int(TF.getscalar(:iter)) == w[3]) || return
                    done[] = true
                    nwfd!(w[1], w[2], w[3]); flush(stdout); exit(0)
                end
            end
        end
        let f0 = TD.STEPFN[], cnt = Ref(0), nat = parse(Int, get(ENV, "NWPROBEAT", "40"))
            TD.STEPFN[] = () -> (cnt[] += 1; cnt[] == nat + 1 && (Trbdf2.NWPROBE[] = parse.(Float64, split(ENV["NWPROBE"], ",")); get(ENV, "NKRYLOV", "0") != "0" && (TD.KRYLOV[] = NewtonKrylov.krylov!)); f0(); nothing)
        end
    end
    if Trbdf2.SDON && (get(ENV, "SDLAG", "0") == "1" || haskey(ENV, "SDFRZG"))     # stalldiag: log / freeze the quantities lagged over one step (the checkpoint's `lag` tuple)
        let f0 = TD.STEPFN[], cnt = Ref(0), saved = Ref{Any}(nothing), grp = split(get(ENV, "SDFRZG", ""), ","), fj = parse(Int, get(ENV, "SDFRZJ", "0"))
            nums(p) = p === nothing ? "-" : join((@sprintf("%s=%.4g", f_, Float64(getfield(p, f_))) for f_ in propertynames(p) if getfield(p, f_) isa Real), ",")
            vecs() = (("GASQ", TX.GASQ[]), ("RAYQ", TX.RAYQ[]), ("CQB", TX.CQB[]), ("YQB", TX.YQB[]), ("CHEMB", TX.CHEMB[]))
            TD.STEPFN[] = function ()
                cnt[] += 1; a_ = Int(TF.getscalar(:ngrs))
                zk_ = collect(Trbdf2.SDZ); ks = [a_ + z - 1 for z in zk_]
                if haskey(ENV, "SDFRZG") && cnt[] == fj       # capture
                    sv = Dict{String,Any}()
                    ("gasq" in grp) && for (nm, v) in vecs()[1:4]; sv[nm] = isempty(v) ? nothing : [v[k] for k in ks]; end
                    ("nlo" in grp) && (sv["NLOPC"] = [k <= length(PM.NLOPC[]) ? deepcopy(PM.NLOPC[][k]) : nothing for k in ks]; sv["NLOST"] = [k <= length(PM.NLOST[]) ? deepcopy(PM.NLOST[][k]) : nothing for k in ks])
                    ("sl" in grp) && (sv["CHEMB"] = isempty(TX.CHEMB[]) ? nothing : [TX.CHEMB[][k] for k in ks]; sv["SLD"] = PM.SLDIAG[] === nothing ? nothing : [(; kHm = PM.SLDIAG[].kHm[k], kH = PM.SLDIAG[].kH[k], kLW = PM.SLDIAG[].kLW[k], fsh = PM.SLDIAG[].fsh[k], tl3 = deepcopy(PM.SLDIAG[].tl3[k])) for k in ks])
                    saved[] = sv; pr("SDFRZ: captured groups %s at step %d of this process\n", join(grp, ","), cnt[])
                elseif haskey(ENV, "SDFRZG") && saved[] !== nothing && cnt[] > fj       # restore
                    sv = saved[]
                    for (nm, v) in vecs(); haskey(sv, nm) && sv[nm] !== nothing && for (i, k) in enumerate(ks); v[k] = sv[nm][i]; end; end
                    haskey(sv, "NLOPC") && for (i, k) in enumerate(ks); PM.NLOPC[][k] = deepcopy(sv["NLOPC"][i]); PM.NLOST[][k] = deepcopy(sv["NLOST"][i]); end
                    if haskey(sv, "SLD") && sv["SLD"] !== nothing
                        d_ = PM.SLDIAG[]; for (i, k) in enumerate(ks); d_.kHm[k] = sv["SLD"][i].kHm; d_.kH[k] = sv["SLD"][i].kH; d_.kLW[k] = sv["SLD"][i].kLW; d_.fsh[k] = sv["SLD"][i].fsh; d_.tl3[k] = deepcopy(sv["SLD"][i].tl3); end
                    end
                end
                if get(ENV, "SDLAG", "0") == "1"
                    er_, tn_ = TF.var(:ern), TF.var(:tn)
                    for (i, k) in enumerate(ks)
                        g_(v) = isempty(v) ? NaN : v[k]
                        pc = k <= length(PM.NLOPC[]) ? PM.NLOPC[][k] : nothing; ps = k <= length(PM.NLOST[]) ? PM.NLOST[][k] : nothing
                        sd_ = PM.SLDIAG[]
                        @printf("SDLAG n %d zone %d | E %.5f T %.2f | GASQ %.4e RAYQ %.4e CQB %.4e YQB %.4e CHEMB %.4e | NLOPC %s | NLOST.p %s | SLD kHm %.4e kH %.4e kLW %.4e\n", cnt[], zk_[i], er_[k], tn_[k], g_(TX.GASQ[]), g_(TX.RAYQ[]), g_(TX.CQB[]), g_(TX.YQB[]), g_(TX.CHEMB[]),
                                pc === nothing ? "-" : nums(pc), ps === nothing ? "-" : nums(ps.p), sd_ === nothing ? NaN : sd_.kHm[k], sd_ === nothing ? NaN : sd_.kH[k], sd_ === nothing ? NaN : sd_.kLW[k])
                    end
                end
                f0(); nothing
            end
        end
    end
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
        FLUSHOUT && flush(stdout)                    # (FLUSHOUT=1: the step's diagnostics reach the output file now)
        NANSCAN && nanscan("after the hydro step", k)
        dT = PM.cpl_chem!(mat)
        NANSCAN && nanscan("after the chemistry", k)
        prof_report()
        GRIDEQP === nothing || grideq_step!()
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
        PM.MGON && SNAPEVERY > 0 && k % SNAPEVERY == 0 && mg_record(k)
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
            tl0 = time(); Λl, Lr, _ = linecool(rr, uu, dd, tt, nh2, k); TLINE[] += time() - tl0
            tcmin = Inf
            for j in 1:nsh
                mat.Λext[a+j-1] = Λl[j]; mat.Tref[a+j-1] = tt[j]; mat.slope[a+j-1] = 0.0
                ue = 1.5 * TF.var(:pgn)[a+j-1] / dd[j]          # thermal energy scale of the cooling time (diagnostic: the line cap of DTDIAG and DTERR's trace)
                Λl[j] != 0 && abs(ue * dd[j] / Λl[j]) < tcmin && (tcmin = abs(ue * dd[j] / Λl[j]); JCMIN[] = j)
            end
            LLINE[] = Lr[end]
            if get(ENV, "DTDIAG", "0") == "1"
                j_ = JCMIN[]; pr("        line cap: t_cool %.3e s at zone %d (m %.4e Msun, T %.1f, n %.3e, e+U %.3e, Lambda %.3e); dtime before %.3e, smax %.3e\n", tcmin, j_, 4π * xmi(a+j_) / Msun, tt[j_], dd[j_] / (1.22PC.mH),
                                TF.var(:egn)[a+j_-1] + PC.UOFF_TH, Λl[j_], TF.getscalar(:dtime), TF.getscalar(:smax))
                LINEIMP && pr("        LINEIMP: zone line models from EP %d, ALI %d, held ALI rates %d (ALI/EP factor clamped in %d); no t_cool cap\n", LISTAT...)
            end
            LINECAPV[] = tcmin                       # (diagnostic only: DTERR sets Δt; LINEIMP has no t_cool cap)
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
                   (RAD ? @sprintf("; E_r: zone %d %.2e, rho: zone %d %.2e", kx[6], dx[6], kx[3], dx[3]) : ""))
            end
        end
        if ENVEUL !== nothing
            envsink!()
            if k % PRINTEVERY == 0
                if TX.SLABFLUX[] && !isempty(TX.SFRATIO)   # SLABFLUX: the interfaces whose optical depth changed by more than x1.5 (zone numbers from the inner boundary)
                    jf = [j for j in a+1:min(b, length(TX.SFRATIO)) if !(1 / 1.5 < TX.SFRATIO[j] < 1.5)]
                    pr("        SLABFLUX: %d interface(s) with dtau/dtau_T outside [2/3, 3/2]%s\n", length(jf), join((@sprintf(" | %d: %.3g", j - a + 1, TX.SFRATIO[j]) for j in jf[1:min(end, 8)]), ""))
                end
                if TX.SUBZONE[] && !isempty(TX.SZOMEGA)    # SUBZONE: the interfaces with a front weight > 0.05 (zone numbers from the inner boundary)
                    jf = [j for j in a+1:min(b, length(TX.SZOMEGA)) if TX.SZOMEGA[j] > 0.05]
                    pr("        SUBZONE: %d front interface(s)%s\n", length(jf), join((@sprintf(" | %d: w %.3f dtau/dtau_T %.3f", j - a + 1, TX.SZOMEGA[j], TX.SZRATIO[j]) for j in jf), ""))
                end
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
                GRIDEQP === nothing || pr("%s", grideq_line())
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
        if TIMING && k % PRINTEVERY == 0
            pr("        timing since last report: step (hydro+radiation) %.1f s, chemistry %.1f s, H2 lines %.1f s, other output %.1f s\n",
               TSTEP[], TCHEM[], TLINE[], time() - tnow - 0.0 + TOTH[] - TCHEM[] - TLINE[])
            REC.on && recperf!(k, tnow)
            TSTEP[] = TCHEM[] = TLINE[] = TOTH[] = 0.0
        else
            TOTH[] += time() - tnow
        end
        TMARK[] = time()
        # CKDUMP=k1:k2[,k3:k4…] (diagnostic, default off): also write a full checkpoint <SAVE>_k<step> at each listed step (restart A/B tests)
        haskey(ENV, "CKDUMP") && haskey(ENV, "SAVE") && any(r_ -> (p_ = parse.(Int, split(r_, ":")); p_[1] <= k <= p_[end]), split(ENV["CKDUMP"], ",")) && checkpoint(ENV["SAVE"] * "_k" * string(k), k; jinc = 1)
        haskey(ENV, "SAVE") && CKEVERY > 0 && k % CKEVERY == 0 && checkpoint(ENV["SAVE"] * "_auto" * string((k ÷ CKEVERY) % 2), k; jinc = 1)
    end)
catch e
    pr("stopped: %s\n", first(sprint(showerror, e), 400))
    get(ENV, "DEBUG", "0") == "1" && Base.display_error(stderr, e, catch_backtrace())
end
Trbdf2.ERESTOL[] > 0 && pr("ERESTOL %g: %d stage solves tested, %d extra Newton iterations, %d left above the tolerance (%d of them at niter); largest accepted measure %.4g L_sun\n",
                           Trbdf2.ERESTOL[], Trbdf2.ERES[:checks], Trbdf2.ERES[:extra], Trbdf2.ERES[:unmet], Trbdf2.ERES[:niter], Trbdf2.ERESMAX[])
pr("run ended: %s; jstep %d jstepe %g timen %.6e tmax %.6e\n", string(RUNEND[]), Int(TF.getscalar(:jstep)), Float64(TF.getscalar(:jstepe)), TF.getscalar(:timen), TF.getscalar(:tmax))
if haskey(ENV, "SAVE")      # never overwrite a checkpoint with a non-finite state (a NaN stop would otherwise be saved and restart as NaN)
    if all(v -> all(isfinite, TF.var(v)[a:b]), (:tn, :dn, :ern, :un))
        checkpoint(ENV["SAVE"], KLAST[]; jinc = RUNEND[] isa Tuple && RUNEND[][2] === :tmax ? 1 : 0); pr("checkpoint %s at step %d\n", ENV["SAVE"], KLAST[])
    else
        pr("checkpoint %s NOT written: the state at step %d is not finite (restart from the last _auto checkpoint)\n", ENV["SAVE"], KLAST[])
    end
end
PM.NLOIT && pr("NLOIT: %d opacity-ratio evaluations (%.1f s), %d slope refreshes; CK expansions %.1f s\n", PM.NLOITN[], PM.NLOITNS[] * 1e-9, PM.NLOITJN[], PM.NLOCNS[] * 1e-9)
NLTE == "hyb" && pr("hybrid NLTE: %d ALI calls\n", NALICALL[])
TD.STEPFN[] === Trbdf2.step! && pr("integrator TR-BDF2, statistics %s\n", string(Trbdf2.STATS))
TD.KRYLOV[] === nothing || pr("%s\n", NewtonKrylov.summary_text())
Trbdf2.ERRCTL[] && pr("DTERR statistics %s; limiting variable (accepted steps) %s; TR-BDF2 %s\n", string(Trbdf2.ESTATS), string(Trbdf2.LIMHIST), string(Trbdf2.STATS))
close(io)
