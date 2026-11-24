# The run itself (statements of the former popiii_m2c.jl, in their original order), evaluated in module PopIIIRun by PopIIIRun.main at run time:
# settings are read from ENV / ARGS here, never at load time. Function definitions live in PopIIIRun.jl.
# default since G2b: θ = 0.55 lets the z = 19 core oscillate at 1e11–1e13 (three-body heating phase) where TR-BDF2 is smooth
get(ENV, "INTEGRATOR", "trbdf2") == "trbdf2" && (TD.STEPFN[] = Trbdf2.step!)
cool = isempty(ARGS) ? "ra04" : ARGS[1]
 nsh = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100
cf = cool == "thin" ? PC.cool_H2_thin : cool == "ra04lines" ? PC.cool_H2_ra04 :
     cool == "lines" ? PC.cool_CIE_ra04 :                     # H₂ lines by transfer (H2Lines), continuum: RA04 CIE
     (ρ, T, y) -> PC.cool_H2_ra04(ρ, T, y) + PC.cool_CIE_ra04(ρ, T, y)
haskey(ENV, "TCMB") && (cf = PC.cool_floor(cf, parse(Float64, ENV["TCMB"])))
# RAD=1 (the stack; from the IC): TITAN's radiation (E_r, F_r, Eddington factors) and the StellarOpacity continuum; the RA04
# CIE escape-factor term is then not used (the continuum is in the transfer), the H₂ lines stay as Λext
RAD = get(ENV, "RAD", "0") == "1"
RAD && (cf = (ρ, T, y) -> zero(T))
# POPIII_OUT: directory for the run log and the snapshot file (default lib/popiii); checkpoints follow CKDIR
OUTDIR = get(ENV, "POPIII_OUT", joinpath(TF.ROOT, "lib", "popiii"))
mkpath(OUTDIR)
io = open(joinpath(OUTDIR, "m2c_$(get(ENV, "IC", "p100") == "p100" ? "" : get(ENV, "IC", "") * "_")$(cool)_n$(nsh)$(get(ENV, "TAG", "")).log"), haskey(ENV, "RESTART") ? "a" : "w")
# (a continued run appends)
ic = get(ENV, "IC", "abn02red")
dm1 = 3e-7Msun
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
        for lm_ in (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0)
            println("$ic: x_e(M = 10^", lm_, " Msun) = ", round(xefin(10^lm_), sigdigits = 3), xefin_flat ? " (XEFLAT)" : " (recombined for one t_ff from XE0)")
        end
        mv1_fin = max(10^lMv_fin[1], mcore_fin); vv1_fin = vfin_(log10(mv1_fin))   # (abn02fin: the first table point, v_fin[1])
        rM1fin = let m1 = mv1_fin * Msun, j = clamp(searchsortedfirst(Mg, m1), 2, length(Mg))
            rg[j-1] + (m1 - Mg[j-1]) / (Mg[j] - Mg[j-1]) * (rg[j] - rg[j-1])
        end
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
CHEM = get(ENV, "CHEM", "netmodernli")
# net / netmodern / netmodernli: NetGas (G5); the RHFC02 PrimordialGas material is gone
CHEM in ("net", "netmodern", "netmodernli") || error("CHEM=$CHEM: the network material net, netmodern or netmodernli")
begin
    # G5: ChemNet network; equilibrium table where the ideal-gas pressure ≥ PSWITCH (EOS validity, RHFC02's 1e4 dyn cm⁻²)
    netw, phot = CHEM == "net" ? (ChemRates.rhfc02(; k5 = k5f), ChemNet.Photo[]) : ChemRates.modern(; li = CHEM == "netmodernli")
    # SL3=1 (with STARLIGHT=1): three-level H atom in the lit network zones (PM.tl3_reactions; rates are zero elsewhere)
    get(ENV, "SL3", "0") == "1" && (netw = PM.tl3_network(netw))
    ChemNet.HM88 && (netw = ChemNet.hm88_network(netw))      # HM88=1: H recombination × Z_w/Z_el (ChemNet.HM88)
    tcmb = parse(Float64, get(ENV, "TCMB", "0"))
    envn = (Trad = tcmb, J = tcmb > 0 ? J_cmb : J_zero, static = true)   # same J in every zone: tabulated photo-rates
    global mat = PM.NetGas([xinit(y) for y in ystore], netw, phot, coolnet; env = envn, pswitch = parse(Float64, get(ENV, "PSWITCH", "1e4")))
end
# artificial viscosity: viscous length l = ql1 r ≈ 2 Δr (the SOD state's ql0 = 1e-4 cm is zero at these scales, M2g);
# Δr/r ≈ 0.04 (core) … 0.12 (envelope) on the geometric mass grid; quadratic term only (cq2), as PORTING_NOTES P7
av = (ql0 = 0.0, ql1 = parse(Float64, get(ENV, "AVQL1", "0.15")), cq1 = 0.0, cq2 = parse(Float64, get(ENV, "AVCQ2", "1.0")))
PS.build!(; rI, rho, T = Tz, u = uI, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = 0.0,
          selfgrav = true, material = mat, av)
TX.AV_TENSION[] = get(ENV, "AVTENSION", "1") == "1"
# 1 (TITAN) | 0 | f in (0,1) | radial | radial,f  (TX.AV_TFRAC)
let p = split(get(ENV, "AVTENSION", "1"), ","), f = tryparse(Float64, p[end])
    TX.AV_TRAD[] = p[1] == "radial"; TX.AV_TFRAC[] = f !== nothing && 0 < f < 1 ? f : 0.0
end
TX.AV_ZONE[] = parse(Float64, get(ENV, "AVZONE", "0"))
# G1: viscous length = AVZONE × zone width (0 = ql1 r̄ above)
TX.REGULAR_CENTRE[] = true
 TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
PS.install_julia!(BlockSolver)
 TD.MSGIO[] = get(ENV, "MSG", "0") == "1" ? stdout : devnull
# (Newton messages also before RAD switches on)
a = Int(TF.getscalar(:ngrs))
 b = Int(TF.getscalar(:ngre))
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
TF.setscalar!(:thet, parse(Float64, get(ENV, "THETA", "0.55")))
# 0.5 lets the centre ring (M2c)
TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.01"))
tff = sqrt(3π / (32Z0.G * ρc))
TF.setscalar!(:dtime, 1e-3tff)
TF.setscalar!(:tmax, 1e3tff)
pr("M2c collapse [IC %s]: cooling %s, %d shells (q = %.4f, dm1 = 3e-7 Msun), R = %.3e cm, rho_c/rho_edge = %.1f, T_edge = %.1f K, p_ext = %.3e\n",
   ic, cool, nsh, q, rI[end], rho[1] / rho[end], Tz[end], TX.PEXT[])
pr("theta = %.2f, C_dyn = %g, k5 = %s, chemistry %s, %s, AV ql1 = %g cq2 = %g, central t_ff0 = %.3e yr\n", TF.getscalar(:thet), TC.DT_DYN[], k5s, CHEM, mat isa PM.NetGas ? @sprintf("equilibrium table where p ≥ %g", mat.pswitch) : @sprintf("rho_eq = %g, eq EOS %s", mat.ρeq, mat.eqtable ? "Saha+SCVH95" : "Saha"), av.ql1, av.cq2, tff / 3.156e7)
# ON98 Fig. 3a (run A), read off: (log n_c, log T_c)
on98 = [(6.3, 2.40), (8.0, 2.60), (10.0, 2.75), (12.0, 3.05), (15.0, 3.30), (16.0, 3.35), (20.0, 3.80)]
pr("%7s %8s %10s %9s %9s %10s %10s %9s %9s\n", "step", "log n_c", "t [yr]", "T_c", "ON98 T_c", "f_H2,c", "u_c [km/s]", "max dT/T", "dt [s]")
nextln = Ref(log10(rho[1] / PC.mH) + 0.25)
SNAPFILE = joinpath(OUTDIR, "snaps_$(get(ENV, "IC", "p100"))_$(cool)_n$(nsh)$(get(ENV, "TAG", "")).jls")
# RUNREC=<dir>: record the run for movies, interactive panels and cost accounting (off by default). Each process segment
# (start step k0, pid) writes snaps_<k0>_<pid>.jls, one serialized snapshot every SNAPEVERY steps, appended and never
# rewritten (a restart opens a new file; readers keep the surviving branch by scanning back in time), and meta_<k0>_<pid>.txt
# (the settings the code reads, git state, threads). perf.csv gets one row per PRINTEVERY steps (needs TIMING=1): the wall-time
# split, Newton/TR-BDF2 counters, ALI calls, chemistry fallbacks, GC time and allocations, the Δt limiter and the core state.
RUNREC = get(ENV, "RUNREC", "")
isempty(RUNREC) || mkpath(RUNREC)
# CONVDIAG=1 (diagnostic, default off; the solution is bit-identical on or off, docs/POPIII_CONVECTION.md M1): with RUNREC, every snapshot carries `cd`, the sub-shock convection
# diagnostics (`convdiag_snap`): TDC intermediates per interface, frozen/equilibrium thermodynamics, N², timescales, the terms of the gas-energy row.
CONVDIAG = get(ENV, "CONVDIAG", "0") == "1"
TX.CONVDIAG[] = CONVDIAG
# Convection = TDC (Stage 2: whenever convection is on, HOT=1; no switch): the turbulent velocity w is the last Newton unknown (TX.TDC_ON, installed below); TDCPAR = A_s,A_d,A_c,alpha_p,alpha_t,alpha_D,alpha_m.
REC = (on = !isempty(RUNREC), io = Ref{Any}(nothing), wall0 = time(), klast = Ref(-1), gct = Ref(Base.gc_num().total_time), gcb = Ref(Base.gc_bytes()),
             dtlim = Ref(""), mstar = Ref(NaN), rstar = Ref(NaN), spn = Ref{Vector{Symbol}}(Symbol[]))
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
# H₂ dissociation with shielding. STARLIGHT=1; SLSHIELD = static (default) | sobolev | none (no shielding) | off (no LW).
if get(ENV, "STARLIGHT", "0") == "1" && mat isa PM.NetGas
    let ph = mat.photos
        PM.STARLIGHT[] = (; pHm = get(ENV, "SLHM", "1") == "1" ? ph[findfirst(p -> p.R === :Hm, ph)] : nothing,
                          pH = nothing,
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
SNAPS = get(ENV, "SNAPS", "0") == "1" ? Any[] : nothing
SNAPNEXT = Ref(log10(rho[1] / PC.mH))
SNAPEVERY = parse(Int, get(ENV, "SNAPEVERY", "0"))
# LINEINC=1: the H₂ levels also see the radiation that falls on the gas from outside the lines (`H2Lines.incident`): the CMB (TCMB)
# and the continuum of the photosphere, B_ν(T_eff) from the sphere of the photosphere's radius (STARLIGHT's, of the last chemistry
# step), diluted. LINEINC=cmb: the CMB alone. Default 0: the lines see only their own radiation, so in the cold outer envelope the
# levels are excited by collisions alone and the line cooling can take the gas below T_CMB. With the incident field the levels are
# pumped as well and the gas term (the net collisional exchange) becomes heating where the radiation is warmer than the
# gas. The ALI rates are net of the absorbed background and can change sign.
LINEINC = get(ENV, "LINEINC", "0")
LINEINCDIAG = get(ENV, "LINEINCDIAG", "0") == "1"
LINEIMP && (PM.LINEIMP[] = true; pr("LINEIMP: H2 line cooling of each zone from its own state (populations at the iterate's T, radiative rates of the last transfer call); no t_cool cap, no IMPCOOL slope\n"))
"(diagnostic, NANSCAN=1) report the first zones whose state is not finite."
NANSCAN = get(ENV, "NANSCAN", "0") == "1"
FLUSHOUT = get(ENV, "FLUSHOUT", "0") == "1"
# (H2Levels is loaded before this module: popiii_m2c.jl includes it into Main, or it is the TitanStack image's)
# (with TITAN_PKG: Main.H2Levels is TitanStack.H2Levels, bound above)
HLV = H2Levels.load()
LTAUCUT = parse(Float64, get(ENV, "LTAUCUT", "10"))
LINECORE = parse(Int, get(ENV, "LINECORE", "1"))
HLL = cool == "lines" ? H2Lines.load() : nothing
NBIN = parse(Int, get(ENV, "NBIN", "30"))
Tstop = parse(Float64, get(ENV, "TSTOP", "6000"))
 nstop = parse(Float64, get(ENV, "NSTOP", "1e18"))
# Checkpoint / restart (S6-A). SAVE=name: checkpoint at the end and every CKEVERY [500] steps (alternating
# <name>_auto0/_auto1); RESTART=name: the set-up above runs as usual (same ENV), then the TITAN state, the material's
# per-zone state, the driver's cooling bookkeeping, the TR-BDF2 counters and the step number are replaced by the
# checkpoint's, so the continuation is the uninterrupted run (checked bit for bit). Files: CKDIR [lib/popiii/ck].
CKDIR = get(ENV, "CKDIR", joinpath(TF.ROOT, "lib", "popiii", "ck"))
"""
THCHAIN=1 (default): the stellar EOS's p, e and MLT coefficients (∇_ad, c_P, δ) inside the Jacobian sweep are evaluated at the zone's values with a
2-partial derivative in (ρ, T), and passed into the outer Dual by the chain rule. Same values and same first derivatives as evaluating
StellarEOS.eos on the outer Dual (whose own nested derivative then carries all of the sweep's partials through every operation), at a
fraction of the cost (it was ~40 % of the run time). THCHAIN=0: the direct evaluation.
"""
THCHAIN = get(ENV, "THCHAIN", "1") == "1"
if HOTON
    DH = parse(Float64, get(ENV, "DH", "2.5e-5"))
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
    SM.TSOLVE[] = T_blend; SM.ZONEOK[] = k -> TF.var(:tn)[k] > PM.HOT[].T1
    # S5's coupling: every step through solve_step! (level factor f_k on ε iterated against the network), the stage ledger, Δt caps on composition
    # and energy; after the step burn! (network along the T, ρ history of the step, the difference to what the residual injected into e) and mix!
    SG[].nuc = true; SG[].dXmax = parse(Float64, get(ENV, "DXMAX", "4e-3")); SG[].efrac = parse(Float64, get(ENV, "EFRAC", "5"))
    SM.REGION_CAP[] = true; SM.install_stage_hook!(SG[])
    push!(TC.DT_HOOKS, () -> min(SG[].dtcap, SM.predict_cap(SG[])))
    SPLIT = (tol = parse(Float64, get(ENV, "SPLITTOL", "3e-3")), maxit = parse(Int, get(ENV, "SPLITIT", "4")), ztol = parse(Float64, get(ENV, "ZSPLITTOL", "3e-3")))
    STEP0 = TD.STEPFN[]
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
end
PM.PSWHYS[] = parse(Float64, get(ENV, "PSWHYS", "1"))
PM.PSWTMIN[] = parse(Float64, get(ENV, "PSWTMIN", "0"))
PM.PSWTEX[] = parse(Float64, get(ENV, "PSWTEX", "0"))
 PM.PSWTAU[] = parse(Float64, get(ENV, "PSWTAU", "0"))
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
ENVEUL = haskey(ENV, "ENVEUL") ? parse.(Float64, split(ENV["ENVEUL"], ",")) : nothing
ENVCENTRE = ENVEUL !== nothing && ENVEUL[1] == 0
# r_in = 0: the mesh runs to a regular centre (ENVCORE=rmin,dlnr,rcb: first interface [R☉], Δln r below rcb); no point mass, wall or sink
ENVCORE = parse.(Float64, split(get(ENV, "ENVCORE", "0.1,0.08,10"), ","))
ENVLIN = parse(Int, get(ENV, "ENVLIN", "1"))
# with a regular centre (ENVEUL=0,...): ENVLIN zones of equal Δr from r = 0 to ENVCORE's first interface (default 1: one central zone)
ENVCONS = get(ENV, "ENVCONS", "0") == "1"
# the ENVEUL remap takes the zone densities from the source's enclosed mass (mass-conserving); default: log-interpolated at the zone centres
ENVSINKTAU = parse(Float64, get(ENV, "ENVSINKTAU", "1e6"))
include(joinpath(@__DIR__, "grideq_driver.jl"))
# GRIDEQ (default off)
# XADV2=1 (default with SL3=1): the abundances cross each interface with the van Leer interface value TITAN uses for e (`advect!`, centred,
# cadv; start-of-step values), explicitly, instead of the donor value. With x_e ~ 0.1–0.5 the first-order scheme carried a different amount of
# ionization energy across each interface than e did: up to ~1 eV per H a step, enough to leave e below u(10 K, x). Falls back to the
# Applied as an explicit correction to the implicit upwind step, at interfaces where less than the donor zone's mass crosses, and per species
# only if it keeps that species ≥ 0 everywhere; `NXADV1` counts the steps where some species lost the correction.
XADV2 = get(ENV, "XADV2", get(ENV, "SL3", "0")) == "1"
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
    pr("GRIDEQ: mode %d (%s), %d zones, weights rho,P,T,E,m,tau %s, gmax %.3g, alph %.3g, kappa %.3g, fv %.3g, fmove %.3g, tau %.3g..%.3g, ctol %.3g, kappa_f %.3g (unused), f_dyn %.3g\n",
       GRIDEQP, GRIDEQP == 1 ? "implicit grid equation" : "prescribed target", nsh, join(GEQW, ","), GEQC[1:10]...)
    GEQGP[1] > 0 && pr("GEQGP (provisional): monitor term w_g %g x (dln rho/dj)^2 x gravity power of the zone / %g L_sun (implicit grid equation only)\n", GEQGP...)
end
JNEXT[] >= 0 && TF.setscalar!(:jstep, JNEXT[])
# restart: the step counter of the uninterrupted run (after the repeated grideq_step!, which ran at the checkpointed step's jstep)
CHEMDIAG = get(ENV, "CHEMDIAG", "0") == "1"
TIMING = get(ENV, "TIMING", "0") == "1"
TMARK = Ref(time())
DTSTOP = parse(Float64, get(ENV, "DTSTOP", "0"))
# stop an unattended run whose Δt has collapsed, before the checkpoints are overwritten
PRINTEVERY = parse(Int, get(ENV, "PRINTEVERY", "500"))
CKEVERY = parse(Int, get(ENV, "CKEVERY", "500"))
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
haskey(ENV, "ZEROSTEP") && ENV["ZEROSTEP"] != "0" && StageCons.configure_zero!()
haskey(ENV, "STAGECONS") && ENV["STAGECONS"] != "0" && (StageCons.configure!(parse(Int, ENV["STAGECONS"])); pr("STAGECONS: stage-consistency diagnostic every %d steps (StageCons.jl)\n", StageCons.EVERY[]))
haskey(ENV, "NWFREEZE") && (TC.NWFREEZE[] = String.(split(ENV["NWFREEZE"], ",")))
# (diagnostic, default off) see TitanCore.NWFREEZE
haskey(ENV, "NWTRACE") && (Trbdf2.NWTRACE[] = open(ENV["NWTRACE"], "a"))
# (diagnostic, default off) Newton trace, see Trbdf2.NWTRACE
haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))
haskey(ENV, "CONV") && TF.setscalar!(:conv, parse(Float64, ENV["CONV"]))
haskey(ENV, "NTRY") && TF.setscalar!(:ntry, parse(Float64, ENV["NTRY"]))
if haskey(ENV, "FAILTFAC"); TD.FAIL_TFAC[] = parse(Float64, ENV["FAILTFAC"]); Trbdf2.TFAC[] = TD.FAIL_TFAC[]; end
haskey(ENV, "ERTHICK") && (TC.ER_THICK[] = parse(Float64, ENV["ERTHICK"]))
# One energy variable, the gas internal energy: the compatible rows (TX.ECOMPAT, below) are the only energy rows (docs/POPIII_ENERGY.md §11, §23); ETOTMON=1 monitors the total energy.
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
end
# Compatible internal-energy rows in every zone, no row switching and no total-energy row (TX.ECOMPAT, docs/POPIII_ENERGY.md §4, §23); implies TIEFRAME
TX.GASROW[] || error("the compatible energy rows need the gas-alone row (GASROW)")
TX.ECOMPAT[] = true
TX.TIEFRAME[] || (TX.TIEFRAME[] = true; pr("ECOMPAT: TIEFRAME switched on\n"))
pr("ECOMPAT: first-law gas row adjoint to the momentum row in every zone, lab-frame radiation energy row\n")
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
# VEFIMP=notie: the same with FLUXTIE off, the moment equations everywhere; nobg: notie without VEF_BG; it / sph: one part (tests).
if get(ENV, "VEFIMP", "0") != "0"
    let v = ENV["VEFIMP"]
        v in ("1", "notie", "nobg", "it", "sph") || error("VEFIMP must be 1, notie, nobg, it or sph")
        TD.VEF_ITER[] = v != "sph"; TD.VEF_SPH[] = v != "it"; TD.VEF_BG[] = v in ("1", "notie")
        if v != "1" && TD.FLUXTIE[][2] > 0
            pr("VEFIMP: FLUXTIE (%g,%g) switched off\n", TD.FLUXTIE[]...)
            TD.FLUXTIE[] = (0.0, 0.0)
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
"""
CPLPROF=1 (diagnostic, default off; also without COUPLED): wall time of matgen (residual + AD Jacobian), matslv, maxdel, update and the
Newton iteration count, and (COUPLED) the local 4-wide Jacobians of PopIIIMaterial.cpl_zone (computed / reused, time); printed every 100 steps.
"""
CPLPROF = get(ENV, "CPLPROF", "0") == "1"
CPLPROFEVERY = parse(Int, get(ENV, "CPLPROFEVERY", "100"))
# (PERF) CPLPROF print interval in steps
if CPLPROF
    for s_ in (:matgen, :matslv, :maxdel, :update, :peq, :timstp)
        let f0 = TD.ROUTINES[s_], s_ = s_
            TD.ROUTINES[s_] = () -> (t0 = time(); b0 = Base.gc_bytes(); try f0() finally PROFT[s_] = get(PROFT, s_, 0.0) + time() - t0; PROFN[s_] = get(PROFN, s_, 0) + 1
                                     PROFB[s_] = get(PROFB, s_, 0.0) + (Base.gc_bytes() - b0) end)
        end
    end
end
PROFSAMPLE = haskey(ENV, "CPLPROFILE") ? parse(Int, ENV["CPLPROFILE"]) : -1
# step at which a 200-step Profile sample starts
PROFSAMPLEN = parse(Int, get(ENV, "CPLPROFILEN", "200"))
# (PERF512) its length in steps
PROFSAMPLE >= 0 && @eval using Profile
# (PERF) CPLPROFALLOC=rate: an allocation profile (Profile.Allocs, sampling rate) over the same steps, summarised by the innermost
# repository frame and by the chain of the three innermost repository frames (bytes of the sampled allocations / rate)
PROFALLOC = parse(Float64, get(ENV, "CPLPROFALLOC", "0"))
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
        TX.CPL_ZONE[] = (k, ρ, T, a, b, xn) -> PM.cpl_zone(mat, k, ρ, T, a, b, xn)
        if PM.SLIMP                                # SLIMP>0: the star light of the coupled zones at the Newton iterate (PM.sl_sweep!, cpl_zone_q)
            PM.SLBLEND[] || error("SLIMP needs SLBLEND=1 (the colour temperature and dilution of the iterate are the blend of PM.sl_sweep!)")
            TX.SLIMP[] = true; TX.SL_SWEEP[] = PM.sl_sweep!; TX.CPL_ZONEQ[] = (k, ρ, T, a, b, q, t, fs, g, dp, xn) -> PM.cpl_zone_q(mat, k, ρ, T, a, b, q, t, fs, g, dp, xn)
            pr("SLIMP: star light of the coupled zones at the Newton iterate: dilution W, colour temperature, Lyman-Werner shielding, Lyα escape (local variables of the zone's chemistry and NLOIT opacity)\n")
        end
        get(ENV, "CPLLOCAL", "1") == "1" && (TX.CPL_LOCAL[] = (k, dt) -> PM.cpl_local!(mat, k, dt))
        TC.CPLH[] = (; newstep = TX.cpl_newstep!, reset = TX.cpl_reset!, update = TX.cpl_update!, dmax = TX.cpl_dmax!)
        if haskey(ENV, "CPLJCHECK")                # Jacobian test of the coupled rows (TitanResidFast.cpl_jcheck), once
            zs = parse.(Int, split(ENV["CPLJCHECK"], ",")); once = Ref(true); g0 = TD.ROUTINES[:matgen]
            TD.ROUTINES[:matgen] = () -> (once[] && (once[] = false; a_ = Int(TF.getscalar(:ngrs)); for z in zs; TX.cpl_jcheck(a_ + z - 1); TX.cplz(a_ + z - 1) && PM.cpl_gcheck(mat, a_ + z - 1); end); g0())
        end
        if haskey(ENV, "LRESPTEST")                # M6 test hook (diagnostic, default off): include the given script once at the first matgen of the run, with `mat` and the loaded modules in Main
            once = Ref(true); g0 = TD.ROUTINES[:matgen]; global LRESP_MAT = mat
            TD.ROUTINES[:matgen] = () -> (once[] && (once[] = false; Base.include(Main, ENV["LRESPTEST"])); g0())
        end
        if haskey(ENV, "SZJCHECK")                 # Jacobian test of the radiation and energy rows (TitanResidFast.rad_jcheck), once
            zs = parse.(Int, split(ENV["SZJCHECK"], ",")); once = Ref(true); g0 = TD.ROUTINES[:matgen]
            TD.ROUTINES[:matgen] = () -> (once[] && (once[] = false; a_ = Int(TF.getscalar(:ngrs)); for z in zs; TX.rad_jcheck(a_ + z - 1); end); g0())
        end
        let f0 = TD.STEPFN[]
            TD.STEPFN[] = () -> (PM.cpl_prepare!(mat); f0(); nothing)
        end
        @printf("COUPLED: x_H+ and x_H2 of the network zones are unknowns 8 and 9 of the hydro Newton (neqn 9, meqn %d)\n", TF.PARAMS[:meqn])
    elseif Int(TF.getscalar(:neqn)) in (9, 10)
        TF.setscalar!(:neqn, 7)                    # (restart from a COUPLED checkpoint)
    end
    COUPLED && !PM.MGON && Int(TF.getscalar(:neqn)) > 9 && TF.setscalar!(:neqn, 9)   # (restart from an MGROUPS checkpoint without groups)
    PM.MGON && mg_install!()
    if TX.CONVECTION[]                              # M3b: the turbulent velocity w as the last Newton unknown (neqn + 1), see TX.TDC_ON
        COUPLED || error("TDC (convection, HOT=1) needs COUPLED=1"); HOTON || error("TDC needs HOT=1 (the convection material)")
        Trbdf2.SCHEME[] === :trbdf2 || error("TDC: DTSCHEME=trbdf2 only (the old w of the other schemes' stages is not combined)")
        isnan(TX.PEXT[]) || error("TDC: PEXT not supported")
        nq_ = Int(TF.getscalar(:neqn)) + 1
        TF.PARAMS[:meqn] >= nq_ || error("TDC: meqn = $(TF.PARAMS[:meqn]) < $nq_")
        TF.setscalar!(:neqn, nq_); TX.TDC_ON[] = true; TX.tdc_alloc!(length(TF.var(:dn)))
        tp_ = parse.(Float64, split(get(ENV, "TDCPAR", "1,1,1,1,0.01,1,0"), ","))      # A_s, A_d, A_c, α_p, α_t, α_D, α_m (B2 best estimate; α_m is not implemented, must be 0)
        length(tp_) == 7 || error("TDCPAR = A_s,A_d,A_c,alpha_p,alpha_t,alpha_D,alpha_m"); tp_[7] == 0 || error("TDC: the eddy viscosity alpha_m is not implemented (TDCPAR 7th = 0)")
        TX.TDCP[] = TX.TDC.TDCPar(; A_s = tp_[1], A_d = tp_[2], A_c = tp_[3], α_p = tp_[4], α_t = tp_[5], α_D = tp_[6], α_m = tp_[7])
        TX.TDC_CEN[] = parse(Float64, get(ENV, "TDCCEN", "0")); atexit(() -> println(TX.buoy_summary()))
        let h0 = TC.CPLH[]
            TC.CPLH[] = merge(h0, (; newstep = () -> (TX.tdc_newstep!(); h0.newstep()), reset = () -> (TX.tdc_reset!(); h0.reset()),     # (M3d: w first, the species' local initial guess uses it; M3e: merge keeps MGROUPS' `predicted`)
                         update = (rhs, a_, b_) -> (h0.update(rhs, a_, b_); TX.tdc_update!(rhs, a_, b_)),
                         dmax = (rhs, a_, b_) -> max(h0.dmax(rhs, a_, b_), TX.tdc_dmax!(rhs, a_, b_))))
        end
        ax_ = CKAUX[]; src_ = "w = 0"
        if ax_ !== nothing && hasproperty(ax_, :tdc) && length(ax_.tdc.WTN) == length(TX.WTN[]) && ax_.tdc.a == Int(TF.getscalar(:ngrs)) && ax_.tdc.b == Int(TF.getscalar(:ngre))
            copyto!(TX.WTN[], ax_.tdc.WTN); copyto!(TX.WTO[], ax_.tdc.WTO); src_ = "the checkpoint's w"
        end
        pr("TDC: turbulent velocity w = unknown %d (neqn %d, meqn %d), A_s/A_d/A_c %g/%g/%g, alpha_p %g, alpha_t %g, alpha_D %g, face weight %g; start from %s\n", nq_, nq_, TF.PARAMS[:meqn], tp_[1], tp_[2], tp_[3], tp_[4], tp_[5], tp_[6], TX.TDC_CEN[], src_)
    end
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
        # TitanResidFast.stage_diag: the conservative total-energy rows of all zones at the converged state of the stage, summed × the stage's
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
            haskey(ENV, "ETOTLOG") && open(ENV["ETOTLOG"], "a") do io_       # ETOTLOG=<file> (diagnostic): one line per converged stage, all attempts: jstep kind dt imb thr G dE fin fout (erg/sr; TDC: E_t included, TX.etot_monitor)
                @printf(io_, "%d %s %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Int(TF.getscalar(:jstep)), kind, dt, m_.imb, m_.thr, m_.G, m_.dE, m_.fin, m_.fout)
            end
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
        @printf("DTERR: error-controlled TR-BDF2 (rtol %g, conv %g, x floor %g, F floor %g cE, E floor %g aT^4, %s norm, PI %s, safety/facmin/facmax %g/%g/%g, Newton-failure factor %g, BE retry %s); TITAN's heuristic and caps are diagnostic only\n",
                Trbdf2.EOPT.rtol, TF.getscalar(:conv), Trbdf2.EOPT.xfl, Trbdf2.EOPT.ffl, Trbdf2.EOPT.efl, Trbdf2.EOPT.maxnorm ? "max" : "rms", Trbdf2.EOPT.pi, Trbdf2.EOPT.safety, Trbdf2.EOPT.facmin, Trbdf2.EOPT.facmax, Trbdf2.EOPT.fail, Trbdf2.EOPT.be)
        let cv = TF.getscalar(:conv), rt = Trbdf2.EOPT.rtol, nk = get(ENV, "NKRYLOV", "0") != "0"
            (cv > rt || cv < 1e-3 * rt) && error(@sprintf("inconsistent tolerances: Newton CONV = %g vs DTERR rtol = %g (require 1e-3 * rtol <= CONV <= rtol)", cv, rt))
            nk && cv < 0.1 * NewtonKrylov.OPT.eta && error(@sprintf("inconsistent tolerances: Newton CONV = %g is tighter than 0.1 * NKRYLOV eta = %g", cv, NewtonKrylov.OPT.eta))
        end
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
    RUNEND[] = TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = output_step)
catch e
    pr("stopped: %s\n", first(sprint(showerror, e), 400))
    LASTOK[] === nothing || (TD.restore!(LASTOK[]); pr("state restored to the last accepted step %d for the checkpoint\n", KLAST[]))
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
pr("ALI line transfer: %d calls\n", NALICALL[])
TD.STEPFN[] === Trbdf2.step! && pr("integrator TR-BDF2, statistics %s\n", string(Trbdf2.STATS))
TD.KRYLOV[] === nothing || pr("%s\n", NewtonKrylov.summary_text())
Trbdf2.ERRCTL[] && pr("DTERR statistics %s; limiting variable (accepted steps) %s; TR-BDF2 %s\n", string(Trbdf2.ESTATS), string(Trbdf2.LIMHIST), string(Trbdf2.STATS))
Trbdf2.ERRCTL[] && pr("Newton true-residual test (converged only if max row |R|/(largest J_b entry of the row) <= CONV): %d convergences refused, largest refused residual %.3e\n", Trbdf2.TRSTATS[:rejected], Trbdf2.TRSTATS[:worst])
close(io)
