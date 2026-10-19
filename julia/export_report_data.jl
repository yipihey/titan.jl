# Data files for the report's interactive figures (report/app/data/*.json).
#   julia --project=env export_report_data.jl bench      benchmark profiles: Fortran 1995-code final dumps vs pure-Julia runs
#   julia --project=env export_report_data.jl qstar NAME[+CONT] [PROFILE.csv] [OUT]   a quasi-star run (lib/qs_runs/NAME_relax_snaps.jls,
#                                                     with the snapshots of its continuation CONT appended)
#   Z0_TABLE=ext julia --project=env export_report_data.jl spectra NAME...   emergent spectra lib/spectra/NAME.jls (spectra_run.jl)
#                                                     → spectra.json, with band luminosities
#   Z0_TABLE=ext julia --project=env export_report_data.jl analytic   simple models vs the W30 non-LTE solution → analytic.json
#   julia --project=env export_report_data.jl tests   Monte Carlo, Balmer series, axisymmetric, He I results → tests.json
using Serialization, Printf
include(joinpath(@__DIR__, "Metrics.jl")); using .Metrics
const OUT = joinpath(@__DIR__, "..", "report", "app", "data"); mkpath(OUT)

js(x::Real) = isfinite(x) ? @sprintf("%.6g", x) : "null"
js(x::Integer) = string(x)
js(x::AbstractString) = "\"" * x * "\""
js(v::AbstractVector) = "[" * join(js.(v), ",") * "]"
js(d::NamedTuple) = "{" * join(["\"$k\":" * js(v) for (k, v) in pairs(d)], ",") * "}"
writejs(name, x) = (open(io -> print(io, js(x)), joinpath(OUT, name), "w"); println("wrote ", joinpath(OUT, name)))

"Fortran final dump CSV → named columns (zone values on rows 1:n-1, interface values on 1:n)."
function readdump(path)
    lines = filter(l -> !startswith(l, "#"), readlines(path))
    hdr = split(lines[1], ","); rows = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    col(n) = [r[findfirst(==(n), hdr)] for r in rows]
    t = parse(Float64, match(r"timen=([-0-9.eE+]+)", readline(path))[1])
    (; t, r = col("r"), u = col("u"), rho = col("rho")[1:end-1], T = col("T")[1:end-1], Trad = col("Trad")[1:end-1])
end

mode = ARGS[1]
if mode == "bench"
    for (name, nsub) in (("SOD", 30), ("SUBC", 40), ("SUPC", 40))
        F = readdump(joinpath(@__DIR__, "..", "reference_runs", name, "$(name)_final.csv"))
        local snaps = deserialize(joinpath(@__DIR__, "..", "lib", "p2_runs", "$(name)_pure", "snaps.jls"))
        J = snaps[end]
        idx = unique(round.(Int, range(1, length(snaps), length = nsub)))
        series = [(t = s.time, r = s.r, rho = s.rho, T = s.T, Trad = Metrics.Trad(s), u = s.u) for s in snaps[idx]]
        writejs("bench_$(name).json", (name = name, fortran = F,
            julia = (t = J.time, r = J.r, u = J.u, rho = J.rho, T = J.T, Trad = Metrics.Trad(J)),
            nsteps_julia = length(snaps), series = series))
    end
elseif mode == "spectra"
    include(joinpath(@__DIR__, "Z0Physics.jl"))
    names = ARGS[2:end]
    LABEL = Dict("star_t0" => "steady star before infall (t = 0)",
                 "A220_late_ext" => "infall, 4.2e8 s, LTE infall + free-fall extension",
                 "W30_late_ext" => "wide domain W30, 4.2e8 s, LTE infall (TITAN to 30 R) + extension",
                 "A220_late_ion" => "infall, 4.2e8 s, infall fully ionized (limit)",
                 "A220_late_noext" => "infall, 4.2e8 s, grid only (to 1.20 R)",
                 "A220_late_ext_nlte" => "infall, 4.2e8 s, non-LTE hydrogen in the infall",
                 "W30_late_ext_nlte" => "wide domain W30, 4.2e8 s, non-LTE hydrogen in the infall",
                 "A220_late_ext10_nlte" => "infall, 4.2e8 s, non-LTE, infall ending at 12 R (test)",
                 "A220_late_bare" => "infall, 4.2e8 s, star only (infall above the shock removed)",
                 "A220_late_hotion" => "infall, 4.2e8 s, cooling layer fully ionized (test)",
                 "R2_late_ext" => "Eulerian run R2, 5.8e7 s, LTE infall + extension",
                 "A220_corona_ext" => "early ramp, 7.9e6 s (Ṁ = 0.003 M⊙/yr at the shock, 0.14 at the grid edge), LTE + extension")
    hP, kB, cc = 6.62607015e-27, 1.380649e-16, 2.99792458e10; keV = 1.602176634e-9
    band(ν, L, ν1, ν2) = (s = 0.0; for i in 1:length(ν)-1
        a, b = max(ν[i], ν1), min(ν[i+1], ν2); b > a || continue
        f(x) = (L[i] > 0 && L[i+1] > 0) ? L[i] * (L[i+1] / L[i])^((log(x) - log(ν[i])) / (log(ν[i+1]) - log(ν[i]))) : L[i] + (L[i+1] - L[i]) * (x - ν[i]) / (ν[i+1] - ν[i]); s += 0.5 * (f(a) + f(b)) * (b - a); end; s)   # power-law interpolation (edges)
    Lat(ν, L, λÅ) = (x = cc / (λÅ * 1e-8); i = argmin(abs.(ν .- x)); L[i])
    runs = []
    for n in names
        d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$n.jls"))
        ν, L = d.ν, d.L
        Ltot = band(ν, L, 0.0, Inf)
        bands = (X_05_2 = band(ν, L, 0.5keV / hP, 2keV / hP), X_2_10 = band(ν, L, 2keV / hP, 10keV / hP),
                 X_10_40 = band(ν, L, 10keV / hP, 40keV / hP), EUV = band(ν, L, 13.6 * 1.602176634e-12 / hP, 0.1keV / hP),
                 FUV = band(ν, L, cc / 2000e-8, cc / 912e-8), opt = band(ν, L, cc / 7000e-8, cc / 3646e-8), total = Ltot)
        Dbalmer = log10(Lat(ν, L, 3664.0) / Lat(ν, L, 3628.0))       # Balmer jump: L_ν just red / just blue of 3646 Å
        R0 = 8.4352e15
        prof = (r = d.atm.r ./ R0, T = d.atm.T, T_titan = d.T_titan, rho = d.atm.ρ,
                xH = haskey(d, :xH) ? d.xH : [Z0Physics.saha(d.atm.ρ[k], d.atm.T[k]).xH for k in eachindex(d.atm.r)])
        tab = Z0Physics.table(); τ = 0.0; kph = length(d.atm.r)          # TITAN-style photosphere: Rosseland τ = 2/3
        for k in length(d.atm.r):-1:2
            τ += Z0Physics.kap(tab, d.atm.ρ[k], d.T_titan[k]) * d.atm.ρ[k] * (d.atm.r[k] - d.atm.r[k-1])
            τ >= 2 / 3 && (kph = k; break)
        end
        rph = d.atm.r[kph]; Teff = (d.Ltitan / (4π * 5.670374419e-5 * rph^2))^0.25
        # non-LTE solutions: the full line spectrum of nlte_spectrum.jl (lib/spectra/linespec_OUT.jls, OUT = n minus "_nlte")
        lsf = joinpath(@__DIR__, "..", "lib", "spectra", "linespec_$(replace(n, "_nlte" => "")).jls")
        linespec = if endswith(n, "_nlte") && isfile(lsf)
            ls = deserialize(lsf)
            (lam = ls.λ, nuLnu = ls.νLν, table = [(name = x.name, lam = x.λ, EW = x.EW, L = x.Lline) for x in ls.lines])
        else
            (lam = Float64[], nuLnu = Float64[], table = [])
        end
        lf = joinpath(@__DIR__, "..", "lib", "spectra", "lines_$n.jls")
        lines = isfile(lf) ? deserialize(lf) : []
        push!(runs, (name = n, lines = lines, linespec = linespec, Lline = haskey(d, :Lline) ? d.Lline : 0.0, label = get(LABEL, n, n), run = d.run, t = d.t, Mdot = d.Mdot * Z0Physics.yr / Z0Physics.Msun,
            Ltitan = d.Ltitan, rph = rph / R0, Teff = Teff, lam = cc ./ ν .* 1e8, E = hP .* ν ./ keV, nuLnu = ν .* L, bands = bands, Dbalmer = Dbalmer,
            nlucy = d.nlucy, ion = haskey(d, :ion) ? d.ion : false, prof = prof))
        @printf("%-18s ∫L = %.4e (TITAN %.4e, ratio %.4f)  X(0.5–2) %.2e  X(2–10) %.2e  X(10–40) %.2e  EUV %.2e  FUV %.3e  opt %.3e  D_B = %.3f  r_ph/R = %.5f  T_eff = %.0f K\n",
                n, Ltot, d.Ltitan, Ltot / d.Ltitan, bands.X_05_2, bands.X_2_10, bands.X_10_40, bands.EUV, bands.FUV, bands.opt, Dbalmer, rph / R0, Teff)
    end
    writejs("spectra.json", (R = 8.4352e15, runs = runs))
elseif mode == "qstar"
    name = ARGS[2]                     # "A+B": run B continued A (QS_FROM); their snapshots are joined
    parts = split(name, "+"); name = String(parts[1])
    snaps = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(parts[1])_relax_snaps.jls"))
    for p in parts[2:end]
        append!(snaps, deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(p)_relax_snaps.jls"))[2:end])
    end
    snaps = filter(s -> all(isfinite, s.r) && all(isfinite, s.T), snaps)
    prof = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, "..", "reference_runs", "qs_profile_f0.1_md0.csv")
    lines = readlines(prof); hdr = Dict(Symbol(split(kv, "=")[1]) => parse(Float64, split(kv, "=")[2]) for kv in split(lstrip(lines[1], ['#', ' '])))
    C = [parse.(Float64, split(l, ",")) for l in lines[3:end]]
    names = split(lines[2], ","); ci(n) = findfirst(==(n), names)
    steady = (r = [c[ci("r")] for c in C], T = [c[ci("T")] for c in C], rho = [c[ci("rho")] for c in C],
              Lconv = [max(0.0, 1 - c[ci("nabla")] / c[ci("nabla_rad")]) for c in C])
    out = length(ARGS) >= 4 ? ARGS[4] : "qstar_$(name).json"
    # time history of every snapshot (shock, photosphere, accretion rate; qstar_summary.jl), profiles of at most 150
    hist = nothing
    if length(snaps) > 2 && length(snaps[end].u) > 0 && any(s -> minimum(s.u) < -1e7, snaps)      # infall runs only
        include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "qstar_summary.jl"))
        M0 = Base.invokelatest(QSum.mass, snaps[1])
        H = [Base.invokelatest(QSum.summary, s, hdr[:R], hdr[:L], M0) for s in snaps]
        hist = (t = [h.t for h in H], rsh = [h.rsh for h in H], Tpk = [h.Tpk for h in H], rph = [h.rph for h in H],
                Teff = [h.Teff for h in H], Lph = [h.Lph for h in H], Mdot_out = [h.Mdot for h in H], Mdot_sh = [h.Mdsh for h in H])
    end
    idx = unique(round.(Int, range(1, length(snaps), length = min(length(snaps), 150))))
    prof = (name = name, R = hdr[:R], L = hdr[:L], M = hdr[:M], steady = steady,
        snaps = [(t = s.t, r = s.r, u = s.u, rho = s.d, T = s.T, pg = s.pg, Er = s.er, Lrad = s.Lrad ./ hdr[:L], Lconv = s.Lconv ./ hdr[:L])
                 for s in snaps[idx]])
    writejs(out, hist === nothing ? prof : merge(prof, (hist = hist,)))
elseif mode == "tests"
    # the decisive tests of the second book (chapter "Testing the comparison") → tests.json: Monte Carlo line profiles
    # (mc_escatter.jl), the Balmer-limit region (balmer_series.jl), axisymmetric profiles (axisym_lines.jl), He I λ10830
    # (he10830.jl); missing inputs are skipped
    sp(n) = joinpath(@__DIR__, "..", "lib", "spectra", n)
    mc = []
    for (tag, label) in (("Halpha_es0_res1", "Hα, no electron scattering"), ("Halpha_es1_res1", "Hα, electron scattering, infall"),
                         ("Halpha_es1_res1_Te1e4", "Hα, electrons at 10⁴ K"), ("Halpha_es1_res1_Te3e4", "Hα, electrons at 3×10⁴ K"),
                         ("Halpha_es1_res1_th90", "Hα, electron scattering, pure outflow"), ("Halpha_es1_res1_th45", "Hα, infall + outflow cones θ_o = 45°"),
                         ("Halpha_es1_res1_th60", "Hα, infall + outflow cones θ_o = 60°"), ("Hbeta_es0_res1", "Hβ, no electron scattering"),
                         ("Hbeta_es1_res1", "Hβ, electron scattering, infall"))
        f = sp("mc_W30_late_ext_$tag.jls"); isfile(f) || continue
        m = deserialize(f)
        push!(mc, (tag = tag, label = label, V = m.V, prof = m.prof, profi = haskey(m, :profi) && m.θo > 0 ? m.profi : [],
                   fwhm = m.fwhm, Wb = m.Wb, Wr = m.Wr, centroid = m.centroid, stats = m.stats))
    end
    bal = []
    for tag in ("n40_b12", "n40_lte", "n80_b12")
        f = sp("balmer_W30_late_ext_$tag.jls"); isfile(f) || continue
        b = deserialize(f); k = 1:3:length(b.λ)
        push!(bal, (tag = tag, lam = b.λ[k], Lnu = b.Lnu[k], Lc = b.Lc[k], L100 = b.L100[k], L30 = b.L30[k], brk = b.brk))
    end
    ax = isfile(sp("axisym_W30_late_ext.jls")) ? [(line = x.line, theta = x.θo, incl = x.incl, v = x.v, F = x.F, EW = x.EW) for x in deserialize(sp("axisym_W30_late_ext.jls"))] : []
    he = []
    for n in ("W30_late_ext", "C30_ext")
        f = sp("he10830_$(n)_s31.0_x1.0.jls"); isfile(f) || continue
        h = deserialize(f); push!(he, (name = n, v = h.v, F = h.F, EW = h.EW, N3 = h.N3))
    end
    writejs("tests.json", (mc = mc, balmer = bal, axisym = ax, he = he))
elseif mode == "analytic"
    # the simple models of the infall's hydrogen (analytic_infall.jl, closed_form.jl) against the full non-LTE solution
    # of W30 → analytic.json: profiles for the second book's "analytic" panel, and the constants of the closed form
    include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
    Z = Z0Physics; NHm = NLTEH; R0 = 8.4352e15; M = 1e6 * Z.Msun; sp(n) = joinpath(@__DIR__, "..", "lib", "spectra", n)
    dn = deserialize(sp("nlte_W30_late_ext.jls"))
    ν32 = (NHm.En[3] - NHm.En[2]) / Z.h; B0 = 2Z.h * ν32^3 / Z.c^2; Bh = Spectrum.Bν(ν32, 7310.0)
    fA(i, l, u) = (st = dn.nl[i]; f = NHm.fosc(l, u); ν = (NHm.En[u] - NHm.En[l]) / Z.h; Blu = 4π * NHm.πe2mc * f / (Z.h * ν);
        Bul = NHm.gn[l] / NHm.gn[u] * Blu; A = 2Z.h * ν^3 / Z.c^2 * Bul; Z.h * ν * st.β[l, u] * (st.n[u] * A - (st.n[l] * Blu - st.n[u] * Bul) * dn.Jhα[i]))
    full = (r = dn.r ./ R0, T = dn.T, x = [st.np / (st.np + sum(st.n)) for st in dn.nl],
            SB = [B0 / (st.n[2] * NHm.gn[3] / (st.n[3] * NHm.gn[2]) - 1) / Bh for st in dn.nl],
            dL = [4π * dn.r[i]^3 * fA(i, 2, 3) for i in eachindex(dn.nl)])       # dL(Hα)/d ln r
    sim(n) = (a = deserialize(sp("analytic_$n.jls")); (r = a.r ./ R0, x = a.xH, dL = 4π .* a.r .^ 3 .* a.PHa, LHa = a.LHa, NHI = a.NHI, τT = a.τT))
    cf = deserialize(sp("closedform_W30.jls"))
    writejs("analytic.json", (Teff = 7310.0, L = 1.4204e44, M = M, R0 = R0, Rst = sqrt(1.4204e44 / (4π * Z.σSB * 7310.0^4)),
        X = 0.76, mH = Z.mH, Msun = Z.Msun, yr = Z.yr, G = Z.G, p2 = cf.p2, p3 = cf.p3, G2 = cf.G2, G3 = cf.G3,
        full = full, closed = (r = cf.r ./ R0, T = cf.T, x = cf.x),
        simple = sim("simple_W30T"), simple_bb = sim("simple_W30T_bb"),
        fuv = [(f = f, s = sim("simple_fuv$f")) for f in ("0.1", "0.01", "0.001")]))
end
