# Reader and exporter for runs recorded with RUNREC=<dir> (popiii_m2c.jl).
#
#   julia --project=<env> runrec.jl <dir> [every | budget=<bytes>] [skip=<field>,<field>...]
#                                                      # writes <dir>/export/{frames.json, frames.f32, track.csv, nuc_v1.csv}
#   every: keep every n-th snapshot; budget=: the smallest stride for which frames.f32 stays below that many bytes;
#   skip=: per-zone fields to leave out of frames.f32.
#
# The run directory holds one snapshot stream per process segment (snaps_<k0>_<pid>.jls). A restart from an earlier checkpoint
# re-runs steps that an abandoned segment had already written, so the surviving history is found by reading all segments in
# start order and keeping, from the newest snapshot backwards, only snapshots strictly earlier in time than the one after them.
# A segment killed mid-write ends in a truncated record, which is dropped.
#
# frames.f32 is Float32, frame after frame, field after field, nz values each (nz varies between frames: rezoning, the remap
# onto the accretion mesh). frames.json has the field names and units, and per frame the step, time, nz, byte offset and the
# scalars of track.csv. Cell-centred quantities; values at interfaces (u, L, Ṁ, F/cE) are averaged onto the cells.
#
# Nuclear burning (runs with HOT=1 record the abundances of the 19 network species in every zone): the mass fractions of the
# nuclei are exported per zone, and the luminosity of each burning process is recomputed from the snapshot with NucNet (the
# five parts add up to the run's own L_nuc). `nuc_v1.csv` caches these for every snapshot of the lineage; the onset marks in
# frames.json ("nuc") are found on that full table, not on the exported subset.
using Serialization, Printf

const Rs = 6.957e10; const Ms = 1.989e33; const Ls = 3.83e33; const yr = 3.156e7; const G = 6.674e-8
const cl = 2.998e10; const arad = 7.5657e-15; const sb = 5.6704e-5; const mH = 1.6726e-24

# the nuclear network, for the burning diagnostics only: the reader and the rest of the export work without it
isdefined(@__MODULE__, :NucNet) || try include(joinpath(@__DIR__, "NucNet.jl")) catch e; @warn "NucNet not loaded: no burning diagnostics" exception = e end

function readstream(f)
    v = Any[]
    open(f) do io
        while !eof(io)
            try push!(v, deserialize(io)) catch; break end
        end
    end
    v
end

"All snapshots of the run's surviving branch, in time order."
function lineage(dir)
    # in write order (attempts restarted from the same step share k0; the newest file holds the surviving branch's end)
    fs = sort(filter(f -> occursin(r"^snaps_\d+_\d+\.jls$", f), readdir(dir)); by = f -> mtime(joinpath(dir, f)))   # RUNREC streams only (a run directory can also hold an old-format SNAPS file)
    all = Any[]; for f in fs; append!(all, readstream(joinpath(dir, f))); end
    isempty(all) && return all
    K = Any[all[end]]
    for x in reverse(all[1:end-1]); x.t < K[end].t && push!(K, x); end
    reverse(K)
end

"Old-format history (SNAPS=1 file, e.g. snaps_abn02z19_lines_n200_fx3.jls) as RUNREC-style snapshots, surviving branch only."
function lineage_old(file)
    A = deserialize(file); K = Any[A[end]]
    for x in reverse(A[1:end-1]); x.t < K[end].t && push!(K, x); end
    K = reverse(K)
    map(K) do x
        n = length(x.ρ)
        Y = vcat([y[1] for y in x.y]', [y[2] for y in x.y]', [y[3] for y in x.y]', [y[4] for y in x.y]')
        (; step = x.step, t = x.t, dt = NaN, r = x.r, u = x.u, m = x.m, ρ = x.ρ, T = x.T, e = x.e, pg = x.pg, er = x.er, fr = x.fr,
           chi = x.chi, kp = fill(NaN, n), fedd = x.fedd, qf = x.qf, eq = x.eq, Λext = x.Λext, spn = [:Hp, :Hm, :H2, :H2p], Y,
           gasq = Float64[], hotY = nothing, Lnuc = NaN, Llines = x.L_lines, nchemfb = 0, kph = 0, Teff = NaN)
    end
end

sp(x, s) = (i = findfirst(==(s), x.spn); i === nothing ? zeros(size(x.Y, 2)) : x.Y[i, :])
avg(v) = 0.5 .* (v[1:end-1] .+ v[2:end])

const FIELDS = [("logr", "log r [R_sun] (cell centre)"), ("logm", "log m [M_sun] (enclosed, cell centre)"),
    ("logrho", "log rho [g cm^-3]"), ("logT", "log T_gas [K]"), ("logTrad", "log T_rad = (E/a)^1/4 [K]"),
    ("u", "u [km/s]"), ("L", "L [L_sun] = 4 pi r^2 F"), ("xHp", "x(H+)"), ("xH2", "x(H2) (per H nucleus)"),
    ("logxHm", "log x(H-)"), ("logtau", "log tau_R: grey continuum optical depth (Rosseland-mean extinction, absorption + scattering), integrated inwards from the outer mesh boundary"), ("Gamma", "kappa L / 4 pi G m c"),
    ("FcE", "F / cE"), ("logkap", "log chi [cm^2 g^-1]"), ("logkapP", "log kappa_P [cm^2 g^-1]"),
    ("Mdot", "Mdot = -4 pi r^2 rho u [M_sun/yr]"), ("eq", "equilibrium (table) zone flag"), ("logLam", "log |line cooling| (Lambda_ext) [erg cm^-3 s^-1]"),
    ("logr_o", "log r of the zone's outer interface [R_sun]"), ("logm_o", "log m at the zone's outer interface [M_sun]"),
    ("logXH", "log mass fraction of 1H nuclei"), ("logXD", "log mass fraction of D"), ("logXHe3", "log mass fraction of 3He"),
    ("logXHe4", "log mass fraction of 4He"), ("logXC", "log mass fraction of C (12C + 13C)"),
    ("logXN", "log mass fraction of N (13N + 14N + 15N)"), ("logXO", "log mass fraction of O (15O + 16O + 17O + 18O)")]
const SCAL = ["step", "t_yr", "dt_s", "logn_c", "T_c", "M_star", "R_star", "r_ph", "T_eff_lum", "T_eff_sl", "L_ph", "Mdot_1e4", "u_min", "L_nuc", "L_lines", "nz", "nchemfb",
    "L_D", "L_pp", "L_He3", "L_CNO", "L_3a", "L_ph_med"]      # luminosity by burning process [L_sun], see nucdiag; L_ph_med: see nucref

# CONVDIAG=1 runs: snapshots carry `cd` (popiii_m2c convdiag_snap). Exported as cv_* fields (cell-centred; interface quantities averaged onto the cells, NaN-aware) and cv_* scalars,
# only if some snapshot has `cd`; a frame without it gets NaN. Runs without CONVDIAG export exactly the fields above.
const CVFIELDS = [("cv_nab", "nabla = dlnT/dlnP used by MLT (after the MLT_HS blend)"), ("cv_nab_raw", "nabla before the MLT_HS blend"), ("cv_hs_w", "MLT_HS blend weight (1 = untouched)"),
    ("cv_nad_used", "nabla_ad that material_thermo gave MLT (1e30 = network zone, gate closed)"), ("cv_nad_eq", "nabla_ad, equilibrium chemistry (table EOS)"), ("cv_nad_fr", "nabla_ad, frozen chemistry (zone's network EOS)"),
    ("cv_nab_rad", "nabla_rad from the radiative flux present"), ("cv_nab_mu", "dln mu/dlnP"), ("cv_nab_L", "Ledoux nabla_ad,fr + phi/delta nabla_mu"),
    ("cv_N2_S", "N^2 Schwarzschild, frozen [s^-2] (subsonic zones only meaningful)"), ("cv_N2_eq", "N^2 Schwarzschild, equilibrium chemistry [s^-2]"), ("cv_N2_L", "N^2 Ledoux, frozen [s^-2]"),
    ("cv_Fconv_Ftot", "F_conv/(F_rad+F_conv), as used"), ("cv_Fpot_Ftot", "gate-free MLT flux / (F_rad + F_conv)"), ("cv_vconv", "MLT velocity as used [km/s]"), ("cv_lmix", "mixing length [R_sun]"),
    ("cv_w_tau", "MLT optical-depth weight"), ("cv_w_inf", "MLT_NOINFALL factor"), ("cv_w_th", "thermodynamic gate (1: nabla_ad,used < 1e29)"), ("cv_mach", "|u|/c_s (zone), c_s = sqrt(5/3 p_gas/rho)"),
    ("cv_branch", "zone branch: 0 network, 1 coupled network, 2 table, 3 hot"), ("cv_logt_settle", "log10 t_settle = H_P/|u| [s]"), ("cv_logt_N", "log10 1/sqrt|N^2_L| [s]"),
    ("cv_logt_turn", "log10 l/v [s]"), ("cv_logt_chem", "log10 chemical time of the coupled species [s] (NaN: not coupled)"),
    ("cv_sup", "MLTBUOY: blob superadiabaticity sup = A/(delta_b |dlnP_hs|) (NaN without MLTBUOY)"), ("cv_gam1b", "M6: Re Gamma_1,b(omega) of the blob (network linear response at omega = 1/t_turn)"), ("cv_gam1_0", "M6: Gamma_1(omega -> 0), the network equilibrium"), ("cv_gam1_inf", "M6: Gamma_1(omega -> inf), frozen composition"), ("cv_wtau", "M6: omega tau_max (tau_max the slowest relaxation time of J = dS/dx)"),
    ("cv_tchem", "M6: tau_max = slowest relaxation time 1/|lambda| of J = dS/dx at the interface [s] (mean of the two zones; 0: no coupled network)"),
    ("cv_e_stor", "gas-energy row: storage [erg/s/g]"), ("cv_e_pdv", "pdV work"), ("cv_e_av", "artificial viscosity work"), ("cv_e_kin", "kinetic-energy correction"),
    ("cv_e_adv", "advection"), ("cv_e_dfe", "energy diffusion"), ("cv_e_rad", "gas-radiation exchange"), ("cv_e_conv", "convective-flux divergence"),
    ("cv_e_chem", "chemistry cooling/heating"), ("cv_e_src", "nuclear source"), ("cv_e_res", "row residual"),
    ("cv_w", "TDC: turbulent velocity w at the interface [km/s] (NaN without TDC)"), ("cv_PtP", "TDC: turbulent pressure P_t/P of the zone"), ("cv_Dmix", "TDC: species mixing diffusivity D = alpha_D l w [cm^2/s]"),
    ("cv_Q", "TDC: realised exchange gas -> turbulent energy Q of the interface CV [erg/s/sr]")]
const CVSCAL = ["cv_rin", "cv_rout", "cv_mlayer", "cv_N2td2", "cv_FconvFtot_max", "cv_FpotFtot_max", "cv_nlayer", "cv_Et", "cv_PtP_ps", "cv_PtP_ps3"]    # the sub-shock unstable layer (convdiag_snap), r in R_star, m in M_sun
hascd(x) = hasproperty(x, :cd) && x.cd !== nothing
nanavg(v) = [(a_ = v[i]; b_ = v[i+1]; isnan(a_) ? b_ : isnan(b_) ? a_ : 0.5 * (a_ + b_)) for i in 1:length(v)-1]

function cvframe(x)
    c = x.cd; n = length(x.ρ); dm = diff(x.m); l10(v) = [isfinite(q) && q > 0 ? log10(q) : NaN for q in v]
    zg(v) = 4π .* v ./ dm      # per steradian -> per gram
    F = Dict("cv_nab" => nanavg(c.nab), "cv_nab_raw" => nanavg(c.nab_raw), "cv_hs_w" => nanavg(c.hs_w), "cv_nad_used" => nanavg(c.nad_used), "cv_nad_eq" => nanavg(c.nad_eq),
        "cv_nad_fr" => nanavg(c.nad_fr), "cv_nab_rad" => nanavg(c.nab_rad), "cv_nab_mu" => nanavg(c.nab_mu), "cv_nab_L" => nanavg(c.nab_L), "cv_N2_S" => nanavg(c.N2_S),
        "cv_N2_eq" => nanavg(c.N2_eq), "cv_N2_L" => nanavg(c.N2_L), "cv_Fconv_Ftot" => nanavg(c.Fconv_over_Ftot), "cv_Fpot_Ftot" => nanavg(c.Fpot_over_Ftot),
        "cv_vconv" => nanavg(c.vconv) ./ 1e5, "cv_lmix" => nanavg(c.lmix) ./ Rs, "cv_w_tau" => nanavg(c.w_tau), "cv_w_inf" => nanavg(c.w_inf), "cv_w_th" => nanavg(c.w_th),
        "cv_mach" => c.mach, "cv_branch" => c.branch, "cv_logt_settle" => l10(nanavg(c.t_settle)), "cv_logt_N" => l10(nanavg(c.t_N)), "cv_logt_turn" => l10(nanavg(c.t_turn)),
        "cv_logt_chem" => l10(c.t_chem))
    for (nm, f) in (("cv_sup", :sup), ("cv_gam1b", :gam1b), ("cv_gam1_0", :gam1_0), ("cv_gam1_inf", :gam1_inf), ("cv_wtau", :wtau), ("cv_tchem", :tchem))     # (records of older runs lack the M6 fields)
        F[nm] = hasproperty(c, f) ? nanavg(getproperty(c, f)) : fill(NaN, n)
    end
    for t in (:stor, :pdv, :av, :kin, :adv, :dfe, :rad, :conv, :chem, :src, :res); F["cv_e_$t"] = zg(getproperty(c, Symbol("e_", t))); end
    if hasproperty(c, :tdc_w)
        F["cv_w"] = nanavg(c.tdc_w) ./ 1e5; F["cv_PtP"] = c.tdc_PtP; F["cv_Dmix"] = nanavg(c.tdc_Dmix); F["cv_Q"] = nanavg(c.tdc_Q)
    else
        for nm in ("cv_w", "cv_PtP", "cv_Dmix", "cv_Q"); F[nm] = fill(NaN, n); end
    end
    F, vcat(Float64.(c.layer[1:7]), hasproperty(c, :tdc_scal) ? Float64.(c.tdc_scal) : fill(NaN, 3))
end

const NUCSPEC = [("logXH", (:p,)), ("logXD", (:d,)), ("logXHe3", (:he3,)), ("logXHe4", (:he4,)), ("logXC", (:c12, :c13)),
    ("logXN", (:n13, :n14, :n15)), ("logXO", (:o15, :o16, :o17, :o18))]
hasnuc(x) = isdefined(@__MODULE__, :NucNet) && x.hotY !== nothing && size(x.hotY, 2) == length(x.ρ)

"""
Luminosity [L_sun] of each burning process in snapshot `x`, from its temperatures, densities and nuclear abundances (zones
above 3e5 K, the run's HOTT1). Net of neutrino losses, like the run's own L_nuc, which the five add up to:
  L_D    burning of deuterium that was not made by p+p in the same zone: d(p,γ)³He beyond the p+p rate, plus d+d and d+³He;
  L_pp   hydrogen burning by the pp chain up to ³He: p+p and pep, and the d(p,γ)³He that follows them;
  L_He3  ³He+³He, ³He+⁴He and the ⁷Be, ⁷Li, ⁸B branches (burning of ³He, whether accreted, made from deuterium or by p+p);
  L_CNO  the CNO cycles;  L_3a  helium burning (3α and α captures).
"""
function nucdiag(x)
    hasnuc(x) || return fill(NaN, 5)
    NN = NucNet; n = length(x.ρ); Y = x.hotY; dm = diff(x.m); c = NN.NA * NN.MeV
    idpg = findfirst(r -> r.name == "d(p,g)he3", NN.NET); ipp = findall(r -> r.name in ("p(p,e+nu)d", "p(p e-,nu)d (pep)"), NN.NET)
    LD = Lpp = L3 = Lc = La = 0.0
    for k in 1:n
        x.T[k] > 3e5 || continue
        Yk = Y[:, k]; kr = NN.rates_k(x.T[k], x.ρ[k], NN.Plasma(x.ρ[k], Yk), NN.ye_of(Yk)); _, _, g = NN.energy(Yk, kr)
        evp = sum(NN.event(NN.NET[i], kr[i], Yk) for i in ipp); evd = NN.event(NN.NET[idpg], kr[idpg], Yk); qd = c * NN.NET[idpg].Q
        epp = sum(c * (NN.NET[i].Q - NN.NET[i].Enu) * NN.event(NN.NET[i], kr[i], Yk) for i in ipp) + qd * min(evp, evd)
        ed = qd * max(evd - evp, 0.0) + g.d
        LD += dm[k] * ed; Lpp += dm[k] * epp; L3 += dm[k] * (g.pp + g.d - epp - ed); Lc += dm[k] * g.cno; La += dm[k] * g.he
    end
    [LD, Lpp, L3, Lc, La] ./ Ls
end

function frame(x, nuc = fill(NaN, 6), cvon = false)
    n = length(x.ρ); r = x.r; m = x.m; rc = sqrt.(r[1:n] .* r[2:n+1] .+ (r[1:n] .== 0) .* r[2:n+1] .^ 2 ./ 4)
    rc[1] = r[1] == 0 ? 0.5r[2] : rc[1]
    mc = avg(m); L = 4π .* r .^ 2 .* x.fr ./ Ls; Lc = avg(L)
    τ = zeros(n + 1); for k in n:-1:1; τ[k] = τ[k+1] + x.chi[k] * x.ρ[k] * (r[k+1] - r[k]); end
    τc = avg(τ)
    Γ = [x.chi[k] * Lc[k] * Ls / (4π * G * max(mc[k], 1e-30) * cl) for k in 1:n]
    fce = avg(x.fr) ./ (cl .* x.er)
    md = -4π .* r .^ 2 .* [x.ρ[min(max(k, 1), n)] for k in 1:n+1] .* x.u .* yr ./ Ms
    l10(v) = [x > 0 ? log10(x) : NaN for x in v]      # non-positive (no radiation field yet, zero line cooling): NaN
    F = Dict("logr" => l10(rc ./ Rs), "logm" => l10(mc ./ Ms), "logrho" => l10(x.ρ), "logT" => l10(x.T),
        "logTrad" => l10((max.(x.er, 0) ./ arad) .^ 0.25), "u" => avg(x.u) ./ 1e5, "L" => Lc, "xHp" => sp(x, :Hp), "xH2" => sp(x, :H2),
        "logxHm" => l10(sp(x, :Hm)), "logtau" => l10(τc), "Gamma" => Γ, "FcE" => fce, "logkap" => l10(x.chi), "logkapP" => l10(x.kp),
        "Mdot" => avg(md), "eq" => Float64.(x.eq), "logLam" => l10(abs.(x.Λext)),
        "logr_o" => l10(r[2:n+1] ./ Rs), "logm_o" => l10(m[2:n+1] ./ Ms))
    for (key, sps) in NUCSPEC      # mass fractions of the nuclei, X = A Y; a nucleus that is absent (Y = 0): NaN
        F[key] = hasnuc(x) ? l10([sum(NucNet.AA[NucNet.IDX[q]] * x.hotY[NucNet.IDX[q], k] for q in sps) for k in 1:n]) : fill(NaN, n)
    end
    # scalars
    js = findfirst(<(-2e5), x.u); Mst = js === nothing ? NaN : m[js] / Ms; Rst = js === nothing ? NaN : r[js] / Rs
    kph = findlast(>(2 / 3), τ); rph = kph === nothing ? NaN : r[kph] / Rs
    Lph = kph === nothing ? NaN : L[kph]; Tl = (kph === nothing || Lph <= 0) ? NaN : (Lph * Ls / (4π * r[kph]^2 * sb))^0.25
    j4 = findfirst(>=(1e4 * Rs), r); md4 = j4 === nothing ? NaN : md[j4]
    S = [x.step, x.t / yr, x.dt, log10(x.ρ[1] / mH), x.T[1], Mst, Rst, rph, Tl, x.Teff, Lph, md4, minimum(x.u) / 1e5,
         x.Lnuc, x.Llines, n, x.nchemfb, nuc...]
    if cvon
        if hascd(x); Fc, Sc = cvframe(x); merge!(F, Fc); append!(S, Sc)
        else append!(S, fill(NaN, length(CVSCAL))) end
    end
    F, S
end

const NUCTHR = 1e-3      # a process counts as burning once it supplies this fraction of the luminosity at the photosphere
const NUCWIN = 100       # ... taken as the median of L_ph over this many snapshots to either side (the photosphere can pulsate)
"""
Burning table of every snapshot of the lineage: (step, t) => [L_D, L_pp, L_He3, L_CNO, L_3a, L_ph, M_star]. Cached in
<out>/nuc_v1.csv, so a refresh computes only the snapshots written since the last one.
"""
function nuctable(L, out)
    f = joinpath(out, "nuc_v1.csv"); C = Dict{Tuple{Int,Float64},Vector{Float64}}()
    if isfile(f)
        for l in Iterators.drop(eachline(f), 1)
            v = split(l, ","); length(v) == 9 || continue
            C[(parse(Int, v[1]), parse(Float64, v[2]))] = parse.(Float64, v[3:9])
        end
    end
    T = Dict{Tuple{Int,Float64},Vector{Float64}}()
    for x in L
        T[(x.step, x.t)] = get(C, (x.step, x.t)) do
            S = frame(x)[2]; vcat(nucdiag(x), S[11], S[6])
        end
    end
    tmp = f * ".tmp"
    open(tmp, "w") do io
        println(io, "step,t_s,L_D,L_pp,L_He3,L_CNO,L_3a,L_ph,M_star")
        for x in L; println(io, x.step, ",", @sprintf("%.17g", x.t), ",", join((@sprintf("%.9g", v) for v in T[(x.step, x.t)]), ",")); end
    end
    mv(tmp, f; force = true)
    T
end

"Reference luminosity of every snapshot of the lineage: the running median of L_ph over ±NUCWIN snapshots."
function nucref(L, T)
    lp = [T[(x.step, x.t)][6] for x in L]; n = length(lp)
    map(1:n) do i
        s = sort(filter(isfinite, lp[max(1, i - NUCWIN):min(n, i + NUCWIN)])); isempty(s) ? NaN : s[(length(s) + 1) ÷ 2]
    end
end

"""
Onset of deuterium burning, of hydrogen burning by the pp chain and of the CNO cycle: the first snapshot from which the
process supplies more than NUCTHR of the reference luminosity `ref` (nucref) in three snapshots in a row. state "before":
that already holds in the first snapshot of the record; "none": not reached. `frac` is the fraction in the last snapshot.
"""
function nucmarks(L, T, ref)
    M = Any[]
    for (key, j) in (("D", 1), ("pp", 2), ("CNO", 4))
        f = [(v = T[(L[i].step, L[i].t)]; ref[i] > 0 ? v[j] / ref[i] : NaN) for i in eachindex(L)]
        on = [isfinite(v) && v > NUCTHR for v in f]
        i = findfirst(i -> all(@view on[i:min(i + 2, end)]), eachindex(on))
        fin = filter(isfinite, f)
        push!(M, (; key, state = i === nothing ? "none" : i == 1 ? "before" : "at", x = L[something(i, 1)],
                    Ms = T[(L[something(i, 1)].step, L[something(i, 1)].t)][7], frac = isempty(fin) ? NaN : fin[end]))
    end
    M
end

function export_frames(dir; every = 1, old = nothing, out = joinpath(dir, "export"), budget = nothing, skip = String[])
    L = old === nothing ? lineage(dir) : lineage_old(old); mkpath(out)
    cvon = any(hascd, L)
    fields = filter(f -> !(f[1] in skip), cvon ? vcat(FIELDS, CVFIELDS) : FIELDS); scal = cvon ? vcat(SCAL, CVSCAL) : SCAL
    if budget !== nothing      # the smallest stride for which frames.f32 stays below `budget` bytes
        every = 1; while sum(4 * length(x.ρ) * length(fields) for x in L[1:every:end]; init = 0) > budget; every += 1; end
    end
    sel = L[1:every:end]; (isempty(sel) || sel[end] !== L[end]) && !isempty(L) && push!(sel, L[end])
    nuc = !isempty(L) && all(hasnuc, L); T = nuc ? nuctable(L, out) : nothing
    ref = nuc ? Dict(objectid(x) => v for (x, v) in zip(L, nucref(L, T))) : nothing
    fr = Any[]; off = 0
    open(joinpath(out, "frames.f32"), "w") do io
        for x in sel
            F, S = nuc ? frame(x, vcat(T[(x.step, x.t)][1:5], ref[objectid(x)]), cvon) : frame(x, fill(NaN, 6), cvon); n = length(x.ρ)
            for (k, _) in fields; write(io, Float32.(get(F, k, fill(NaN, n)))); end
            push!(fr, (S, n, off)); off += 4 * n * length(fields)
        end
    end
    open(joinpath(out, "track.csv"), "w") do io
        println(io, join(scal, ","))
        for (S, _, _) in fr; println(io, join((@sprintf("%.12g", v) for v in S), ",")); end    # 12 digits: Δt of 1e4 s at t = 2e7 yr
    end
    js(v) = isfinite(v) ? @sprintf("%.12g", v) : "null"
    open(joinpath(out, "frames.json"), "w") do io
        print(io, "{\"fields\":[", join(("[\"$k\",\"$d\"]" for (k, d) in fields), ","), "],\"scalars\":[", join(("\"$s\"" for s in scal), ","), "],")
        if nuc      # onset marks, found on every snapshot of the lineage
            print(io, "\"nuc\":{\"thr\":$NUCTHR,\"marks\":[", join(("{\"key\":\"$(m.key)\",\"state\":\"$(m.state)\",\"step\":$(m.x.step),\"t_yr\":" *
                js(m.x.t / yr) * ",\"M_star\":" * js(m.Ms) * ",\"frac\":" * js(m.frac) * "}" for m in nucmarks(L, T, nucref(L, T))), ","), "]},")
        end
        print(io, "\"frames\":[", join(("{\"n\":$n,\"off\":$o,\"s\":[" * join(js.(S), ",") * "]}" for (S, n, o) in fr), ","))
        print(io, "]}")
    end
    println("exported $(length(sel)) of $(length(L)) snapshots (steps $(L[1].step)..$(L[end].step)) to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # runrec.jl <dir> [every | budget=<bytes>] [skip=<field>,...]   or   runrec.jl --old <snaps file> <outdir> [every]
    if ARGS[1] == "--old"
        export_frames(ARGS[3]; old = ARGS[2], out = ARGS[3], every = length(ARGS) > 3 ? parse(Int, ARGS[4]) : 1)
    else
        kw = Dict{Symbol,Any}()
        for a in ARGS[2:end]
            if startswith(a, "budget="); kw[:budget] = round(Int, parse(Float64, a[8:end]))
            elseif startswith(a, "skip="); kw[:skip] = String.(split(a[6:end], ","))
            else kw[:every] = parse(Int, a) end
        end
        export_frames(ARGS[1]; kw...)
    end
end
