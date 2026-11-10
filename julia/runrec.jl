# Reader and exporter for runs recorded with RUNREC=<dir> (popiii_m2c.jl).
#
#   julia --project=<env> runrec.jl <dir> [every]     # writes <dir>/export/{frames.json, frames.f32, track.csv}
#
# The run directory holds one snapshot stream per process segment (snaps_<k0>_<pid>.jls). A restart from an earlier checkpoint
# re-runs steps that an abandoned segment had already written, so the surviving history is found by reading all segments in
# start order and keeping, from the newest snapshot backwards, only snapshots strictly earlier in time than the one after them.
# A segment killed mid-write ends in a truncated record, which is dropped.
#
# frames.f32 is Float32, frame after frame, field after field, nz values each (nz varies between frames: rezoning, the remap
# onto the accretion mesh). frames.json has the field names and units, and per frame the step, time, nz, byte offset and the
# scalars of track.csv. Cell-centred quantities; values at interfaces (u, L, Ṁ, F/cE) are averaged onto the cells.
using Serialization, Printf

const Rs = 6.957e10; const Ms = 1.989e33; const Ls = 3.83e33; const yr = 3.156e7; const G = 6.674e-8
const cl = 2.998e10; const arad = 7.5657e-15; const sb = 5.6704e-5; const mH = 1.6726e-24

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
    fs = sort(filter(f -> startswith(f, "snaps_") && endswith(f, ".jls"), readdir(dir)); by = f -> mtime(joinpath(dir, f)))
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
    ("logr_o", "log r of the zone's outer interface [R_sun]"), ("logm_o", "log m at the zone's outer interface [M_sun]")]
const SCAL = ["step", "t_yr", "dt_s", "logn_c", "T_c", "M_star", "R_star", "r_ph", "T_eff_lum", "T_eff_sl", "L_ph", "Mdot_1e4", "u_min", "L_nuc", "L_lines", "nz", "nchemfb"]

function frame(x)
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
    # scalars
    js = findfirst(<(-2e5), x.u); Mst = js === nothing ? NaN : m[js] / Ms; Rst = js === nothing ? NaN : r[js] / Rs
    kph = findlast(>(2 / 3), τ); rph = kph === nothing ? NaN : r[kph] / Rs
    Lph = kph === nothing ? NaN : L[kph]; Tl = (kph === nothing || Lph <= 0) ? NaN : (Lph * Ls / (4π * r[kph]^2 * sb))^0.25
    j4 = findfirst(>=(1e4 * Rs), r); md4 = j4 === nothing ? NaN : md[j4]
    S = [x.step, x.t / yr, x.dt, log10(x.ρ[1] / mH), x.T[1], Mst, Rst, rph, Tl, x.Teff, Lph, md4, minimum(x.u) / 1e5,
         x.Lnuc, x.Llines, n, x.nchemfb]
    F, S
end

function export_frames(dir; every = 1, old = nothing, out = joinpath(dir, "export"))
    L = old === nothing ? lineage(dir) : lineage_old(old); mkpath(out)
    sel = L[1:every:end]; (isempty(sel) || sel[end] !== L[end]) && !isempty(L) && push!(sel, L[end])
    fr = Any[]; off = 0
    open(joinpath(out, "frames.f32"), "w") do io
        for x in sel
            F, S = frame(x); n = length(x.ρ)
            for (k, _) in FIELDS; write(io, Float32.(F[k])); end
            push!(fr, (S, n, off)); off += 4 * n * length(FIELDS)
        end
    end
    open(joinpath(out, "track.csv"), "w") do io
        println(io, join(SCAL, ","))
        for (S, _, _) in fr; println(io, join((@sprintf("%.12g", v) for v in S), ",")); end    # 12 digits: Δt of 1e4 s at t = 2e7 yr
    end
    js(v) = isfinite(v) ? @sprintf("%.12g", v) : "null"
    open(joinpath(out, "frames.json"), "w") do io
        print(io, "{\"fields\":[", join(("[\"$k\",\"$d\"]" for (k, d) in FIELDS), ","), "],\"scalars\":[", join(("\"$s\"" for s in SCAL), ","), "],\"frames\":[")
        print(io, join(("{\"n\":$n,\"off\":$o,\"s\":[" * join(js.(S), ",") * "]}" for (S, n, o) in fr), ","))
        print(io, "]}")
    end
    println("exported $(length(sel)) of $(length(L)) snapshots (steps $(L[1].step)..$(L[end].step)) to $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # runrec.jl <dir> [every]   or   runrec.jl --old <snaps file> <outdir> [every]
    if ARGS[1] == "--old"
        export_frames(ARGS[3]; old = ARGS[2], out = ARGS[3], every = length(ARGS) > 3 ? parse(Int, ARGS[4]) : 1)
    else
        export_frames(ARGS[1]; every = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1)
    end
end
