# photfluct analysis: per-snapshot photosphere / shock diagnostics of a RUNREC directory
# (photfluct study, 16 Nov 2026; run from scripts/photfluct with the JLENV project: julia ana.jl <run dir> <out.csv> [minstep])
#   julia ana.jl <dir> <out.csv> [minstep]
using Serialization, Printf
const Rs = 6.957e10; const Ms = 1.989e33; const Ls = 3.83e33; const yr = 3.156e7; const G = 6.674e-8
const cl = 2.998e10; const arad = 7.5657e-15; const sb = 5.6704e-5; const kB = 1.380649e-16; const mH = 1.6726e-24
isdefined(Main, :NucNet) || (@eval Main module NucNet end)   # (skip the network: runrec only needs it for burning)
include(joinpath(@__DIR__, "..", "..", "julia", "runrec.jl"))

function diag(x)
    n = length(x.ρ); r = x.r; m = x.m
    τ = zeros(n + 1); for k in n:-1:1; τ[k] = τ[k+1] + x.chi[k] * x.ρ[k] * (r[k+1] - r[k]); end
    L = 4π .* r .^ 2 .* x.fr ./ Ls
    kph = findlast(>(2 / 3), τ)                    # interface index (runrec convention)
    js = findfirst(<(-2e5), x.u)                    # first interface with infall > 2 km/s (shock)
    Lph = L[kph]; rph = r[kph]
    zp = min(kph, n)                                # zone just outside the photospheric interface
    zi = max(kph - 1, 1)                            # zone inside
    Trad(k) = (max(x.er[k], 0) / arad)^0.25
    Tl = Lph > 0 ? (Lph * Ls / (4π * rph^2 * sb))^0.25 : NaN
    je = something(findfirst(>=(1.5rph), r), n + 1); Lem = L[je]; Lk2 = L[min(kph + 2, n + 1)]
    Tem = Lem > 0 ? (Lem * Ls / (4π * rph^2 * sb))^0.25 : NaN
    # accretion rate just outside the shock and accretion luminosity
    ja = min(js + 1, n + 1); mdot = -4π * r[ja]^2 * x.ρ[min(ja, n)] * x.u[ja]
    M = m[js]; R = r[js]; Lacc = G * M * mdot / R / Ls
    # post-shock zone (inside the shock interface): cooling, diffusion, sound crossing, advection
    k = js - 1; dr = r[k+1] - r[k]
    e = x.e[k]; ρ = x.ρ[k]; T = x.T[k]
    Λ = x.kp[k] * ρ * cl * (arad * T^4 - x.er[k])  # erg/cm3/s net emission
    tcool = abs(Λ) > 0 ? ρ * e / abs(Λ) : Inf
    τz = x.chi[k] * ρ * dr; tdiff = max(τz, 1) * dr / cl
    cs = sqrt(max(x.pg[k], 0) / ρ); tsound = dr / cs
    upre = abs(x.u[js]); tadv = dr / max(upre, 1)
    # pre-shock zone
    k2 = js; dr2 = r[k2+1] - r[k2]; τz2 = x.chi[k2] * x.ρ[k2] * dr2
    (; step = x.step, t = x.t, dt = x.dt, n, kph, js, rph = rph / Rs, Lph, Tl, Lem, Lk2, Tem, Trad_ph = Trad(zp), Tgas_ph = x.T[zp], Trad_in = Trad(zi),
       L_out = L[end], L_sh = L[js], L_sh2 = L[min(js + 2, n + 1)], Lmax = maximum(L), Rst = R / Rs, Mst = M / Ms, mdot = mdot * yr / Ms, Lacc,
       tacc = M / mdot, tcool, tdiff, tsound, tadv, τz, τz2, dr = dr / Rs, Tps = T, ρps = ρ, upre = upre / 1e5, τsh = τ[js])
end

dir = ARGS[1]; out = ARGS[2]; mins = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0
L = lineage(dir)
ok_(x) = (findfirst(<(-2e5), x.u) !== nothing && findfirst(<(-2e5), x.u) > 1)
rows = filter(!isnothing, [(try diag(x) catch; nothing end) for x in L if x.step >= mins && ok_(x)])
open(out, "w") do io
    ks = keys(rows[1]); println(io, join(ks, ","))
    for rw in rows; println(io, join((@sprintf("%.10g", Float64(v)) for v in values(rw)), ",")); end
end
println("wrote ", length(rows), " rows to ", out)
