# Where the non-ideality criterion lands along run profiles (checkpoint states): s = ((Γ/G0)^4 + (D/D0)^4 + (P/P0)^4)^(1/4), w = blend weight.
ENV["HM88"] = "1"
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl")); include(joinpath(@__DIR__, "mats_hm.jl"))
function readstate(path)
    io = IOBuffer(read(`gzip -dc $path`)); String(read(io, 8)) == "TITANST1" || error("format")
    n = read(io, Int64); st = Dict{String,Any}()
    for _ in 1:n
        ln = read(io, Int64); name = String(read(io, ln)); code = read(io, UInt8); nd = read(io, Int64); dims = [read(io, Int64) for _ in 1:nd]
        len = prod(dims); T = code == 0x01 ? Float64 : code == 0x02 ? Int32 : NTuple{80,UInt8}
        a = reinterpret(T, read(io, sizeof(T) * len)); st[name] = reshape(collect(a), dims...)
    end
    st
end
const G0, D0, P0 = 0.3, 0.1, 0.01
sval(nd) = ((nd.Γ / G0)^4 + (nd.D / D0)^4 + (nd.Pneu / P0)^4)^0.25
wval(s) = (x = clamp(log10(3s) / log10(3.0), 0.0, 1.0); x * x * (3 - 2x))
for (lab, path) in ((ARGS[i], ARGS[i+1]) for i in 1:2:length(ARGS))
    st = readstate(path); ρ = vec(st["dn"]); T = vec(st["tn"]); r = vec(st["rn"])
    a = findfirst(>(0), ρ); b = findlast(>(0), ρ)
    println("== ", lab, " (", path, "): zones ", a, "..", b)
    first_w0 = nothing; nblend = 0; ntab = 0
    for k in a:b
        (ρ[k] > 0 && T[k] > 0) || continue
        T[k] > 3e5 && (ntab += 1; continue)        # HOT branch (stellar EOS) above 3e5 K
        nd = nonideal(ρ[k], min(T[k], 9.9e5)); s = sval(nd); w = wval(s)
        w >= 1 && (ntab += 1); 0 < w < 1 && (nblend += 1)
        if w < 1 && first_w0 === nothing
            first_w0 = k
        end
    end
    # print a compact profile: every 4th zone between the HOT edge and the photosphere-ish
    for k in a:max(1, (b - a) ÷ 40):b
        (ρ[k] > 0 && T[k] > 0) || continue
        T[k] > 3e5 && (@printf("  k %3d r %.3e  log ρ %6.2f log T %5.2f  HOT\n", k, r[k] / 6.957e10, log10(ρ[k]), log10(T[k])); continue)
        nd = nonideal(ρ[k], T[k]); s = sval(nd)
        @printf("  k %3d r %.3e R☉ log ρ %6.2f log T %5.2f  Γ %.2e D %.2e P_neu %.2e  s %.2f w %.2f\n", k, r[k] / 6.957e10, log10(ρ[k]), log10(T[k]), nd.Γ, nd.D, nd.Pneu, s, wval(s))
    end
    @printf("  -> zones on the table (w = 1, or HOT) %d, in the blend %d, outermost blended/table zone boundary at k = %s\n", ntab, nblend, string(first_w0))
end
