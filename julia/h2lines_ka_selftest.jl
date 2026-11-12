# Self-test of H2LinesKA.jl (the KernelAbstractions formal solution of the H₂ lines) on a chosen backend, without any run data:
# a synthetic accreting-protostar envelope (opaque star, accretion shock, free-fall envelope, sub-thermal sources), solved
#   (1) on KernelAbstractions' CPU backend in Float64  — the reference (round-off equal to H2Lines.transfer_line_fast!), and
#   (2) on the chosen backend in Float32 (and in Float64 where the device has it).
#
#   julia -t 4 --project=<environment with KernelAbstractions and CUDA or Metal> h2lines_ka_selftest.jl metal|cuda|cpu [zones] [lines]
#
# Pass: Float32 errors of J̄ and Λ* around 1e-6 or below and net cooling within ~1e-6 of its peak (Float32 round-off); the
# script prints PASS/FAIL against 1e-4. On Apple silicon: `] add KernelAbstractions Metal` in a fresh project, then `metal`.
using Printf
include(joinpath(@__DIR__, "H2LinesKA.jl"))
using KernelAbstractions
const KA = KernelAbstractions
name = length(ARGS) >= 1 ? lowercase(ARGS[1]) : "cpu"
n = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 80
nl = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 48
name == "cuda" && @eval using CUDA
name == "metal" && @eval using Metal
backend() = name == "cuda" ? CUDA.CUDABackend() : name == "metal" ? Metal.MetalBackend() : KA.CPU()

# ---- the synthetic envelope (cgs)
const cl = 2.99792458e10; const hP = 6.62607015e-27; const kB = 1.380649e-16; const mH = 1.6726e-24; const GN = 6.674e-8
Bnu(ν, T) = 2hP * ν^3 / cl^2 / expm1(hP * ν / (kB * T))
function problem(n, nl)
    rs = 3e12                                                       # accretion shock = photosphere
    r = [0.0; exp.(range(log(5e11), log(3e17), length = n))]        # n + 1 interfaces, the first at the centre
    kq = argmin(abs.(log.(r[2:end] ./ rs))) + 1                     # the interface nearest the shock
    rc = 0.5 .* (r[1:end-1] .+ r[2:end])
    M = 10 * 1.989e33
    u = [x > r[kq] ? -sqrt(2GN * M / x) / (1 + (x / 3e16)^2) - 2e4 : 0.0 for x in r]
    T = [x > r[kq] ? max(300.0, 2500.0 * (x / rs)^-0.4) : 6000.0 for x in rc]
    vth = sqrt.(2kB .* T ./ (2mH))
    nH2 = [x > r[kq] ? 1e17 * (x / rs)^-1.5 : 0.0 for x in rc]      # no H₂ in the star: zones without opacity
    lines = Any[]
    for i in 1:nl
        f = (i - 1) / max(nl - 1, 1)
        ν0 = cl * (350.0 * (5000.0 / 350.0)^f)                      # 350 … 5000 cm⁻¹
        El = 170.0 * (1 + 60f); A = 1e-11 * 1e5^f                   # lower-level energy [K], Einstein A [1/s]; 10^(i mod 5 − 2) spreads the optical depths
        κ0 = [10.0^(mod(i, 5) - 2) * nH2[k] * exp(-El / T[k]) / (1 + T[k] / 170) * A * cl^2 / (8π * ν0^2) * (-expm1(-hP * ν0 / (kB * T[k]))) for k in 1:n]
        Sz = [Bnu(ν0, T[k] * (0.55 + 0.45 / (1 + (rc[k] / 1e15)^0.5))) for k in 1:n]     # sub-thermal outside
        act = findall(>(0), κ0)
        vspan = 3maximum(vth[act]) + maximum(max(abs(u[k]), abs(u[k+1])) for k in act)
        nb = clamp(ceil(Int, 2vspan / (0.5minimum(vth[act]))) + 1, 30, 200)
        push!(lines, (; κ0, ν0, xs = range(-vspan, vspan, length = nb), lam = true, Sz, bg = nothing))
    end
    τc = [sum(L.κ0[k] / (sqrt(π) * vth[k] / cl * L.ν0) * (r[k+1] - r[k]) for k in 1:n) for L in lines]
    (; r, u, vth, lines, kq, τc)
end
withbg(P) = [(; L.κ0, L.ν0, L.xs, L.lam, L.Sz, bg = (Bc = Bnu(L.ν0, 54.5), Bs = Bnu(L.ν0, 6000.0), kq = P.kq)) for L in P.lines]

relcol(A, B) = maximum(maximum(abs, A[:, l] .- B[:, l]) / max(maximum(abs, B[:, l]), floatmin()) for l in axes(B, 2))
function netcool(F, r)                                              # per zone, summed over the lines: 4π (r²F)|out − (r²F)|in
    n = length(r) - 1
    [sum(4π * (r[k+1]^2 * F[k+1, l] - r[k]^2 * F[k, l]) for l in axes(F, 2)) for k in 1:n]
end
function run(be, T, P, lines)
    G = H2LinesKA.raygeom(P.r, P.u)
    ws = H2LinesKA.Workspace(be, T, P.r, P.u, P.vth, G; nlmax = length(lines), nbmax = 200)
    res = H2LinesKA.formal!(ws, lines)                              # compiles
    t = minimum(@elapsed(H2LinesKA.formal!(ws, lines)) for _ in 1:3)
    res = H2LinesKA.formal!(ws, lines)
    mb = H2LinesKA.devbytes(ws) / 2^20
    H2LinesKA.free!(ws)
    (; F = res[1], J = res[2], Λ = res[3], t, mb)
end

P = problem(n, nl)
@printf("H2LinesKA self-test: %d zones, %d lines, bins %d–%d, line-centre optical depths %.1e … %.1e, %d threads, backend %s\n", n, nl,
        minimum(length(L.xs) for L in P.lines), maximum(length(L.xs) for L in P.lines), minimum(P.τc), maximum(P.τc), Threads.nthreads(), name)
ok = true
for (tag, lines) in (("no incident field", P.lines), ("CMB + photosphere ", withbg(P)))
    ref = run(KA.CPU(), Float64, P, lines)
    Nr = netcool(ref.F, P.r); mx = maximum(abs, Nr)
    @printf("  %s  reference (KA CPU Float64): %.3f s per batch\n", tag, ref.t)
    for T in (name == "metal" ? (Float32,) : (Float32, Float64))
        o = run(backend(), T, P, lines)
        eJ = relcol(o.J, ref.J); eΛ = relcol(o.Λ, ref.Λ); eF = relcol(o.F, ref.F); eN = maximum(abs, netcool(o.F, P.r) .- Nr) / mx
        tol = T === Float32 ? 1e-4 : 1e-10
        pass = all(isfinite, o.F) && max(eJ, eΛ, eN) < tol
        global ok &= pass
        @printf("  %s  %-5s %s: %.4f s per batch, device %.0f MB | J̄ %.1e  Λ* %.1e  flux %.1e  net cooling/peak %.1e  %s\n", tag, name, T, o.t,
                o.mb, eJ, eΛ, eF, eN, pass ? "PASS" : "FAIL")
    end
end
println(ok ? "PASS" : "FAIL")
