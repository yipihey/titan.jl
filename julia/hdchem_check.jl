# Offline checks of the HDCHEM block of PopIIIMaterial.jl: rates and cooling fits against ChemistryKernels.jl / EmissionKernels.jl,
# the chain step against BigFloat linear algebra, the cooling function.
#   julia --project=<env with ForwardDiff> hdchem_check.jl      (ENV CKPATH, EKPATH: the two packages; default ~/Projects/...)
# The block is cut out of PopIIIMaterial.jl and evaluated in a stub module, so the driver is not loaded.
using Printf, ForwardDiff, LinearAlgebra, Random
push!(LOAD_PATH, get(ENV, "CKPATH", joinpath(homedir(), "Projects", "ChemistryKernels.jl")))
push!(LOAD_PATH, get(ENV, "EKPATH", joinpath(homedir(), "Projects", "EmissionKernels.jl")))
import ChemistryKernels; const CK = ChemistryKernels
include(joinpath(@__DIR__, "PrimordialChem.jl"))
module T
using Printf, ForwardDiff
const PC = Main.PrimordialChem
module CNn; const mH = 1.6726e-24; μH() = mH * (1 + 4 / 12); full(x, net) = x; end
struct NetGas end
const STARLIGHT = Ref{Any}(nothing); const SLDIAG = Ref{Any}(nothing)
chem_core!(m) = 0.0
module TF; getscalar(s) = 0.0; end
module TM; v(s) = Float64[]; end
src = read(joinpath(@__DIR__, "PopIIIMaterial.jl"), String)
i0 = findfirst("# ------------------------------------------------------------------------------------------ deuterium and HD (HDCHEM=1)", src)
i0 === nothing && error("HDCHEM block not found in PopIIIMaterial.jl")
src = replace(src[first(i0):end], "end # module" => "")
include_string(@__MODULE__, src)
end
Ts = 10 .^ range(1, 6, length = 501)
rel(a, b) = (a == b) ? 0.0 : abs(a - b) / max(abs(a), abs(b))
for (nm, f, g) in (("k50", T.hd_k50, CK.k50), ("k51", T.hd_k51, CK.k51), ("k52", T.hd_k52, CK.k52), ("k53", T.hd_k53, CK.k53), ("k54", T.hd_k54, CK.k54),
                   ("k55", T.hd_k55, CK.k55), ("k56", T.hd_k56, CK.k56))
    e = [rel(f(t), g(t)) for t in Ts]; i = argmax(e)
    @printf("%s  max rel diff %.2e at T = %.3g   values at 30/100/300/1000/1e4 K: %.3e %.3e %.3e %.3e %.3e\n", nm, e[i], Ts[i], f(30.0), f(100.0), f(300.0), f(1e3), f(1e4))
end
Hl = isdefined(CK, :HDlte) ? CK.HDlte : getfield(Base.require(Main, :EmissionKernels), :HDlte)
Hw = isdefined(CK, :HDlow) ? CK.HDlow : getfield(Base.require(Main, :EmissionKernels), :HDlow)
for (nm, f, g) in (("HDlte", T.hd_lte, Hl), ("HDlow", T.hd_low, Hw))
    Tc = 10 .^ range(0.5, 5, length = 451); e = [rel(f(t), g(t)) for t in Tc]; i = argmax(e)
    @printf("%s  max rel diff %.2e at T = %.3g   values at 30/60/100/300/1000 K: %.3e %.3e %.3e %.3e %.3e\n", nm, e[i], Tc[i], f(30.0), f(60.0), f(100.0), f(300.0), f(1e3))
end
# chain step against BigFloat: p (I − dt Q) = p0, Q generator (rows sum to zero)
Random.seed!(7); worst = 0.0; worstc = 0.0; worsts = 0.0
for it in 1:20000
    q = Tuple(10.0 .^ (rand(6) .* 40 .- 30) .* (rand(6) .> 0.15))
    p0r = rand(3) .* (rand(3) .> 0.2); sum(p0r) == 0 && (p0r = [1.0, 0, 0]); p0 = Tuple(p0r ./ sum(p0r))
    dt = 10.0^(rand() * 30 - 10)
    q12, q21, q23, q32, q13, q31 = big.(q)
    Q = [-(q12 + q13) q12 q13; q21 -(q21 + q23) q23; q31 q32 -(q31 + q32)]
    A = I - big(dt) * Q
    pref = (collect(big.(p0))' / A)'
    p = T.hd_step(q, p0, 1 / dt)
    global worst = max(worst, maximum(abs.(collect(p) .- Float64.(pref)))); global worstc = max(worstc, abs(sum(p) - 1))
    # relative accuracy of the components that are not negligible
    for i in 1:3; pref[i] > 1e-30 && (global worsts = max(worsts, Float64(abs(p[i] - pref[i]) / pref[i]))); end
end
@printf("chain step, 20000 random cases (rates 1e-30..1e10 1/s with zeros, dt 1e-10..1e20 s): max abs err %.2e, max rel err of any component > 1e-30: %.2e, max |sum-1| %.2e\n", worst, worsts, worstc)
# stationary state = limit of long steps; and backward Euler converges to the exact exponential for small steps
q = (3e-13, 2e-12, 5e-14, 1e-15, 4e-16, 2e-17); ps = T.hd_step(q, (1.0, 0.0, 0.0), 0.0); pl = T.hd_step(q, (0.2, 0.3, 0.5), 1e-40)
@printf("stationary %s ; step with dt = 1e40 s %s\n", string(ps), string(pl))
Q = [-(q[1] + q[5]) q[1] q[5]; q[2] -(q[2] + q[3]) q[3]; q[6] q[4] -(q[6] + q[4])]
for ns in (1, 10, 100, 1000)
    tt = 3e12; p = (1.0, 0.0, 0.0); for i in 1:ns; p = T.hd_step(q, p, ns / tt); end
    pe = ([1.0 0 0] * exp(Q * tt))
    @printf("  %4d backward-Euler steps over 3e12 s: HD fraction %.6e (exact %.6e)\n", ns, p[3], pe[3])
end
# cooling: values, sign below T_CMB, dual numbers
T.HDTCMB[] = 54.5
d = Dict(:H => 0.999, :Hp => 1e-7, :e => 1e-7, :H2 => 5e-4, :Hm => 1e-12)
for (n, t) in ((1e4, 300.0), (1e4, 100.0), (1e4, 54.5), (1e4, 30.0), (1e8, 300.0), (1e2, 300.0))
    ρ = n * T.CNn.μH(); c = T.hd_cool(ρ, t, d, 2.5e-5); dc = ForwardDiff.derivative(x -> T.hd_cool(ρ, x, d, 2.5e-5), t)
    ek = n * 2.5e-5 * (Hl(t) / (1 + Hl(t) / (n * 0.999 * Hw(t))))
    @printf("  n %.0e T %6.1f: HD cooling %.4e erg/cm3/s (without the CMB term, from EmissionKernels fits: %.4e), dΛ/dT %.3e\n", n, t, c, ek, dc)
end
dd = Dict(k => ForwardDiff.Dual(v, 1.0) for (k, v) in d); println("dual composition: ", T.hd_cool(1e4 * T.CNn.μH(), ForwardDiff.Dual(200.0, 1.0), dd, 2.5e-5))
println("rates at 300 K, nH 1e4: ", T.hd_rates(300.0, 1e4, d, 0.0))
