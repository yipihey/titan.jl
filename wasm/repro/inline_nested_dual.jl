using WasmTarget, ForwardDiff
const kB = 1.380649e-16; const h = 6.62607015e-27; const me = 9.1093837e-28; const eV = 1.602176634e-12
const mH = 1.6735575e-24; const mHe = 6.6464731e-24
const χH = 13.598434eV; const χHe1 = 24.587389eV; const χHe2 = 54.417765eV
λe3(T) = (2π * me * kB * T / h^2)^1.5
_rates(T) = (L = λe3(T); kT = kB * T; (L * exp(-χH / kT), 4L * exp(-χHe1 / kT), L * exp(-χHe2 / kT)))
function polish(f, x0)
    x = x0 + zero(f(x0))
    for _ in 1:2
        fx = f(x); x = x - fx / ForwardDiff.derivative(f, x)
    end
    x
end
@inline function charge(ne, nH, nHe, SH, S1, S2)
    xH = SH / (ne + SH); r1 = S1 / ne; r2 = S2 / ne
    f0 = 1 / (1 + r1 + r1 * r2); y1 = r1 * f0; y2 = r1 * r2 * f0
    return nH * xH + nHe * (y1 + 2y2)
end
const NE0 = 4.538237247137191e11
function k(ρ, ncap)
    nH = 0.76ρ / mH; nHe = 0.24ρ / mHe; SH, S1, S2 = _rates(1e4)
    if ncap == 2
        return exp(polish(lne -> log(charge(exp(lne), nH, nHe, 1.0e15, 1.0e4, 1.0)) - lne, log(NE0)))
    elseif ncap == 3
        return exp(polish(lne -> log(charge(exp(lne), nH, nHe, SH, 1.0e4, 1.0)) - lne, log(NE0)))
    else
        return exp(polish(lne -> log(charge(exp(lne), nH, nHe, SH, S1, S2)) - lne, log(NE0)))
    end
end
x2(ρ::Float64)::Float64 = ForwardDiff.derivative(r -> k(r, 2), ρ)
x3(ρ::Float64)::Float64 = ForwardDiff.derivative(r -> k(r, 3), ρ)
x5(ρ::Float64)::Float64 = ForwardDiff.derivative(r -> k(r, 5), ρ)
# capture a single tuple instead of 5 scalars
function kt(ρ)
    nH = 0.76ρ / mH; nHe = 0.24ρ / mHe; S = _rates(1e4)
    exp(polish(lne -> log(charge(exp(lne), nH, nHe, S[1], S[2], S[3])) - lne, log(NE0)))
end
xt(ρ::Float64)::Float64 = ForwardDiff.derivative(kt, ρ)
for (f, x) in ((x2, 1e-12), (x3, 1e-12), (x5, 1e-12), (xt, 1e-12))
    print(nameof(f), " native=", f(x), "  ")
    try; b = compile_multi([(f, (Float64,))]); println("OK ", length(b))
    catch e; println("FAIL ", first(sprint(showerror, e), 150)); end
end
