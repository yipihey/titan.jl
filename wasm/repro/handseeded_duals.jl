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
function charge(ne, nH, nHe, SH, S1, S2)
    xH = SH / (ne + SH); r1 = S1 / ne; r2 = S2 / ne
    f0 = 1 / (1 + r1 + r1 * r2); y1 = r1 * f0; y2 = r1 * r2 * f0
    return nH * xH + nHe * (y1 + 2y2)
end
const NE0 = 4.538237247137191e11
function ne_of(ρ, T)
    nH = 0.76ρ / mH; nHe = 0.24ρ / mHe; SH, S1, S2 = _rates(T)
    g(lne) = log(charge(exp(lne), nH, nHe, SH, S1, S2)) - lne
    exp(polish(g, log(NE0)))
end
struct MyTag end
ForwardDiff.:≺(::Type{MyTag}, ::Type{<:ForwardDiff.Tag}) = true
ForwardDiff.:≺(::Type{<:ForwardDiff.Tag}, ::Type{MyTag}) = false
function y2p(ρ::Float64)::Float64      # explicit 2-partial Duals with a custom tag, as eos_rhoT does
    D = ForwardDiff.Dual{MyTag}
    T = ρ * 1e16; ρd = D(ρ, ρ, 0.0); Td = D(T, 0.0, T)
    ForwardDiff.partials(ne_of(ρd, Td), 2)
end
y1p(ρ::Float64)::Float64 = ForwardDiff.derivative(r -> ne_of(r, r * 1e16), ρ)   # T Dual too, 1 partial
yg(ρ::Float64)::Float64 = ForwardDiff.gradient(v -> ne_of(v[1], v[2]), [ρ, 1e4])[2]

fdummy(x) = x
function z1(ρ::Float64)::Float64      # tag from a function, like gradient's
    D = ForwardDiff.Dual{typeof(ForwardDiff.Tag(fdummy, Float64))}
    T = ρ * 1e16; ρd = D(ρ, ρ, 0.0); Td = D(T, 0.0, T)
    ForwardDiff.partials(ne_of(ρd, Td), 2)
end
function z2(ρ::Float64)::Float64      # Partials constructed explicitly
    D = ForwardDiff.Dual{typeof(ForwardDiff.Tag(fdummy, Float64))}
    T = ρ * 1e16; ρd = D(ρ, ForwardDiff.Partials((ρ, 0.0))); Td = D(T, ForwardDiff.Partials((0.0, T)))
    ForwardDiff.partials(ne_of(ρd, Td), 2)
end
function z3(ρ::Float64)::Float64      # gradient of ne_of(exp(lr), exp(lT))
    ForwardDiff.gradient(v -> ne_of(exp(v[1]), exp(v[2])), [log(ρ), log(ρ * 1e16)])[2]
end
for (f, x) in ((z1, 1e-12), (z2, 1e-12), (z3, 1e-12))
    print(nameof(f), " native=", f(x), "  ")
    try; b = compile_multi([(f, (Float64,))]); println("OK ", length(b))
    catch e; println("FAIL ", first(sprint(showerror, e), 150)); end
end
