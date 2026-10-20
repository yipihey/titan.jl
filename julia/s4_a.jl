# S4 test a: with the nuclear network off, StarGas (Z0 opacity table) reproduces StellarGas bit-for-bit in a TITAN rad-hydro run.
include(joinpath(@__DIR__, "s4_star.jl"))
nst = parse(Int, get(ENV, "NSTEPS", "25"))
function run(kind)
    st = build_star(; M = 100Msun, R = 6.5Rsun, n = 40, nenv = 12, nuc = false, opac = :z0, kind)
    TF.setscalar!(:dtime, 1e3)
    TD.run_julia!(maxsteps = nst)
    (; t = TF.getscalar(:timen), (Symbol(s) => copy(TF.var(s)) for s in (:rn, :dn, :tn, :ern, :frn, :un, :pgn, :egn))...)
end
A = run(:stellar); B = run(:star)
dmax = 0.0; ident = true
for s in (:rn, :dn, :tn, :ern, :frn, :un, :pgn, :egn)
    d = maximum(abs.(A[s] .- B[s]) ./ (abs.(A[s]) .+ 1e-300)); global dmax = max(dmax, d); global ident &= A[s] == B[s]
    @printf("  %-4s max rel. diff %.3e\n", s, d)
end
@printf("%d steps, t = %.6e s (both: %.6e); bit-identical: %s; max relative difference %.2e\n", nst, A.t, B.t, ident, dmax)
println(dmax <= 1e-14 ? "PASS" : "FAIL")
