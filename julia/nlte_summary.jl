# Summary numbers of a non-LTE solution (nlte_run.jl): column densities, Thomson depth, Hα source function, energy.
#   julia --project=env nlte_summary.jl OUT [R=cm]   (lib/spectra/nlte_OUT.jls; R: the star's radius, the unit of r/R,
#   default 8.4352e15 cm, run W30)
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const NH = NLTEH; const Z0 = Z0Physics; const SP = Spectrum
o = ARGS[1]
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$o.jls"))
sp = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$(d.name).jls"))
R0 = length(ARGS) > 1 ? parse(Float64, split(ARGS[2], "=")[2]) : 8.4352e15; r = d.r
dr = [q < length(r) ? r[q+1] - r[q] : 0.0 for q in eachindex(r)]     # segment to the next point outward
col(f) = sum(f(q) * dr[q] for q in eachindex(r))
nHt(q) = d.nl[q].np + sum(d.nl[q].n)
s_lte(q) = Z0.saha(sp.atm.ρ[d.sel[q]], d.T[q])
@printf("%s: cscale %.2g, caseB %s, X-rays %s\n", o, d.cscale, d.caseB, d.xray)
@printf("  N_H %.3e; neutral N_H I: non-LTE %.3e, LTE %.3e\n", col(nHt), col(q -> sum(d.nl[q].n)), col(q -> s_lte(q).nH0))
for (lo, hi) in ((0.99, 1.21), (1.21, 5.0), (5.0, 200.0))
    m = findall(q -> lo <= r[q] / R0 < hi, eachindex(r)); isempty(m) && continue
    @printf("    r/R %5.2f–%6.1f: N_H I non-LTE %.3e, LTE %.3e (of N_H %.3e)\n", lo, hi, sum(sum(d.nl[q].n) * dr[q] for q in m),
            sum(s_lte(q).nH0 * dr[q] for q in m), sum(nHt(q) * dr[q] for q in m))
end
@printf("  Thomson depth: non-LTE %.3f, LTE %.3f\n", col(q -> d.nl[q].ne * Z0.σT), col(q -> s_lte(q).ne * Z0.σT))
ν32 = (NH.En[3] - NH.En[2]) / Z0.h; λ = Z0.c / ν32
println("  Hα in the infall (r/R, x_H, S_L/B(T), S_L/J_c, τ_Sobolev radial-tangential average K/(v/r)):")
for q in eachindex(r)
    (q == 1 || q % 10 == 0 || q == length(r)) || continue
    st = d.nl[q]; T = d.T[q]
    SL = 2Z0.h * ν32^3 / Z0.c^2 / (st.n[2] * NH.gn[3] / (st.n[3] * NH.gn[2]) - 1)
    K = NH.πe2mc * NH.fosc(2, 3) * λ * (st.n[2] - NH.gn[2] / NH.gn[3] * st.n[3])
    @printf("    %8.4f  x %.4f  S/B %9.3e  S/J %9.3e  τ %9.3e   (LTE τ %9.3e)\n", r[q] / R0, st.np / nHt(q), SL / SP.Bν(ν32, T), SL / d.Jhα[q],
            K / (d.vabs[q] / r[q]), let s = s_lte(q), pH = Z0.level_pops(T, s)[1]; NH.πe2mc * NH.fosc(2, 3) * λ * pH[2] * (1 - exp(-Z0.h * ν32 / (Z0.kB * T))) / (d.vabs[q] / r[q]) end)
end
Lc = sum(0.5 * (d.L[i] + d.L[i+1]) * (d.ν[i+1] - d.ν[i]) for i in 1:length(d.ν)-1)
Ll = sum(d.nl[q].P * d.dV[q] for q in eachindex(r))
@printf("  energy: continuum %.4e + lines %.3e = %.4e erg/s; TITAN L_rad(out) %.4e (ratio %.4f)\n", Lc, Ll, Lc + Ll, sp.Ltitan, (Lc + Ll) / sp.Ltitan)
