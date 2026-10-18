# Per-line energy check of nlte_spectrum.jl: the ray-integrated line luminosity (from the emergent profile) against the
# local Sobolev net line power of the same non-LTE solution, Σ hν β [n_u A − (n_l B_lu − n_u B_ul) J_c] dV, and its gross
# (emission-only) part. Last line: totals over the listed lines and the SE sum over all 66 transitions.
#   Z0_TABLE=ext TITAN_NOLIB=1 julia -t 2 --project=env nlte_linepower.jl NAME
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const SP = Spectrum; const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]
d = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "$name.jls")); dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
ls = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "linespec_$name.jls"))
atm = merge(d.atm, (T = dn.Tstar,)); inf = dn.sel
nlfull = Any[nothing for _ in eachindex(atm.r)]; for (q, j) in enumerate(inf); nlfull[j] = dn.nl[q]; end
op = NH.NOpac(SP.Opac(atm), nlfull); _, R = SP.rays(atm.r; ncore = 12); R = SP.weights!(atm, R)
h = Z0.h; c = Z0.c
tot = [0.0, 0.0, 0.0]
for x in ls.lines
    l, u = x.l, x.u; ν = c / (x.λ * 1e-8); f = NH.fosc(l, u)
    _, J, _ = SP.monochromatic(atm, op, R, ν)
    Blu = 4π * NH.πe2mc * f / (h * ν); Bul = NH.gn[l] / NH.gn[u] * Blu; A = 2h * ν^3 / c^2 * Bul
    Pn = 0.0; Pg = 0.0
    for (q, j) in enumerate(inf)
        st = dn.nl[q]; β = st.β[l, u]
        Pg += h * ν * β * st.n[u] * A * dn.dV[q]
        Pn += h * ν * β * (st.n[u] * A - (st.n[l] * Blu - st.n[u] * Bul) * J[j]) * dn.dV[q]
    end
    tot .+= (Pn, Pg, x.Lline)
    @printf("%-5s Sobolev net %.3e  gross %.3e  ray %.3e  ray/net %.2f\n", x.name, Pn, Pg, x.Lline, x.Lline / Pn)
end
@printf("total Sobolev net %.3e gross %.3e ray %.3e (SE P sum %.3e)\n", tot..., sum(dn.nl[q].P * dn.dV[q] for q in eachindex(dn.nl)))
