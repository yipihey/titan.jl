# Variational gravity forces on a checkpoint state: the node-mass force F_j = g_j m̄_j (= −∂W/∂r_j at fixed masses for W = −Σ_j G M_j m̄_j/r_j) against the force of the zone-centred
# energy W_c = −Σ_k G m_k M_k^c / r_k^c (M_k^c = M_k + ½ m_k, r_k^c = ½(r_k + r_k+1)): F_c,j = −½(m_{j-1} g^c_{j-1} + m_j g^c_j), g^c_k = G M_k^c/(r_k^c)².
#   TITAN_NOLIB=1 julia --project=env gravity_ratio_check.jl path/to/ck/m2c_nw0.tstate.gz
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
const TD = TitanDriver; const TF = TD.TF
TF.load_state!(ARGS[1])
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); G4 = 4π * TF.getscalar(:cgrav)
r = TF.var(:rn); d = TF.var(:dn); xm = TF.var(:xmn); dv = TF.var(:dvoln)
m(k) = d[k] * dv[k]
gc(k) = (rc = 0.5 * (r[k] + r[k+1]); G4 * (xm[k] + 0.5 * m(k)) / rc^2)
@printf("%5s %12s %12s %9s\n", "zone", "F_node", "F_centred", "ratio")
rat = Float64[]
for j in a+1:b
    Fn = G4 * xm[j] / r[j]^2 * 0.5 * (m(j-1) + m(j)); Fc = 0.5 * (m(j-1) * gc(j-1) + m(j) * gc(j)); push!(rat, Fc / Fn)
    (j - a + 1) % 5 == 0 && @printf("%5d %12.4e %12.4e %9.5f\n", j - a + 1, Fn, Fc, Fc / Fn)
end
@printf("ratio F_centred/F_node over interfaces 2..%d: min %.5f median %.5f max %.5f; zones 2-26 (star): min %.5f max %.5f; zones 27-50: min %.5f max %.5f\n", b - a + 1, minimum(rat), sort(rat)[end÷2], maximum(rat),
        minimum(rat[1:25]), maximum(rat[1:25]), minimum(rat[26:49]), maximum(rat[26:49]))
