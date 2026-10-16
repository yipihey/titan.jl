# Summarise a quasi-star infall run from its snapshots (qstar_run.jl relax): per snapshot the shock radius and the
# post-shock peak temperature, the photosphere (Rosseland τ = 2/3 from the outer boundary) and its effective temperature,
# the emergent luminosity, the mass flux at the outer boundary and just outside the shock, and the mass on the grid.
#   TITAN_NOLIB=1 Z0_TABLE=ext julia --project=env qstar_analyze.jl NAME [stride]
using Serialization, Printf
include(joinpath(@__DIR__, "Z0Physics.jl")); const Z0 = Z0Physics
name = ARGS[1]; stride = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(name)_relax_snaps.jls"))
meta = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "setup_meta.jls"))
tab = Z0.table(); σ = 5.670374e-5; Msun = 1.989e33; yr = 3.156e7
R0 = meta.R; L0 = S[1].Lrad[end-1]
mass(s) = sum(4π / 3 * (s.r[k+1]^3 - s.r[k]^3) * s.d[k] for k in eachindex(s.d))
M0 = mass(S[1])
function summary(s)
    nz = length(s.d)
    # photosphere: integrate Rosseland τ inward from the outer boundary
    τ = 0.0; kph = 1
    for k in nz:-1:1
        τ += Z0.kap(tab, s.d[k], s.T[k]) * s.d[k] * (s.r[k+1] - s.r[k])
        τ >= 2 / 3 && (kph = k; break)
    end
    rph = s.r[kph + 1]; Lph = s.Lrad[kph+1] + s.Lconv[kph+1]
    Teff = (Lph / (4π * σ * rph^2))^0.25
    # shock: largest velocity jump between adjacent interfaces (infall u < 0 above, settling below)
    du = [s.u[k+1] - s.u[k] for k in 1:length(s.u)-1]
    ks = argmin(du)                       # most negative du: u drops outward into the infall
    Tpk = maximum(s.T[max(1, ks - 3):min(nz, ks + 3)])
    Mdot = -4π * s.r[end]^2 * s.d[end] * s.u[end]                     # at the outer boundary
    j = min(ks + 4, length(s.u) - 1)                                  # just outside the shock
    Mdsh = -4π * s.r[j]^2 * 0.5 * (s.d[j-1] + s.d[j]) * s.u[j]
    (t = s.t, dt = s.dt, rsh = s.r[ks+1] / R0, Tpk, rph = rph / R0, Teff, Lph = Lph / L0, Lout = s.Lrad[end] / L0,
     Mdot = Mdot / (Msun / yr), Mdsh = Mdsh / (Msun / yr), dM = (mass(s) - M0) / Msun, uin = s.u[end])
end
@printf("%-11s %-10s %-9s %-10s %-9s %-7s %-7s %-7s %-10s %-10s %-10s\n", "t [s]", "dt", "r_sh/R", "T_peak", "r_ph/R", "T_eff", "L_ph/L0",
        "L_out", "Mdot_out", "Mdot_sh", "ΔM [Msun]")
for i in 1:stride:length(S)
    q = summary(S[i])
    @printf("%-11.4e %-10.3e %-9.6f %-10.3e %-9.6f %-7.0f %-7.4f %-7.4f %-10.3e %-10.3e %-10.3e\n", q.t, q.dt, q.rsh, q.Tpk, q.rph, q.Teff,
            q.Lph, q.Lout, q.Mdot, q.Mdsh, q.dM)
end
