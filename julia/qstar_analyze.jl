# Summarise a quasi-star infall run from its snapshots (qstar_run.jl relax): per snapshot the shock radius and the
# post-shock peak temperature, the photosphere (Rosseland τ = 2/3 from the outer boundary) and its effective temperature,
# the emergent luminosity, the mass flux at the outer boundary and just outside the shock, and the mass on the grid.
#   TITAN_NOLIB=1 Z0_TABLE=ext julia --project=env qstar_analyze.jl NAME [stride]
using Serialization, Printf
include(joinpath(@__DIR__, "Z0Physics.jl")); const Z0 = Z0Physics
name = ARGS[1]; stride = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1
S = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", "$(name)_relax_snaps.jls"))
# QS_META=NAME: the set-up the run started from (lib/qs_runs/setup_meta_NAME.jls); default: the last set-up made
meta = deserialize(joinpath(@__DIR__, "..", "lib", "qs_runs", haskey(ENV, "QS_META") ? "setup_meta_$(ENV["QS_META"]).jls" : "setup_meta.jls"))
include(joinpath(@__DIR__, "qstar_summary.jl"))
# L0: luminosity of the steady star at t = 0 (the set-up), the same for every run and continuation (it was the first
# snapshot's outer L, which for a continuation such as A220d or W30d is not the t = 0 value)
R0 = meta.R; L0 = meta.L; M0 = QSum.mass(S[1])
summary(s) = QSum.summary(s, R0, L0, M0)
@printf("%-11s %-10s %-9s %-10s %-9s %-7s %-7s %-7s %-10s %-10s %-10s\n", "t [s]", "dt", "r_sh/R", "T_peak", "r_ph/R", "T_eff", "L_ph/L0",
        "L_out", "Mdot_out", "Mdot_sh", "ΔM [Msun]")
for i in 1:stride:length(S)
    q = summary(S[i])
    @printf("%-11.4e %-10.3e %-9.6f %-10.3e %-9.6f %-7.0f %-7.4f %-7.4f %-10.3e %-10.3e %-10.3e\n", q.t, q.dt, q.rsh, q.Tpk, q.rph, q.Teff,
            q.Lph, q.Lout, q.Mdot, q.Mdsh, q.dM)
end
