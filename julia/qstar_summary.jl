# Per-snapshot summary of a quasi-star infall run (snapshots of qstar_run.jl relax), shared by qstar_analyze.jl and
# export_report_data.jl. Needs Z0Physics (with the same opacity table as the run) loaded at top level.
module QSum
import ..Z0Physics
const Z0 = Z0Physics
const σ = 5.670374e-5; const Msun = 1.989e33; const yr = 3.156e7
"Mass on the grid (g)."
mass(s) = sum(4π / 3 * (s.r[k+1]^3 - s.r[k]^3) * s.d[k] for k in eachindex(s.d))
"""
summary(s, R0, L0, M0): shock radius r_sh/R0 (largest outward velocity drop) and the peak temperature within 3 zones of
it; photosphere r_ph/R0 (Rosseland τ = 2/3 integrated inward from the outer boundary), L there (/L0) and
T_eff = (L/4πσr_ph²)^¼; L_rad at the outer boundary (/L0); mass flux (M⊙/yr) at the outer boundary and 4 interfaces
outside the shock; mass gained on the grid (M⊙).
"""
function summary(s, R0, L0, M0, tab = Z0.table())
    nz = length(s.d)
    τ = 0.0; kph = 1
    for k in nz:-1:1
        τ += Z0.kap(tab, s.d[k], s.T[k]) * s.d[k] * (s.r[k+1] - s.r[k])
        τ >= 2 / 3 && (kph = k; break)
    end
    rph = s.r[kph+1]; Lph = s.Lrad[kph+1] + s.Lconv[kph+1]
    Teff = (Lph / (4π * σ * rph^2))^0.25
    du = [s.u[k+1] - s.u[k] for k in 1:length(s.u)-1]
    ks = argmin(du)
    Tpk = maximum(s.T[max(1, ks - 3):min(nz, ks + 3)])
    Mdot = -4π * s.r[end]^2 * s.d[end] * s.u[end]
    j = min(ks + 4, length(s.u) - 1)
    Mdsh = -4π * s.r[j]^2 * 0.5 * (s.d[j-1] + s.d[j]) * s.u[j]
    (t = s.t, dt = s.dt, rsh = s.r[ks+1] / R0, Tpk, rph = rph / R0, Teff, Lph = Lph / L0, Lout = s.Lrad[end] / L0,
     Mdot = Mdot / (Msun / yr), Mdsh = Mdsh / (Msun / yr), dM = (mass(s) - M0) / Msun)
end
end # module
