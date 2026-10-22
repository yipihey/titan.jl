# S7c-1: E_s of the S6b precursor boundary condition at the stall state (M ~ 32 Msun, Mdot ~ 2.3e-2): scan over L* and R.
include(joinpath(@__DIR__, "s4_common.jl"))
include(joinpath(@__DIR__, "s6_star.jl"))
include(joinpath(@__DIR__, "s6_setup.jl"))
include(joinpath(@__DIR__, "s6b_setup.jl"))
include(joinpath(@__DIR__, "s6c_precursor.jl"))
M = parse(Float64, get(ENV, "M", "32.58")) * Msun; Md = parse(Float64, get(ENV, "MD", "2.3e-2")) * Msun / YR
for R_ in (350.0, 400.0, 450.0), L_ in (5e4, 1e5, 2e5, 3e5, 4e5, 6.8e5)
    R = R_ * Rsun; Lin = L_ * Lsun
    p = precursor_thick(R, M, Md, Lin)
    @printf("R=%4.0f L*=%8.2e Ts=%7.0f (Lucy %7.0f) tau=%8.3f hR/GMR=%.4f w=%.2f Linf/Lin=%.3f Tph=%.0f rph/R=%.2f\n", R_, L_, p.Ts, (p.Es_lucy/Z0Physics.arad)^0.25, p.tau, p.hR/(GRAV*M/R), p.w, p.Linf/(L_*Lsun+GRAV*M*Md/R), p.Tph, p.rph/R)
end
