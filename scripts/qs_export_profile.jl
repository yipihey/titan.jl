# Copy of QuasiStarInfall/scripts/export_profile.jl (run it from that project: it needs its src/ and opacity table).
# Export a steady quasi-star envelope (hydrostatic, given Ṁ boundary) as CSV for TITAN initial states.
#   julia --project=. scripts/export_profile.jl fBH Mdot_Msun_yr out.csv
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Printf
fBH = parse(Float64, ARGS[1]); md = parse(Float64, ARGS[2]); out = ARGS[3]
tab = load_table(); M = 1e6Msun; L = LEdd_es(M)
p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
s = solve_R(p; Mdot = md * Msun / yr, lRmin = 15.0, lRmax = 16.8, nscan = 37)
pr = s.prof
open(out, "w") do io
    println(io, "# M=$(M) MBH=$(fBH*M) L=$(L) R=$(s.R) Mdot=$(md*Msun/yr) Tb=$(s.bc.Tb) Pb=$(s.bc.Pb) alpha=$(p.αMLT)")
    println(io, "r,m,P,T,rho,nabla,nabla_ad,nabla_rad,v,cs,kappa,Gamma,beta,HP")
    for i in eachindex(pr[:r])
        println(io, join((@sprintf("%.10e", pr[k][i]) for k in (:r,:m,:P,:T,:ρ,:∇,:∇ad,:∇rad,:v,:cs,:κ,:Γ,:β,:HP)), ","))
    end
end
@printf("R = %.5e  points = %d  Tb = %.1f\n", s.R, length(pr[:r]), s.bc.Tb)
