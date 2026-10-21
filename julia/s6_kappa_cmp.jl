# S6a: Rosseland opacity of the two Z = 0 tables in the photospheric regime of the protostar (StellarOpacity: S3; Z0Physics table: P7), X = 0.75.
#   TITAN_NOLIB=1 julia --project=env s6_kappa_cmp.jl
include(joinpath(@__DIR__, "s4_common.jl"))
tab = Z0.table()
println("   T[K]   rho=1e-10: stellar  z0        rho=1e-9: stellar  z0        rho=1e-8: stellar  z0")
for T in (3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 9000, 1e4, 1.2e4, 1.5e4, 2e4)
    row = String[]
    for ρ in (1e-10, 1e-9, 1e-8)
        push!(row, @sprintf("%9.3e %9.3e", StellarOpacity.kappa_R(ρ, T, 0.75), Z0.kap(tab, ρ, T)))
    end
    println(@sprintf("%8.0f   ", T), join(row, "   "))
end
