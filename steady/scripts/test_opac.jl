include("../src/constants.jl"); include("../src/eos.jl"); include("../src/opacity.jl")
let
println("H− bf at 0.85 μm: ", σHm_bf(0.85), "  [expect ≈4e-17]")
println("H− ff at 1.6μm, 6300K: ", kHm_ff(1.6, 6300.0), " cm^4/dyn [Gray: ~ few e-26]")
println("κ_R high T (1e7 K, ρ=1e-8): ", kappa_ross(1e-8, 1e7), " [→ κ_es = 0.2(1+X) = $(0.2*(1+0.76)) plus ff]")
println("κ_R 1e8 K: ", kappa_ross(1e-8, 1e8))
for T in (4000.,5000.,6241.,8000.,1e4,1.5e4,2e4,3e4,5e4,1e5,3e5,1e6)
  println("T=$T  κR(ρ=1e-10)=$(round(kappa_ross(1e-10,T),sigdigits=3))  κR(ρ=1e-8)=$(round(kappa_ross(1e-8,T),sigdigits=3))")
end
end
