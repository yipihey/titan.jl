include("../src/constants.jl"); include("../src/eos.jl")
let
t = thermo(1e17, 1e7); println("hot gas-dominated: ∇ad=$(round(t.∇ad,digits=4)) Γ1=$(round(t.Γ1,digits=4)) β=$(round(t.β,digits=4))  [expect ~0.4, 5/3]")
t = thermo(Prad(1e6)*1.001, 1e6); println("rad-dominated: ∇ad=$(round(t.∇ad,digits=4)) Γ1=$(round(t.Γ1,digits=4)) β=$(t.β)  [expect ~0.25, 4/3]")
# Eddington standard-model check: Γ1 = (32-24β-3β²)/(24-21β) for fully ionized gas+rad
for β in (0.5, 0.1, 0.01)
  T = 1e7; Pr = Prad(T); P = Pr/(1-β); t = thermo(P,T)
  Γth = (32-24β-3β^2)/(24-21β)
  println("β=$β: Γ1=$(round(t.Γ1,digits=5)) analytic=$(round(Γth,digits=5))")
end
t = thermo(1e3, 3000.0); println("neutral: ∇ad=$(round(t.∇ad,digits=4)) Γ1=$(round(t.Γ1,digits=4))")
for T in (8e3,1e4,1.3e4,2e4,4e4,1e5)
  P = 30.0 + Prad(T); t = thermo(P, T); s=saha(t.ρ,T); println("Pgas≈30 T=$T: ρ=$(round(t.ρ,sigdigits=3)) xH=$(round(s.xH,digits=3)) ∇ad=$(round(t.∇ad,digits=3)) Γ1=$(round(t.Γ1,digits=3)) β=$(round(t.β,digits=3))")
end
end
