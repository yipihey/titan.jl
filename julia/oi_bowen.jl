# Order-of-magnitude O I λ8446 from Lyβ fluorescence (Bowen 1947) in the infall of a non-LTE solution, with the O I
# λ1025.77 oscillator strength f_O and the oxygen abundance as explicit parameters (the model itself has no metals).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia --project=env oi_bowen.jl NAME [fO=0.01] [OH=4.9e-4] [overlap=0.5]
# Lyβ is optically thick: every decay n=3 → 1 is a scattering, at a rate n₃A₃₁ per volume. O I absorbs a fraction
#   p = (n_OI f_O)/(n_HI f_Lyβ) × overlap
# of these flights (n_OI/n_HI ≈ O/H, both neutral together by charge exchange; overlap: the 0.05 Å = 15 km/s offset of
# the two lines against their Doppler widths), and each absorption gives one λ8446 photon (cascade 3d → 3p → 3s → 2p).
# Upper limit in spirit: it ignores Lyβ destruction competing with O I and O I absorption of other photons.
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "Spectrum.jl")); include(joinpath(@__DIR__, "NLTEH.jl"))
using Serialization, Printf
const Z0 = Z0Physics; const NH = NLTEH
name = ARGS[1]; opt = Dict(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in ARGS[2:end])
fO = get(opt, :fO, 0.01); OH = get(opt, :OH, 4.9e-4); ov = get(opt, :overlap, 0.5)
dn = deserialize(joinpath(@__DIR__, "..", "lib", "spectra", "nlte_$name.jls"))
ν31 = (NH.En[3] - NH.En[1]) / Z0.h; f31 = NH.fosc(1, 3)
A31 = 2Z0.h * ν31^3 / Z0.c^2 * (NH.gn[1] / NH.gn[3]) * 4π * NH.πe2mc * f31 / (Z0.h * ν31)
ν32 = (NH.En[3] - NH.En[2]) / Z0.h; f32 = NH.fosc(2, 3)
A32 = 2Z0.h * ν32^3 / Z0.c^2 * (NH.gn[2] / NH.gn[3]) * 4π * NH.πe2mc * f32 / (Z0.h * ν32)
flights = sum(st.n[3] * A31 * dV for (st, dV) in zip(dn.nl, dn.dV))                     # Lyβ scatterings per s
Nha = sum(Z0.h * ν32 * st.β[2, 3] * st.n[3] * A32 * dV for (st, dV) in zip(dn.nl, dn.dV)) / (Z0.h * ν32)   # Hα photons escaping per s (gross)
p = OH * fO / f31 * ov
N8446 = flights * p; L8446 = N8446 * Z0.h * Z0.c / 8446e-8
@printf("%s: A31 = %.3e s⁻¹, f(Lyβ) = %.4f; Lyβ scatterings %.3e s⁻¹, Hα photons %.3e s⁻¹ (ratio %.1f)\n", name, A31, f31, flights, Nha, flights / Nha)
@printf("O/H = %.2e, f_O = %.3f, overlap %.2f → p = %.2e per flight; L(O I 8446) ≈ %.2e erg/s, photon ratio to Hα %.3e\n", OH, fO, ov, p, L8446, N8446 / Nha)
