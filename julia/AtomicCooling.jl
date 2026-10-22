"""
    AtomicCooling

Atomic cooling of metal-free gas — hydrogen and helium — for the Pop III collapse (the sink
boundary region and the protostellar precursor reach 6e3–4e5 K with no atomic cooling in the
network, POPIII_NOTES S6b).  Self-contained (no dependencies; the same fits as
EmissionKernels.jl, whose `cooling_rate_total` + `cooling_rate_helium` this reproduces).

Sources, per coefficient: Cen (1992, ApJS 78, 341) Table 1, as compiled by Anninos, Zhang,
Abel & Norman (1997, New Astron. 2, 209) Table 1 (AZAN97); collisional ionisation rates k1,
k3, k5 from Abel, Anninos, Zhang & Norman (1997, New Astron. 2, 181) Table 3 (Janev et al.
1987 fits); case-B recombination coefficients k2, k4, k6 from Hui & Gnedin (1997, MNRAS 292,
27); bremsstrahlung Gaunt factor as in the Enzo/grackle implementation of this set.
[read]: the fits are transcribed from EmissionKernels.jl (cooling_atomic.jl, rates_atomic.jl),
which cites the same primary sources; they have NOT been re-checked here against the papers.

Channels (erg cm⁻³ s⁻¹; physical number densities):
  H:  collisional excitation (Lyα)  ceHI·n_e·n_HI;  collisional ionisation ciHI·n_e·n_HI
      (= 13.6 eV × k1); recombination (kinetic energy of the captured electron, case B)
      reHII·n_e·n_HII; bremsstrahlung brem·n_e·n_HII.
  He: He I 2³S excitation ceHeI·n_e²·n_HeII; He II excitation ceHeII·n_e·n_HeII;
      ionisation ciHeI·n_e·n_HeI (24.59 eV × k3), ciHeII·n_e·n_HeII (54.42 eV × k5),
      ciHeIS·n_e²·n_HeII (2³S); recombination reHeII1 (radiative) + reHeII2 (dielectronic)
      ×n_e·n_HeII, reHeIII·n_e·n_HeIII; bremsstrahlung brem·n_e·(n_HeII + 4 n_HeIII).

**Energy bookkeeping — read before plugging into ChemNet/NetGas.**  The Cen/AZAN set is
written for codes whose energy variable is the THERMAL energy: the "collisional
ionisation cooling" ciX is the kinetic energy spent on ionisation, and the recombination
photon's χ is not an energy loss of the gas.  In ChemNet the ionisation energy is part of
u by construction (species carry E0), so there (i) ciHI/ciHe must be omitted — the
Newton step already removes χ from the thermal energy when the composition ionises — and
(ii) a radiative recombination removes χ + ⟨E_kin⟩ from u (the photon carries both),
i.e. the loss is χ·α_rad·n_e·n_X⁺ + reX.  `cooling_atomic(…; ionization_in_u = true)`
does exactly that for hydrogen; helium is inert in ChemNet today (no He⁺ species), so for
helium the thermal-energy convention is used in both modes.  In optically thick gas
(J_ν = B_ν) these radiative terms are balanced by absorption and the transfer, not this
function, decides the net exchange (M3c): use the adapter only where the lines/continuum
are thin, or pass `escape` factors.
"""
module AtomicCooling

const kB = 1.3806504e-16; const eV = 1.602176634e-12
const χH = 13.598eV; const χHeI = 24.587eV; const χHeII = 54.418eV
const _LOG_DHUGE = 69.07755278982137          # log(1e30): guard of the Boltzmann exponents
const yHe = 1 / 12                            # n_He/n_H (RHFC02 App. A1.2, as PrimordialChem)
const mH = 1.6735575e-24
const μH = mH * (1 + 4yHe)

_expg(x) = exp(-min(_LOG_DHUGE, x))

# ------------------------------------------------------------------ collisional excitation (Cen 1992)
"H I collisional excitation (Lyα) [erg cm³ s⁻¹] × n_e n_HI."
ceHI(T) = 7.5e-19 * _expg(118348 / T) / (1 + sqrt(T / 1e5))
"He I 2³S collisional excitation [erg cm⁶ s⁻¹] × n_e² n_HeII."
ceHeI(T) = 9.1e-27 * _expg(13179 / T) * T^-0.1687 / (1 + sqrt(T / 1e5))
"He II collisional excitation [erg cm³ s⁻¹] × n_e n_HeII."
ceHeII(T) = 5.54e-17 * _expg(473638 / T) * T^-0.3970 / (1 + sqrt(T / 1e5))

# ------------------------------------------------------------------ collisional ionisation rates (AAZN97 Table 3)
_poly(x, c) = (s = 0.0; for i in length(c):-1:1; s = s * x + c[i]; end; s)
"H + e → H⁺ + 2e [cm³ s⁻¹] (Janev et al. 1987 fit; T in eV inside)."
k1(T) = T > 0 ? exp(_poly(log(T / 11605.0), (-32.71396786375, 13.53655609057, -5.739328757388, 1.563154982022,
    -0.2877056004391, 0.03482559773736999, -0.00263197617559, 0.0001119543953861, -2.039149852002e-6))) : 0.0
"He + e → He⁺ + 2e [cm³ s⁻¹]; the fit is used above 0.8 eV only (AAZN97)."
k3(T) = T / 11605.0 <= 0.8 ? 0.0 : exp(_poly(log(T / 11605.0), (-44.09864886561001, 23.91596563469, -10.75323019821,
    3.058038757198, -0.5685118909884001, 0.06795391233790001, -0.005009056101857001, 0.0002067236157507, -3.649161410833e-6)))
"He⁺ + e → He²⁺ + 2e [cm³ s⁻¹]; above 0.8 eV only."
k5(T) = T / 11605.0 <= 0.8 ? 0.0 : exp(_poly(log(T / 11605.0), (-68.71040990212001, 43.93347632635, -18.48066993568,
    4.701626486759002, -0.7692466334492, 0.08113042097303, -0.005324020628287001, 0.0001975705312221, -3.165581065665e-6)))
"Collisional ionisation cooling [erg cm³ s⁻¹]: χ × k (thermal-energy convention)."
ciHI(T) = 2.18e-11 * k1(T)
ciHeI(T) = 3.94e-11 * k3(T)
ciHeII(T) = 8.72e-11 * k5(T)
"He I 2³S collisional ionisation [erg cm⁶ s⁻¹] × n_e² n_HeII."
ciHeIS(T) = 5.01e-27 * T^-0.1687 / (1 + sqrt(T / 1e5)) * _expg(55338 / T)

# ------------------------------------------------------------------ recombination (Hui & Gnedin 1997, case B)
"H⁺ + e → H + γ, case B [cm³ s⁻¹]."
k2(T) = 4.881357e-6 * T^-1.5 * (1 + 1.14813e2 * T^-0.407)^-2.242
"He⁺ + e → He + γ, case B [cm³ s⁻¹] (radiative part)."
k4(T) = 1.26e-14 * (5.7067e5 / T)^0.75
"He²⁺ + e → He⁺ + γ, case B [cm³ s⁻¹]."
k6(T) = 7.8155e-5 * T^-1.5 * (1 + 2.0189e2 * T^-0.407)^-2.242
"He⁺ dielectronic recombination [cm³ s⁻¹]: the Aldrovandi & Péquignot (1973) form as used by Cen (1992) [memory: the
coefficient 1.9e-3 and the exponents are from memory — verify against AAZN97 Table 3 (k4) before production use]."
k4d(T) = 1.9e-3 * T^-1.5 * _expg(470000 / T) * (1 + 0.3 * _expg(94000 / T))
"Recombination cooling [erg cm³ s⁻¹] (kinetic energy of the captured electron), case B."
function reHII(T); λ = 2 * 157807 / T; 3.435e-30 * T * λ^1.970 / (1 + (λ / 2.25)^0.376)^3.720; end
function reHeII1(T); λ = 2 * 285335 / T; 1.26e-14 * kB * T * λ^0.75; end
reHeII2(T) = 1.24e-13 * T^-1.5 * _expg(470000 / T) * (1 + 0.3 * _expg(94000 / T))
function reHeIII(T); λ = 2 * 631515 / T; 8 * 3.435e-30 * T * λ^1.970 / (1 + (λ / 2.25)^0.376)^3.720; end

# ------------------------------------------------------------------ bremsstrahlung
"Free–free [erg cm³ s⁻¹] × n_e (n_HII + n_HeII + 4 n_HeIII); Gaunt factor 1.1 + 0.34 e^{−(5.5 − log T)²/3}."
brem(T) = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T))^2 / 3))
"Free–free for ions of charge Z: 1.43e-27 Z² √T ḡ(T/Z²) (the integrated Gaunt factor depends on kT/Z² Ry)."
brem_Z(T, Z) = 1.43e-27 * Z^2 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T / Z^2))^2 / 3))

# ------------------------------------------------------------------ helium ionisation state
"""
Helium stages (x_HeI, x_HeII, x_HeIII as fractions of He) in collisional ionisation equilibrium
(`:cie`, independent of n_e) or Saha (`:saha`, with the electron density `ne`; the collisional
terms are added so that the ratio is Saha + coronal, as in a J = B environment).
"""
function he_fractions(T, ne; mode = :cie)
    r1 = k3(T) / (k4(T) + k4d(T)); r2 = k5(T) / k6(T)
    if mode === :saha
        nQ = (2π * 9.1093837e-28 * kB * T / 6.62607015e-27^2)^1.5
        r1 += 4 * nQ * exp(-χHeI / (kB * T)) / max(ne, 1e-300)
        r2 += nQ * exp(-χHeII / (kB * T)) / max(ne, 1e-300)
    end
    d = 1 + r1 + r1 * r2
    (1 / d, r1 / d, r1 * r2 / d)
end

# ------------------------------------------------------------------ the cooling function
"""
    channels(T; nHI, nHII, nHeI, nHeII, nHeIII, ne) -> NamedTuple

All atomic channels [erg cm⁻³ s⁻¹] in the thermal-energy convention (Cen/AZAN), and the
totals `H`, `He`, `total`.
"""
function channels(T; nHI, nHII, nHeI = 0.0, nHeII = 0.0, nHeIII = 0.0, ne)
    H_exc = ceHI(T) * ne * nHI
    H_ci  = ciHI(T) * ne * nHI
    H_rec = reHII(T) * ne * nHII
    He_exc = (ceHeI(T) * ne + ceHeII(T)) * ne * nHeII
    He_ci  = ciHeI(T) * ne * nHeI + (ciHeII(T) + ciHeIS(T) * ne) * ne * nHeII
    He_rec = (reHeII1(T) + reHeII2(T)) * ne * nHeII + reHeIII(T) * ne * nHeIII
    ff = brem(T) * ne * (nHII + nHeII) + brem_Z(T, 2) * ne * nHeIII
    H = H_exc + H_ci + H_rec; He = He_exc + He_ci + He_rec
    (; H_exc, H_ci, H_rec, He_exc, He_ci, He_rec, brem = ff, H, He, total = H + He + ff)
end

"""
    cooling_atomic(nH, T, d; he = :cie, ionization_in_u = true) -> Λ [erg cm⁻³ s⁻¹]

Adapter for ChemNet/NetGas: `d` is the composition Dict of `ChemNet.full` (number fractions
x_i = n_i/n_H with keys :H, :Hp, :e, :He [, :Hep, :Hepp when a He network exists]).  Helium
stages from `d` if present, else `he_fractions(T, n_e; mode = he)` on the He budget `d[:He]`
(or yHe).  `ionization_in_u = true` (ChemNet: χ in u): no ciHI; hydrogen radiative
recombination removes χ_H + ⟨E_kin⟩ (k2·n_e·n_H⁺·χ_H + reHII); helium in the thermal
convention.  `false`: the plain Cen/AZAN sum (thermal-energy codes).
"""
function cooling_atomic(nH, T, d; he = :cie, ionization_in_u = true)
    ne = d[:e] * nH; nHI = d[:H] * nH; nHII = d[:Hp] * nH
    nHe = get(d, :He, yHe) * nH
    if haskey(d, :Hep)
        nHeII = d[:Hep] * nH; nHeIII = get(d, :Hepp, 0.0) * nH; nHeI = max(nHe - nHeII - nHeIII, 0.0)
    else
        f1, f2, f3 = he_fractions(T, ne; mode = he)
        nHeI, nHeII, nHeIII = f1 * nHe, f2 * nHe, f3 * nHe
    end
    c = channels(T; nHI, nHII, nHeI, nHeII, nHeIII, ne)
    if ionization_in_u
        return c.H_exc + c.H_rec + χH * k2(T) * ne * nHII + c.He + c.brem
    end
    c.total
end
"Λ(ρ, T, d) for `PopIIIMaterial.NetGas` (n_H = ρ/μ_H)."
cooling_atomic_rho(ρ, T, d; kw...) = cooling_atomic(ρ / μH, T, d; kw...)

end # module
