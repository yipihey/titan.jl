//! Physical constants (cgs) and the composition, identical to `constants.jl` and `eos.jl` of the
//! Julia code. Values are CODATA 2018 / IAU 2015 nominal; ionization energies in erg.

/// Gravitational constant (cm³ g⁻¹ s⁻²).
pub const G: f64 = 6.67430e-8;
/// Speed of light (cm s⁻¹).
pub const C: f64 = 2.99792458e10;
/// Boltzmann constant (erg K⁻¹).
pub const KB: f64 = 1.380649e-16;
/// Planck constant (erg s).
pub const H: f64 = 6.62607015e-27;
/// Electron mass (g).
pub const ME: f64 = 9.1093837e-28;
/// Hydrogen atom mass (g).
pub const MH: f64 = 1.6735575e-24;
/// Helium-4 atom mass (g).
pub const MHE: f64 = 6.6464731e-24;
/// Stefan–Boltzmann constant (erg cm⁻² s⁻¹ K⁻⁴).
pub const SIGMA_SB: f64 = 5.670374e-5;
/// Radiation constant a = 4σ/c (erg cm⁻³ K⁻⁴).
pub const ARAD: f64 = 4.0 * SIGMA_SB / C;
/// Thomson cross section (cm²).
pub const SIGMA_T: f64 = 6.6524587e-25;
/// Electron volt (erg).
pub const EV: f64 = 1.602176634e-12;
/// Solar mass (g).
pub const MSUN: f64 = 1.98847e33;
/// Solar radius (cm).
pub const RSUN: f64 = 6.957e10;
/// Solar luminosity (erg s⁻¹).
pub const LSUN: f64 = 3.828e33;
/// Julian year (s).
pub const YR: f64 = 3.15576e7;

/// Ionization energy of H I (erg).
pub const CHI_H: f64 = 13.598434 * EV;
/// First ionization energy of He (erg).
pub const CHI_HE1: f64 = 24.587389 * EV;
/// Second ionization energy of He (erg).
pub const CHI_HE2: f64 = 54.417765 * EV;
/// Binding energy of H⁻ (erg).
pub const CHI_HM: f64 = 0.754195 * EV;

/// Hydrogen mass fraction X (Z = 0; as in the MESA-QUEST inlist).
pub const X: f64 = 0.76;
/// Helium mass fraction Y.
pub const Y: f64 = 0.24;

/// Eddington luminosity for electron scattering, κ_es = 0.2(1 + X):
/// `L_Edd = 4π c G M / (0.2 (1 + X))` (erg s⁻¹, M in g). This is the default envelope luminosity.
pub fn l_edd(m: f64) -> f64 {
    4.0 * std::f64::consts::PI * C * G * m / (0.2 * (1.0 + X))
}
