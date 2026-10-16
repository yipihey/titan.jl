//! Numerical core of the report: a steady quasi-star envelope with infall, ported from the Julia
//! package `QuasiStarInfall.jl` (which remains the specification; `tests/` pins this port to it).
//!
//! * [`constants`]: cgs constants and composition (X = 0.76, Y = 0.24, Z = 0).
//! * [`eos`]: Saha ionization (eq-saha), pressure (eq-pressure), entropy (eq-entropy) and the exact
//!   thermodynamic derivatives (eq-thermo).
//! * [`opacity`]: the tabulated Rosseland (eq-rosseland) and Planck means.
//! * [`envelope`]: structure equations (eq-structure), MLT (eq-mlt-cubic, eq-mlt-grad) and the
//!   Dormand–Prince inward integration.
//! * [`boundary`]: photospheric (eq-bc-phot) and infall (eq-bc-infall, eq-lucy) boundaries and the
//!   shooting on R.
//!
//! Every public item names the equation of the text it implements; the chapters link here with
//! the `{{< src >}}` shortcode.

pub mod boundary;
pub mod constants;
mod dual;
pub mod envelope;
pub mod eos;
pub mod opacity;

#[cfg(any(target_arch = "wasm32", doc))]
pub mod wasm;
