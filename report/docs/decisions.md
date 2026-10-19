# Decisions and conventions

Record here what cannot be derived from the code or the git history: sign
conventions, clocks and units, normalisations, why a method was chosen over
another, and the numbers that were verified when a claim entered the text.
One dated entry per decision; link the chapter and the test.

## 2026-10-18 — a second book, `spectra/`

- The short book on the spectra of the fed quasi-star lives in `../spectra/` as its own Quarto book. CI renders it
  after this one and copies `spectra/_site` to `_site/spectra/`, so it is served at `…/titan.jl/spectra/` and its
  widgets load this book's laboratory code and data with relative paths (`../app/js/widgets.js`, `../app/data/`).
  Locally: render both, then `cp -r ../spectra/_site _site/spectra`.
- New widget `analytic` (`app/js/panels/analytic.js`, data `analytic.json` from `export_report_data.jl analytic`)
  evaluates the closed-form ionization of the infall (second book, sec-sm-closed); its state keys are av, amd, afu,
  atr, alg. The closed form is checked against the full non-LTE solution in `julia/closed_form.jl` (neutral column
  8.8e23 vs 8.0e23 cm⁻², Thomson depth 4.98 vs 5.06 on the same points), not by a Rust test.
- The simple models use the emergent spectrum of the TITAN star without infall (`star_t0`) as the stellar intensity,
  because its far-UV excess (about 700× a 7310 K blackbody at 1000 Å) is what ionizes the infall; with a blackbody the infall
  stays neutral. Both books flag that excess as the weak point (photosphere without line opacity).

## 2026-10-19 — the five tests (second book, chapter 6)

- New chapter `spectra/06-tests.qmd` (sec-tests) and widget `tests` (`app/js/panels/tests.js`, data `tests.json` from
  `export_report_data.jl tests`; state keys tv, tma, tb, tl, tc). Views: mc, mcincl, balmer, axisym, he, cool.
- The outflow is kinematic: polar cones at +v_ff with the spherical solution's densities, populations and source
  functions (`axisym_lines.jl`, `mc_escatter.jl theta=`). Checks: θ_o = 0 reproduces the spherical Hα (EW −810 vs −805 Å);
  θ_o = 90° is its exact mirror. The equivalent widths of mixed geometries are indicative only (cross-flow resonance
  coupling is not in the source functions).
- The cooler star C30 (f_BH = 0.01, R = 1.3477e16 cm) has no Unsöld–Lucy correction: the photosphere sits on the H
  recombination front and the correction diverged (`ltarget=rad` and `lrmin` tried). Its non-LTE solution keeps 59 % of
  TITAN's luminosity; the book states this and compares ratios and EWs. The closed form, with the bare C30 star, gives the
  same neutral column (3.08e24 vs 3.07e24 cm⁻²) without fitting.
- The O I λ1025.77 oscillator strength could not be verified from a reachable source; `oi_bowen.jl` keeps it as an
  explicit parameter (f_O = 0.01) and the book says L(λ8446) scales with it.
- Balmer series: levels above n = 12 bracketed by b_n = b_12 and LTE; lines above n = 40 bounded by `balmer_bound.jl`
  (strip between the edge and the n = 40 line at 0.18 of the continuum → 4.31, black → 4.63).
