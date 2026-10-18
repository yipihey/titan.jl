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
