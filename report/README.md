# Quasi-stars under infall — the book

Source of https://yipihey.github.io/titan.jl/, a Quarto book in the
[quadrant](https://github.com/yipihey/quadrant) format: equations link to the Rust functions that evaluate them
(`core/`, compiled to WebAssembly for the figures and the laboratory), quantitative claims link to the tests that pin
them, and every experiment opens in the browser laboratory (`app/`) from a permalink.

- `chapters/` — the text; `chapters/references.bib` (every entry verified against ADS).
- `core/` — `qs_core`, the steady quasi-star model (Saha EOS, tabulated opacities, MLT envelope, infall boundary),
  ported from the Julia reference implementation in `../steady/` and tested against values exported from it
  (`core/data/reference.json`, written by `../steady/scripts/export_reference.jl`).
- `app/` — the laboratory and chapter widgets; `app/data/` holds the opacity table and the TITAN run output replayed
  in the figures (written by `../julia/export_report_data.jl`).

Build (Rust with the wasm32 target, wasm-pack, Quarto ≥ 1.6):

```sh
make test          # cargo test --release: the validation suite behind the ✓ badges
make wasm api check site
python3 scripts/serve.py 8790 _site
```

The site is rebuilt and deployed by `.github/workflows/publish.yml` at the repository root on every push to `main`
that touches `report/`.
