# Julia → WebAssembly with WasmTarget.jl: an investigation (16 Oct 2026)

**Question.** Could the book's browser figures run *our Julia code itself* rather than the Rust port in `report/core`, so
that what readers run is literally the code that produced the results?

**Tool.** [WasmTarget.jl](https://github.com/GroupTherapyOrg/WasmTarget.jl) v0.5.3 (Apache-2.0, registered, "experimental").
It reads Julia's inferred IR (`code_typed`) and emits WasmGC bytecode directly: no LLVM, no runtime, no linear memory.
Needs Julia ≥ 1.12 (tested here with 1.12.1; our code otherwise runs on 1.11.5) and `wasm-tools` for validation. Runs in
Node 22 / current browsers (WasmGC).

**Target.** The steady quasi-star model, `steady/src` (Saha EOS with ForwardDiff derivatives, opacity table interpolation,
MLT, Dormand–Prince envelope integration, infall boundary, root scan). Everything below is **[run]** unless marked.

## Results

| piece (original source unless noted) | compiles | runs in Node | vs native Julia |
|---|---|---|---|
| Saha `n_e`, `P_gas`, entropy (Float64 path) | yes | yes | bit-identical |
| `rho_PT` (Newton polish with nested Duals) | after patch 1 | yes | bit-identical |
| `eos_rhoT` (hand-seeded `EOSTag` Duals) | **no** (see bug C) | – | – |
| `eos_rhoT_ref` (three `ForwardDiff.gradient`s) | yes | yes | bit-identical (∇_ad, Γ₁) |
| `thermo` (→ ∇_ad, δ, c_P, Γ₁) | after patch 2 | yes | bit-identical to the patched native code; ≈1e-11 from the original `eos_rhoT` path |
| opacity table (filled at run time from JS) + `kap` | yes | yes | bit-identical |
| `mlt`, `rhs`, `outer_bc`, `infall_tau` | yes | yes | bit-identical |
| `integrate_inward` | after patches 3–4 | **traps** ("illegal cast") | – |
| same integrator without the `save`/profile `Dict` path (`probe.jl`, `Probe3`) | yes | yes | bit-identical m(r_in) in 5 of 6 test radii; one differs by 1.5e-11 relative |
| scan + bisection for all roots (`Roots2`) | yes | yes | roots identical to all digits except one at 2e-12 (inside the noise of the shooting function, see below) |

Module size: 0.5 MB unoptimized (`wasm-opt` not run here; upstream reports ~85 % reduction). Compile time: 1.5–2.5 min.

**Speed** (one full root scan, 105 trial radii + bisections): native Julia 3.4–4.7 s, Wasm in Node 4.2–5.8 s. The Rust port
does the same in 0.36–0.77 s in Wasm. The gap to Rust is algorithmic, not Wasm: the Julia source bisects the charge balance
200 times per Saha call, and patch 2 replaces the single-pass dual EOS by three gradient calls (13 M allocations per scan
instead of 1 M). With those fixed in the Julia source, the Julia-compiled Wasm should come close to native Julia.

**Bit-identity.** Elementary functions are bit-identical (`exp`, `log`, `log1p`, `^` with real exponent, `sqrt`), and so are
the EOS, opacity, MLT and boundary values tested. Along a 100-step integration, one of six test radii gave m(r_in) differing
by 1.5e-11. The steady roots sit where m(r_in) − M_BH changes sign at the 1e-11 level (the shooting function is itself that
noisy near a root: native values at R = 8.26785524971/72/73/74 e15 alternate in sign), so a root can move by ~2e-12. That is
far below every tolerance used in the book. The source of the 1.5e-11 difference was not located.

## Source changes needed (all exact rewrites; applied as patches in `qsw.jl`, the original files are untouched)

1. `eos.jl`: drop `@inline` on `_fractions`. With it, the Newton polish of the charge balance (a closure called with nested
   Duals) does not compile (bug B).
2. `eos.jl` `thermo`: derivatives from `ForwardDiff.gradient` in (ln ρ, ln T) (`eos_rhoT_jac` in `qsw.jl`) instead of the
   hand-seeded `EOSTag` Duals of `eos_rhoT` (bug C). Equal to rounding, but slower.
3. `envelope.jl`: the Dormand–Prince stage loop indexes a heterogeneous tuple of tuples (`DP_a[s]`, mixed lengths, `Int` and
   `Float64`) with a run-time stage number, and copies `ks[1:s-1]`. That is dynamic dispatch plus a `memmove`; WasmTarget
   lowers neither. Replaced by padded `NTuple{6,Float64}` coefficients and an explicit loop with the same summation order.
   This also helps natively: 2.5× fewer allocations, same results bit for bit.
4. `envelope.jl`: `DP_b`, `DP_bs` (and `DP_c`) contain `Int` zeros; converted to `Float64` (exact).
5. `integrate_inward` still traps with the optional profile (`Dict{Symbol,Vector{Float64}}` behind `save = false`); a copy
   without it works. Not narrowed further.

None of these change results; 3–4 are improvements we may want in `steady/src` anyway.

## WasmTarget problems found (minimal reproductions in `repro/`)

- **A. Silent wrong result (most serious):** assigning to a global `Ref` (`CNT[] = 0`) is dropped when the global's initial
  value is baked into the module. `repro/ref_reset_miscompile.jl`: native `3 3 | 3 3`, Wasm `3 3 | 8 11`. Workaround: no
  mutable globals for state that is reset; pass state explicitly. This contradicts the package's "correct-or-loud" design
  goal and is worth reporting upstream.
- **B.** A function marked `@inline` and called with nested ForwardDiff Duals inside a closure fails to compile
  ("unsupported method", `repro/inline_nested_dual.jl`).
- **C.** Duals constructed by hand (`Dual{Tag}(v, p1, p2)`) passed into code that nests a further `ForwardDiff.derivative`
  fail to compile, whereas the same computation through `ForwardDiff.gradient` compiles (`repro/handseeded_duals.jl`).
- **D.** Compile order matters: `compile_multi([table_init, kap_R, …])` hits an internal `StackImbalanceError` while
  emitting `Matrix{Float64}` construction; other orders compile (`repro/order_stack_imbalance.jl`).
- **E.** Run-time "illegal cast" trap in `integrate_inward` (item 5 above).
- Constants: any global array that holds data at compile time is baked into the module (a 490 k-element table refuses to
  compile; reset globals before `compile_multi`). Large data must be loaded at run time (done here through an exported setter).

## TITAN (not attempted; assessment **[read, inferred]**)

The TITAN residual is a much harder target: state lives in `Dict{Symbol,Array}` (`TitanFortran.var`), which is not
type-stable; the Z0 material relies on nested Duals with a hand-made tag (bugs B/C); `BlockSolver` uses StaticArrays matrices
(`SMatrix` is out of WasmTarget's stated scope). A browser TITAN would need the state refactored into a concrete struct
first, which we would want for GPU work anyway.

## Verdict

- **Feasible for the steady model now**, with four exact source rewrites, a loader for the table, and care with globals.
  What runs in the browser would then be the Julia functions of `steady/src` (patched as above), not a re-implementation,
  and they agree with native Julia bit for bit in almost every quantity tested.
- **Not yet trustworthy enough to replace the tested Rust core unchecked:** one silent miscompile (A) and an unexplained
  trap (E) in an experimental 0.5.x compiler. If used, the same reference tests that pin the Rust port should pin the
  Julia-compiled Wasm (the Node harness here does that comparison).
- A reasonable path: make `steady/src` type-stable (items 3–5), keep the Rust core as the published engine, and add the
  Julia-compiled module as a second engine checked against it in CI, until WasmTarget matures.

## Reproduce

```sh
cd wasm
julia +1.12 --project=. -e 'using Pkg; Pkg.instantiate()'     # WasmTarget 0.5.3, ForwardDiff
julia +1.12 --project=. build_probe.jl      # native reference values, then out/probe.wasm (~2.5 min)
node run4.mjs                                # the same functions in Node (WasmGC), with timings
node run5.mjs                                # shooting values near a root, to compare with julia cmp.jl
```
Needs `wasm-tools` on the PATH, and the v2 opacity table at `steady/output/kaptable_Z0_v2.jls`
(`steady/scripts/build_kap.jl`).
