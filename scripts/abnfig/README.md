# abnfig — ABN02 Fig. 2–style PDF figures of a run

`abnfig.py` makes a two-page, vector, US-Letter PDF for one run of the 1D Pop III protostar code, styled after
Fig. 2 of Abel, Bryan & Norman 2002 (Science 295, 93): profiles at several output times overlaid, shared axes
(r on the left, enclosed mass on the right), one colour and dash pattern per time.

```
python3 -I scripts/abnfig/abnfig.py RUN [--panel P] [--out PATH] [--times auto|t1,t2,...] [--xaxis abn|r|m] [--n N]
```

- `--panel P`: explorer export directory (default: the session panel directory hard-coded in the script). Reads
  `P/RUN/frames.json`, `P/RUN/frames.b64.txt`, `P/RUN/perf.csv` (if present) and `P/runs.json` (for the description).
- `--out`: output PDF (default `P/RUN/abnfig.pdf`).
- `--times`: `auto` (default), or a comma-separated list of t − t0 in yr; the nearest frames are used.
- `--xaxis`: `abn` (default; left column vs r, right column vs M_enc, as in ABN02), `r` or `m` (all panels vs one
  coordinate; the redundant M_enc(r) panel then becomes r(M)).
- `--n`: number of automatically chosen times (default 7, at most 8).

Uses only numpy + matplotlib (Anaconda python). Takes ~3–5 s per run.

## Page 1 — radial structure
Left column vs radius [R☉] (top axis in AU): **A** n_H = X_H ρ/m_H, **B** M_enc, **C** L = 4πr²F, **D** grey
Rosseland τ from the outer boundary (τ = 2/3 marked), **E** Ṁ = −4πr²ρv_r.
Right column vs enclosed mass [M☉]: **F** T_gas (thick) and T_rad (thin), **G** v_r (thick) and −c_s (thin),
**H** x(H₂), **I** x(H⁺), **J** x(H⁻) (per H nucleus).
On every curve a filled circle marks the protostar radius R* (the code's accretion shock) and an open square
the photosphere r_ph (τ_R = 2/3). Legend, marker key and a short caption sit below the panels.
c_s = (γkT/μm_H)^½ with γ = 5/3 and μ from x(H₂), x(H⁺) and neutral He — an estimate, since the run's EOS
(e.g. He ionisation) is not exported. Non-positive L or Ṁ are masked on the log axes.

## Page 2 — histories and summary
vs t − t0 (log): **A** M*, **B** T_eff (from L_ph and r_ph, and the run's own), **C** R* and r_ph, **D** Ṁ through
10⁴ R☉, **E** L_ph, L_nuc, L_lines (the latter two converted from erg/s), **F** −u_min, **G** n_c, **H** T_c.
Dots in the curve colours mark the page-1 times. Dense perf.csv columns (mstar, rstar, umin, log_nc, Tc) are used
when perf.csv has more rows than the frames. Below: a table at each plotted time (t − t0, M*, R*, r_ph, T_eff,
L_ph, Ṁ, log n_c, T_c, peak infall speed and its radius).

## Choice of times
- t0 = time of the first frame with a protostar: finite M_star > 0 **and** log n_c ≥ 20. The density condition is
  needed because the code's M_star (mass inside the first interface with v_r < −2 km/s) is also finite during
  the pre-protostellar collapse (e.g. M_star = 0.39 M☉ at n_c = 10^14.6 in abnfin6).
- Up to two pre-protostar frames (first, and the one farthest in profile/n_c) are added in greys, labelled by
  log n_c and their (negative) t − t0. A run without a protostar is labelled by t − t_first and log n_c.
- The protostar frames are picked by greedy farthest-point sampling starting from the first and last frame, with a
  distance that combines the RMS difference of the log ρ, log T and soft-log v_r profiles (on a common log r grid,
  so zone counts may differ), the difference in log M*, and the difference in log(t − t0). This gives roughly
  log spacing in time/M* while skipping near-duplicate profiles (e.g. stalled phases).
- Colours: matplotlib `plasma` from 0 to 0.8 (dark → light = early → late; colour-blind safe and monotone in
  greyscale), combined with eight distinct dash patterns so the curves remain separable in black-and-white print.

The title block gives the run key, the description from runs.json, mesh size (zone-count range), radial
extent, frame/step range, status, generation date and the commit of this script (the run's own code commit is
not in the export).
