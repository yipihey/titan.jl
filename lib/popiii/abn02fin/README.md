# ABN02 Fig. 2: the final output, from the vector paths

Source: Abel, Bryan & Norman 2002, Science 295, 93 (arXiv:astro-ph/0112088). Fig. 2 is on page 4 of the arXiv PDF.
The final output is "finally 200 years later (dotted line with circles) at z = 18.181164". This README covers
`abn02_final_fig2.csv`, `calib.json` and `abn02_final_fig2_overlay.png`. They are made by `scripts/abn02_fig2/`:

```
pdftocairo -svg -f 4 -l 4 abn02.pdf p4.svg ; pdftotext -bbox -f 4 -l 4 abn02.pdf bbox.html
python3 -I scripts/abn02_fig2/extract.py p4.svg bbox.html lib/popiii/abn02fin            # CSV + calib.json
pdftoppm -r 300 -f 4 -l 4 -png abn02.pdf p4r300
python3 -I scripts/abn02_fig2/overlay.py p4r300-4.png 300 lib/popiii/abn02fin/calib.json lib/popiii/abn02fin/abn02_final_fig2.csv lib/popiii/abn02fin/abn02_final_fig2_overlay.png
python3 -I scripts/abn02_fig2/crosscheck.py lib/popiii/abn02fin/abn02_final_fig2.csv     # M(r) check, headline numbers
```

## Method

1. **Paths.** The figure in the PDF is vector graphics. Every stroked `<path>` in the SVG is parsed (M/L/C commands)
   and mapped through its `transform` to page points (y measured downward), the same frame that `pdftotext -bbox` uses.
2. **Line styles.** The caption gives the order: z = 19 solid, 9 Myr dotted+circles, 0.3 Myr dashed, 3e4 yr long-dashed,
   3e3 yr dot-dashed, 1.5e3 yr solid, final dotted+circles. In the SVG this becomes blue solid; red dasharray
   6.66,13.32 with red circles; black 33.3,33.3; green 66.6,33.3; blue 66.6,33.3,20,33.3; red solid; and
   **black dasharray 6.66,13.32 with black circles**. The 9 Myr and final curves share the dotted+circles style but differ in
   colour, red versus black, and the central densities increase monotonically along this sequence.
   The final curves are these SVG paths (numbered by line of p4.svg):
   - A: 1037, n(r); 1113, the final-time DM line.
   - B: 2005, M_enc(r).
   - C: 1308, f_H2(M), and 1369, f_HI(M).
   - D: 1806, T(M).
   - E: 1579, v_r(M). The purple dotted path 1641 and the filled purple symbols are −c_s.

   Every circle marker lies on its polyline to within 0.003 pt.
3. **Axes.** The axis frames are the closed black rectangles. Each axis uses its major-tick path, the 11.65-wide
   segments that start on the frame edge. Tick labels are the `10` and superscript words in `bbox.html`. Each label is
   assigned to its nearest major tick, and log10(value), or the value on the linear v axis, is fitted as a linear function
   of the page coordinate. Two checks follow. Every major tick must map to an integer, and every minor tick to
   log10(k) for k = 2..9 within its decade.

| axis | labelled ticks | fit residual (dex) | major ticks max dev | minor ticks max dev |
|---|---|---|---|---|
| A, B x: r [pc] | 1e-6, 1e-4, 1e-2, 1, 1e2 | ≤ 0.0033 | 0.004 (9) | 0.006 (74) |
| A y: n [cm^-3] | 1e-2 … 1e12 (8 labels) | ≤ 0.0053 | 0.0065 (16) | 0.009 (144) |
| B y: M_enc [M☉] | 1e-3, 1e-1, 10, 1e3, 1e5 | < 1e-5 | 0.0054 (9) | 0.0057 (73) |
| C/D/E x: M_enc [M☉] | 0.1, 10, 1e3 | ≤ 5e-5 | 0.007 (7) | 0.007 (56) |
| C y: mass fraction | 1e-4, 1e-2, 1 | ≤ 0.0021 | 0.0037 (5) | 0.0037 (37) |
| D y: T [K] | 1e2, 1e3 (only two labels) | 0 (2 points) | – | 0.0014 (11 minor ticks verify) |
| E y: v [km/s], linear | −3, −2, −1, 0 | ≤ 0.0014 km/s | 0.0014 | – |

So the extraction is accurate to about 0.01 dex, which is the pdfTeX/cairo coordinate quantisation. The decorative a.u.
axis at the top of B is not consistent with the pc axis: its log(au/pc) offset varies by 0.02–0.08 dex across the panel
from the expected 5.314. It is not used.

## CSV

`abn02_final_fig2.csv` is in long format: `curve, x, y, page_x_pt, page_y_pt`, with one row per path vertex. The curves are:

| curve | x | y |
|---|---|---|
| `n_r` | log r [pc] | log n [cm⁻³] |
| `dm_r` | log r [pc] | log ρ_DM [GeV cm⁻³] |
| `M_r` | log r [pc] | log M_enc [M☉] |
| `fH2_M`, `fHI_M` | log M [M☉] | log mass fraction |
| `T_M` | log M [M☉] | log T [K] |
| `vr_M` | log M [M☉] | v_r [km/s] |
| `mcs_M` | log M [M☉] | −c_s [km/s] (dotted path) |
| `mcs_M_symbols` | log M [M☉] | −c_s [km/s] (filled symbols) |

The −c_s curve leaves the frame below −3 km/s for M ≈ 7–60 M☉ and at its right end. Those vertices sit at −3,
because the path is clipped there.

**Panel C.** The two curves are the mass fractions of molecular and of neutral atomic hydrogen. Their sum is
f_H2 + f_HI = 0.756 ± 0.005 at every mass from 0.01 to 2e4 M☉, which is the hydrogen mass fraction X. The z = 19
curve is a flat 0.76. The PI, T. Abel, confirmed in Nov 2026 that the second curve is neutral H. An earlier suggestion
that it is H⁺ was withdrawn. **The figure contains no ionisation fraction.** The PI gives the ionised fraction of the
outer cloud as x(H⁺) = 2e-4 (pers. comm., Nov 2026). Nothing in the figure constrains it in the dense core.

## The final-output profile (from `crosscheck.py`)

- **Centre.** n_c = 10^13.52 cm⁻³ at the innermost point, r = 1e-6 pc. n is 1 dex below n_c at r = 10^-4.25 pc
  (12 au), where M_enc = 0.29 M☉.
- **Radii.** M_enc = 1 M☉ at r = 2.0e-4 pc (41 au). M_enc = 1e3 M☉ at r = 1.69 pc. Panel B ends at
  M = 10^4.632 = 4.3e4 M☉ at r = 201 pc. Panels C, D and E cover M = 0.0102 … 2.0e4 M☉ (10^4.30).
- **Temperature.** T ranges from 164 K (at ~1e3 M☉) to 1118 K (at ~1e4 M☉). The core is at 530–570 K.
- **Composition.** f_H2 = 0.74 and f_HI = 0.013 in the core, so the core is fully molecular. Outside, f_H2 falls to
  1e-4 at 1.5e4 M☉.
- **Velocity.** v_r is at its minimum, −2.29 km/s, at M = 0.63 M☉. The infall is supersonic around 0.3–1 M☉:
  |v_r|/c_s reaches 1.22 at 0.5 M☉. It is 0.21, 0.55 and 0.98 at 0.01, 0.1 and 1 M☉.
- **Dark matter.** At the final time, log ρ_DM = 7.56 GeV cm⁻³ at 7e-3 pc and −0.71 at 200 pc.

## M(r) cross-check: panel A integrated against panel B

Panel A is integrated as ρ = n μ m_H. The table gives Δ log M = (integral − panel B) in dex.

| μ | 1e-5 pc | 1e-4 | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 | mean / rms over B |
|---|---|---|---|---|---|---|---|---|---|
| 1.22 (atomic; the abn02* default) | −0.41 | −0.27 | −0.12 | −0.06 | −0.03 | −0.04 | −0.03 | −0.04 | −0.12 / 0.18 |
| μ(M) = 1/(f_HI + f_H2/2 + Y/4) from panel C | −0.14 | −0.04 | −0.06 | −0.05 | −0.02 | −0.03 | −0.03 | −0.04 | −0.05 / 0.06 |
| 1 (n = ρ/m_H) | −0.49 | −0.36 | −0.21 | −0.14 | −0.11 | −0.12 | −0.12 | −0.12 | −0.21 / 0.25 |
| 2.27 (fully molecular) | −0.14 | −0.01 | +0.15 | +0.21 | +0.24 | +0.23 | +0.24 | +0.23 | +0.15 / 0.19 |

The best match is the composition-dependent μ: 2.26 in the core, 1.63 at 1 M☉ and 1.23 at 1e3 M☉. It is within
0.03–0.06 dex of panel B everywhere except the steep innermost point, where M ~ 1e-2 M☉ and panel B is poorly resolved.
A constant μ = 1.22 underestimates the core mass by a factor of 1.9–2.6. For that reason **IC=abn02fin takes M(r)
from panel B directly**, so no μ is assumed. Inside panel B's first point, r = 4.8e-6 pc and M = 1e-3 M☉, the density
is uniform. That density is ρ/m_H = 8.9e13 = 10^13.95, which is μ_c = 2.71 times panel A's n_c. A fully molecular
core would give 2.26. The 0.08 dex excess is the innermost-point uncertainty of panel B, whose curve rises 2 dex per
0.5 dex in r there.

## IC=abn02fin (julia/popiii_m2c.jl) and the smoke test

`IC=abn02fin` reads the CSV at run time, or `ABNFINCSV` if that is set. It builds the IC from these pieces:

- **Mass.** M(r) comes from panel B (pchip in log–log) and the gas mass is MCLOUD, default 2e4 M☉. 2e4 M☉ is the upper
  end of panels C, D and E (10^4.30). Panel B reaches 4.3e4 M☉, and beyond 2e4 the T, f and v tables are extrapolated.
- **Temperature.** T(M) from panel D.
- **Composition.** x_H2 = f_H2 / (2 (f_H2 + f_HI)) per H nucleus, from panel C.
- **Velocity.** v_r(M) from panel E, with v ∝ r inside 0.01 M☉.
- **Ionisation (15 Nov 2026, T. Abel).** 2e-4 (XE0) is the ionised fraction of the outer cloud. Inside, x_e is that value
  recombined for one local free-fall time at the IC's density and temperature, so the dense regions do not start with an
  artificially high x_e (and H₂ formation through H⁻): x = x0 / (1 + k2(T) n_H x0 t_ff), with t_ff = sqrt(3π/(32 G ρ)),
  k2 the network's H⁺ + e rate (AAZN97 case A), n_H = ρ/μ_H and ρ from panel B's dM/dr. The profile is made
  non-decreasing outward (only the innermost ~0.01 M☉ needed it). XEFLAT=1 restores the constant XE0 of the first runs
  (abnfin1–6, z100fin). He is neutral. Values in the IC (printed at start-up, run abnfin7):

  | M (M☉) | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 | 1e3 | 1e4 |
  |---|---|---|---|---|---|---|---|---|
  | x_e | 3.0e-11 | 3.3e-11 | 5.4e-11 | 4.5e-10 | 4.7e-9 | 4.8e-8 | 5.4e-6 | 1.3e-4 |

- **Dark matter.** TX.EXTMASS comes from the final-time DM line, as with DM=1 for the other abn02* ICs.

ICEUL=1 accepts it. All other ICs are unchanged: an IC=abn02z19 ICEUL start (5 steps) gives an identical checkpoint
to popiii-s6 79d4dc1.

Smoke test, z100ic3 settings with ICEUL=1, NZ=128 and 4 threads (runs/abnfin*). The fixed mesh has 96 zones from
0.1 R☉ to 63.9 pc.

- **Settings unchanged (abnfin4).** The run stops at step 95, at t = 0.65 yr and n_c = 10^14.7. Step 94 is accepted
  with T = 0 in zone 2, and step 95 then fails 60 times. The dense core stays in the network: with PSWTMIN=4500, no zone
  below 4500 K uses the equilibrium table.
- **With PSWTMIN=0 PSWITCH=1, table in the core (abnfin6).** 1500 steps run with 70 Newton failures and 40 error
  rejections. A second core forms at t = 6.69 yr (central t_ff = 5.5 yr), with n_c 10^23.3 and T_c 1.07e5 K by the end.
  M* reaches 0.43 M☉ at R* = 73 R☉ by t = 15.9 yr, and dt is about 2e4 s late in the run.
