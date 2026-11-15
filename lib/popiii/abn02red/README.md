# ABN02 Fig. 2: the second output ("nine Myrs later"), from the vector paths

This is the same pipeline as `lib/popiii/abn02fin`: method, axis calibration and residuals are in `../abn02fin/README.md`,
and the axis maps in `../abn02fin/calib.json`. Here the curves are the **red dotted lines with red circles** (stroke
rgb(100%,0,0), dasharray 6.66,13.32). The caption's second time is "nine Myrs later (dotted lines with circles)"; the
other red line, 1.5e3 yr, is solid. These are SVG paths 981 (A), 1949 (B), 1210 and 1254 (C), 1756 (D) and 1529 (E).
Every red circle marker lies on its polyline within 0.003 pt. The figure has no DM line and no −c_s at this time.

```
python3 -I scripts/abn02_fig2/extract.py p4.svg bbox.html <final outdir> lib/popiii/abn02red   # final CSV unchanged (byte-identical)
python3 -I scripts/abn02_fig2/overlay.py p4r300-4.png 300 lib/popiii/abn02fin/calib.json lib/popiii/abn02red/abn02_red_fig2.csv \
        lib/popiii/abn02red/abn02_red_fig2_overlay.png "second output (9 Myr later, red dotted)"
python3 -I scripts/abn02_fig2/crosscheck_red.py lib/popiii/abn02red/abn02_red_fig2.csv
```

`abn02_red_fig2.csv` uses the same long format and curve names as the final CSV (`n_r`, `M_r`, `fH2_M`, `fHI_M`, `T_M`, `vr_M`).

## Profile

- **Central plateau.** n = 10^5.811 = 6.5e5 cm⁻³, exactly flat in the path data from r = 1e-6 pc to the vertex at
  log r = −2.034 (r_p = 9.2e-3 pc). The next vertex is 0.013 dex lower. Panel B gives M(r_p) = 0.0905 M☉.
  T = 253 K is flat within 0.01 dex to M = 13 M☉. f_H2 = 1.40e-3 is flat within 0.01 dex to 0.6 M☉.
  f_HI = 0.751 is constant to 5.8e3 M☉.
- **Radii.** M = 1 M☉ at 0.021 pc, 10 M☉ at 0.048 pc, 100 M☉ at 0.13 pc, 1e3 M☉ at 1.79 pc.
  1e9 R☉ is 22.55 pc, where panel B gives 7.14e3 M☉ (pchip; 7.16e3 linear).
- **Ranges.** Panel B covers M = 10^-2.351 … 10^4.632 M☉ (r = 10^-2.498 … 201 pc). Panels C, D and E cover 0.0102 … 2.0e4 M☉.
- **Temperature, composition, velocity.** T runs 165–1118 K. f_H2 falls from 1.4e-3 to 1e-4.
  v_r is −0.015 km/s at 0.01 M☉, −0.06 at 1 M☉, −0.38 at 100 M☉, with a minimum of −1.79 km/s at 3.3e3 M☉.
  The flow is subsonic everywhere: the isothermal c_s is 1.3 km/s at 253 K.

## M(r): panel A integrated against panel B

| μ | 10^-2.5 pc | 1e-2 | 10^-1.5 | 0.1 | 1 | 10 | 100 | mean / rms (dex) |
|---|---|---|---|---|---|---|---|---|
| 1.22 | −0.24 | −0.13 | −0.14 | −0.09 | −0.04 | −0.02 | −0.03 | −0.07 / 0.08 |
| μ(M) from panel C (1.229 at the centre) | −0.23 | −0.13 | −0.13 | −0.09 | −0.04 | −0.02 | −0.03 | −0.07 / 0.08 |
| 1 | −0.32 | −0.22 | −0.22 | −0.18 | −0.13 | −0.11 | −0.12 | −0.16 / 0.16 |

The gas here is atomic, so μ ≈ 1.23 and the choice of μ hardly matters. Integrated panel A still falls 0.13 dex below
panel B inside 0.1 pc, and 0.02–0.04 dex below it outside. A uniform sphere of the plateau density n_A μ m_H out to r_p
holds 0.065 M☉; panel B gives 0.0905 M☉ there.

## IC=abn02red (julia/popiii_m2c.jl; T. Abel, 15 Nov 2026)

The PI asked for this IC because the outer profile has stopped changing (outside about 0.1 pc it coincides with the final
output) and because starting from the final output gives a strong initial transient. It uses the abn02fin code, with
these differences:

- **Core.** Density, temperature and composition are constant from r_p to r = 0. r_p is the end of the exactly flat part
  of the n(r) path. The core density is panel B's mean density inside r_p, 1.85e-18 g cm⁻³, so the mass stays
  continuous with panel B at r_p. That is n = 10^5.95 at μ = 1.229, 0.14 dex above panel A's plateau (see the M(r)
  table). T, f_H2, f_HI and x_e inside are their values at M(r_p); T(M) and f(M) are flat there anyway.
  v_r is homologous (∝ r), from v_r(M(r_p)) at r_p to 0 at the centre.
- **Mass and tables.** M(r) comes from panel B. T(M) from D, f_H2 and f_HI from C, v_r(M) from E.
- **x_e.** The abn02fin recipe: 2e-4 in the outer cloud (T. Abel), recombined for one local free-fall time at the IC
  density and T, non-decreasing outward. XEFLAT=1 gives the constant. He is neutral.
- **Dark matter: `DMLINE=z19` (default) or `final`.** The figure shows DM only at z = 19 and at the final time. In time
  the final line is the closer one: this output is 9 Myr after z = 19, and the final output follows about 0.33 Myr later.
  But inside 0.1 pc the final-time DM has contracted along with the collapsed baryons. Inside r_p it reaches
  10^7.3–7.56 GeV cm⁻³, 20–35 times the core gas density, which is 3.2 M☉ of DM against 0.09 M☉ of gas. Run abnred6
  (DMLINE=final) shows the consequence: the core collapses on the DM free-fall time, 8e3 yr, against 4.9e4 yr for the
  gas. By t = 1.4e4 yr (step 15) n_c has gone from 10^5.81 to 10^6.57 and T_c from 254 to 505 K, nearly adiabatically.
  The default is therefore the z = 19 line (blue solid, SVG path 980, `dm_z19_r` in the CSV). It starts at 0.087 pc
  at 10^4.38 GeV cm⁻³, about 2 % of the core gas density. Between 0.1 and 1 pc it lies 0.3–0.5 dex below the final
  line; further out the two agree. Inside 0.087 pc it is held at 10^4.38, the pchip's constant continuation, which gives
  a uniform DM core with no cusp.
- **DMCORE** (default on for abn02red only). Inside the DM table's first radius, 1e-4 pc, the DM mass falls as r³.
  Without it, EXTMASS returns the constant M_DM(1e-4 pc) there, which acts as a central point mass. That is 4e-6 M☉
  for the final-time line and 2e-8 M☉ at z = 19. Run abnred5, with the ENVLIN mesh reaching 0.01 R☉, shows the effect:
  zone 1 reached −6.6 km/s and 2400 K in two steps, and step 5 then failed 60 times. abn02fin and abn02z19 still have the
  point mass, so that DMCORE=0 stays bit-identical.
- **MCLOUD.** The default is 2e4 M☉, the upper end of panels C, D and E. The smoke test uses 7137 M☉, panel B's mass
  inside 1e9 R☉.

## Comparison with the hand read-off IC=abn02 (same profile, "9 Myr later")

Each entry is hand minus vector.

- **log n(r)** at log r = −3, −1.57, −1.28, −1.08, −0.84, −0.59, −0.35, 0.06 and 0.87: +0.02, +0.05, −0.28, −0.28,
  −0.43, −0.56, −0.43, 0.00 and −0.10 dex. The hand table puts the plateau edge and the steep part about 0.2 dex in r
  too far inside.
- **log T(M):** within 0.02 dex, except +0.095 at 10^3.4 M☉ in the steep rise.
- **log f_H2:** −0.02 to −0.22 dex. The hand values are systematically low.
- **v_r:** −0.38 to +0.34 km/s, with the largest differences in the 10^2.5–10^3.7 M☉ wiggles.

## Smoke test (runs/abnred*, 15 Nov 2026)

Settings: the z300acc stack (run.env and phaseB.env) with `IC=abn02red MCLOUD=7137`, the `PSWTMIN=0 PSWITCH=1` line
dropped, NZ=128, 8 threads. The mesh is ENVLIN=10 linear zones from 0 to 0.1 R☉ plus log zones out to 1e9 R☉
(22.54 pc), 300 zones in all.

- **abnred5 (DM point mass, before DMCORE).** Zone 1 reached −6.6 km/s and 2400 K within two steps of about 400 s;
  step 5 then failed 60 times.
- **abnred6 (DMLINE=final, DMCORE).** The core collapses on the DM free-fall time: n_c went from 10^5.81 to 10^6.57
  and T_c from 254 to 505 K by t = 1.4e4 yr. Stopped at step 16.
- **abnred7 and abnred8 (defaults: DMLINE=z19, DMCORE; abnred8 adds PRINTEVERY=10).** 40 steps to t = 8.7e4 yr,
  which is 1.8 central t_ff0. There were 2 Newton failures (steps 12 and 40), each retried at dt/4, and no decoupled
  zones. dt grew from 1.5e9 s to 3.9e11 s. Each step took 62–76 s, of which GC was 11–15 s, and allocated
  150–230 GB. z300acc took 7.3 s per step with 8 GB.

  Profiles compared with the IC (log ρ, log T):

  | t | core (r < r_p) | r_p … 0.1 pc | beyond 0.1 pc |
  |---|---|---|---|
  | 7.2 kyr | +0.09, +0.06 | | ≤ 0.04 |
  | 13 kyr | −0.08, −0.06 | +0.11 at the plateau edge | ≤ 0.04 |
  | 31 kyr | −0.02, −0.04 | +0.14 | ≤ 0.05 |
  | 87 kyr | +0.38, +0.09 (n_c 2.2e6, T_c 315 K) | +0.32 | mean −0.02 (max 0.16 for 0.1–1 pc) |

  The first two rows are a ±0.09 dex adiabatic "breathing" of the uniform core on its sound-crossing time,
  r_p/c_s ≈ 7 kyr. The density step at r_p starts it: the core is panel B's mean density, about 15 % above the
  envelope just outside. After that the collapse continues; ABN02's next output, 0.3 Myr later, has n_c = 10^7.7.
  There is no strong initial transient.
