# Primordial-gas reaction rates and H₂ cooling: source-referenced extraction

Compiled 2026-10-19 for the RHFC02-faithful chemistry module. No code was written.

## Conventions

**Tags**
- **[paper]**: copied from the source. Every formula below was read from the PDF page image, not only from the pdftotext output, unless marked otherwise.
- **[interp]**: my reading, conversion or check.
- **[num]**: a numerical check I ran with `python3 -I /tmp/ratecheck/check.py`, a scratch script outside every project tree.
- **[OCR-unverified]**: taken from the hidden ADS OCR text layer of `PSS83/full.pdf` (see §0). It cannot be checked against an image.

**Page numbers** are the printed page numbers of the local PDFs:
- **AAZN97**: arXiv preprint pages 1–30.
- **GP98**: article pages 1–18. The A&A page is probably 402 + p [interp, not checked].
- **GA08**: arXiv pages 1–26.
- **RA04**: preprint pages 1–17.
- **RA05**: chapter pages 1–51.
- **RHFC02**: preprint pages 1–34.

**Temperature units**
- T is in K unless a formula says T_eV or T_e (eV).
- 1 eV = 11605 K, the conversion used by AAZN97 and the code.

**Source abbreviations**
- **AAZN97**: Abel, Anninos, Zhang & Norman 1997, New Astron. 2, 181. RHFC02 cites it as "AANZ98" (see §A.0).
- **GP98**: Galli & Palla 1998.
- **GA08**: Glover & Abel 2008.
- **RA04**: Ripamonti & Abel 2004.
- **RA05**: Ripamonti & Abel 2005.
- **RHFC02**: Ripamonti, Haardt, Ferrara & Colpi 2002.
- **PSS83**: Palla, Salpeter & Stahler 1983.

**Code**
- **CK**: `/Users/tabel/Projects/ChemistryKernels.jl/src`.
- **EK**: `/Users/tabel/Projects/EmissionKernels.jl/src`.

---

## 0. Status of PSS83

### 0a. Verified from the scan (2026-10-19)
Tom supplied the full ADS scan (pp. 632–641), kept in `~/Projects/popiii_lit/PSS83/scan/`. Table 1 on p. 634 was read
from a 300 dpi render (`table1.png`). The [OCR-unverified] readings below are now **[paper]**:

| # | reaction | rate (as printed) | reference |
|---|---|---|---|
| k1 | H⁺ + e⁻ → H + hν | 1.88×10⁻¹⁰ T^−0.644 cm³ s⁻¹ | Spitzer 1956; Hutchins 1976 |
| k2 | H + e⁻ → H⁻ + hν | 1.83×10⁻¹⁸ **T** cm³ s⁻¹ (exponent 1; RA05's T^0.88 is not PSS83's) | Doughty & Fraser 1964; Hutchins 1976 |
| k3 | H + H⁻ → H₂ + e⁻ | 1.35×10⁻⁹ cm³ s⁻¹ | Schmeltekopf et al. 1967 |
| k4 | H+H+H → H₂+H | 5.5×10⁻²⁹ T⁻¹ cm⁶ s⁻¹ | Jacobs et al. 1967 |
| k5 | H₂+H → 3H | 6.5×10⁻⁷ T^−1/2 exp(−52 000/T)[1 − exp(−6000/T)] **cm³** s⁻¹ (RA05's cm⁶ is a typo) | Jacobs et al. 1967 |
| k6, k7 | with H₂ as third body | k4/8, k5/8 | Jacobs et al. 1967 |
| k8 | H + e⁻ → H⁺ + 2e⁻ | 5.8×10⁻¹¹ T^1/2 exp(−158 000/T) cm³ s⁻¹ | Lotz 1967 |
| k9 | H + H → H⁺ + e⁻ + H | 1.7×10⁻⁴ k8 | Drawin 1969 |

PrimordialChem's k4_pss, k5_pss and k9_pss agree with the scan exactly. k1, k2, k3 are not used: AAZN97 k2, k7 and k8
are used instead.

Also confirmed on pp. 633–635:
- H₂ levels: v = 0–2, J ≤ 20, with ortho:para fixed at 3:1.
- Thick lines use the escape prescription of eq. (6).
- The basic model starts at T₀ = 100 K, n₀ = 2.9 cm⁻³, f_H₂ = 6.3×10⁻⁶ and f_e = 2.8×10⁻⁵.

The brief says PSS83 is not available. The local file `/Users/tabel/Projects/popiii_lit/PSS83/full.pdf` (identical to `try.pdf`) has **one image page**, ApJ 271, p. 632: abstract and introduction, with no rates.

However, its hidden text layer is the ADS OCR of **all of pp. 632–641**. `pdftotext -layout full.pdf -` shows Table 1 on p. 634.

The main tables below rely only on the secondary sources, as instructed. The OCR is used only as a labelled cross-check: it agrees with RA05 Table 1.1 wherever it is legible.

**[OCR-unverified]** PSS83 p. 634, Table 1, "Reactions and Rate Constants" (garbled OCR, my reading):

| PSS83 # | reaction | OCR text → reading | ref. column |
|---|---|---|---|
| k1 | H⁺ + e → H + hν | "= 1.88X10“" (rest lost) | Spitzer 1956; Hutchins 1976 |
| k2 | H + e⁻ → H⁻ + hν | "k2 = 1.83 X 10~lsTcm3 s“" → 1.83×10⁻¹⁸ T^0.88 cm³ s⁻¹ | Doughty & Fraser 1964; Hutchins 1976 |
| k3 | H⁻ + H → H₂ + e | "1.35x 10 9 cm3 s 1" → 1.35×10⁻⁹ cm³ s⁻¹ | Schmeltekopf et al. 1967 |
| k4 | H+H+H → H₂+H | "5.5XlO" r- cm s~" with exponents "29 1 6 1" → 5.5×10⁻²⁹ T⁻¹ cm⁶ s⁻¹ | Jacobs et al. 1967 |
| k5 | H₂+H → 3H | "6.5 X 10~77’~1/2 exp(−52,000/7) X [1 − exp(−6000/7)] cm3 s~1" → 6.5×10⁻⁷ T^(−1/2) exp(−52000/T)[1 − exp(−6000/T)] **cm³ s⁻¹** | Jacobs et al. 1967 |
| k6, k7 | (H+H+H₂ → 2H₂; H₂+H₂ → 2H+H₂) | label "k_6" only; the values are lost | Jacobs et al. 1967 (×2) |
| k8 | H + e⁻ → H⁺ + 2e⁻ | "5.8X 10~ 7 / exp(−158,000/7) cm s“ 3 1" with "1, 1 2" → 5.8×10⁻¹¹ T^(1/2) exp(−158000/T) cm³ s⁻¹ | Lotz 1967 |
| k9 | H + H → H⁺ + e⁻ + H | "k9 = \JX 10 /c8" with "-4" → 1.7×10⁻⁴ k8 | Drawin 1969 |

**[OCR-unverified]** PSS83 p. 634 text:

> "The dissociation rates for H₂ were derived, using the principle of detailed balance, from the measured recombination rates of Jacobs, Giedt, and Cohen (1967). Over the relevant temperature range, our dissociation rates agree within a factor of 5 with the measurements of Breshears and Bird (1973)."

This is the formation/dissociation detailed-balance relation asked about in item C. It is not discussed in GA08. My numerical check is in §C.3.

**[OCR-unverified]** PSS83 p. 633: "We keep track of the energy lost or gained through dissociation, ionization, or recombination by appropriately defining ε for each species." This is an internal-energy (binding-energy) bookkeeping, the same as RHFC02 eq. A8. There is no n_cr suppression.

---

## A. RHFC02 Table A1 set ("faithful")

### A.0 What RHFC02 says [paper, RHFC02 p. 29, Appendix A1.2, Table A1]

```
Table A1. Chemical reactions and reaction rates
  H+ + e- -> H + γ        AANZ98, k2
  H + e-  -> H- + γ       AANZ98, k7
  H + H-  -> H2 + e-      AANZ98, k8
  H+H+H   -> H2 + H       PSS83, k4
  H2 + H  -> H+H+H        PSS83, k5
  H+H+H2  -> H2 + H2      PSS83, k6
  H2 + H2 -> H + H + H2   PSS83, k7
  H + e-  -> H+ + e- + e- AANZ98, k1
  H+H     -> H+ + H + e-  PSS83, k9
"AANZ98 stands for Abel, Anninos, Norman & Zhang (1998) and PSS83 stands for Palla, Salpeter & Stahler (1983).
Note that the H recombination rate (first reaction in the list) is appropriate for case A (case B could be
more appropriate, but results are not affected by the difference, as can be seen from the comparison to ON98)."
```

Further notes on RHFC02 [paper, p. 29–30]:
- The species list includes He⁰, He⁺ and He⁺⁺. The He number abundance is 1/12 of H (X = 0.75, Y = 0.25). However, **Table A1 lists no He reaction**, and no H₂⁺ reaction although H₂⁺ is listed as a species. This is an internal inconsistency in RHFC02.
- The text says "(see Table 1)"; it means Table A1.
- Above ρ_eq ≡ 5×10⁻¹¹ g cm⁻³ the network is switched off and equilibrium chemistry is used:
  - Saha equations for p ≤ 10⁴ dyn cm⁻²;
  - SCVH95 tables for p ≥ 10⁶;
  - interpolation between the two.
- Chemical energy enters through the internal energy, not as a heating rate (eqs. A6–A8, p. 30): `u = u_th + u_chem`, `u_chem = −N_tot Σ_i y_i χ_i` (χ_i = binding energy of species i).
- §3.2 (p. 9): "4.48 eV per dissociated molecule are lost by the gas thermal energy".

**Citation [interp].** The RHFC02 reference list reads "Abel T., Anninos P., Norman M. L., Zhang Y. 1998, ApJ, 518". This is presumably ApJ 508, 518, which uses the AAZN97 network. **The k-numbers k1, k2, k7, k8 match the AAZN97 appendix tables exactly**, so AAZN97 is the primary source.

**AAZN97 internal numbering trap [paper].** AAZN97 numbers its processes in two places, and the two schemes differ for processes 2–6:
- **§3 text** (pp. 6–7): (1) collisional ionization of H; (2) of He; (3) of He⁺; (4) radiative recombination to H; (5, 6) recombination to He.
- **Appendix tables** (pp. 24–26): (1) H ionization; **(2) H recombination**; (3) He ionization; (4) He⁺ recombination; (5) He⁺ ionization; (6) He⁺⁺ recombination.

RHFC02's "k2 = recombination" follows the appendix numbering. From (7) on, text and tables agree.

### A.1 Summary table

The detailed formulas follow in A.2–A.10.

| # | reaction | RHFC02 source | rate as given in the determining source | determined from | CK code (function) | code vs paper |
|---|---|---|---|---|---|---|
| 1 | H⁺+e→H+γ (case A) | AANZ98 k2 | AAZN97 k2 (A.2) | **primary** | `k2` = **case B**, Hui & Gnedin 97 | **differs**: case B is 0.81× case A at 100 K, 0.75× at 10³ K, 0.62× at 10⁴ K [num] |
| 2 | H+e→H⁻+γ | AANZ98 k7 | AAZN97 k7 (A.3) | **primary** | `k7` = Stancil et al. 98: 3.0e-16 (T/300)^0.95 e^(−T/9320) | different fit. Within 2% for T ≤ 10³ K; +4% at 3000 K; −19% at 10⁴ K [num] |
| 3 | H+H⁻→H₂+e | AANZ98 k8 | AAZN97 k8 (A.4) | **primary** | `k8` = Kreckel et al. 2010 | **differs**: 3.1× AAZN97 at 100 K, 1.8× at 10³ K, 0.56× at 10⁴ K [num] |
| 4 | H+e→H⁺+2e | AANZ98 k1 | AAZN97 k1 (A.5) | **primary** | `k1` same Janev fit (more digits) | equal. AAZN97 sets it to zero below 1000 K; CK evaluates the fit (≈ 0 anyway) |
| 5 | H+H+H→H₂+H | PSS83 k4 | 5.5×10⁻²⁹ T⁻¹ cm⁶ s⁻¹ (A.6) | **secondary only** (RA04 eq. 31, RA05 k8, OCR agrees) | `k22` = ABN02: 1.3e-32 (T/300)^−0.38 (T ≤ 300), (T/300)^−1 (T > 300) | **differs**: CK is 14.1× lower than PSS83 for T > 300 K; 14.0× at 300 K [num] |
| 6 | H₂+H→3H | PSS83 k5 | (a) original 6.5×10⁻⁷ T^(−1/2) e^(−52000/T)(1−e^(−6000/T)); (b) RA04-modified 2.2×10⁻⁹ T^0.2 e^(−51800/T)(1−e^(−6000/T)) cm³ s⁻¹ (A.7) | **secondary only** | `k13` = Dove & Mandy 86 via AAZN97 k13 (functional form ≠ AAZN97 print). **Zero for T ≤ 0.3 eV = 3481 K** | **differs by orders of magnitude**: below 3481 K CK is 0 where PSS83 is 5×10⁻²⁰ at 2000 K. At 5000 K CK/PSS83 = 1/450; at 8000 K 1/130 [num] |
| 7 | H+H+H₂→2H₂ | PSS83 k6 | k6 = k4/8 (A.8) | **secondary only** (RA05 k10 = k8/8; RA04 n_f; GA08 k31 = k30/8 "following PSS83") | **absent** | missing reaction |
| 8 | H₂+H₂→2H+H₂ | PSS83 k7 | k7 = k5/8 (A.9) | **secondary only** (RA05 k11 = k9/8; RA04 n_f) | **absent** | missing reaction |
| 9 | H+H→H⁺+H+e | PSS83 k9 | 1.7×10⁻⁴ × 5.8×10⁻¹¹ T^(1/2) e^(−158000/T) cm³ s⁻¹ (A.10) | **secondary only** (RA05 k13 = 1.7e-4 k12) | `k57` = Lenzuni et al. 91: 1.2e-17 T^1.2 e^(−1.578e5/T) | different fit: 0.35× at 3000 K, 0.49× at 5000 K, 0.78× at 10⁴ K, 1.26× at 2×10⁴ K [num] |

### A.2 H⁺ + e⁻ → H + γ, case A (AAZN97 k2) [paper, AAZN97 p. 24, Table 1 (2); text p. 7, §3.2 "(4)"]

```
k2 = exp[−28.6130338 − 0.72411256 ln(T) − 2.02604473×10⁻² ln(T)² − 2.38086188×10⁻³ ln(T)³
         − 3.21260521×10⁻⁴ ln(T)⁴ − 1.42150291×10⁻⁵ ln(T)⁵ + 4.98910892×10⁻⁶ ln(T)⁶
         + 5.75561414×10⁻⁷ ln(T)⁷ − 1.85676704×10⁻⁸ ln(T)⁸ − 3.07113524×10⁻⁹ ln(T)⁹]  cm³ s⁻¹
T in eV. Source: "Our fit to data given by Ferland et al. (1992)", the sum of α_n for n = 1..1000, i.e. case A.
Text: "accurate for temperatures from 1K to 10⁹ K".
```

**Other sources for case A**
- **GA08 k13,A** [paper, GA08 p. 24, Table A1, ref. 13 = Ferland et al. 1992]: `k13,A = 1.269×10⁻¹³ (315614/T)^1.503 [1.0 + (604625/T)^0.470]^−1.923`. Agrees with AAZN97 k2 to within 1–3% over 100–3×10⁴ K [num].
- **RA05 Table 1.1 k1** [paper, RA05 p. 22], attributed to T97: `k1 ≃ 1.88×10⁻¹⁰ T^−0.644 cm³ s⁻¹`. This is 12% above AAZN97 at 100 K, 20% at 10⁴ K and 33% at 3×10⁴ K [num]. PSS83 k1 is the same, 1.88×10⁻¹⁰ T^−0.644 [paper, scan §0a].
- **GP98 H1** [paper, GP98 p. 2]: GP98 uses the cosmological effective rate `R_c2 = 8.76×10⁻¹¹ (1+z)^−0.58 cm³ s⁻¹` (recombination to excited states, Grachev & Dubrovich 1991). It is not a function of T and is not a substitute.

### A.3 H + e⁻ → H⁻ + γ (AAZN97 k7) [paper, AAZN97 p. 25, Table 2 (7); text p. 7]

```
k7 = 1.429×10⁻¹⁸ T^0.7620 T^(0.1523 log10 T) T^(−3.274×10⁻² log10² T)  cm³ s⁻¹      for T ≤ 6000 K
   = 3.802×10⁻¹⁷ T^(0.1998 log10 T) dex(4.0415×10⁻⁵ log10⁶ T − 5.447×10⁻³ log10⁴ T) cm³ s⁻¹   otherwise
"This work from data by Wishart (1979)"
```

**[interp]** AAZN97's appendix says "temperatures are in eV unless stated otherwise". The branch condition here is in K, and **T is in K** in this formula. Evidence: it is identical (to < 0.01%) to GA08 k1 with T in K [num].

Text accuracy claims: "within a few percent for 1K < T < 100K, about one percent for 100K < T < 10⁴ K, and better than 10% for T > 10⁴ K".

**Other sources**
- **GA08 k1** [paper, GA08 p. 24, ref. 1 = Wishart 1979]: `k1 = dex[−17.845 + 0.762 log T + 0.1523 (log T)² − 0.03274 (log T)³]` for T ≤ 6000 K; `= dex[−16.4199 + 0.1998 (log T)² − 5.447×10⁻³ (log T)⁴ + 4.0415×10⁻⁵ (log T)⁶]` for T > 6000 K. Algebraically the same as AAZN97, since log10(1.429e-18) = −17.845.
- **GP98 H3** [paper, GP98 p. 3]: `1.4×10⁻¹⁸ Tg^0.928 exp(−Tg/16200)` (fit to de Jong 1972). Within about 6% below 3000 K and 10–12% above at 6000–10⁴ K [num].
- **RA05 k2** [paper, p. 22]: `1.83×10⁻¹⁸ T^0.88`. PSS83 k2 [paper, scan §0a] is `1.83×10⁻¹⁸ T` (exponent 1). Close below 1000 K; 1.4× at 6000 K and 1.7× at 10⁴ K [num].

### A.4 H + H⁻ → H₂ + e⁻ (AAZN97 k8) [paper, AAZN97 p. 25, Table 2 (8); text p. 7]

```
T > 0.1 eV : k8 = exp[−20.06913897 + 0.22898 ln(T) + 3.5998377×10⁻² ln(T)² − 4.55512×10⁻³ ln(T)³
                      − 3.10511544×10⁻⁴ ln(T)⁴ + 1.0732940×10⁻⁴ ln(T)⁵ − 8.36671960×10⁻⁶ ln(T)⁶
                      + 2.23830623×10⁻⁷ ln(T)⁷]  cm³ s⁻¹          (T in eV)
T < 0.1 eV : k8 = 1.428×10⁻⁹ cm³ s⁻¹
"Our Integration of data from Janev et al. (1987)"
```

The text says "about 10% higher than the constant values given by Bieniek (1980) and de Jong (1972)", and "Launay et al. (1991) ... 1.49×10⁻⁹ [at 300 K], 4% above".

**Other sources**
- **GP98 H5** [paper, p. 3]: 1.5×10⁻⁹ (Tg ≤ 300); 4.0×10⁻⁹ Tg^−0.17 (Tg > 300) (LDZ = Launay et al. 1991). GP98 p. 3 comments that "Abel et al. (1997) favor a rate increasing with temperature for T > 1000 K".
- **GA08 k2** [paper, p. 3, eqs. 5–7]: default 1.3×10⁻⁹ (Schmeltekopf et al. 1967). Bounds are 5.0×10⁻⁹ and 0.65×10⁻⁹ (Glover, Savin & Jappsen 2006).
- **RA05 k3** [paper, p. 22]: 1.3×10⁻⁹. **PSS83 k3** [paper, scan §0a]: 1.35×10⁻⁹.

### A.5 H + e⁻ → H⁺ + 2e⁻ (AAZN97 k1) [paper, AAZN97 p. 24, Table 1 (1); text p. 6]

```
k1 = exp[−32.71396786 + 13.536556 ln(T) − 5.73932875 ln(T)² + 1.56315498 ln(T)³ − 0.2877056 ln(T)⁴
         + 3.48255977×10⁻² ln(T)⁵ − 2.63197617×10⁻³ ln(T)⁶ + 1.11954395×10⁻⁴ ln(T)⁷
         − 2.03914985×10⁻⁶ ln(T)⁸]  cm³ s⁻¹       (T in eV; Janev et al. 1987, 2.1.5)
Text: "For temperatures below 1000K ∼ 0.086eV we take it to be zero since it is smaller than 10⁻⁶⁰ cm³ s⁻¹ ... The ionization threshold is 13.6eV."
```

**Other sources**
- **GA08 k12**: the same fit [paper, p. 24, ref. 12].
- **RA05 k12** [paper, p. 22] and PSS83 k8 [paper, scan §0a]: `5.8×10⁻¹¹ T^(1/2) e^(−158000/T) cm³ s⁻¹` (Lotz 1967). This is 0.19–0.003× AAZN97 at 3000–5000 K (both negligible there) and 1.1–1.35× at 10⁴–2×10⁴ K [num].

### A.6 H + H + H → H₂ + H (PSS83 k4)

| source | formula | location |
|---|---|---|
| RA04 eq. (31) [paper] | `k4 = 5.5 × 10⁻²⁹ T⁻¹ cm⁶ s⁻¹` ("We use the PSS83 formation rate (k4)") | RA04 p. 11, §3.2.2 |
| RA05 Table 1.1 k8 [paper] | `k8 ≃ 5.5 × 10⁻²⁹ T⁻¹ cm⁶ s⁻¹`; caption: reactions 8–13 from "Palla et al. 1983" | RA05 p. 22 |
| PSS83 k4 [OCR-unverified] | 5.5×10⁻²⁹ T⁻¹ cm⁶ s⁻¹, Jacobs et al. 1967 | PSS83 p. 634 |

There is no stated validity range. All three agree.

### A.7 H₂ + H → 3H (PSS83 k5); RA04's modification

| version | formula | location |
|---|---|---|
| **PSS83 original, as quoted by RA05** [paper] | `k9 ≃ 6.5 × 10⁻⁷ T^(−1/2) e^(−52000/T) × (1 − e^(−6000/T))` **printed "cm⁶ s⁻¹"** | RA05 p. 22, Table 1.1 k9 |
| PSS83 [OCR-unverified] | same constants; unit **cm³ s⁻¹** | PSS83 p. 634 |
| **RA04 modified** [paper] | `k5 = 2.2 × 10⁻⁹ T^0.2 e^(−B5/T) (1 − e^(−C5/T)) cm³ s⁻¹, B5 = 51800 K, C5 = 6000 K` | RA04 p. 11, eq. (32) |

RA04's justification [paper, RA04 p. 11]: "We use the PSS83 formation rate (k4), but we modify the dissociation rate (k5) in order to better approximate the 'high density' (n ≳ 10⁹ cm⁻³) results of Martin, Schwarz & Mandy (1996)."

**Discrepancies**
- **RA05 unit typo.** The cm⁶ s⁻¹ in RA05 is a typo: a two-body dissociation rate is cm³ s⁻¹. RA04 and the PSS83 OCR both have cm³ s⁻¹.
- **Size of RA04's modification [num].**
  - The ratio (RA04 k5)/(PSS83 k5) is 0.36 at 300 K, 0.52 at 1000 K, 0.76 at 2000 K, 0.98 at 3000 K, 1.37 at 5000 K and 1.87 at 8000 K.
  - So the modification is within a factor of 2 over 1000–8000 K.
- **Which version RHFC02 used.** RHFC02 attributes the rate to "PSS83, k5". The faithful choice is therefore (a), the original. RA04 eq. 32 is the version RA04 used in its own 1D test with the modified R02 code.

**Code: CK `k13`** (`rates_h2.jl`) = "Abel et al. (1997)" Dove & Mandy (1986):
```
Tev ≤ 0.3: 0
else      1.0670825e-10 · Tev^2.012 / ( exp(4.463/Tev) · (1 + 0.2472 Tev)^3.512 )
```

**AAZN97 k13 as printed** [paper, AAZN97 p. 26, Table 3 (13), Dove and Mandy 1986]:
```
k13 = 1.067 × 10⁻¹⁰ T^2.012 × exp(−(4.463/T)(1 + 0.2472 T)^3.512) cm³ s⁻¹        (T in eV)
```
- **Form mismatch.** AAZN97 puts (1+0.2472T)^3.512 inside the exponent; CK divides by it. These are different functions.
  - At 5000 K the AAZN97 print gives 7.5×10⁻¹⁸ and CK gives 4.4×10⁻¹⁶ (59×) [num].
  - The CK form is the one that agrees with GA08's v = 0 rate (GA08 k9, Mac Low & Shull 1986: 5.2×10⁻¹⁶ at 5000 K).
  - So the AAZN97 print is probably a typesetting error [interp].
- **Neither is PSS83.** Both are low-density (v = 0) rates. CK also zeroes the rate below 3481 K.
  - Consequence: there is **no H-impact H₂ dissociation in CK below 3481 K**.
  - CK's `_threebody_h2_equilibrium` (network_step.jl l. 46–57) then returns "all molecular" whenever k13 = 0.
  - RHFC02/PSS83 collapse has H₂ dissociating from ≈1550 K (RHFC02 §3.2), so this is a qualitative difference.

### A.8 H + H + H₂ → 2H₂ (PSS83 k6)

| source | formula | location |
|---|---|---|
| RA05 Table 1.1 [paper] | `k10 ≃ k8/8` (k8 = PSS83 k4) | RA05 p. 22 |
| RA04 [paper] | implicit: third-body density `n_f = n_p (1 − 15 f_H2/16)` multiplies **both** formation and dissociation. `n_f = n_H + n_H2/8`, i.e. H₂ is 1/8 as efficient a third body as H [interp, consistent with RA04 extraction §5.1] | RA04 p. 11, eqs. (29)–(30) |
| GA08 §2.1.5 [paper] | "To fix the rate of reaction 31, we follow Palla, Salpeter, & Stahler (1983) and assume that k31 = k30/8." | GA08 p. 4 |
| PSS83 [OCR-unverified] | the value is lost in the OCR | p. 634 |

RA04 eqs. (29)–(30), verbatim [paper, p. 11]:
```
F ≡ df_H2/dt = n_f [2 n k4 (1 − f_H2)² − k5 f_H2]                        (29)   (A9 writes 2 n_p k4)
n_f = n_p (1 − 15 f_H2/16)                                                (30)
n_p = ρX/m_H = n_H + 2n_H2 ;  f_H2 ≡ 2 n_H2/(n_H + 2n_H2)
f0(ρ,T) = 1 + (k5/(4 n_p k4)) [1 − (1 + 8 n_p k4/k5)^(1/2)]              (33)  (equilibrium, F = 0)
```

He is not counted as a third body in RA04. CK: absent.

### A.9 H₂ + H₂ → 2H + H₂ (PSS83 k7)

| source | formula | location |
|---|---|---|
| RA05 Table 1.1 [paper] | `k11 ≃ k9/8` (k9 = PSS83 k5) | RA05 p. 22 |
| RA04 [paper] | implicit via n_f (as A.8): H₂-impact dissociation = (RA04 k5)/8 | RA04 p. 11 |
| PSS83 [OCR-unverified] | the value is lost in the OCR | p. 634 |

CK: absent.

For comparison only (not in RHFC02): GA08 k10 [paper, p. 24]:
- v = 0 (ref. 9 = Martin, Keogh & Mandy 1998): `5.996×10⁻³⁰ T^4.1881 / (1.0 + 6.761×10⁻⁶ T)^5.6881 × exp(−54657.4/T)`.
- LTE (ref. 10 = Shapiro & Kang 1987): `1.3×10⁻⁹ exp(−53300/T)`.

### A.10 H + H → H⁺ + e⁻ + H (PSS83 k9)

| source | formula | location |
|---|---|---|
| RA05 Table 1.1 [paper] | `k13 ≃ 1.7 × 10⁻⁴ k12`, `k12 ≃ 5.8 × 10⁻¹¹ T^(1/2) e^(−158000/T) cm³ s⁻¹` | RA05 p. 22 |
| PSS83 [OCR-unverified] | k9 = 1.7×10⁻⁴ k8 (Drawin 1969); k8 = 5.8×10⁻¹¹ T^(1/2) exp(−158000/T) (Lotz 1967) | p. 634 |

**PSS83 ↔ RA05 numbering [interp].** PSS83 k4…k9 = RA05 k8…k13, an offset of +4. This is confirmed by PSS83 k5 ↔ RA05 k9 and PSS83 k8 ↔ RA05 k12 in the OCR.

**RHFC02 deliberately mixes sources.** It uses AANZ98 k1 for H + e ionization, not PSS83 k8. Yet the PSS83 k9 that it does use is defined as a multiple of PSS83 k8, not of AANZ98 k1. A faithful implementation must therefore use **1.7×10⁻⁴ × (5.8×10⁻¹¹ T^½ e^(−158000/T))**, not 1.7×10⁻⁴ × AAZN97 k1 [interp].

CK `k57` (Lenzuni, Chernoff & Salpeter 1991): `1.2e-17 T^1.2 exp(−1.578e5/T)`. Ratios are in A.1.

---

## B. Other H₂-related channels from AAZN97 (their appendix numbering)

All entries are [paper], AAZN97 pp. 25–26 (Tables 2–3). T is in eV unless "where T is in K".

| AAZN97 # | reaction | AAZN97 rate (verbatim) | AAZN97 source | other sources [paper] | CK |
|---|---|---|---|---|---|
| (9) | H + H⁺ → H₂⁺ + γ | `3.833×10⁻¹⁶ × T^1.8 cm³ s⁻¹` for T < 0.577 eV; `5.81×10⁻¹⁶ × (0.20651 T)^(−0.2891×log(0.20651×T)) cm³ s⁻¹` for T ≥ 0.577 eV | Shapiro & Kang 1987 | GP98 H8 = GA08 k3: `dex[−19.38 − 1.523 log T + 1.118 (log T)² − 0.1269 (log T)³]`, 1 ≤ T ≤ 32000 K (Ramaker & Peek 1976; Stancil et al. 1993). RA05 k5: `1.85×10⁻²³ T^1.8` (= AAZN97 low branch in K: 3.833e-16 × 11605^−1.8 = 1.85e-23 [num]) | `k9` = Latif et al. 2015, `dex[−18.20 − 3.194 lT + 1.786 lT² − 0.2072 lT³]` (clamped at 3.2e4 K; power law below 30 K). Within 8% of AAZN97 for 100–5000 K; −16% at 6700 K; −11% at 10⁴ K [num] |
| (10) | H₂⁺ + H → H₂ + H⁺ | `k10 = (6.4 ± 1.2)×10⁻¹⁰ cm³ s⁻¹` | Karpas et al. 1979 | GP98 H10 6.4×10⁻¹⁰; GA08 k4 6.4×10⁻¹⁰; RA05 k6 6.4×10⁻¹⁰ | `k10 = 6.0e-10`: **6% below every paper** (within the ±1.2 error) |
| (11) | H₂ + H⁺ → H₂⁺ + H | `ln(k11) = −24.24914687 + 3.40082444 ln(T) − 3.89800396 ln(T)² + 2.04558782 ln(T)³ − 0.541618285 ln(T)⁴ + 8.41077503×10⁻² ln(T)⁵ − 7.87902615×10⁻³ ln(T)⁶ + 4.13839842×10⁻⁴ ln(T)⁷ − 9.36345888×10⁻⁶ ln(T)⁸ cm³ s⁻¹` | "This work" (Janev 1987 data) | Text p. 8 also quotes Donahue & Shull: `k11 = 6.4×10⁻¹⁰ exp(−2.65 eV / T_eV)`. GP98 H15: `3.0×10⁻¹⁰ exp(−21050/Tg)` (Tg ≤ 10⁴), `1.5×10⁻¹⁰ exp(−14000/Tg)` (Tg > 10⁴) (Holliday et al. 1971); GP98 p. 5: the Janev-based AAZN97 rate "shows a dramatic drop at T < 10⁴ K ... should be disregarded". GA08 k7 (Savin et al. 2004), Table A1 p. 24: `[−3.3232183e-7 + 3.3735382e-7 log T − 1.4491368e-7 (log T)² + 3.4172805e-8 (log T)³ − 4.7813720e-9 (log T)⁴ + 3.9731542e-10 (log T)⁵ − 1.8171411e-11 (log T)⁶ + 3.5311932e-13 (log T)⁷] × exp(−21237.15/T)`. GA08 eq. 11 alt.: `2.4×10⁻⁹ exp(−21200/T)` (Shapiro & Kang). GA08 eq. 12 quotes AAZN97 k11 | `k11` = Savin 2004 with **natural** ln T, zero for T ≤ 0.3 eV |
| (14) | H⁻ + e⁻ → H + 2e⁻ | `k14 = exp[−18.01849334 + 2.3608522 ln(T) − 0.28274430 ln(T)² + 1.62331664×10⁻² ln(T)³ − 3.36501203×10⁻² ln(T)⁴ + 1.17832978×10⁻² ln(T)⁵ − 1.65619470×10⁻³ ln(T)⁶ + 1.06827520×10⁻⁴ ln(T)⁷ − 2.63128581×10⁻⁶ ln(T)⁸] cm³ s⁻¹` | Janev et al. 1987 (7.1.1); text: "net energy loss ... 0.755eV" | GA08 k14 identical (ref. 12) | `k14` identical; zero for T ≤ 0.04 eV (464 K) |
| (16) | H⁻ + H⁺ → 2H | `k16 = 7×10⁻⁸ (T/100 K)^(−1/2) cm³ s⁻¹, where T is in K` | Dalgarno & Lepp 1987 | GP98 H7: `5.7×10⁻⁶ Tg^−0.5 + 6.3×10⁻⁸ − 9.2×10⁻¹¹ Tg^0.5 + 4.4×10⁻¹³ Tg` (fit by PAMS; MAP). GA08 default k5 eq. (8): `2.4×10⁻⁶ T^−0.5 (1.0 + T/20000)` (Croft et al. 1999); eq. (9) = GP98 form (Moseley et al. 1970); eq. (10) `7.0×10⁻⁷ T^−0.5` (Dalgarno & Lepp; = AAZN97 k16); GA08 notes the last "may be erroneously small, owing to typographical errors in Dalgarno & Lepp (1987)" (Glover, Savin & Jappsen 2006) | `k16` = Croft (GA08 default): **3.4× AAZN97 at 100 K, 3.6× at 10³ K, 5.1× at 10⁴ K** [num] |
| (17) | H⁻ + H⁺ → H₂⁺ + e⁻ | `2.291×10⁻¹⁰ T^−0.4 cm³ s⁻¹` for T < 1.719 eV; `8.4258×10⁻¹⁰ T^−1.4 exp(−1.301/T) cm³ s⁻¹` otherwise | modified Shapiro & Kang | GP98 H6 = GA08 k16: `6.9×10⁻⁹ T^−0.35` (T ≤ 8000 K), `9.6×10⁻⁷ T^−0.90` (T > 8000 K) (Poulaert et al. 1978) | `k17`: `1.0e-8 T^−0.4` (T ≤ 10⁴ K), `4.0e-4 T^−1.4 e^(−15100/T)` (T > 10⁴ K). These are the AAZN97 branches converted to K (9.68e-9, 4.13e-4, 15098 K [interp]) but rounded, **with the break at 10⁴ K instead of 1.719 eV = 19950 K** |
| (18) | H₂⁺ + e⁻ → 2H | `1.0×10⁻⁸` (T < 617 K); `1.32×10⁻⁶ T^−0.76` (T > 617 K), T in K | fit to Schneider et al. 1994 | GA08 k6 identical; GP98 H11 `2.0×10⁻⁷ Tg^−0.5` | `k18` identical |
| (19) | H₂⁺ + H⁻ → H₂ + H | `k19 = 5×10⁻⁷ (100 K/T)^(1/2)`, T in K | Dalgarno & Lepp 1987 | GA08 k21 `1.4×10⁻⁷ (T/300)^−0.5` | `k19` identical to AAZN97 |
| (12) | H₂ + e⁻ → 2H + e⁻ | `k12 = 5.6×10⁻¹¹ T^(1/2) exp(−102,124 K/T)`, T in K | Donahue & Shull 1991 | GP98 H17 `4.4×10⁻¹⁰ Tg^0.35 exp(−102000/Tg)`; GA08 k8 v=0 `4.49×10⁻⁹ T^0.11 exp(−101858/T)`, LTE `1.91×10⁻⁹ T^0.136 exp(−53407.1/T)` (Trevisan & Tennyson 2002a) | `k12 = 4.4886e-9 T^0.109127 e^(−101858/T)` (TT02 v=0). 3.8× AAZN97 at 3000 K; 2.2× at 10⁴ K [num] |
| (15) | H⁻ + H → 2H + e⁻ | T > 0.1 eV: `exp[−20.37260896 + 1.13944933 ln T − 0.14210135 ln² T + 8.4644554e-3 ln³ T − 1.4327641e-3 ln⁴ T + 2.0122503e-4 ln⁵ T + 8.6639632e-5 ln⁶ T − 2.5850097e-5 ln⁷ T + 2.4555012e-6 ln⁸ T − 8.0683825e-8 ln⁹ T]`; T < 0.1 eV: `2.5634×10⁻⁹ × T^1.78186` | modified Janev | GA08 k15 identical | `k15` identical, except the low-T prefactor 2.56e-9 vs 2.5634e-9 (0.13%) |

### B.1 Case B recombination

**AAZN97 gives no case-B H rate.** Its k2 is a case-A fit (sum over n = 1..1000). RHFC02 likewise used case A.

The case-B fits in the local sources are [paper, GA08 p. 24, Table A1 k13, ref. 13 = Ferland et al. 1992]:
```
k13,A = 1.269×10⁻¹³ (315614/T)^1.503 [1.0 + (604625/T)^0.470]^−1.923     Case A
k13,B = 2.753×10⁻¹⁴ (315614/T)^1.500 [1.0 + (115188/T)^0.407]^−2.242     Case B
```

CK `k2` (labelled "Hui & Gnedin (1997) CaseB"): `4.881357e-6 · T^−1.5 · (1 + 114.813 · T^−0.407)^−2.242`.
- This is **algebraically identical to GA08 k13,B**: 2.753e-14 × 315614^1.5 = 4.8813e-6, and 115188^0.407 = 114.8 [num].
- So CK implements case B. A RHFC02-faithful run needs case A (AAZN97 k2 or GA08 k13,A).

### B.2 Physical threshold of H₂ + H⁺ ⇌ H₂⁺ + H [interp]

- **Energetics.** The reaction is endothermic by D(H₂) − D(H₂⁺) = 4.48 − 2.65 = 1.83 eV. AAZN97 p. 17 states "1.83 eV for process 10" (the reverse, exothermic direction).
- **Consistent thresholds.** The Savin fit's exp(−21237.15/T) is 1.830 eV and GP98's exp(−21050/T) is 1.81 eV. Both are consistent.
- **Inconsistent quote.** The Donahue & Shull form quoted by AAZN97, exp(−2.65 eV/T), uses the H₂⁺ binding energy and is not the reaction threshold.

### B.3 Notation and numerical checks [num]

- **GA08 k7 must use ln T, not log10 T.** GA08 Table A1 k7 prints "log T". Evaluating it with log10 gives a **negative** rate (−3.0×10⁻¹⁸ at 1000 K). With ln T it gives 1.30×10⁻¹⁹ at 1000 K and 2.91×10⁻¹¹ at 10⁴ K, consistent with GP98 H15 (2.2×10⁻¹⁹ and 3.7×10⁻¹¹). So GA08's "log" in k7 means ln. CK uses ln (correct).
- **AAZN97 k9 high branch also uses ln.** Its "log(0.20651 T)" is continuous with the low branch at 0.577 eV only with ln: 1.57×10⁻¹⁶ vs 1.43×10⁻¹⁶. With log10 the high branch gives 3.29×10⁻¹⁶.
- **AAZN97 k11 below 10⁴ K.** The Janev-based fit is 1.0×10⁻⁵⁰ at 1000 K, 1.6×10⁻¹⁹ at 3000 K, 2.3×10⁻¹⁴ at 5000 K and 1.6×10⁻¹¹ at 10⁴ K. Savin is 1.3×10⁻¹⁹, 2.2×10⁻¹³, 3.8×10⁻¹² and 2.9×10⁻¹¹ at the same temperatures. These match the "dramatic drop" GP98 describes.

---

## C. Three-body H₂ formation rates

### C.1 Every three-body formation rate in GA08 [paper, GA08 p. 4, §2.1.5; p. 19, §3.3.3; p. 26, Table A1]

| GA08 # | reaction | formula (verbatim) | GA08 attribution | note |
|---|---|---|---|---|
| 30 (default) | H+H+H → H₂+H | `k30 = 1.14×10⁻³¹ T^−0.38 cm³ s⁻¹` (T ≤ 300 K); `3.9×10⁻³⁰ T^−1.0 cm³ s⁻¹` (T > 300 K) — eq. (21) | "taken from Abel, Bryan, & Norman (2002) and partially based on Orel (1987) is the lowest of the values we have found in the literature" | **printed unit cm³ s⁻¹ is a typo for cm⁶ s⁻¹** |
| 30 (alt.) | H+H+H → H₂+H | `k30 = 1.44×10⁻²⁶ T^−1.54 cm³ s⁻¹` — eq. (22) | "from the recent paper of Flower & Harris (2007) has the highest value at low temperatures of any of the rates we have found" | same unit typo |
| 31 | H+H+H₂ → 2H₂ | `k31 = k30/8` | "we follow Palla, Salpeter, & Stahler (1983)" | — |
| 32 | H+H+He → H₂+He | `k32 = 6.9×10⁻³² T^−0.4` (Table A1 header says cm³ s⁻¹; it is a three-body rate, so cm⁶ s⁻¹) | ref. 24 = Walkauskas & Kaufman (1975) | not discussed in the text |

Further notes:
- GA08 §2.1.5 says the rates are "highly uncertain, as previously discussed in Glover (2007)". Glover 2007 is not available locally. The GA08 text quotes **no other numerical k30**. The PSS83 value is never quoted in GA08; GA08 credits PSS83 only for the /8 scaling.
- §3.3.3 (p. 19): with FH07 versus ABN02, final temperatures at n = 10¹² cm⁻³ differ by about 65%. FH07 heats faster for 5×10⁶ < n < 10⁸ (greater three-body formation heating) and is cooler above 10⁸.
- §3.3.4 "maximal" model uses ABN02; "minimal" model uses FH07.

Comparison with PSS83 and CK [num], in cm⁶ s⁻¹:

| T [K] | PSS83 k4 = 5.5e-29/T | ABN02 (GA08 eq. 21 = CK `k22`) | FH07 (GA08 eq. 22) |
|---|---|---|---|
| 300 | 1.83e-31 | 1.30e-32 | 2.21e-30 |
| 1000 | 5.50e-32 | 3.90e-33 | 3.45e-31 |
| 2000 | 2.75e-32 | 1.95e-33 | 1.19e-31 |
| 5000 | 1.10e-32 | 7.80e-34 | 2.90e-32 |
| 8000 | 6.88e-33 | 4.87e-34 | 1.40e-32 |

- CK `k22` = 1.3e-32 (T/300)^−0.38 / (T/300)^−1 is algebraically GA08 eq. 21 (ABN02): 1.14e-31 × 300^−0.38 = 1.304e-32 and 3.9e-30/300 = 1.3e-32.
- PSS83 / ABN02 = 14.1 for T > 300 K.
- FH07 / PSS83 = 6.3 at 1000 K, 4.3 at 2000 K, 2.0 at 8000 K.

### C.2 Collisional dissociation in GA08 and its density dependence [paper, GA08 p. 4, eqs. 13–17; p. 24, Table A1]

```
log k_i = (n/n_cr)/(1 + n/n_cr) · log k_i,LTE + 1/(1 + n/n_cr) · log k_i,v=0           (13)
1/n_cr  = x_H/n_cr,H + x_H2/n_cr,H2 + x_He/n_cr,He                                    (14)
x_H = n_H/n, x_H2 = 2n_H2/n, x_He = n_He/n, n = number density of hydrogen nuclei
n_cr,H  = dex[3.0 − 0.416 log T4 − 0.327 (log T4)²]                                   (15)
n_cr,H2 = dex[4.845 − 1.3 log T4 + 1.62 (log T4)²]                                    (16)
n_cr,He = dex[5.0792 {1.0 − 1.23×10⁻⁵ (T − 2000)}]                                    (17)
T4 = T/10000 K
```

Sources of the critical densities [paper]:
- n_cr,H is from Lepp & Shull (1983) "but has been decreased by an order of magnitude, as recommended by Martin, Schwarz & Mandy (1996)".
- n_cr,H2 is from Shapiro & Kang (1987).
- n_cr,He is from Dove et al. (1987).

H₂ + H → 3H dissociation (GA08 Table A1 k9):
- v = 0: `6.67×10⁻¹² T^0.5 exp[−(1 + 63593/T)]` (ref. 7, Mac Low & Shull 1986).
- LTE: `3.52×10⁻⁹ exp(−43900/T)` (ref. 8, Lepp & Shull 1983).

He-impact dissociation (k11, Dove et al. 1987): v = 0 `dex[−27.029 + 3.801 log T − 29487/T]`; LTE `dex[−2.729 − 1.75 log T − 23474/T]`.

### C.3 Detailed balance between formation and dissociation

**What the sources say**
- **GA08**: does **not** discuss detailed balance between k30 and k9. The words do not occur in GA08. It pairs the ABN02 or FH07 formation rate with the v = 0/LTE dissociation interpolation of eq. 13. These were never constructed as a detailed-balance pair.
- **RA04**: no explicit statement either. Its "equilibrium" f0 (eq. 33) is the kinetic steady state F = 0 of PSS83 k4 with modified k5.
- **PSS83** [OCR-unverified, p. 634]: dissociation rates were "derived, using the principle of detailed balance, from the measured recombination rates of Jacobs, Giedt, and Cohen (1967)".

**My check [num].** The equilibrium constant is
```
K(T) = n_H²/n_H2 = 16/(Z_rot Z_vib) · (π m_H k T/h²)^{3/2} · exp(−D0/kT)
```
with:
- D0 = 4.478 eV;
- rigid rotor θ_rot = 85.3 K, with ortho/para weights 3/1;
- harmonic Z_vib, θ_vib = 6332 K.

The table gives (k_diss/k_form)/K. A value of 1 means exact detailed balance.

| T [K] | K(T) [cm⁻³] | PSS83 k5/k4 (RA05 form) | RA04 k5/k4 | LS83-LTE k9 / PSS83 k4 | LS83-LTE / ABN02 | LS83-LTE / FH07 |
|---|---|---|---|---|---|---|
| 800 | 7.7e-5 | 0.26 | 0.12 | 979 | 1.4e4 | 138 |
| 1000 | 3.8e1 | 0.26 | 0.13 | 145 | 2042 | 23 |
| 1500 | 1.5e9 | 0.26 | 0.17 | 12 | 171 | 2.4 |
| 2000 | 1.0e13 | 0.25 | 0.19 | 3.7 | 53 | 0.86 |
| 3000 | 6.6e16 | 0.25 | 0.25 | 1.3 | 18 | 0.37 |
| 5000 | 7.1e19 | 0.25 | 0.34 | 0.69 | 9.8 | 0.26 |
| 8000 | 3.4e21 | 0.25 | 0.47 | 0.63 | 8.9 | 0.31 |

**Reading of the table [interp]**
- **PSS83 is a detailed-balance pair.** PSS83 k5/k4 = 0.25 K(T) with a constant factor over 800–8000 K. That confirms PSS83 derived k5 from k4 by detailed balance. The constant factor 1/4 is a statistical-weight convention, possibly the H electron-spin degeneracy.
- **RA04's modification breaks the T-scaling.** The modified k5 is no longer a detailed-balance partner of k4 (0.12→0.47 across the range).
- **CK's pair is not in balance.** CK's k13/k22 (DM86 v = 0 vs ABN02) is 0 below 3481 K. Above that it is 0.004 K at 4000 K, 0.013 K at 6000 K and 0.027 K at 8000 K.

---

## D. H₂ line cooling

### D.1 GP98 low-density limit [paper, GP98 p. 17, eq. A.7]

```
log Λ_H2[n(H) → 0] = −103.0 + 97.59 log Tg − 48.05 (log Tg)² + 10.80 (log Tg)³ − 0.9032 (log Tg)⁴    (A.7)
over the range 10 K ≤ Tg ≤ 10⁴ K.  Units: erg cm³ s⁻¹ (log = log10).
```

- **What it is.** H–H₂ collisions only. Ortho:para fixed at 3:1. Mandy & Martin (1993) rates above 600 K, Forrey et al. (1997) below.
- **Factor of 2.** This is a rate per H₂ molecule per atomic-H collider: the volumetric rate is `n(H) n(H₂) Λ`.
  - There is no ambiguity about molecule vs nucleus.
  - n(H) is the **atomic** H density (GP98 eq. A.6 "n(H)"; RA05 writes n_H).

### D.2 Hollenbach & McKee (1979) LTE fit

The text is identical in RA05 Eq. (1.23), p. 18 and RA04 Eq. (7), p. 4 [paper; both PDF-verified in the existing extractions]:

```
Λ_H2,LTE(T) = [9.5×10⁻²² T3^3.76 / (1 + 0.12 T3^2.1)] · exp[−(0.13/T3)³]
              + 3×10⁻²⁴ exp(−0.51/T3) + 6.7×10⁻¹⁹ exp(−5.86/T3) + 1.6×10⁻¹⁸ exp(−11.7/T3)     erg s⁻¹
T3 ≡ T/1000 K. "the first row ... accounts for rotational cooling, while the second row accounts for the
first two vibrational terms" (RA05).
```

**Factor of 2 [paper + interp]**
- Λ_LTE is per **H₂ molecule**. Check: H(T3 = 1) = 2.77×10⁻²¹ erg s⁻¹ [num], matching the GP98 per-molecule curve in RA04 Fig. 1.
- **RA04 eq. (6)** writes `L_lines,HM = (X f_H2/m_H) · H(T3)` per unit mass. X f_H2/m_H is the number of **H nuclei** in H₂ per gram, i.e. 2× the number of molecules.
- RA04 eq. (8) uses X f_H2/(2 m_H).
- So eq. 6, and through it eqs. 36 and A1, double-counts by a factor of 2 (RA04 extraction FLAG A).
- GP98 calls Λ(LTE) "the LTE cooling function given by Hollenbach & McKee (1979)" without restating it.

### D.3 Combination formula

**GP98** [paper, p. 17, eqs. A.5–A.6]:
```
Λ_H2[n(H), Tg] = Λ(LTE) / (1 + [n^cr/n(H)])                    (A.5)
n^cr/n(H) = Λ_H2(LTE) / Λ_H2[n(H) → 0]                          (A.6)
```

- **Unit inconsistency.** GP98 says Λ_H2 is "in erg cm³ s⁻¹" in the A.5 sentence. But Λ(LTE) (per molecule, erg s⁻¹) and Λ[n(H)→0] (erg cm³ s⁻¹) have different units, so A.6 as printed is dimensionally inconsistent.
- **Consistent reading [interp].** n^cr/n(H) = Λ_LTE / (n(H) Λ_{n→0}), and the result is per molecule in erg s⁻¹. This is exactly how RA05 writes it.

**RA05** [paper, p. 17, Eq. (1.22), PDF-verified]:
```
Λ_H2(ρ,T) = Λ_H2,LTE(T) / [1 + Λ_H2,LTE(T) / (n_H Λ_H2,ρ→0(T))]
"high and low density limits ... (which apply at n ≳ 10⁴ cm⁻³ and at n ≲ 10² cm⁻³, respectively)"
```

RA05 Eq. (1.24) [paper, p. 18] is GP98 A.7 with "10 K ≤ T ≤ 10⁴ K"; T in K and Λ in erg s⁻¹ cm³.

**GA08** [paper, p. 8, eqs. 37–39]:
```
Λ_H2,n→0 = Σ_k Λ_H2,k n_k ,  k = H, H2, He, H⁺, e⁻          (37)   per H2 molecule, erg s⁻¹
Λ_H2,LTE = Σ_{i,j>i} A_ji E_ji f_j                             (38)   per H2 molecule
Λ_H2 = Λ_H2,LTE / (1 + Λ_H2,LTE/Λ_H2,n→0)                      (39)
Footnote 1 (p. 5): "This cooling rate has units of erg cm³ s⁻¹, as do all of the other cooling rates quoted in this paper, unless indicated otherwise."
```

GA08 gives **no fit** for Λ_LTE; eq. 38 is a level sum. GA08 eq. (26) quotes GP98 as `Λ_H2,GP = dex[−103.0 + 97.59 log T − 48.05 (log T)² + 10.80 (log T)³ − 0.9032 (log T)⁴]`.

### D.4 GA08 rate coefficients (low-density limit, per H₂ molecule per collider, erg cm³ s⁻¹; T3 = T/1000 K; log = log10) [paper]

**Separate ortho/para fits**
- **H₂–H** (Wrathmall & Flower 2007), eq. 27, `log Λ = Σ_{i=0}^7 a_i (log T3)^i`. Coefficients in Tables 1 (ortho) and 2 (para), p. 6; reproduced below.
  - Ortho fit to within 2% for 100 < T < 6000 K. Para fit to within 3% over the same range.
  - Extrapolations below 100 K: ortho, eq. 28: `Λ_oH2,H = 5.09×10⁻²⁷ T3^(1/2) exp(−852.5/T)`; para, eq. 29: `Λ_pH2,H = 8.16×10⁻²⁶ T3^(1/2) exp(−509.85/T)`.
  - Combination, eq. 30: `Λ_H2,H = x_o/(x_o+x_p) Λ_oH2,H + x_p/(x_o+x_p) Λ_pH2,H`.
- **H₂–H₂**, eqs. 31–32, Tables 3 and 4, 100 < T < 6000 K, within 2%:
  `Λ_H2,H2 = x_p² Λ_pH2,pH2 + x_p x_o Λ_pH2,oH2 + x_o x_p Λ_oH2,pH2 + x_o² Λ_oH2,oH2`.
- **H₂–He**, eq. 33, Table 5. T < 6000 K, within 1%.
- **H₂–H⁺**, eq. 34, Table 6. 10 < T < 10⁴ K, within 2% (Gerlich 1990; Krstić 2002).
  - Ortho-para transitions J = 0 ↔ 1 by protons, eq. 35: `Λ_H2,H⁺,o⇔p = 4.76×10⁻²⁴ [9 exp(−170.5/T) x_p − x_o]`.
- **H₂–e⁻**, eq. 36: `log Λ = log[exp(−x/kT)] × Σ a_i (log T3)^i`, with x = 509.85k (para) and 845k (ortho). Table 7; within 10% for 10 < T < 10⁴ K.

**Table 8 (p. 9), 3:1 ortho:para, form `log Λ_H2 = Σ_{i=0}^5 a_i (log T3)^i` (eq. 40)**

| collider | T range [K] | a0 | a1 | a2 | a3 | a4 | a5 |
|---|---|---|---|---|---|---|---|
| H | 10 < T ≤ 100 | −16.818342 | 37.383713 | 58.145166 | 48.656103 | 20.159831 | 3.8479610 |
| H | 100 < T ≤ 1000 | −24.311209 | 3.5692468 | −11.332860 | −27.850082 | −21.328264 | −4.2519023 |
| H | 1000 < T ≤ 6000 | −24.311209 | 4.6450521 | −3.7209846 | 5.9369081 | −5.5108047 | 1.5538288 |
| H₂ | 100 < T ≤ 6000 | −23.962112 | 2.09433740 | −0.77151436 | 0.43693353 | −0.14913216 | −0.033638326 |
| He | 10 < T ≤ 6000 | −23.689237 | 2.1892372 | −0.81520438 | 0.29036281 | −0.16596184 | 0.19191375 |
| H⁺ | 10 < T ≤ 10000 | −21.716699 | 1.3865783 | −0.37915285 | 0.11453688 | −0.23214154 | 0.058538864 |
| e⁻ | 10 < T ≤ 200 | −34.286155 | −48.537163 | −77.121176 | −51.352459 | −15.169160 | −0.98120322 |
| e⁻ | 200 < T ≤ 10000 | −22.190316 | 1.5728955 | −0.21335100 | 0.96149759 | −0.91023195 | 0.13749749 |

"Note that the H₂-H⁺ rate quoted here does not include the effects of collisional transitions from J = 0 to J = 1 or vice versa" (add eq. 35 with x_p = 0.25, x_o = 0.75).

**Tables 1–7, verbatim [paper, pp. 6–8]**
- **T1, ortho-H₂–H**:
  - 100–1000 K: −24.330855, 4.4404496, −4.0460989, −1.1390725, 9.8094223, 8.6273872, 0.0, 0.0.
  - 1000–6000 K: −24.329086, 4.6105087, −3.9505350, 12.363818, −32.403165, 48.853562, −38.542008, 12.066770.
- **T2, para-H₂–H**:
  - 100–1000 K: −24.216387, 3.3237480, −11.642384, −35.553366, −35.105689, −10.922078, 0, 0.
  - 1000–6000 K: −24.216387, 4.2046488, −1.3155285, −1.6552763, 4.1780102, −0.56949697, −3.3824407, 1.0904027.
- **T3, para-H₂ excited by H₂**:
  - by para: −23.889798, 1.8550774, −0.55593388, 0.28429361, −0.20581113, 0.13112378.
  - by ortho: −23.748534, 1.76676480, −0.58634325, 0.31074159, −0.17455629, 0.18530758.
- **T4, ortho-H₂ excited by H₂**:
  - by para: −24.126177, 2.3258217, −1.0082491, 0.54823768, −0.33679759, 0.20771406.
  - by ortho: −24.020047, 2.2687566, −1.0200304, 0.83561432, −0.40772247, 0.096025713.
- **T5, He**:
  - para: −23.489029, 1.8210825, −0.59110559, 0.42280623, −0.30171138, 0.12872839.
  - ortho: −23.7749, 2.40654, −1.23449, 0.739874, −0.258940, 0.120573.
- **T6, H⁺**:
  - para: −21.757160, 1.3998367, −0.37209530, 0.061554519, −0.37238286, 0.23314157.
  - ortho: −21.706641, 1.3901283, −0.34993699, 0.075402398, −0.23170723, 0.068938876.
- **T7, e⁻**:
  - para T ≤ 10³ K: −22.817869, 0.95653474, 0.79283462, 0.56811779, 0.27895033, 0.056049813.
  - para T > 10³ K: −22.817869, 0.66916141, 7.1191428, −11.176835, 7.0467275, −1.6471816.
  - ortho: −21.703215, 0.76059565, 0.50644890, 0.050371349, −0.10372467, −0.035709409.

**Factor of 2 (GA08).**
- Every Λ_k is per H₂ molecule per collider, so the volumetric rate is n_H2 · Λ_H2 (eq. 37 sums Λ_k n_k).
- The H₂–H₂ term in eq. 37 uses n_k = n_H2 (molecules), not 2n_H2 [interp from eq. 32's quadratic-in-x form].
- GA08's chemistry uses x_H2 = 2n_H2/n (nuclei) for n_cr (eq. 14). Do not mix the two conventions.

**GA08 vs GP98 [paper, p. 8–9].** "the Galli & Palla cooling rate provides significantly more cooling at T < 1000 K than the newer Wrathmall & Flower cooling rate, with the rates differing most significantly at temperatures 300 ≲ T ≲ 500 K, where the Galli & Palla rate provides almost five times more cooling" [num: 6.0×/4.0× at 300/500 K; 1.2× at 1000 K].

### D.5 RA04's opacity factor (used by CK/EK as `ih2optical`) [paper, RA04 eq. 22, p. 9]

```
L_lines,thick(T) = L_lines,thin(T) · min(1, (n/n0)^−β),   n0 = 8×10⁹ cm⁻³ (ρ0 = 1.34×10⁻¹⁴ g cm⁻³), β = 0.45
n ≡ ρ/m_H (total mass density / m_H, not n_H). Eqs. 36/A1 print "max" — typo for min.
```

RHFC02 itself used full line transfer (J ≤ 20, v ≤ 2, LTE populations, Turner et al. 1977 A-values; p. 31), not any fit.

### D.6 Code vs papers (EK `cooling_h2.jl`, `emission.jl`) [num]

**Assembly** (`emission.jl` l. 39–57, 140–152):
```
Λ = n_H2 · fudge · [ LTE/(1 + LTE/Σ_k Λ_k n_k)  at T  −  same at T_CMB ]
```
This is the GA08 eq. 39 form, per molecule × molecule density. CK passes n_H2 = yH2I/2 and n_He = yHeI/4, so the normalisation is correct.

Additions not in GP98/GA08:
- **CMB subtraction.** Λ(T) − Λ(T_CMB), an Enzo convention.
- **Range clamps.** All fits are clamped to [10, 10⁴] K, extrapolating GA08's 6000 K limits for the H, H₂ and He colliders to 10⁴ K.

| code fn | label in code | matches paper? | numbers [num] |
|---|---|---|---|
| `GAH2` | GA08 | **= GA08 Table 8 H₂** | — |
| `GAHe` | GA08 | **= GA08 Table 8 He** | — |
| `GAHI` | "Lique 2015" | not in any local paper | vs GA08 Table 8 H: 1.0× at 100 K, 1.2× at 300 K, 1.6× at 500 K, 1.7× at 1000 K, 1.6× at 2000 K, 1.0× at 8000 K. Vs GP98: 0.43× at 100 K, 0.21× at 300 K, 1.4× at 1000 K, 1.9× at 2000 K. Returns 0 below 100 K (GA08 has a 10–100 K branch) |
| `GAHp` | "Glover & Abel 2008" | **does NOT match GA08 Table 8 H⁺**: code a0..a5 = −22.089523, 1.5714711, 0.015391166, −0.23619985, −0.51002221, 0.32168730 vs GA08 −21.716699, 1.3865783, −0.37915285, 0.11453688, −0.23214154, 0.058538864 | code/GA08 = 0.44 at 100 K, 0.42 at 1000 K, 0.64 at 5000 K, 0.68 at 8000 K |
| `GAel` | "Glover & Abel 2008" | **does NOT match GA08**: code uses two degree-8 branches split at 500 K and is zero below 100 K. GA08 uses degree 5, split at 200 K, valid from 10 K | code/GA08 = 0.68 at 100 K, 0.28 at 200 K, 0.18–0.20 at 300–1000 K, 0.42 at 3000 K, 0.68 at 8000 K |
| `H2LTE` | "Glover & Abel 2008" | **GA08 gives no LTE fit.** This is a different fit from HM79 (RA05 1.23 / RA04 7). Below 100 K: `7e-27 T^1.5 e^(−512/T)` | code/HM79 = 1.15 at 100 K, 0.92 at 200–300 K, 1.07 at 500 K, 0.94 at 1000 K, 0.99 at 2000 K, 0.92 at 3000 K, **0.72 at 5000 K, 0.55 at 8000 K** |
| `fudge` | (RA04) | uses **n_H** (H nuclei). RA04 uses n = ρ/m_H = n_H/X | Above 8×10⁹ the code's factor is X^−0.45 = 1.14× RA04 (less suppression); the onset shifts by 1/X = 1.33 |
| file header | "Galli-Palla (2008) H₂" (`EmissionKernels.jl` l. 17) | mis-citation: GP is 1998; GA is 2008 | — |

---

## E. Chemical heating and cooling energies

### E.1 Energies per reaction

| process | energy | source [paper] | note |
|---|---|---|---|
| H₂ binding / dissociation | **4.48 eV** per molecule | RHFC02 §3.2 p. 9: "4.48 eV per dissociated molecule are lost by the gas thermal energy". RA04 §3.1 p. 10: χ_H2 = 4.48 eV, eq. 25 `Ė_chem ≃ −ṅ_H2 χ_H2`. RA05 p. 39: "each H2 dissociation absorbs 4.48 eV of thermal energy" | RA04/RHFC02 apply it to formation (heating) too, with no suppression. RHFC02 does it through u_chem (eq. A8) |
| H⁻ + H → H₂ + e (AAZN97 8) | **3.53 eV** | AAZN97 p. 17 §6.1: "Most of the reaction enthalpy (1.83 eV for process 10 and 3.53 eV for process 8) goes rather into ro-vibrational excitation" | [interp] D(H₂) − EA(H) = 4.48 − 0.754 = 3.72 eV. 3.53 is AAZN97's number |
| H₂⁺ + H → H₂ + H⁺ (AAZN97 10) | **1.83 eV** | AAZN97 p. 17 | = D(H₂) − D(H₂⁺) |
| AAZN97 formation cooling | `Λ_H2 Formation = n_H (3.53 n_H⁻ k8 + 1.38 n_H2⁺) eV cm⁻³ s⁻¹` (PDF-verified) | AAZN97 p. 17 | **Typo in AAZN97**: "1.38" where the text says 1.83, and the k10 factor is missing. Intended: n_H(3.53 n_H⁻ k8 + 1.83 n_H2⁺ k10) [interp] |
| H⁻ + e → H + 2e (AAZN97 14) | 0.755 eV | AAZN97 p. 8 | threshold |
| H ionization (k1) | 13.6 eV | AAZN97 p. 6 | — |
| Solomon photodissociation | 0.4 eV per pair: `Γ27 = 6.4×10⁻¹³ n_H2 k27 erg cm⁻³ s⁻¹` (eq. 16) | AAZN97 p. 17 | not needed here |
| GA08 | **no energies given** | GA08 §2.3.9 p. 11: chemical cooling from reactions 12, 17, 18 (collisional ionization), 8–11 (collisional dissociation), 7 (charge transfer) and 13, 19, 20 (recombination); "chemical heating arising from the formation of H₂ via reactions 2 and 4"; "Further details ... can be found in Glover & Jappsen (2007)" | Glover & Jappsen 2007 is not available locally. GA08 §3.3.3 mentions "three-body H₂ formation heating" (reaction 30 heats) without a value |
| PSS83 | via ε per species | [OCR-unverified] p. 633 | as RHFC02 |

### E.2 Density-dependent suppression of formation heating

- **AAZN97** [paper, p. 17]: "These reactions, however, are important heat sources at high densities (n ≥ 10⁸ 1/cm³) where collisional de-excitation can transform most of the excitation energy into kinetic (thermal) energy. The density dependent heating rates are given by Hollenbach and McKee (1979, hereafter HM)."
  - This is qualitative only. **No formula and no n_cr is given**, and HM79 is not available locally.
- **GA08**: its n_cr (eqs. 14–17, §C.2) interpolates **dissociation rates** between the v = 0 and LTE limits. It is not a heating-efficiency factor. GA08 states no suppression formula for heating.
- **RHFC02 / RA04 / PSS83**: full thermalisation, no suppression. For RHFC02 this holds through u_chem; RHFC02 runs equilibrium chemistry above ρ_eq = 5×10⁻¹¹ g cm⁻³ anyway.
- **GP98 / RA05**: no chemical heating formula. RA05 extraction §2.6 says "NOT PRESENT".

### E.3 Code (CK `fast.jl` l. 82–100, `_h2_chem_heat`)

```
fform = k8·n_H⁻·n_HI (3.53 eV); gform = k10·n_H2⁺·n_HI (1.83 eV); tform = k22·n_HI³ (4.48 eV)
diss  = (k13·n_HI + k11·n_HII + k12·n_e)·n_H2 (4.48 eV absorbed)
fcr   = nH/(nH + 1e8)
edot  = eV·[fcr·(3.53 fform + 1.83 gform + 4.48 tform) − 4.48 diss]
```

Discrepancies with the sources:
1. **Unsourced suppression factor.** `fcr = n_H/(n_H + 10⁸ cm⁻³)` is attributed in the comment to "Palla, Salpeter & Stahler 1983 / Omukai 2000".
   - No local source gives this form or the constant n_cr = 10⁸.
   - AAZN97 only says "n ≥ 10⁸" qualitatively; PSS83 and RHFC02 thermalise fully.
   - A RHFC02-faithful run should use fcr = 1. RHFC02's network region (ρ < 5×10⁻¹¹ g cm⁻³) includes 10⁸–10¹³ cm⁻³ where fcr ≈ 1 anyway. Below 10⁸ it halves the H⁻ heating.
2. **Charge-exchange energy.** The k11 channel (H₂ + H⁺ → H₂⁺ + H) is charged 4.48 eV. Its endothermicity is 1.83 eV (B.2) [interp].
3. **Missing three-body channels in heating.** The three-body heating uses only H+H+H (`k22`, the ABN02 rate, 14× below PSS83). The H+H+H₂ formation channel (PSS83 k6) and the H₂+H₂ dissociation channel (PSS83 k7) are absent from both the network and the heating.
4. **Missing helium term.** Neither the network nor the heating has the GA08 k32 (H+H+He) term.

---

## F. Coverage

**Fully determined from primary sources (local PDFs of the originating paper)**
- RHFC02 #1 H recombination, case A (AAZN97 k2).
- RHFC02 #2 H⁻ formation (AAZN97 k7).
- RHFC02 #3 H₂ via H⁻ (AAZN97 k8).
- RHFC02 #4 H collisional ionization (AAZN97 k1).
- All of section B (AAZN97 k9–k19).
- Case-B H recombination (GA08 k13,B).
- GA08's own k30 values (as quoted; the primary papers ABN02 and FH07 were not checked).
- GP98 A.5–A.7.
- GA08 eqs. 26–40 and Tables 1–8.
- AAZN97 chemical energies.

**Determined only from secondary sources**
- **PSS83 k4** (RA04, RA05; unanimous).
- **PSS83 k5**: RA05 only for the original. RA04 gives a *modified* version.
- **PSS83 k6 = k4/8 and k7 = k5/8**: RA05; RA04 via n_f; GA08 for k6.
- **PSS83 k9 = 1.7×10⁻⁴ k8**: RA05 only.
- **HM79 LTE fit**: RA04, RA05. GP98 only names it.
- **PSS83 detailed-balance statement**: OCR only.
- The PSS83 OCR side-channel agrees with RA05 on k2, k4, k5, k8 and k9. It shows the RA05 k9 unit "cm⁶ s⁻¹" to be a RA05 typo.
