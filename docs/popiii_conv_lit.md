# TDC literature extraction (Kuhfuss 1986, Smolec & Moskalik 2008, Paxton 2019, Jermyn 2023)

All equations below were read from full text (arXiv PDFs via pdftotext; Kuhfuss 1986 from the ADS scan, read as page images;
MESA defaults from github.com/MESAHub/mesa `star/defaults/controls.defaults`, main branch, fetched 2026-11-22).
Tags: [FT] full text read; [derived] my own algebra from quoted equations; [garbled] PDF text extraction ambiguous.
NOTE: the Jermyn et al. 2023 paper (MESA VI) is arXiv:**2208.03651**, NOT 2208.07903 (that ID is an unrelated computer-vision paper).

## 0. Citations
- Kuhfuss R., 1986, A&A 160, 116-120, "A model for time-dependent turbulent convection". ADS 1986A&A...160..116K. [FT]
- Smolec R., Moskalik P., 2008, Acta Astron. 58, 193, "Convective hydrocodes for radial stellar pulsation. Physical and numerical formulation". arXiv:0809.1979. [FT]
- Paxton B. et al., 2019, ApJS 243, 10 (MESA V), arXiv:1903.01426, Sect. 2.1 (RSP), Tables 1, 3, 4. [FT]
- Jermyn A.S. et al., 2023, ApJS 265, 15 (MESA VI), arXiv:2208.03651 (v2, 30 Dec 2022), Sect. 3 "Time-Dependent Local Convection". [FT; the arXiv v2 text, journal version not compared]
- Wuchterl & Feuchtinger 1998 (A&A 340, 419; radiative-loss term) and Kuhfuss 1987 (MLT limit, alphas) are only cited second-hand [not verified].
- Farag et al. 2026, arXiv:2603.15766 (TDC_alpha_M eddy-viscosity implementation in MESA): seen only as a link in controls.defaults [not verified].

## 1. Kuhfuss 1986 (A&A 160, 116) [FT]
Symbols: varpi = <v'^2>/2 specific turbulent KE (RSP: e_t); Lambda = typical length (mixing length); mu = turbulent viscosity; <.> spherical average.
Closures (eq numbers as in paper):
- (19a) varpi = < v'^2/2 >. (19b) d_t varpi = (source) - (1/<rho>)<sigma_ij D_i v'_j> + (1/<rho>)<sigma~_ij D_i<v>_j> - (1/<rho>) div(j_t - <v' sigma~>).
- (19c) (source) := - < v' (1/rho grad p)' >.  (20) j_t := < v' (<rho> v'^2/2)' > (flux of turbulent energy).
- (22) j_w := <v'(<rho> w)'> (enthalpy flux).  (23b) turbulent Reynolds stress modelled as viscous stress with mu; trace part fixes p_t.
- (24) dissipation: (1/<rho>)<sigma_ij D_i v_j> ~ varpi / tau_Dis,  tau_Dis := Lambda / (varpi^{1/2} c_D).  (c_D dimensionless, order 1)
- (26) turbulent pressure p_t := (2/3)<rho> varpi  (remark: Stellingwerf overestimated it by factor 3/2).
- (30) velocity scale  dv := [<v'^2>/3]^{1/2} = sqrt(2/3) varpi^{1/2}.
- (31) mu = alpha_mu varpi^{1/2} Lambda.  (32) Re = (dv Lambda)/nu = sqrt(2/3)(1/alpha_mu)(mu/nu) >> 1.
- (34) j_i = -alpha_i Lambda varpi^{1/2} <rho> grad<c_i>; (35) alpha_i = sqrt(2/3) alpha_MLT (all alpha_i equal by sum j_i=0).
- (36) j_t = -alpha_t Lambda varpi^{1/2} <rho> grad varpi  (turbulent-energy flux = overshoot transport).
- (39)-(41) entropy flux <v'(<rho>s')> = -alpha_s Lambda varpi^{1/2}<rho> grad<s>; enthalpy flux
  (41b) j_w = -C varpi^{1/2} Lambda <rho> { <T> grad<s> + sum' <psi_i> grad<c_i> }, with (41a) C := alpha_s = alpha_i.
- (43) source = C varpi^{1/2} Lambda <alpha><T> (grad p)^2 /(<rho> p) * { nabla - nabla_ad + (composition term) }, alpha = -(1/rho)(drho/dT)_{p,c}
  (Schwarzschild-Ledoux criterion (45) follows: source>0 iff nabla > nabla_ad - composition term).
- Model equations (25a-e): 25d:
  d_t varpi = (source) - c_D varpi^{3/2}/Lambda - (1/<rho>) div j_t - (1/<rho>) p_t div<v> + (4/3) mu r^2 ( d_r (<v_r>/r) )^2 ;
  25c momentum: d_t<v> = g - grad(p+p_t)/<rho> + e_r (4/3)(1/(<rho> r^3)) d_r( r^4 mu <rho> d_r(<v_r>/r) ) ; 25e total energy.
- Four free parameters: C = sqrt(2/3) alpha_ML (convection drive), c_D (dissipation), alpha_t (overshoot, "order of C"), alpha_mu ("one or two orders smaller than C", citing Rodi 1979).
- Radiative losses of eddies are NOT in the original model (it has no D_r); gamma in (53) "V - V' = gamma(V' - V_ad)" is the MLT eddy-excess parameter, used only to relate to MLT cubic.
- Steady local limit (j_t = 0, no radial motion), eq (46): dv = [ 2 C A^2 (d_r p)^2 <alpha><T> (V - V_ad) / (3 c_D <rho> p) ]^{1/2}; with A = H_p, <alpha> ~ 1/<T> (ideal gas): (47) dv = [ 2 C p (V-V_ad)/(3 c_D <rho>) ]^{1/2}  (as printed: 2Cp / (3 c_D <rho>)).
  Matches Bohm-Vitense dv if (48) c_D = (32/9)(1/C); matches the MLT convective flux (incl. factor 1/2) if (49) c_D = (128/9)(1/C).
- Temperature-gradient equation (energy conservation, using (49)): (51) f^{3/2} = (8x/9)(1-f); MLT (Kippenhahn cubic) (52) f_MLT = 1 - (9/(8x))(sqrt(f_MLT) + x^2 - x)^3,
  with f := (nabla - nabla_ad)/(nabla_rad - nabla_ad) (50a), U (50b), x := U/(nabla_rad - nabla_ad) (50c). Ratio f_MLT/f never > ~1.8.
- Stability: eq (25d) damps varpi to 0 when V < V_ad (source negative) -- "source term can damp when layer becomes stable" (Kuhfuss's stated advantage over Stellingwerf 1982).
- Fig. 1 (p.120): plot of solutions f of cubic (51), f_MLT of cubic (52), and ratio f_MLT/f vs x in [-1, 11]. Fully reproducible (see Sect. 6).

## 2. Smolec & Moskalik 2008 (Acta Astron. 58, 193) [FT]
Equations (Sect. 2) as printed (R radius, U velocity, V=1/rho, E gas+rad energy):
- (1) dU/dt = -(1/rho) d(P+P_t)/dR + U_q - G M_R/R^2
- (2) dE/dt + P dV/dt = -(1/rho)(1/R^2) d[R^2(F_r+F_c)]/dR - C
- (3) de_t/dt + P_t dV/dt = -(1/rho)(1/R^2) d(R^2 F_t)/dR + E_q + C        (4) sum of (2),(3) = total energy.
- (6) P_t = alpha_p rho e_t.   (10) mu_Q = alpha_m rho Lambda e_t^{1/2}, Lambda = alpha H_P.
- (7) U_q = (1/(rho R^3)) d/dR( (4/3) mu_Q R^3 d(U/R)/dR ) = (4 pi/R) dX/dM ;  (8) E_q = (4/(3 rho)) mu_Q R^2 (d(U/R)/dR)^2 ; (9) X = (16/3) pi mu_Q R^6 rho d(U/R)/dM.
- (11) C = S - D - D_r.
- (12) D = alpha_d e_t^{3/2}/Lambda.
- (13) D_r = 4 sigma gamma_r^2 T^3 V^2 e_t /( alpha^2 c_P kappa H_P^2 ).
- (14) S = (T P Q /(c_P H_P)) Pi ,  Q = (dV/dT)_P.
- (15) Pi = alpha alpha_s e_t^{1/2} c_P Y.   (16) Y = -(H_P/c_P) ds/dR = nabla - nabla_a.
- (17) F_c = (alpha_c/alpha_s) rho T Pi = alpha alpha_c rho T c_P e_t^{1/2} Y.
- (18) F_t = -alpha_t rho Lambda e_t^{1/2} de_t/dR  ("sole cause of overshooting").
- (19) optional Wuchterl-Feuchtinger enthalpy-flux limiter FL replacing (15).
- (20) F_r = -(4 sigma/3) 4 pi R^2 (1/kappa) dT^4/dM (diffusion).
- Source S keeps sign of Y: negative (damping) in convectively stable zones (S ~ Y, Sect. 5.3). Variant "S ~ Y_+" (non-negative, "PP model") and Florida-Budapest code are discussed; "NN" = negative allowed. Sect. 5.3: PP leaves large e_t in stable zones -> too much eddy-viscous damping; NN drops e_t by ~25 decades below the convection zone (Fig. 9a).
- Numerics: e_t = 0 is an absorbing solution (all terms power of e_t); fix: D -> (eq. 25, footnote/eq.28) D = alpha_d (e_t^{3/2} - e_0^{3/2})/Lambda with a small floor e_0 (e0=1 erg/g in the code; results independent for e0 in 1e-4..1e6 erg/g; e_t reset to e0 if below). e_t=0 at inner and outermost boundary (e_{t,N}=0).
- Table 1 parameters (alpha values are in units of standard values): alpha (-), alpha_m (-), alpha_p 2/3, alpha_s (1/2)sqrt(2/3), alpha_d (8/3)sqrt(2/3), alpha_c (1/2)sqrt(2/3), alpha_t (-), gamma_r 2 sqrt3.
  Text: "Neglecting radiative losses, local static solution of eqs (1)-(3) may be reduced to MLT if alpha_c = alpha_s = (1/2) sqrt(2/3) and alpha_d = (8/3) sqrt(2/3) (Kuhfuss 1987). With radiative losses included, exact MLT cannot be reproduced. Following W&F 1998: gamma_r = 2 sqrt3. These are 'standard values' ... reference values only." alpha_p and alpha_c introduced by Yecko et al. 1998 (original: alpha_p = 2/3, alpha_c = alpha_s).
- Table 2 (test sets, units of standard values): A: alpha=1.5, alpha_m=0.20, alpha_s=alpha_c=alpha_d=1, alpha_p=0, alpha_t=0, gamma_r=0; B: alpha_m=0.25, alpha_p=1, alpha_t=0.01.
- Tests: Fig. 1 instability strips (sets A,B); Fig. 2 static 4.5 Msun model: Y, nabla_ad, F_c, F_t, e_t vs zone showing overshoot into Y<0; Fig. 3 linear work integrals; Figs 4-6 nonlinear limit cycles; Figs 9-10 NN vs PP e_t profiles/work integrals. All need the pulsation hydrocode/model builder.

## 3. Paxton et al. 2019 (ApJS 243, 10), Sect. 2.1 [FT]
- Eqs (1)-(4) are the Smolec & Moskalik equations in mass coordinate (A=4 pi r^2): (3) De_t/Dt + p_t DV/Dt = -d(A F_t)/dm + q + C ; C = S - D - D_r.
- Table 1 definitions: F_c = alpha alpha_c rho T c_p e_t^{1/2} Y_sag ; F_t = -(alpha alpha_t rho H_p e_t) d e_t/dr [garbled in extraction; Smolec eq.18 has e_t^{1/2} and Lambda=alpha H_p]; p_t = alpha_p rho e_t ; q(visc) = alpha_m rho alpha H_p e_t^{1/2};
  S = alpha alpha_s e_t^{1/2} (T p Q/H_p) Y_sag ; D = alpha_d e_t^{3/2}/(alpha H_p) [garbled; consistent with Smolec 12]; D_r = (4 sigma gamma_r^2/alpha^2)(T^3 V^2/(c_p kappa H_p^2)) e_t.
- Table 3 (parameter = base value x MESA control): alpha 1 x RSP_alfa; alpha_m 1 x RSP_alfam; alpha_s (1/2)sqrt(2/3) x RSP_alfas; alpha_c (1/2)sqrt(2/3) x RSP_alfac; alpha_d (8/3)sqrt(2/3) x RSP_alfad; alpha_p 2/3 x RSP_alfap; alpha_t 1 x RSP_alfat; gamma_r 2 sqrt3 x RSP_gammar.
- Text: "If radiative cooling and turbulent pressure are neglected, the time-independent version of Kuhfuss reduces to standard MLT provided base values are used for alpha_s, alpha_c, alpha_d (controls set to 1)." Suggested starting values: alpha_t ~ 0.01, alpha_m <~ 1, alpha <~ 2.
- Table 4 sets (RSP_alfa/alfam/alfas/alfac/alfad/alfap/alfat/gammar): A 1.5/0.25/1/1/1/0/0/0; B 1.5/0.50/1/1/1/0/0/1; C 1.5/0.40/1/1/1/1/0.01/0; D 1.5/0.70/1/1/1/1/0.01/1.
- Tests (RSP): Fig. 4 growth rate Gamma = 2(E_k,max^{i+1}-E_k,max^i)/(sum) vs LNA predictions (RR Lyrae, M=0.65, L=45, Teff=7100); Fig. 5 mode selection; Figs 6-10 parameter sensitivity of light curves. Need RSP.
- MESA RSP defaults (controls.defaults, main, 2026-11): RSP_alfa=1.5, RSP_alfac=1, RSP_alfas=1, RSP_alfad=1, RSP_alfap=0, RSP_alfat=0, RSP_alfam=0.25, RSP_gammar=0 (all in units of standard values above; comment: "RSP uses P_turb = RSP_alfap*(2/3)*rho*e_t in both momentum and turbulent energy eqs"). Other RSP: RSP_cq=4.0, RSP_zsh=0.1, RSP_nz=150, RSP_nz_outer=40, RSP_T_anchor=11d3.

## 4. Jermyn et al. 2023 (ApJS 265, 15), Sect. 3 [FT]
Variable: e_t = w^2 ("w = turbulent velocity"); steady link w = sqrt(3/2) v_c (v_c = MLT convective velocity). Y := nabla - nabla_ad (5) (Ledoux gradient nabla_L replaces nabla_ad with Ledoux criterion). Local limit: L_t = 0 in (6), and in the solver.
- (2) De_t/Dt + alpha_Pt P_t D rho^-1/Dt = q + C - dL_t/dm ,  P_t = (2/3) rho e_t. (q dropped from MESA energy equation.)
- (3) L_t = -A alpha alpha_t rho h e_t^{1/2} de_t/dr  (h = P/(rho g)); [turbulent transport -- NOT used: "for simplicity we currently set L_t = 0"].
- (4) C = alpha e_t^{1/2} (T P Q/(h sqrt6)) Y - alpha_D (8 sqrt2/(3 sqrt3)) e_t^{3/2}/(alpha h) - (48 sigma alpha_r/alpha^2)(T^3/(rho^2 c_P kappa h^2)) e_t.
  i.e. S with alpha_s=1/sqrt6; D with alpha_d = alpha_D (8/3)sqrt(2/3); D_r with gamma_r^2/12 -> alpha_r [derived: 4 gamma_r^2 = 48 alpha_r, so RSP standard gamma_r=2 sqrt3 <-> alpha_r = 1].
- (6) L = L_rad + L_conv ; (7) L_conv = 4 pi r^2 alpha rho c_P T (w/sqrt6) Y  [alpha_c = 1/sqrt6 folded in]; MLT-consistent version (24): L_conv = 4 pi r^2 alpha (Gamma/(1+Gamma)) rho c_P T (w/sqrt6) Y, Gamma from the Cox MLT call (held fixed during iteration); implemented when include_mlt_corr_to_TDC = .true. (default). In effect Y -> Y Gamma/(1+Gamma) in source/flux/velocity.
- Rewriting with e_t = w^2 (8): 2 Dw/Dt = xi0 + xi1 w + xi2 w^2 (divided by w; branch w=0 excluded), with (9)-(11)
  xi0 = (alpha/(h sqrt6)) c_P T nabla_ad Y ;
  xi1 = -[ (4 sigma T^3/(rho^2 c_P kappa)) (2 sqrt(3 alpha_r)... /(alpha h))^2 [garbled; from (4) it is 48 sigma alpha_r T^3/(alpha^2 rho^2 c_P kappa h^2)] + alpha_Pt (2/3) rho D rho^-1/Dt ] ;
  xi2 = -(8 sqrt2/(3 sqrt3)) alpha_D/(alpha h).
  Coefficients frozen at end-of-step values; w solved at end of step from start-of-step w (implicit, unconditionally stable).
  Discriminant J^2 := xi1^2 - 4 xi0 xi2.
- J^2 > 0 (convectively unstable) (12): w = -(1/(2 xi2)) [ J tanh((lambda + J dt)/4) + xi1 ], lambda fixed by initial w; dt->inf: w -> -(J+xi1)/(2 xi2) independent of w_0 (plateau).
- J^2 < 0 (stable) (13): w = (1/(2 xi2)) [ |J| tan((lambda + |J| dt)/4) - xi1 ]; reaches w=0 at finite dt* ~ 1/J; then stays w=0 (valid solution of (2)). Solver must check whether first root occurs inside the step and then set w=0 at end of step.
- Start-up: Y<0 and w_0=0 => no pre-existing turbulence, w=0 for whole step, L = L_rad, (16) linear in Y. Y>0 with w_0=0: NOTE the paper only states the w>0 solution is used when J^2>0; it does not describe a seed value; exact start from w_0=0 exactly is a fixed point of (2) (all terms ~ w^k) [not verified how MESA seeds; check star/private/tdc*.f90 source].
- Y is NOT advanced in time; the solver inverts L(Y) = L_solver ((14)-(16)): R(Y) = (L_solver - L0 nabla_L) - (L0 + c0 w) Y = 0, L0 = (16 pi a c/3) G m T^4/(kappa P), c0 = 4 pi r^2 (alpha/sqrt6) rho T c_P (17). Sign of Y fixed by sign of (L_solver - L0 nabla_L) (L0+c0 w >0). Uses Z = ln|Y| in [-100,100]; dw/dY >= 0.
  Y>0: monotonic -> bisection in Z (Delta Z<=1) then Newton-Raphson (autodiff partials).
  Y<0, w_0>0: up to 3 roots (Fig 1); rule: choose the root with smallest |Y| (slowest-decaying convection), found via special points Y0 (smallest |Y| with w=0) and Y1 (dR/dY=0, w>0).
- Complementarity: no explicit complementarity formulation; piecewise: w>=0 enforced by branch logic (tan solution clipped to w=0 once root passed; w=0 branch if w_0=0, Y<0).
- Spatial transport: none (L_t=0, local limit). Overshoot only via MESA's separate overshoot schemes; paper warns (and controls.defaults MLT_option comment repeats) that combining TDC with overshoot/other mixing models "might lead to physically inconsistent solutions". Eddy viscosity q dropped.
- Steady limit/MLT (Sect. 3.5): alpha_r=0, Dp/Dt=0 => xi1=0, vc^2 = (1/6)(J^2/xi2^2)... (31) vc^2 = (1/6) J^2/xi2^2 = -2 xi0/(3 xi2) = alpha^2 c_P T nabla_ad Y/(8 alpha_D); with Cox f1=1/8, f2=1/2 and alpha_D=1, identity (35) 1/(8 alpha_D f1) = (Gamma3-1) c_V/(chi_rho c_P nabla_ad)=1 gives equality with MLT (25)-(28) provided Y -> Y Gamma/(1+Gamma). Agreement to ~7 decimals in Y. The correction is needed even for alpha_r>0; without it vc scaling in the inefficient limit differs (Fig. 2).
- Parameters (paper): alpha = 2 (mixing_length_alpha default 2), alpha_D = 1, alpha_Pt = 0, alpha_r = 0 ("By default alpha_r = 0 ... approximation that changes physical contents", cites Kuhfuss 1986/87, W&F 1998); alpha_c, alpha_s fixed 1/sqrt6; alpha_t only in (3) (unused, L_t=0).
- Current MESA control names (controls.defaults main; old names in parentheses): TDC_alpha_D = 1 (alpha_TDC_DAMP), TDC_alpha_R = 0 (alpha_TDC_DAMPR), TDC_alpha_Pt = 0 (alpha_TDC_PtdVdt), TDC_alpha_C = 1, TDC_alpha_S = 1, TDC_alpha_M = 0 (eddy viscous damping, analogous to RSP_alfam; RSP default 0.25), MLT_option = 'TDC' (default in main), mlt_Pturb_factor = 0 (P_turb = factor*rho*vc^2/3; "=3 equals TDC_alpha_Pt=1"), include_mlt_corr_to_TDC = .true., TDC_include_eturb_in_energy_equation = .false., steps_before_use_TDC = 0, TDC_num_{innermost,outermost}_cells_forced_nonturbulent = 0, use_rsp_form_of_scale_height = .false. In main, TDC_alpha_C/S are scalings on top of the 1/sqrt6 constants (per controls.defaults text "prefactor on convective flux/source"); the 2022 paper does not contain them.

## 5. Side-by-side table
| Term | Kuhfuss 1986 | RSP (Smolec&Moskalik 2008 / Paxton 2019) | MESA TDC (Jermyn 2023 / controls.defaults) |
|---|---|---|---|
| Variable | varpi = <v'^2>/2 (19a) | e_t | e_t = w^2 |
| Source S | C varpi^{1/2} Lambda <alpha><T>(grad p)^2/(<rho>p)(nabla-nabla_ad) (43); C = sqrt(2/3) alpha_ML | alpha alpha_s e_t^{1/2} (T P Q/H_P) Y (14,15); alpha_s std (1/2)sqrt(2/3), control RSP_alfas=1 | alpha e_t^{1/2} (T P Q/(h sqrt6)) Y (4); TDC_alpha_S=1; Y->Y Gamma/(1+Gamma) |
| Dissipation D | c_D varpi^{3/2}/Lambda (24,25d); c_D=32/(9C) or 128/(9C) | alpha_d e_t^{3/2}/(alpha H_P) (12) (floor e_0: e_t^{3/2}-e_0^{3/2}); std alpha_d=(8/3)sqrt(2/3); RSP_alfad=1 | alpha_D (8 sqrt2/(3 sqrt3)) e_t^{3/2}/(alpha h); TDC_alpha_D=1 |
| Radiative damping D_r | none in original (added by W&F 1998) | 4 sigma gamma_r^2 T^3 V^2 e_t/(alpha^2 c_P kappa H_P^2) (13); std gamma_r = 2 sqrt3; RSP_gammar=0 (default) | 48 sigma alpha_r T^3 e_t/(alpha^2 rho^2 c_P kappa h^2); TDC_alpha_R=0 (alpha_r = gamma_r^2/12 [derived]) |
| Convective flux F_c | enthalpy flux j_w (41b), C | alpha alpha_c rho T c_P e_t^{1/2} Y (17); std alpha_c=(1/2)sqrt(2/3); RSP_alfac=1 | 4 pi r^2 alpha (Gamma/(1+Gamma)) rho c_P T w Y/sqrt6 (24); TDC_alpha_C=1 |
| Turbulent pressure P_t | (2/3)<rho> varpi (26), in momentum (25c) | alpha_p rho e_t (6); std alpha_p=2/3 (so P_t=(2/3)rho e_t); RSP_alfap=0 default | (2/3) rho e_t; in e_t eq. via alpha_Pt (TDC_alpha_Pt=0); momentum via mlt_Pturb_factor=0 |
| Turbulent energy flux / overshoot | j_t = -alpha_t Lambda varpi^{1/2} <rho> grad varpi (36) | F_t = -alpha_t rho Lambda e_t^{1/2} d e_t/dr (18); RSP_alfat=0 (0.01 suggested) | L_t (3) defined but set to 0; none (local) |
| Eddy viscosity | mu = alpha_mu varpi^{1/2} Lambda (31); terms in 25c,25d | mu_Q = alpha_m rho Lambda e_t^{1/2} (10), U_q, E_q (7,8); RSP_alfam=0.25 | q dropped (paper); TDC_alpha_M=0 default (later feature, Farag 2026 [not verified]) |
| P_t dV work | -(1/rho)p_t div<v> (25d) | + P_t dV/dt on LHS of (3) | alpha_Pt P_t D rho^-1/Dt in (2); alpha_Pt=0 |
| Lambda | Lambda (A) | alpha H_P, RSP_alfa=1.5 | alpha h, h=P/(rho g); alpha=mixing_length_alpha (default 2) |

## 6. Steady-state limits and MLT relation
- Kuhfuss: local, j_t=0, no radial motion: (46)-(47). MLT recovered for c_D=32/(9C) [BV velocity] or 128/(9C) [BV flux]; cubic (51) vs MLT cubic (52) differ by <=1.8x in f (Fig 1) because no "horizontal" eddy-cooling.
- RSP: alpha_c=alpha_s=(1/2)sqrt(2/3), alpha_d=(8/3)sqrt(2/3), gamma_r=0 (and alpha_p=alpha_t=0) -> standard MLT (Kuhfuss 1987, via S&M p.6). gamma_r=2 sqrt3 (W&F 1998 value) is part of the "standard set" but exact MLT is lost when D_r on. Brief check [derived]: S=D, with T P Q = g H_p delta, gives e_t = alpha^2 (alpha_s/alpha_d) g delta H_p Y, v_c^2=(2/3)e_t = alpha^2 (1/8) g delta H_p Y, i.e. MLT f1 = 1/8 for alpha_s/alpha_d = 3/16 (= standard ratio). Then F_c=alpha alpha_c rho T c_P e_t^{1/2} Y gives f2=1/2 with alpha_c=1/(2)sqrt(2/3).
- MESA: alpha_D=1 <=> alpha_d standard; constants (alpha_s=alpha_c=1/sqrt6=(1/2)sqrt(2/3) fixed); alpha_r=0; plus Gamma/(1+Gamma) factor => Cox MLT (f1=1/8, f2=1/2) in long-step limit (7 decimals). RSP itself (without the factor) is only MLT in the efficient limit; Kuhfuss-form differs in the inefficient limit.

## 7. Reproducible without MESA
1. Kuhfuss 1986 Fig. 1: solve cubics (51) f^{3/2}=(8x/9)(1-f) and (52) f_MLT = 1-(9/(8x))(sqrt f_MLT + x^2 - x)^3 for x in [-1,11]; plot f, f_MLT and f_MLT/f (peak ~1.8 at small x, ->1 at large x).
2. Jermyn Fig. 1 (toy): solve R(Y)=0 with R=(L_solver - L0 nabla_L) - (L0+c0 w)Y, w from (13) with L_solver=0.1, L0=2, nabla_ad=1 (L0 nabla_L term: nabla_L as printed [not verified], c0=20, xi1=0, xi2=-1, w_initial=1, |J|=-Y, dt=5; shows 3 roots (w~0.2,Y~-0.35), (w~0.05,Y~-0.75), (w=0,Y~-0.95).
3. Jermyn Fig. 2: Y and v_c vs (nabla_rad - nabla_ad) (1e-7..10), MLT (Cox) vs TDC at long dt vs TDC without Gamma/(1+Gamma). Inputs of the stellar zone are not stated in the text [not verified; Gamma needs the Cox cubic]. Algebraic check available: eqs (28)-(35).
4. Pointwise ODE tests from Jermyn eq. (2)/(8) with frozen coefficients: tanh growth to plateau -(J+xi1)/(2 xi2) (J^2>0); tan decay reaching w=0 at finite dt* (J^2<0). Analytic, no stellar model needed.
5. Smolec & Moskalik Fig. 9/10 (NN vs PP), RSP Figs 4-5 (growth rate vs LNA), Jermyn Figs 3-4 (He-flash WD) need the full codes/MESA: not reproducible standalone.
