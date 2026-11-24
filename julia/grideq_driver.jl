# The GRIDEQ moving mesh: the settings (read from ENV when this file is evaluated by the run, popiii_run_body.jl, after envmenc/gm_tau) of the
# functions in PopIIIRun.jl (geq_*, grideq_step!, grideq_line: part of the TitanStack image). Uses the driver's globals a, b, nsh.
# GridEq.jl (the module) is loaded with the other module files.
"""
GRIDEQ=mode (with ENVEUL=0,… (regular centre), lgrid = 2; default off): one global, scale-free equidistribution rule for the whole
mesh (GridEq.jl): every zone between the first interface above the centre and the outer edge (both fixed in radius; with ICMLOG the first interface
is instead Lagrangian, its grid row holds its enclosed mass: TitanResidFast.GEQ_LAGFIRST, and it counts in the fmove limit) covers the same
∫ M d ln r, M = sqrt(1 + Σ_q w_q (∂ln q/∂ln r)²) over q = ρ, P, T, E_r, m and τ (clamped to [τlo, τhi]), with a smooth gradient cap
gmax and Dorfi–Drury smoothing α (neighbours within (1 + α)/α). No markers, thresholds or zone budgets.

  mode 1 (preferred): TITAN's implicit grid equation over the whole mesh (TitanResidFast.GRIDEQ), time-smoothed per zone with
          γ_j = κ t_cross,j / Δt, t_cross = r_c/(|u| + c_s) the zone's own signal-crossing time; the τ term (GEQW entry 6, wτ) is formed from the Newton iterate (τ = Σχρ Δr from the outer zone, ln τ smoothly confined to [ln τlo, ln τhi]) and differentiated; the speed floor of the relaxation time too;
          Δt limited so that no interface moves more than fmove of its smaller neighbouring zone per step.
  mode 0 (fallback): the explicit target of the same rule (GridEq.target); every interface relaxes to it in max(3 Δt, κ t_cross),
          at most fmove of the zone it moves into per step and |ṙ| ≤ fv (|u| + c_s); no zone shrinks below half its narrower
          neighbour. Mesh motion never limits Δt in this mode.

GEQN=N (with ICEUL): the IC is built on N zones (even; default: envlr's count). GEQW=wρ,wP,wT,wE,wm,wτ,wu,wχ (default 1,1,1,1,1,1,0,0; wτ: (Δ ln τ)² of the Newton iterate within τlo..τhi, mode 1; wχ: (∂ln(χρ)/∂ln r)² of the Newton iterate, mode 1 only; wu: the term wu ((r ∂u/∂r)/(|u| + c_s))², mode 1 only); GEQC=gmax,alph,κ,fv,fmove,τlo,τhi,ctol,κf,f_dyn (default 100,1.5,30,0.5,0.5,0.1,1e4,0.3,0,0.2; mode 1: κf is no longer used in τ; f_dyn: the grid relaxation time of every zone is smin(κ t_cross, f_dyn/√(Gρ)) of the Newton iterate, see TitanResidFast.grideq_rows!).

Feature tracking (κf > 0): the features of M move with the velocity v = −(∂M/∂t)(∂M/∂ln r)/((∂M/∂ln r)² + (2M)²) (in ln r per s;
the pattern speed of ∂M/∂t + v ∂M/∂ln r = 0, damped where M is flat: scale-free), estimated from M of this step and of the last
one; each zone then relaxes in max(min(κ t_cross, κf Δln r/|v|), Δln r t_cross/fv) (|v| the largest of the zone and its
neighbours; the floor keeps the grid speed below about fv (|u| + c_s)), so that the mesh keeps up with a moving photosphere or
shock without being fast everywhere.
With ICEUL the IC is built on the rule's mesh of the IC. No state beyond the mesh: a checkpoint restarts with the same setting.

GEQDEC=pr,pm[,a] (default 6,6,0.4): guaranteed resolution floor as a MIXTURE of separately normalised monitors. Physics monitor M_p = sqrt(1 + G²) (all terms above, G² capped at gmax²) and
floor monitor M_f = sqrt(ω_r² + ω_m² (∂ln m/∂ln r)²), ω_x = p_x/ln 10 (pr zones per decade of r, pm per decade of enclosed mass; one-sided stencil from the pinned first interface with ICMLOG).
The mesh equidistributes M = a (I_p/I_f) M_f + (1-a) M_p, I_x = Σ_j M_x,j Δξ_j over the zones, Δξ_j = 1/ν_j ≈ Δln r_j (implicit rows: x̂_j/M_j, TitanResidFast grideq_rows!, from the Newton iterate, differentiated;
I_p, I_f are global sums, dense dependence only through J·v; IC mesh: GridEq.density, I = trapezoid ∫M dln r). a = 0 gives the plain physics rows. Equidistribution gives M_j Δξ_j = I_p/N, hence
a (I_p/I_f) M_f,j Δξ_j <= I_p/N, i.e. Δln r_j <= I_f/(a N M_f,j) <= I_f/(a N ω_r) and, at fixed ∂ln m/∂ln r, Δln m_j <= I_f/(a N ω_m):
Δlog10 r <= I_f/(a N pr) and Δlog10 m <= I_f/(a N pm) (I_f >= ω_r Δln r_total: the floor is a guaranteed share a of the N zones, independent of how strong the physics features are; bounds hold up to
the (1+α)/α neighbour smoothing of the xnt row and the O(1) stencil approximations). I_f, I_p are printed in the grid-eq line.
The 1 in M_p and wm still act in M_p.
GEQX=wx[,xfl] (default 0: off): a composition term wx Σ_s (∂ ln(x_s + xfl)/∂ ln r)² over every species of the network (abundances
per H nucleus, ChemNet.full; xfl default 1e-3), frozen per step like the τ term (mode 1) and in the explicit target (mode 0). A
species far below xfl contributes nothing, so only changes of composition that carry a share ≳ xfl of the nuclei (ionization and
H2 dissociation fronts) are resolved; no list of species, no thresholds per front.
GEQRAMP=κ0,n[,f[,fvr]] (default off; fvr: the ramped speed floor never allows more than fvr (|u| + c_s), default unlimited): κ of the time smoothing (and with it the speed floor) rises geometrically from κ0 to GEQC's κ over
the first n steps after the start of the process, if the mesh is then far from the rule's, max_j(M_j Δln r_j) > f × median (default
f 10; a fixed-mesh checkpoint with its front in one zone: 26.5 for zdiagB's S5; an equidistributed mesh: ~3.4). A restart from a
fixed-mesh checkpoint needs it: κ 30 would take ~30 crossing times of the front to adapt the mesh, while Δt stays at the fixed
mesh's and the fixed-mesh run fails within ~30 steps. The next segment of the same run (mesh already adapted) does not ramp again.
"""
GRIDEQP = let v_ = get(ENV, "GRIDEQ", "off"); v_ in ("", "off") ? nothing : parse(Int, v_) end   # GRIDEQ=off: the fixed mesh (the stack default is 1)
GEQW = geqv("GEQW", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])   # entry 8 wχ: ln(χρ) of the Newton iterate (mode 1; with wτ = 0 and GEQX off the monitor has no lagged term)
GEQC = geqv("GEQC", [100.0, 1.5, 30.0, 0.5, 0.5, 0.1, 1e4, 0.3, 0.0, 0.2])
GEQX = geqv("GEQX", [0.0, 1e-3])
# GEQGP=w_g,L_ref [L_sun] (provisional, default 0 = off; implicit grid equation only): density jump x gravity power of the zone in the monitor
# (TitanResidFast.grideq_rows!, docs/POPIII_ENERGY.md §14.6)
GEQGP = geqv("GEQGP", [0.0, 100.0])
GEQDEC = geqv("GEQDEC", [6.0, 6.0, 0.4])   # zones per decade in r, in enclosed mass; a, the floor monitor's share of the zones (see above)
GEQFL = ((GEQDEC[1] / log(10))^2, (GEQDEC[2] / log(10))^2)   # (ω_r², ω_m²)
GEQRAMP = haskey(ENV, "GEQRAMP") ? parse.(Float64, split(ENV["GEQRAMP"], ",")) : nothing
