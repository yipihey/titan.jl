# S7b: coupled protostellar feedback — the pieces that do not depend on a TITAN state. Loaded into Main after ChemNet, AtomicHHe, ContinuumOpacity, RadTransfer.
#   - exchange records (key = value text files, atomic rename) and waiting
#   - the source spectrum of the star+precursor (S6b record) and its attenuation through the analytic free-fall column between the stellar surface and the sink radius r_s
#     (the region that is not on the envelope grid): `inner_column`, a photon-conserving equilibrium sweep (RadTransfer.ion_sweep) of H, He, He⁺ in the prescribed flow
#     ρ ∝ r^{-3/2}, u ∝ r^{-1/2} normalised to the density and velocity of the envelope at r_s (continuity at the grid edge)
#   - the field on the envelope grid: `envfields`, ion_sweep on the hydro state (equilibrium ionisation of the lit gas, where the front crosses many zones per call), then RadTransfer.transfer on the
#     updated composition (other continua, diffuse field, LW, radiation acceleration)
#   - semi-analytic estimates (OI02 criteria, Lyα pressure bound)
module S7b
using Printf
import ..ChemNet, ..RadTransfer, ..AtomicHHe
const CN = ChemNet; const RT = RadTransfer
const G = 6.67430e-8; const Msun = 1.98847e33; const Rsun = 6.957e10; const AU = 1.495978707e13; const YR = 3.15576e7; const Lsun = 3.828e33
const σSB = 5.670374419e-5; const c = CN.c; const hP = CN.hP; const eV = CN.eV; const kB = CN.kB; const mp = CN.mp; const σT = 6.6524587e-25

# exchange records: s7b_exchange.jl (module S7bX, no dependencies; loaded first)
const write_rec = Main.S7bX.write_rec; const read_rec = Main.S7bX.read_rec; const wait_rec = Main.S7bX.wait_rec

# ------------------------------------------------------------------ source spectrum
"""
    star_source(g, rec; mode = :star) -> RT.Source

Spectrum of the radiation that leaves the stellar surface. `mode = :star`: blackbody of the stellar surface (T_eff,* from L*, R*) plus, for the shock luminosity that does not
go into the ionisation of the accreted gas (L_∞ − L*), a blackbody at the precursor photosphere temperature T_ph; `:photosphere`: only the blackbody (T_ph, L_∞) of the precursor photosphere (S6b).
`rec` keys: Lstar, R, Linf, Tph [cgs; L in erg/s]. No Z = 0 atmosphere corrections.
"""
function star_source(g::RT.FreqGrid, rec; mode = :star)
    Ls = max(rec["Lstar"], 1e-3 * Lsun); R = rec["R"]; Linf = max(rec["Linf"], Ls); Tph = max(rec["Tph"], 3000.0)
    if mode === :photosphere
        return RT.blackbody_source(g, Tph, Linf; Rstar = R)
    end
    Teff = (Ls / (4π * R^2 * σSB))^0.25
    s1 = RT.blackbody_source(g, Teff, Ls; Rstar = R)
    Lp = Linf - Ls
    Lp <= 1e-6 * Ls && return s1
    s2 = RT.blackbody_source(g, Tph, Lp)
    RT.Source(s1.L .+ s2.L, R, Ls + Lp)
end
"Ionising photon rate (≥ 13.6 eV) of a power-per-bin vector [s⁻¹]."
qion(ia::RT.IonAtoms, L::AbstractVector) = sum(L[f] / (hP * ia.g.ν[f]) for f in ia.f0:length(L); init = 0.0)

# ------------------------------------------------------------------ analytic free-fall column inside the sink radius
"""
    inner_column(ia, src, M, Rin, rs, rho_s, u_s; nz = 300, T = 1.5e4, passes = 6, reemit = nothing) -> NamedTuple

The photons of `src` (power per bin of the star+precursor, emitted at the radius `Rin` of the stellar surface) cross the gas between `Rin` and the sink radius `rs`, which is not on the envelope grid:
ρ(r) = ρ_s (r/r_s)^{-3/2}, inflow speed u_s (r/r_s)^{-1/2} (the free-fall shape, normalised to the envelope's density `rho_s` [g cm⁻³] and speed `u_s` [cm/s, positive inward] at r_s), H, He neutral
where not photoionised (the same opacities as on the envelope grid: ground-level bound–free of H I, He I, He II, `ion_sweep` equilibrium at gas temperature `T`), including the ionised flux that the inflow carries inwards (OI02 eq. 5).
The power absorbed in the column is re-emitted as a thermal continuum (steady state: everything absorbed leaves as recombination, free–free and degraded line radiation) — the caller adds it as a
blackbody of the precursor photosphere (`reemit = (g, Tph)` → added to `Lout`). Returns `Lout` (power per bin entering the envelope at r_s) and the diagnostics.
"""
function inner_column(ia::RT.IonAtoms, src::RT.Source, M, Rin, rs, rho_s, u_s; nz = 300, T = 1.5e4, passes = 6, reemit = nothing)
    Qs = qion(ia, src.L); Lsrc = sum(src.L)
    Mdot = 4π * rs^2 * rho_s * abs(u_s)
    αH = RT.ion_alpha(ia, T)[1]
    Qcrit = αH * (Mdot / mp)^2 / (8π * G * M) * 1.0            # OI02 eq. 8 (α Ṁ²/(8π G m_p² M)), Ṁ/m_p written with the proton mass
    Qcrit = αH * Mdot^2 / (8π * G * mp^2 * M)
    if rs <= 1.05 * Rin
        return (; Lout = copy(src.L), Qstar = Qs, Qesc = Qs, RII = Rin, Qcrit = Qcrit, Mdot = Mdot, nz = 0, Labs = 0.0, ledger = nothing, xedge = 0.0, nHin = 0.0)
    end
    r = Rin .* exp.(range(0, log(rs / Rin), length = nz + 1)); rc = sqrt.(r[1:end-1] .* r[2:end])
    rho = rho_s .* (rc ./ rs) .^ (-1.5); nH = rho ./ CN.μH()
    s = RT.ion_sweep(ia, src.L, r, fill(T, nz), nH; mdotH = fill(Mdot / CN.μH(), nz), passes = passes)
    kf = findlast(>=(0.5), s.x); RII = kf === nothing ? Rin : (kf == nz ? rs : exp(log(rc[kf]) + (s.x[kf] - 0.5) / (s.x[kf] - s.x[kf+1]) * (log(rc[kf+1]) - log(rc[kf]))))
    Lout = copy(view(s.Ld, nz + 1, :))
    Labs = Lsrc - sum(Lout)
    if reemit !== nothing
        g, Tph = reemit
        Labs > 0 && (Lout .+= RT.blackbody_source(g, Tph, Labs).L)
    end
    (; Lout, Qstar = Qs, Qesc = qion(ia, s.Ld[nz+1, :]), RII, Qcrit, Mdot, nz, Labs, ledger = s.ledger, xedge = s.x[end], nHin = nH[1])
end

# ------------------------------------------------------------------ fields on the envelope grid
"""
    envfields(g, ia, ab, net, photos, r, T, rho, mdotH, xnet, Lsrc; opts, ...) -> NamedTuple

One feedback evaluation on the hydro grid (n zones; r has n + 1 interfaces; `xnet[k]` the free-species vector of the ChemNet state of zone k). (1) The other opacities from the lagged state
(`RT.fill_opacity!`), (2) `ion_sweep` of the incoming spectrum `Lsrc` through the zones with those opacities: the equilibrium ionisation of H, He, He⁺ of the lit gas and the exactly photon-conserving
direct beam (the zone's state is raised to the sweep's wherever the sweep ionises more than the network has), (3) `RT.transfer` on the composition updated this way (direct beam as (2), the diffuse field,
H₂ Lyman–Werner rates, radiation acceleration). Returns `J` (n × nf, bin-wise mean intensity; + `Jcmb` if given), `kLW`, `arad` [cm s⁻²], `heat_other` [erg cm⁻³ s⁻¹] and diagnostics.
"""
function envfields(g::RT.FreqGrid, ia::RT.IonAtoms, ab, net, photos, r, T, rho, mdotH, xnet, Lsrc; opts = RT.Opts(diffuse = false), Tsweep_min = 8000.0, Jcmb = nothing, passes = 4, sweep = true, perf = nothing)
    n = length(T); nf = RT.nfreq(g); nH = rho ./ CN.μH()
    med = RT.medium_from_net(r, T, nH, xnet, net)
    X = med.X
    sec(name, f) = perf === nothing ? f() : perf(name, f)
    local sw = nothing
    xH0 = copy(get(X, :Hp, zeros(n)))
    op = RT.new_opacity(n, nf, length(ab))
    sec(:rt_opacity, () -> RT.fill_opacity!(op, g, ab, med, opts))
    function do_sweep(kext)
        Tsw = max.(T, Tsweep_min)
        sw = sec(:rt_sweep, () -> RT.ion_sweep(ia, Lsrc, r, Tsw, nH; kext = kext, mdotH = mdotH, passes = passes,
                                                x0 = (x = xH0, a1 = ia.yHe > 0 ? get(X, :Hep, zeros(n)) ./ ia.yHe : zeros(n), a2 = ia.yHe > 0 ? get(X, :Hepp, zeros(n)) ./ ia.yHe : zeros(n))))
        Hp = get!(X, :Hp, zeros(n)); Hn = get!(X, :H, zeros(n)); Ee = get!(X, :e, zeros(n)); Hep = get!(X, :Hep, zeros(n)); Hepp = get!(X, :Hepp, zeros(n)); Hen = get!(X, :He, zeros(n))
        for k in 1:n
            dx = max(sw.x[k] - Hp[k], 0.0)
            if dx > 0
                # the zone is ionised to the sweep's level: nuclei in neutral form = 1 − x; H₂ cannot hold more than 2 % of them (so that the neutral H atoms of the sweep stay in the medium)
                Hp[k] += dx; Ee[k] += dx; avail = max(1 - Hp[k], 0.0)
                H2k = get(X, :H2, nothing); h2n = H2k === nothing ? 0.0 : min(H2k[k], 0.01 * avail); H2k === nothing || (H2k[k] = h2n)
                Hn[k] = max(avail - 2h2n - (haskey(X, :Hm) ? X[:Hm][k] : 0.0) - (haskey(X, :H2p) ? 2 * X[:H2p][k] : 0.0), 0.0)
            end
            if ia.yHe > 0
                h1 = max(Hep[k], ia.yHe * sw.a1[k]); h2 = max(Hepp[k], ia.yHe * sw.a2[k])
                if h1 + h2 > ia.yHe; h1 = max(ia.yHe - h2, 0.0); end
                Ee[k] += (h1 - Hep[k]) + 2 * (h2 - Hepp[k]); Hep[k] = h1; Hepp[k] = h2; Hen[k] = max(ia.yHe - h1 - h2, 0.0)
            end
        end
        sw
    end
    if sweep && sum(Lsrc) > 0
        sw = do_sweep(op.kother .+ op.ssc)
        if any(k -> sw.x[k] > xH0[k] + 1e-3, 1:n)
            # the other opacities (Rayleigh of the neutral H and H₂ that are gone, …) depend on the composition: second pass with the opacities of the updated medium
            op = RT.new_opacity(n, nf, length(ab)); sec(:rt_opacity, () -> RT.fill_opacity!(op, g, ab, med, opts))
            sw = do_sweep(op.kother .+ op.ssc)
        end
    end
    rf = sec(:rt_transfer, () -> RT.transfer(g, ab, RT.Source(Lsrc, 0.0, sum(Lsrc)), med, opts; rho = rho))
    J = copy(rf.J); Jcmb === nothing || (J .+= reshape(Jcmb, 1, :))
    # net heating of the continua that the network does not carry (CIA, H⁻ and H free–free, …): absorption of the field minus the thermal emission of the same opacity, 4π ∫ k (J − B_ν(T)) dν
    # (one-sided absorption would heat the gas without the emission that balances it at J → B)
    emis = zeros(n)
    for k in 1:n, f in 1:nf; emis[k] += 4π * op.kother[k, f] * RT.Bν(g.ν[f], T[k]) * g.dν[f]; end
    heat_net = rf.heat_other .- emis
    (; J, kLW = rf.kLW, arad = rf.arad, heat_other = rf.heat_other, heat_net, rf, sw, med, xH0)
end

# ------------------------------------------------------------------ analytic estimates
"OI02 critical ionising rates for the state (M, Mdot, Rin): Q_crit^FF (eq. 8, free fall) and the radiation-force Q_crit (eq. 29, R_in ≪ R_II,crit) [s⁻¹]."
function oi02_qcrit(M, Mdot, Rin; T = 1e4, hν = 13.6 * eV, α = nothing)
    α === nothing && (α = 2.59e-13 * (1e4 / T)^0.7)
    QFF = α * Mdot^2 / (8π * G * mp^2 * M)
    RH = α * hν / (4π * c * G * mp^2)
    Qrad = (c / hν) * G * M^2 / (Rin * RH)
    (QFF = QFF, Qrad = Qrad, RH = RH)
end
"Eddington ratio for electron scattering of fully ionised H+He (X = 0.75): κ = 0.2 (1 + X)."
Γedd(L, M; X = 0.75) = 0.2 * (1 + X) * L / (4π * c * G * M)

"""
Upper bound of the force of trapped Lyα photons from the H II region (case B: 0.68 Lyα per recombination ≈ 0.68 Q): F ≤ (L_Lyα/c) N_sc, N_sc = min(τ_0-limited, 1/ε); ε = n_p q_{2p→2s}/A_{21}
(collisional mixing 2p → 2s and two-photon decay, q = 4.7e-4 cm³ s⁻¹ × (T/1e4)^{0.5}) at the proton density `np` of the H II gas. Returns (L_Lyα, N_sc, F [dyn]).
"""
function lya_force_bound(Q, np; T = 1.5e4, Nsc_cap = 1e8)
    LLy = 0.68 * Q * 10.2 * eV
    ε = np * 4.7e-4 * sqrt(T / 1e4) / 6.265e8
    Nsc = min(1 / max(ε, 1e-30), Nsc_cap)
    (LLy, Nsc, LLy / c * Nsc)
end
end # module
