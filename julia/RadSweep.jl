# S7b: equilibrium inside-out photoionisation sweep for H, He, He⁺ and the interface radiation acceleration. Included at the end of the module RadTransfer
# (uses its FreqGrid, PhotoQuad, quad_milne, absorbers, volume, CN, hP).

"""
Atomic data of the sweep (`ion_atoms`): bin-averaged ground-level cross sections of H I, He I, He II, ionisation energies, and the recombination coefficients to the
excited levels (= case B, the Milne coefficient of the network's `excited` photo pair at J = 0) tabulated in log T; `dr(T)` adds the dielectronic recombination of He⁺.
"""
struct IonAtoms
    g::FreqGrid
    sig::NTuple{3,Vector{Float64}}; Eion::NTuple{3,Float64}
    tab::Matrix{Float64}; lT0::Float64; dlT::Float64
    yHe::Float64; dr::Function
    f0::Int                                              # first bin above the H edge
end
function ion_atoms(photos, g::FreqGrid; yHe = CN.yHe, dr = T -> 0.0)
    syms = (:H, :He, :Hep); sig = Vector{Float64}[]; Ei = Float64[]; tab = zeros(301, 3)
    for (i, s) in enumerate(syms)
        gr = [p for p in photos if p.R === s && length(p.levels) == 1]; ex = [p for p in photos if p.R === s && length(p.levels) > 1]
        (length(gr) == 1 && length(ex) == 1) || error("ion_atoms: need one ground and one excited Photo for $s")
        a = absorbers(gr, g)[1]; push!(sig, a.σ[1, :]); push!(Ei, a.Eion)
        pq = PhotoQuad(ex[1]; n = 160)
        for j in 1:301; tab[j, i] = quad_milne(pq, 10^(1.0 + 0.02 * (j - 1)))[1]; end
    end
    f0 = findfirst(i -> g.edges[i] >= CN.χH / hP * (1 - 1e-9), 1:nfreq(g))
    IonAtoms(g, (sig[1], sig[2], sig[3]), (Ei[1], Ei[2], Ei[3]), tab, 1.0, 0.02, yHe, dr, f0)
end
"Recombination coefficients (α_H, α_He = excited + dielectronic, α_He⁺) at T [cm³ s⁻¹] (log-linear interpolation of the table)."
function ion_alpha(ia::IonAtoms, T)
    x = clamp((log10(T) - ia.lT0) / ia.dlT, 0.0, 299.999); j = floor(Int, x) + 1; w = x - (j - 1)
    a(i) = exp((1 - w) * log(max(ia.tab[j, i], 1e-300)) + w * log(max(ia.tab[j+1, i], 1e-300)))
    (a(1), a(2) + ia.dr(T), a(3))
end

"Root of a function increasing in u on [lo, hi] (Illinois); the end point is returned when F does not change sign."
function root_inc(F, lo, hi; rtol = 1e-12, maxit = 80)
    flo = F(lo); flo >= 0 && return lo
    fhi = F(hi); fhi <= 0 && return hi
    side = 0; mid = hi
    for _ in 1:maxit
        mid = (lo * fhi - hi * flo) / (fhi - flo)
        (mid <= lo || mid >= hi) && (mid = 0.5 * (lo + hi))
        fm = F(mid)
        fm == 0 && return mid
        if fm < 0; lo = mid; flo = fm; side == -1 && (fhi *= 0.5); side = -1
        else hi = mid; fhi = fm; side == 1 && (flo *= 0.5); side = 1; end
        hi - lo < rtol * max(abs(hi), 1e-300) && return 0.5 * (lo + hi)
    end
    mid
end

"""
    ion_sweep(ia, Lsrc, r, T, nH; kext = nothing, mdotH = nothing, x0 = nothing, Tfun = nothing, passes = 3, tol = 1e-9, ne_extra = nothing) -> NamedTuple

Steady photoionisation equilibrium of H, He, He⁺ solved zone by zone from the source outwards with exact photon conservation (the generalisation of `HII.sweep!`
to three absorbers; equilibrium instead of the time-dependent ODE, i.e. the R-type front that crosses many zones per call is found directly). `Lsrc` is the power per
bin entering the first zone [erg/s], `r` the n + 1 interfaces, `T`, `nH` per zone (H nuclei cm⁻³), `kext[k, f]` the other absorption [cm⁻¹] (Thomson, H⁻, Rayleigh, …: absorbs
photons without ionising). In each zone (H⁺ fraction x, He⁰/He⁺/He²⁺ fractions a0, a1, a2) the absorbed photon numbers A_i of the species, with the optical depth τ_f = Δr Σ n_i σ_if
evaluated at the zone's own state, balance the case-B recombinations: A_H = α_H n_e n_p V + Ṁ_H (x_k − x_{k+1}) (the last term the ionised flux that the inflow carries inwards, optional,
`mdotH` = H nuclei per second through the zone, neutral gas arrives from the outer neighbour; iterated `passes` times), A_He⁰ = (α_He + α_dr) n_e n_He⁺ V, A_He⁺ = α_He⁺ n_e n_He²⁺ V;
the transmitted power is L_in e^{−τ}. Recombination photons are absorbed on the spot (case B; the He recombination photons that would ionise H are neglected).
`Tfun(G, x, nH, T)` (optional) returns the gas temperature for the photoheating rate G [erg cm⁻³ s⁻¹] of the zone (iterated with the state).
Returns `(x, a1, a2, T, Ld (n+1 × nf), Jdir (n × nf), A (n × 3 absorbed photons/s), G, ledger)`.
"""
function ion_sweep(ia::IonAtoms, Lsrc::AbstractVector, r::AbstractVector, T0::AbstractVector, nH::AbstractVector; kext = nothing, mdotH = nothing, x0 = nothing,
                   Tfun = nothing, passes = 3, tol = 1e-9, ne_extra = nothing)
    g = ia.g; nf = nfreq(g); n = length(nH); f0 = ia.f0; yHe = ia.yHe
    hν = hP .* g.ν
    xs_ = x0 === nothing ? zeros(n) : copy(x0.x); a1s = x0 === nothing ? zeros(n) : copy(x0.a1); a2s = x0 === nothing ? zeros(n) : copy(x0.a2)
    Ts = collect(Float64, T0); Ld = zeros(n + 1, nf); Jdir = zeros(n, nf); A = zeros(n, 3); Gv = zeros(n); Nadv = zeros(n)
    σH = ia.sig[1]; σHe = ia.sig[2]; σHp = ia.sig[3]
    for pass in 1:(mdotH === nothing ? 1 : passes)
        Lin = collect(Float64, Lsrc); Ld[1, :] .= Lin
        xold = copy(xs_)
        for k in 1:n
            V = volume(r, k); Δr = r[k+1] - r[k]; nHk = nH[k]; T = Ts[k]
            Qf = [Lin[f] / hν[f] for f in f0:nf]                       # photons/s in the ionising bins
            kx = kext === nothing ? nothing : view(kext, k, :)
            xin = mdotH === nothing ? 0.0 : (pass > 1 ? xold[min(k + 1, n)] : 0.0)
            MH = mdotH === nothing ? 0.0 : mdotH[k]
            nex = ne_extra === nothing ? 0.0 : ne_extra[k]
            uH = 1.0 - xs_[k]; a1 = a1s[k]; a2 = a2s[k]; a0 = max(1.0 - a1 - a2, 0.0)
            Gk = 0.0
            absorb(uH_, a0_, a1_) = begin
                AH = 0.0; A0 = 0.0; A1 = 0.0
                @inbounds for j in 1:nf-f0+1
                    f = f0 + j - 1
                    kH = nHk * uH_ * σH[f]; k0 = nHk * yHe * a0_ * σHe[f]; k1 = nHk * yHe * a1_ * σHp[f]; ko = kx === nothing ? 0.0 : kx[f]
                    kt = kH + k0 + k1 + ko; kt <= 0 && continue
                    ab = Qf[j] * (-expm1(-kt * Δr)) / kt
                    AH += ab * kH; A0 += ab * k0; A1 += ab * k1
                end
                (AH, A0, A1)
            end
            nefun(uH_, a1_, a2_) = nHk * ((1 - uH_) + yHe * (a1_ + 2a2_)) + nex
            for tit in 1:(Tfun === nothing ? 1 : 5)
                αH, αHe, αHp = ion_alpha(ia, T)
                for it in 1:300
                    uHo, a0o, a2o = uH, a0, a2
                    FH(u) = (ne = nefun(u, a1, a2); absorb(u, a0, a1)[1] - αH * ne * nHk * (1 - u) * V - MH * ((1 - u) - xin))
                    uH = root_inc(FH, 1e-16, 1.0)
                    if yHe > 0
                        F1(a0_) = (a1_ = 1 - a2 - a0_; ne = nefun(uH, a1_, a2); absorb(uH, a0_, a1_)[2] - αHe * ne * nHk * yHe * a1_ * V)
                        a0 = root_inc(F1, 0.0, max(1 - a2, 1e-16)); a1 = max(1 - a2 - a0, 0.0)
                        F2(a2_) = (a1_ = 1 - a0 - a2_; ne = nefun(uH, a1_, a2_); -(absorb(uH, a0, a1_)[3] - αHp * ne * nHk * yHe * a2_ * V))
                        a2 = root_inc(F2, 0.0, max(1 - a0, 1e-16)); a1 = max(1 - a0 - a2, 0.0)
                    else
                        a0 = 1.0; a1 = 0.0; a2 = 0.0
                    end
                    abs(uH - uHo) <= tol * max(uH, 1e-12) && abs(a0 - a0o) <= tol * max(a0, 1e-12) + 1e-14 && abs(a2 - a2o) <= tol + 1e-14 && break
                end
                if Tfun !== nothing
                    G = 0.0
                    @inbounds for j in 1:nf-f0+1
                        f = f0 + j - 1
                        kH = nHk * uH * σH[f]; k0 = nHk * yHe * a0 * σHe[f]; k1 = nHk * yHe * a1 * σHp[f]; ko = kx === nothing ? 0.0 : kx[f]
                        kt = kH + k0 + k1 + ko; kt <= 0 && continue
                        ab = Qf[j] * (-expm1(-kt * Δr)) / kt
                        G += ab * (kH * (hν[f] - ia.Eion[1]) + k0 * (hν[f] - ia.Eion[2]) + k1 * (hν[f] - ia.Eion[3]) + ko * hν[f])
                    end
                    Gk = G / V
                    Tn = Tfun(Gk, 1 - uH, nHk, T); abs(Tn / T - 1) < 1e-4 && (T = Tn; break); T = sqrt(T * Tn)
                end
            end
            Ts[k] = T; xs_[k] = 1 - uH; a1s[k] = a1; a2s[k] = a2
            AH, A0, A1 = absorb(uH, a0, a1)
            A[k, 1] = AH; A[k, 2] = A0; A[k, 3] = A1; Gv[k] = Gk
            Nadv[k] = MH * ((1 - uH) - xin)
            for f in 1:nf
                kt = (f >= f0 ? nHk * uH * σH[f] + nHk * yHe * a0 * σHe[f] + nHk * yHe * a1 * σHp[f] : 0.0) + (kext === nothing ? 0.0 : kext[k, f])
                τ = kt * Δr
                Jdir[k, f] = (τ > 1e-300 ? Lin[f] * (-expm1(-τ)) / τ : Lin[f]) * Δr / (4π * V) / g.dν[f]
                Lin[f] *= exp(-τ); Ld[k+1, f] = Lin[f]
            end
        end
    end
    Nsrc = sum(Lsrc[f] / hν[f] for f in f0:nf); Nesc = sum(Ld[n+1, f] / hν[f] for f in f0:nf)
    Nrec = 0.0; Nion = sum(A)
    for k in 1:n
        V = volume(r, k); αH, αHe, αHp = ion_alpha(ia, Ts[k]); nHk = nH[k]
        ne = nHk * (xs_[k] + yHe * (a1s[k] + 2a2s[k])) + (ne_extra === nothing ? 0.0 : ne_extra[k])
        Nrec += (αH * ne * nHk * xs_[k] + αHe * ne * nHk * yHe * a1s[k] + αHp * ne * nHk * yHe * a2s[k]) * V
    end
    (x = xs_, a1 = a1s, a2 = a2s, T = Ts, Ld = Ld, Jdir = Jdir, A = A, G = Gv,
     ledger = (Nsrc = Nsrc, Nesc = Nesc, Nion = Nion, Nother = Nsrc - Nesc - Nion, Nrec = Nrec, Nadv = sum(Nadv)))
end

"""
    interface_accel(a, m) -> Vector (n + 1)

Per-interface acceleration for `TitanResidFast.AEXT` from zone accelerations `a` and zone masses `m` [g]: interface j carries half of the force a_{j−1} m_{j−1} and half of a_j m_j on
the interface mass ½(m_{j−1} + m_j) (the end interfaces take the adjacent zone's acceleration); the force is conserved, Σ a_j M_j = Σ a_k m_k up to the end interfaces.
"""
function interface_accel(a::AbstractVector, m::AbstractVector)
    n = length(a); out = zeros(n + 1)
    out[1] = a[1]; out[n+1] = a[n]
    for j in 2:n
        M = 0.5 * (m[j-1] + m[j]); out[j] = 0.5 * (a[j-1] * m[j-1] + a[j] * m[j]) / M
    end
    out
end
