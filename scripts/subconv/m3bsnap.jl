# M3b: per-snapshot summary of a CONVDIAG=1 run directory (snaps_*.jls, the first process segment): step, dt, R*, the convective layer, F_c/F_tot, w, P_t/P, w in the supersonic inflow, Ṁ.
# Usage: julia m3bsnap.jl <rundir> [every]    (Float64 only; reads `cd` as written by popiii_m2c convdiag_snap)
using Serialization, Printf
d = ARGS[1]; ev = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
fs = sort(filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(d)))
nz(v) = [isnan(q) ? 0.0 : q for q in v]
@printf("%6s %9s %7s %8s | %9s %6s %6s | %9s %5s | %9s %9s %9s | %9s\n", "step", "dt", "R*/Rsun", "Mdot[Msun/yr]", "maxFc/Ft", "@int", "maxw/cs", "w_max", "@int", "Pt/P ps", "maxPt/P", "w_sup", "E_t")
for f in fs
  open(joinpath(d, f)) do io
    while !eof(io)
      s = deserialize(io)
      s.step % ev == 0 || s.step <= 2 || continue
      c = s.cd; c === nothing && continue
      n = length(s.ρ); u = s.u; r = s.r
      js = findfirst(<(-2e5), u); Rst = js === nothing ? NaN : r[js] / 6.957e10
      # accretion rate through the outermost-inflow interface just outside R*: Mdot = 4π r² ρ |u| at interface js+3 (ρ of the zone above), Msun/yr
      Md = NaN
      if js !== nothing && js + 3 <= n
        jj = js + 3; Md = 4π * r[jj]^2 * s.ρ[jj] * abs(u[jj]) * 3.156e7 / 1.989e33
      end
      fo = nz(c.Fconv_over_Ftot); jm = argmax(fo)
      w = hasproperty(c, :tdc_w) ? nz(c.tdc_w) : zeros(n + 1); jw = argmax(w)
      cs = sqrt.(5 / 3 .* s.pg ./ s.ρ); csi = [j == 1 || j == n + 1 ? NaN : 0.5 * (cs[j-1] + cs[j]) for j in 1:n+1]
      wcs = maximum(nz(w ./ csi))
      mach = nz(c.mach_i); sup = [(mach[j] > 1 && u[j] < 0) ? w[j] : 0.0 for j in 1:n+1]
      ptp = hasproperty(c, :tdc_PtP) ? nz(c.tdc_PtP) : zeros(n)
      sc = hasproperty(c, :tdc_scal) ? c.tdc_scal : [NaN, NaN, NaN]
      @printf("%6d %9.3e %7.1f %8.2e | %9.3g %6d %6.3g | %9.3e %5d | %9.3g %9.3g %9.3e | %9.3e\n", s.step, s.dt, Rst, Md, fo[jm], jm, wcs, w[jw], jw, sc[2], sc[3], maximum(sup), sc[1])
    end
  end
end
