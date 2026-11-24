# M6: per-snapshot summary of a CONVDIAG=1 TDC pilot: step, dt, max F_c/F_tot and its interface, w there, sup, Re Γ₁,b, Γ₁(0), Γ₁(∞), ωτ_max (M6 runs only) at the interface of the largest sup.
# Usage: julia lresp_pilot.jl <rundir> [every]
using Serialization, Printf
d = ARGS[1]; ev = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
fs = sort(filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(d))); nz(v) = [isnan(q) ? 0.0 : q for q in v]
@printf("%6s %9s | %9s %5s %9s | %9s %5s | %8s %8s %8s %9s\n", "step", "dt", "maxFc/Ft", "@int", "w[cm/s]", "max sup", "@int", "Γ1b", "Γ1(0)", "Γ1(inf)", "ωτmax")
for f in fs
  open(joinpath(d, f)) do io
    while !eof(io)
      s = deserialize(io); (s.step % ev == 0 || s.step <= 2) || continue
      c = s.cd; c === nothing && continue
      fo = nz(c.Fconv_over_Ftot); jm = argmax(fo); w = hasproperty(c, :tdc_w) ? nz(c.tdc_w) : zeros(length(fo))
      sp = [isnan(q) ? -Inf : q for q in c.sup]; js = argmax(sp)
      g(nm) = hasproperty(c, nm) ? c[nm][js] : NaN
      @printf("%6d %9.3e | %9.3g %5d %9.3e | %9.3g %5d | %8.4f %8.4f %8.4f %9.2e\n", s.step, s.dt, fo[jm], jm, w[jm], sp[js], js, g(:gam1b), g(:gam1_0), g(:gam1_inf), g(:wtau))
    end
  end
end
