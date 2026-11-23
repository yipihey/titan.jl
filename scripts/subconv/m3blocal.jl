# M3b (e): the local MLT limit inside the full code. For a TDC run with CONVDIAG: at the given snapshot step, every interface with a positive
# gate-free MLT flux of the same state (cd.fpot, 'fpot' = the Kippenhahn–Weigert flux of the buoyancy `sup`, ℓ, U of that interface) and a hydrostatic, subsonic
# state (mach_i < 0.1): F_c(TDC)/F_MLT, w/w_MLT, with F_c/F_tot and the superadiabaticity. Usage: julia m3blocal.jl <rundir> <step>
using Serialization, Printf
d = ARGS[1]; st = parse(Int, ARGS[2])
f = sort(filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(d)))[1]
open(joinpath(d, f)) do io
  while !eof(io)
    s = deserialize(io); s.step == st || continue; c = s.cd; n = length(s.ρ)
    @printf("step %d: interface, F_c [erg/s/cm2], F_MLT(same sup) , ratio, w [cm/s], w_MLT, ratio, sup, U, F_c/F_tot, mach\n", s.step)
    rows = Tuple[]
    for j in 2:n
      fp = c.fpot[j]; (isfinite(fp) && fp > 0) || continue
      c.mach_i[j] < 0.1 || continue
      push!(rows, (j, c.fconv[j], fp, c.fconv[j] / fp, c.tdc_w[j], c.vpot[j], c.tdc_w[j] / c.vpot[j], c.sup[j], c.tdc_U[j], c.Fconv_over_Ftot[j], c.mach_i[j]))
    end
    for r in rows; @printf("%4d %11.4e %11.4e %8.5f %11.4e %11.4e %8.5f %10.3e %10.3e %8.4f %6.3f\n", r...); end
    isempty(rows) && println("(no subsonic interface with F_MLT > 0)")
    fr = [r[4] for r in rows if r[10] > 0.01]
    isempty(fr) || @printf("interfaces with F_c/F_tot > 0.01: %d, F_c/F_MLT min %.5f median %.5f max %.5f\n", length(fr), minimum(fr), sort(fr)[(length(fr)+1)÷2], maximum(fr))
  end
end
