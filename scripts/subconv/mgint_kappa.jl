# M3e: the eddy opacity under MGROUPS (group Rosseland-type mean, cd.kap) against the grey avchi (cd.kap_grey), F_rad as the group sum, and the terms of the gas-energy row
# at the interfaces/zones with turbulence. Usage: julia mgint_kappa.jl <rundir> <step>   (CONVDIAG=1 run, snapshot `step`)
using Serialization, Printf
d = ARGS[1]; st = parse(Int, ARGS[2])
for f in sort(filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(d)))
  open(joinpath(d, f)) do io
    while !eof(io)
      s = deserialize(io); s.step == st || continue; c = s.cd
      n = length(s.ρ); w = c.tdc_w
      @printf("step %d: interfaces with w > 1e3 cm/s\n  k      w[cm/s]   kap(group)  kap(grey)  ratio   Fconv/Ftot  frad[erg/s/cm2]  sup\n", s.step)
      for j in 1:n+1
        (isnan(w[j]) || w[j] < 1e3) && continue
        @printf("%4d %11.3e %11.4e %11.4e %7.4f %10.3g %14.5e %10.3e\n", j, w[j], c.kap[j], c.kap_grey[j], c.kap[j] / c.kap_grey[j], c.Fconv_over_Ftot[j], c.frad[j], c.sup[j])
      end
      nm = (:e_stor, :e_pdv, :e_av, :e_kin, :e_adv, :e_dfe, :e_rad, :e_conv, :e_chem, :e_src, :e_tdc); mx = 0.0; zw = 0
      for j in 1:n
        isnan(c.e_res[j]) && continue
        sm = sum(c[k][j] for k in nm if !isnan(c[k][j])); tm = maximum(abs(c[k][j]) for k in nm if !isnan(c[k][j]))
        r = abs(sm - c.e_res[j]) / max(tm, 1e-300); r > mx && (mx = r; zw = j)
      end
      @printf("gas-energy row: max |sum of the 11 terms - e_res| / (largest term) over zones = %.3e (zone %d)\n", mx, zw)
    end
  end
end
