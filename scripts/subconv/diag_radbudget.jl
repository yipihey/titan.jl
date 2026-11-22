# zone radiative exchange (e_rad, total = 4 pi x per-sr) versus the L_rad profile, zie1 k655 step 700: julia diag_radbudget.jl <snaps.jls>
using Serialization, Printf
v = Any[]; open(ARGS[1]) do io; while !eof(io); push!(v, deserialize(io)); end; end
s = v[findlast(x -> x.step == 700, v)]; cd = s.cd
for k in 39:42
  mz = 4pi/3 * s.ρ[k] * (s.r[k+1]^3 - s.r[k]^3)
  # e_rad is per sr: total = 4pi*e_rad ; sign: positive = cooling
  @printf("zone %d mass %.3f Msun  e_rad(total) %.3e erg/s = %.1f Lsun (positive = net cooling)  e_adv %.3e = %.1f Lsun\n", k, mz/1.989e33, 4pi*cd.e_rad[k], 4pi*cd.e_rad[k]/3.828e33, 4pi*cd.e_adv[k], 4pi*cd.e_adv[k]/3.828e33)
end
for j in 40:43; @printf("j=%d L_rad=%.1f Lsun\n", j, 4pi*s.r[j]^2*cd.frad[j]/3.828e33); end
