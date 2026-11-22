import sys
p = sys.argv[1]; s = open(p).read()
old = "    TX.GASROW[] = true\n    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lum = ENVCENTRE ? nothing : Lin, lgrid = 2,"
assert s.count(old) == 1
new = ("    TX.GASROW[] = true\n"
       "    # per-zone material on the new indices before build! fills the EOS (as in the ENVMERGE rebuild: else a new zone far out inherits the\n"
       "    # table-EOS flag of the old zone with its index: \"rho below the table\" when the source has more zones in the star than the target\n"
       "    # has in all, e.g. a refined stage A3 onto 64 zones); set again below on the rebuilt grid\n"
       "    for k in eachindex(mat.x); j = clamp(k - a + 1, 1, n); mat.x[k] = copy(xz[j]); mat.eqzone[k] = eqz[j]; end\n"
       "    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr, lum = ENVCENTRE ? nothing : Lin, lgrid = 2,")
open(sys.argv[2], "w").write(s.replace(old, new))
