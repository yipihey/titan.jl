# T, u, log rho at fixed radii over selected frames of a RUNREC export:  fixr.py <exportdir> <logr,...> <frame indices,... | first:N>
import json, sys, numpy as np
d=sys.argv[1]; lr=[float(x) for x in sys.argv[2].split(",")]
m=json.load(open(d+"/frames.json")); raw=np.fromfile(d+"/frames.f32",dtype=np.float32); nm=[f[0] for f in m["fields"]]
F=[]
for fr in m["frames"]:
    n,o=fr["n"],fr["off"]//4; a=raw[o:o+n*len(nm)].reshape(len(nm),n).astype(float); F.append(({k:a[i] for i,k in enumerate(nm)}, dict(zip(m["scalars"],fr["s"]))))
sel=sys.argv[3]; idx=list(range(int(sel[6:]))) if sel.startswith("first:") else [int(x) for x in sel.split(",")]
idx=[i if i>=0 else len(F)+i for i in idx if -len(F)<=i<len(F)]
print(len(F),"frames; log r =",lr); t0=F[idx[0]][1]["t_yr"]
for i in idx:
    f,s=F[i]; g=lambda q:[np.interp(x,f["logr"],f[q]) for x in lr]
    print(f"step {int(s['step']):6d} t-t0={s['t_yr']-t0:10.1f} yr  nz {len(f['logr'])}")
    print("     T: "+" ".join(f"{10**x:8.1f}" for x in g("logT"))); print("     u: "+" ".join(f"{x:8.3f}" for x in g("u"))); print("  lrho: "+" ".join(f"{x:8.3f}" for x in g("logrho")))
