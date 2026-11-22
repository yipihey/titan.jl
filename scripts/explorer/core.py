# core.py <exportdir> <logn_c,...> : central structure at the frames nearest given central densities
import sys, json, numpy as np
ex=sys.argv[1]; meta=json.load(open(ex+"/frames.json")); raw=np.fromfile(ex+"/frames.f32",dtype=np.float32)
names=[f[0] for f in meta["fields"]]; sc=meta["scalars"]; iln=sc.index("logn_c")
ln=np.array([fr["s"][iln] for fr in meta["frames"]],dtype=float)
for n in map(float,sys.argv[2].split(",")):
    i=int(np.argmin(abs(ln-n))); fr=meta["frames"][i]; nz=fr["n"]; o=fr["off"]//4
    a=raw[o:o+nz*len(names)].reshape(len(names),nz); F={k:a[j] for j,k in enumerate(names)}
    s=dict(zip(sc,fr["s"]))
    print(f"--- step {int(s['step'])} log n_c {s['logn_c']:.3f} T_c {s['T_c']:.1f} L_lines {s['L_lines']:.3e}")
    print("  zone  logm    logr   logrho   T      Trad    logtau  xH2    L[Lsun]   logLam   FcE     eq  u")
    tau=F["logtau"]
    idx=sorted(set([0,2,5]+[int(k) for k in np.linspace(8,nz-1,26)]))
    for k in idx:
        print(f"  {k:4d} {F['logm'][k]:7.3f} {F['logr'][k]:7.3f} {F['logrho'][k]:7.3f} {10**F['logT'][k]:8.1f} {10**F['logTrad'][k]:8.1f} {tau[k]:7.2f} {F['xH2'][k]:6.3f} {F['L'][k]:10.3e} {F['logLam'][k]:7.2f} {F['FcE'][k]:8.1e} {int(F['eq'][k]):2d} {F['u'][k]:7.3f}")
