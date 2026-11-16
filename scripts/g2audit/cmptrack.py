import sys, numpy as np
def load(p):
    rows = []
    for l in open(p):
        s = l.split()
        if len(s) > 10 and s[0].isdigit():
            rows.append(dict(step=int(s[0]), dt=float(s[2]), M=float(s[3]), R=float(s[4]), L=float(s[5]), Te=float(s[6]), Tsl=float(s[7]),
                             on=('YES' in l), dlgT=float(s[10])))
    return rows
if __name__ == '__main__':
    smax = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
    for p in sys.argv[1:3]:
        r = [x for x in load(p) if x['step'] <= smax]
        st = np.array([x['step'] for x in r]); dt = np.array([x['dt'] for x in r]); L = np.array([x['L'] for x in r]); Te = np.array([x['Te'] for x in r])
        late = st >= 3040
        last = r[-1]
        print(f"{p}: snaps {len(r)} steps {st[0]}-{st[-1]}; dt median {np.median(dt):.3g} (from 3040: median {np.median(dt[late]):.3g}, min {dt[late].min():.3g}, "
              f"frac >= 1e4 s {np.mean(dt[late] >= 1e4):.2f}); <log L_ph> {np.mean(np.log10(L[late])):.3f} (std {np.std(np.log10(L[late])):.3f}); "
              f"<T_eff,L> {np.mean(Te[late]):.0f} (std {np.std(Te[late]):.0f}); at step {last['step']}: M* {last['M']:.5f} R* {last['R']:.3f}; "
              f"cliff on interface {sum(x['on'] for x in r)}/{len(r)}; <dlogT cliff> {np.mean([x['dlgT'] for x in r]):.3f}")
