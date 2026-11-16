import json, base64, numpy as np, sys
d = '/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/panel/zdiagB/'
J = json.load(open(d + 'frames.json'))
raw = base64.b64decode(open(d + 'frames.b64.txt').read().strip())
A = np.frombuffer(raw, dtype=np.float32)
F = [f[0] for f in J['fields']]


def get(fr, name):
    n = fr["n"]; off = fr["off"] // 4; k = F.index(name)
    return A[off + k * n:off + (k + 1) * n]


if __name__ == '__main__':
    print(len(J['frames']), J['scalars'])
    sel = [int(a) for a in sys.argv[1:]] or [len(J['frames']) - 1]
    for fi in sel:
        fr = J['frames'][fi]
        print(fi, dict(zip(J['scalars'], fr['s'])))
        lr = get(fr, 'logr'); rho = get(fr, 'logrho'); T = get(fr, 'logT'); Tr = get(fr, 'logTrad')
        xp = get(fr, 'xHp'); x2 = get(fr, 'xH2'); tau = get(fr, 'logtau'); xm = get(fr, 'logxHm')
        n = fr['n']
        # cliff: largest |dlogT| between neighbours with T > 3000 K
        dT = np.abs(np.diff(T)); dT[(T[:-1] < 3.4)] = 0
        i = int(np.argmax(dT))
        for j in range(max(0, i - 6), min(n, i + 8)):
            print(j, f"logr {lr[j]:.3f} logrho {rho[j]:.3f} logT {T[j]:.3f} logTr {Tr[j]:.3f} xHp {xp[j]:.3e} xH2 {x2[j]:.3e} logxHm {xm[j]:.2f} logtau {tau[j]:.2f}")
