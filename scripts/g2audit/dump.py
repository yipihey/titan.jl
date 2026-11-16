import sys
sys.path.insert(0, '/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2')
import frames as F
fi, z0, z1, path = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
fr = F.J['frames'][fi]
cols = ('logrho', 'logT', 'logTrad', 'xHp', 'xH2', 'logxHm', 'logtau')
with open(path, 'w') as f:
    for k in range(z0, z1 + 1):
        f.write(' '.join([str(k)] + [repr(float(F.get(fr, c)[k])) for c in cols]) + '\n')
