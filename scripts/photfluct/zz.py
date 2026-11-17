# period-2 (zig-zag) amplitude over a step window: rms of x_k - (x_{k-1}+x_{k+1})/2 in percent, and fraction of sign alternations
import sys, numpy as np
s0, s1 = int(sys.argv[1]), int(sys.argv[2])
for run in sys.argv[3:]:
    d = np.genfromtxt('hist_%s.csv' % run, delimiter=',', names=True); d = d[(d['step'] >= s0) & (d['step'] <= s1)]
    out = []
    for c in ['Lph', 'Lem', 'Trad_ph', 'Tgas_ph']:
        x = np.log(np.maximum(d[c], 1e-30)); z = x[1:-1] - 0.5 * (x[:-2] + x[2:])
        a = np.mean(np.sign(z[1:]) != np.sign(z[:-1]))
        out.append('%s zz %.3f%% alt %.2f' % (c, 100 * np.sqrt(np.mean(z ** 2)), a))
    print('%-7s n %d dt med %.3g | ' % (run, len(d), np.median(d['dt'])) + ' | '.join(out) + ' | M %.5f R %.2f Lem mean %.3g' % (d['Mst'][-1], d['Rst'][-1], d['Lem'].mean()))
