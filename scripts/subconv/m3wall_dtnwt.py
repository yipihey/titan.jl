import re,sys,statistics as st
# per jstep: dt of the accepted step from the last trapezoidal attempt (dtime = gamma*dt, gamma = 2-sqrt(2))
g=2-2**0.5; cur=None; last={}; fails=[]; nat=0
for l in open(sys.argv[1]):
    m=re.match(r'NWT head att (\d+) jstep (\d+) thet (\S+) dtime (\S+)',l)
    if m: cur=(int(m.group(2)),float(m.group(3)),float(m.group(4)))
    m=re.match(r'NWT end att (\d+) ok (\w+) iter (\d+)',l)
    if m:
        nat+=1
        if m.group(2)=='false': fails.append(cur[0])
        elif cur[1]==0.5: last[cur[0]]=cur[2]/g
lo=int(sys.argv[2]) if len(sys.argv)>2 else 0; hi=int(sys.argv[3]) if len(sys.argv)>3 else 10**9
js=sorted(j for j in last if lo<=j<=hi); dts=[last[j] for j in js]
nf=[j for j in fails if lo<=j<=hi]
print('steps %d (%d..%d): sum dt %.3e s, mean %.3e, median %.3e, min %.3e, max %.3e | Newton-failed attempts %d (%.1f per 100 steps)'%(len(js),js[0],js[-1],sum(dts),st.mean(dts),st.median(dts),min(dts),max(dts),len(nf),100*len(nf)/len(js)))
