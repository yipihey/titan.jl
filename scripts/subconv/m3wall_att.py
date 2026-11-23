import re,sys
cur=None
for l in open(sys.argv[1]):
    if l.startswith('NWT head'):
        m=re.search(r'att (\d+) jstep (\d+) thet (\S+) dtime (\S+)',l); cur=m.groups()
    elif l.startswith('NWT end'):
        m=re.search(r'ok (\w+) iter (\d+)',l); print(cur[0],cur[1],cur[2],'%.3g'%float(cur[3]),m.groups())
