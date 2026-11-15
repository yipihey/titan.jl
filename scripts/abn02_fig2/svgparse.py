import re, sys, json
def parse(fn):
    s = open(fn).read()
    body = s[s.index('<g id="surface1">'):]
    out = []
    for m in re.finditer(r'<path ([^>]*)/>', body):
        attrs = m.group(1)
        st = re.search(r'style="([^"]*)"', attrs); d = re.search(r' d="([^"]*)"', attrs); tr = re.search(r'transform="matrix\(([^)]*)\)"', attrs)
        if d is None: continue
        style = dict(kv.split(':', 1) for kv in st.group(1).strip().strip(';').split(';') if ':' in kv) if st else {}
        T = [float(x) for x in tr.group(1).split(',')] if tr else [1, 0, 0, 1, 0, 0]
        toks = re.findall(r'[MLCZmlcz]|-?[0-9.]+(?:e-?[0-9]+)?', d.group(1))
        subs = []; cur = []; cmd = None; nums = []
        def flush():
            pass
        i = 0
        pts_cmds = []
        while i < len(toks):
            t = toks[i]
            if t in 'MLCZ':
                cmd = t; i += 1
                if t == 'Z': pts_cmds.append(('Z', None))
                continue
            if cmd in ('M', 'L'):
                x, y = float(toks[i]), float(toks[i+1]); i += 2; pts_cmds.append((cmd, (x, y)))
            elif cmd == 'C':
                xs = [float(toks[i+k]) for k in range(6)]; i += 6; pts_cmds.append(('C', (xs[4], xs[5]), xs))
            else:
                raise ValueError(cmd)
        def tf(p):
            a, b, c, dd, e, f = T
            return (a * p[0] + c * p[1] + e, b * p[0] + dd * p[1] + f)
        sp = []
        for pc in pts_cmds:
            if pc[0] == 'M':
                if cur: sp.append(cur)
                cur = [tf(pc[1])]
            elif pc[0] in ('L', 'C'):
                cur.append(tf(pc[1]))
        if cur: sp.append(cur)
        ncurve = sum(1 for pc in pts_cmds if pc[0] == 'C')
        out.append(dict(style=style, subpaths=sp, ncurve=ncurve, line=s[:m.start() + s.index('<g id="surface1">')].count('\n') + 1))
    return out
if __name__ == '__main__':
    P = parse(sys.argv[1])
    for p in P:
        allp = [q for sp in p['subpaths'] for q in sp]
        xs = [q[0] for q in allp]; ys = [q[1] for q in allp]
        st = p['style']
        print(p['line'], st.get('stroke', '-')[:22], st.get('fill', '-')[:22], st.get('stroke-dasharray', ''), st.get('stroke-width', ''), len(p['subpaths']), len(allp), p['ncurve'],
              '%.1f-%.1f x %.1f-%.1f' % (min(xs), max(xs), min(ys), max(ys)))
