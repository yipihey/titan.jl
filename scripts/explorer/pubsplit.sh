#!/bin/bash
# Split the PUBLISH lines of a watch5 output ($1) into pubA.json (runs.json first, then files, <= 45 MB), pubB.json, pubC.json ... (<= 45 MB each)
S=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad
grep '^PUBLISH:' "$1" | sed 's/^PUBLISH: //' | tr ' ' '\n' | grep -v '^$' > $S/publist.txt
rm -f $S/pubC.json $S/pubD.json 2>/dev/null
cd $S/panel && python3 -I -c '
import os,json
S="/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad"
L=[l.strip() for l in open(S+"/publist.txt") if l.strip() and os.path.exists(l.strip())]   # (a run just created may not have every file yet)
# a data file and its metadata (run/spec.json + run/spec.b64.txt, run/frames.json + run/frames.b64.txt) always go in the same publish:
# split across two publishes, the page could read new metadata with old data in between
units=[]; idx={}
for f in L:
    key=f.split(".")[0]
    if key in idx: units[idx[key]].append(f)
    else: idx[key]=len(units); units.append([f])
groups=[[]]; size=[0]
for u in units:
    s=sum(os.path.getsize(f) for f in u)
    if size[-1]+s>45e6 and groups[-1]: groups.append([]); size.append(0)
    groups[-1].extend(u); size[-1]+=s
names="ABCDEFG"
for i in range(max(2,len(groups))):
    g=groups[i] if i<len(groups) else []
    json.dump({f:f for f in g},open(S+"/pub%s.json"%names[i],"w"))
print(" ".join("%s %.1f MB %d files"%(names[i],size[i]/1e6,len(groups[i])) for i in range(len(groups))))'
