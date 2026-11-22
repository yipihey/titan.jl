#!/bin/bash
# Refresh the explorer data: runs.json for all runs, frames + perf for every run whose chain is running (live process, perf.csv
# updated in the last 30 min). Prints the published paths of all files to send with the next Artifact publish (PUBLISH: lines).
S=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad; PY=/home/tabel/.local/share/anaconda3/bin/python3
$PY -I $S/panel_manifest.py > $S/panel_manifest.log 2>&1 || { echo "manifest failed"; cat $S/panel_manifest.log; exit 1; }
RUNNING=$($PY -I -c "import json;print(' '.join(r['key'] for r in json.load(open('$S/panel/runs.json'))['runs'] if r['status']=='running'))")
# CHANGED: runs that are not running any more but whose perf.csv is newer than their explorer frames (finished or failed in the
# last 24 h, e.g. a test that ended between refreshes): their frames are rebuilt and published too, so every recent run is visible
CHANGED=$($PY -I -c "
import json,os,time
S='$S'; D='/zpool/nvme/data/tabel_scratch/titan_popiii/runs'
for r in json.load(open(S+'/panel/runs.json'))['runs']:
    k=r['key']; p=os.path.join(D,k,'perf.csv'); f=os.path.join(S,'panel',k,'frames.json')
    if r['status']=='running' or not os.path.exists(p) or not any(n.startswith('snaps_') for n in os.listdir(os.path.join(D,k))): continue
    m=os.path.getmtime(p)
    if time.time()-m<86400 and (not os.path.exists(f) or os.path.getmtime(f)<m): print(k)
")
for r in $RUNNING $CHANGED; do
  case $r in c128) B=6000000;; t_xu1|t_xu2|t_ec|t_gt2) B=2000000;; *) B=3500000;; esac
  RUN=$r BUDGET=$B bash $S/panel_refresh2.sh > $S/panel_refresh2_$r.log 2>&1 || echo "refresh of $r failed (see $S/panel_refresh2_$r.log)"
done
# spectra: one background queue (4 threads, CPUs 0-3) adds the new snapshots of the running runs to their caches; started here
# unless the previous one is still busy (its PID in speclog/queue.pid)
qp=$S/speclog/queue.pid
if ! { [ -f $qp ] && kill -0 $(cat $qp) 2>/dev/null; }; then
  nohup bash $S/spec_queue.sh 24-27 0 $RUNNING > /dev/null 2>&1 < /dev/null & echo $! > $qp
fi
# spectra (snapspec.jl, filled by $S/spec_queue.sh): copy the export of every running run that has one
D=/zpool/nvme/data/tabel_scratch/titan_popiii/runs; SPECR=""
for r in $RUNNING; do
  [ -f $D/$r/export/spec.json ] && [ -f $D/$r/export/spec.f32 ] || continue
  cp $D/$r/export/spec.json $S/panel/$r/spec.json.tmp && mv $S/panel/$r/spec.json.tmp $S/panel/$r/spec.json
  base64 -w0 $D/$r/export/spec.f32 > $S/panel/$r/spec.b64.txt.tmp && mv $S/panel/$r/spec.b64.txt.tmp $S/panel/$r/spec.b64.txt; SPECR="$SPECR $r"
done
# optical depth by process (snapspec.jl export, spectau.*): copy where both files exist
TAUR=""
for r in $RUNNING; do
  [ -f $D/$r/export/spectau.json ] && [ -f $D/$r/export/spectau.i16 ] || continue
  cp $D/$r/export/spectau.json $S/panel/$r/spectau.json.tmp && mv $S/panel/$r/spectau.json.tmp $S/panel/$r/spectau.json
  base64 -w0 $D/$r/export/spectau.i16 > $S/panel/$r/spectau.b64.txt.tmp && mv $S/panel/$r/spectau.b64.txt.tmp $S/panel/$r/spectau.b64.txt; TAUR="$TAUR $r"
done
# full spectra with lines (snapfull.jl, filled by $S/specfull_queue.sh): copy the export of every running run that has one
# keep one light queue going (4 threads, CPUs 112-115, SPECMAX=2 new spectra per run and pass, the newest first), unless one runs
fq=$S/speclog/specfull_queue.pid; fqp=$(head -1 $fq 2>/dev/null | awk '{print $2}')
if ! { [ -n "$fqp" ] && kill -0 $fqp 2>/dev/null; }; then
  SPECMAX=2 setsid nohup bash $S/specfull_queue_light.sh $RUNNING > $S/speclog/specfull_queue.out 2>&1 < /dev/null &
  sleep 1; ps -C bash -o pid=,args= | awk -v f=$S/specfull_queue_light.sh '$3==f {print "queue", $1}' > $fq
fi
FULLR=""
for r in $RUNNING; do
  [ -f $D/$r/export/specfull.json ] && [ -f $D/$r/export/specfull.f32 ] || continue
  cp $D/$r/export/specfull.json $S/panel/$r/specfull.json.tmp && mv $S/panel/$r/specfull.json.tmp $S/panel/$r/specfull.json
  base64 -w0 $D/$r/export/specfull.f32 > $S/panel/$r/specfull.b64.txt.tmp && mv $S/panel/$r/specfull.b64.txt.tmp $S/panel/$r/specfull.b64.txt; FULLR="$FULLR $r"
done
# printable ABN02-style figure (scripts/abnfig) of every running run with frames
for r in $RUNNING; do [ -f $S/panel/$r/frames.json ] && timeout 120 $PY -I /home/tabel/Projects/titan.jl/scripts/abnfig/abnfig.py $r --panel $S/panel --out $S/panel/$r/abnfig.pdf > $S/speclog/abnfig_$r.log 2>&1; done
# photons ahead of the photosphere: ionfs_post.jl (frozen tree of ionfs-v2 3499e3a) over the new snapshots of every running run
for r in $RUNNING; do ls $D/$r/snaps_*.jls >/dev/null 2>&1 || continue
  timeout 600 numactl --physcpubind=96-99 --membind=12 julia -t 4 --project=$S/jlenv $S/ionfs_tree/julia/ionfs_post.jl $D/$r > $S/speclog/ionfs_post_$r.log 2>&1
  [ -f $D/$r/ionfs_post.csv ] && cp $D/$r/ionfs_post.csv $S/panel/$r/ionfs_post.csv; done
echo "PUBLISH: runs.json"; for r in $RUNNING; do [ -f $S/panel/$r/ionfs_post.csv ] && echo "PUBLISH: $r/ionfs_post.csv"; done; for r in $RUNNING; do [ -f $S/panel/$r/abnfig.pdf ] && echo "PUBLISH: $r/abnfig.pdf"; done; for r in $RUNNING $CHANGED; do [ -f $S/panel/$r/frames.json ] && echo "PUBLISH: $r/frames.json $r/frames.b64.txt $r/perf.csv"; done
for r in $SPECR; do echo "PUBLISH: $r/spec.json $r/spec.b64.txt"; done
for r in $TAUR; do echo "PUBLISH: $r/spectau.json $r/spectau.b64.txt"; done
for r in $FULLR; do echo "PUBLISH: $r/specfull.json $r/specfull.b64.txt"; done
du -cb $S/panel/runs.json $(for r in $RUNNING $CHANGED; do echo $S/panel/$r/frames.json $S/panel/$r/frames.b64.txt $S/panel/$r/perf.csv; done) $(for r in $SPECR; do echo $S/panel/$r/spec.json $S/panel/$r/spec.b64.txt; done) $(for r in $TAUR; do echo $S/panel/$r/spectau.json $S/panel/$r/spectau.b64.txt; done) $(for r in $FULLR; do echo $S/panel/$r/specfull.json $S/panel/$r/specfull.b64.txt; done) 2>/dev/null | tail -1
