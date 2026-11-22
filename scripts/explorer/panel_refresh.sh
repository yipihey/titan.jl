#!/bin/bash
# Re-export the run (RUN, default full2) for the panel (frame stride chosen so the base64 file stays under ~11 MB), copy perf.csv.
# Spectra (snapspec.jl): one background job per run fills <run>/spec and writes export/spec.*; whatever export exists is copied.
S=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad; RUN=${RUN:-full2}; R=/zpool/nvme/data/tabel_scratch/titan_popiii/runs/$RUN; mkdir -p $S/panel/$RUN $S/speclog
n=$(ls $R/snaps_*.jls | xargs -n1 du -b | awk '{s+=$1} END{print s}')
# stride chosen by the exporter so that frames.f32 stays below 10.4 MB (base64: 13.9 MB; the limit per published file is 16 MB)
cd ~/Projects/titan.jl/julia && julia --project=$S/jlenv runrec.jl $R budget=10400000 skip=logkap,logkapP,eq,logLam | tail -1
cp $R/export/frames.json $S/panel/$RUN/; base64 -w0 $R/export/frames.f32 > $S/panel/$RUN/frames.b64.txt; cp $R/perf.csv $S/panel/$RUN/
if [ "${NOSPEC:-0}" != 1 ]; then
  pf=$S/speclog/$RUN.pid
  if ! { [ -f $pf ] && kill -0 $(cat $pf) 2>/dev/null; }; then
    last=$(tail -1 $R/perf.csv | cut -d, -f1); smod=160; while (( last / smod > 300 )); do smod=$((smod*2)); done
    TAG=_spec_$RUN setsid nohup julia -t ${SPECTHREADS:-8} --project=$S/jlenv snapspec.jl $R $smod >> $S/speclog/$RUN.out 2>&1 &
    echo $! > $pf
  fi
fi
[ -f $R/export/spec.json ] && { cp $R/export/spec.json $S/panel/$RUN/; base64 -w0 $R/export/spec.f32 > $S/panel/$RUN/spec.b64.txt; }
ls -la $S/panel/$RUN
