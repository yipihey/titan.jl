#!/bin/bash
# Copy of panel_refresh.sh for the auto-refreshing explorer: frame budget from BUDGET (bytes of frames.f32), julia pinned, no spectra.
S=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad; RUN=${RUN:?}; R=/zpool/nvme/data/tabel_scratch/titan_popiii/runs/$RUN; mkdir -p $S/panel/$RUN
cd ~/Projects/titan.jl/julia && numactl --physcpubind=${CPUS:-62-63} julia --startup-file=no --project=$S/jlenv runrec.jl $R budget=${BUDGET:-3500000} skip=logkap,logkapP,eq,logLam | tail -1
cp $R/export/frames.json $S/panel/$RUN/frames.json.tmp && mv $S/panel/$RUN/frames.json.tmp $S/panel/$RUN/frames.json
base64 -w0 $R/export/frames.f32 > $S/panel/$RUN/frames.b64.txt.tmp && mv $S/panel/$RUN/frames.b64.txt.tmp $S/panel/$RUN/frames.b64.txt
cp $R/perf.csv $S/panel/$RUN/perf.csv
