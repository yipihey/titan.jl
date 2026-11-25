#!/bin/bash
# ws.sh <step> : wait (max 590 s) until the A run log reaches <step>, then print the dt trajectory every 5th snapshot
P=/zpool/nvme/data/tabel_scratch/titan_popiii/proto/stage2c
timeout 590 bash -c "until [ \$(grep 'LEDGER k' $P/A/m2c_abn02red_lines_n128_zie1.log | tail -1 | sed 's/.*LEDGER k \([0-9]*\).*/\1/') -ge $1 ]; do sleep 20; done"
cd /home/tabel/Projects/titan.jl/.claude/worktrees/proto
julia --startup-file=no scripts/proto/stage2c_table.jl $P/A /tmp/tA.tsv >/dev/null
cut -f1-3,6-8,11,14 /tmp/tA.tsv | column -t | awk 'NR>1 && NR%5==0' | tail -${2:-6}
grep 'LEDGER k' $P/A/m2c_abn02red_lines_n128_zie1.log | tail -1 | cut -c1-30; date
