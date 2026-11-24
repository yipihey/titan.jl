#!/bin/bash
# T3 (nightly): from the initial condition to M* >= GATE_MSTOP (34) Msun at 100 zones (case z100, as zie1) and 298 zones (case z298, as zie3), in segments of
# GATE_SEG steps (default 400) that each restart the previous segment's final checkpoint; resumable (rerun continues after the last completed segment,
# procedure.log, as scripts/subconv/m4run.sh). The two cases run on two slots in parallel. Then runrec export of each case and the comparison with the
# production export(s): --ref-export DIR (default zie1 -> $DATA/subconv/m4/exports/B0, zie3 -> .../B0z3), via scripts/subconv/m4extra.jl + m4compare.py.
# Metrics t3.<case>.*: M*, R*, Mdot(1e4 Rsun), u_min, T_c, T_eff, L_ph, s_ps, Mdot_sh, F_c/F_tot ... AT FIXED M* (1,3,10,30 Msun) and their relative
# difference to the reference at the same M*; plus run-level counts (steps, Newton failures, wall). Sourced by gate.sh.

t3_case() {   # case cores
  local cs=$1 cores=$2 R=$OUT/t3/$1 SEG=${GATE_SEG:-400} MSTOP=${GATE_MSTOP:-34} NSEGMAX=${GATE_NSEGMAX:-60}
  mkdir -p $R/ck; local log=$R/procedure.log; lg() { echo "$(date '+%F %T') $*" >> $log; }
  [ -f $R/eqeos_he_xe.jls ] || cp ${EQEOS:-$DATA/runs/zie1/eqeos_he_xe.jls} $R/
  local n=$(grep -c "segment done" $log 2>/dev/null); n=${n:-0}; local ck=""
  [ "$n" -gt 0 ] && ck=$(grep 'segment done' $log | tail -1 | awk '{print $5}')
  if [ -n "$ck" ] && python3 -c "import sys; sys.exit(0 if 'reached' in open('$log').read() else 1)"; then lg "already reached M* >= $MSTOP"; return; fi
  lg "lane start: case $cs from ${ck:-the IC}, code $(git -C $WT rev-parse --short HEAD), config $CFG, cores $cores"
  while [ $n -lt $NSEGMAX ]; do
    n=$((n+1)); local save=${cs}S$n out=$R/seg_${cs}S$n.out
    # run_pop writes run.out per call; keep each segment's own
    run_pop $R $cores $WT $SEG $CFG $cs CONVDIAG=1 ETOTMON=1 ETOTLOG=$R/etot.log NWTRACE=$R/nwt.txt ${ck:+RESTART=$ck} SAVE=$save $GATE_T3_ENV
    mv $R/run.out $out; cat $R/wall.txt >> $R/wall_all.txt
    local end=$(grep -E 'run ended|stopped' $out | tail -1 | cut -c1-200)
    local ms=$(grep -E '^ *hot: t' $out | tail -1 | sed -E 's/.*m\* ([0-9.eE+-]+) Msun.*/\1/')
    if ! echo "$end" | grep -q ':maxsteps'; then lg "segment FAILED $save: $end (M* $ms)"; return 1; fi
    lg "segment done $save : $end (M* $ms)"; ck=$save
    python3 -c "import sys; sys.exit(0 if float('$ms')>=$MSTOP else 1)" 2>/dev/null && { lg "reached M* $ms >= $MSTOP"; break; }
  done
}

tier_t3() {
  local D=$OUT/t3; mkdir -p $D; : > $D/notes.txt; rm -f $D/*.kv; pool_init $D
  local cases=${GATE_T3_CASES:-"z100 z298"} cs
  for cs in $cases; do pool_launch t3_case $cs; done
  pool_wait
  for cs in $cases; do
    local R=$D/$cs
    [ -f $R/snaps_*.jls ] 2>/dev/null || ls $R/snaps_*.jls > /dev/null 2>&1 || { echo "no snapshots for $cs" >> $D/notes.txt; continue; }
    mkdir -p $D/exports/$cs
    ( cd $WT/julia && numactl --physcpubind=${SLOTS[0]} --localalloc julia -t 2 --project=$JLENV $WT/scripts/subconv/m4extra.jl $R $D/exports/$cs ) > $D/export_$cs.log 2>&1
    local ref=${REFEXPORT[$cs]:-$DATA/subconv/m4/exports/$([ $cs = z100 ] && echo B0 || echo B0z3)} smin=$([ $cs = z100 ] && echo 290 || echo 775)
    python3 $PROTO/t3_metrics.py $D/exports/$cs t3.$cs $ref $smin > $D/$cs.kv 2> $D/$cs.metrics.err
    ls -v $R/seg_*.out | xargs cat > $R/run.out; awk '{s+=$1} END{print s}' $R/wall_all.txt > $R/wall.txt
    python3 $PROTO/parse_run.py $R t3.$cs.run | sed 's/^t3\.\([^ ]*\)\.run\./t3.\1.run_/' >> $D/$cs.kv
    python3 $PROTO/ledger.py $R/etot.log t3.$cs.etot >> $D/$cs.kv
    # the human-readable m4compare tables: new case against the reference export (copied in as case REF), first case = reference
    mkdir -p $D/exports_cmp; ln -sfn $D/exports/$cs $D/exports_cmp/$cs; ln -sfn $ref $D/exports_cmp/REF
    python3 $WT/scripts/subconv/m4compare.py REF $cs --exports $D/exports_cmp --out $D --tag $cs --ref REF --smin $smin > $D/m4compare_$cs.log 2>&1 \
      || echo "m4compare failed for $cs (see t3/m4compare_$cs.log)" >> $D/notes.txt
  done
  python3 $PROTO/finalize.py $D T3
}
