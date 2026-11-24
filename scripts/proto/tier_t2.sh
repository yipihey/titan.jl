#!/bin/bash
# T2 (1-2 h): GATE_PILOT_STEPS (default 200) steps from every library state, one 8-core slot each, at most 4 at once. Sourced by gate.sh.
# Metrics per state <n> (keys pilot.<n>.*): Newton-failure / error rejections, Newton iteration statistics, dt min/median, simulated time, wall s/step,
# F_c/F_tot and P_t/P (CONVDIAG snapshots; P_t only with TDC), R*, Mdot, the ETOTMON ledger (|sum imb|/sum thr, max |G|).

t2_pilot() {  # state cores
  local line=$1 cores=$2; local n=$(st_field "$line" 1) cs=$(st_field "$line" 5); local R=$OUT/t2/pilot_$n
  rm -rf $R; stage_ck $(st_field "$line" 2) $(st_field "$line" 3) $R || { kv $R.kv pilot.$n.stage_failed 1; return; }
  run_pop $R $cores $WT ${GATE_PILOT_STEPS:-200} $CFG $cs RESTART=$(st_field "$line" 3) CONVDIAG=1 ETOTMON=1 ETOTLOG=$R/etot.log NWTRACE=$R/nwt.txt $GATE_PILOT_ENV
  python3 $PROTO/parse_run.py $R pilot.$n > $R.kv
  python3 $PROTO/ledger.py $R/etot.log pilot.$n.etot >> $R.kv
  julia --project=$JLENV $PROTO/pilot_snap.jl $R pilot.$n >> $R.kv 2> $R.snap.err
  rm -rf $R/ck
}

tier_t2() {
  local D=$OUT/t2; rm -rf $D; mkdir -p $D; : > $D/notes.txt; pool_init $D
  local line
  while IFS= read -r line; do [ -n "$line" ] && pool_launch t2_pilot "$line"; done < <(state_lines)
  pool_wait
  python3 $PROTO/finalize.py $D T2
}
