#!/bin/bash
# T1 (~30 min): per library state NKCHECK (grey and MGROUPS=4), ETOTMON ledger (40 steps), checkpoint round trip; once: bit-identity of the IC start
# against --bitref. Jobs run on the core slots of GATE_CORES (at most 4 at once). Sourced by gate.sh.

t1_nk() {   # state mode cores    mode = grey | mg
  local line=$1 mode=$2 cores=$3; local n=$(st_field "$line" 1) cs=$(st_field "$line" 5) nkz=$(st_field "$line" 6) dt=$(st_field "$line" 7)
  local ckdir=$(st_field "$line" 2) ckn=$(st_field "$line" 3) mgck=$(st_field "$line" 8)
  local R=$OUT/t1/nk_${n}_$mode extra="" it
  rm -rf $R $R.kv
  if [ $mode = mg ]; then
    extra="MGROUPS=4"
    # A grey checkpoint has no group fields and no equilibrium partition (mgprod docs 24.4: the first step needs ~20 failed attempts, NKCHECK is never reached
    # from it). The MG probe therefore runs from a checkpoint that carries the group fields (optional 8th column of states.txt), or, if none is listed, from the
    # grey state advanced GATE_MGWARM (default 20) steps with MGROUPS=4 (SAVE=mgw).
    if [ -n "$mgck" ] && [ "$mgck" != "-" ]; then ckn=$mgck
    else
      local W=$OUT/t1/mgwarm_$n; rm -rf $W; stage_ck $ckdir $ckn $W || { kv $R.kv t1.nk.$n.mg.stage_failed 1; return; }
      run_pop $W $cores $WT ${GATE_MGWARM:-20} $CFG $cs RESTART=$ckn MGROUPS=4 SAVE=mgw
      if [ ! -f $W/ck/m2c_mgw.aux ]; then kv $R.kv t1.nk.$n.mg.warm_failed 1; echo "nk $n mg: the MGROUPS=4 warm-up from $ckn wrote no checkpoint (t1/mgwarm_$n/run.out)" >> $OUT/t1/notes.txt; return; fi
      ckdir=$W/ck; ckn=mgw
    fi
  fi
  for it in 2 1; do                      # NKCHECK=2 as in the docs; a state that converges before iteration 2 only gets iteration 1 (recorded)
    rm -rf $R; stage_ck $ckdir $ckn $R || { kv $R.kv t1.nk.$n.$mode.stage_failed 1; return; }
    run_pop $R $cores $WT 3 $CFG $cs RESTART=$ckn NWPROBE=$dt NWPROBEAT=0 NKRYLOV=1 NKCHECK=$it NKZ=$nkz $extra
    grep -aq "^NKCHECK eq" $R/run.out && break
    cp $R/run.out $R.attempt_it$it.out 2>/dev/null      # keep the log of an attempt that gave no listing (the retry removes the run dir)
  done
  python3 $PROTO/nkparse.py $R/run.out nk.$n.$mode $mode > $R.kv; kv $R.kv nk.$n.$mode.iteration $it
  [ -f $R/run.out ] && rm -rf $R/ck
}

t1_etot() {  # state cores
  local line=$1 cores=$2; local n=$(st_field "$line" 1) cs=$(st_field "$line" 5); local R=$OUT/t1/etot_$n
  rm -rf $R; stage_ck $(st_field "$line" 2) $(st_field "$line" 3) $R || return
  run_pop $R $cores $WT ${GATE_ETOT_STEPS:-40} $CFG $cs RESTART=$(st_field "$line" 3) ETOTMON=1 ETOTLOG=$R/etot.log
  python3 $PROTO/ledger.py $R/etot.log etot.$n > $R.kv; python3 $PROTO/parse_run.py $R t1.etot_$n >> $R.kv
  rm -rf $R/ck
}

t1_roundtrip() {  # state cores: continue 15 steps with a checkpoint written at k0+5, restart from it, 10 steps: snapshots of steps > k0+5 bit-identical
  local line=$1 cores=$2; local n=$(st_field "$line" 1) cs=$(st_field "$line" 5) ckn=$(st_field "$line" 3); local k0=$(ck_step $ckn) k5
  k5=$((k0+5)); local A=$OUT/t1/rt_${n}_A B=$OUT/t1/rt_${n}_B
  rm -rf $A $B; stage_ck $(st_field "$line" 2) $ckn $A || return
  run_pop $A $cores $WT 15 $CFG $cs RESTART=$ckn SAVE=rtA CKDUMP=$k5:$k5
  if [ ! -f $A/ck/m2c_rtA_k$k5.aux ]; then kv $A.kv t1.roundtrip.$n.failed 1; kv $A.kv t1.roundtrip.$n.no_checkpoint 1; return; fi
  mkdir -p $B/ck; cp $A/ck/m2c_rtA_k$k5.* $B/ck/
  run_pop $B $cores $WT 10 $CFG $cs RESTART=rtA_k$k5
  ( cd $WT/julia && julia --project=$JLENV $PROTO/cmp_snaps.jl $A $B $k5 t1.roundtrip.$n ) > $A.kv 2> $A.cmp.err
  local nd=$(awk '/fields_diff/{print $2}' $A.kv) nc=$(awk '/snaps_compared/{print $2}' $A.kv)
  kv $A.kv t1.roundtrip.$n.failed $([ "${nd:-1}" = 0 ] && [ "${nc:-0}" -gt 0 ] && echo 0 || echo 1)
  rm -rf $A/ck $B/ck
}

t1_bitid() {  # code-tree-of-reference cores: IC start, MGROUPS unset, 40 steps, both trees on the same case; compares every field of every snapshot
  local ref=$1 cores=$2 cs=${GATE_BITCASE:-z100} n=${GATE_BITSTEPS:-40}
  local A=$OUT/t1/bitid_new B=$OUT/t1/bitid_ref
  rm -rf $A $B
  run_pop $A $cores $WT $n $CFG $cs &
  local pa=$!; local c2=${SLOTS[1]:-$cores}
  run_pop $B $c2 $ref $n $CFG $cs; wait $pa
  ( cd $WT/julia && julia --project=$JLENV $PROTO/cmp_snaps.jl $A $B 0 t1.bitid ) > $OUT/t1/bitid.kv 2> $OUT/t1/bitid.err
  local nd=$(awk '/fields_diff/{print $2}' $OUT/t1/bitid.kv) nc=$(awk '/snaps_compared/{print $2}' $OUT/t1/bitid.kv)
  kv $OUT/t1/bitid.kv t1.bitid.failed $([ "${nd:-1}" = 0 ] && [ "${nc:-0}" -gt 0 ] && echo 0 || echo 1)
  rm -rf $A/ck $B/ck
}

tier_t1() {
  local D=$OUT/t1; rm -rf $D; mkdir -p $D; : > $D/notes.txt; pool_init $D
  local line
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    case " ${GATE_T1_JOBS:-nk etot rt} " in *" nk "*) pool_launch t1_nk "$line" grey; pool_launch t1_nk "$line" mg;; esac
    case " ${GATE_T1_JOBS:-nk etot rt} " in *" etot "*) pool_launch t1_etot "$line";; esac
    case " ${GATE_T1_JOBS:-nk etot rt} " in *" rt "*) pool_launch t1_roundtrip "$line";; esac
  done < <(state_lines)
  if [ -n "${BITREF:-}" ]; then
    # the reference tree: a directory (worktree) or a commit/branch, checked out detached under $OUT (removed afterwards)
    local refdir=$BITREF made=0
    if [ ! -d "$BITREF/julia" ]; then refdir=$OUT/bitref_wt; git -C $WT worktree add --detach $refdir $BITREF > $D/bitref.log 2>&1 && made=1 || { echo "SKIP bitid: cannot check out $BITREF" >> $D/notes.txt; refdir=; }; fi
    if [ -n "$refdir" ]; then
      pool_wait; echo "bit-identity reference: $BITREF ($(git -C $refdir rev-parse --short HEAD 2>/dev/null)) vs $(git -C $WT rev-parse --short HEAD), $GATE_CORES" >> $D/notes.txt
      t1_bitid $refdir ${SLOTS[0]}
      [ $made = 1 ] && git -C $WT worktree remove --force $refdir
    fi
  else echo "SKIP bitid: no --bitref given" >> $D/notes.txt; fi
  pool_wait
  python3 $PROTO/finalize.py $D T1
}
