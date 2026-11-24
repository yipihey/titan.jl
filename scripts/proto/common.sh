#!/bin/bash
# Shared helpers of the prototype gate (sourced by gate.sh, make_library.sh). See docs/POPIII_PROTO_GATE.md.
DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
PROTO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WT=${CODE:-$(cd "$PROTO/../.." && pwd)}            # the tree under test (julia/ is run from here)
JLENV=${JLENV:-$DATA/subconv/jlenv}
GATE_CORES=${GATE_CORES:-"120-127"}               # core slots, 8 cores each; T1/T2 use at most 4 at once (fewer slots = more serial)
GATE_THREADS=${GATE_THREADS:-8}
PY=${PY:-python3}
export TITAN_NOLIB=${TITAN_NOLIB:-1}

die() { echo "gate: $*" >&2; exit 2; }
stamp() { date '+%F %T'; }

# gate_config <file> : the KEY=VALUE words of a config file ('#' comments, '@include <file>' relative to the file)
gate_config() {
  local f=$1 d; d=$(dirname "$f"); [ -f "$f" ] || die "config $f not found"
  while IFS= read -r l; do
    l=${l%%#*}
    case "$l" in "@include "*) gate_config "$d/${l#@include }"; continue;; esac
    for w in $l; do echo "$w"; done
  done < "$f"
}
# case_words <case> : mesh/IC words of cases/<case>.env
case_words() { gate_config "$PROTO/cases/$1.env"; }

# states: one line per state of states.txt (filtered by GATE_STATES=a,b,c)
state_lines() {
  grep -v '^\s*#' "$PROTO/states.txt" | awk 'NF>=7' | while read -r l; do
    n=${l%% *}
    if [ -n "${GATE_STATES:-}" ]; then case ",$GATE_STATES," in *",$n,"*) ;; *) continue;; esac; fi
    echo "$l"
  done
}
# st_field <line> <1-based column>
st_field() { echo "$1" | awk -v c=$2 '{print $c}'; }
ck_step() { echo "$1" | sed -E 's/.*_k([0-9]+)$/\1/'; }   # step of a checkpoint name (…_kNNN)

# stage_ck <ckdir> <ckname> <rundir>: copy the checkpoint files into the run's own ck dir (never written in place)
stage_ck() { mkdir -p "$3/ck"; for e in aux tstate.gz; do [ -f "$3/ck/m2c_$2.$e" ] || cp "$1/m2c_$2.$e" "$3/ck/" || return 1; done; }

# run_pop <rundir> <cores> <code-tree> <maxsteps> <config-file> <case> [extra env words...]
# One popiii_m2c.jl process (as m3brun.sh / mgint_run.sh): the config's words, then the case's, then the extras (later wins), SNAPEVERY=5.
# Writes <rundir>/run.out and <rundir>/wall.txt (wall seconds of the process incl. start-up). Returns julia's exit status.
run_pop() {
  local R=$1 cores=$2 code=$3 max=$4 cfg=$5 cs=$6; shift 6
  mkdir -p "$R/ck"; [ -f "$R/eqeos_he_xe.jls" ] || cp "${EQEOS:-$DATA/runs/zie1/eqeos_he_xe.jls}" "$R/"
  local words=(); mapfile -t words < <(gate_config "$cfg"; case_words "$cs")
  local s=$(date +%s.%N) rc
  ( cd "$code/julia" && env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$(basename $R) \
      "${words[@]}" "$@" MAXSTEPS=$max SNAPEVERY=5 \
      ${GATE_TIMEOUT:+timeout $GATE_TIMEOUT} numactl --physcpubind=$cores --localalloc julia -t $GATE_THREADS --project=$JLENV popiii_m2c.jl lines 128 \
      > $R/run.out 2>&1 < /dev/null )
  rc=$?
  echo "$(echo "$(date +%s.%N) - $s" | bc)" > $R/wall.txt
  return $rc
}

# ---- slot pool: pool_launch <function> <args...> runs "<function> <args...> <cores>" in the background on a free slot (at most 4 slots)
POOLD=
pool_init() { POOLD=$1/.slots; rm -rf "$POOLD"; mkdir -p "$POOLD"; read -ra SLOTS <<< "$GATE_CORES"; NS=${#SLOTS[@]}; [ $NS -gt 4 ] && NS=4; [ $NS -ge 1 ] || die "no cores"; }
pool_launch() {
  while true; do
    for ((i=0;i<NS;i++)); do
      if mkdir "$POOLD/slot$i" 2>/dev/null; then ( "$@" "${SLOTS[$i]}"; rmdir "$POOLD/slot$i" ) & return; fi
    done; sleep 3
  done
}
pool_wait() { wait; }

# kv <file> <key> <value> : append a metric (flat key -> number); non-numeric values are rejected by finalize.py
kv() { echo "$2 $3" >> "$1"; }
