#!/bin/bash
# run_prof.sh <tree> <script.jl> [args...]: one of the g2 scratch scripts with the G2 table settings, pinned to 104-111
G=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2
TREE=$1; S=$2; shift 2
cd $G
H2PF=${H2PF:-1} SCVHXE=1 EQEOS_CACHE=$G/eqeos_he_xe.jls TREE=$TREE timeout 3000 numactl --physcpubind=104-111 --membind=13 \
    julia -t ${NT:-1} --startup-file=no --project=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/jlenv $S "$@"
