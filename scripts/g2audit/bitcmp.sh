#!/bin/bash
# bitcmp.sh <ckA_prefix> <ckB_prefix>: byte comparison of the decompressed TITAN state and of the .aux
A=$1; B=$2
if cmp -s <(gzip -dc $A.tstate.gz) <(gzip -dc $B.tstate.gz); then echo "tstate: identical ($(gzip -dc $A.tstate.gz | wc -c) bytes)"; else echo "tstate: DIFFER"; fi
if cmp -s $A.aux $B.aux; then echo "aux: identical ($(wc -c < $A.aux) bytes)"; else echo "aux: DIFFER"; fi
echo "arrays in state: $(gzip -dc $A.tstate.gz | head -c 16 | tail -c 8 | od -An -t d8)"
