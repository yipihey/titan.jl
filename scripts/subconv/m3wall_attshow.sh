#!/bin/bash
# att_show.sh file attno : condensed trace
awk "/NWT head att $2 /,/NWT end att $2 /" $1 | grep -E "^NKR|^NWT" | sed -E 's/rel residual ([0-9.e+-]{6})[0-9]*/rel \1/; s/\|delta0\|_W ([0-9.]{5})[0-9]*/d0 \1/; s/\|delta\|_W ([0-9.]{5})[0-9]*/d \1/; s/krylov its/k/; s/drift 0.0//; s/ base .* \| / | /; s/converged //' | cut -c1-${3:-260}
