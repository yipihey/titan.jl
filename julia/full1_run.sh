#!/bin/zsh
# FULL-1: physical ABN02 z = 19 chain with feedback (collapse -> formation -> coupled -> ms). Usage: full1_run.sh [stage ...]   (default: all four; resumable)
# E2E_REUSE = work dir of the no-feedback chain (E2E-A) whose collapse/formation snapshots are reused (verified on 240 formation steps, see notes FULL-1).
# Two Julia processes of E2E_THREADS [4] threads in the coupled stage, one in the others.
cd "$(dirname "$0")"
export E2E_FB=1 E2E_PERF=full1 E2E_THREADS=${E2E_THREADS:-4}
export E2E_REUSE=${E2E_REUSE:-/Users/tabel/Projects/titan.jl/.claude/worktrees/agent-a09b72bed2ba32d23/lib/popiii/e2e_work/abn02z19_nofb}
exec julia --project=env e2e_chain.jl "$@"
