#!/bin/bash
# Branch/worktree cleanup after consolidation into the prototype branch.
# For every LOCAL branch fully merged into <target>: tag archive/<branch> at its head, remove its worktree (if any), delete the
# branch. Never touches: remote branches, the protected list, branches whose worktree has a running process (cwd inside it),
# locked worktrees, worktrees with uncommitted changes, and the target itself.
# Usage: cleanup_branches.sh <target> [--do]      (default: dry run, prints the plan)
set -u
TARGET=$1; DO=${2:-}
REPO=$(git rev-parse --show-toplevel 2>/dev/null) || { echo "run inside the repository"; exit 1; }
PROTECT=" $TARGET popiii-subconv subconv-mg tdc-henyey proto-gate iedef logmass mgprod mgsmooth prodpkg pkgstart gasheat popiii-combined popiii-s6 ${PROTECT_EXTRA:-} "
# cwd of every process -> set of busy worktree paths
busy=$(for p in /proc/[0-9]*; do readlink $p/cwd 2>/dev/null; done | sort -u)
declare -A WT LOCK
while read -r line; do
  case $line in
    worktree\ *) cur=${line#worktree };;
    branch\ *) b=${line#branch refs/heads/}; WT[$b]=$cur;;
    locked*) LOCK[$cur]=1;;
  esac
done < <(git -C "$REPO" worktree list --porcelain)
n_del=0; n_keep=0
for b in $(git -C "$REPO" branch --merged "$TARGET" --format='%(refname:short)'); do
  reason=""
  [[ "$PROTECT" == *" $b "* ]] && reason="protected"
  w=${WT[$b]:-}
  if [ -z "$reason" ] && [ -n "$w" ]; then
    [ -n "${LOCK[$w]:-}" ] && reason="worktree locked"
    echo "$busy" | grep -q "^$w\(/\|$\)" && reason="process running in worktree"
    [ -z "$reason" ] && [ -n "$(git -C "$w" status --porcelain 2>/dev/null | grep -v '^??')" ] && reason="uncommitted changes"
  fi
  if [ -n "$reason" ]; then echo "KEEP   $b ($reason)"; n_keep=$((n_keep+1)); continue; fi
  echo "DELETE $b ${w:+(worktree $w)} -> tag archive/$b $(git -C "$REPO" rev-parse --short "$b")"
  n_del=$((n_del+1))
  if [ "$DO" = "--do" ]; then
    git -C "$REPO" tag -f "archive/$b" "$b" >/dev/null
    [ -n "$w" ] && git -C "$REPO" worktree remove "$w"
    git -C "$REPO" branch -d "$b" >/dev/null
  fi
done
echo "merged into $TARGET: $n_del to delete, $n_keep kept${DO:+ (done)}"
echo "not merged into $TARGET (left alone): $(git -C "$REPO" branch --no-merged "$TARGET" | wc -l)"
