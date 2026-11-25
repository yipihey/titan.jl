#!/bin/bash
# mkimage.sh <worktree> [grey|mg|name]   (name: print the image base name (<branch>-<sha>[-dirty-<h>]) only)   (generalises prod_env/build_prod.sh; docs/POPIII_DEVTOOLS.md)
# Builds the TitanStack image (docs/POPIII_PRECOMPILE.md) of <worktree> into
#   $DATA/images/<branch>-<sha>[-dirty-<hash of git diff + untracked files of julia/>][-mg]/{env,depot}
# and prints, on stdout, the three exports a run needs (all other output goes to stderr, so `eval "$(mkimage.sh WT)"` works):
#   export JLENV=...  export JULIA_DEPOT_PATH=...  export TITAN_PKG=1
# Variants: grey (default, MGROUPS unset) and mg (MGROUPS=4: the image freezes the conditional JuliaSolver include, so MG runs need their own).
# Reuse: an existing image with the same name is a no-op (checked with Base.isprecompiled, ~5 s). A dirty julia/ gets the -dirty-<hash> name.
# The image is built from the same words a gate run uses: the Stage-1 defaults of popiii_stack.jl + the tree's scripts/proto/run_zie1.env
# (zie3 differs only in ENVEUL/ENVCORE, which are run-time reads; a run that falls back prints "falling back to the include path").
# A tree that is not a git checkout (git-archive scratch tree) carries its sha in <tree>/.titan_sha and is named tree-<sha>.
# Env: DATA, MKIMAGE_NODE (NUMA node of the build, default 3), MKIMAGE_THREADS (8), BASEENV (env dir whose Project.toml/Manifest.toml are copied).
set -u
WT=$(cd "${1:?usage: mkimage.sh <worktree> [grey|mg]}" && pwd); VAR=${2:-grey}
DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
BASEENV=${BASEENV:-$DATA/prod_env/env}          # Project.toml + Manifest.toml of the production TitanStack env (never modified; copied)
NODE=${MKIMAGE_NODE:-3}; TH=${MKIMAGE_THREADS:-8}
case $VAR in grey|name) SUF=;; mg) SUF=-mg;; *) echo "variant must be grey or mg" >&2; exit 2;; esac
if git -C $WT rev-parse --git-dir >/dev/null 2>&1 && [ -e $WT/.git ]; then
  SHA=$(git -C $WT rev-parse --short HEAD); BR=$(git -C $WT branch --show-current); BR=${BR:-detached}
  if [ -n "$(git -C $WT status --porcelain -- julia | head -1)" ]; then
    H=$( { git -C $WT diff HEAD -- julia; git -C $WT ls-files -o --exclude-standard julia | sort | while read -r f; do echo "$f"; cat "$WT/$f"; done; } | sha1sum | cut -c1-8)
    DIRTY=-dirty-$H
  else DIRTY=; fi
else SHA=$(cat $WT/.titan_sha); BR=tree; DIRTY=; fi
NAME=${BR//\//_}-$SHA$DIRTY$SUF; IMG=$DATA/images/$NAME
[ $VAR = name ] && { echo $NAME; exit 0; }
export_lines() {
  echo "export JLENV=$IMG/env"; echo "export JULIA_DEPOT_PATH=$IMG/depot:/home/tabel/.julia"; echo "export TITAN_PKG=1"
}
mkdir -p $IMG/env $IMG/depot
exec 9>$IMG/.lock; flock 9                       # one build per image; a second caller waits and then reuses
T0=$(date +%s)
check() {   # the image's cache is valid for the tree as it is now
  [ -f $IMG/built.ok ] || return 1
  JULIA_DEPOT_PATH=$IMG/depot:/home/tabel/.julia numactl --cpunodebind=$NODE --localalloc \
    julia --project=$IMG/env -e 'id = Base.identify_package("TitanStack"); exit(id !== nothing && Base.isprecompiled(id) ? 0 : 1)' >/dev/null 2>&1
}
if check; then
  echo "mkimage: $NAME reused ($(( $(date +%s) - T0 )) s check; built in $(cat $IMG/built.ok) s)" >&2
  export_lines; exit 0
fi
cp $BASEENV/Project.toml $IMG/env/Project.toml
sed -E "/^\[\[deps\.TitanStack\]\]/,/^$/ s#^path = .*#path = \"$WT/julia/TitanStack\"#" $BASEENV/Manifest.toml > $IMG/env/Manifest.toml
grep -q "path = \"$WT/julia/TitanStack\"" $IMG/env/Manifest.toml || { echo "mkimage: cannot point the Manifest at $WT" >&2; exit 3; }
WORDS=$WT/scripts/proto/run_zie1.env; [ -f $WORDS ] || WORDS=$DATA/cklib/zie1/run.env
W=$(grep -v '^#' $WORDS | tr '\n' ' '); [ $VAR = mg ] && W="$W MGROUPS=4"
[ -f $IMG/built.ok ] && mv $IMG/built.ok $IMG/built.ok.old     # (never rm)
echo "mkimage: building $NAME (node $NODE, $TH threads) from $WT" >&2
( cd $WT/julia && env $W JULIA_DEPOT_PATH=$IMG/depot:/home/tabel/.julia numactl --cpunodebind=$NODE --localalloc \
    julia -t $TH --project=$IMG/env TitanStack/build.jl ) > $IMG/build.log 2>&1; rc=$?
tail -1 $IMG/build.log >&2
[ $rc = 0 ] || { echo "mkimage: build FAILED (rc $rc), see $IMG/build.log" >&2; exit 4; }
echo $(( $(date +%s) - T0 )) > $IMG/built.ok
echo "mkimage: built $NAME in $(cat $IMG/built.ok) s" >&2
export_lines
