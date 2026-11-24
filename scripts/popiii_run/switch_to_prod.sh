#!/bin/bash
# Switch a running Pop III production run (driven by iceulrun.sh) to the prodpkg tree + precompiled image, at a segment boundary.
# Usage: switch_to_prod.sh <run> "<exact argv of the old driver loop>"        (dry run: lists, changes nothing)
#        GO=1 switch_to_prod.sh <run> "<argv>"                                (stops the old driver loop only, waits for the segment, relaunches)
# The old driver's argv is read from `ps -eo pid=,args=`; it is matched by EQUALITY of the whole string, never by pattern. The segment's
# julia process is the (only) child of the old driver; it is not signalled, it finishes its segment (setsid'd, output file already open).
# The new driver is launched with TITAN_PKG=1 and the prod depot in its environment (iceulrun.sh passes the environment on to julia through
# `env`; nothing is written to the run directory's launch.env / run.env). Its first segment restarts the finished segment's SAVE.
set -u
RUN=$1; OLDARGV=$2; GO=${GO:-0}
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; R=$DATA/runs/$RUN; E=$DATA/prod_env
PROD=/home/tabel/Projects/titan.jl/.claude/worktrees/prod
drv=$(ps -eo pid=,args= | awk -v a="$OLDARGV" -v me1=$$ -v me2=$PPID '{p=$1; sub(/^ *[0-9]+ /,""); if ($0==a && p!=me1 && p!=me2) print p}')
[ "$(echo $drv | wc -w)" = 1 ] || { echo "expected exactly one driver with argv '$OLDARGV', got: '$drv'"; exit 1; }
seg=$(ps -o pid= --ppid $drv | tr -d ' ')
[ "$(echo $seg | wc -w)" = 1 ] || { echo "expected exactly one segment child of driver $drv, got: '$seg'"; exit 1; }
SAVE=$(tr '\0' '\n' < /proc/$seg/environ | sed -n 's/^SAVE=//p')
echo "run $RUN: old driver pid $drv, segment pid $seg (SAVE=$SAVE, $(tr '\0' ' ' < /proc/$seg/cmdline | cut -c1-120))"
echo "plan: kill -TERM $drv; wait for $seg to exit and ck/m2c_$SAVE.aux; launch $PROD/scripts/popiii_run/iceulrun.sh $RUN $SAVE with JLENV=$E/env JULIA_DEPOT_PATH=$E/depot:\$HOME/.julia TITAN_PKG=1"
[ "$GO" = 1 ] || { echo "(dry run)"; exit 0; }
kill -TERM $drv
while kill -0 $seg 2>/dev/null; do sleep 20; done
grep -E 'run ended|stopped' $R/seg_$SAVE.out | tail -1 | grep -q ':maxsteps' || { echo "segment $SAVE did not end with :maxsteps; not relaunching"; exit 1; }
[ -f $R/ck/m2c_$SAVE.aux ] || { echo "no checkpoint $SAVE; not relaunching"; exit 1; }
echo "$(date '+%F %T') $RUN: switched to prodpkg (TitanStack image, TITAN_PKG=1, MGROUPS unset): old driver $drv stopped, segment $SAVE finished, relaunch from $SAVE" >> $R/procedure.log
cd $PROD
JLENV=$E/env JULIA_DEPOT_PATH=$E/depot:$HOME/.julia TITAN_PKG=1 setsid nohup bash scripts/popiii_run/iceulrun.sh $RUN $SAVE > /dev/null 2>&1 < /dev/null &
sleep 5; tail -n 2 $R/procedure.log | cut -c1-200
