#!/bin/bash
#---------------------------------------------------------------------------
#  scaling.sh : measure SMP (OpenMP) parallel scaling and append one row per
#               run to a TSV file.
#
#  Usage
#    tools/scaling/scaling.sh [options]
#
#  Typical invocations
#    # First, a smoke test (small cases only, 1 and max threads, once each)
#    tools/scaling/scaling.sh --tier quick --repeats 1 --threads "1 112"
#
#    # Production (assuming a 112-core machine)
#    tools/scaling/scaling.sh --threads "1 2 4 8 14 28 56 84 112" \
#                             --backends "eigen amgcl" --repeats 3
#
#    # If it stops part-way through, simply run the same command again.
#    # Any (case, backend, threads, repeat) already in the TSV is skipped.
#
#  Options
#    --exe PATH        yuzen executable          (default build/bin/yuzen)
#    --cases "a b"     case names (case column of manifest.tsv); default from --tier
#    --tier quick|full|all                        (default all)
#    --threads "1 2 4" list of thread counts      (default 1 2 4 8 16 32 64)
#    --backends "..."  list of auto|eigen|amgcl (default eigen amgcl)
#    --repeats N       times to measure each condition (default 3)
#    --out FILE        TSV to append to            (default tools/scaling/results/scaling-<host>.tsv)
#    --workdir DIR     where to run the cases      (default runs/scaling)
#    --bind close|spread|none  OMP_PROC_BIND      (default close)
#    --timeout SEC     limit for one run (0 = none) (default 0)
#    --keep            do not delete the per-run output directory
#    --dry-run         only print what would be run
#
#  For the result format see tools/scaling/README.md. To summarize:
#    python3 tools/scaling/summarize.py <TSV>
#
#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0
#---------------------------------------------------------------------------
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
MANIFEST="$ROOT/cases/scaling/manifest.tsv"
CASEDIR="$ROOT/cases/scaling"

EXE="$ROOT/build/bin/yuzen"
CASES=""
TIER="all"
THREADS="1 2 4 8 16 32 64"
BACKENDS="eigen amgcl"
REPEATS=3
OUT=""
WORKDIR="$ROOT/runs/scaling"
BIND="close"
TIMEOUT=0
KEEP=0
DRY=0

die() { echo "scaling.sh: $*" >&2; exit 2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --exe)      EXE=$2; shift 2;;
    --cases)    CASES=$2; shift 2;;
    --tier)     TIER=$2; shift 2;;
    --threads)  THREADS=$2; shift 2;;
    --backends) BACKENDS=$2; shift 2;;
    --repeats)  REPEATS=$2; shift 2;;
    --out)      OUT=$2; shift 2;;
    --workdir)  WORKDIR=$2; shift 2;;
    --bind)     BIND=$2; shift 2;;
    --timeout)  TIMEOUT=$2; shift 2;;
    --keep)     KEEP=1; shift;;
    --dry-run)  DRY=1; shift;;
    -h|--help)  sed -n '2,50p' "$0"; exit 0;;
    *) die "unknown option: $1";;
  esac
done

[ -x "$EXE" ] || die "executable not found: $EXE  (pass --exe, or build it first)"
#  Each run cds into its working directory, so a relative path would be lost.
case "$EXE" in /*) ;; *) EXE=$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE");; esac
[ -f "$MANIFEST" ] || die "manifest not found: $MANIFEST  (run python3 tools/scaling/make_cases.py first)"

HOST=$(hostname -s 2>/dev/null || echo unknown)
[ -n "$OUT" ] || OUT="$ROOT/tools/scaling/results/scaling-$HOST.tsv"
mkdir -p "$(dirname "$OUT")" "$WORKDIR"

#---------------------------------------------------------------------------
#  Decide which cases to run
#---------------------------------------------------------------------------
if [ -z "$CASES" ]; then
  if [ "$TIER" = "all" ]; then
    CASES=$(awk -F'\t' 'NR>1{print $1}' "$MANIFEST")
  else
    CASES=$(awk -F'\t' -v t="$TIER" 'NR>1 && $8==t{print $1}' "$MANIFEST")
  fi
fi
[ -n "$CASES" ] || die "no cases selected (tier=$TIER)"

meta() {  # meta <case> <column-number>
  awk -F'\t' -v c="$1" -v n="$2" 'NR>1 && $1==c{print $n}' "$MANIFEST"
}

#---------------------------------------------------------------------------
#  Record the environment. Written once, as # lines at the top of the TSV.
#  This is what later tells you which machine and which build the numbers
#  came from.
#---------------------------------------------------------------------------
cpu=$(awk -F: '/^model name/{gsub(/^ +/,"",$2); print $2; exit}' /proc/cpuinfo 2>/dev/null)
[ -n "${cpu:-}" ] || cpu=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)
sockets=$(lscpu 2>/dev/null | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2}')
corespersocket=$(lscpu 2>/dev/null | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2); print $2}')
threadspercore=$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2}')
numanodes=$(lscpu 2>/dev/null | awk -F: '/^NUMA node\(s\)/{gsub(/ /,"",$2); print $2}')
memkb=$(awk '/^MemTotal/{print $2}' /proc/meminfo 2>/dev/null)
#  The version, as --version prints it (Yuzen-1.0.0).
ver=$("$EXE" --version 2>/dev/null | awk '/Yuzen-[0-9]/{print $1; exit}')
[ -n "${ver:-}" ] || ver=unknown

if [ ! -s "$OUT" ]; then
  {
    echo "# Yuzen SMP scaling measurement"
    echo "# format-version	1"
    echo "# generated-by	tools/scaling/scaling.sh"
    echo "# date	$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# host	$HOST"
    echo "# cpu	${cpu:-unknown}"
    echo "# sockets	${sockets:-unknown}"
    echo "# cores-per-socket	${corespersocket:-unknown}"
    echo "# threads-per-core	${threadspercore:-unknown}"
    echo "# numa-nodes	${numanodes:-unknown}"
    echo "# memory-kb	${memkb:-unknown}"
    echo "# os	$(uname -sr)"
    echo "# yuzen-version	${ver:-unknown}"
    #  Write the executable as a path relative to the repository. The TSV is
    #  meant to be handed to other people, so do not leave an absolute path on
    #  a shared machine in it (these can contain project codes or user IDs).
    #  If it lies outside the repository, keep the file name only.
    exe_rel="${EXE#"$ROOT"/}"
    case "$exe_rel" in
        /*) exe_rel="$(basename "$EXE")" ;;
    esac
    echo "# executable	$exe_rel"
    echo "# omp-proc-bind	$BIND"
    echo "#"
    echo "# One row per run. wall_s is the time-step loop only (JobInfo.txt"
    echo "# 'Time loop [s]'), which excludes grid construction and start-up."
    printf "case\tdim\tnx\tny\tnodes\tbackend\tthreads\trepeat\tsteps\twall_s\tnk_total\tbcg_total\tstatus\n"
  } > "$OUT"
  echo "scaling.sh: created $OUT"
else
  echo "scaling.sh: appending to $OUT (already-measured rows are skipped)"
fi

#  Look up rows that have already been measured (so a run can be resumed)
done_key() {  # done_key <case> <backend> <threads> <repeat>
  awk -F'\t' -v a="$1" -v b="$2" -v c="$3" -v d="$4" \
    '!/^#/ && $1==a && $6==b && $7==c && $8==d {found=1} END{exit !found}' "$OUT"
}

#---------------------------------------------------------------------------
#  1 run
#---------------------------------------------------------------------------
run_one() {  # run_one <case> <backend> <threads> <repeat>
  local c=$1 be=$2 th=$3 rep=$4
  local dim nx ny nodes steps d rc j
  dim=$(meta "$c" 2); nx=$(meta "$c" 3); ny=$(meta "$c" 4)
  nodes=$(meta "$c" 5); steps=$(meta "$c" 7)

  d="$WORKDIR/${c}_${be}_t${th}_r${rep}"
  rm -rf "$d"; mkdir -p "$d" || return 1

  #  The backend is selected by inserting a Backend key into the input.
  #  None of the shipped cases sets Backend, so it is enough to add it at the
  #  start of the LinearSolver block (just before PreconditionType).
  #  It must go into both the height and the concentration equation, or one
  #  side stays on auto.
  if [ "$be" = "auto" ]; then
    cp "$CASEDIR/$c.toml" "$d/case.toml"
  else
    sed "s/^PreconditionType = \"ILU\"/Backend = \"$be\"\nPreconditionType = \"ILU\"/" \
        "$CASEDIR/$c.toml" > "$d/case.toml"
  fi

  (
    cd "$d" || exit 1
    export OMP_NUM_THREADS="$th"
    if [ "$BIND" != "none" ]; then
      export OMP_PROC_BIND="$BIND"
      export OMP_PLACES=cores
    fi
    #  Nested parallelism or dynamic thread adjustment would make the
    #  measurements erratic.
    export OMP_DYNAMIC=false
    export OMP_NESTED=false
    #  Each run uses the thread count it is given: the solver's cap on the
    #  threads of a small grid is off, since measurements like these are what
    #  the cap is chosen from.
    export YUZEN_THREAD_CAP=0
    if [ "$TIMEOUT" -gt 0 ] 2>/dev/null; then
      timeout "$TIMEOUT" "$EXE" case.toml > stdout.txt 2>&1
    else
      "$EXE" case.toml > stdout.txt 2>&1
    fi
  )
  rc=$?

  j=$(ls -d "$d"/*/ 2>/dev/null | head -1)
  local wall="" nsteps="" nk="" bcg="" status="ok" real_be=""
  if [ -n "$j" ] && [ -f "$j/JobInfo.txt" ]; then
    wall=$(awk -F'\t'   '/^Time loop/         {print $2}' "$j/JobInfo.txt")
    nsteps=$(awk -F'\t' '/^Time steps/        {print $2}' "$j/JobInfo.txt")
    nk=$(awk -F'\t'     '/^Newton iterations/ {print $2}' "$j/JobInfo.txt")
    bcg=$(awk -F'\t'    '/^Linear iterations/ {print $2}' "$j/JobInfo.txt")
    real_be=$(awk -F'\t' '/^Linear solver <height>/{gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2}' "$j/JobInfo.txt")
    if [ -s "$j/FatalErrors.txt" ]; then status="fatal"; fi
    if [ -s "$j/Warnings.txt" ] && [ "$status" = "ok" ]; then status="warned"; fi
  fi
  if [ $rc -ne 0 ]; then status="exit$rc"; fi
  if [ -z "$wall" ]; then status="nodata"; wall="NA"; nsteps="NA"; nk="NA"; bcg="NA"; fi
  #  If the run did not take the requested number of steps, it is not comparable.
  if [ "$nsteps" != "NA" ] && [ "$nsteps" != "$steps" ] && [ "$status" = "ok" ]; then
    status="steps$nsteps"
  fi
  #  If something other than auto was asked for but a different backend was
  #  chosen (not built, etc.), record the one that was actually used.
  if [ -n "$real_be" ] && [ "$be" != "auto" ] && [ "$real_be" != "$be" ]; then
    status="${status}/used:$real_be"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$c" "$dim" "$nx" "$ny" "$nodes" "$be" "$th" "$rep" \
    "$nsteps" "$wall" "$nk" "$bcg" "$status" >> "$OUT"

  printf "  %-16s %-7s th=%-4s r=%s  %8s s  %s\n" "$c" "$be" "$th" "$rep" "$wall" "$status"

  [ $KEEP -eq 1 ] || rm -rf "$d"
  return 0
}

#---------------------------------------------------------------------------
#  Main body. The loops run from the outside in: case, backend, thread count,
#  repeat. Putting the repeats innermost would run the same condition back to
#  back and bias the cache and frequency state, so the repeats are put
#  outermost here and each pass goes through the whole set once.
#---------------------------------------------------------------------------
total=0
for c in $CASES; do
  for be in $BACKENDS; do
    for th in $THREADS; do
      total=$((total + REPEATS))
    done
  done
done
echo "scaling.sh: $total runs planned (cases: $(echo $CASES | wc -w), backends: $BACKENDS, threads: $THREADS, repeats: $REPEATS)"

if [ $DRY -eq 1 ]; then
  for rep in $(seq 1 "$REPEATS"); do
    for c in $CASES; do
      for be in $BACKENDS; do
        for th in $THREADS; do
          if done_key "$c" "$be" "$th" "$rep"; then
            echo "  skip  $c $be th=$th r=$rep"
          else
            echo "  run   $c $be th=$th r=$rep"
          fi
        done
      done
    done
  done
  exit 0
fi

n=0
for rep in $(seq 1 "$REPEATS"); do
  for c in $CASES; do
    for be in $BACKENDS; do
      for th in $THREADS; do
        n=$((n + 1))
        if done_key "$c" "$be" "$th" "$rep"; then
          printf "  [%d/%d] skip (already measured) %s %s th=%s r=%s\n" "$n" "$total" "$c" "$be" "$th" "$rep"
          continue
        fi
        printf "  [%d/%d]\n" "$n" "$total"
        run_one "$c" "$be" "$th" "$rep"
      done
    done
  done
done

echo
echo "scaling.sh: done. results in $OUT"
echo "scaling.sh: summarize with  python3 tools/scaling/summarize.py $OUT"
