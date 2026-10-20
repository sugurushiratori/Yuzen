#!/bin/sh
# usage: bench.sh <input> [threads...]
#   Run the same case with varying thread counts and report wall-clock time.
#   The input may be either .toml or .inp (run_case.sh preserves the extension).
#   Note that specifying more threads than physical cores is extremely slow.
#
#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INP=$1; shift
[ $# -eq 0 ] && set -- 1 2
for t in "$@"; do
  D="$ROOT/runs/bench_t$t"
  s=$(date +%s.%N)
  # Through sh, so that a checkout which has not kept the executable bit
  # still works.
  sh "$ROOT/tools/run_case.sh" "$INP" "$D" "$t" >/dev/null 2>&1
  e=$(date +%s.%N)
  printf "threads=%-3s  %6.1f s\n" "$t" "$(echo "$e-$s"|bc)"
done
