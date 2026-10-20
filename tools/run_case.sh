#!/bin/sh
# usage: run_case.sh <input> <outdir> [threads]
#
#   Copy the input file into the working directory and run it there.
#   Yuzen determines the input format from the extension, so the copy
#   **keeps the extension**.
#   (A copy named .inp would send a .toml file to the .inp parser, which
#    stops on a syntax error.)
#   Any .dat files referenced by Restart are carried along too.
#
#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0
set -e
INP=$1; OUT=$2; TH=${3:-1}
[ -n "$INP" ] && [ -n "$OUT" ] || { echo "usage: $0 <input> <outdir> [threads]" >&2; exit 2; }
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRCDIR=$(cd "$(dirname "$INP")" && pwd)
EXT=${INP##*.}
rm -rf "$OUT"; mkdir -p "$OUT"
cp "$INP" "$OUT/case.$EXT"
# pick up the .dat files referenced by Restart
for d in $(grep -oE '"[^"]+\.dat"' "$INP" | tr -d '"'); do
  [ -e "$SRCDIR/$d" ] && cp "$SRCDIR/$d" "$OUT/"
done
cd "$OUT"
OMP_NUM_THREADS=$TH "$ROOT/build/bin/yuzen" "case.$EXT" > stdout.txt 2>&1
