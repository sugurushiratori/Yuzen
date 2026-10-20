#!/usr/bin/env python3
"""Compare two result directories numerically.

usage: compare.py <runA> <runB> [--job <result-dir-name>]

Matches up History.dat (scalar time series) and the final field*.dat
(spatial distribution) and reports the largest relative error per column.
Intended as a regression check, to confirm that a change to the code has
not changed the numbers.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import sys, os, glob


def read_table(path):
    with open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append([float(x) for x in line.split("\t")])
    return [h.strip() for h in head], rows


def reldiff(a, b):
    d = abs(a - b)
    s = max(abs(a), abs(b))
    return d if s < 1e-300 else d / s


def compare_table(pa, pb, label):
    """Compare column by column as (largest difference) / (largest absolute
    value in that column). Normalizing by the column scale avoids the
    relative error blowing up on entries close to zero."""
    ha, ra = read_table(pa)
    hb, rb = read_table(pb)
    if ha != hb:
        print(f"  [{label}] column layout differs")
        return None
    if len(ra) != len(rb):
        print(f"  [{label}] row count differs: {len(ra)} vs {len(rb)}")
        n = min(len(ra), len(rb))
    else:
        n = len(ra)
    # Iteration counts and residuals are diagnostics of how the solve
    # proceeded, not physical quantities. They may legitimately change, so
    # they are reported but excluded from the pass/fail decision.
    DIAG = {"NK_iter", "NK_resid", "BCG_iter", "BCG_resid"}
    worst = 0.0
    worst_col = ""
    diag_worst = 0.0
    ncol = min(len(ha), min(len(ra[0]), len(rb[0]))) if n else 0
    for j in range(ncol):
        scale = 0.0
        d = 0.0
        for i in range(n):
            scale = max(scale, abs(ra[i][j]), abs(rb[i][j]))
            d = max(d, abs(ra[i][j] - rb[i][j]))
        if scale > 0.0:
            r = d / scale
            if ha[j] in DIAG:
                diag_worst = max(diag_worst, r)
            elif r > worst:
                worst, worst_col = r, ha[j]
    msg = f"  [{label}] rows={n} max diff/column scale={worst:.3e}"
    msg += f"  ({worst_col})" if worst > 0 else "  (exact match)"
    if diag_worst > 0:
        msg += f"   [diagnostics (iterations/residual) diff={diag_worst:.1e}]"
    print(msg)
    return worst


def find_job(d):
    """The result directory inside a run directory.

    Each case sets JobName to its own file name, so the directory is named
    after the case rather than being one fixed name. Find it instead of
    assuming it: the one subdirectory holding a History.dat.
    """
    if not os.path.isdir(d):
        return None
    cand = [n for n in sorted(os.listdir(d))
            if os.path.isfile(os.path.join(d, n, "History.dat"))]
    if len(cand) == 1:
        return cand[0]
    if not cand:
        return None
    print("  !! %s holds more than one result directory (%s); "
          "name one with --job" % (d, ", ".join(cand)))
    return None


def main():
    a, b = sys.argv[1], sys.argv[2]
    job = None
    if "--job" in sys.argv:
        job = sys.argv[sys.argv.index("--job") + 1]
    else:
        job = find_job(a) or find_job(b)
        if job is None:
            print("  !! no result directory found under %s or %s" % (a, b))
            return 2
    da, db = os.path.join(a, job), os.path.join(b, job)
    print(f"comparing: {a}  vs  {b}")

    # Do not report "zero difference" when the directory or the target
    # files are missing. Passing over this silently makes a comparison
    # against a non-existent reference look like a pass.
    for d in (da, db):
        if not os.path.isdir(d):
            print(f"  !! directory missing: {d}")
            return 2

    worst = 0.0
    compared = 0
    for name in ("History.dat", "CalcInfo.dat"):
        pa, pb = os.path.join(da, name), os.path.join(db, name)
        if os.path.exists(pa) and os.path.exists(pb):
            w = compare_table(pa, pb, name)
            compared += 1
            if w is not None:
                worst = max(worst, w)
    fa = sorted(glob.glob(os.path.join(da, "field*.dat")))
    fb = sorted(glob.glob(os.path.join(db, "field*.dat")))
    if fa and fb:
        w = compare_table(fa[-1], fb[-1], os.path.basename(fa[-1]))
        compared += 1
        if w is not None:
            worst = max(worst, w)

    if compared == 0:
        print("  !! no comparable files at all")
        return 2

    print(f"  => overall max diff/column scale = {worst:.3e}  ({compared} files)")
    return 0 if worst == 0.0 else (1 if worst > 1e-10 else 0)


if __name__ == "__main__":
    sys.exit(main())
