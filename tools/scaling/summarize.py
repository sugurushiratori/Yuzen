#!/usr/bin/env python3
"""Read the TSV written by scaling.sh and tabulate speed-up and parallel efficiency.

    python3 tools/scaling/summarize.py results/scaling-host.tsv
    python3 tools/scaling/summarize.py results/*.tsv --format latex
    python3 tools/scaling/summarize.py results/*.tsv --format csv > summary.csv

Output formats
    text   table for reading in a terminal (default)
    latex  ready to paste into the manual's yzntable
    csv    for loading into a spreadsheet

Speed-up S(p) = T(1)/T(p), parallel efficiency E(p) = S(p)/p. The baseline is
the single-thread value for the same case and the same backend, taking the
**minimum** where there are repeats. The minimum rather than the mean is used
because the scatter in the measurements is one-sided noise from being
interrupted by other work (it cannot stray toward faster). The scatter itself
is reported in the spread column (max/min - 1), so do not trust a row where
that is large.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import math
import signal
import os
import sys
from collections import OrderedDict

NUM_COLS = ("threads", "repeat", "steps", "nk_total", "bcg_total")


def read_tsv(paths):
    """Read several TSVs and return (rows, meta). meta holds the # header per file."""
    rows = []
    meta = OrderedDict()
    for p in paths:
        m = {}
        head = None
        for line in open(p):
            line = line.rstrip("\n")
            if line.startswith("#"):
                bits = line[1:].strip().split("\t", 1)
                if len(bits) == 2:
                    m[bits[0].strip()] = bits[1].strip()
                continue
            if not line.strip():
                continue
            f = line.split("\t")
            if head is None:
                head = f
                continue
            r = dict(zip(head, f))
            r["_file"] = os.path.basename(p)
            rows.append(r)
        meta[os.path.basename(p)] = m
    return rows, meta


def to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def collect(rows, include_bad=False):
    """(case, backend, threads) -> {wall: [..], steps, nodes, dim, nk, bcg, bad}"""
    out = OrderedDict()
    for r in rows:
        ok = r.get("status", "") == "ok"
        if not ok and not include_bad:
            continue
        w = to_float(r.get("wall_s"))
        if w is None:
            continue
        key = (r["case"], r["backend"], int(r["threads"]))
        e = out.setdefault(key, {"wall": [], "nk": set(), "bcg": set(),
                                 "nodes": r.get("nodes"), "dim": r.get("dim"),
                                 "steps": r.get("steps"), "bad": 0})
        e["wall"].append(w)
        e["nk"].add(r.get("nk_total"))
        e["bcg"].add(r.get("bcg_total"))
        if not ok:
            e["bad"] += 1
    return out


def summarize(data):
    """Assemble the table rows. One row = (case, backend, threads)."""
    cases = []
    for (c, be, th) in data:
        if (c, be) not in cases:
            cases.append((c, be))

    table = []
    for (c, be) in cases:
        base = data.get((c, be, 1))
        t1 = min(base["wall"]) if base else None
        threads = sorted(th for (cc, bb, th) in data if cc == c and bb == be)
        for th in threads:
            e = data[(c, be, th)]
            t = min(e["wall"])
            spread = (max(e["wall"]) / t - 1.0) if t > 0 else float("nan")
            s = (t1 / t) if (t1 and t > 0) else None
            table.append({
                "case": c, "backend": be, "threads": th,
                "dim": e["dim"], "nodes": e["nodes"], "steps": e["steps"],
                "n": len(e["wall"]),
                "wall": t, "spread": spread,
                "speedup": s,
                "eff": (s / th) if s else None,
                "s_per_step": t / float(e["steps"]) if e["steps"] not in (None, "NA") else None,
                "nk": sorted(e["nk"])[0] if e["nk"] else "",
                "bcg": sorted(e["bcg"])[0] if e["bcg"] else "",
                "work_varies": len(e["nk"]) > 1 or len(e["bcg"]) > 1,
                "bad": e["bad"],
            })
    return table


def amdahl(table):
    """Fit Amdahl's parallel fraction f by least squares, per case and backend.

    S(p) = 1/((1-f) + f/p)  <=>  1/S = (1-f) + f*(1/p)
    Regress 1/S as a linear function of 1/p. Where there are only two points,
    the fit is no more than an order-of-magnitude guide.

    *** The fit breaks down if the saturated, slowing-down region is included ***
    Amdahl's law can express "adding more brings no further gain", but not
    "adding more makes it slower". If the regression includes the range where
    adding threads genuinely slows the run down (crossing NUMA nodes, too small
    a problem), the slope points the wrong way and f sticks at 0. In the
    112-core measurements the small 1-D cases behave this way: fitting all the
    points gives f = 0.000, whereas refitting over 14 threads and below gives
    0.4 to 0.75.
    Read a row with f = 0 as "does not follow Amdahl's form over this range",
    not as "nothing is parallel". The measurements themselves are the S and E
    columns.
    """
    out = OrderedDict()
    groups = OrderedDict()
    for r in table:
        if r["speedup"] is None or r["threads"] < 2:
            continue
        groups.setdefault((r["case"], r["backend"]), []).append(r)
    for k, rs in groups.items():
        xs = [1.0 / r["threads"] for r in rs]
        ys = [1.0 / r["speedup"] for r in rs]
        n = len(xs)
        if n < 2:
            continue
        mx = sum(xs) / n
        my = sum(ys) / n
        den = sum((x - mx) ** 2 for x in xs)
        if den == 0:
            continue
        slope = sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / den
        out[k] = max(0.0, min(1.0, slope))   # slope = f
    return out


def fmt(v, spec="%.3g", na="--"):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return na
    return spec % v


def print_text(table, meta, fits):
    for fn, m in meta.items():
        print("# %s" % fn)
        for k in ("host", "cpu", "sockets", "cores-per-socket",
                  "threads-per-core", "numa-nodes", "yuzen-version",
                  "omp-proc-bind", "date"):
            if k in m:
                print("#   %-18s %s" % (k, m[k]))
    print()
    hdr = ("case", "size", "backend", "th", "n", "wall[s]", "s/step",
           "S", "E", "spread", "NK", "BCG")
    print("%-16s %-11s %-7s %4s %3s %10s %10s %7s %7s %7s %8s %9s" % hdr)
    print("-" * 108)
    last = None
    for r in table:
        if last is not None and (r["case"], r["backend"]) != last:
            print()
        last = (r["case"], r["backend"])
        size = ("1-D N=%s" % r["nodes"]) if r["dim"] == "1" else \
               ("2-D %sx%s" % (int(math.isqrt(int(r["nodes"]))),
                               int(math.isqrt(int(r["nodes"])))))
        flag = ""
        if r["work_varies"]:
            flag += "  <- work differs between runs"
        if r["bad"]:
            flag += "  <- %d non-ok run(s) included" % r["bad"]
        print("%-16s %-11s %-7s %4d %3d %10s %10s %7s %7s %7s %8s %9s%s" % (
            r["case"], size, r["backend"], r["threads"], r["n"],
            fmt(r["wall"], "%.4g"), fmt(r["s_per_step"], "%.4g"),
            fmt(r["speedup"], "%.2f"), fmt(r["eff"], "%.2f"),
            fmt(r["spread"], "%.3f"),
            r["nk"], r["bcg"], flag))
    #  Some backends change their iteration count with the thread count
    #  (amgcl's preconditioner depends on the number of threads, so rounding
    #   differences shift the iterations by 1-2). That means the numerator and
    #  the denominator of the speed-up did different amounts of work, so print
    #  a note wherever the difference is large.
    groups = OrderedDict()
    for r in table:
        groups.setdefault((r["case"], r["backend"]), []).append(r)
    notes = []
    for (c, be), rs in groups.items():
        vals = [to_float(r["bcg"]) for r in rs if to_float(r["bcg"])]
        if len(set(vals)) > 1:
            lo, hi = min(vals), max(vals)
            notes.append("  %-16s %-7s  linear iterations vary with thread "
                         "count: %g .. %g (%.2f %%)" %
                         (c, be, lo, hi, 100.0 * (hi / lo - 1.0)))
    if notes:
        print()
        print("Work is not identical across thread counts:")
        for n in notes:
            print(n)

    if fits:
        print()
        print("Amdahl fit (parallel fraction f, from all measured thread counts):")
        print("  NOTE: f = 0.000 means the fit does not hold over this range")
        print("        (the run slows down past its optimum, which Amdahl's law")
        print("        cannot express), not that nothing is parallel.")
        for (c, be), f in fits.items():
            smax = (1.0 / (1.0 - f)) if f < 1 else float("inf")
            print("  %-16s %-7s  f = %.3f   S(inf) = %s" %
                  (c, be, f, ("%.1f" % smax) if smax != float("inf") else "inf"))


def print_csv(table):
    cols = ("case", "dim", "nodes", "backend", "threads", "n", "steps",
            "wall", "s_per_step", "speedup", "eff", "spread", "nk", "bcg", "bad")
    print(",".join(cols))
    for r in table:
        print(",".join("" if r.get(c) is None else str(r.get(c, "")) for c in cols))


def print_latex(table):
    """Ready to paste into the yzntable of docs/manual."""
    print(r"\begin{yzntable}{Wall-clock time, speed-up $S$ and parallel"
          r" efficiency $E = S/p$ against thread count $p$."
          r"\label{tab:scaling}}{@{}L{0.17\textwidth} L{0.13\textwidth}"
          r" L{0.10\textwidth} d d d d@{}}")
    print(r"\textsf{\bfseries Case} & \textsf{\bfseries Size} &"
          r" \textsf{\bfseries Backend} & \textsf{\bfseries $p$} &"
          r" \textsf{\bfseries $T$ [s]} & \textsf{\bfseries $S$} &"
          r" \textsf{\bfseries $E$} \\ \midrule")
    last = None
    for r in table:
        size = ("1-D, $N_x=%s$" % r["nodes"]) if r["dim"] == "1" else \
               ("2-D, $%d^2$" % int(math.isqrt(int(r["nodes"]))))
        c = r["case"].replace("_", r"\_")
        if (r["case"], r["backend"]) == last:
            cell_case, cell_size, cell_be = "", "", ""
        else:
            if last is not None and r["case"] != last[0]:
                print(r"\addlinespace")
            cell_case = r"\file{%s}" % c
            cell_size = size
            cell_be = "\\val{\"%s\"}" % r["backend"]
        last = (r["case"], r["backend"])
        print("%s & %s & %s & %d & %s & %s & %s \\\\" % (
            cell_case, cell_size, cell_be, r["threads"],
            fmt(r["wall"], "%.4g"), fmt(r["speedup"], "%.2f"),
            fmt(r["eff"], "%.2f")))
    print(r"\end{yzntable}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="+")
    ap.add_argument("--format", choices=("text", "csv", "latex"), default="text")
    ap.add_argument("--include-bad", action="store_true",
                    help="include runs whose status is not ok (discarded by default)")
    args = ap.parse_args()

    rows, meta = read_tsv(args.tsv)
    if not rows:
        sys.exit("no data rows found")
    data = collect(rows, args.include_bad)
    if not data:
        sys.exit("no usable rows (all runs had status != ok?)  use --include-bad to see them")
    table = summarize(data)

    if args.format == "csv":
        print_csv(table)
    elif args.format == "latex":
        print_latex(table)
    else:
        print_text(table, meta, amdahl(table))


if __name__ == "__main__":
    #  So that `| head` does not make it crash
    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except (AttributeError, ValueError):
        pass
    main()
