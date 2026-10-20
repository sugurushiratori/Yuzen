#!/usr/bin/env python3
"""How the run time goes with the thread count, under variants of the OpenMP
runtime, its settings and the thread counts by phase (YUZEN_GRAIN).

    python3 tools/thread_variants.py --exe build-amg/bin/yuzen --variants variants.tsv \\
            --workdir runs/tv --threads "1 2 4 8 16 28 42 56" --resume --budget 1500

The variants file has one variant a line, tab separated:

    name    exe    VAR=value VAR=value ...    [threads]

"exe" is "-" for --exe, else the solver to run (a second build); the third
field is the environment of the variant, on top of the one every run has
(OMP_PROC_BIND=close, OMP_PLACES=cores, YUZEN_PROFILE=1, YUZEN_THREAD_CAP=0,
OMP_WAIT_POLICY and YUZEN_GRAIN unset; "-VAR" unsets VAR); the fourth, if there, the thread
counts of the variant instead of --threads (those of --threads it has), for
a variant that cannot differ from the first at some counts (one thread).
Lines starting with # are skipped.
For example:

    base     -                     (nothing)
    active   -                     OMP_WAIT_POLICY=ACTIVE
    iomp     -                     LD_PRELOAD=/path/to/libiomp5.so
    grain    -                     YUZEN_GRAIN=lin=1000,amg=1000

Each case below is run with each backend, at each thread count, under each
variant, --repeats times (all variants once before any is repeated). The
runs of profile.tsv are not repeated with --resume; --budget works as in
tools/case_runtimes.py (exit code 3: runs left for the next call). A run
over --timeout is recorded with rc "timeout" and not repeated.

The call that completes the set writes summary.txt: for each case and
backend, the loop time [s] (least over the repeats) against the thread count
for each variant, and

    best     the least time and its thread count
    rise     the largest t(T) / min over T' < T of t(T'): 1.00 when the time
             never grows with the thread count, 1.30 when it grows by 30 %
             somewhere
    vs base  the best time of the variant against the best of the first
             variant of the file

A count a variant leaves out is shown, and counted, with the time of the
first variant that runs the same solver at every count.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from case_runtimes import host_info          # noqa: E402
from profile_threads import run, PARTS       # noqa: E402

#  (case, steps or None for the case as shipped, backends): the cases whose
#  time grew with the thread count within a socket, the large ones and CZ
CASES = [
    ("leveling_1d_reg", None, ("eigen", "amgcl")),
    ("leveling_2d_reg", 100, ("eigen", "amgcl")),
    ("scaling/bench1d_n08000", 100, ("eigen", "amgcl")),
    ("scaling/bench2d_n0128", 30, ("eigen", "amgcl")),
    ("scaling/bench2d_n0256", 10, ("eigen", "amgcl")),
    ("scaling/bench2d_n0512", 4, ("eigen", "amgcl")),
    ("scaling/bench2d_n1024", 2, ("amgcl",)),
    ("concz_2d", 20, ("eigen", "amgcl")),
]
COLS = ["variant", "rep", "case", "backend", "threads", "used", "rc", "wall_s", "loop_s", "parts"]


def read_variants(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            v = line.split("\t")
            name, exe = v[0].strip(), (v[1].strip() if len(v) > 1 else "-")
            env = {}
            for tok in (v[2].split() if len(v) > 2 else []):
                if tok.startswith("-") and "=" not in tok:
                    env[tok[1:]] = None             # unset it
                    continue
                k, _, val = tok.partition("=")
                env[k] = val
            thr = [int(t) for t in v[3].split()] if len(v) > 3 and v[3].strip() else None
            out.append((name, exe, env, thr))
    return out


def read_rows(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        cols = None
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.rstrip("\n").split("\t")
            if cols is None:
                cols = v
                continue
            if len(v) == len(cols):
                rows.append(dict(zip(cols, v)))
    return rows


def summary(rows, variants, threads):
    out = []
    names = [v[0] for v in variants]
    keys = []
    for case, _, bes in CASES:
        for be in bes:
            keys.append((case, be))
    overall = []
    for case, be in keys:
        sel = [d for d in rows if d["case"] == case and d["backend"] == be and d["rc"] == "0" and d["loop_s"]]
        if not sel:
            continue
        out.append("")
        out.append("%s, %s: loop time [s] by thread count (least of the repeats)" % (case, be))
        out.append("%-14s" % "variant" + "".join("%8d" % t for t in threads)
                   + "   best      @T   rise  vs base")
        base_best = None
        for name in names:
            t_of = {}
            for d in sel:
                if d["variant"] != name:
                    continue
                t, v = int(d["threads"]), float(d["loop_s"])
                t_of[t] = min(v, t_of.get(t, v))
            if not t_of:
                continue
            #  A count the variant leaves out (fourth field) is taken from the
            #  first variant when both run the same solver: the same run
            v = [x for x in variants if x[0] == name][0]
            full = [x[0] for x in variants if x[1] == v[1] and x[3] is None]
            if v[3] is not None and full:
                for d in sel:
                    t = int(d["threads"])
                    if d["variant"] == full[0] and t not in v[3]:
                        t_of[t] = min(float(d["loop_s"]), t_of.get(t, float(d["loop_s"])))
            ts = [t for t in threads if t in t_of]
            best_t = min(ts, key=lambda t: t_of[t])
            best = t_of[best_t]
            rise, low = 1.0, None
            for t in ts:
                if low is not None:
                    rise = max(rise, t_of[t] / low)
                low = t_of[t] if low is None else min(low, t_of[t])
            if base_best is None:
                base_best = best
            vs = base_best / best if best else 0.0
            out.append("%-14s" % name + "".join("%8s" % ("%.3g" % t_of[t] if t in t_of else "-") for t in threads)
                       + "  %6.3g  %4d  %5.2f  %5.2f" % (best, best_t, rise, vs))
            overall.append((case, be, name, best, best_t, rise, vs, t_of.get(max(ts)) / best))
    #  The parts of each run at the largest thread count, base and every variant
    out.append("")
    out.append("Parts [s] at the largest thread count measured (first repeat)")
    for case, be in keys:
        sel = [d for d in rows if d["case"] == case and d["backend"] == be and d["rc"] == "0" and d["rep"] == "0"]
        if not sel:
            continue
        tmax = max(int(d["threads"]) for d in sel)
        out.append("")
        out.append("%s, %s, %d threads" % (case, be, tmax))
        line = "%-12s" % "part" + "".join("%11s" % n[:10] for n in names)
        out.append(line)
        per = {}
        for d in sel:
            if int(d["threads"]) == tmax:
                per[d["variant"]] = dict(p.split("=") for p in d["parts"].split(";") if "=" in p)
                per[d["variant"]]["loop"] = "0," + d["loop_s"]
        for part in ["loop"] + PARTS:
            vals = [per.get(n, {}).get(part, "").split(",")[-1] for n in names]
            if not any(vals):
                continue
            out.append("%-12s" % part + "".join("%11s" % ("%.3g" % float(v) if v else "-") for v in vals))
    #  One line per variant: how many case-backend pairs rise, and by how much
    out.append("")
    out.append("By variant: pairs whose time grows with the thread count (rise > 1.05), the worst rise,")
    out.append("the geometric mean of the best time against the base, and of t(largest T) / best")
    out.append("%-14s %6s %8s %10s %10s" % ("variant", "rising", "worst", "best/base", "tmax/best"))
    for name in names:
        r = [o for o in overall if o[2] == name]
        if not r:
            continue
        rising = sum(1 for o in r if o[5] > 1.05)
        worst = max(o[5] for o in r)
        gm = 1.0
        gt = 1.0
        for o in r:
            gm *= o[6]
            gt *= o[7]
        gm **= 1.0 / len(r)
        gt **= 1.0 / len(r)
        out.append("%-14s %3d/%-2d %8.2f %10.2f %10.2f" % (name, rising, len(r), worst, gm, gt))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", default="1 2 4 8 16 28 42 56")
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--only", action="append", help="this case only (repeatable)")
    ap.add_argument("--only-variant", action="append", help="this variant only (repeatable)")
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--budget", type=float, default=None)
    args = ap.parse_args()

    t_start = time.perf_counter()
    work = os.path.abspath(args.workdir)
    os.makedirs(work, exist_ok=True)
    cores = sorted(os.sched_getaffinity(0))
    threads = [int(t) for t in args.threads.split()]
    variants = read_variants(args.variants)
    if not variants:
        print("no variants in %s" % args.variants)
        return 2
    for name, exe, env, _ in variants:
        e = os.path.abspath(args.exe if exe == "-" else exe)
        if not os.access(e, os.X_OK):
            print("variant %s: %s is not there" % (name, e))
            return 2
        if "LD_PRELOAD" in env and not os.path.exists(env["LD_PRELOAD"]):
            print("variant %s: LD_PRELOAD %s is not there" % (name, env["LD_PRELOAD"]))
            return 2
    path = os.path.join(work, "profile.tsv")
    rows = read_rows(path) if args.resume else []
    done = {(d["variant"], d["rep"], d["case"], d["backend"], d["threads"]) for d in rows}
    info = host_info(os.path.abspath(args.exe))
    if rows:
        tsv = open(path, "a")
        tsv.write("# resumed\t%s\n" % info["date"])
    else:
        tsv = open(path, "w")
        for k, v in info.items():
            tsv.write("# %s\t%s\n" % (k, v))
        for name, exe, env, thr in variants:
            tsv.write("# variant\t%s\t%s\t%s\t%s\n" % (name, exe, " ".join("-%s" % k if v is None else "%s=%s" % (k, v) for k, v in env.items()),
                                                      " ".join(map(str, thr)) if thr else ""))
        tsv.write("\t".join(COLS) + "\n")
    tsv.flush()
    print("host %s, %s, %d cores allowed; %s" % (info["host"], info["cpu"], len(cores), info["yuzen"]),
          flush=True)

    left_over, too_long = 0, []
    for rep in range(args.repeats):
        for name, exe, env, vthr in variants:
            if args.only_variant and name not in args.only_variant:
                continue
            e = os.path.abspath(args.exe if exe == "-" else exe)
            for case, steps, backends in CASES:
                if args.only and case not in args.only:
                    continue
                for be in backends:
                    for t in threads:
                        if vthr is not None and t not in vthr:
                            continue
                        key = (name, str(rep), case, be, str(t))
                        if t > len(cores) or key in done:
                            continue
                        timeout = args.timeout
                        if args.budget is not None:
                            left = args.budget - (time.perf_counter() - t_start)
                            if left < 5.0:
                                left_over += 1
                                continue
                            timeout = min(timeout, left)
                        started = time.perf_counter() - t_start
                        rundir = os.path.join(work, "run_%s_%s_%s_%d" % (name, case.replace("/", "_"), be, t))
                        r = run(e, case, steps, be, t, cores[:t], rundir, timeout, "default", "default", env)
                        shutil.rmtree(rundir, ignore_errors=True)
                        if r is None and timeout >= args.timeout:
                            #  Over --timeout itself, not cut by the budget: a
                            #  failed run, not one for the next call
                            r = {"rc": "timeout", "wall": timeout, "loop": "", "used": "", "parts": {}}
                        if r is None:
                            left_over += 1
                            if args.budget is not None and started < 5.0:
                                too_long.append(key)
                            continue
                        parts = ";".join("%s=%d,%.6f" % (k, c, s) for k, (c, s) in r["parts"].items())
                        row = {"variant": name, "rep": str(rep), "case": case, "backend": be, "threads": str(t),
                               "used": r["used"], "rc": str(r["rc"]), "wall_s": "%.3f" % r["wall"],
                               "loop_s": r["loop"], "parts": parts}
                        tsv.write("\t".join(row[c] for c in COLS) + "\n")
                        tsv.flush()
                        rows.append(row)
                        print("%-10s r%d %-24s %-6s thr %-3d loop %8s s  wall %8.2f s  rc %s"
                              % (name, rep, case, be, t, r["loop"], r["wall"], r["rc"]), flush=True)
    tsv.close()
    if too_long:
        print("%s did not finish within --budget, which it had to itself: raise --budget" % too_long)
        return 2
    if left_over:
        print("budget spent: %d run(s) left for the next call (--resume)" % left_over)
        return 3
    lines = ["host %s; %s; %d cores allowed" % (info["host"], info["cpu"], len(cores)), info["yuzen"]]
    lines += summary(rows, variants, threads)
    with open(os.path.join(work, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 1 if any(r["rc"] != "0" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
