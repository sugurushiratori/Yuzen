#!/usr/bin/env python3
"""Where the time of a run goes, at each thread count (YUZEN_PROFILE).

    python3 tools/profile_threads.py --exe build-amg/bin/yuzen --workdir runs/prof
    python3 tools/profile_threads.py --exe build-amg/bin/yuzen --workdir runs/prof \\
            --threads "1 2 4 8 16 28 56 112" --resume --budget 1500

Each case of the set below is run with each backend at each thread count,
with YUZEN_PROFILE=1, so that the solver reports at exit the wall-clock time
spent in its named parts (yuzen/src/util/Profile.h):

    step          one time step, all of it (TimeAdvance)
    output        the field and history output between steps
    h.assemble    the film-height matrix A (SetMatrixA)
    h.jacobian    its Jacobian and the matrix of the implicit step
    h.rhs         the explicit part and the right-hand side
    h.linsolve    the linear solve of the film height
    h.couple      what each Newton pass does besides (the well-mixed
                  concentration trial among it)
    h.after       the end of the step (the concentration of the step)
    c.solve       a concentration solve, all of it; c.linsolve its linear solve
    velocity      the velocities
    lin.prepare   copying and scaling the matrix for the backend
    lin.precond   building the preconditioner (ILUT factorization, AMG setup)
    lin.iterate   the BiCGStab iterations

Times of nested parts are each inclusive (h.linsolve contains lin.*). The
threads are bound to cores 0..T-1 (OMP_PROC_BIND=close, OMP_PLACES=cores).
Each run is made with every OMP_WAIT_POLICY of --wait (default: left unset,
the runtime's own, with which GCC spins a while before it sleeps; PASSIVE,
which ctest sets, makes the threads sleep between parallel regions, and each
region then costs a wake-up of tens of microseconds), and with every value
of --opt: "default" leaves YUZEN_OPT unset, "none" turns the algorithms of
util/Opt.h off, to compare against the ones they replaced. The solver's cap on the
thread count of small grids (util/ThreadCap.h) is turned off
(YUZEN_THREAD_CAP=0), so that each run uses the count it is given.
The cases are cut to a fixed number of steps so that the set fits a short
batch job; the work of a run does not depend on the thread count.

--resume and --budget work as in tools/case_runtimes.py: runs already in
profile.tsv are not repeated, no run is started (or let finish) past the
budget, and a call that leaves runs undone exits with code 3. The call that
completes the set writes summary.txt: for each case and backend, the seconds
of each part against the thread count, and the speed-up of the step.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from case_runtimes import prepare, host_info   # noqa: E402

#  (case, steps or None for the case as shipped, backends)
CASES = [
    ("leveling_1d_reg", None, ("eigen", "amgcl")),
    ("leveling_2d_reg", 100, ("eigen", "amgcl")),
    ("scaling/bench1d_n08000", 100, ("eigen", "amgcl")),
    ("scaling/bench2d_n0128", 30, ("eigen", "amgcl")),
    ("scaling/bench2d_n0256", 10, ("eigen", "amgcl")),
    ("scaling/bench2d_n0512", 4, ("eigen", "amgcl")),
    ("scaling/bench2d_n1024", 2, ("eigen", "amgcl")),
    ("concz_2d", 20, ("eigen", "amgcl")),
]
#  Too slow with every algorithm of util/Opt.h off (the old ILUT on a large
#  2-D film: minutes a step); their "none" is left out
SKIP_NONE = {("scaling/bench2d_n0512", "eigen"), ("scaling/bench2d_n1024", "eigen")}
PARTS = ["step", "output", "h.assemble", "h.jacobian", "h.rhs", "h.linsolve", "h.couple", "h.after",
         "c.solve", "c.linsolve", "velocity", "lin.prepare", "lin.precond", "lin.iterate"]


def run(exe, case, steps, backend, threads, cores, rundir, timeout, wait, opt="default", extra_env=None):
    name = prepare(case, rundir, backend)
    path = os.path.join(rundir, name)
    if steps is not None:
        with open(path, encoding="utf-8") as f:
            text = f.read()
        text, n = re.subn(r"^MaxSteps = .*$", "MaxSteps = %d" % steps, text, count=1, flags=re.M)
        if n != 1:
            text = re.sub(r"^\[Numerics\]\n", "[Numerics]\nMaxSteps = %d\n" % steps, text, count=1, flags=re.M)
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PROC_BIND="close", OMP_PLACES="cores",
               YUZEN_PROFILE="1", YUZEN_THREAD_CAP="0")
    env.pop("OMP_WAIT_POLICY", None)
    if wait != "default":
        env["OMP_WAIT_POLICY"] = wait
    env.pop("YUZEN_OPT", None)
    if opt != "default":
        env["YUZEN_OPT"] = opt
    for k, v in (extra_env or {}).items():
        if v is None:
            env.pop(k, None)        # "-VAR" in a variant: unset
        else:
            env[k] = v
    t0 = time.perf_counter()
    try:
        r = subprocess.run([exe, name], cwd=rundir, capture_output=True, text=True, encoding="utf-8", errors="replace", env=env,
                           timeout=timeout, preexec_fn=lambda: os.sched_setaffinity(0, cores))
        rc, out = r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return None
    wall = time.perf_counter() - t0
    parts = {}
    for m in re.finditer(r"^\[profile\] (\S+)\s+(\d+)\s+([0-9.eE+-]+)$", out, re.M):
        parts[m.group(1)] = (int(m.group(2)), float(m.group(3)))
    loop, used = "", ""
    try:
        job = os.path.join(rundir, os.path.basename(case))
        with open(os.path.join(job, "JobInfo.txt"), errors="replace", encoding="utf-8") as f:
            ji = f.read()
        m = re.search(r"^Time loop \[s\] :\s*(\S+)", ji, re.M)
        loop = m.group(1) if m else ""
        m = re.search(r"^Linear solver <height> :\s*(\S+)", ji, re.M)
        used = m.group(1) if m else ""
    except OSError:
        pass
    return {"rc": rc, "wall": wall, "loop": loop, "used": used, "parts": parts}


def read_done(path):
    done, rows = set(), []
    if not os.path.exists(path):
        return done, rows
    with open(path, encoding="utf-8") as f:
        cols = None
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.rstrip("\n").split("\t")
            if cols is None:
                cols = v
                continue
            if len(v) != len(cols):
                continue
            d = dict(zip(cols, v))
            rows.append(d)
            d.setdefault("opt", "default")
            done.add((d["case"], d["backend"], int(d["threads"]), d["wait"], d["opt"]))
    return done, rows


def summary(rows, threads):
    out = []
    keys = []
    for d in rows:
        k = (d["case"], d["backend"], d["wait"], d.get("opt", "default"))
        if k not in keys:
            keys.append(k)
    for case, be, wait, opt in keys:
        sel = {int(d["threads"]): d for d in rows
               if d["case"] == case and d["backend"] == be and d["wait"] == wait
               and d.get("opt", "default") == opt}
        ts = [t for t in threads if t in sel]
        if not ts:
            continue
        out.append("")
        out.append("%s, %s (used %s), OMP_WAIT_POLICY %s, YUZEN_OPT %s: seconds by thread count"
                   % (case, be, sel[ts[0]]["used"], wait, opt))
        out.append("%-12s" % "part" + "".join("%9d" % t for t in ts))
        parts = {t: dict(p.split("=") for p in sel[t]["parts"].split(";") if "=" in p) for t in ts}
        for name in ["loop"] + PARTS:
            vals = []
            for t in ts:
                if name == "loop":
                    vals.append(sel[t]["loop_s"])
                else:
                    vals.append(parts[t].get(name, "").split(",")[-1] if name in parts[t] else "")
            if not any(vals):
                continue
            out.append("%-12s" % name + "".join("%9s" % ("%.3g" % float(v) if v else "-") for v in vals))
        base = float(sel[ts[0]]["loop_s"] or 0)
        out.append("%-12s" % "speed-up" + "".join(
            "%9s" % ("%.2f" % (base / float(sel[t]["loop_s"])) if base and sel[t]["loop_s"] else "-") for t in ts))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", default="1 2 4 8 16 28 56 112")
    ap.add_argument("--only", action="append", help="this case only (repeatable)")
    ap.add_argument("--wait", default="default",
                    help='OMP_WAIT_POLICY values to run each with; "default" leaves it unset'
                         ' (default "default"; "default PASSIVE" runs both)')
    ap.add_argument("--opt", default="default",
                    help='YUZEN_OPT values to run each with; "default" leaves it unset'
                         ' (default "default"; "default none" also runs every algorithm of'
                         ' util/Opt.h off)')
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--budget", type=float, default=None)
    args = ap.parse_args()

    t_start = time.perf_counter()
    exe = os.path.abspath(args.exe)
    work = os.path.abspath(args.workdir)
    os.makedirs(work, exist_ok=True)
    cores = sorted(os.sched_getaffinity(0))
    threads = [int(t) for t in args.threads.split()]
    path = os.path.join(work, "profile.tsv")
    cols = ["case", "backend", "threads", "wait", "opt", "used", "rc", "wall_s", "loop_s", "parts"]
    waits = args.wait.split()
    opts = args.opt.split()
    done, rows = read_done(path) if args.resume else (set(), [])
    info = host_info(exe)
    if rows:
        tsv = open(path, "a", encoding="utf-8")
        tsv.write("# resumed\t%s\n" % info["date"])
    else:
        tsv = open(path, "w", encoding="utf-8")
        for k, v in info.items():
            tsv.write("# %s\t%s\n" % (k, v))
        tsv.write("\t".join(cols) + "\n")
    tsv.flush()
    print("host %s, %s, %d cores allowed; %s" % (info["host"], info["cpu"], len(cores), info["yuzen"]),
          flush=True)

    left_over, too_long = 0, []
    for opt, wait in [(o, w) for o in opts for w in waits]:
      for case, steps, backends in CASES:
        if args.only and case not in args.only:
            continue
        for be in backends:
            for t in threads:
                if t > len(cores) or (case, be, t, wait, opt) in done:
                    continue
                if opt == "none" and (case, be) in SKIP_NONE:
                    continue
                timeout = args.timeout
                if args.budget is not None:
                    left = args.budget - (time.perf_counter() - t_start)
                    if left < 5.0:
                        left_over += 1
                        continue
                    timeout = min(timeout, left)
                started = time.perf_counter() - t_start
                rundir = os.path.join(work, "run_%s_%s_%d_%s_%s" % (case.replace("/", "_"), be, t, wait,
                                                                     opt.replace(",", "_")))
                r = run(exe, case, steps, be, t, cores[:t], rundir, timeout, wait, opt)
                shutil.rmtree(rundir, ignore_errors=True)
                if r is None:
                    left_over += 1
                    if args.budget is not None and started < 5.0:
                        too_long.append((case, be, t, wait, opt))
                    continue
                parts = ";".join("%s=%d,%.6f" % (k, c, s) for k, (c, s) in r["parts"].items())
                row = {"case": case, "backend": be, "threads": str(t), "wait": wait, "opt": opt,
                       "used": r["used"], "rc": str(r["rc"]),
                       "wall_s": "%.3f" % r["wall"], "loop_s": r["loop"], "parts": parts}
                tsv.write("\t".join(row[c] for c in cols) + "\n")
                tsv.flush()
                rows.append(row)
                print("%-24s %-6s thr %-3d %-8s %-8s loop %8s s  wall %8.2f s  rc %d"
                      % (case, be, t, wait, opt, r["loop"], r["wall"], r["rc"]), flush=True)
    tsv.close()
    if too_long:
        print("%s did not finish within --budget, which it had to itself: raise --budget" % too_long)
        return 2
    if left_over:
        print("budget spent: %d run(s) left for the next call (--resume)" % left_over)
        return 3
    lines = ["host %s; %s; %d cores allowed" % (info["host"], info["cpu"], len(cores)), info["yuzen"]]
    lines += summary(rows, threads)
    with open(os.path.join(work, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 1 if any(r["rc"] != "0" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
