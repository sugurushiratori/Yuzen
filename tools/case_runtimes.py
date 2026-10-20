#!/usr/bin/env python3
"""Measure the run time of every shipped case, for the manual's tables.

    python3 tools/case_runtimes.py --exe build/bin/yuzen --workdir runs/rt
    python3 tools/case_runtimes.py --exe build/bin/yuzen --workdir runs/rt --repeats 3 --jobs 28
    python3 tools/case_runtimes.py --exe build-amg/bin/yuzen --workdir runs/rt-be \\
            --backends-only --threads "1 2 4 8 16 28 56 112"

What it measures:

  cases     Every case under cases/ and cases/eq/ (28), single-threaded, with
            Backend left at "auto": the CPU time of the run (user + system of
            the child process, from wait4), the wall-clock time, the number
            of steps (JobInfo.txt), how it ended (EndTime / Steady state /
            failed) and the number of warnings (lines of Warnings.txt). This
            is tab:cases-runtimes of the manual (appendix "The example cases").
  backends  With --backends: leveling_1d_reg and leveling_2d_reg with
            Backend = "eigen" and "amgcl" at each thread count of --threads,
            wall-clock time. This is tab:linsolve-timing (chapter "Linear
            solvers"). Needs a build with -DYUZEN_WITH_AMGCL=ON; a build
            without it falls back to eigen, which the script detects from
            JobInfo.txt and marks.

Each measurement is repeated --repeats times; the table reports the minimum
(the run least disturbed by the rest of the machine) and the TSV keeps every
repeat. Independent single-threaded runs may run side by side with --jobs;
each is pinned to a core of its own (os.sched_setaffinity), spread over the
cores the script itself is allowed to use. Runs with several threads are
always run one at a time, on cores 0..T-1 of that set. OMP_WAIT_POLICY is
left unset (removed from the environment if set), since no run shares its
cores.

Output, in --workdir:
  runtimes.tsv        one row per run, with the host description at the top
  summary.txt         the tables as text
  table_cases.tex     the rows of tab:cases-runtimes, ready to paste
  table_backends.tex  the rows of tab:linsolve-timing (with --backends)

The run directories are deleted after each run unless --keep is given.

Measuring in pieces (a batch queue with a short time limit): with --resume
the runs already in runtimes.tsv are not repeated and new ones are appended;
with --budget SECONDS no run is started past that much time from the start,
and a run still going when the budget is spent is stopped and left for the
next call (it is not recorded). A call that leaves runs undone exits with
code 3 and writes no tables; the call that completes the set writes them
from every row of runtimes.tsv. So

    python3 tools/case_runtimes.py ... --resume --budget 1500

repeated until it exits 0 gives the same tables as one long call.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import concurrent.futures
import datetime
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = os.path.join(ROOT, "cases")

#  The order of tab:cases-runtimes.
ORDER = ["smoke_1d", "leveling_1d_reg", "periodic_1d_reg", "leveling_1d_logout", "tutorial_1d", "allphys_1d_reg",
         "allphys_1d_neg_reg", "concz_1d_reg", "nozzle_1d_offset", "nozzle_1d_fixed",
         "nozzle_1d_rotating", "allphys_2d_reg", "leveling_2d_reg", "leveling_2d_tight_reg",
         "leveling_2d_fine_reg", "concz_2d", "stress_all", "stress_dt", "stress_coupling",
         "stress_evap", "stress_wetting", "wet_coupled", "wet_static_5", "wet_static_15", "wet_static",
         "eq/eq_run", "eq/spread_n500", "eq/eq_n1000", "eq/eq_n2000", "eq/eq_n4000"]
BACKEND_CASES = ["leveling_1d_reg", "leveling_2d_reg"]


def all_cases():
    found = []
    for sub in ("", "eq"):
        d = os.path.join(CASES, sub)
        for f in sorted(os.listdir(d)):
            if f.endswith(".toml"):
                found.append((sub + "/" if sub else "") + f[:-5])
    rank = {c: i for i, c in enumerate(ORDER)}
    return sorted(found, key=lambda c: (rank.get(c, len(ORDER)), c))


def host_info(exe):
    info = {"host": socket.gethostname(), "date": datetime.datetime.now().isoformat(timespec="seconds"),
            "platform": platform.platform(), "python": platform.python_version()}
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as f:
            m = re.search(r"^model name\s*:\s*(.*)$", f.read(), re.M)
            info["cpu"] = m.group(1).strip() if m else "?"
    except OSError:
        info["cpu"] = platform.processor() or "?"
    info["cores_allowed"] = len(os.sched_getaffinity(0))
    try:
        info["lscpu"] = "; ".join(l.strip() for l in subprocess.run(
            ["lscpu"], capture_output=True, text=True, encoding="utf-8", errors="replace").stdout.splitlines()
            if re.match(r"(Socket|Core|Thread|NUMA node\(s\)|CPU max MHz|Model name)", l))
    except OSError:
        pass
    v = subprocess.run([exe, "--version"], capture_output=True, text=True, encoding="utf-8", errors="replace", cwd=ROOT)
    info["yuzen"] = " ".join((v.stdout + v.stderr).split())[:200]
    return info


def prepare(case, rundir, backend=None):
    """Copy the case and the .dat files it refers to into rundir; with
    backend, set Backend of both linear solvers to it."""
    src = os.path.join(CASES, case + ".toml")
    srcdir = os.path.dirname(src)
    if os.path.exists(rundir):
        shutil.rmtree(rundir)
    os.makedirs(rundir)
    with open(src, encoding="utf-8") as f:
        text = f.read()
    for d in set(re.findall(r'"([^"]+\.dat)"', text)):
        if os.path.exists(os.path.join(srcdir, d)):
            shutil.copy(os.path.join(srcdir, d), rundir)
    if backend is not None:
        for blk in ("Numerics.HeightEquation.LinearSolver", "Numerics.ConcEquation.LinearSolver"):
            pat = re.compile(r"(^\[" + re.escape(blk) + r"\]\n(?:(?!\[).*\n)*?)Backend = .*\n", re.M)
            if pat.search(text):
                text = pat.sub(r'\1Backend = "%s"\n' % backend, text, count=1)
            else:
                text = re.sub(r"^\[" + re.escape(blk) + r"\]\n", r'[%s]\nBackend = "%s"\n' % (blk, backend),
                              text, count=1, flags=re.M)
    name = os.path.basename(case) + ".toml"
    with open(os.path.join(rundir, name), "w", encoding="utf-8") as f:
        f.write(text)
    return name


def run_one(exe, case, rundir, threads, cores, backend=None, keep=False, timeout=None):
    name = prepare(case, rundir, backend)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PROC_BIND="close", OMP_PLACES="cores")
    # The runtime's own wait policy: every run has its cores to itself, and
    # PASSIVE puts a sleep and a wake-up into every parallel loop (25 us at
    # 8 threads, 258 us at 112 on a Xeon Max 9480, against 1.7 and 24 us).
    env.pop("OMP_WAIT_POLICY", None)
    cores = list(cores)

    def pin():
        os.sched_setaffinity(0, cores)
    t0 = time.perf_counter()
    with open(os.path.join(rundir, "stdout.txt"), "w", encoding="utf-8") as out:
        p = subprocess.Popen([exe, name], cwd=rundir, stdout=out, stderr=subprocess.STDOUT, env=env,
                             preexec_fn=pin)
        deadline = None if timeout is None else t0 + timeout
        while True:
            pid, status, ru = os.wait4(p.pid, os.WNOHANG)
            if pid:
                break
            if deadline and time.perf_counter() > deadline:
                p.kill()
                pid, status, ru = os.wait4(p.pid, 0)
                break
            time.sleep(0.01)
    wall = time.perf_counter() - t0
    rc = os.waitstatus_to_exitcode(status)
    with open(os.path.join(rundir, "stdout.txt"), errors="replace", encoding="utf-8") as f:
        log = f.read()
    if "is more threads than this run can use" in log:
        end = "refused"         # above the solver's thread limit (util/ThreadCap.h)
    elif "Time reached to the EndTime" in log:
        end = "EndTime"
    elif "Steady state" in log:
        end = "Steady state"
    else:
        end = "failed (exit %d)" % rc
    job = os.path.join(rundir, os.path.basename(case))
    steps, warns, used, tused = "", 0, "", ""
    try:
        with open(os.path.join(job, "JobInfo.txt"), errors="replace", encoding="utf-8") as f:
            ji = f.read()
        m = re.search(r"^Time steps :\s*(\d+)", ji, re.M)
        steps = m.group(1) if m else ""
        m = re.search(r"^Linear solver <height> :\s*(\S+)", ji, re.M)
        used = m.group(1) if m else ""
        m = re.search(r"^OMP_MAX_THREADS :\s*(\d+)", ji, re.M)
        tused = m.group(1) if m else ""
    except OSError:
        pass
    try:
        with open(os.path.join(job, "Warnings.txt"), errors="replace", encoding="utf-8") as f:
            warns = sum(1 for l in f if l.strip())
    except OSError:
        pass
    if not keep:
        shutil.rmtree(rundir, ignore_errors=True)
    return {"case": case, "threads": threads, "threads_used": tused, "backend": backend or "auto", "used": used,
            "cpu": ru.ru_utime + ru.ru_stime, "wall": wall, "steps": steps, "end": end,
            "warnings": warns, "rc": rc}


def fmt_time(t):
    return "%.2f" % t if t < 0.95 else "%.1f" % t


def tex_file(case):
    return "\\file{%s}" % case.replace("_", "\\_")


def read_rows(path):
    """The rows of an existing runtimes.tsv as run_one's dicts, with their
    kind and repeat, and the yuzen line of its header."""
    rows, version = [], ""
    if not os.path.exists(path):
        return rows, version
    cols = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("# yuzen\t"):
                version = line.split("\t", 1)[1]
                continue
            if line.startswith("#") or not line:
                continue
            v = line.split("\t")
            if cols is None:
                cols = v
                continue
            if len(v) != len(cols):
                continue        # a row cut short by a job killed while writing it
            d = dict(zip(cols, v))
            rows.append({"kind": d["kind"], "repeat": int(d["repeat"]), "case": d["case"],
                         "threads": int(d["threads"]), "backend": d["backend"], "used": d["used"],
                         "cpu": float(d["cpu_s"]), "wall": float(d["wall_s"]), "steps": d["steps"],
                         "end": d["end"], "warnings": int(d["warnings"]), "rc": int(d["rc"]),
                         "threads_used": d.get("threads_used", "")})
    return rows, version


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", default=os.path.join(ROOT, "build", "bin", "yuzen"))
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=1,
                    help="single-threaded runs side by side, each on a core of its own (default 1)")
    ap.add_argument("--only", action="append", help="this case only (repeatable)")
    ap.add_argument("--backends", action="store_true", help="also measure tab:linsolve-timing")
    ap.add_argument("--threads", default="1 2", help='thread counts for --backends (default "1 2"; counts'
                    ' above the solver\'s limit for the case are refused by it and left out of the table)')
    ap.add_argument("--timeout", type=float, default=3600.0, help="seconds per run (default 3600)")
    ap.add_argument("--keep", action="store_true", help="keep the run directories")
    ap.add_argument("--backends-only", action="store_true", help="measure tab:linsolve-timing only")
    ap.add_argument("--resume", action="store_true",
                    help="keep the runs already in runtimes.tsv and append the missing ones")
    ap.add_argument("--budget", type=float, default=None,
                    help="seconds from the start after which no run is started or let finish;"
                         " exit code 3 if runs are left")
    args = ap.parse_args()

    t_start = time.perf_counter()
    exe = os.path.abspath(args.exe)
    os.makedirs(args.workdir, exist_ok=True)
    work = os.path.abspath(args.workdir)
    cores = sorted(os.sched_getaffinity(0))
    info = host_info(exe)
    if args.backends_only:
        args.backends = True
    cases = [] if args.backends_only else (args.only or all_cases())
    jobs = max(1, min(args.jobs, len(cores)))
    threads = [int(t) for t in args.threads.split()] if args.backends else []
    path = os.path.join(work, "runtimes.tsv")
    cols = ["kind", "case", "repeat", "threads", "backend", "used", "cpu_s", "wall_s", "steps",
            "end", "warnings", "rc", "threads_used"]
    prev, prev_version = read_rows(path) if args.resume else ([], "")
    if prev and prev_version and prev_version != info["yuzen"]:
        sys.exit("runtimes.tsv was written by %r, this binary is %r: use another --workdir"
                 % (prev_version, info["yuzen"]))
    if prev:
        tsv = open(path, "a", encoding="utf-8")
        tsv.write("# resumed\t%s on %s\n" % (info["date"], info["host"]))
    else:
        tsv = open(path, "w", encoding="utf-8")
        for k, v in info.items():
            tsv.write("# %s\t%s\n" % (k, v))
        tsv.write("# repeats\t%d\n# jobs\t%d\n" % (args.repeats, jobs))
        tsv.write("\t".join(cols) + "\n")
    tsv.flush()
    lock = threading.Lock()
    rows = list(prev)
    done = {(r["kind"], r["case"], r["repeat"], r["threads"], r["backend"]) for r in prev}
    deferred, too_long = [], []

    def record(kind, rep, r):
        with lock:
            tsv.write("\t".join(str(x) for x in [kind, r["case"], rep, r["threads"], r["backend"],
                                                 r["used"], "%.3f" % r["cpu"], "%.3f" % r["wall"],
                                                 r["steps"], r["end"], r["warnings"], r["rc"],
                                                 r["threads_used"]]) + "\n")
            tsv.flush()
            rows.append(dict(r, kind=kind, repeat=rep))
            print("%-9s %-24s rep %d  thr %-3d %-6s cpu %8.2f s  wall %8.2f s  %s  warnings %d"
                  % (kind, r["case"], rep, r["threads"], r["backend"], r["cpu"], r["wall"],
                     r["end"], r["warnings"]), flush=True)

    def measure(kind, case, rep, t, run_cores, backend, rundir):
        """One run, unless it is done already or the budget is spent."""
        key = (kind, case, rep, t, backend or "auto")
        if key in done:
            return
        timeout = args.timeout
        if args.budget is not None:
            left = args.budget - (time.perf_counter() - t_start)
            if left < 5.0:
                with lock:
                    deferred.append(key)
                return
            timeout = min(timeout, left)
        started = time.perf_counter() - t_start
        r = run_one(exe, case, rundir, t, run_cores, backend=backend, keep=args.keep, timeout=timeout)
        if timeout < args.timeout and r["wall"] >= timeout - 0.5 and r["end"].startswith("failed"):
            #  stopped by the budget, not failed: left for the next call --
            #  unless it had the whole budget, when no call would finish it
            with lock:
                deferred.append(key)
                if started < 5.0 * jobs:
                    too_long.append(key)
            return
        record(kind, rep, r)

    print("host %s, %s, %d cores allowed; %s" % (info["host"], info["cpu"], len(cores), info["yuzen"]))
    if prev:
        print("resuming: %d run(s) already in %s" % (len(prev), path))

    # --- the cases, single-threaded -------------------------------------
    free = list(cores[:jobs])
    flock = threading.Lock()

    def task(case, rep):
        with flock:
            core = free.pop()
        try:
            measure("case", case, rep, 1, [core], None,
                    os.path.join(work, "run_%s_%d" % (case.replace("/", "_"), rep)))
        finally:
            with flock:
                free.append(core)

    #  The long cases first, so that with --jobs they do not end up running
    #  last on their own.
    LONG = ["concz_2d", "leveling_2d_fine_reg", "leveling_2d_tight_reg", "leveling_2d_reg",
            "allphys_2d_reg", "eq/eq_n4000", "wet_static_5", "wet_coupled", "wet_static"]
    first = sorted(cases, key=lambda c: (LONG.index(c) if c in LONG else len(LONG)))
    order = [(c, rep) for rep in range(args.repeats) for c in first]
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        for f in [ex.submit(task, c, rep) for c, rep in order]:
            f.result()

    # --- the backends, at each thread count, one run at a time ------------
    for case in BACKEND_CASES:
        for t in threads:
            if t > len(cores):
                print("skipping %d threads: only %d cores allowed" % (t, len(cores)))
                continue
            for be in ("eigen", "amgcl"):
                for rep in range(args.repeats):
                    measure("backend", case, rep, t, cores[:t], be,
                            os.path.join(work, "run_be_%s_%s_%d_%d" % (case, be, t, rep)))
    tsv.close()

    if too_long:
        print("\n%s did not finish within --budget %g s, which it had to itself: raise --budget"
              % (" ".join("%s/%s/%d threads" % (k[1], k[4], k[3]) for k in too_long), args.budget))
        return 2
    if deferred:
        print("\nbudget spent: %d run(s) left for the next call (--resume); no tables written"
              % len(deferred))
        return 3

    # --- the tables, from every row of runtimes.tsv -----------------------
    res = {c: [r for r in rows if r["kind"] == "case" and r["case"] == c and r["repeat"] < args.repeats]
           for c in cases}
    lines = ["Cases, single-threaded, Backend auto: minimum CPU time over %d repeat(s)" % args.repeats,
             "%-24s %9s %9s %7s  %-14s %s" % ("case", "cpu [s]", "wall [s]", "steps", "end", "warnings")]
    tex = []
    for c in cases:
        rs = res[c]
        best = min(rs, key=lambda r: r["cpu"])
        same = all(r["end"] == best["end"] and r["warnings"] == best["warnings"] and
                   r["steps"] == best["steps"] for r in rs)
        lines.append("%-24s %9.2f %9.2f %7s  %-14s %d%s" % (c, best["cpu"], best["wall"], best["steps"],
                                                           best["end"], best["warnings"],
                                                           "" if same else "   (repeats differ!)"))
        tex.append("%-28s & %6s & %s & %4d \\\\" % (tex_file(c), fmt_time(best["cpu"]),
                                                     best["end"], best["warnings"]))
    if cases:
        with open(os.path.join(work, "table_cases.tex"), "w", encoding="utf-8") as f:
            f.write("\n".join(tex) + "\n")
    else:
        lines = []

    if args.backends:
        bres = {}
        refused = sorted({(r["case"], r["threads"]) for r in rows
                          if r["kind"] == "backend" and r["end"] == "refused"})
        for r in rows:
            if r["kind"] == "backend" and r["repeat"] < args.repeats and r["end"] != "refused":
                bres.setdefault((r["case"], r["threads"], r["backend"]), []).append(r)
        lines += ["" if lines else "", "Backends: minimum wall-clock time over %d repeat(s)" % args.repeats,
                  "%-18s %7s %10s %10s" % ("case", "threads", "eigen [s]", "amgcl [s]"),
                  "(threads: as asked for; in parentheses the count used, where the solver's"
                  " cap for small grids lowered it)"]
        if refused:
            lines.append("refused (more threads than the run can use): " +
                         ", ".join("%s at %d" % k for k in refused))
        tex = []
        for case in BACKEND_CASES:
            shown = [t for t in threads if (case, t, "eigen") in bres and (case, t, "amgcl") in bres]
            for t in shown:
                e = min(r["wall"] for r in bres[(case, t, "eigen")])
                a = min(r["wall"] for r in bres[(case, t, "amgcl")])
                used = {r["used"] for r in bres[(case, t, "amgcl")]}
                flag = "" if used <= {"amgcl", ""} else "   (amgcl not in this build: ran %s)" % ",".join(used)
                tu = {r["threads_used"] for r in bres[(case, t, "eigen")] + bres[(case, t, "amgcl")]} - {""}
                tshow = "%d" % t if tu <= {str(t)} else "%d (%s)" % (t, ",".join(sorted(tu)))
                lines.append("%-18s %7s %10.2f %10.2f%s" % (case, tshow, e, a, flag))
                be = "\\textbf{%s s}" % fmt_time(e) if e < a else "%s s" % fmt_time(e)
                ba = "\\textbf{%s s}" % fmt_time(a) if a < e else "%s s" % fmt_time(a)
                tex.append("%-40s & %d & %s & %s \\\\" % (tex_file(case) if t == shown[0] else "", t, be, ba))
        with open(os.path.join(work, "table_backends.tex"), "w", encoding="utf-8") as f:
            f.write("\n".join(tex) + "\n")

    hdr = ["host %s; %s; %d cores allowed" % (info["host"], info["cpu"], len(cores)),
           info.get("lscpu", ""), info["yuzen"], ""]
    with open(os.path.join(work, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(hdr + lines) + "\n")
    print("\n".join([""] + lines))
    bad = [c for c in cases if any(r["end"].startswith("failed") for r in res[c])]
    if bad:
        print("\nFAILED: " + " ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
