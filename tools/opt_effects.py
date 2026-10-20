#!/usr/bin/env python3
"""The effect of each algorithm of YUZEN_OPT on its own.

    python3 tools/opt_effects.py --exe build-amg/bin/yuzen --workdir runs/opt
    python3 tools/opt_effects.py --exe build-amg/bin/yuzen --workdir runs/opt \\
            --threads "1 8 28 56 112" --repeats 2 --resume --budget 1500

Every case of the set below is run with each backend it lists, at each
thread count, once with every algorithm of util/Opt.h off ("base",
YUZEN_OPT=none: the algorithms they replaced), once with each on its own,
and once with all of them on, as a run is by default ("default"):

    regions    YUZEN_OPT=none,regions    fewer parallel regions, serial loops parallel
    assembly   YUZEN_OPT=none,assembly   face properties shared by A and its Jacobian
    ew         YUZEN_OPT=none,ew         Eisenstat-Walker linear tolerance (Newton only)
    concpc     YUZEN_OPT=none,concpc     light preconditioner of the concentration
    bilu       YUZEN_OPT=none,bilu       block-Jacobi ILUT (eigen only)
    amgcoarse  YUZEN_OPT=none,amgcoarse  coarsest level 500 instead of 3000 (amgcl only)
    numa       YUZEN_OPT=none,numa       first-touch placement of the large arrays
    interleave numactl --interleave=all with YUZEN_OPT=none (skipped without numactl)
    extrap     YUZEN_OPT=none,extrap     Newton start extrapolated (Newton only)
    default    YUZEN_OPT unset           all of them

(util/Opt.h has what each one does). Combinations are not run. Each run
has YUZEN_PROFILE=1, the thread cap off (YUZEN_THREAD_CAP=0), the threads
bound to cores 0..T-1 and OMP_WAIT_POLICY unset; the field output is
switched on (NeutralFile), so that the final state of a variant can be
compared with that of the base run at the same thread count: dH and dC are
the largest differences in the film height and the concentration over the
largest value of the base run. The time compared is that of the steps
(the profile's "step", which leaves the output out); the loop time of
JobInfo.txt is kept as well. With --repeats each run is repeated and the
fastest kept.

--resume and --budget work as in tools/case_runtimes.py (exit code 3: runs
left for the next call). The call that completes the set writes
summary.txt: for each case, backend and thread count, every variant's step
time, its speed-up over the base, the Newton and linear iterations and dH,
dC; then each variant's best time over the thread counts against the base's
best.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from case_runtimes import prepare, host_info   # noqa: E402

GENERIC = ["regions", "assembly", "ew", "concpc", "numa", "interleave", "extrap"]
AMG = ["amgcoarse"]
#  (case, steps or None for the case as shipped, Newton or not, {backend: variants})
CASES = [
    ("leveling_1d_reg", None, True, {"eigen": GENERIC}),
    ("scaling/bench1d_n08000", 100, True, {"eigen": GENERIC}),
    ("leveling_2d_reg", 100, True, {"amgcl": GENERIC + AMG, "eigen": ["bilu", "concpc"]}),
    ("scaling/bench2d_n0128", 30, True, {"amgcl": GENERIC + AMG, "eigen": ["bilu", "concpc"]}),
    ("scaling/bench2d_n0256", 10, True, {"amgcl": GENERIC + AMG, "eigen": ["bilu"]}),
    ("scaling/bench2d_n0512", 4, True, {"amgcl": GENERIC + AMG}),
    ("concz_2d", 20, False, {"amgcl": GENERIC + AMG, "eigen": ["bilu", "concpc"]}),
]
NEWTON_ONLY = {"ew", "extrap"}


def variants(newton, vs):
    return ["base"] + [v for v in vs if newton or v not in NEWTON_ONLY] + ["default"]


def set_case(path, steps):
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if steps is not None:
        text, n = re.subn(r"^MaxSteps = .*$", "MaxSteps = %d" % steps, text, count=1, flags=re.M)
        if n != 1:
            text = re.sub(r"^\[Numerics\]\n", "[Numerics]\nMaxSteps = %d\n" % steps, text, count=1, flags=re.M)
    #  the final state, for dH and dC (written after the loop)
    text = re.sub(r"^NeutralFile = false", "NeutralFile = true", text, flags=re.M)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def final_field(job):
    files = sorted(glob.glob(os.path.join(job, "field*.dat")))
    return files[-1] if files else None


def read_cols(path):
    """The HH and CC columns of a field .dat file."""
    with open(path, encoding="utf-8") as f:
        head = f.readline().split()
        if head and head[0] == "#":
            head = head[1:]         # the ConcZ files start their header with "# "
        idx = {n.split("[")[0]: i for i, n in enumerate(head)}
        cols = {k: [] for k in ("HH", "CC") if k in idx}
        for line in f:
            v = line.split()
            if len(v) != len(head):
                continue
            for k in cols:
                cols[k].append(float(v[idx[k]]))
    return cols


def rel_diff(a, b):
    if a is None or b is None or len(a) != len(b) or not a:
        return ""
    scale = max(abs(x) for x in a) or 1.0
    return "%.3e" % (max(abs(x - y) for x, y in zip(a, b)) / scale)


def run(exe, case, steps, backend, threads, cores, rundir, timeout, variant, numactl):
    name = prepare(case, rundir, backend)
    set_case(os.path.join(rundir, name), steps)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PROC_BIND="close", OMP_PLACES="cores",
               YUZEN_PROFILE="1", YUZEN_THREAD_CAP="0")
    for k in ("OMP_WAIT_POLICY", "YUZEN_OPT"):
        env.pop(k, None)
    cmd = [exe, name]
    if variant == "interleave":
        cmd = [numactl, "--interleave=all"] + cmd
        env["YUZEN_OPT"] = "none"
    elif variant == "base":
        env["YUZEN_OPT"] = "none"
    elif variant != "default":
        env["YUZEN_OPT"] = "none," + variant
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True, encoding="utf-8", errors="replace", env=env,
                           timeout=timeout, preexec_fn=lambda: os.sched_setaffinity(0, cores))
        rc, out = r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return None
    wall = time.perf_counter() - t0
    parts = {}
    for m in re.finditer(r"^\[profile\] (\S+)\s+(\d+)\s+([0-9.eE+-]+)$", out, re.M):
        parts[m.group(1)] = (int(m.group(2)), float(m.group(3)))
    job = os.path.join(rundir, os.path.basename(case))
    info = {}
    try:
        with open(os.path.join(job, "JobInfo.txt"), errors="replace", encoding="utf-8") as f:
            ji = f.read()
        for key, pat in (("loop", r"^Time loop \[s\] :\s*(\S+)"), ("steps", r"^Time steps :\s*(\S+)"),
                         ("nk", r"^Newton iterations :\s*(\S+)"), ("lin", r"^Linear iterations :\s*(\S+)"),
                         ("opt", r"^YUZEN_OPT :\s*(\S+)")):
            m = re.search(pat, ji, re.M)
            info[key] = m.group(1) if m else ""
    except OSError:
        pass
    ff = final_field(job)
    return {"rc": rc, "wall": wall, "parts": parts, "info": info,
            "field": read_cols(ff) if (ff and rc == 0) else None}


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
            done.add((d["case"], d["backend"], int(d["threads"]), d["variant"], int(d["repeat"])))
    return done, rows


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def summary(rows, threads):
    out = []
    best = {}
    for d in rows:
        if d["rc"] != "0" or fnum(d["step_s"]) is None:
            continue
        k = (d["case"], d["backend"], int(d["threads"]), d["variant"])
        if k not in best or float(d["step_s"]) < float(best[k]["step_s"]):
            best[k] = d
    pairs = []
    for case, steps, newton, bes in CASES:
        for be, vs in bes.items():
            pairs.append((case, be, variants(newton, vs)))
    for case, be, vs in pairs:
        ts = [t for t in threads if (case, be, t, "base") in best]
        if not ts:
            continue
        out.append("")
        out.append("== %s, %s" % (case, be))
        for t in ts:
            b = best[(case, be, t, "base")]
            out.append("  %d thread(s)" % t)
            out.append("    %-10s %9s %8s %7s %8s %10s %10s" % ("variant", "step [s]", "speed-up", "Newton",
                                                              "linear", "dH", "dC"))
            for v in vs:
                d = best.get((case, be, t, v))
                if d is None:
                    out.append("    %-10s %9s" % (v, "-"))
                    continue
                sp = float(b["step_s"]) / float(d["step_s"]) if float(d["step_s"]) > 0 else 0
                out.append("    %-10s %9.3f %8.3f %7s %8s %10s %10s"
                           % (v, float(d["step_s"]), sp, d["nk"], d["lin"], d["dH"] or "-", d["dC"] or "-"))
        out.append("  best over the thread counts (time, at threads; speed-up over the base's best)")
        bb = min((float(best[(case, be, t, "base")]["step_s"]), t) for t in ts)
        for v in vs:
            cand = [(float(best[(case, be, t, v)]["step_s"]), t) for t in ts if (case, be, t, v) in best]
            if not cand:
                continue
            m = min(cand)
            out.append("    %-10s %9.3f s at %3d   %6.3f" % (v, m[0], m[1], bb[0] / m[0]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", default="1 8 28 56 112")
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--only", action="append", help="this case only (repeatable)")
    ap.add_argument("--variant", action="append", help="this variant only (repeatable; base always runs)")
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--budget", type=float, default=None)
    args = ap.parse_args()

    t_start = time.perf_counter()
    exe = os.path.abspath(args.exe)
    work = os.path.abspath(args.workdir)
    os.makedirs(os.path.join(work, "ref"), exist_ok=True)
    cores = sorted(os.sched_getaffinity(0))
    threads = [int(t) for t in args.threads.split()]
    numactl = shutil.which("numactl")
    path = os.path.join(work, "effects.tsv")
    cols = ["case", "backend", "threads", "variant", "repeat", "opt", "rc", "wall_s", "loop_s", "step_s",
            "steps", "nk", "lin", "dH", "dC", "parts"]
    done, rows = read_done(path) if args.resume else (set(), [])
    info = host_info(exe)
    if rows:
        tsv = open(path, "a", encoding="utf-8")
        tsv.write("# resumed\t%s\n" % info["date"])
    else:
        tsv = open(path, "w", encoding="utf-8")
        for k, v in info.items():
            tsv.write("# %s\t%s\n" % (k, v))
        tsv.write("# numactl\t%s\n" % (numactl or "none"))
        tsv.write("\t".join(cols) + "\n")
    tsv.flush()
    print("host %s, %s, %d cores allowed; %s; numactl %s"
          % (info["host"], info["cpu"], len(cores), info["yuzen"], numactl or "none"), flush=True)

    left_over, too_long = 0, []
    for rep in range(args.repeats):
        for case, steps, newton, bes in CASES:
            if args.only and case not in args.only:
                continue
            for be, vs in bes.items():
                for t in threads:
                    if t > len(cores):
                        continue
                    for v in variants(newton, vs):
                        if args.variant and v != "base" and v not in args.variant:
                            continue
                        if v == "interleave" and not numactl:
                            continue
                        if (case, be, t, v, rep) in done:
                            continue
                        ref = os.path.join(work, "ref", "%s_%s_%d.dat" % (case.replace("/", "_"), be, t))
                        if v != "base" and not os.path.exists(ref):
                            continue        # the base of this thread count first (it runs before)
                        timeout = args.timeout
                        if args.budget is not None:
                            left = args.budget - (time.perf_counter() - t_start)
                            if left < 5.0:
                                left_over += 1
                                continue
                            timeout = min(timeout, left)
                        started = time.perf_counter() - t_start
                        rundir = os.path.join(work, "run_%s_%s_%d_%s" % (case.replace("/", "_"), be, t, v))
                        r = run(exe, case, steps, be, t, cores[:t], rundir, timeout, v, numactl)
                        dH = dC = ""
                        if r is not None and v == "base" and r["field"] is not None and not os.path.exists(ref):
                            ff = final_field(os.path.join(rundir, os.path.basename(case)))
                            shutil.copy(ff, ref)
                        if r is not None and r["field"] is not None and os.path.exists(ref):
                            rf = read_cols(ref)
                            dH = rel_diff(rf.get("HH"), r["field"].get("HH"))
                            dC = rel_diff(rf.get("CC"), r["field"].get("CC"))
                        shutil.rmtree(rundir, ignore_errors=True)
                        if r is None:
                            left_over += 1
                            if args.budget is not None and started < 5.0:
                                too_long.append((case, be, t, v))
                            continue
                        step = r["parts"].get("step", (0, None))[1]
                        parts = ";".join("%s=%d,%.6f" % (k, c, s) for k, (c, s) in r["parts"].items())
                        row = {"case": case, "backend": be, "threads": str(t), "variant": v, "repeat": str(rep),
                               "opt": r["info"].get("opt", ""), "rc": str(r["rc"]), "wall_s": "%.3f" % r["wall"],
                               "loop_s": r["info"].get("loop", ""),
                               "step_s": "" if step is None else "%.6f" % step,
                               "steps": r["info"].get("steps", ""), "nk": r["info"].get("nk", ""),
                               "lin": r["info"].get("lin", ""), "dH": dH, "dC": dC, "parts": parts}
                        tsv.write("\t".join(row[c] for c in cols) + "\n")
                        tsv.flush()
                        rows.append(row)
                        print("%-24s %-6s thr %-3d %-10s rep %d step %9s s  nk %6s lin %7s  dH %10s dC %10s rc %d"
                              % (case, be, t, v, rep, row["step_s"], row["nk"], row["lin"], dH or "-",
                                 dC or "-", r["rc"]), flush=True)
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
