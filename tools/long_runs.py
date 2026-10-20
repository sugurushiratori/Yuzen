#!/usr/bin/env python3
"""Long runs and large grids: the balances over many steps, and the memory.

    python3 tools/long_runs.py --exe build/bin/yuzen --workdir runs/long
    python3 tools/long_runs.py --exe build/bin/yuzen --workdir runs/long --only wm1d_fixed

Each run is a shipped case with a few keys changed, run single-threaded:

  wm1d_fixed   allphys_1d_reg, dt fixed at 1e-4 s, to 50 s: 500,001 steps
  wm1d_long    allphys_1d_reg, adaptive dt, to 1000 s
  cz1d_fixed   concz_1d_reg, dt fixed at 1e-4 s, to 20 s: 200,001 steps
  wm2d_long    allphys_2d_reg, adaptive dt, to 30 s
  pvd_many     allphys_1d_reg on 20 cells, dt 1e-4 s, a field file every
               step, to 3 s: 30,001 field outputs
  wm2d_1024    scaling/bench2d_n1024 (1024^2 cells), its 4 steps
  wm1d_1e6     scaling/bench1d_n131072 at Nx = 10^6, 5 steps
  cz2d_200     concz_2d at 200 x 200 x 20, 3 steps

For each it reports the steps taken, the end time, the largest |VolErr| and
|ResinVolErr| of History.dat over the run, the peak resident memory, and,
on Linux, the resident memory a quarter of the way through and at the end
(sampled from /proc every 2 s), which do not differ unless memory grows
with the step count. The large grids need about 4 GB (wm2d_1024), 2 GB and
1.6 GB; the whole set takes about an hour on one core.

Writes long_runs.tsv and summary.txt in --workdir. The run directories are
deleted after each run unless --keep is given.

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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = os.path.join(ROOT, "cases")

FIXED = [(r"^InitialDt = .*$", "InitialDt = 1e-4"),
         (r"(\[Numerics\.TimeStepControl\]\n)Enable = true", r"\1Enable = false")]


def field_every(v):
    return (r"(\[Output\.FieldData\]\n)TimeInterval = .*", r"\g<1>TimeInterval = %s" % v)


def key(k, v):
    return (r"^%s = .*$" % k, "%s = %s" % (k, v))


RUNS = [
    ("wm1d_fixed", "allphys_1d_reg", FIXED + [field_every(5.0), key("EndTime", "50.0")]),
    ("wm1d_long", "allphys_1d_reg", [field_every(100.0), key("EndTime", "1000.0")]),
    ("cz1d_fixed", "concz_1d_reg", FIXED + [field_every(2.0), key("EndTime", "20.0")]),
    ("wm2d_long", "allphys_2d_reg", [field_every(3.0), key("EndTime", "30.0")]),
    ("pvd_many", "allphys_1d_reg", FIXED + [field_every(1e-4), key("Nx", "20"), key("EndTime", "3.0")]),
    ("wm2d_1024", "scaling/bench2d_n1024", []),
    ("wm1d_1e6", "scaling/bench1d_n131072", [key("Nx", "1000000"), key("MaxSteps", "5")]),
    ("cz2d_200", "concz_2d", [key("Nx", "200"), key("Ny", "200"), key("Nz", "20"), key("MaxSteps", "3")]),
]


def make_case(src, subs, job):
    s = open(os.path.join(CASES, src + ".toml")).read()
    for pat, rep in subs + [key("JobName", '"%s"' % job)]:
        s, n = re.subn(pat, rep, s, count=1, flags=re.M)
        if n != 1:
            sys.exit("%s: no line matches %s" % (src, pat))
    return s


def rss_kb(pid):
    try:
        with open("/proc/%d/status" % pid) as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except OSError:
        pass
    return None


def history(path):
    """Largest |VolErr|, |ResinVolErr| and the last time of a History.dat."""
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        names = [h.split(" ")[0] for h in head]
        iv, ir = names.index("VolErr"), names.index("ResinVolErr")
        ev = er = 0.0
        t = 0.0
        for line in f:
            p = line.split("\t")
            if len(p) < len(names):
                continue
            t = float(p[0])
            ev = max(ev, abs(float(p[iv])))
            er = max(er, abs(float(p[ir])))
    return ev, er, t


def jobinfo(path, label):
    with open(path) as f:
        for line in f:
            if line.startswith(label):
                return line.split(":", 1)[1].strip()
    return ""


def run_one(exe, workdir, name, src, subs, keep):
    d = os.path.join(workdir, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    with open(os.path.join(d, "case.toml"), "w") as f:
        f.write(make_case(src, subs, name))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as log:
        p = subprocess.Popen([exe, "case.toml"], cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
        samples = []
        while True:
            pid, st, usage = os.wait4(p.pid, os.WNOHANG)
            if pid == p.pid:
                break
            r = rss_kb(p.pid)
            if r is not None:
                samples.append(r)
            time.sleep(2.0)
        p.returncode = status = os.waitstatus_to_exitcode(st)
    wall = time.time() - t0
    #  ru_maxrss: kilobytes on Linux, bytes on macOS.
    peak = usage.ru_maxrss / (1024.0 if sys.platform == "darwin" else 1.0)
    res = os.path.join(d, name)
    ev, er, t = history(os.path.join(res, "History.dat"))
    steps = jobinfo(os.path.join(res, "JobInfo.txt"), "Time steps")
    nvts = len([x for x in os.listdir(res) if x.endswith(".vts")])
    q = samples[len(samples) // 4] if samples else None
    e = samples[-1] if samples else None
    if not keep:
        shutil.rmtree(d, ignore_errors=True)
    return dict(run=name, case=src, rc=status, steps=steps, t_end="%.6g" % t, vol_err="%.2e" % ev,
                resin_err="%.2e" % er, rss_peak_MB=fmt(peak), rss_q_MB=fmt(q), rss_end_MB=fmt(e),
                fields=nvts, wall_s="%.0f" % wall)


def fmt(kb):
    return "" if kb is None else "%.1f" % (kb / 1024.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", action="append", help="run only this one (repeatable)")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    exe = os.path.abspath(a.exe)
    os.makedirs(a.workdir, exist_ok=True)
    cols = ["run", "case", "rc", "steps", "t_end", "vol_err", "resin_err", "rss_peak_MB",
            "rss_q_MB", "rss_end_MB", "fields", "wall_s"]
    rows = []
    for name, src, subs in RUNS:
        if a.only and name not in a.only:
            continue
        r = run_one(exe, a.workdir, name, src, subs, a.keep)
        rows.append(r)
        print("\t".join(str(r[c]) for c in cols), flush=True)
    with open(os.path.join(a.workdir, "long_runs.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(a.workdir, "summary.txt"), "w") as f:
        w = [max(len(c), *(len(str(r[c])) for r in rows)) for c in cols]
        f.write("  ".join(c.ljust(n) for c, n in zip(cols, w)) + "\n")
        for r in rows:
            f.write("  ".join(str(r[c]).ljust(n) for c, n in zip(cols, w)) + "\n")
    bad = [r["run"] for r in rows if r["rc"] != 0]
    if bad:
        print("failed:", " ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
