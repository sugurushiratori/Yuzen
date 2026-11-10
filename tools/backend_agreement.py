#!/usr/bin/env python3
"""Run representative cases with each linear-solver backend and compare.

    python3 tools/backend_agreement.py --exe build/bin/yuzen --workdir /tmp/ba

The same input is run twice, once with Backend = "eigen" and once with
Backend = "amgcl" for both the film-height and the concentration equation,
on one thread, and the two results are compared: every physical column of
the last field file and of History.dat, relative to the largest value of
the column, must agree to TOL; the conservation errors (VolErr,
ResinVolErr), which are rounding, must both stay below ABS_BALANCE. The two
backends stop their iterations at different points inside the same
tolerance, so the answers differ by about that tolerance and no more; a
backend that solved a different system, dropped a boundary row or lost
solute would show far more.

The cases cover the well-mixed model in one and two dimensions with every
term but the disjoining pressure (allphys_*), and the CZ model in one and
two dimensions. Needs a build with AMGCL (-DYUZEN_WITH_AMGCL=ON); ctest
registers it as backend_agreement only then.

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

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(os.path.dirname(HERE), "cases")

#  (case, EndTime or None for the case's own); concz_2d is cut short, it is
#  the expensive one
PLAN = [("leveling_1d_reg", None),
        ("allphys_1d_reg", None),
        ("concz_1d_reg", None),
        ("allphys_2d_reg", None),
        ("concz_2d", "5e-4")]

TOL = 1e-7              # physical columns, relative to the column's scale
ABS_BALANCE = 1e-11     # VolErr, ResinVolErr: absolute, both runs
BALANCE = {"VolErr", "ResinVolErr"}
#  How the solve went, not what it found: reported, not judged
DIAG = {"NK_iter", "NK_resid", "BCG_iter", "BCG_resid"}


def with_backend(text, backend, end_time):
    """The case with Backend set in both LinearSolver blocks (added to each
    block, or the block created) and, if given, a shorter EndTime."""
    for blk in ("Numerics.HeightEquation.LinearSolver", "Numerics.ConcEquation.LinearSolver"):
        head = "[%s]" % blk
        if re.search(r"^%s\s*$" % re.escape(head), text, flags=re.M):
            text = re.sub(r"^(%s\s*\n)(?:Backend = [^\n]*\n)?" % re.escape(head),
                          lambda m: m.group(1) + 'Backend = "%s"\n' % backend,
                          text, count=1, flags=re.M)
        else:
            text = text.rstrip("\n") + '\n\n%s\nBackend = "%s"\n' % (head, backend)
    if end_time:
        text, n = re.subn(r"^EndTime = [^\n]*$", "EndTime = " + end_time, text, count=1, flags=re.M)
        assert n == 1
    return text


def run(exe, d, name, text):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name + ".toml"), "w", encoding="utf-8") as f:
        f.write(text)
    for dat in set(re.findall(r'"([^"]+\.dat)"', text)):
        if os.path.exists(os.path.join(CASES, dat)):
            shutil.copy(os.path.join(CASES, dat), d)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    r = subprocess.run([exe, name + ".toml"], cwd=d, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=env, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError("%s in %s ended rc=%d:\n%s" % (name, d, r.returncode, (r.stdout + r.stderr)[-1500:]))
    jobs = [p for p in glob.glob(os.path.join(d, "*")) if os.path.isfile(os.path.join(p, "History.dat"))]
    assert len(jobs) == 1, jobs
    return jobs[0]


def table(path):
    with open(path, encoding="utf-8") as f:
        head = [h.split("[")[0].strip() for h in f.readline().lstrip("#").rstrip("\n").split("\t")]
        rows = [[float(x) for x in line.split()] for line in f if line.strip()]
    return head, rows


def compare(pa, pb):
    """{column: (difference over scale, largest |value| of either run)}"""
    ha, ra = table(pa)
    hb, rb = table(pb)
    if ha != hb or len(ra) != len(rb):
        return None
    out = {}
    for j, name in enumerate(ha):
        if not name:
            continue
        a = [r[j] for r in ra]
        b = [r[j] for r in rb]
        scale = max(max(abs(x) for x in a), max(abs(x) for x in b))
        d = max(abs(x - y) for x, y in zip(a, b))
        out[name] = (d / scale if scale > 0 else 0.0, scale)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", action="append", help="this case only (repeatable)")
    args = ap.parse_args()
    args.exe = os.path.abspath(args.exe)
    bad = 0
    for case, end in PLAN:
        if args.only and case not in args.only:
            continue
        src = open(os.path.join(CASES, case + ".toml"), encoding="utf-8").read()
        res = {}
        try:
            for be in ("eigen", "amgcl"):
                d = os.path.join(args.workdir, case, be)
                shutil.rmtree(d, ignore_errors=True)
                res[be] = run(args.exe, d, case, with_backend(src, be, end))
        except RuntimeError as e:
            print("FAIL %-16s %s" % (case, e))
            bad += 1
            continue
        problems, worst, worst_col = [], 0.0, ""
        fa = sorted(glob.glob(os.path.join(res["eigen"], "field*.dat")))
        fb = sorted(glob.glob(os.path.join(res["amgcl"], "field*.dat")))
        if [os.path.basename(p) for p in fa] != [os.path.basename(p) for p in fb]:
            problems.append("different field files written")
        for label, pa, pb in (("field", fa[-1], fb[-1]),
                              ("History", os.path.join(res["eigen"], "History.dat"),
                               os.path.join(res["amgcl"], "History.dat"))):
            c = compare(pa, pb)
            if c is None:
                problems.append("%s: layout or row count differs" % label)
                continue
            for col, (rel, scale) in c.items():
                if col in DIAG or col in ("Time", "time", "Step"):
                    continue
                if col in BALANCE:
                    if scale > ABS_BALANCE:
                        problems.append("%s %s reaches %.2e" % (label, col, scale))
                    continue
                if rel > worst:
                    worst, worst_col = rel, "%s %s" % (label, col)
                if rel > TOL:
                    problems.append("%s %s differs by %.2e of its scale" % (label, col, rel))
        if problems:
            print("FAIL %-16s %s" % (case, "; ".join(problems)))
            bad += 1
        else:
            print("ok   %-16s largest difference %.1e (%s)" % (case, worst, worst_col))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
