#!/usr/bin/env python3
"""The well-mixed and the CZ model on the same evaporating film, at three
Peclet numbers.

    python3 tools/wm_cz_evaporation.py --exe build/bin/yuzen --workdir runs/wmcz --run

Six runs, all made from cases/concz_1d_reg.toml: a one-
dimensional film on a substrate step, 5 um thick on the low side and 10 um
on the high side (top flat), resin fraction 0.5 at the start, evaporating at
E0 = 1 um/s (Raoult), capillary leveling on, solutal Marangoni off (with it, the
stratified film is unstable in the model at Pe >= 1; see the manual,
"What has not been verified"), for 1.5 s. The diffusivity D is constant and takes three
values, D = 1e-10, 1e-11, 1e-12 m^2/s, so that the Peclet number
Pe = E0 h0 / D, h0 = 10 um, is 0.1, 1 and 10. Each D is run with the CZ
model and with the WM model on the same grid (Nx = 120, Nz = 121).

With --run the cases are written to --workdir and run there (one thread
each, two at a time; about ten minutes on two cores); without it only the
extraction is done, from runs already in --workdir. A summary is printed
(evaporated volume and surface concentration at 1.5 s, and the ratio of the
evaporated volumes); with --out the time series and the profiles are
written there as wmcz_*.dat. Needs NumPy; not part of ctest.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(os.path.dirname(HERE), "cases")
DATA = None                                        # --out
PE = [("pe01", 1e-10, "0.1"), ("pe1", 1e-11, "1"), ("pe10", 1e-12, "10")]
E0 = 1e-6


def sub(s, pat, rep):
    s2, n = re.subn(pat, rep, s, count=1, flags=re.M)
    assert n == 1, pat
    return s2


def make_case(base, tag, D, model):
    s = base
    s = sub(s, r'^JobName = .*$', 'JobName = "%s_%s"' % (tag, model.lower()))
    s = sub(s, r'^Model = "CZ"$', 'Model = "%s"' % model)
    s = sub(s, r'^Nx = 60$', 'Nx = 120')
    s = sub(s, r'^Nz = 11$', 'Nz = 121')
    s = sub(s, r'^PrecursorRatio = 0\.5$', 'PrecursorRatio = 0.001')
    i = s.index('[Initial.SurfaceShape.EllipseBump]')
    j = s.index('[Substrate]')
    s = s[:i] + s[i:j].replace('Base = 0.0', 'Base = 0.5', 1) + s[j:]
    s = sub(s, r'^(\[Physics\.SolutoCapillarity\]\n)Enable = true', r'\g<1>Enable = false')
    s = sub(s, r'^(\[Physics\.Evaporation\]\nEnable = true\n)Coef = 1e-7', r'\g<1>Coef = %g' % E0)
    s = sub(s, r'^(\[Fluid\.Diffusivity\]\n)BaseValue = 1\.0e-9\nBaseConc = 0\.5\nExponent = -7\.0',
            r'\g<1>BaseValue = %g\nBaseConc = 0.5\nExponent = 0.0' % D)
    s = sub(s, r'^EndTime = 0\.5', 'EndTime = 1.5')
    s = sub(s, r'^(\[Output\.FieldData\]\n)TimeInterval = 1e-1', r'\g<1>TimeInterval = 0.05')
    return s


def run_all(exe, wd):
    with open(os.path.join(CASES, "concz_1d_reg.toml"), encoding="utf-8") as f:
        base = f.read()
    os.makedirs(wd, exist_ok=True)
    names = []
    for tag, D, _ in PE:
        for model in ("CZ", "WM"):
            name = "%s_%s" % (tag, model.lower())
            with open(os.path.join(wd, name + ".toml"), "w", encoding="utf-8") as f:
                f.write(make_case(base, tag, D, model))
            names.append(name)
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def one(name):
        with open(os.path.join(wd, name + ".log"), "w", encoding="utf-8") as log:
            r = subprocess.run([exe, name + ".toml"], cwd=wd, stdout=log, stderr=subprocess.STDOUT, env=env)
        return name, r.returncode
    with ThreadPoolExecutor(2) as ex:
        for name, rc in ex.map(one, names):
            print("%-8s rc=%d" % (name, rc))
            assert rc == 0, name


def rd(p):
    with open(p, encoding="utf-8") as f:
        head = [c.split("[")[0].strip() for c in f.readline().lstrip("#").strip().split("\t") if c.strip()]
    return dict(zip(head, np.loadtxt(p, skiprows=1).T))


def field(path):
    """x [mm], h [um], surface c, depth-averaged c, substrate c, surface
    height and substrate height [um] of one field file."""
    d = rd(path)
    if "zz" in d:                                  # CZ: x fastest, then zeta levels
        nx = len(np.unique(d["xx"]))
        nz = len(d["CC"]) // nx
        C = d["CC"].reshape(nz, nx)
        cbar = 0.5 * (C[1:] + C[:-1]).mean(axis=0)   # levels equally spaced in zeta
        return d["xx"][:nx], d["HH"][:nx], C[-1], cbar, C[0], d["HS"][:nx], d["SS"][:nx]
    return d["xx"], d["HH"], d["CC"], d["CC"], d["CC"], d["HS"], d["SS"]


def write(name, head, rows):
    if DATA is None:
        return
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, name), "w", encoding="utf-8") as f:
        for h in head:
            f.write("# %s\n" % h)
        for r in rows:
            f.write("\t".join("%.8g" % v for v in r) + "\n")


def extract(wd):
    flat = 0.3                                     # x > 0.3 mm: the flat 10 um side
    out_t, summary = {}, []
    for tag, D, pe in PE:
        for m in ("cz", "wm"):
            job = os.path.join(wd, "%s_%s" % (tag, m))
            hist = rd(os.path.join(job, "History.dat"))
            t, ev = hist["Time"], hist["EvaporatedVol"]
            rate = np.gradient(ev, t) / 1e-3 * 1e6      # m^3/m/s over the 1 mm domain -> um/s
            fs = sorted(glob.glob(os.path.join(job, "field*.dat")))
            ts, cs = [], []
            for k, p in enumerate(fs):
                x, h, cs_, cb, c0, _, _ = field(p)
                ts.append(0.05 * k)
                cs.append(cs_[x > flat].mean())
            out_t[(tag, m)] = (t, rate, np.array(ts), np.array(cs))
            x, h, cs_, cb, c0, hs, ss = field(fs[30])    # t = 1.5 s
            write("wmcz_%s_%s_x.dat" % (tag, m),
                  ["%s %s, t = 1.5 s, Pe = %s" % (tag, m.upper(), pe),
                   "x[mm] h[um] c_surface c_depthavg c_substrate surface[um] substrate[um]"],
                  zip(x, h, cs_, cb, c0, hs, ss))
            write("wmcz_%s_%s_t.dat" % (tag, m),
                  ["%s %s, Pe = %s: evaporation rate averaged over the 1 mm domain" % (tag, m.upper(), pe),
                   "t[s] rate[um/s]"], zip(t, rate))
            write("wmcz_%s_%s_c.dat" % (tag, m),
                  ["%s %s, Pe = %s: surface concentration over x > %g mm (WM: its only value)" % (tag, m.upper(), pe, flat),
                   "t[s] c_surface"], zip(ts, cs))
            summary.append((tag, m, ev[-1], cs[-1]))
            if (tag, m) == ("pe1", "wm"):
                x, h, _, _, _, hs, ss = field(fs[0])
                write("wmcz_init_x.dat", ["t = 0 (every run)", "x[mm] h[um] surface[um] substrate[um]"], zip(x, h, hs, ss))
    print("tag  model  evaporated[m^3/m]  c_surface(1.5 s)")
    for tag, m, ev, c in summary:
        print("%-5s %-3s %.5e %.4f" % (tag, m, ev, c))
    for tag, _, pe in PE:
        e = {m: ev for t_, m, ev, c in summary if t_ == tag}
        print("Pe = %-3s evaporated CZ/WM = %.3f" % (pe, e["cz"] / e["wm"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", default=os.path.join(os.path.dirname(HERE), "build", "bin", "yuzen"))
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--run", action="store_true", help="write and run the six cases first")
    ap.add_argument("--out", help="write the data files here")
    a = ap.parse_args()
    global DATA
    DATA = a.out and os.path.abspath(a.out)
    wd = os.path.abspath(a.workdir)
    if a.run:
        run_all(os.path.abspath(a.exe), wd)
    extract(wd)


if __name__ == "__main__":
    main()
