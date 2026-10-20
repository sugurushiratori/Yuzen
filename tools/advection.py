#!/usr/bin/env python3
"""The advection of the concentration by each convection scheme.

    python3 tools/advection.py --exe build/bin/yuzen --workdir /tmp/adv
    python3 tools/advection.py --exe build/bin/yuzen --workdir /tmp/adv --all
    python3 tools/advection.py --exe build/bin/yuzen --workdir /tmp/adv --oblique

A flat film of uniform thickness, periodic in x, is carried along by
in-plane gravity at the uniform velocity rho g h^2 / 3mu (the viscosity does
not depend on the concentration, and nothing but the flow moves the solute:
the diffusivity is 1e-13 m^2/s, a cell Peclet number of 1e6, and there is no
Marangoni flux and no evaporation). Its answer is the initial concentration
translated, and two profiles are put through it with the scheme of
[Numerics.Convection] Concentration:

    front_wm   a step, c = 0.7 on the middle half of the domain and 0.3
               elsewhere, over one period at a Courant number of 0.4: the
               extremes of c over the run against the bounds 0.3 and 0.7.
               The limited schemes (Koren, VanLeer) must stay inside them to
               rounding; for the others the over- and undershoot is reported.
    front_cz   the same with ConcZ (4 cells through the film; each level is
               carried at its own speed, 1.5 times the mean at the surface).
    smooth_wm  c = 0.5 + 0.2 sin(2 pi x / L) over one period on 32, 64 and
               128 cells at a Courant number of 0.05 (the time error of
               Crank-Nicolson is then far below the grid's): the largest
               error against the translated profile, and its order in dx.

    courant    the Courant check of the limited schemes: front_wm and
               front_cz with Koren over a tenth of a period at several
               Courant numbers u dt / dx, and the warnings each raises in
               Warnings.txt. The concentration is warned about above 1 (WM:
               1.5, not 0.4; CZ, whose surface moves at 1.5 u: 1.5, not
               0.4), the film thickness, whose fluxes move at dQ/dh = 3 u,
               above a third (0.4 warns, 0.3 does not).

    oblique_wm the asymptotic order of the well-mixed concentration in 2-D,
               with diffusion: a square periodic film carried by gravity in
               the direction (2, 1)/sqrt(5), a sine along the diagonal,
               c = 0.5 + 0.2 sin(k (x + y)), oblique to the flow and to both
               axes, so that neither the field nor its motion is separable,
               and D = 5e-7 m^2/s: c = 0.5 + 0.2 exp(-D |k|^2 t)
               sin(k.(x - u t)). 16, 32 and 64 cells a side, with 20, 40 and
               80 steps over one decay time (0.025 s), dt refined with dx;
               the largest error and its order on each pair of grids (see
               OBLIQUE for the expected values).

The film is made by a first run of one step, whose field00000.dat is
edited to carry the concentration and read back through Restart.

By default (as ctest runs it) the fixed schemes KK, UTOPIA and Central4 and
the limited ones Koren and VanLeer are run, and what is checked is: the
limited schemes bounded, every scheme's smooth order at least that of its
kind in the max norm (3 for KK and UTOPIA, 2 for Central4; the limited
ones clip the extrema of the sine to first order there, and converge in the
max norm at 1.56 (Koren) and 1.41 (VanLeer), checked at 1.4 and 1.3; in the
mean of |error| the order is shown as well), and for oblique_wm, with the
default scheme UTOPIA and the limited Koren, the order on the finer pair of
grids between 1.9 and 2.1 and that on the coarser pair within 0.1 of it.
--all adds QUICK, Upwind1 and Upwind2, reported only.

Only the Python standard library is used, so the test runs wherever ctest
does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import math
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from concz_vs_wm import edit_toml               # noqa: E402
from equivalence import arr, fixed_steps, run, read_field   # noqa: E402
from analytic import case, quiet, wave_film     # noqa: E402

L = 1.0e-3                        # the domain, m (periodic)
H0 = 5.0e-6                       # the film, m
RHO, MU = 1000.0, 0.1
U = 1.0e-2                        # the velocity, m/s
G = 3.0 * MU * U / (RHO * H0 ** 2)
PERIOD = L / U                    # 0.1 s

#  (scheme, minimum smooth order in the max norm) -- None: reported only
CHECKED = [("KK", 2.8), ("UTOPIA", 2.8), ("Central4", 1.8), ("Koren", 1.4), ("VanLeer", 1.3)]
REPORTED = [("QUICK", None), ("Upwind1", None), ("Upwind2", None)]
LIMITED = ("Koren", "VanLeer")


def film(model, nx, dt, steps, interval, scheme):
    s = quiet(case("periodic_1d_reg"))
    s = wave_film(s, H0 / 10e-6, 0.0, L)
    s = edit_toml(s, {}, {
        "Grid": {"Nx": str(nx), "Xmin": repr(-L / 2), "Xmax": repr(L / 2),
                 "Nz": "4" if model == "CZ" else "1"},
        "Physics.ConcField": {"Model": '"%s"' % model, "EnableConvection": "true", "EnableDiffusion": "true"},
        "Physics.LaplacePressure": {"Enable": "true"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(G, 0.0)},
        "Fluid.Density": {"Value": repr(RHO)},
        "Fluid.Viscosity": {"BaseValue": repr(MU), "Exponent": "0.0"},
        "Fluid.Diffusivity": {"BaseValue": "1.0e-13", "Exponent": "0.0"},
        "Numerics.Convection": {"Concentration": '"%s"' % scheme},
    })
    return fixed_steps(s, dt, steps, interval)


def write_restart(src_dat, dst_dat, cfun):
    """Copy a field*.dat with its concentration column replaced by cfun(x)
    (in 2-D by cfun(x, y))."""
    with open(src_dat) as f:
        lines = f.read().split("\n")
    head = lines[0].lstrip("#").split()
    names = [c.split("[")[0] for c in head]
    ix, ic = names.index("xx"), names.index("CC")
    iy = names.index("yy") if "yy" in names else None
    out = [lines[0]]
    for ln in lines[1:]:
        v = ln.split("\t") if "\t" in ln else ln.split()
        if len([t for t in v if t.strip()]) != len(names):
            out.append(ln)
            continue
        v = [t for t in v if t.strip()]
        if iy is None:
            v[ic] = "%.17e" % cfun(float(v[ix]) * 1e-3)
        else:
            v[ic] = "%.17e" % cfun(float(v[ix]) * 1e-3, float(v[iy]) * 1e-3)
        out.append("\t".join(v))
    with open(dst_dat, "w") as f:
        f.write("\n".join(out))


def advect(workdir, exe, tag, model, scheme, nx, courant, periods, cfun):
    """Run the film from the concentration cfun(x) for `periods` periods;
    returns the list of (time, field) written every tenth of a period."""
    dx = L / nx
    dt = courant * dx / U
    steps = int(round(periods * PERIOD / dt))
    dt = periods * PERIOD / steps
    per_out = max(1, steps // 10)
    interval = per_out * dt
    d = os.path.join(workdir, "%s_%s_%s_n%d" % (tag, model, scheme, nx))
    #  the film, from a run of one step
    rc, log = run(exe, os.path.join(d, "init"), "a", film(model, nx, dt, 1, dt, scheme), 300)
    if rc != 0:
        raise RuntimeError("initial run failed:\n" + log[-600:])
    write_restart(os.path.join(d, "init", "a", "field00000.dat"), os.path.join(d, "start.dat"), cfun)
    s = film(model, nx, dt, steps, interval, scheme)
    s = edit_toml(s, {}, {"Restart": {"Enable": "true", "ASCII": "true", "FileName": '"../start.dat"',
                                      "RestartTime": "0", "FileIndex": "0"}})
    rc, log = run(exe, os.path.join(d, "run"), "a", s, 900)
    if rc != 0:
        raise RuntimeError("run failed:\n" + log[-600:])
    out = []
    for p in sorted(glob.glob(os.path.join(d, "run", "a", "field*.dat"))):
        n = int(re.search(r"(\d+)\.dat$", p).group(1))
        out.append((n * interval, read_field(p)))
    return out, steps


def interior(f, nx_, ny_, nz_):
    """Rows of the real cells: in x the cell centers (1 .. nx_-2), and for
    ConcZ the levels 1 .. nz_-2 of each."""
    rows = []
    for k in (range(1, nz_ - 1) if nz_ > 1 else [0]):
        for i in range(1, nx_ - 1):
            rows.append(k * nx_ * ny_ + i)
    return rows


def step_c(x):
    return 0.7 if abs(x) < L / 4 else 0.3


def smooth_c(x):
    return 0.5 + 0.2 * math.sin(2 * math.pi * x / L)


def front(workdir, exe, model, scheme):
    res, _ = advect(workdir, exe, "front", model, scheme, 100, 0.4, 1.0, step_c)
    cmin, cmax = 1.0, 0.0
    for t, (f, nx_, ny_, nz_) in res:
        for r in interior(f, nx_, ny_, nz_):
            cmin = min(cmin, f["CC"][r])
            cmax = max(cmax, f["CC"][r])
    return cmin, cmax


def smooth(workdir, exe, scheme, grids=(32, 64, 128)):
    errs = []
    for nx in grids:
        res, _ = advect(workdir, exe, "smooth", "WM", scheme, nx, 0.05, 1.0, smooth_c)
        t, (f, nx_, ny_, nz_) = res[-1]
        e, e1, n = 0.0, 0.0, 0
        for r in interior(f, nx_, ny_, nz_):
            x = f["xx"][r] * 1e-3
            #  one period: the profile is back where it started
            d = abs(f["CC"][r] - smooth_c(x - U * t))
            e = max(e, d)
            e1 += d
            n += 1
        errs.append((e / 0.2, e1 / n / 0.2))
    return errs


#  oblique_wm: the asymptotic order of the well-mixed concentration in 2-D.
#  With the film flat and the viscosity and the diffusivity constant, the
#  concentration equation is c_t + u.grad c = D lap c, and the error of the
#  scheme follows from its modified equation. The diffusion, a central
#  difference in each direction, under-damps the wave by
#  D (kx^4 + ky^4) dx^2/12; UTOPIA (and Koren where it does not limit)
#  over-damps it by (|ux| kx^4 + |uy| ky^4) dx^3/12; Crank-Nicolson, with dt
#  refined with dx, adds an error in dt^2, here about 1 % of the others. So
#
#      e / a = T exp(-D |k|^2 T) [D (kx^4 + ky^4) dx^2/12
#                                 - (|ux| kx^4 + |uy| ky^4) dx^3/12]:
#
#  second order, with a third-order part of the other sign whose share, for
#  kx = ky, is r = (|ux| + |uy|) dx / 2D, 0.67 times the cell Peclet number
#  u dx / D. The observed order is 2 + log2[(1 - r) / (1 - r/2)], below 2
#  until r is small -- the cancellation that keeps the order of c low in
#  wm_all_2d of tools/mms.py on the grids 2-D can afford, and the reason
#  that case takes a tenth of the diffusivity of 1-D. Here the cell Peclet
#  number is 0.125, 0.0625 and 0.031 on 16, 32 and 64 cells (r = 0.084,
#  0.042, 0.021), and the formula gives orders of 1.94 and 1.97. UTOPIA
#  converges at 1.91 and 1.97, its error within 3 % of the formula's; Koren,
#  which clips the extrema of the sine to first order, at 1.88 and 1.92.
#  Each must converge at 2 on the finer pair (OBLIQUE_BAND), with the order
#  on the coarser pair within OBLIQUE_SETTLED of it, i.e. in its asymptotic
#  range. The convection is not lost behind the diffusion: with first-order
#  upwinding the orders are 0.72 and 0.87, with KK (three times UTOPIA's
#  dissipation) 1.76 and 1.90. (The film-height equation is not checked
#  here -- the film stays flat to rounding; its order in 2-D is that of the
#  manufactured solutions.)
OBL_L = 1.0e-3                    # the square, m (periodic)
OBL_U = 1.0e-3                    # the velocity, m/s
OBL_DIR = (2.0 / math.sqrt(5.0), 1.0 / math.sqrt(5.0))
OBL_D = 5.0e-7                    # m^2/s
OBL_A = 0.2
OBL_T = 0.025                     # s; exp(-D |k|^2 T) = 0.37
OBL_GRIDS = ((16, 20), (32, 40), (64, 80))   # (cells a side, steps)
OBLIQUE = ("UTOPIA", "Koren")
OBLIQUE_BAND = (1.9, 2.1)
OBLIQUE_SETTLED = 0.1


def oblique_c(x, y, t):
    k = 2 * math.pi / OBL_L
    ux, uy = OBL_U * OBL_DIR[0], OBL_U * OBL_DIR[1]
    return 0.5 + OBL_A * math.exp(-2 * OBL_D * k * k * t) * math.sin(k * (x - ux * t) + k * (y - uy * t))


def oblique_predicted(n):
    """The leading error of the modified equation on n cells a side,
    relative to a (above)."""
    k, dx = 2 * math.pi / OBL_L, OBL_L / n
    rate = (OBL_D * 2 * k ** 4 * dx ** 2 / 12
            - OBL_U * (OBL_DIR[0] + OBL_DIR[1]) * k ** 4 * dx ** 3 / 12)
    return OBL_T * math.exp(-2 * OBL_D * k * k * OBL_T) * rate


def film_oblique(n, dt, steps, interval, scheme):
    g = 3.0 * MU * OBL_U / (RHO * H0 ** 2)
    s = quiet(case("periodic_1d_reg"), twoD=True)
    s = wave_film(s, H0 / 10e-6, 0.0, OBL_L, OBL_L)
    s = edit_toml(s, {}, {
        "Grid": {"Dimension": '"xyh"', "Nx": str(n), "Ny": str(n), "Nz": "1",
                 "Xmin": repr(-OBL_L / 2), "Xmax": repr(OBL_L / 2),
                 "Ymin": repr(-OBL_L / 2), "Ymax": repr(OBL_L / 2)},
        "Physics.ConcField": {"Model": '"WM"', "EnableConvection": "true", "EnableDiffusion": "true"},
        #  (the film stays flat, and without the capillary pressure its
        #  equation is cheaper to solve)
        "Physics.LaplacePressure": {"Enable": "false"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(g * OBL_DIR[0], g * OBL_DIR[1])},
        "Fluid.Density": {"Value": repr(RHO)},
        "Fluid.Viscosity": {"BaseValue": repr(MU), "Exponent": "0.0"},
        "Fluid.Diffusivity": {"BaseValue": repr(OBL_D), "Exponent": "0.0"},
        "Numerics.Convection": {"Concentration": '"%s"' % scheme},
    })
    return fixed_steps(s, dt, steps, interval)


def oblique(workdir, exe, scheme):
    """The largest error of c (relative to a) and of h (relative to the
    film) at the end, on each grid of OBL_GRIDS."""
    out = []
    for n, steps in OBL_GRIDS:
        dt = OBL_T / steps
        d = os.path.join(workdir, "oblique_WM_%s_n%d" % (scheme, n))
        rc, log = run(exe, os.path.join(d, "init"), "a", film_oblique(n, dt, 1, dt, scheme), 300)
        if rc != 0:
            raise RuntimeError("initial run failed:\n" + log[-600:])
        write_restart(os.path.join(d, "init", "a", "field00000.dat"), os.path.join(d, "start.dat"),
                      lambda x, y: oblique_c(x, y, 0.0))
        s = film_oblique(n, dt, steps, OBL_T, scheme)
        s = edit_toml(s, {}, {"Restart": {"Enable": "true", "ASCII": "true", "FileName": '"../start.dat"',
                                          "RestartTime": "0", "FileIndex": "0"}})
        rc, log = run(exe, os.path.join(d, "run"), "a", s, 900)
        if rc != 0:
            raise RuntimeError("run failed:\n" + log[-600:])
        f, nx_, ny_, _ = read_field(sorted(glob.glob(os.path.join(d, "run", "a", "field*.dat")))[-1])
        ec = eh = 0.0
        for j in range(1, ny_ - 1):
            for i in range(1, nx_ - 1):
                r = j * nx_ + i
                ec = max(ec, abs(f["CC"][r] - oblique_c(f["xx"][r] * 1e-3, f["yy"][r] * 1e-3, OBL_T)))
                eh = max(eh, abs(f["HH"][r] * 1e-6 - H0))
        out.append((ec / OBL_A, eh / H0))
    return out


#  (model, Courant number of u, the keys warned about)
COURANT = [("WM", 1.5, {"Concentration", "Height"}), ("WM", 0.4, {"Height"}),
           ("WM", 0.3, set()), ("CZ", 1.5, {"Concentration"}), ("CZ", 0.4, set())]


def courant_warnings(workdir, exe, model, courant):
    """The Numerics.Convection keys a short Koren run of front_<model> at
    this Courant number is warned about, and its Warnings.txt."""
    tag = "courant%g" % courant
    advect(workdir, exe, tag, model, "Koren", 100, courant, 0.1, step_c)
    path = os.path.join(workdir, "%s_%s_Koren_n100" % (tag, model), "run", "a", "Warnings.txt")
    with open(path) as f:
        txt = f.read()
    return {k for k in ("Concentration", "Height")
            if 'Courant number above 1' in txt and 'Numerics.Convection.%s = "Koren"' % k in txt}, txt


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--all", action="store_true", help="also QUICK, Upwind1, Upwind2 (reported only)")
    ap.add_argument("--scheme", action="append", help="only these schemes (repeatable)")
    ap.add_argument("--no-cz", action="store_true", help="leave out front_cz")
    ap.add_argument("--oblique", action="store_true", help="only oblique_wm")
    args = ap.parse_args()
    exe = os.path.abspath(args.exe)
    shutil.rmtree(args.workdir, ignore_errors=True)
    plan = CHECKED + (REPORTED if args.all else [])
    if args.scheme:
        plan = [p for p in plan if p[0] in args.scheme]
    if args.oblique:
        plan = []
    bad = 0
    if plan:
        print("%-9s %-24s %-24s %s" % ("scheme", "front_wm min/max", "front_cz min/max",
                                      "smooth_wm error 32/64/128 (orders)"))
    for scheme, lo in plan:
        try:
            fw = front(args.workdir, exe, "WM", scheme)
            fc = None if args.no_cz else front(args.workdir, exe, "CZ", scheme)
            es = smooth(args.workdir, exe, scheme)
        except RuntimeError as e:
            print("FAIL %-9s %s" % (scheme, e))
            bad += 1
            continue
        po = [math.log(es[i][0] / es[i + 1][0]) / math.log(2.0) for i in range(len(es) - 1)]
        p1 = [math.log(es[i][1] / es[i + 1][1]) / math.log(2.0) for i in range(len(es) - 1)]
        ok = True
        tol = 1e-12
        if scheme in LIMITED:
            for mm in (fw, fc):
                if mm is not None and (mm[0] < 0.3 - tol or mm[1] > 0.7 + tol):
                    ok = False
        if lo is not None and po[-1] < lo:
            ok = False
        bad += not ok

        def over(mm):
            if mm is None:
                return "-"
            return "%.6f/%.6f" % mm
        print("%s %-9s %-24s %-24s max %s (%s)  mean %s (%s)" % (
            "ok  " if ok else "FAIL", scheme, over(fw), over(fc),
            "/".join("%.2e" % e[0] for e in es), ", ".join("%.2f" % p for p in po),
            "/".join("%.2e" % e[1] for e in es), ", ".join("%.2f" % p for p in p1)))
    for scheme in OBLIQUE:
        if args.scheme and scheme not in args.scheme:
            continue
        try:
            eo = oblique(args.workdir, exe, scheme)
        except RuntimeError as e:
            print("FAIL oblique_wm %s: %s" % (scheme, e))
            bad += 1
            continue
        po = [math.log(eo[i][0] / eo[i + 1][0]) / math.log(OBL_GRIDS[i + 1][0] / OBL_GRIDS[i][0])
              for i in range(len(eo) - 1)]
        ok = OBLIQUE_BAND[0] <= po[-1] <= OBLIQUE_BAND[1] and abs(po[-2] - po[-1]) <= OBLIQUE_SETTLED
        bad += not ok
        print("%s oblique_wm %-7s error %s (orders %s; %s times the modified equation's), h off by %.0e"
              % ("ok  " if ok else "FAIL", scheme, "/".join("%.2e" % e[0] for e in eo),
                 ", ".join("%.2f" % p for p in po),
                 "/".join("%.2f" % (e[0] / oblique_predicted(n)) for e, (n, _) in zip(eo, OBL_GRIDS)),
                 max(e[1] for e in eo)))
    if not args.oblique and (not args.scheme or "Koren" in args.scheme):
        for model, courant, want in COURANT:
            if model == "CZ" and args.no_cz:
                continue
            try:
                got, txt = courant_warnings(args.workdir, exe, model, courant)
            except RuntimeError as e:
                print("FAIL courant %s %g: %s" % (model, courant, e))
                bad += 1
                continue
            ok = got == want
            bad += not ok
            print("%s courant  %s at %.2f: warned about %s (expected %s)" % (
                "ok  " if ok else "FAIL", model, courant, sorted(got) or "nothing", sorted(want) or "nothing"))
            if not ok:
                print(txt[-600:])
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
