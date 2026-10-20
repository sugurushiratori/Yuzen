#!/usr/bin/env python3
"""Compare the ConcZ model against well-mixed in the limit of a
concentration that is uniform through the film thickness.

    python3 tools/concz_vs_wm.py --exe build/bin/yuzen
    python3 tools/concz_vs_wm.py --exe build/bin/yuzen --nx 120 --nz 21

ConcZ is a model that carries concentration and viscosity profiles
through the thickness, so in the limit where those profiles are absent
it is solving the same problem as well-mixed. Here both are run under
conditions that produce that limit, and the time evolution of the film
thickness is compared.

Constructing the limit
    Fluid.Viscosity.Exponent = 0    removes the concentration dependence
                                    of the viscosity. Concentration then
                                    no longer enters the thickness
                                    equation, so h is the same whether or
                                    not there is a profile in z
    Physics.Evaporation      = off  keeps the concentration static
    Physics.SolutoCapillarity= off  no flow driven by concentration
                                    gradients

How to run the comparison (get this wrong and the numbers are meaningless)
    - Fix the time step. With adaptive stepping left on, the step
      histories of the two models diverge and the output times no
      longer line up. You then compare "the film thickness at different
      times", see a difference of a few percent that does not shrink
      when the grid is refined, and mistake it for a discretization
      difference
    - Take an output interval well below EndTime. If it is too coarse,
      the only common time is the initial condition and the difference
      looks like 0
    - Turn off the steady-state test. If one run is cut short before the
      other there is nothing to compare

Two-dimensional (pseudo three-dimensional) comparison
    Passing --ny compares with Grid.Dimension = "xyh". In two dimensions
    the initial film thickness is uniform, so with nothing driving it
    both models stay at rest, the difference looks like 0, and nothing
    is verified. A bump is placed on the substrate so that capillary
    pressure drives a flow (the substrate does not touch the
    concentration, so the limit is preserved).
    The x-y and y-zeta cross terms only come into play here, so this is
    the case to use when checking the pseudo three-dimensional
    implementation.

Expected result
    In either dimension the two agree to the printed precision (of order
    1e-6 in relative terms). Both models are in conservative form and use
    the same positivity preserving mobility (Eq. 23 of Diez & Kondic
    2002), but the fluxes are constructed by independent implementations,
    so agreement supports the validity of both. Volume is conserved
    exactly in both.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
BASE_CASE = os.path.join(ROOT, "cases", "concz_1d_reg.toml")


def edit_toml(src, top=None, sect=None):
    """Rewrite only the keys of the relevant TOML sections.

    Rewriting with a regular expression matches greedily across section
    boundaries and corrupts the file. Instead, scan line by line and keep
    track of which section we are currently in.
    """
    top = top or {}
    sect = sect or {}
    out, cur, seen = [], None, set()
    for line in src.split("\n"):
        st = line.strip()
        if st.startswith("[") and st.endswith("]"):
            cur = st[1:-1]
            out.append(line)
            continue
        if "=" in st and not st.startswith("#"):
            key = st.split("=", 1)[0].strip()
            if cur is None and key in top:
                out.append("%s = %s" % (key, top[key]))
                seen.add(key)
                continue
            if cur in sect and key in sect[cur]:
                out.append("%s = %s" % (key, sect[cur][key]))
                seen.add(cur + "." + key)
                continue
        out.append(line)
    #  Keys absent from the original case are added just after the
    #  relevant section
    for sname, kv in sect.items():
        add = [(k, v) for k, v in kv.items() if (sname + "." + k) not in seen]
        if not add:
            continue
        res, placed = [], False
        for line in out:
            res.append(line)
            if not placed and line.strip() == "[" + sname + "]":
                for k, v in add:
                    res.append("%s = %s" % (k, v))
                placed = True
        if not placed:
            res.append("[" + sname + "]")
            for k, v in add:
                res.append("%s = %s" % (k, v))
        out = res
    return "\n".join(out)


def make_case(path, job, model, nx, nz, dt, endtime, out_interval, ny=None):
    """Passing ny builds a two-dimensional (pseudo three-dimensional) case.

    In two dimensions the initial film thickness is uniform, so with
    nothing driving it both models stay at rest, the difference looks
    like 0, and nothing is verified. A bump is placed on the substrate so
    that capillary pressure drives a flow. The substrate does not touch
    the concentration, so the limit in which concentration does not
    affect h is preserved.
    """
    src = open(BASE_CASE, encoding="utf-8").read()
    sect = {
        "Grid": {"Nx": str(nx), "Nz": str(nz)},
        "Physics.ConcField": {"Model": '"%s"' % model},
        "Numerics": {"InitialDt": dt, "EndTime": endtime,
                     "MaxSteps": "100000000", "SteadyStateCriterion": "-1"},
        #  Fix the time step; the comparison depends on it
        "Numerics.TimeStepControl": {"Enable": "false"},
        #  Take the limit in which concentration does not enter the
        #  thickness equation
        "Fluid.Viscosity": {"Exponent": "0.0"},
        "Physics.Evaporation": {"Enable": "false"},
        "Physics.SolutoCapillarity": {"Enable": "false"},
        "Output.GeneralHistory": {"TimeInterval": out_interval},
        "Output.FieldData": {"TimeInterval": out_interval},
    }
    if ny is not None:
        sect["Grid"].update({"Dimension": '"xyh"', "Ny": str(ny),
                             "Ymin": "-0.5e-3", "Ymax": "0.5e-3"})
        sect["Substrate.EllipseBump"] = {
            "CreateType": '"add"', "Height": "0.3",
            "Center": "[0.0, 0.0]", "Width": "[0.4e-3, 0.4e-3]",
            "SlopeWidth": "0.1e-3"}
    open(path, "w", encoding="utf-8").write(edit_toml(src, {"JobName": '"%s"' % job}, sect))


def read_history(path):
    #  Column labels carry their unit, as in "maxH [m]". Key the rows on the
    #  bare name so a change of unit in the header does not break the lookup.
    rows, head = [], None
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if head is None:
            head = [c.split("[")[0].strip() for c in f]
            continue
        if len(f) < 2:
            continue
        rows.append(dict(zip(head, [float(x) for x in f])))
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="yuzen executable")
    ap.add_argument("--nx", type=int, default=60)
    ap.add_argument("--ny", type=int, default=None,
                    help="if given, compare in two dimensions (pseudo 3-D)")
    ap.add_argument("--nz", type=int, default=11,
                    help="number of through-thickness points on the ConcZ side")
    ap.add_argument("--dt", default="1e-6")
    ap.add_argument("--endtime", default="2e-3")
    ap.add_argument("--interval", default="2e-4", help="output interval")
    ap.add_argument("--workdir", default=None,
                    help="default is ./concz_vs_wm")
    ap.add_argument("--keep", action="store_true",
                    help="do not delete the run output")
    ap.add_argument("--tol", type=float, default=None,
                    help="upper bound on the relative difference in maxH; "
                         "exit code 1 if exceeded")
    ap.add_argument("--vol-tol", type=float, default=None,
                    help="upper bound on the relative change in volume; "
                         "exit code 1 if exceeded")
    ap.add_argument("--vol-match-tol", type=float, default=1e-6,
                    help="upper bound on the relative difference between the "
                         "Volume the two models write for the same film "
                         "(default 1e-6, the seven digits of History.dat); "
                         "exit code 1 if exceeded")
    args = ap.parse_args()

    exe = os.path.abspath(args.exe)
    if not os.path.exists(exe):
        sys.exit("not found: %s" % exe)

    wd = os.path.abspath(args.workdir or os.path.join(os.getcwd(), "concz_vs_wm"))
    os.makedirs(wd, exist_ok=True)

    runs = {}
    for model, nz in (("WM", 1), ("CZ", args.nz)):
        job = "cmp_%s" % model.lower()
        toml = os.path.join(wd, job + ".toml")
        make_case(toml, job, model, args.nx, nz, args.dt,
                  args.endtime, args.interval, args.ny)
        outdir = os.path.join(wd, job)
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        r = subprocess.run([exe, toml], cwd=wd,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            sys.stdout.write(r.stdout.decode("utf-8", "replace")[-2000:])
            sys.exit("%s run failed (exit %d)" % (model, r.returncode))
        runs[model] = read_history(os.path.join(outdir, "History.dat"))

    wm, cz = runs["WM"], runs["CZ"]
    n = min(len(wm), len(cz))
    if n < 2:
        sys.exit("not enough output rows; lower --interval")

    print("ConcZ vs well-mixed, uniform-concentration limit")
    if args.ny is None:
        print("  grid      Nx = %d,  ConcZ Nz = %d  (1-D film)" % (args.nx, args.nz))
    else:
        print("  grid      Nx = %d, Ny = %d,  ConcZ Nz = %d  (2-D film)"
              % (args.nx, args.ny, args.nz))
    print("  time step %s (fixed),  end %s" % (args.dt, args.endtime))
    print()
    print("%-13s %-15s %-15s %-11s" % ("t", "WM maxH", "CZ maxH", "rel.diff"))
    print("-" * 56)
    worst = 0.0
    for i in range(n):
        ta, tb = wm[i]["Time"], cz[i]["Time"]
        if abs(ta - tb) > 1e-9 * max(1.0, abs(ta)):
            print("!! output times differ: %.6e vs %.6e" % (ta, tb))
            print("   (adaptive stepping left on, or intervals differ)")
            sys.exit(2)
        a, b = wm[i]["maxH"], cz[i]["maxH"]
        d = abs(a - b) / abs(a) if a else 0.0
        worst = max(worst, d)
        print("%-13.4e %-15.7e %-15.7e %-11.3e" % (ta, a, b, d))
    print("-" * 56)
    print("worst relative difference in maxH: %.3e" % worst)

    #  Check volume conservation as well. Both models are in conservative
    #  form, so with no flux through the boundaries in this setup the
    #  volume must stay constant.
    print()
    print("%-9s %-17s %-17s %-11s" % ("model", "V(0)", "V(end)", "rel.change"))
    print("-" * 56)
    volerr = {}
    for name, rows in (("WM", wm), ("ConcZ", cz)):
        v0, v1 = rows[0]["Volume"], rows[-1]["Volume"]
        rel = abs(v1 - v0)/abs(v0) if v0 else float("nan")
        volerr[name] = rel
        print("%-9s %-17.9e %-17.9e %-11.3e" % (name, v0, v1, rel))

    #  And the two report the same quantity: the volume of the film in the
    #  same units (a 1-D film per meter of width). A ConcZ in 1-D that took
    #  the Ymin..Ymax of the input as the width would write (Ymax - Ymin) / L0
    #  times the well-mixed volume (1.0e-8 against 8.7e-12 m^3 for this film).
    vmatch = 0.0
    for i in range(n):
        a, b = wm[i]["Volume"], cz[i]["Volume"]
        vmatch = max(vmatch, abs(a - b) / abs(a) if a else float("inf"))
    print()
    print("worst relative difference in Volume between the models: %.3e" % vmatch)

    print()
    print("The two models solve the same problem in this limit, by different")
    print("routes: well-mixed advances h from the Diez & Kondic positivity")
    print("preserving mobility (Eq. 23), ConcZ from a flux integrated through")
    print("the thickness. Agreement to the printed precision is the check.")

    rc = 0
    if args.tol is not None and worst > args.tol:
        print()
        print("FAIL: maxH differs by %.3e, above the allowed %.3e"
              % (worst, args.tol))
        rc = 1
    if not vmatch <= args.vol_match_tol:
        print()
        print("FAIL: the Volume of the two models differs by %.3e, above the allowed %.3e"
              % (vmatch, args.vol_match_tol))
        rc = 1
    if args.vol_tol is not None:
        for name, rel in volerr.items():
            if rel > args.vol_tol:
                print()
                print("FAIL: %s volume changed by %.3e, above the allowed %.3e"
                      % (name, rel, args.vol_tol))
                rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
