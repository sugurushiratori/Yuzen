#!/usr/bin/env python3
"""Check the features that no other test switches on.

    python3 tools/features.py --exe build/bin/yuzen --workdir /tmp/ft
    python3 tools/features.py --exe build/bin/yuzen --workdir /tmp/ft --only probe_wm_1d

    probe_wm_1d, probe_wm_2d, probe_cz_1d
                 what Probe00.dat writes at the end of the run is what the
                 field file of the same instant holds at the nearest point
                 (free surface, film, substrate, concentration). The check
                 found the probe landing on the first point within a cell of
                 the one asked for, not the nearest, and ConcZ comparing the
                 meters of the input with its dimensionless grid; corrected. For ConcZ the
                 level through the film (the zeta of the header) must be
                 the one nearest Z as well.
    steady_wm, steady_cz
                 a leveling film with SteadyStateCriterion stops with
                 "Steady state" long before EndTime, flat
    dry_wm, dry_cz
                 a drying film with EndConcType "mean" and EndConcValue 0.6
                 stops with "Dry out" in the step in which the mean
                 concentration passes 0.6
    restart_edited_wm, restart_edited_cz
                 a field file edited after the run wrote it: the restart
                 takes it as it stands, warns that it is not the state saved
                 beside it, and an unedited copy is continued without that
                 warning
    restart_pinned
                 a pinned film pushed against its wall and dried, past a
                 tenth of its thickest film at the wall: a restart from it
                 goes on, the line of the check of the pinned edges as a
                 warning (the check, which is of the initial condition,
                 refused it; found by tools/proptest.py, corrected)
    maxsteps     Numerics.MaxSteps stops the run after exactly that many
                 steps, and says so
    adaptive_dt  the time-step control, from a step far too small: the
                 step must grow ("dt+") and the answer must agree with a
                 run at a fixed small step to within the error the control
                 is asked for
    adaptive_clock_wm_1d, adaptive_clock_wm_2d, adaptive_clock_cz_1d,
    adaptive_clock_cz_2d
                 the time-step control from a step far too long (cut back,
                 "dt-") and then lengthened every third step ("dt+"): the
                 time the clock reached is the sum of the steps the
                 equations were solved with (JobInfo.txt), to rounding. The
                 check found the control lengthening the step before the
                 clock had advanced, so the clock added the lengthened step,
                 which had not been solved, and ran ahead (by 0.37 % at the
                 end of leveling_1d_reg); corrected
    adaptive_source
                 a nozzle whose flow rate relaxes as 1 - exp(-t/tau), the
                 step adapted: the volume supplied up to each History.dat
                 row, divided by the integral of the relaxation up to the
                 time the row reports, is the same on every row (the full
                 rate times the footprint). With the clock ahead the source
                 was evaluated, and the time reported, at the wrong times
    coupling_cz  the film-concentration passes of ConcZ with Marangoni that
                 do not settle: the stratified manufactured film of
                 tools/mms.py (cz_strat_1d) on 32 cells in a step of
                 0.02 s: at a fixed step the run must stop and name the
                 passes as the cause; under a step control loose enough to
                 accept that step otherwise, the step must be cut (one
                 warning, counted in JobInfo.txt), the run end normally,
                 and its error at 0.02 s stay below 0.5 of the film's
                 amplitude (0.17; with the step kept as it stood, 2.1)
    evap_precursor_wm, evap_precursor_cz
                 leveling_1d_reg (concz_1d_reg) evaporating, a film over a
                 substrate step that thins below its precursor at the foot
                 of the step: no solvent may condense there, so
                 EvaporatedVol never falls. Found the precursor correction
                 of the rate negative below the precursor, which condensed
                 solvent onto the film; corrected (0 there)
    stagnation_wm
                 a Newton correction that stagnates far from convergence
                 (the first steps of stress_wetting): at a fixed step the run
                 must stop and name the stagnation as the cause; under step
                 control the step must be cut (one warning, the steps
                 counted in JobInfo.txt), the run end normally, and |VolErr|
                 stay at rounding. Found taken with a warning; corrected

Only the Python standard library is used, so the test runs wherever ctest
does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import math
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from concz_vs_wm import edit_toml                       # noqa: E402
from equivalence import fixed_steps, run, read_field  # noqa: E402
from analytic import case, quiet, wave_film, closed_box  # noqa: E402
import mms                                                 # noqa: E402

L = 0.25e-3


def wave(model, twoD):
    """A leveling sinusoid of 30 % on a periodic film."""
    s = quiet(case("concz_1d_reg" if model == "CZ" else "leveling_1d_reg"), twoD=twoD)
    s = edit_toml(s, {}, {"Grid": {"Nx": "32", "Xmin": repr(-L / 2), "Xmax": repr(L / 2)},
                          "Physics.ConcField": {"Model": '"%s"' % model}})
    if model == "CZ":
        s = edit_toml(s, {}, {"Grid": {"Nz": "8"}})
    if twoD:
        s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": "24",
                                       "Ymin": repr(-0.75 * L / 2), "Ymax": repr(0.75 * L / 2)}})
    return wave_film(s, 0.5, 0.15, L, 0.75 * L if twoD else None)


def last_lines(path, n=1):
    rows = [l.split() for l in open(path) if l.strip() and not l.lstrip().startswith("#")
            and re.match(r"\s*[-+0-9.]", l)]
    return rows[-n:]


# ----------------------------------------------------------------------
def probe(model, twoD):
    def check(workdir, exe):
        s = wave(model, twoD)
        px, py = 0.037e-3, -0.021e-3
        s = fixed_steps(s, 2e-5, 100, 1e-3)
        s = edit_toml(s, {}, {"Output.ProbeList": {"TimeInterval": "1e-3"}})
        s = s.rstrip() + "\n\n[[Output.ProbeList.Probe]]\nX = %r\nY = %r\nZ = 2.0e-6\n" % (px, py)
        d = os.path.join(workdir, "probe_%s_%s" % (model, "2d" if twoD else "1d"))
        rc, log = run(exe, d, "a", s, 600)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        res = os.path.join(d, "a")
        pf = os.path.join(res, "Probe00.dat")
        if not os.path.exists(pf):
            return None, "no Probe00.dat"
        t, hs, cc, ux, uy, hh, ss = [float(v) for v in last_lines(pf)[0][:7]]
        f, nx, ny, nz = read_field(sorted(os.path.join(res, n) for n in os.listdir(res)
                                          if n.startswith("field") and n.endswith(".dat"))[-1])
        #  The nearest point of the field file (coordinates in mm)
        best, r = 1e30, 0
        for k in range(len(f["xx"])):
            dx = f["xx"][k] - px * 1e3
            dy = (f["yy"][k] - py * 1e3) if "yy" in f else 0.0
            dd = dx * dx + dy * dy
            if dd < best:
                best, r = dd, k
        err = max(abs(hs - f["HS"][r]), abs(hh - f["HH"][r]), abs(ss - f["SS"][r])) / f["HH"][r]
        err = max(err, abs(cc - f["CC"][r]))
        #  The header names the point the data are from (to the six digits
        #  it is printed with), and that point is on the side of the probe:
        #  the film is symmetric about x = 0, so the data alone would not
        #  tell the probe from its mirror image.
        head = open(pf).readline()
        m = re.search(r"\[mm\] = \(([-+0-9.eE]+)(?:, ([-+0-9.eE]+))?", head)
        if not m:
            return None, "no position in the header of Probe00.dat: " + head.strip()
        where = [(float(m.group(1)), f["xx"][r])]
        if twoD:
            where.append((float(m.group(2)), f["yy"][r]))
        for got, data in where:
            if abs(got - data) > 1e-3 * (L * 1e3 / 32):
                return None, "the header gives %s, the data are those of %s" % (head.strip(), (f["xx"][r], f["yy"][r] if twoD else None))
        if model == "CZ":
            #  The level through the film: the cell whose zeta is nearest
            #  Z over the thickness of the probe's column at the start
            #  (the concentration of this film is the same at every level,
            #  so only the header tells which cell was taken).
            m = re.search(r"zeta = ([-+0-9.eE]+)", open(pf).readline())
            if not m:
                return None, "no zeta in the header of Probe00.dat"
            want = 2.0 / float(last_lines(pf, 10 ** 6)[0][5])     # Z / h, both in um
            levels = [f["zz"][k] / f["HH"][k] for k in range(len(f["xx"]))
                      if f["xx"][k] == f["xx"][r] and 0.0 < f["zz"][k] < f["HH"][k]]
            near = min(levels, key=lambda z: abs(z - want))
            err = max(err, abs(float(m.group(1)) - near))
        return err, None
    return check


def stops(model, kind):
    def check(workdir, exe):
        s = wave(model, False)
        if kind == "steady":
            s = edit_toml(s, {}, {"Numerics": {"InitialDt": "1e-3", "EndTime": "100.0",
                                               "MaxSteps": "100000",
                                               "SteadyStateCriterion": "1e-3", "SteadyStateSteps": "10"},
                                  "Numerics.TimeStepControl": {"Enable": "true"},
                                  "Output.FieldData": {"TimeInterval": "10.0", "NeutralFile": "true",
                                                       "VTKFile": "false"},
                                  "Output.GeneralHistory": {"TimeInterval": "1e-2"}})
        else:
            s = edit_toml(s, {}, {"Physics.Evaporation": {
                "Enable": "true", "Coef": "2e-5", "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"',
                "EndConcType": '"mean"', "EndConcValue": "0.6"},
                "Fluid.Diffusivity": {"BaseValue": "1e-8", "Exponent": "0.0"}})
            s = fixed_steps(s, 1e-3, 100000, 1.0)
            s = edit_toml(s, {}, {"Output.GeneralHistory": {"TimeInterval": "1e-3"}})
        d = os.path.join(workdir, "%s_%s" % (kind, model))
        rc, log = run(exe, d, "a", s, 900)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        want = "Steady state" if kind == "steady" else "Dry out"
        if want not in log:
            return None, "the run did not report %r" % want
        hist = last_lines(os.path.join(d, "a", "History.dat"), 2)
        if kind == "steady":
            t = float(hist[-1][0])
            if not t < 50.0:
                return None, "it ran to t = %g" % t
            f, nx, ny, nz = read_field(sorted(os.path.join(d, "a", n) for n in os.listdir(os.path.join(d, "a"))
                                              if n.startswith("field") and n.endswith(".dat"))[-1])
            hh = f["HH"][1:-1] if model == "WM" else f["HH"][nx + 1:2 * nx - 1]
            return (max(hh) - min(hh)) / 3.0, None     # of the 3 um the wave spanned
        #  AveC is column 4; the step that stopped must be the first past 0.6
        prev, last = float(hist[0][3]), float(hist[1][3])
        return (0.0 if (prev <= 0.6 < last) else 1.0), None
    return check


def restart_edited(model):
    """A field file edited after it was written: the restart must take it as
    it stands, say that it is not the state saved beside it, and still read
    the rest (time, step history) from the state file."""
    def check(workdir, exe):
        d = os.path.join(workdir, "restart_edited_" + model)
        shutil.rmtree(d, ignore_errors=True)
        s = fixed_steps(wave(model, False), 1e-4, 20, 1e-3)
        rc, log = run(exe, d, "a", s, 600)
        if rc != 0:
            return None, "source run failed: " + log[-400:]
        t1 = "1e-3"     # field00001: 10 fixed steps of 1e-4 s
        #  The copy as written (u) and one with a single film value raised by 1 % (e)
        for sub, edit in (("u", False), ("e", True)):
            os.makedirs(os.path.join(d, sub))
            shutil.copy(os.path.join(d, "a", "field00001.state"), os.path.join(d, sub))
            lines = open(os.path.join(d, "a", "field00001.dat")).read().split("\n")
            if edit:
                #  every row of that column (ConcZ writes one per level)
                head = [c.lstrip("#").split("[")[0].strip() for c in lines[0].split("\t")]
                xs = sorted(set(l.split("\t")[0] for l in lines[1:] if l.strip()), key=float)
                x0 = xs[len(xs) // 2]       # a column in the middle, not on a wall
                for k in range(1, len(lines)):
                    row = lines[k].split("\t")
                    if row[0] != x0:
                        continue
                    for name in ("HH", "HS"):
                        i = head.index(name)
                        row[i] = repr(float(row[i]) * 1.01)
                    lines[k] = "\t".join(row)
            open(os.path.join(d, sub, "field00001.dat"), "w").write("\n".join(lines))
        out = {}
        for sub in ("u", "e"):
            r = edit_toml(s, {}, {"Restart": {"Enable": "true", "ASCII": "true",
                                              "FileName": '"%s/field00001.dat"' % sub,
                                              "RestartTime": t1, "FileIndex": "2"},
                                  "Numerics": {"MaxSteps": "10"}})
            rc, log = run(exe, d, "r" + sub, r, 600)
            if rc != 0:
                return None, "restart from %s failed: %s" % (sub, log[-400:])
            warn = open(os.path.join(d, "r" + sub, "Warnings.txt")).read()
            out[sub] = (warn, read_field(os.path.join(d, "r" + sub, "field00002.dat"))[0])
        if "is not the state saved beside it" in out["u"][0]:
            return None, "the unedited copy was reported as edited"
        if "is not the state saved beside it" not in out["e"][0]:
            return None, "the edited field file was not reported"
        if out["u"][1]["HH"] == out["e"][1]["HH"]:
            return None, "the edit did not reach the continued run"
        return 0.0, None
    return check


def restart_pinned(workdir, exe):
    """A pinned film pushed against its wall and dried (tools/
    restart_roundtrip.py, PINNED_DRY), its film at the wall past a tenth of
    the thickest by 0.4 s: the check of the pinned edges refuses the same
    film as an initial condition, but a restart from it goes on, with the
    line as a warning (tools/proptest.py seed 1206 is such a restart)."""
    import restart_roundtrip as rr
    d = os.path.join(workdir, "restart_pinned")
    shutil.rmtree(d, ignore_errors=True)
    v = rr.CHAIN["wm1d_pinned_dry"]
    base = rr.toml_edit(case(v["case"]), {k: dict(x) for k, x in v["sect"].items()}, v["append"])
    s = edit_toml(base, {}, {"Numerics": {"EndTime": "0.4"}})
    rc, log = run(exe, d, "a", s, 600)
    if rc != 0:
        return None, "source run failed: " + log[-400:]
    r = edit_toml(s, {}, {"Restart": {"Enable": "true", "ASCII": "false",
                                      "FileName": '"a/field00004.vts"',
                                      "RestartTime": "0.4", "FileIndex": "5"},
                          "Numerics": {"EndTime": "0.45"}})
    rc, log = run(exe, d, "r", r, 600)
    if rc != 0:
        return None, "the restart was refused: " + log[-400:]
    warn = open(os.path.join(d, "r", "Warnings.txt")).read()
    if "Pinning Xmin: the restarted film reaches" not in warn or "taken as it is" not in warn:
        return None, "no warning about the film at the pinned wall: " + warn[-300:]
    return 0.0, None


def maxsteps(workdir, exe):
    s = wave("WM", False)
    s = edit_toml(s, {}, {"Numerics": {"InitialDt": "1e-5", "EndTime": "100.0", "MaxSteps": "37",
                                       "SteadyStateCriterion": "-1"},
                          "Numerics.TimeStepControl": {"Enable": "false"}})
    d = os.path.join(workdir, "maxsteps")
    rc, log = run(exe, d, "a", s, 600)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    if "MaxSteps" not in log:
        return None, "the run did not say it stopped at MaxSteps"
    t = float(last_lines(os.path.join(d, "a", "History.dat"))[0][0])
    return abs(t - 37e-5) / 1e-5, None


def adaptive_dt(workdir, exe):
    s = wave("WM", False)
    T = 0.02
    a = edit_toml(s, {}, {"Numerics": {"InitialDt": "1e-9", "EndTime": repr(T), "MaxSteps": "1000000",
                                       "SteadyStateCriterion": "-1"},
                          "Numerics.TimeStepControl": {"Enable": "true", "Criterion": "1e-4"},
                          "Output.FieldData": {"TimeInterval": repr(T), "NeutralFile": "true",
                                               "VTKFile": "false"}})
    d = os.path.join(workdir, "adaptive_dt")
    rc, log = run(exe, d, "a", a, 900)
    if rc != 0:
        return None, "adaptive run failed: " + log[-400:]
    grew = len(re.findall(r"dt\+ \(", log))
    if grew < 10:
        return None, "the step grew only %d times" % grew
    b = fixed_steps(s, T / 2000, 2000, T)
    rc, log2 = run(exe, d, "b", b, 900)
    if rc != 0:
        return None, "reference run failed: " + log2[-400:]
    fa = read_field(sorted(os.path.join(d, "a", n) for n in os.listdir(os.path.join(d, "a"))
                           if n.startswith("field") and n.endswith(".dat"))[-1])[0]
    fb = read_field(sorted(os.path.join(d, "b", n) for n in os.listdir(os.path.join(d, "b"))
                           if n.startswith("field") and n.endswith(".dat"))[-1])[0]
    ta = float(last_lines(os.path.join(d, "a", "History.dat"))[0][0])
    if abs(ta - T) > 1e-3 * T:
        return None, "the adaptive run ended at t = %g, not %g" % (ta, T)
    change = max(abs(x - 5.0) for x in fb["HH"][1:-1])
    return max(abs(x - y) for x, y in zip(fa["HH"][1:-1], fb["HH"][1:-1])) / change, None


def jobinfo_value(path, title):
    for line in open(path):
        if line.startswith(title + " :"):
            return float(line.split(":", 1)[1])
    raise ValueError("%s: no '%s'" % (path, title))


def adaptive_clock(model, twoD):
    #  In 2-D a step far too long is cut back and then lengthened again in
    #  one run; in 1-D it is not lengthened again, so a second run starts
    #  from a step far too short.
    starts = [("2e-3", "1e-3")] if twoD else [("2e-3", "1e-3"), ("1e-9", "1e-4")]

    def check(workdir, exe):
        T = 0.02
        worst, logs = 0.0, ""
        for k, (dt0, crit) in enumerate(starts):
            s = edit_toml(wave(model, twoD), {}, {
                "Numerics": {"InitialDt": dt0, "EndTime": repr(T), "MaxSteps": "1000000",
                             "SteadyStateCriterion": "-1"},
                "Numerics.TimeStepControl": {"Enable": "true", "Criterion": crit,
                                             "KeepSteps": "3", "ChangeFactor": "0.75",
                                             "MaxDtFactor": "1e10", "MinDt": "1e-15"},
                "Output.FieldData": {"TimeInterval": repr(T), "NeutralFile": "true",
                                     "VTKFile": "false"}})
            d = os.path.join(workdir, "adaptive_clock_%s_%s" % (model, "2d" if twoD else "1d"))
            rc, log = run(exe, d, "a%d" % k, s, 900)
            if rc != 0:
                return None, "run failed: " + log[-400:]
            logs += log
            ji = os.path.join(d, "a%d" % k, "JobInfo.txt")
            try:
                t0 = jobinfo_value(ji, "Time at the start of the loop")
                t1 = jobinfo_value(ji, "Time at the end of the loop")
                sdt = jobinfo_value(ji, "Sum of the steps solved")
            except ValueError as e:
                return None, str(e)
            worst = max(worst, abs((t1 - t0) - sdt) / (t1 - t0))
        if "dt+" not in logs or "dt-" not in logs:
            return None, "the step was not both cut back and lengthened"
        return worst, None
    return check


def adaptive_source(workdir, exe):
    tau = 0.01
    T = 0.05
    s = closed_box(False, nozzle=True)
    s = fixed_steps(s, 1e-7, 1000000, T)
    s = edit_toml(s, {}, {"Physics.Source": {"TimeConstant": repr(tau)},
                          "Numerics": {"EndTime": repr(T)},
                          "Numerics.TimeStepControl": {"Enable": "true", "KeepSteps": "3",
                                                       "ChangeFactor": "0.75", "Criterion": "1e-2",
                                                       "MaxDtFactor": "1e10", "MinDt": "1e-15"},
                          "Output.GeneralHistory": {"TimeInterval": repr(tau / 20)}})
    d = os.path.join(workdir, "adaptive_source")
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    if "dt+" not in log:
        return None, "the step was never lengthened"
    hist = [l.split("\t") for l in open(os.path.join(d, "a", "History.dat"))]
    k = [c.split("[")[0].strip() for c in hist[0]].index("SuppliedVol")
    ratios = []
    for r in hist[1:]:
        t = float(r[0])
        if t < 2 * tau:
            continue        # the supply is still small: the 7 digits of the file would show
        ratios.append(float(r[k]) / (t - tau * (1.0 - math.exp(-t / tau))))
    if len(ratios) < 10:
        return None, "only %d History.dat rows past 2 tau" % len(ratios)
    return (max(ratios) - min(ratios)) / ratios[-1], None


def coupling_cz(workdir, exe):
    #  The step is twice the largest that settles on this state (1e-2 s does,
    #  in four passes; ch-verification, "The coupling passes of CZ"); the
    #  passes stop at kMaxCouplingPasses = 8. A step whose passes do not
    #  settle has failed: at a fixed step the run stops and says why ...
    wd = os.path.join(workdir, "coupling_cz")
    terms = mms.cz_strat_terms(False, {32: 8})
    try:
        mms.run_mms(wd, exe, "CZ", False, 32, 0.02, dict(terms), tag="c")
        return None, "the run at a fixed step went on past passes that did not settle"
    except RuntimeError as e:
        why = str(e)
    if "did not settle in 8 coupling passes" not in why or "cannot be retried at a shorter one" not in why:
        return None, "the run did not stop for the passes: " + why[-400:]
    #  ... and under step control the step is cut and solved again. The
    #  criterion of the step control is loose enough to accept the step of
    #  0.02 s, so that only the passes cut it. Kept as it stood, that step
    #  left an error of 2.1 (of the film's amplitude, as mms.py measures it);
    #  cut, 0.17. What is returned is that error.
    try:
        err = mms.run_mms(wd, exe, "CZ", False, 32, 0.02, dict(terms), tag="a", after={
            "Numerics": {"MaxSteps": "1000", "EndTime": repr(mms.T_END)},
            "Numerics.TimeStepControl": {"Enable": "true", "Criterion": "1.0", "KeepSteps": "3",
                                         "ChangeFactor": "0.75", "MaxDtFactor": "1e10", "MinDt": "1e-12"}})[0]
    except RuntimeError as e:
        return None, "run failed: " + str(e)[-400:]
    res = os.path.join(wd, "CZa_1d_n32_dt0.02", "a")
    warn = open(os.path.join(res, "Warnings.txt")).read()
    if warn.count("did not settle in 8 coupling passes") != 1 or "the step is cut" not in warn:
        return None, "not one warning that the passes did not settle and the step was cut"
    ji = open(os.path.join(res, "JobInfo.txt")).read()
    m = re.search(r"Steps failed: the passes did not settle :\s*(\d+)", ji)
    if not m or int(m.group(1)) < 1:
        return None, "steps failed for the passes: %s, not 1 or more" % (m.group(1) if m else "?")
    if not re.search(r"Termination Type :\s*normal", ji):
        return None, "the run did not end normally"
    return err, None


def evap_precursor(model):
    def check(workdir, exe):
        #  leveling_1d_reg: a film draped over a substrate step, the
        #  precursor at half its scale, evaporating (Raoult). A dip opens
        #  below the step where the film thins under the precursor; no
        #  solvent may condense there, so EvaporatedVol never falls.
        s = case("concz_1d_reg" if model == "CZ" else "leveling_1d_reg")
        s = fixed_steps(s, 1e-3, 200, 0.1)
        s = edit_toml(s, {}, {"Output.GeneralHistory": {"TimeInterval": "1e-3"}})
        d = os.path.join(workdir, "evap_precursor_" + model)
        rc, log = run(exe, d, "a", s, 900)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        hist = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "a", "History.dat"))]
        names = [c.split("[")[0].strip() for c in hist[0]]
        iv, ie, ih = names.index("Volume"), names.index("EvaporatedVol"), names.index("minH")
        rows = [[float(r[k]) for k in (iv, ie, ih)] for r in hist[1:] if len(r) == len(names)]
        if not min(r[2] for r in rows) < 5e-6:
            return None, "the film never went below the precursor (minH %g)" % min(r[2] for r in rows)
        v0 = rows[0][0]
        #  The largest fall of EvaporatedVol from one row to the next,
        #  relative to the initial volume.
        return max([0.0] + [(a[1] - b[1]) / v0 for a, b in zip(rows, rows[1:])]), None
    return check


def stagnation_wm(workdir, exe):
    #  stress_wetting: a slab meeting a thick precursor film, far from
    #  equilibrium. At its InitialDt (1e-6 s) the Newton correction of the
    #  first step stagnates at 2.7 times the first correction, and at the
    #  shorter steps tried after it at 0.068 to 25 times: far from convergence, so the step has failed. At a fixed
    #  step the run stops and says why ...
    s = case("stress_wetting")
    d = os.path.join(workdir, "stagnation_wm")
    f = edit_toml(s, {}, {"Numerics": {"MaxSteps": "5"},
                          "Numerics.TimeStepControl": {"Enable": "false"},
                          "Output.FieldData": {"VTKFile": "false"}})
    rc, log = run(exe, d, "f", f, 300)
    if rc == 0:
        return None, "the run at a fixed step went on past a correction that stagnated"
    if "stagnated far from convergence" not in log or "cannot be retried at a shorter one" not in log:
        return None, "the run did not stop for the stagnation: " + log[-400:]
    #  ... and under step control the step is cut until it converges (thirteen
    #  cuts for the stagnation, the first warned), and the run goes on. The
    #  steps the Newton iteration stagnated in were rejected by the test on
    #  the truncation error before (the fields are the same to the bit); a
    #  correction that stagnates at its floor is still accepted
    #  (yuzen_tests). What is returned is the largest |VolErr|.
    a = edit_toml(s, {}, {"Numerics": {"EndTime": "2e-8"},
                          "Output.FieldData": {"VTKFile": "false"}})
    rc, log = run(exe, d, "a", a, 300)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    res = os.path.join(d, "a")
    warn = open(os.path.join(res, "Warnings.txt")).read()
    if warn.count("stagnated far from convergence") != 1 or "the step is cut" not in warn:
        return None, "not one warning that the correction stagnated and the step was cut"
    ji = open(os.path.join(res, "JobInfo.txt")).read()
    m = re.search(r"Steps failed: the Newton correction stagnated :\s*(\d+)", ji)
    if not m or int(m.group(1)) < 1:
        return None, "steps failed on a stagnation: %s, not 1 or more" % (m.group(1) if m else "?")
    if not re.search(r"Termination Type :\s*normal", ji):
        return None, "the run did not end normally"
    hist = [l.split("\t") for l in open(os.path.join(res, "History.dat"))]
    k = [c.split("[")[0].strip() for c in hist[0]].index("VolErr")
    return max(abs(float(r[k])) for r in hist[1:]), None


#  The probe file carries seven significant digits
CHECKS = [("probe_wm_1d", probe("WM", False), 1e-6),
          ("probe_wm_2d", probe("WM", True), 1e-6),
          ("probe_cz_1d", probe("CZ", False), 1e-6),
          ("steady_wm", stops("WM", "steady"), 1e-3),
          ("steady_cz", stops("CZ", "steady"), 1e-3),
          ("dry_wm", stops("WM", "dry"), 0.0),
          ("dry_cz", stops("CZ", "dry"), 0.0),
          ("restart_edited_wm", restart_edited("WM"), 0.0),
          ("restart_edited_cz", restart_edited("CZ"), 0.0),
          ("restart_pinned", restart_pinned, 0.0),
          ("maxsteps", maxsteps, 1e-6),
          ("adaptive_dt", adaptive_dt, 1e-2),
          ("adaptive_clock_wm_1d", adaptive_clock("WM", False), 1e-12),
          ("adaptive_clock_wm_2d", adaptive_clock("WM", True), 1e-12),
          ("adaptive_clock_cz_1d", adaptive_clock("CZ", False), 1e-12),
          ("adaptive_clock_cz_2d", adaptive_clock("CZ", True), 1e-12),
          ("adaptive_source", adaptive_source, 1e-3),
          ("coupling_cz", coupling_cz, 0.5),
          ("evap_precursor_wm", evap_precursor("WM"), 1e-15),
          ("evap_precursor_cz", evap_precursor("CZ"), 1e-15),
          ("stagnation_wm", stagnation_wm, 1e-12)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", action="append")
    args = ap.parse_args()
    exe = os.path.abspath(args.exe)
    shutil.rmtree(args.workdir, ignore_errors=True)
    bad = 0
    for name, fn, tol in CHECKS:
        if args.only and name not in args.only:
            continue
        v, why = fn(args.workdir, exe)
        if why:
            print("FAIL %-12s %s" % (name, why))
            bad += 1
            continue
        ok = v <= tol
        bad += not ok
        print("%s %-12s %.2e (allowed %g)" % ("ok  " if ok else "FAIL", name, v, tol))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
