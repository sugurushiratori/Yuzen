#!/usr/bin/env python3
"""Check that a restart starts from the state it was given.

    python3 tools/restart_roundtrip.py --exe build/bin/yuzen --workdir /tmp/rt

For each model -- well-mixed 1-D and 2-D, ConcZ 1-D and 2-D -- a shipped case
is run for a few steps, its last field is used as a restart, both from the
binary .vts and from the ASCII field*.dat, and the state the restarted run
writes at the restart time is compared with the state it was given. Also
checked: that the restarted run takes steps (its clock was set in seconds,
not in units of the time scale), and that its .pvd points at the files it
wrote, stamped with the time it restarted from. ConcZ is restarted a third
time, from the .vts rewritten without the film thickness ("old"), so the
film has to be taken back from the point heights.

Before the restarts, the .vts and the field*.dat the source run wrote of the
same state are checked against each other: the .vts must hold together as
ParaView reads it (each array where its offset says, as long as it should
be, the blocks back to back up to the end of the appended data), and carry
the concentration, film, substrate, velocity and point positions of the
.dat (for the well-mixed model, the velocity is also the sum of the
per-term velocities). The restart reads the blocks in order and never the
offsets, and nothing read the velocities, so none of that was checked.

During development every one of these failed for ConcZ: the time was not
divided by T0, so the run started past EndTime; the film came back thinner
by (nz-1)/nz; the ASCII flag was ignored; and after a restart of either
model the .pvd listed field00000.vts onwards. Nothing in the unit tests reached any of it.

These restarts are given the field file alone (no field*.state beside it),
which is how a field file without its state file (for example a hand-made
one) is read.

    python3 tools/restart_roundtrip.py --chain --exe build/bin/yuzen --workdir /tmp/rt

compares runs stopped and restarted again and again with the same run made
straight through (see the comment above CHAIN below): with the state file,
to the bit, in every field file, every row of History.dat, CalcInfo.dat and
the probe files, the steady-state and dry-out stops and the iteration counts.

Only the Python standard library is used, so the test runs wherever ctest
does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import os
import re
import shutil
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(os.path.dirname(HERE), "cases")

#  (case, steps before the restart). Few steps: the point is the restart,
#  not the physics, and the 2-D cases are the expensive ones.
PLAN = [("leveling_1d_reg", 40),
        ("leveling_2d_reg", 6),
        ("concz_1d_reg", 40),
        ("concz_2d", 6)]


def read_field(path):
    """A field*.dat as {column name: [values]}."""
    with open(path) as f:
        head = f.readline().lstrip("#").strip().split("\t")
        names = [c.split("[")[0].strip() for c in head if c.strip()]
        cols = [[] for _ in names]
        for line in f:
            v = line.split()
            if len(v) != len(names):
                continue
            for k, x in enumerate(v):
                cols[k].append(float(x))
    return dict(zip(names, cols))


def edit(src, steps, job, end_time=None):
    src = re.sub(r"^JobName = .*$", 'JobName = "%s"' % job, src, flags=re.M)
    if re.search(r"^MaxSteps = ", src, flags=re.M):
        src = re.sub(r"^MaxSteps = .*$", "MaxSteps = %d" % steps, src, flags=re.M)
    else:
        src = re.sub(r"^\[Numerics\]$", "[Numerics]\nMaxSteps = %d" % steps, src, flags=re.M)
    if end_time is not None:
        src = re.sub(r"^EndTime = .*$", "EndTime = %s" % end_time, src, flags=re.M)
    return src


def to_old_vts(src, dst):
    """Rewrite a ConcZ .vts in the layout written before FilmThickness and
    Substrate were added (be91f56^, OutputVTK::Write(Field_CZ&)): no
    thickness arrays, and the film only in the point heights, every layer kz
    of the nz (ghosts included) at z = s + h kz / nz. The reader has to take
    the film back from those heights."""
    raw = open(src, "rb").read()
    mark = raw.index(b"   _") + 4
    head = raw[:mark].decode()
    nx, ny, nz = (int(re.search(r'%s="(\d+)"' % k, head).group(1)) for k in ("Nx", "Ny", "Nz"))
    npt = nx * ny * nz
    pos = mark
    blocks = []
    for count in (npt, 3 * npt, 3 * npt, npt, npt):     # Conc, Velocity, Coordinate, h, s
        size = struct.unpack_from("<I", raw, pos)[0]
        if size != 8 * count:
            raise ValueError("%s: block of %d bytes, expected %d" % (src, size, 8 * count))
        blocks.append(list(struct.unpack_from("<%dd" % count, raw, pos + 4)))
        pos += 4 + size
    conc, vel, coord, hh, ss = blocks
    for kz in range(nz):
        for col in range(nx * ny):
            p = kz * nx * ny + col
            coord[3 * p + 2] = ss[p] + hh[p] * kz / nz
    head = re.sub(r'<DataArray Name="(FilmThickness|Substrate)"[^\n]*\n', "", head)
    if "FilmThickness" in head or "Substrate" in head:
        raise ValueError("%s: thickness arrays not removed" % src)
    out = head.encode()
    for b in (conc, vel, coord):
        out += struct.pack("<I", 8 * len(b)) + struct.pack("<%dd" % len(b), *b)
    out += b"\n</AppendedData>\n</VTKFile>\n"
    with open(dst, "wb") as f:
        f.write(out)


def read_vts(path):
    """The appended arrays of a .vts, found through the offsets its header
    gives them, as {name: [values]}, and the grid size. Also checks that the
    file holds together the way ParaView reads it: each block is as long as
    its array (points x components), the blocks follow one another with no
    gap or overlap, and the last one ends where the appended data does.
    (The restart reads the blocks one after the other and never looks at the
    offsets, so a wrong offset, or a block of the wrong length at the end,
    reached no test.)"""
    raw = open(path, "rb").read()
    mark = raw.index(b"   _") + 4
    head = raw[:mark].decode()
    ext = [int(v) for v in re.search(r'WholeExtent="([^"]+)"', head).group(1).split()]
    nx, ny, nz = ext[1], ext[3], ext[5]
    npt = nx * ny * nz
    arrays, blocks, problems = {}, [], []
    for tag in re.findall(r"<DataArray ([^>]*)/>", head):
        at = dict(re.findall(r'(\w+)="([^"]*)"', tag))
        ncomp = int(at.get("NumberOfComponents", "1"))
        off = int(at["offset"])
        if mark + off + 4 > len(raw):
            problems.append("%s: offset %d is past the end of the file" % (at["Name"], off))
            continue
        size = struct.unpack_from("<I", raw, mark + off)[0]
        if size != 8 * npt * ncomp:
            problems.append("%s: block of %d bytes, %d expected" % (at["Name"], size, 8 * npt * ncomp))
            continue
        arrays[at["Name"]] = list(struct.unpack_from("<%dd" % (npt * ncomp), raw, mark + off + 4))
        blocks.append((off, 4 + size, at["Name"]))
    blocks.sort()
    pos = 0
    for off, length, name in blocks:
        if off != pos:
            problems.append("%s at offset %d, the block before it ends at %d" % (name, off, pos))
        pos = off + length
    end = raw.index(b"\n</AppendedData>", mark) - mark
    if pos != end:
        problems.append("the blocks end at %d, the appended data at %d" % (pos, end))
    return arrays, (nx, ny, nz), problems


def check_vts(vts, ref, tol=1e-13):
    """The .vts written with a field*.dat must carry the same state: the
    concentration, film, substrate, velocity and point positions of the .dat
    at the points both have (the well-mixed files have the same points; a
    ConcZ .vts has every cell, ghosts included, and its real cell centers are
    among the points of the .dat). For the well-mixed model the velocity is
    also the sum of the per-term velocities written beside it."""
    arr, (nx, ny, nz), problems = read_vts(vts)
    if problems:
        return problems

    def cmp(what, a, b):
        scale = max([abs(x) for x in b] + [1e-300])
        err = max(abs(x - y) for x, y in zip(a, b)) / scale
        if err > tol:
            problems.append("%s: .vts and .dat differ by %.2e" % (what, err))

    vel, xyz = arr["Velocity"], arr["Coordinate"]
    if "zz" not in ref:                               # well-mixed: nx*ny points, 2 layers
        n = nx * ny
        if len(ref["CC"]) != n:
            return ["%d points in the .dat, %d in a layer of the .vts" % (len(ref["CC"]), n)]
        cmp("Conc", arr["Conc"][:n], ref["CC"])
        cmp("Height", arr["Height"][:n], ref["HS"])
        cmp("Substrate", arr["Substrate"][:n], ref["SS"])
        cmp("Velocity x", vel[0:3 * n:3], ref["Ux"])
        cmp("Velocity y", vel[1:3 * n:3], ref.get("Uy", [0.0] * n))
        cmp("x", xyz[0:3 * n:3], ref["xx"])
        cmp("y", xyz[1:3 * n:3], ref.get("yy", [0.0] * n))
        cmp("z of the top layer", xyz[3 * n + 2::3], ref["HS"])
        for c in (0, 1):
            tot = [sum(arr[k][3 * p + c] for k in ("V_Ca", "V_Ma", "V_Bo", "V_Pg", "V_Sp"))
                   for p in range(2 * n)]
            cmp("Velocity %s against the sum of the terms" % "xy"[c], tot, vel[c::3])
        return problems

    #  ConcZ: the .dat walks wall face, cell centers, wall face along each axis
    #  (z outermost, then y, then x); the .vts every cell, ghosts included.
    xs = sorted(set(ref["xx"]))
    ys = sorted(set(ref["yy"])) if "yy" in ref else [0.0]
    ndx, ndy = len(xs), len(ys)
    ndz = len(ref["CC"]) // (ndx * ndy)
    gx, gy, gz = (nx - (ndx - 2)) // 2, (ny - (ndy - 2)) // 2 if ny > 1 else 0, (nz - (ndz - 2)) // 2
    a, b, c = [], [], []
    for kr in range(ndz - 2):
        for jr in range(ndy - 2 if ny > 1 else 1):
            for ir in range(ndx - 2):
                p = (kr + gz) * nx * ny + (jr + gy) * nx + (ir + gx)
                r = (kr + 1) * ndx * ndy + ((jr + 1) if ny > 1 else 0) * ndx + (ir + 1)
                a.append(p); b.append(r)
    for name, col, comp in (("Conc", "CC", None), ("FilmThickness", "HH", None),
                            ("Substrate", "SS", None), ("Velocity x", "Ux", 0),
                            ("Velocity y", "Uy", 1), ("Velocity z", "Uz", 2),
                            ("x", "xx", 0), ("y", "yy", 1), ("z", "zz", 2)):
        if col not in ref:
            continue
        if comp is None:
            got = [arr[name][p] for p in a]
        else:
            got = [(vel if name.startswith("Velocity") else xyz)[3 * p + comp] for p in a]
        cmp(name, got, [ref[col][r] for r in b])
    return problems


def run(exe, d, name, src, timeout):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name + ".toml"), "w") as f:
        f.write(src)
    try:
        r = subprocess.run([exe, name + ".toml"], cwd=d, capture_output=True,
                           text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return -1, "timed out after %g s" % timeout
    return r.returncode, r.stdout + r.stderr


def pvd_entries(d):
    p = glob.glob(os.path.join(d, "*.pvd"))
    if not p:
        return []
    return re.findall(r'timestep="([^"]+)" file="([^"]+)"', open(p[0]).read())


def maxrel(a, b):
    scale = max(abs(x) for x in a) or 1.0
    return max(abs(x - y) for x, y in zip(a, b)) / scale


#--------------------------------------------------------------------------
#  --chain: a run interrupted and restarted again and again must be the run
#  it would have been without the interruptions.
#
#  Each case is run once straight through ("U") and once as a chain of
#  parts ("C"): every part but the last is stopped by Numerics.MaxSteps at a
#  step that is no output time, and the next part restarts from the last
#  field file the part before wrote (the one written as it stopped), from
#  the .vts and the field*.dat in turn, into the same result directory, as
#  a user resuming a job would. One case restarts once from an earlier file
#  instead (a job killed some way past its last useful output), so that the
#  part after it is done twice. With the field*.state written beside each
#  field file, everything is compared to the bit: every field file U wrote
#  (the .dat as text, the arrays of the .vts byte for byte), every row of
#  History.dat, CalcInfo.dat (Newton iterations, residuals, dt) and the
#  probe files, the step at which the steady-state test ends the run, and
#  the step, Newton and linear iteration counts of JobInfo.txt summed over
#  the parts. C may only add the outputs written at the moments it stopped.
#
#  Without the state file (a field file written by hand, or one copied
#  without it) a restart takes the field alone. Two fixed-step cases check
#  that this still holds the output schedule (no burst of files after the
#  restart) and stays within a stated distance of U: the field file holds
#  the state to its digits, and the micrometers it is written in round
#  once on the way back.
#--------------------------------------------------------------------------
def toml_edit(src, sect, append=""):
    """Set keys of the given tables, by scanning line by line (a regular
    expression across tables corrupts the file); "" is the root table."""
    out, cur, seen = [], "", set()
    for line in src.split("\n"):
        st = line.strip()
        if st.startswith("[") and st.endswith("]"):
            cur = st.strip("[]")
            out.append(line)
            continue
        if "=" in st and not st.startswith("#"):
            key = st.split("=", 1)[0].strip()
            if cur in sect and key in sect[cur]:
                out.append("%s = %s" % (key, sect[cur][key]))
                seen.add(cur + "." + key)
                continue
        out.append(line)
    for name, kv in sect.items():
        add = ["%s = %s" % (k, v) for k, v in kv.items() if name + "." + k not in seen]
        if not add:
            continue
        if not name:
            out = add + out
        elif "[%s]" % name in out:
            i = out.index("[%s]" % name)
            out = out[:i + 1] + add + out[i + 1:]
        else:
            out += ["", "[%s]" % name] + add
    return "\n".join(out) + append


PROBE = "\n\n[[Output.ProbeList.Probe]]\nX = 0.0\nY = 0.0\nZ = 2.0e-6\n"


def chain_case(case, sect, end, fint, oint, parts, back=None, state=True, tol=0.0, append=""):
    """parts: the MaxSteps of every part but the last; back: the part after
    which the restart is taken from the file before the last; append: text
    added to the end of the input (tables of arrays). fint None: the field
    output follows a list (Output.FieldData.Times or LogSchedule, given in
    sect or append), and the case's TimeInterval is taken out."""
    fd = {"NeutralFile": "true", "VTKFile": "true"}
    if fint is not None:
        fd["TimeInterval"] = fint
    s = {"Output.FieldData": fd,
         "Output.GeneralHistory": {"TimeInterval": oint},
         "Output.CalcInfo": {"TimeInterval": oint},
         "Output.ProbeList": {"TimeInterval": oint},
         "Output": {"Overwrite": "true"},
         "Numerics": {"EndTime": end, "MaxSteps": "4000000000"}}
    for k, v in sect.items():
        s.setdefault(k, {}).update(v)
    return dict(case=case, sect=s, parts=parts, back=back, state=state, tol=tol, append=append,
                listed=fint is None)


LOG_SCHEDULE = ("\n\n[Output.FieldData.LogSchedule]\nStart = 1e-5\nEnd = 0.04\nPerDecade = 10\n"
                "ExactTimes = %s\n")

PINNED_DRY = {
    "Grid": {"Nx": 60},
    "Numerics.TimeStepControl": {"Enable": "false"},
    "Boundary.XminBC": {"Type": '"Pinning"'},
    "Initial.Concentration": {"Value": "0.4"},
    "Initial.SurfaceShape": {"PrecursorRatio": "0.02"},
    "Initial.SurfaceShape.EllipseBump": {"Height": "0.0"},
    "Substrate.EllipseBump": {"Height": "0.0"},
    "Physics.HorizontalGravity": {"Acceleration": "[-100.0, 0]"},
    "Physics.Evaporation": {"Coef": "4e-6"}}
PINNED_DRY_EDGE = ("\n\n[[Initial.SurfaceShape.EdgeBump]]\nOrder = 1\nCreateType = \"add\"\n"
                   "Base = 0.5\nHeight = 1.0\nLocation = \"Xmin\"\nAlpha = 1.0\n"
                   "InnerLength = 0.05e-3\nOuterLength = 0.15e-3\n")

CHAIN = {
    #  adaptive dt, evaporation, every body force, Marangoni, WM
    #  concentration (UTOPIA: the explicit half from the last Newton pass),
    #  height Koren (frozen limiter weights)
    "wm1d": chain_case("allphys_1d_reg", {"Grid": {"Nx": 60}}, "0.03", "0.01", "0.002", [97, 89]),
    #  fixed dt, limited concentration and height schemes (the explicit half
    #  rebuilt at commit), the steady-state test ending the run mid-way, and
    #  a restart from an earlier file than the last
    "wm1d_steady": chain_case("allphys_1d_reg", {
        "Grid": {"Nx": 60}, "Numerics.TimeStepControl": {"Enable": "false"},
        "Numerics.Convection": {"Concentration": '"Koren"', "Height": '"VanLeer"'},
        "Numerics": {"SteadyStateCriterion": "1e9", "SteadyStateSteps": "45"}},
        "0.1", "0.01", "0.002", [17, 19], back=1),
    #  the nozzle on its arm: rotation ramped up, held, ramped down
    "nozzle": chain_case("nozzle_1d_rotating", {
        "Grid": {"Nx": 60},
        "Physics.Source": {"RotationRate": "20.0", "RampTime": "0.003", "PlateauTime": "0.003"}},
        "0.009", "0.003", "0.001", [251, 307]),
    #  the same with the liquid supplied at its own concentration
    #  (ConcType = "supplied"): c under the nozzle moves toward Conc, and
    #  SuppliedResin goes on from the state
    "nozzle_supplied": chain_case("nozzle_1d_rotating", {
        "Grid": {"Nx": 60},
        "Physics.Source": {"RotationRate": "20.0", "RampTime": "0.003", "PlateauTime": "0.003",
                           "ConcType": '"supplied"', "Conc": "0.2"}},
        "0.009", "0.003", "0.001", [251, 307]),
    "wm2d": chain_case("allphys_2d_reg", {"Grid": {"Nx": 16, "Ny": 16}}, "0.03", "0.01", "0.005", [13, 11]),
    #  open edges: fed through a FreshInlet, drained through an Outlet,
    #  adaptive dt. What crossed them is summed step by step from the face
    #  fluxes the step put through them (OutflowedVol, OutflowedResin): the
    #  sums go on from the state, and the balance is the straight run's
    "wm1d_open": chain_case("allphys_1d_reg", {
        "Grid": {"Nx": 60}, "Boundary.XminBC": {"Type": '"FreshInlet"'},
        "Boundary.XmaxBC": {"Type": '"Outlet"'}}, "0.03", "0.01", "0.002", [97, 89]),
    "wm2d_open": chain_case("allphys_2d_reg", {
        "Grid": {"Nx": 16, "Ny": 16}, "Boundary.XminBC": {"Type": '"FreshInlet"'},
        "Boundary.XmaxBC": {"Type": '"Outlet"'}, "Boundary.YminBC": {"Type": '"Outlet"'}},
        "0.03", "0.01", "0.005", [13, 11]),
    #  CZ: predictor-corrector (c of the step before), Marangoni coupling
    #  passes, evaporation, adaptive dt
    "cz1d": chain_case("concz_1d_reg", {}, "0.03", "0.01", "0.002", [23, 41]),
    "cz2d": chain_case("concz_2d", {"Grid": {"Nx": 12, "Ny": 12, "Nz": 6},
                                    "Numerics": {"SteadyStateCriterion": "0"}},
                       "6e-4", "2e-4", "1e-4", [11, 13]),
    #  stopped (MaxSteps) at the very step at which the steady-state test
    #  ends the straight run: the next part must take no step
    "cz2d_steady": chain_case("concz_2d", {"Grid": {"Nx": 12, "Ny": 12, "Nz": 6}},
                              "6e-4", "2e-4", "1e-4", [11]),
    #  stopped (MaxSteps) at the very step at which the dry-out test ends
    #  the straight run: the next part must take no step
    #  (tools/proptest.py). A flat film well above its precursor,
    #  so that it dries, fixed dt.
    "wm1d_dryout": chain_case("allphys_1d_reg", {
        "Grid": {"Nx": 60}, "Numerics.TimeStepControl": {"Enable": "false"},
        "Physics.Evaporation": {"Coef": "1e-4", "EndConcValue": "0.6"},
        "Initial.SurfaceShape": {"PrecursorRatio": "0.1"},
        "Initial.SurfaceShape.EllipseBump": {"Base": "0.5", "Height": "0.0"}},
        "0.2", "0.01", "0.002", [56]),
    #  a pinned edge (Xmin) with the film pushed against it by the in-plane
    #  gravity and dried by a strong evaporation: its film at the wall
    #  grows past a tenth of the thickest film (0.14 at 0.4 s), which the check of the
    #  pinned edges, made on the initial condition, must not refuse on a
    #  restart (tools/proptest.py seed 1206); fixed dt, then adaptive
    "wm1d_pinned_dry": chain_case("allphys_1d_reg", PINNED_DRY, "0.6", "0.1", "0.02", [140, 330],
                                  append=PINNED_DRY_EDGE),
    "wm1d_pinned_dry_adapt": chain_case("allphys_1d_reg", dict(PINNED_DRY, **{
        "Numerics.TimeStepControl": {"Enable": "true"}}), "0.6", "0.1", "0.02", [700, 330],
        append=PINNED_DRY_EDGE),
    #  field output at listed times (Output.FieldData.LogSchedule): the
    #  list goes on where the part before left it, each time landed on
    #  (ExactTimes; the step after a landing goes back to the step before
    #  it) or written at the step closest to it; adaptive dt
    "wm1d_listed": chain_case("allphys_1d_reg", {"Grid": {"Nx": 60}}, "0.03", None, "0.002", [97, 89],
                              append=LOG_SCHEDULE % "true"),
    "wm1d_listed_near": chain_case("allphys_1d_reg", {"Grid": {"Nx": 60}}, "0.03", None, "0.002", [97, 89],
                                   append=LOG_SCHEDULE % "false"),
    "cz1d_listed": chain_case("concz_1d_reg", {}, "0.03", None, "0.002", [23, 41],
                              append=LOG_SCHEDULE % "true"),
    #  the driving potential written (WritePotential): one more column of
    #  field*.dat, read past by the ASCII restart, and one more array of
    #  the .vts after the coordinates
    "wm1d_potential": chain_case("allphys_1d_reg", {"Grid": {"Nx": 60},
                                 "Output.FieldData": {"WritePotential": "true"}},
                                 "0.03", "0.01", "0.002", [97, 89]),
    "wm2d_potential": chain_case("allphys_2d_reg", {"Grid": {"Nx": 16, "Ny": 16},
                                 "Output.FieldData": {"WritePotential": "true"}},
                                 "0.03", "0.01", "0.005", [13, 11]),
    #  Output.Energy: Energy.dat continued to the bit, its dissipated
    #  energy and the balance of the interval across the restart from the
    #  state file
    "wm1d_energy": chain_case("allphys_1d_reg", {"Grid": {"Nx": 60}}, "0.03", "0.01", "0.002", [97, 89],
                              append="\n[Output.Energy]\nEnable = true\n"),
    "wm2d_energy": chain_case("allphys_2d_reg", {"Grid": {"Nx": 16, "Ny": 16}}, "0.03", "0.01", "0.005", [13, 11],
                              append="\n[Output.Energy]\nEnable = true\n"),
    #  ... and without the state file: from the first listed time after
    #  the restart
    "wm1d_listed_nostate": chain_case("allphys_1d_reg", {
        "Grid": {"Nx": 60}, "Numerics.TimeStepControl": {"Enable": "false"}},
        "0.04", None, "0.002", [13, 17], state=False, tol=1e-11, append=LOG_SCHEDULE % "true"),
    #  no state file: the field alone, to the digits of the file
    "wm1d_nostate": chain_case("allphys_1d_reg", {
        "Grid": {"Nx": 60}, "Numerics.TimeStepControl": {"Enable": "false"}},
        "0.04", "0.01", "0.002", [13, 17], state=False, tol=1e-11),
    "cz1d_nostate": chain_case("concz_1d_reg", {"Numerics.TimeStepControl": {"Enable": "false"}},
        "0.04", "0.01", "0.002", [13, 17], state=False, tol=1e-7),
}

JOBINFO_SUMS = ("Time steps", "Newton iterations", "Linear iterations",
                "Film-concentration passes taken again", "Steps failed: the passes did not settle",
                "Steps failed: the Newton correction stagnated")


def jobinfo_counts(path):
    d = {}
    for line in open(path):
        k, _, v = line.partition(":")
        if k.strip() in JOBINFO_SUMS:
            d[k.strip()] = int(v.split()[0])
    return d


def timeline_rows(path):
    """{time string: row} of a History/CalcInfo/Probe file, and the times in
    the order written (a time twice means an output written twice)."""
    rows, order = {}, []
    for line in open(path):
        tok = line.split("\t", 1)[0].strip()
        try:
            float(tok)
        except ValueError:
            continue
        rows[tok] = line.rstrip("\n")
        order.append(tok)
    return rows, order


def vts_payload(path):
    raw = open(path, "rb").read()
    return raw[raw.index(b"   _") + 4:]


def run_chain(exe, workdir, name, v, timeout):
    """Run U and C; return the list of problems and a one-line summary."""
    base = open(os.path.join(CASES, v["case"] + ".toml")).read()
    d = os.path.join(workdir, "chain", name)

    def src(job, maxsteps, restart):
        s = {k: dict(x) for k, x in v["sect"].items()}
        s.setdefault("", {})["JobName"] = '"%s"' % job
        if maxsteps:
            s["Numerics"]["MaxSteps"] = str(maxsteps)
        if restart:
            s["Restart"] = restart
        b = base
        if v.get("listed"):
            b = re.sub(r"(\[Output\.FieldData\]\n)TimeInterval = [^\n]*\n", r"\1", b, count=1)
        return toml_edit(b, s, PROBE + v.get("append", ""))

    problems = []
    rc, log = run(exe, os.path.join(d, "U"), "in", src("run", None, None), timeout)
    if rc != 0:
        return ["straight run: rc=%d %s" % (rc, log[-400:])], ""
    ud = os.path.join(d, "U", "run")

    #  The chain. With the state file every part writes into the same
    #  directory; without it each part has its own (as before, a restart
    #  without the state starts its History.dat afresh), merged by time.
    cdir = os.path.join(d, "C")
    os.makedirs(cdir, exist_ok=True)
    extra = set()           # times C may write and U not: the stops
    counts = {}
    segdirs = []
    restart = None
    nparts = len(v["parts"]) + 1
    for k in range(nparts):
        job = "run" if v["state"] else "run%d" % k
        last = k == nparts - 1
        rc, log = run(exe, cdir, "in%d" % k, src(job, None if last else v["parts"][k], restart), timeout)
        if rc != 0:
            return problems + ["part %d: rc=%d %s" % (k, rc, log[-400:])], ""
        rd = os.path.join(cdir, job)
        segdirs.append(rd)
        if k > 0 and v["state"] and ("continuing exactly" not in log or "taken to the bit" not in log):
            problems.append("part %d did not continue exactly from the state file" % k)
        if counts is not None:
            for key, val in jobinfo_counts(os.path.join(rd, "JobInfo.txt")).items():
                counts[key] = counts.get(key, 0) + val
        if last:
            break
        ent = pvd_entries(rd)
        extra.add(ent[-1][0])
        pick = ent[-1]
        if v["back"] is not None and k == v["back"]:
            pick = ent[-2]
            extra.discard(ent[-1][0])
        t, f = pick
        idx = int(re.search(r"(\d+)\.vts$", f).group(1))
        ascii_ = k % 2 == 1
        #  "/" whatever the OS: the name goes into a TOML basic string, where
        #  the backslash of a Windows path starts an escape ("run" and
        #  "field00002.vts" joined by one read as "run", a form feed, "ield...").
        fname = job + "/" + f[:-4] + (".dat" if ascii_ else ".vts")
        if v["state"]:
            restart = {"Enable": "true", "ASCII": "true" if ascii_ else "false",
                       "FileName": '"%s"' % fname, "RestartTime": t, "FileIndex": str(idx + 1)}
        else:
            #  the field file alone, in a directory of its own
            sd = os.path.join(cdir, "start%d" % k)
            os.makedirs(sd, exist_ok=True)
            shutil.copy(os.path.join(cdir, fname), sd)
            restart = {"Enable": "true", "ASCII": "true" if ascii_ else "false",
                       "FileName": '"start%d/%s"' % (k, os.path.basename(fname)),
                       "RestartTime": t, "FileIndex": str(idx + 1)}
            extra.add(t)        # written again by the restarted part
        if v["back"] is not None and k == v["back"]:
            counts = None       # a part done twice: the sums mean nothing

    def is_stop(tk):
        """A time C may have written and U not: one of the stops, written
        to six digits in a .pvd and to seven in the timeline files."""
        return any(abs(float(tk) - float(e)) <= 5e-6 * abs(float(e)) for e in extra)

    #  Field files, by time
    ue = pvd_entries(ud)
    ce = {}
    for rd in segdirs:
        for t, f in pvd_entries(rd):
            ce.setdefault(t, os.path.join(rd, f))
    ut = [t for t, _ in ue]
    missing = [t for t in ut if t not in ce and float(t) >= 0]
    surplus = [t for t in ce if t not in ut and not is_stop(t)]
    if missing:
        problems.append("field times of U missing in C: %s" % missing[:5])
    if surplus:
        problems.append("C wrote field files U did not, at %s" % surplus[:5])
    worst = 0.0
    for t, f in ue:
        if t not in ce:
            continue
        a, b = read_field(os.path.join(ud, f[:-4] + ".dat")), read_field(ce[t][:-4] + ".dat")
        if v["state"]:
            if open(os.path.join(ud, f[:-4] + ".dat")).read() != open(ce[t][:-4] + ".dat").read():
                problems.append("t = %s: field*.dat differs" % t)
            if vts_payload(os.path.join(ud, f)) != vts_payload(ce[t]):
                problems.append("t = %s: .vts arrays differ" % t)
        for col in a:
            if col in b and any(a[col]):
                worst = max(worst, maxrel(a[col], b[col]))
    if not v["state"] and worst > v["tol"]:
        problems.append("fields differ by %.2e, more than %.0e" % (worst, v["tol"]))

    #  History, CalcInfo, probes: every row of U, the same, and nothing else
    #  but the stops
    hist = {}
    for fn in sorted(os.path.basename(p) for p in glob.glob(os.path.join(ud, "*.dat"))
                     if not os.path.basename(p).startswith("field")):
        ur, _ = timeline_rows(os.path.join(ud, fn))
        cr, corder = {}, []
        for rd in segdirs:
            r, o = timeline_rows(os.path.join(rd, fn))
            for tk in o:
                if tk not in cr:
                    corder.append(tk)
                cr.setdefault(tk, r[tk])
            if v["state"]:
                break       # one directory, one file
        if v["state"] and len(corder) != len(set(corder)):
            problems.append("%s: a time written twice" % fn)
        diff = [tk for tk in ur if tk in cr and cr[tk] != ur[tk]]
        miss = [tk for tk in ur if tk not in cr]
        surp = [tk for tk in cr if tk not in ur and not is_stop(tk)]
        if v["state"] and diff:
            problems.append("%s: %d rows differ, first at t = %s" % (fn, len(diff), diff[0]))
        if miss:
            problems.append("%s: rows of U missing at %s" % (fn, miss[:3]))
        if surp:
            problems.append("%s: rows U does not have, at %s" % (fn, surp[:3]))
        if fn == "History.dat" and diff:
            #  how far the totals are off without the state file
            cols = ur[diff[-1]].split("\t")
            got = cr[diff[-1]].split("\t")
            hist = {k: (cols[k], got[k]) for k in range(1, len(cols)) if cols[k] != got[k]}

    if v["state"] and counts is not None:
        ucnt = jobinfo_counts(os.path.join(ud, "JobInfo.txt"))
        for key, val in ucnt.items():
            if counts.get(key) != val:
                problems.append("JobInfo %s: %s over the parts, %d straight" % (key, counts.get(key), val))
    summary = "%d field times, %s" % (len(ue), "bit-identical" if v["state"] else "max rel. diff %.1e" % worst)
    if hist:
        summary += "; History columns differing at the end: %d" % len(hist)
    return problems, summary


def main_chain(args):
    exe = os.path.abspath(args.exe)
    shutil.rmtree(os.path.join(args.workdir, "chain"), ignore_errors=True)
    bad = 0
    for name, v in CHAIN.items():
        if args.only and name not in args.only:
            continue
        problems, summary = run_chain(exe, args.workdir, name, v, args.run_timeout)
        if problems:
            print("FAIL %-14s %s" % (name, "; ".join(problems)))
            bad += 1
        else:
            print("ok   %-14s %s" % (name, summary))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--run-timeout", type=float, default=600,
                    help="seconds allowed for any one run of the solver")
    ap.add_argument("--tol", type=float, default=1e-12,
                    help="relative tolerance on thickness, height and concentration")
    ap.add_argument("--chain", action="store_true",
                    help="restart a run again and again and compare it with one run straight through")
    ap.add_argument("--only", action="append", help="with --chain: run only this case (repeatable)")
    args = ap.parse_args()
    if args.chain:
        return main_chain(args)
    exe = os.path.abspath(args.exe)
    shutil.rmtree(args.workdir, ignore_errors=True)

    bad = 0
    for case, steps in PLAN:
        base = open(os.path.join(CASES, case + ".toml")).read()
        d0 = os.path.join(args.workdir, case, "source")
        rc, log = run(exe, d0, "src", edit(base, steps, "src"), args.run_timeout)
        fields = sorted(glob.glob(os.path.join(d0, "src", "field*.dat")))
        if rc != 0 or not fields:
            print("FAIL %-18s source run: rc=%d\n%s" % (case, rc, log[-800:]))
            bad += 1
            continue
        last = fields[-1]
        t_src = pvd_entries(os.path.join(d0, "src"))[-1][0]
        idx = int(re.search(r"(\d+)\.dat$", last).group(1))
        ref = read_field(last)

        #  The two files the run wrote of the same state must agree.
        problems = check_vts(last[:-4] + ".vts", ref)
        if problems:
            print("FAIL %-18s .vts against .dat: %s" % (case, "; ".join(problems)))
            bad += 1
        else:
            print("ok   %-18s .vts and .dat of %s agree" % (case, os.path.basename(last)))

        #  "old": a ConcZ .vts without the thickness array, the film to be
        #  recovered from the point heights.
        kinds = ["vts", "dat"] + (["old"] if case.startswith("concz") else [])
        for kind in kinds:
            ascii_ = kind == "dat"
            ext = ".dat" if ascii_ else ".vts"
            d = os.path.join(args.workdir, case, kind)
            os.makedirs(d, exist_ok=True)
            if kind == "old":
                to_old_vts(last[:-4] + ext, os.path.join(d, "start" + ext))
            else:
                shutil.copy(last[:-4] + ext, os.path.join(d, "start" + ext))
            src = edit(base, 3, "rst", end_time="1.0")
            blk = ("[Restart]\nEnable = true\nASCII = %s\nFileName = \"start%s\"\n"
                   "RestartTime = %s\nFileIndex = %d\n"
                   % ("true" if ascii_ else "false", ext, t_src, idx + 1))
            if re.search(r"^\[Restart\]$", src, flags=re.M):
                src = re.sub(r"^\[Restart\]\n(?:[A-Za-z].*\n)*", blk, src, count=1, flags=re.M)
            else:
                src += "\n" + blk
            rc, log = run(exe, d, "rst", src, args.run_timeout)
            outs = sorted(glob.glob(os.path.join(d, "rst", "field*.dat")))
            label = "%-18s %s" % (case, kind)
            if rc != 0 or not outs:
                print("FAIL %s restart run: rc=%d\n%s" % (label, rc, log[-800:]))
                bad += 1
                continue

            got = read_field(outs[0])
            problems = []
            if os.path.basename(outs[0]) != "field%05d.dat" % (idx + 1):
                problems.append("first file %s" % os.path.basename(outs[0]))
            for col in ("HH", "SS", "HS", "CC") + (("zz",) if "zz" in ref else ()):
                if len(got[col]) != len(ref[col]):
                    problems.append("%s: %d rows, expected %d" % (col, len(got[col]), len(ref[col])))
                    continue
                if any(ref[col]) and maxrel(ref[col], got[col]) > args.tol:
                    problems.append("%s differs by %.2e" % (col, maxrel(ref[col], got[col])))
            ent = pvd_entries(os.path.join(d, "rst"))
            if not ent or ent[0][1] != "field%05d.vts" % (idx + 1):
                problems.append(".pvd starts at %s" % (ent[0][1] if ent else "nothing"))
            elif abs(float(ent[0][0]) - float(t_src)) > 1e-9 * max(1.0, abs(float(t_src))):
                problems.append(".pvd time %s, restarted from %s" % (ent[0][0], t_src))
            if len(ent) < 2 or float(ent[-1][0]) <= float(t_src):
                problems.append("the restarted run did not advance")

            if problems:
                print("FAIL %s: %s" % (label, "; ".join(problems)))
                bad += 1
            else:
                print("ok   %s  restarted at %s from %s" % (label, t_src, os.path.basename(last)))

    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
