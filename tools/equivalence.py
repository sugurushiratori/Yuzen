#!/usr/bin/env python3
"""Check that the solver gives the same answer to the same problem posed twice.

    python3 tools/equivalence.py --exe build/bin/yuzen --workdir /tmp/eq
    python3 tools/equivalence.py --exe build/bin/yuzen --workdir /tmp/eq --only wm2d_transpose
    python3 tools/equivalence.py --exe build/bin/yuzen --workdir /tmp/eq --scheme Koren

(--scheme, --conc and --height put a [Numerics.Convection] block into every
run; see tools/mms.py.)

Each check runs a problem and a copy of it that differs only in how it is
written down, and compares the two results point by point:

    wm2d_transpose   x and y exchanged: grid, shapes, gravity, rotation axis
                     and edges. h(x, y) must come back as h'(y, x), Ux as Uy'.
    wm2d_mirror_x    x -> -x. With the transpose this covers y -> -y as well.
    wm2d_nozzle_transpose  the transpose with a fixed nozzle off the center.
    wm2d_nozzle_thin  the same transpose on cells 5.6 times
                     longer one way than the other.
    wm2d_nozzle_mirror_x   the mirror with a fixed nozzle off the center.
    wm2d_spin_transpose, wm2d_spin_mirror_x  the transpose and the mirror
                     with a strong rotation and UTOPIA for the film height:
                     the upwind side of the centrifugal flux of a fixed
                     scheme (the default, Koren, chooses it elsewhere).
    wm_1d_vs_2d      a 1-D run against a 2-D run of the same problem made
                     uniform in y (periodic in y, shapes that do not vary in
                     y). The 1-D and 2-D well-mixed models are separate code.
    wm_1d_vs_2d_koren  the same with the Koren scheme for the concentration
                     (a limited scheme rebuilds the explicit half of each
                     step after the evaporation is taken out of it).
    wm_io_1d_vs_2d   the same with a FreshInlet and an Outlet in x, y periodic
    wm_io_1d_vs_2d_sym  and with y symmetric (the inlet/outlet corners).
    wm_io_1d_vs_2d_ymin  the 1-D inlet run against a 2-D one flowing along
                     y from an inlet on Ymin, x periodic.
    wm_io_1d_vs_2d_ymax_sym  and along -y from an inlet on Ymax, x symmetric.
    wm1d_io_mirror   the 1-D inlet run against its mirror, the inlet at Xmax.
    wm1d_shift       a periodic 1-D film moved round the domain by a whole
                     number of cells: the result must move with it.
    wm2d_shift       the same in 2-D, moved in x and y at once.
    wm1d_scales      the same physical problem made dimensionless with other
                     scales: TotalHeight 2.5 times larger (so another H0) and
                     a custom length scale in place of the capillary length.
    cz1d_hscale      ConcZ with another H0.
    cz1d_lscale      ConcZ with another L0.
    wm1d_nozzle_mirror  a nozzle off the center of a flat 1-D film against
                     the same nozzle on the other side.
    cz1d_mirror      ConcZ, x -> -x.
    cz1d_shift       ConcZ, a periodic film moved round the domain.
    cz2d_transpose   ConcZ, x and y exchanged.
    cz2d_dry_transpose, cz2d_dry_mirror_x  a thin ConcZ film drying with
                     Marangoni, transposed and mirrored, to 1e-9 of the
                     change: the concentration varies through the film and
                     in the plane, and the film and the solute are coupled
                     by the passes of the Marangoni stress (a
                     non-conservative form differs by 6e-8 to 8e-8).
    cz_1d_vs_2d      ConcZ 1-D against ConcZ 2-D, uniform in y.
    cz2d_deficit_transpose  a ConcZ column next to a wall restarted with
                     less than no solute and filled from its neighbors,
                     against the transpose.

and six against a known answer, the other model or another run:

    body_velocity    on a flat film the in-plane gravity and centrifugal
                     velocities are rho g h^2/(3 mu) and rho w^2 r h^2/(3 mu),
                     along Acceleration and away from the axis (1-D and 2-D).
    cz_solute        a drying flat ConcZ film keeps its solute.
    cz_solute_slope  a drying ConcZ film on a sloping substrate, flowing in
                     the plane, with Marangoni (1-D and 2-D), keeps its
                     solute to 1e-12 (it lost 3 to 25 % of what the
                     evaporated fraction should have concentrated).
    cz_uniform       a uniform concentration in a ConcZ film that flows and
                     levels (1-D and 2-D, no evaporation) stays uniform to
                     1e-13 (it drifted by ~1e-6).
    cz_wm_advection  with the solute uniform through the thickness, a strong
                     in-plane flow changes the depth-mean concentration the
                     same way in ConcZ as in well-mixed.
    wm_retry         a step rejected by the time-step control and retried
                     at a smaller dt, against that step taken at that dt.

Any term that treats one direction differently from the other, an edge or
corner that is handled on one side only, an index slip in a stencil or an
output writer, shows up as a difference far above rounding. What is left
when everything is right is rounding: the two runs assemble their matrices
in a different order, and the iterative solvers stop at their tolerances on
different iterates.

The difference is reported against the change of the field over the run,
since that is what a wrong term would disturb -- a difference of 1e-12 of h
means little if h moved by 1e-12. Fields that do not change (the substrate)
are reported against their size.

The time step is fixed so that both runs visit the same times; with adaptive
stepping, rounding in the error estimate could pick different steps and the
comparison would measure that instead.

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
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(os.path.dirname(HERE), "cases")
sys.path.insert(0, HERE)
from concz_vs_wm import edit_toml  # noqa: E402  (section-aware key rewriting)


# ----------------------------------------------------------------------
#  Input files
# ----------------------------------------------------------------------
def case(name):
    with open(os.path.join(CASES, name + ".toml")) as f:
        return f.read()


def drop_section(src, name):
    """Remove the table [name] and its keys (sub-tables are kept)."""
    out, skip = [], False
    for line in src.split("\n"):
        st = line.strip()
        if st.startswith("[") and st.endswith("]"):
            skip = (st[1:-1] == name)
        if not skip:
            out.append(line)
    return "\n".join(out)


def add_section(src, name, kv, after):
    """Insert [[name]] (one more element of an array of tables) just before
    the table that follows [after]."""
    body = "[[%s]]\n%s\n" % (name, "\n".join("%s = %s" % x for x in kv.items()))
    lines = src.split("\n")
    seen = False
    for n, line in enumerate(lines):
        st = line.strip()
        if st == "[%s]" % after:
            seen = True
            continue
        if seen and st.startswith("[") and st.endswith("]"):
            lines.insert(n, body)
            return "\n".join(lines)
    return src + "\n" + body


def arr(*v):
    return "[%s]" % ", ".join(repr(float(x)) for x in v)


def fixed_steps(src, dt, steps, interval):
    """Fixed dt, a set number of steps, a field written every `interval` s.

    Both equations are put on the "eigen" backend. "auto" would pick "amgcl"
    for the 2-D runs of an AMGCL build and "eigen" for the 1-D ones, and a
    comparison of the two would then measure the backends' tolerances rather
    than the discretization.
    """
    return edit_toml(src, {}, {
        "Numerics.HeightEquation.LinearSolver": {"Backend": '"eigen"'},
        "Numerics.ConcEquation.LinearSolver": {"Backend": '"eigen"'},
        "Numerics": {"InitialDt": repr(dt), "EndTime": "1.0e3",
                     "MaxSteps": str(steps), "SteadyStateCriterion": "-1"},
        "Numerics.TimeStepControl": {"Enable": "false"},
        "Output.FieldData": {"TimeInterval": repr(interval), "NeutralFile": "true",
                             "VTKFile": "false"},
        "Output.GeneralHistory": {"TimeInterval": repr(interval)},
        "Output.CalcInfo": {"TimeInterval": repr(interval)},
        "Output.ProbeList": {"TimeInterval": repr(interval)},
    })


def bcs(xmin, xmax, ymin, ymax):
    return {"Boundary.XminBC": {"Type": '"%s"' % xmin},
            "Boundary.XmaxBC": {"Type": '"%s"' % xmax},
            "Boundary.YminBC": {"Type": '"%s"' % ymin},
            "Boundary.YmaxBC": {"Type": '"%s"' % ymax}}


# ----------------------------------------------------------------------
#  Field files and point maps
# ----------------------------------------------------------------------
def read_field(path):
    """A field*.dat as (column dict, nx, ny, nz); rows run x, then y, then z."""
    with open(path) as f:
        head = f.readline().lstrip("#").split()
        names = [c.split("[")[0] for c in head]
        cols = {n: [] for n in names}
        for line in f:
            v = line.split()
            if len(v) != len(names):
                continue
            for n, x in zip(names, v):
                cols[n].append(float(x))
    nx = len(set(cols["xx"]))
    ny = len(set(cols["yy"])) if "yy" in cols else 1
    nz = len(cols["xx"]) // (nx * ny)
    return cols, nx, ny, nz


def last_field(d):
    fs = sorted(glob.glob(os.path.join(d, "field*.dat")))
    return fs[0], fs[-1]


def compare(A, B, A0, point_map, col_map, centers_only=False):
    """max |a - T(b)| over the points of A, relative to how far A moved.

    point_map(i, j, k, dims) gives the point of B that point (i, j, k) of A
    corresponds to; col_map[name] = (name in B, sign). A field that did not
    move over the run (the substrate) is measured against its size instead,
    and one that is zero in both runs reports 0.
    """
    ca, nx, ny, nz = A
    cb, mx, my, mz = B
    c0 = A0[0]
    res = {}
    for name, (bname, sign) in col_map.items():
        worst, change, size = 0.0, 0.0, 0.0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if centers_only and (i in (0, nx - 1) or (ny > 1 and j in (0, ny - 1))):
                        continue
                    r = (k * ny + j) * nx + i
                    i2, j2, k2 = point_map(i, j, k, (nx, ny, nz))
                    r2 = (k2 * my + j2) * mx + i2
                    a, b = ca[name][r], sign * cb[bname][r2]
                    worst = max(worst, abs(a - b))
                    change = max(change, abs(a - c0[name][r]))
                    size = max(size, abs(a), abs(b))
        if change > 1e-12 * size:
            res[name] = (worst / change, "change")
        elif size > 0:
            res[name] = (worst / size, "size")
        else:
            res[name] = (0.0, "zero")
    return res


def coords_agree(A, B, point_map, xy_map, tol=1e-9):
    """The mapped points must sit where the transformation says (in mm).
    xy_map gives the expected (x, y) in B, y None where it is not checked."""
    ca, nx, ny, nz = A
    cb, mx, my, mz = B
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                r = (k * ny + j) * nx + i
                i2, j2, k2 = point_map(i, j, k, (nx, ny, nz))
                r2 = (k2 * my + j2) * mx + i2
                x, y = ca["xx"][r], ca.get("yy", [0.0] * len(ca["xx"]))[r]
                X, Y = xy_map(x, y)
                x2, y2 = cb["xx"][r2], cb.get("yy", [0.0] * len(cb["xx"]))[r2]
                if abs(X - x2) > tol or (Y is not None and abs(Y - y2) > tol):
                    return "point (%d,%d,%d) at x = %g, y = %g maps to x = %g, y = %g, " \
                           "expected x = %g%s" % (i, j, k, x, y, x2, y2, X,
                                                  "" if Y is None else ", y = %g" % Y)
    return None


# ----------------------------------------------------------------------
#  The checks
# ----------------------------------------------------------------------
WM_COLS = ["HS", "HH", "SS", "CC", "Ux", "Uy"]
CZ_COLS = ["zz", "CC", "Ux", "Uy", "Uz", "HS", "HH", "SS"]


def same(names):
    return {n: (n, 1) for n in names}


def transposed(names):
    m = same(names)
    if "Ux" in m:
        m["Ux"], m["Uy"] = ("Uy", 1), ("Ux", 1)
    return m


def mirrored_x(names):
    m = same(names)
    m["Ux"] = ("Ux", -1)
    return m


def wm2d_base():
    """allphys_2d_reg on a grid with dx != dy and x, y ranges that differ, so
    that a slip between the two axes cannot hide behind a square domain."""
    s = case("allphys_2d_reg")
    return edit_toml(s, {}, {"Grid": {"Nx": "40", "Ny": "24",
                                      "Xmin": "-0.5e-3", "Xmax": "0.5e-3",
                                      "Ymin": "-0.3e-3", "Ymax": "0.45e-3"}})


def wm2d_transpose():
    a = fixed_steps(wm2d_base(), 2e-4, 40, 4e-3)
    b = edit_toml(a, {}, {
        "Grid": {"Nx": "24", "Ny": "40", "Xmin": "-0.3e-3", "Xmax": "0.45e-3",
                 "Ymin": "-0.5e-3", "Ymax": "0.5e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(0.2, 0.3), "Width": arr(0.5, 1.0)},
        "Substrate.EllipseBump": {"Center": arr(0.2, 0.3), "Width": arr(0.5, 1.0)},
        "Physics.HorizontalGravity": {"Acceleration": arr(0.2, 0.4903)},
        "Physics.Centrifugal": {"RotationAxis": arr(0.02e-3, 0.05e-3)},
    })
    a = edit_toml(a, {}, {"Physics.Centrifugal": {"RotationAxis": arr(0.05e-3, 0.02e-3)}})
    return (a, b, lambda i, j, k, d: (j, i, k), lambda x, y: (y, x),
            transposed(WM_COLS), False)


def wm2d_mirror_x():
    a = fixed_steps(wm2d_base(), 2e-4, 40, 4e-3)
    a = edit_toml(a, {}, {"Grid": {"Xmin": "-0.5e-3", "Xmax": "0.4e-3", "Nx": "36"},
                          "Physics.Centrifugal": {"RotationAxis": arr(0.05e-3, 0.02e-3)}})
    b = edit_toml(a, {}, {
        "Grid": {"Xmin": "-0.4e-3", "Xmax": "0.5e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(-0.3, 0.2)},
        "Substrate.EllipseBump": {"Center": arr(-0.3, 0.2)},
        "Physics.HorizontalGravity": {"Acceleration": arr(-0.4903, 0.2)},
        "Physics.Centrifugal": {"RotationAxis": arr(-0.05e-3, 0.02e-3)},
    })
    return (a, b, lambda i, j, k, d: (d[0] - 1 - i, j, k), lambda x, y: (-x, y),
            mirrored_x(WM_COLS), False)


def wm2d_nozzle_transpose():
    """wm2d_transpose with a fixed nozzle off the center (and off the
    diagonal), which none of the other checks has in 2-D."""
    a, b, pmap, xymap, cols, centers = wm2d_transpose()
    src = {"Enable": "true", "Radius": "0.12e-3", "InletRate": "0.05",
           "TimeConstant": "0.0", "ArmLength": "0.0", "SubRegions": "4"}
    a = edit_toml(a, {}, {"Physics.Source": dict(src, ArmCenter=arr(0.1e-3, -0.05e-3))})
    b = edit_toml(b, {}, {"Physics.Source": dict(src, ArmCenter=arr(-0.05e-3, 0.1e-3))})
    return a, b, pmap, xymap, cols, centers


def wm2d_nozzle_thin():
    """wm2d_nozzle_transpose on cells 5.6 times longer in one direction than
    in the other (60 x 8), the nozzle's area found from 10 x 10 points a
    cell: the cells the nozzle covers are found by a test on the distance
    from the nozzle in units of dx and of dy, which on near-square cells
    cannot miss one whatever it does."""
    a, b, pmap, xymap, cols, centers = wm2d_nozzle_transpose()
    a = edit_toml(a, {}, {"Grid": {"Nx": "60", "Ny": "8"}, "Physics.Source": {"SubRegions": "10"}})
    b = edit_toml(b, {}, {"Grid": {"Nx": "8", "Ny": "60"}, "Physics.Source": {"SubRegions": "10"}})
    return a, b, pmap, xymap, cols, centers


def wm2d_nozzle_mirror_x():
    """wm2d_mirror_x with a fixed nozzle off the center."""
    a, b, pmap, xymap, cols, centers = wm2d_mirror_x()
    src = {"Enable": "true", "Radius": "0.12e-3", "InletRate": "0.05",
           "TimeConstant": "0.0", "ArmLength": "0.0", "SubRegions": "4"}
    a = edit_toml(a, {}, {"Physics.Source": dict(src, ArmCenter=arr(0.1e-3, -0.05e-3))})
    b = edit_toml(b, {}, {"Physics.Source": dict(src, ArmCenter=arr(-0.1e-3, -0.05e-3))})
    return a, b, pmap, xymap, cols, centers


def spin_utopia(make):
    """make() with a rotation that dominates the in-plane forces and the
    fixed third-order upwinding (UTOPIA) for the film height. The default
    scheme (Koren) takes the upwind side of each face from the drive on the
    face itself; the fixed schemes take it from the mean of the two cell
    centers (CentrifugalKernel::Drive), which no other check reaches in 2-D.
    An upwind side picked wrongly is an error of third order, below every
    tolerance of a comparison with a known answer, but it breaks the
    symmetry of the problem outright."""
    def f():
        a, b, pmap, xymap, cols, centers = make()
        cv = '\n[Numerics.Convection]\nHeight = "UTOPIA"\n'
        a, b = (edit_toml(s, {}, {"Physics.Centrifugal": {"RotationRate": "60.0"}}) + cv for s in (a, b))
        return a, b, pmap, xymap, cols, centers
    return f


def uniform_in_y(src, ny=4):
    """A 1-D case as a 2-D one that does not vary in y."""
    return edit_toml(src, {}, dict({
        "Grid": {"Dimension": '"xyh"', "Ny": str(ny),
                 "Ymin": "-0.1e-3", "Ymax": "0.1e-3"},
        #  A width of 1e30 m takes y out of the radius exactly: (y/W)^2 is
        #  below the last bit of the x part everywhere in the domain.
        "Initial.SurfaceShape.EllipseBump": {"Width": arr(1.0, 1e30)},
        "Substrate.EllipseBump": {"Width": arr(1.0, 1e30)},
    }, **bcs("Symmetry", "Symmetry", "Periodic", "Periodic")))


def wm_1d_vs_2d():
    s = edit_toml(case("allphys_1d_reg"), {}, {
        "Grid": {"Nx": "100"},
        #  The centrifugal force has a y component wherever y != axis, so a
        #  run with it cannot be uniform in y; it is covered by the 2-D checks.
        "Physics.Centrifugal": {"Enable": "false"},
    })
    a = fixed_steps(s, 2e-4, 50, 5e-3)
    b = uniform_in_y(a)
    cols = same(["HS", "HH", "SS", "CC", "Ux"])
    return (a, b, lambda i, j, k, d: (i, 1, k), lambda x, y: (x, None),
            cols, False)


def wm_1d_vs_2d_koren():
    """wm_1d_vs_2d with the Koren scheme for the concentration. A limited
    scheme rebuilds the explicit half of the next step after the step (the
    limiter of c^n), when the evaporation of the step has been taken out of
    the operator: only this path reaches that rebuild with the solvent
    evaporating, and the 1-D and 2-D models build it separately."""
    a, b, pmap, xymap, cols, centers = wm_1d_vs_2d()
    kv = {"Numerics.Convection": {"Concentration": '"Koren"'}}
    return (edit_toml(a, {}, kv), edit_toml(b, {}, kv), pmap, xymap, cols, centers)


def io_1d():
    """allphys_1d_reg fed through a "FreshInlet" at Xmin and drained through
    an "Outlet" at Xmax, the film carried across by in-plane gravity. The
    solvent evaporates, so the film inside becomes richer in solute than
    what the inlet supplies (0.5, the initial value), and a wrong inlet
    concentration shows in the interior. The film itself (not only its
    surface) thickens across the inlet -- the rim of a bump outside the
    domain, about 0.18 mm wide, is centered on Xmin -- so the thickness the
    inlet holds depends on which points it is read from."""
    s = edit_toml(case("allphys_1d_reg"), {}, dict({
        "Grid": {"Nx": "100"},
        "Initial.SurfaceShape.EllipseBump": {"CreateType": '"add"', "Center": arr(-0.5005, 0.0),
                                             "SlopeWidth": "4e-4"},
        "Physics.Centrifugal": {"Enable": "false"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(2.0e3, 0.0)},
    }, **bcs("FreshInlet", "Outlet", "Symmetry", "Symmetry")))
    return fixed_steps(s, 2e-4, 50, 5e-3)


def io_1d_vs_2d(yedge):
    """wm_1d_vs_2d with a "FreshInlet" and an "Outlet" in x, the film carried
    across by in-plane gravity, and y periodic or symmetric: the edges and
    the corners they make with y must not make the 2-D run differ from the
    1-D one."""
    def make():
        a = io_1d()
        b = edit_toml(uniform_in_y(a), {}, bcs("FreshInlet", "Outlet", yedge, yedge))
        return (a, b, lambda i, j, k, d: (i, 1 if yedge == "Periodic" else 2, k),
                lambda x, y: (x, None), same(["HS", "HH", "SS", "CC", "Ux"]), False)
    return make


def io_1d_vs_2d_y(inlet, xedge):
    """The 1-D inlet/outlet run against a 2-D one in which the film flows
    along y instead, uniform in x: from an inlet on Ymin (the 1-D problem
    transposed) or on Ymax (transposed and mirrored), with x periodic or
    symmetric. The checks above have the inlet on Xmin only; these put it
    on the other edges, where it reads its boundary values through other
    indices, and make the corners it forms with x."""
    def make():
        a = io_1d()
        up = (inlet == "Ymin")
        #  1-D x -> 2-D y (Ymin inlet) or -y (Ymax inlet)
        cy = 0.5 if up else -0.5
        b = edit_toml(a, {}, dict({
            "Grid": {"Dimension": '"xyh"', "Nx": "4", "Ny": "100",
                     "Xmin": "-0.1e-3", "Xmax": "0.1e-3",
                     "Ymin": "-0.5e-3", "Ymax": "0.5e-3"},
            "Initial.SurfaceShape.EllipseBump": {"Center": arr(0.0, -0.5005 if up else 0.5005),
                                                 "Width": arr(1e30, 1.0)},
            "Substrate.EllipseBump": {"Center": arr(0.0, cy), "Width": arr(1e30, 1.0)},
            "Physics.HorizontalGravity": {"Acceleration": arr(0.0, 2.0e3 if up else -2.0e3)},
        }, **bcs(xedge, xedge, "FreshInlet" if up else "Outlet",
                 "Outlet" if up else "FreshInlet")))
        col = 1 if xedge == "Periodic" else 2
        if up:
            pmap = lambda i, j, k, d: (col, i, k)
            cols = same(["HS", "HH", "SS", "CC"])
            cols["Ux"] = ("Uy", 1)
        else:
            pmap = lambda i, j, k, d: (col, d[0] - 1 - i, k)
            cols = same(["HS", "HH", "SS", "CC"])
            cols["Ux"] = ("Uy", -1)
        return (a, b, pmap, None, cols, False)
    return make


def wm1d_io_mirror():
    """The 1-D inlet/outlet run against its mirror image, the inlet at Xmax.
    The domain is not symmetric about 0, so the image is not the same
    problem again. Xmax held its inlet at the mean of its two ghost points
    rather than of the two points either side of the face."""
    a = edit_toml(io_1d(), {}, {"Grid": {"Xmin": "-0.5e-3", "Xmax": "0.4e-3", "Nx": "90"}})
    b = edit_toml(a, {}, dict({
        "Grid": {"Xmin": "-0.4e-3", "Xmax": "0.5e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(0.5005, 0.0)},
        "Substrate.EllipseBump": {"Center": arr(-0.5, 0.0)},
        "Physics.HorizontalGravity": {"Acceleration": arr(-2.0e3, 0.0)},
    }, **bcs("Outlet", "FreshInlet", "Symmetry", "Symmetry")))
    return (a, b, lambda i, j, k, d: (d[0] - 1 - i, j, k), lambda x, y: (-x, None),
            mirrored_x(["HS", "HH", "SS", "CC", "Ux"]), False)


def periodic_film(src, twoD):
    """A film and substrate made of waves that fit the domain, periodic edges."""
    for sec in ("Initial.SurfaceShape.EllipseBump", "Substrate.EllipseBump"):
        src = drop_section(src, sec)
    L = 1.0e-3
    src = add_section(src, "Initial.SurfaceShape.PeriodicWave", {
        "Order": "1", "CreateType": '"stack"', "Base": "0.25", "Amplitude": "0.15",
        "Phase": arr(0.13e-3, 0.07e-3), "Wavelength": arr(L, L)}, "Initial.SurfaceShape")
    src = add_section(src, "Initial.SurfaceShape.PeriodicWave", {
        "Order": "2", "CreateType": '"stack"', "Base": "0.0", "Amplitude": "0.05",
        "Phase": arr(-0.21e-3, 0.3e-3), "Wavelength": arr(L / 2, L / 2)}, "Initial.SurfaceShape")
    src = add_section(src, "Substrate.PeriodicWave", {
        "CreateType": '"add"', "Base": "0.2", "Amplitude": "0.1",
        "Phase": arr(0.31e-3, -0.11e-3), "Wavelength": arr(L, L)}, "Substrate")
    return edit_toml(src, {}, bcs("Periodic", "Periodic",
                                  "Periodic" if twoD else "Symmetry",
                                  "Periodic" if twoD else "Symmetry"))


def shifted(src, sx, sy):
    """The same film, moved by (sx, sy) meters."""
    def mv(m):
        x, y = (float(v) for v in m.group(1).split(","))
        return "Phase = " + arr(x + sx, y + sy)
    return re.sub(r"^Phase = \[([^\]]*)\]$", mv, src, flags=re.M)


def wm1d_shift():
    s = edit_toml(case("periodic_1d_reg"), {}, {"Grid": {"Nx": "80"},
                  "Physics.Centrifugal": {"Enable": "false"}})
    a = fixed_steps(periodic_film(s, False), 2e-4, 50, 5e-3)
    k, dx = 13, 1.0e-3 / 80
    b = shifted(a, k * dx, 0.0)
    return (a, b, lambda i, j, kk, d: ((i - 1 + k) % (d[0] - 2) + 1, j, kk),
            None, same(["HS", "HH", "SS", "CC", "Ux"]), True)


def wm2d_shift():
    s = edit_toml(case("allphys_2d_reg"), {}, {
        "Grid": {"Nx": "32", "Ny": "20", "Xmin": "-0.5e-3", "Xmax": "0.5e-3",
                 "Ymin": "-0.5e-3", "Ymax": "0.5e-3"},
        "Physics.Centrifugal": {"Enable": "false"}})
    a = fixed_steps(periodic_film(s, True), 2e-4, 30, 3e-3)
    kx, ky = 5, 7
    b = shifted(a, kx * 1.0e-3 / 32, ky * 1.0e-3 / 20)
    return (a, b, lambda i, j, k, d: ((i - 1 + kx) % (d[0] - 2) + 1,
                                      (j - 1 + ky) % (d[1] - 2) + 1, k),
            None, same(WM_COLS), True)


def cz2d_base():
    s = case("concz_2d")
    return edit_toml(s, {}, {
        "Grid": {"Nx": "30", "Ny": "20", "Nz": "8",
                 "Xmin": "-0.5e-3", "Xmax": "0.5e-3",
                 "Ymin": "-0.35e-3", "Ymax": "0.45e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(0.08e-3, 0.03e-3),
                                             "Width": arr(0.6e-3, 0.45e-3)},
        "Substrate.EllipseBump": {"Center": arr(-0.05e-3, 0.06e-3),
                                  "Width": arr(0.5e-3, 0.35e-3)},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(3.0, 1.2)},
    })


def cz2d_transpose():
    a = fixed_steps(cz2d_base(), 1e-5, 30, 1e-4)
    b = edit_toml(a, {}, {
        "Grid": {"Nx": "20", "Ny": "30", "Xmin": "-0.35e-3", "Xmax": "0.45e-3",
                 "Ymin": "-0.5e-3", "Ymax": "0.5e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(0.03e-3, 0.08e-3),
                                             "Width": arr(0.45e-3, 0.6e-3)},
        "Substrate.EllipseBump": {"Center": arr(0.06e-3, -0.05e-3),
                                  "Width": arr(0.35e-3, 0.5e-3)},
        "Physics.HorizontalGravity": {"Acceleration": arr(1.2, 3.0)},
    })
    return (a, b, lambda i, j, k, d: (j, i, k), lambda x, y: (y, x),
            transposed(CZ_COLS), False)


def cz2d_dry_src():
    """cz2d_base, thin and drying: 0.5 um over 1 mm, so that the diffusion
    through the film is stiff at the steps of cz2d_dry, with a thinner
    precursor (the evaporation rate is cut off below 20 times its
    thickness) and a diffusivity at which the surface concentration rises
    by a few percent, unevenly with the film's shape, so that the
    Marangoni stress it drives moves the film."""
    return edit_toml(cz2d_base(), {}, {
        "Initial.SurfaceShape": {"TotalHeight": "0.5e-6", "PrecursorRatio": "0.1"},
        "Substrate": {"TotalHeight": "0.5e-6"},
        "Physics.Evaporation": {"Coef": "3e-6"},
        "Fluid.Diffusivity": {"BaseValue": "1e-8"}})


def cz2d_dry():
    """cz2d_dry_src in 30 steps of 1e-2 s. Transposed and mirrored, it
    agrees to rounding, 4e-11 of the change; a non-conservative form of the
    concentration equation differs by 8e-8 and 6e-8."""
    return fixed_steps(cz2d_dry_src(), 1e-2, 30, 0.1)


def cz2d_dry_transpose():
    a = cz2d_dry()
    b = edit_toml(a, {}, {
        "Grid": {"Nx": "20", "Ny": "30", "Xmin": "-0.35e-3", "Xmax": "0.45e-3",
                 "Ymin": "-0.5e-3", "Ymax": "0.5e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(0.03e-3, 0.08e-3),
                                             "Width": arr(0.45e-3, 0.6e-3)},
        "Substrate.EllipseBump": {"Center": arr(0.06e-3, -0.05e-3),
                                  "Width": arr(0.35e-3, 0.5e-3)},
        "Physics.HorizontalGravity": {"Acceleration": arr(1.2, 3.0)},
    })
    return (a, b, lambda i, j, k, d: (j, i, k), lambda x, y: (y, x),
            transposed(CZ_COLS), False)


def cz2d_dry_mirror_x():
    a = cz2d_dry()
    b = edit_toml(a, {}, {
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(-0.08e-3, 0.03e-3)},
        "Substrate.EllipseBump": {"Center": arr(0.05e-3, 0.06e-3)},
        "Physics.HorizontalGravity": {"Acceleration": arr(-3.0, 1.2)},
    })
    return (a, b, lambda i, j, k, d: (d[0] - 1 - i, j, k), lambda x, y: (-x, y),
            mirrored_x(CZ_COLS), False)


def cz_1d_vs_2d():
    s = edit_toml(case("concz_1d_reg"), {}, {
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(3.0, 0.0)}})
    a = fixed_steps(s, 1e-4, 40, 1e-3)
    b = uniform_in_y(a)
    cols = same(["zz", "CC", "Ux", "Uz", "HS", "HH", "SS"])
    return (a, b, lambda i, j, k, d: (i, 1, k), lambda x, y: (x, None),
            cols, False)


def rescaled_height(src, f):
    """The same film written with TotalHeight f times larger: the scale H0 the
    model is made dimensionless with changes, the physical film does not."""
    def val(key, sec):
        m = re.search(r"^\[%s\]\n(?:[^\[].*\n)*?%s = (\S+)" % (re.escape(sec), key), src, flags=re.M)
        return float(m.group(1))
    T = val("TotalHeight", "Initial.SurfaceShape")
    hs = val("PrecursorRatio", "Initial.SurfaceShape")
    b = val("Height", "Initial.SurfaceShape.EllipseBump")
    #  h = T (shape (1 - h*) + h*) with shape = Height inside the bump and 0
    #  outside; keep h and h* T fixed.
    hs2 = hs / f
    b2 = b * (1 - hs) / (f * (1 - hs2))
    return edit_toml(src, {}, {
        "Initial.SurfaceShape": {"TotalHeight": repr(T * f), "PrecursorRatio": repr(hs2)},
        "Initial.SurfaceShape.EllipseBump": {"Height": repr(b2)}})


def custom_length(src, L):
    """The same problem made dimensionless with a length scale of L meters
    instead of the capillary length."""
    assert "LengthScaleType" not in src
    return re.sub(r"^(JobName = .*)$", r'\1\nLengthScaleType = "custom"\nLengthScale = %r' % L,
                  src, count=1, flags=re.M)


def wm1d_scales():
    s = edit_toml(case("allphys_1d_reg"), {}, {"Grid": {"Nx": "100"}})
    a = fixed_steps(s, 2e-4, 50, 5e-3)
    b = custom_length(rescaled_height(a, 2.5), 0.7e-3)
    return (a, b, lambda i, j, k, d: (i, j, k), lambda x, y: (x, None),
            same(["HS", "HH", "SS", "CC", "Ux"]), False)


def cz1d_hscale():
    a = fixed_steps(case("concz_1d_reg"), 1e-4, 40, 1e-3)
    b = rescaled_height(a, 2.5)
    return (a, b, lambda i, j, k, d: (i, j, k), lambda x, y: (x, None),
            same(["zz", "CC", "Ux", "Uz", "HS", "HH", "SS"]), False)


def cz1d_lscale():
    a = fixed_steps(case("concz_1d_reg"), 1e-4, 40, 1e-3)
    b = custom_length(a, 0.7e-3)
    return (a, b, lambda i, j, k, d: (i, j, k), lambda x, y: (x, None),
            same(["zz", "CC", "Ux", "Uz", "HS", "HH", "SS"]), False)


def cz1d_mirror():
    s = edit_toml(case("concz_1d_reg"), {}, {
        "Grid": {"Xmin": "-0.5e-3", "Xmax": "0.4e-3"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(3.0, 0.0)}})
    a = fixed_steps(s, 1e-4, 40, 1e-3)
    b = edit_toml(a, {}, {
        "Grid": {"Xmin": "-0.4e-3", "Xmax": "0.5e-3"},
        "Initial.SurfaceShape.EllipseBump": {"Center": arr(-0.5, 0.0)},
        "Substrate.EllipseBump": {"Center": arr(-0.5, 0.0)},
        "Physics.HorizontalGravity": {"Acceleration": arr(-3.0, 0.0)}})
    return (a, b, lambda i, j, k, d: (d[0] - 1 - i, j, k), lambda x, y: (-x, None),
            mirrored_x(["zz", "CC", "Ux", "Uz", "HS", "HH", "SS"]), False)


def wm1d_nozzle_mirror():
    """A nozzle off the center of a flat film, against the same nozzle on
    the other side. The footprint fraction of each cell was measured at the
    left end of each sub-division of the cell, which fed the film half a
    sub-division to the right of the nozzle."""
    s = edit_toml(case("nozzle_1d_offset"), {}, {
        "Grid": {"Nx": "100"},
        "Initial.SurfaceShape.EllipseBump": {"Height": "0.0"},
        "Substrate.EllipseBump": {"Height": "0.0"},
        "Physics.HorizontalGravity": {"Enable": "false"},
        "Physics.Centrifugal": {"Enable": "false"},
        "Physics.Source": {"SubRegions": "4", "TimeConstant": "0.0"}})
    a = fixed_steps(s, 2e-4, 40, 4e-3)
    b = edit_toml(a, {}, {"Physics.Source": {"ArmCenter": arr(-395.25e-3, 0.0)}})
    return (a, b, lambda i, j, k, d: (d[0] - 1 - i, j, k), lambda x, y: (-x, None),
            mirrored_x(["HS", "HH", "SS", "CC", "Ux"]), False)


def cz1d_shift():
    s = edit_toml(case("concz_1d_reg"), {}, {"Grid": {"Nx": "60"}})
    a = fixed_steps(periodic_film(s, False), 1e-4, 40, 1e-3)
    k = 11
    b = shifted(a, k * 1.0e-3 / 60, 0.0)
    return (a, b, lambda i, j, kk, d: ((i - 1 + k) % (d[0] - 2) + 1, j, kk),
            None, same(["zz", "CC", "Ux", "Uz", "HS", "HH", "SS"]), True)


# ----------------------------------------------------------------------
#  Two checks against a known answer. They sit here because they share the
#  machinery with the ones above: a body-force velocity must be the one the
#  film-height equation moves the film with, wherever it is read, and no
#  solute may be created at the free surface.
# ----------------------------------------------------------------------
def flat_film(src, extra):
    """A flat film of 5.5 um on a flat substrate, precursor 0.1."""
    return edit_toml(src, {}, dict({
        "Initial.SurfaceShape": {"PrecursorRatio": "0.1"},
        "Initial.SurfaceShape.EllipseBump": {"Base": "0.5", "Height": "0.0"},
        "Substrate.EllipseBump": {"Height": "0.0"},
        "Physics.SolutoCapillarity": {"Enable": "false"},
        "Physics.VerticalGravity": {"Enable": "false"},
        "Physics.Evaporation": {"Enable": "false"}}, **extra))


def body_velocity(exe, workdir, timeout):
    """On a flat film the in-plane gravity and centrifugal velocities are
    closed-form: u = rho g h^2 / (3 mu) and rho w^2 x h^2 / (3 mu), mean over
    the depth, pointing along Acceleration and away from the axis."""
    rho, mu, g, w = 1000.0, 0.1, 4.903, 41.89
    worst = 0.0
    for dim in ("xh", "xyh"):
        for term in ("HorizontalGravity", "Centrifugal"):
            other = "Centrifugal" if term == "HorizontalGravity" else "HorizontalGravity"
            s = flat_film(case("allphys_1d_reg"), {
                "Grid": {"Nx": "40"},
                "Physics.%s" % other: {"Enable": "false"},
                "Physics.%s" % term: ({"Enable": "true", "Acceleration": arr(g, 0.6 * g)}
                                      if term == "HorizontalGravity" else
                                      {"Enable": "true", "RotationRate": repr(w),
                                       "RotationAxis": arr(0.0, 0.0)})})
            if dim == "xyh":
                s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": "30",
                                               "Ymin": "-0.3e-3", "Ymax": "0.3e-3"}})
            s = fixed_steps(s, 1e-4, 1, 1e-4)
            d = os.path.join(workdir, "body_%s_%s" % (dim, term))
            rc, log = run(exe, d, "a", s, timeout)
            if rc != 0:
                return None, "run failed: " + log[-400:]
            f, nx, ny, nz = read_field(last_field(os.path.join(d, "a"))[0])
            #  The cell centers away from the walls
            for j in (range(3, ny - 3) if ny > 1 else [0]):
                for i in range(3, nx - 3):
                    r = j * nx + i
                    h = f["HH"][r] * 1e-6
                    x = f["xx"][r] * 1e-3
                    y = f["yy"][r] * 1e-3 if "yy" in f else 0.0
                    if term == "HorizontalGravity":
                        ex, ey = rho * g * h * h / (3 * mu), rho * 0.6 * g * h * h / (3 * mu)
                    else:
                        ex, ey = rho * w * w * x * h * h / (3 * mu), rho * w * w * y * h * h / (3 * mu)
                    worst = max(worst, abs(f["Ux"][r] - ex) / abs(ex) if ex else abs(f["Ux"][r]))
                    if ny > 1 and ey:
                        worst = max(worst, abs(f["Uy"][r] - ey) / abs(ey))
    return worst, None


def cz_solute(exe, workdir, timeout):
    """A drying flat ConcZ film: the solvent leaves, the solute stays. The
    depth-integrated solute h * mean(c) must hold while h falls by 0.4 %."""
    s = flat_film(case("concz_1d_reg"), {"Grid": {"Nx": "20"},
                  "Physics.HorizontalGravity": {"Enable": "false"},
                  "Physics.Evaporation": {"Enable": "true", "Coef": "1e-6"}})
    s = fixed_steps(s, 1e-4, 1000, 0.05)
    d = os.path.join(workdir, "cz_solute")
    rc, log = run(exe, d, "a", s, timeout)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    first, last = last_field(os.path.join(d, "a"))
    out = []
    for p in (first, last):
        f, nx, ny, nz = read_field(p)
        i = nx // 2
        h = f["HH"][i]
        cm = sum(f["CC"][k * nx + i] for k in range(1, nz - 1)) / (nz - 2)
        out.append((h, h * cm))
    (h0, m0), (h1, m1) = out
    if (h0 - h1) / h0 < 1e-3:
        return None, "the film did not dry (h %g -> %g)" % (h0, h1)
    return abs(m1 - m0) / (h0 - h1) / (m0 / h0), None


def cz_columns(path):
    """The solute of a ConcZ field file, the sum over the columns of h times
    the mean of c over the cells through the film (the points on the edges,
    the substrate and the surface, which the file carries too, left out),
    and the concentrations of every point."""
    f, nx, ny, nz = read_field(path)
    js = range(1, ny - 1) if ny > 1 else [0]
    tot = 0.0
    for j in js:
        for i in range(1, nx - 1):
            col = [f["CC"][(k * ny + j) * nx + i] for k in range(1, nz - 1)]
            tot += f["HH"][j * nx + i] * sum(col) / len(col)
    return tot, f["CC"]


def cz_slope_1d():
    """concz_1d_reg on a sloping substrate, flowing (capillary, in-plane
    gravity, Marangoni), with a diffusivity high enough that the solute
    spreads along the film over the run: the conormal part of the diffusive
    flux through the sloping surface and substrate, which leaks when a
    surface condition leaves it out, is then as large as the vertical
    one."""
    s = edit_toml(case("concz_1d_reg"), {}, {
        "Grid": {"Nx": "40", "Nz": "8"},
        "Initial.SurfaceShape": {"PrecursorRatio": "0.01"},
        "Initial.SurfaceShape.EllipseBump": {"Base": "0.5", "Height": "0.3"},
        "Fluid.Diffusivity": {"BaseValue": "1.0e-7"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(3.0, 0.0)}})
    return s


def cz_solute_slope(exe, workdir, timeout):
    """The solute of a drying ConcZ film that is not flat, in 1-D and 2-D:
    |S_end / S_start - 1|. Periodic edges: the in-plane gravity passes
    through a Symmetry wall with the liquid it carries (an open finding,
    ch-notverified), and the solute would go with it."""
    worst = 0.0
    per = bcs("Periodic", "Periodic", "Periodic", "Periodic")
    for tag, s in (("1d", fixed_steps(edit_toml(cz_slope_1d(), {}, dict(per, **{
                                          "Physics.Evaporation": {"Coef": "5e-5"}})), 2e-4, 40, 4e-3)),
                   ("2d", edit_toml(cz2d_dry(), {}, per))):
        d = os.path.join(workdir, "cz_solute_slope")
        rc, log = run(exe, d, tag, s, timeout)
        if rc != 0:
            return None, "run %s failed: %s" % (tag, log[-400:])
        first, last = last_field(os.path.join(d, tag))
        (s0, _), (s1, _) = cz_columns(first), cz_columns(last)
        v0, v1 = (sum(read_field(p)[0]["HH"]) for p in (first, last))
        if (v0 - v1) / v0 < 1e-3:
            return None, "the %s film did not dry" % tag
        worst = max(worst, abs(s1 / s0 - 1))
    return worst, None


def cz_uniform(exe, workdir, timeout):
    """A uniform concentration in a ConcZ film that is not flat, flows and
    levels, without evaporation (1-D and 2-D): max |c - c0| over every
    point of every field file. The 2-D one is cz2d_dry_src without the
    evaporation, whose steps make the diffusion through the film stiff:
    rows in non-conservative form sum to zero only to the rounding of terms
    that large, and the concentration drifts (by 1e-6 in tools/proptest.py
    seed 2104)."""
    worst = 0.0
    wet = {"Physics.Evaporation": {"Enable": "false"}}
    for tag, s in (("1d", fixed_steps(edit_toml(cz_slope_1d(), {}, wet), 2e-4, 40, 4e-3)),
                   ("2d", fixed_steps(edit_toml(cz2d_dry_src(), {}, wet), 1e-2, 30, 0.1))):
        d = os.path.join(workdir, "cz_uniform")
        rc, log = run(exe, d, tag, s, timeout)
        if rc != 0:
            return None, "run %s failed: %s" % (tag, log[-400:])
        fs = sorted(glob.glob(os.path.join(d, tag, "field*.dat")))
        c0 = read_field(fs[0])[0]["CC"][0]
        moved = 0.0
        for p in fs:
            f = read_field(p)[0]
            worst = max(worst, max(abs(c - c0) for c in f["CC"]))
            moved = max(moved, max(abs(a - b) for a, b in zip(f["HH"], read_field(fs[0])[0]["HH"])))
        if moved < 1e-6 * max(read_field(fs[0])[0]["HH"]):
            return None, "the %s film did not move" % tag
    return worst, None


def wm_retry(exe, workdir, timeout):
    """A step that the time-step control rejected and retried with a smaller
    dt, against the same step taken at that dt from the start. The retry
    used to build the explicit half of its Crank-Nicolson update from the
    rejected attempt's answer instead of the start of the step, so the two
    differed by 3e-5 of the change; they must agree to rounding."""
    base = edit_toml(case("leveling_1d_reg"), {}, {
        "Numerics": {"InitialDt": "1e-3", "MaxSteps": "1", "SteadyStateCriterion": "-1"},
        "Numerics.TimeStepControl": {"Criterion": "1e-5"},
        "Numerics.HeightEquation.LinearSolver": {"Backend": '"eigen"'},
        "Output.FieldData": {"TimeInterval": "1.0", "NeutralFile": "true",
                             "VTKFile": "false"}})
    rc, log = run(exe, workdir, "a", base, timeout)
    if rc != 0:
        return None, "adaptive run: rc=%d\n%s" % (rc, log[-1200:])
    k = len(re.findall(r"dt- \(", log))
    if k == 0:
        return None, "the first step was not rejected, so no retry was tested"
    dt = 1e-3
    for _ in range(k):          # the same products the solver forms
        dt *= 0.75
    fixed = edit_toml(base, {}, {"Numerics": {"InitialDt": repr(dt)},
                                 "Numerics.TimeStepControl": {"Enable": "false"}})
    rc, log = run(exe, workdir, "b", fixed, timeout)
    if rc != 0:
        return None, "fixed-dt run: rc=%d\n%s" % (rc, log[-1200:])
    a0, a1 = last_field(os.path.join(workdir, "a"))
    b0, b1 = last_field(os.path.join(workdir, "b"))
    if os.path.basename(a1) != os.path.basename(b1) or os.path.basename(a1) == "field00000.dat":
        return None, "the two runs did not reach the same output"
    diff = compare(read_field(a1), read_field(b1), read_field(a0),
                   lambda i, j, kk, d: (i, j, kk), same(["HS"]))
    return diff["HS"][0], None


def cz2d_deficit_transpose(exe, workdir, timeout):
    """A ConcZ column with less than no solute, next to a wall, filled from
    the columns beside it (CZ_HeightEquation_2D::RedistributeDeficit): the
    2-D film of cz2d_transpose and its transpose run one step, every layer
    of the column at x index 10 on the first row in y (in the transpose,
    the same column on the first row in x) is set to c = -0.05 in the field
    file, and the two are continued from it. The wall leaves the column one
    neighbor in y and two in x, so a fault in either direction makes the
    two runs part. The result is the largest difference relative to the
    change, as for the checks of CHECKS; the run before the edit must
    report no solute moved. (Mutation testing: the neighbors in y read
    from the wrong side passed every other check.)"""
    a, b, pmap, xymap, cols, _ = cz2d_transpose()
    res = []
    for tag, src, pick in (("a", a, lambda xs, ys: (xs[10], ys[1])),
                           ("b", b, lambda xs, ys: (xs[1], ys[10]))):
        first = fixed_steps(src, 1e-5, 1, 1e-5)
        rc, log = run(exe, workdir, tag + "0", first, timeout)
        if rc != 0:
            return None, "run %s0 failed: %s" % (tag, log[-400:])
        #  Before the edit nothing fell below 0: JobInfo.txt says 0 and the
        #  log has no summary line of the moves (mutation testing: counting
        #  every step as one that moved solute changed nothing else)
        ji = open(os.path.join(workdir, tag + "0", "JobInfo.txt")).read()
        m = re.search(r"Solute moved to c < 0 \[-\] :\s*(.*)", ji)
        if not m or m.group(1).strip() != "0" or "Solute moved to c < 0 from" in log:
            return None, "a run without c < 0 reports solute moved: %r" % (m.group(1) if m else None)
        lines = open(os.path.join(workdir, tag + "0", "field00001.dat")).read().split("\n")
        hdr = [h.split("[")[0].strip() for h in lines[0].lstrip("#").split("\t")]
        ix, iy, ic = hdr.index("xx"), hdr.index("yy"), hdr.index("CC")
        rows = [l.split("\t") for l in lines[1:] if l.strip()]
        x0, y0 = pick(sorted({float(r[ix]) for r in rows}), sorted({float(r[iy]) for r in rows}))
        for n, l in enumerate(lines[1:], 1):
            t = l.split("\t")
            if l.strip() and float(t[ix]) == x0 and float(t[iy]) == y0:
                t[ic] = "-0.05"
                lines[n] = "\t".join(t)
        os.makedirs(os.path.join(workdir, tag + "_start"), exist_ok=True)
        with open(os.path.join(workdir, tag + "_start", "field00001.dat"), "w") as f:
            f.write("\n".join(lines))
        cont = fixed_steps(src, 1e-5, 10, 1e-4)
        cont = re.sub(r"^\[Restart\]\n(?:[A-Za-z].*\n)*",
                      '[Restart]\nEnable = true\nASCII = true\nFileName = "%s_start/field00001.dat"\n'
                      "RestartTime = 1e-5\nFileIndex = 2\n" % tag, cont, count=1, flags=re.M)
        rc, log = run(exe, workdir, tag, cont, timeout)
        if rc != 0:
            return None, "run %s failed: %s" % (tag, log[-400:])
        fs = sorted(glob.glob(os.path.join(workdir, tag, "field*.dat")))
        res.append((read_field(os.path.join(workdir, tag + "_start", "field00001.dat")), read_field(fs[-1])))
    (a0, a1), (b0, b1) = res
    why = coords_agree(a1, b1, pmap, xymap)
    if why:
        return None, why
    diff = compare(a1, b1, a0, pmap, cols)
    return max(v for v, _ in diff.values()), None


def cz_wm_advection(exe, workdir, timeout):
    """The solute is carried by the flow the same way in both models when
    it is uniform through the thickness. A strong in-plane gravity is
    switched on over a drying flat film whose evaporation varies with x
    (SpatialFunc "negPow"), and the change it makes to the depth-mean
    concentration is compared between ConcZ and well-mixed; a vertical
    diffusivity of 1e-7 keeps the ConcZ profile flat. The viscosity does
    not depend on c and Marangoni is off, so the film is the same in both.
    The result is the largest difference relative to the largest change."""
    base = flat_film(case("concz_1d_reg"), {})
    common = {"Grid": {"Nx": "50"}, "Fluid.Viscosity": {"Exponent": "0.0"},
              "Physics.Evaporation": {"Enable": "true", "Coef": "1e-5",
                                      "SpatialFunc": '"negPow"', "NegPowRmax": "0.6e-3"},
              "Fluid.Diffusivity": {"BaseValue": "1e-7", "Exponent": "0.0"}}
    c = {}
    for model in ("CZ", "WM"):
        for g in ("0.0", "4903.0"):
            s = edit_toml(base, {}, dict(common, **{
                "Physics.ConcField": {"Model": '"%s"' % model},
                "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": "[%s, 0.0]" % g}}))
            s = fixed_steps(s, 1e-5, 500, 0.005)
            d = os.path.join(workdir, "%s_g%s" % (model, g))
            rc, log = run(exe, d, "a", s, timeout)
            if rc != 0:
                return None, "run failed: " + log[-400:]
            f, nx, ny, nz = read_field(last_field(os.path.join(d, "a"))[1])
            c[(model, g)] = ([sum(f["CC"][k * nx + i] for k in range(1, nz - 1)) / (nz - 2)
                              for i in range(nx)] if model == "CZ" else f["CC"])
    dcz = [a - b for a, b in zip(c[("CZ", "4903.0")], c[("CZ", "0.0")])]
    dwm = [a - b for a, b in zip(c[("WM", "4903.0")], c[("WM", "0.0")])]
    big = max(abs(v) for v in dwm)
    if big == 0:
        return None, "the flow changed nothing"
    #  The four cells next to each wall are left out: there the liquid the
    #  in-plane gravity pushes through a Symmetry edge brings in the
    #  mirrored ghost, and the two models' different stencils at the wall
    #  (advective with the ghost for ConcZ, the well-mixed form for WM) give
    #  answers that differ by tens of percent without either being the
    #  right one.
    return max(abs(a - b) for a, b in zip(dcz[5:-5], dwm[5:-5])) / big, None


CHECKS = [("wm2d_transpose", wm2d_transpose),
          ("wm2d_mirror_x", wm2d_mirror_x),
          ("wm2d_nozzle_transpose", wm2d_nozzle_transpose),
          ("wm2d_nozzle_thin", wm2d_nozzle_thin),
          ("wm2d_nozzle_mirror_x", wm2d_nozzle_mirror_x),
          ("wm2d_spin_transpose", spin_utopia(wm2d_transpose)),
          ("wm2d_spin_mirror_x", spin_utopia(wm2d_mirror_x)),
          ("wm_1d_vs_2d", wm_1d_vs_2d),
          ("wm_1d_vs_2d_koren", wm_1d_vs_2d_koren),
          ("wm_io_1d_vs_2d", io_1d_vs_2d("Periodic")),
          ("wm_io_1d_vs_2d_sym", io_1d_vs_2d("Symmetry")),
          ("wm_io_1d_vs_2d_ymin", io_1d_vs_2d_y("Ymin", "Periodic")),
          ("wm_io_1d_vs_2d_ymax_sym", io_1d_vs_2d_y("Ymax", "Symmetry")),
          ("wm1d_io_mirror", wm1d_io_mirror),
          ("wm1d_shift", wm1d_shift),
          ("wm2d_shift", wm2d_shift),
          ("wm1d_scales", wm1d_scales),
          ("cz1d_hscale", cz1d_hscale),
          ("cz1d_lscale", cz1d_lscale),
          ("wm1d_nozzle_mirror", wm1d_nozzle_mirror),
          ("cz1d_mirror", cz1d_mirror),
          ("cz1d_shift", cz1d_shift),
          ("cz2d_transpose", cz2d_transpose),
          ("cz2d_dry_transpose", cz2d_dry_transpose),
          ("cz2d_dry_mirror_x", cz2d_dry_mirror_x),
          ("cz_1d_vs_2d", cz_1d_vs_2d)]


# ----------------------------------------------------------------------
#  [Numerics.Convection] for every run of this script (--scheme / --conc /
#  --height); left empty, the inputs keep the models' defaults. A case that
#  sets the block itself keeps its own.
CONVECTION = {}

#  Checks held to less than --tol (relative to the change over the run)
TIGHT = {"cz2d_dry_transpose": 1e-9, "cz2d_dry_mirror_x": 1e-9}


def with_convection(src):
    if not CONVECTION or re.search(r"^\[Numerics\.Convection\]", src, flags=re.M):
        return src
    cz = re.search(r'^Model = "CZ"', src, flags=re.M) is not None
    kv = ["%s = \"%s\"" % (k, v) for k, v in sorted(CONVECTION.items())
          if v and not (cz and k == "Height")]
    return src + ("\n[Numerics.Convection]\n%s\n" % "\n".join(kv) if kv else "")


def run(exe, d, name, src, timeout):
    os.makedirs(d, exist_ok=True)
    src = with_convection(src)
    src = re.sub(r"^JobName = .*$", 'JobName = "%s"' % name, src, count=1, flags=re.M)
    with open(os.path.join(d, name + ".toml"), "w") as f:
        f.write(src)
    try:
        r = subprocess.run([exe, name + ".toml"], cwd=d, capture_output=True,
                           text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return -1, "timed out after %g s" % timeout
    return r.returncode, r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", action="append", help="run only this check (repeatable)")
    ap.add_argument("--tol", type=float, default=1e-7,
                    help="largest difference allowed, relative to the change over the run")
    ap.add_argument("--run-timeout", type=float, default=600)
    ap.add_argument("--scheme", help="Numerics.Convection: this scheme for the concentration"
                    " and (well-mixed) the film height in every run")
    ap.add_argument("--conc", help="Numerics.Convection.Concentration alone")
    ap.add_argument("--height", help="Numerics.Convection.Height alone (well-mixed runs)")
    args = ap.parse_args()
    if args.scheme:
        CONVECTION["Concentration"] = CONVECTION["Height"] = args.scheme
    if args.conc:
        CONVECTION["Concentration"] = args.conc
    if args.height:
        CONVECTION["Height"] = args.height
    exe = os.path.abspath(args.exe)
    shutil.rmtree(args.workdir, ignore_errors=True)

    bad = 0
    for name, make in CHECKS:
        if args.only and name not in args.only:
            continue
        a, b, pmap, xymap, cols, centers_only = make()
        res = []
        for tag, src in (("a", a), ("b", b)):
            d = os.path.join(args.workdir, name)
            rc, log = run(exe, d, tag, src, args.run_timeout)
            if rc != 0:
                print("FAIL %-15s run %s: rc=%d\n%s" % (name, tag, rc, log[-1200:]))
                res = None
                break
            first, last = last_field(os.path.join(d, tag))
            res.append((read_field(first), read_field(last), last))
        if res is None:
            bad += 1
            continue
        (a0, a1, fa), (b0, b1, fb) = res
        if os.path.basename(fa) != os.path.basename(fb) or os.path.basename(fa) == "field00000.dat":
            print("FAIL %-15s the two runs did not reach the same output (%s, %s)"
                  % (name, os.path.basename(fa), os.path.basename(fb)))
            bad += 1
            continue
        if xymap is not None:
            why = coords_agree(a1, b1, pmap, xymap)
            if why:
                print("FAIL %-15s %s" % (name, why))
                bad += 1
                continue
        diff = compare(a1, b1, a0, pmap, cols, centers_only)
        worst = max(v for v, _ in diff.values())
        ok = worst <= min(args.tol, TIGHT.get(name, args.tol))
        bad += not ok
        print("%s %-15s worst %.1e   %s" % ("ok  " if ok else "FAIL", name, worst,
              "  ".join("%s %.0e" % (n, v) for n, (v, _) in diff.items())))
    #  The known answers. body_velocity: relative error of the velocity.
    #  cz_solute: the solute gained or lost, as a fraction of what the
    #  solvent lost -- 1 would mean the solute went with it.
    for name, fn, tol in (("body_velocity", body_velocity, 1e-10),
                          ("cz_solute", cz_solute, 1e-10),
                          ("cz_solute_slope", cz_solute_slope, 1e-12),
                          ("cz_uniform", cz_uniform, 1e-13),
                          ("cz_wm_advection", cz_wm_advection, 3e-2),
                          ("cz2d_deficit_transpose", cz2d_deficit_transpose, 1e-7),
                          ("wm_retry", wm_retry, 1e-8)):
        if args.only and name not in args.only:
            continue
        err, why = fn(exe, os.path.join(args.workdir, name), args.run_timeout)
        if why:
            print("FAIL %-15s %s" % (name, why))
            bad += 1
            continue
        ok = err <= tol
        bad += not ok
        print("%s %-15s error %.1e (allowed %.0e)" % ("ok  " if ok else "FAIL", name, err, tol))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
