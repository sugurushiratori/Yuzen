#!/usr/bin/env python3
"""Property-based testing of the solver on random valid inputs.

    python3 tools/proptest.py --exe build/bin/yuzen --workdir /tmp/pt
    python3 tools/proptest.py --exe build/bin/yuzen --workdir /tmp/pt --count 300 --seed 1000
    python3 tools/proptest.py --exe build/bin/yuzen --workdir /tmp/pt --only 1234 --keep
    python3 tools/proptest.py --exe build/bin/yuzen --workdir /tmp/pt --quick   # the ctest slice

Each seed gives one input, drawn at random but valid and physically sensible:
the model (well-mixed or ConcZ), 1-D or 2-D on a small grid (16 to 48
cells; ConcZ with 4 to 10 layers), a random subset of the physics terms
(capillary pressure, Marangoni, the two gravities, the centrifugal term,
evaporation, the disjoining pressure, a nozzle) with parameters spread over
orders of magnitude around those of the shipped cases, a boundary condition
per edge (Symmetry, a Periodic pair, Pinning with the edge film it needs,
FreshInlet, Outlet; ConcZ only the first two), one to three initial shapes
and a substrate, convection schemes, a fixed or an adaptive time step, for a
few dozen steps. The same seed gives the same input on every machine (only
the standard library's random.Random is used).

Every run is held to these properties:

  ends      the run ends with exit code 0, or 1 with a stop the solver
            means (a divergence it detects, a time step below MinDt); never
            a signal, a sanitizer report, a hang, or an exit 1 without a
            reason
  finite    no NaN or inf in field*.dat, History.dat or CalcInfo.dat
  film      the film thickness of every written field stays positive in
            the interior (above half the precursor film with the disjoining
            pressure, whose minimum it is)
  bounds    0 <= c <= 1 with a bounded scheme (Koren, VanLeer, Upwind1) for
            the concentration while the Courant number stays at most 1
  files     every row of History.dat and CalcInfo.dat has the columns
            of its header
  clock     the time the run reached is the sum of the steps the equations
            were solved with (JobInfo.txt), to 1e-12 of the time run: the
            time-step control lengthened the step before the clock had
            advanced, so the clock added a step that had not been solved
            and ran ahead (0.37 % at the end of leveling_1d_reg)
  volume    |VolErr| <= 1e-9 without evaporation, whatever the edges
            (VolErr counts what a nozzle supplies, in 1-D what the
            centrifugal term flings sideways, and what crosses a
            FreshInlet or an Outlet as the step's face fluxes put it
            through)
  solute    |ResinVolErr| <= 1e-11 for the well-mixed model, whatever the
            edges, whose balance of h c keeps it to rounding (it counts
            the solute a nozzle supplies, SuppliedResin, and the solute
            flux through the open edges; with the advection
            on, since without it the model does not carry the solute the
            flow moves; 1e-2 when the
            concentration equation was solved for c, which kept it only
            to the accuracy of the scheme: found, corrected); for ConcZ,
            the solute sum h mean(c) of the first and the last field file
            to 1e-13 with closed edges (the only ones it takes),
            evaporation or not (with the 1-D centrifugal term,
            which takes solute out of the field files, ResinVolErr); both
            also for a film dried out to c = 1
  uniform   without evaporation: the uniform initial concentration stays
            uniform, max |c - c0| over the last field written <= 1e-12 for
            ConcZ and 1e-10 for the well-mixed model (the discrete
            geometric conservation law; the well-mixed solute is formed
            from the face fluxes once more, so a uniform c moves by their
            rounding relative to the film, eps dt sum|F| / h: up to 3e-12
            over seeds 1000-1199, on a thin film next to a Pinning edge)

"Closed" is Symmetry, Pinning or Periodic on every edge (a Periodic edge
is not drawn with the centrifugal term, which the solver refuses: a rotating
film is not periodic). The 1-D centrifugal term flings the film into the y
that is not solved for: found uncounted in VolErr and ResinVolErr, which
count it as outflowed since. A wall carries no flux of any term: the in-plane body
forces and the Marangoni stress took liquid through a Symmetry or a Pinning
wall (VolErr -5.9e-5 over 0.2 s), and the 1-D volume took the ghosts into
the edge cells, so that a pinned film seemed to lose ~5e-5 of its volume;
ConcZ let the body forces through a Symmetry wall as well (seed 1015);
all found by this tool, and asserted since.

and, for a subset of the seeds:

  mirror / transpose  a 2-D configuration against itself mirrored in x, or
            with x and y exchanged (everything that has a direction moves
            with it), fixed dt: the fields agree to 1e-7 of their change
            (the velocities and a film that hardly moves are measured as
            compare() says)
  threads   one thread against two (eigen backend): every output file the
            same to the bit
  restart   the run stopped twice and restarted from its field files and
            state files (tools/restart_roundtrip.py --chain): the run it
            would have been straight through, to the bit (a pinned edge with
            evaporation was left out, seeds 1035 and 1206: corrected, see
            subsets())

A violation prints the seed; --only SEED --keep reproduces it and leaves the
input and outputs in the work directory. --quick runs a fixed slice of 20
small seeds with every property, and the seeds that found a defect
(REGRESSION_SEEDS), in well under a minute (the ctest entry `proptest`).

Only the Python standard library is used.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

SCHEMES = ["UTOPIA", "KK", "Koren", "VanLeer", "Central4", "QUICK", "Upwind1", "Upwind2", "Central2"]
BOUNDED = {"Koren", "VanLeer", "Upwind1"}
#  The edges that let nothing through. (A pinned film is closed too: its
#  volume is the sum of the cells, not a mean of face values that would take
#  the ghosts in, -h_first dx/2 at a pinned edge; FieldInfo.cpp.)
CLOSED = {"Symmetry", "Periodic", "Pinning"}
G = 9.80665

#  The stops a run may end with that are the solver's own (exit code 1)
INTENDED = ["Diverged:", "Under flow: the requested step fell below",
            "cannot be retried at a shorter one"]

VOL_TOL = 1e-9
SOLUTE_TOL = 1e-11
#  (1e-11 until seeds 1054 and 2105: the shift that puts the linear
#  solver's residual back into the ConcZ solute; 2.4e-14 over 500 seeds since)
CZ_SOLUTE_TOL = 1e-13
CZ_UNIFORM_TOL = 1e-12
WM_UNIFORM_TOL = 1e-10
EQUIV_TOL = 1e-7


# ---------------------------------------------------------------------------
#  The input
# ---------------------------------------------------------------------------
def num(x):
    """A float as TOML, exactly (repr round-trips)."""
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, str):
        return '"%s"' % x
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(num(v) for v in x) + "]"
    r = repr(float(x))
    return r if ("e" in r or "." in r or "inf" in r or "nan" in r) else r + ".0"


class Input:
    """Tables in order; an array of tables is a list of dicts."""

    def __init__(self):
        self.t = {"": {}}

    def set(self, table, **kv):
        self.t.setdefault(table, {})
        self.t[table].update(kv)

    def add(self, table, **kv):
        self.t.setdefault(table, [])
        self.t[table].append(dict(kv))

    def text(self):
        out = []
        for name, body in self.t.items():
            if isinstance(body, list):
                for b in body:
                    out.append("\n[[%s]]" % name)
                    out += ["%s = %s" % (k, num(v)) for k, v in b.items()]
            else:
                if name:
                    out.append("\n[%s]" % name)
                out += ["%s = %s" % (k, num(v)) for k, v in body.items()]
        return "\n".join(out) + "\n"


def logu(rng, a, b):
    """log-uniform in [10^a, 10^b]"""
    return 10.0 ** rng.uniform(a, b)


def generate(seed, quick=False):
    """One random valid input. Returns (Input, meta)."""
    rng = random.Random(seed)
    m = {"seed": seed}
    cz = rng.random() < 0.3
    twoD = rng.random() < (0.35 if cz else 0.45)
    m["model"], m["twoD"] = ("CZ" if cz else "WM"), twoD
    if quick:
        nx = rng.randint(12, 20) if twoD else rng.randint(16, 32)
        ny = rng.randint(10, 16) if twoD else 1
    elif twoD:
        nx, ny = (rng.randint(10, 20), rng.randint(10, 20)) if cz else (rng.randint(12, 32), rng.randint(12, 32))
    else:
        nx, ny = rng.randint(16, 48), 1
    nz = (rng.randint(4, 6 if (twoD or quick) else 10)) if cz else 1
    lx = logu(rng, -3.7, -2.5)                 # 0.2 .. 3 mm
    ly = lx * logu(rng, -0.3, 0.3) if twoD else 2.0
    x0 = -lx * rng.uniform(0.0, 1.0)
    y0 = -ly * rng.uniform(0.0, 1.0) if twoD else -1.0
    m.update(nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, x0=x0, y0=y0)

    I = Input()
    I.set("", SolverName="Yuzen-1.0", InputFormatVersion="1.0", JobName="run")
    I.set("Grid", Dimension="xyh" if twoD else "xh", Xmin=x0, Xmax=x0 + lx,
          Ymin=y0, Ymax=y0 + ly, Nx=nx, Ny=ny, Nz=nz)

    # --- the fluid
    rho = logu(rng, 2.9, 3.3)
    sig = logu(rng, -2.0, -1.2)
    mu = logu(rng, -3.0, 0.0)
    c0 = rng.uniform(0.05, 0.95)
    m["c0"] = c0
    I.set("Initial.Concentration", Type="uniform", Value=c0)
    H0 = logu(rng, -6.5, -4.7)
    m["H0"] = H0

    # --- edges
    if cz:
        kinds = ["Symmetry", "Periodic"]
    else:
        kinds = ["Symmetry", "Symmetry", "Periodic", "Pinning", "FreshInlet", "Outlet"]
    bc = {}
    for pair in (("Xmin", "Xmax"), ("Ymin", "Ymax")) if twoD else (("Xmin", "Xmax"),):
        if rng.random() < 0.3:
            bc[pair[0]] = bc[pair[1]] = "Periodic"
        else:
            for e in pair:
                bc[e] = rng.choice([k for k in kinds if k != "Periodic"])
    if twoD and any(bc[e] == "Pinning" for e in ("Xmin", "Xmax")):
        #  The film can be brought down to a pinned wall on one axis only:
        #  EdgeBumps on two adjacent edges add up, and each leaves its full
        #  height along the other's wall.
        for e in ("Ymin", "Ymax"):
            if bc[e] == "Pinning":
                bc[e] = "Symmetry"
    if twoD:
        #  A Pinning corner next to an inlet or an outlet is not implemented
        #  (ch-ref-boundary-output): the run aborts, by design.
        for a in ("Xmin", "Xmax"):
            for b in ("Ymin", "Ymax"):
                if {bc[a], bc[b]} & {"Pinning"} and {bc[a], bc[b]} & {"FreshInlet", "Outlet"}:
                    bc[b] = "Symmetry" if bc[b] != "Pinning" else bc[b]
                    if bc[a] in ("FreshInlet", "Outlet"):
                        bc[a] = "Symmetry"
    else:
        bc["Ymin"] = bc["Ymax"] = "Symmetry"
    #  One pinned edge at most: the EdgeBumps of two add up, and each leaves
    #  its far-field height at the other's wall.
    pins = [e for e in ("Xmin", "Xmax", "Ymin", "Ymax") if bc[e] == "Pinning"]
    for e in pins[1:]:
        bc[e] = "Symmetry"
    m["bc"] = bc
    for e in ("Xmin", "Xmax", "Ymin", "Ymax"):
        I.set("Boundary.%sBC" % e, Type=bc[e])
    pinned = [e for e in bc if bc[e] == "Pinning" and (twoD or e[0] == "X")]

    # --- the physics
    disj = (not cz) and not pinned and rng.random() < 0.2
    #  A pinned edge needs the film to reach zero at the wall, and the film is
    #  never thinner than the precursor: a thin one (ch-ref-boundary-output)
    hstar = (logu(rng, -2.5, -1.5) if pinned else logu(rng, -2.0, -0.7) if disj
             else rng.uniform(0.1, 0.6))
    m["hstar"], m["disj"] = hstar, disj
    I.set("Initial.SurfaceShape", TotalHeight=H0, PrecursorRatio=hstar)

    # the initial shapes
    def cx():
        return x0 + lx * rng.uniform(0.2, 0.8)

    def cy():
        return y0 + ly * rng.uniform(0.2, 0.8) if twoD else 0.0

    shapes = []
    for e in pinned:
        #  The film falls to the wall over OuterLength; the check before the first
        #  step extrapolates it from the first two cells, so OuterLength spans
        #  several of them.
        L, n = (lx, nx) if e[0] == "X" else (ly, ny)
        lout = max(0.15 * L, 6.0 * L / n)
        I.add("Initial.SurfaceShape.EdgeBump", Order=1, CreateType="add", Base=0.5, Height=1.0,
              Location=e, Alpha=1.0, InnerLength=lout / 3.0, OuterLength=lout)
        shapes.append("EdgeBump")
    if not pinned:
        for k in range(rng.randint(1, 3)):
            kind = rng.choice(["EllipseBump", "RectBump", "PeriodicWave", "MexicanHat", "SmoothConvex"])
            order = k + 1
            ct = rng.choice(["add", "stack"])
            w = [lx * rng.uniform(0.2, 0.9), ly * rng.uniform(0.2, 0.9) if twoD else 1.0]
            if kind == "EllipseBump":
                I.add("Initial.SurfaceShape.EllipseBump", Order=order, CreateType=ct, Base=0.0,
                      Height=rng.uniform(0.2, 1.0), Center=[cx(), cy()], Width=w,
                      SlopeWidth=rng.uniform(0.05, 0.3))
            elif kind == "RectBump":
                I.add("Initial.SurfaceShape.RectBump", Order=order, CreateType=ct, Base=0.0,
                      Height=rng.uniform(0.2, 1.0), Center=[cx(), cy()], Width=w,
                      SlopeWidth=[lx * rng.uniform(0.3, 1.0), ly * rng.uniform(0.3, 1.0) if twoD else 1.0])
            elif kind == "PeriodicWave":
                #  whole wavelengths in a periodic direction
                nwx = rng.randint(1, 2)
                I.add("Initial.SurfaceShape.PeriodicWave", Order=order, CreateType=ct, Base=0.5,
                      Amplitude=rng.uniform(0.05, 0.4),
                      Wavelength=[lx / nwx, ly / rng.randint(1, 2) if twoD else 1e3],
                      Phase=[cx(), cy()])
            elif kind == "MexicanHat":
                I.add("Initial.SurfaceShape.MexicanHat", Order=order, CreateType=ct, Base=0.2,
                      Height=rng.uniform(0.2, 0.7), Center=[cx(), cy()],
                      Width=[0.5 * w[0], 0.5 * w[1]], Beta=rng.uniform(0.5, 3.0))
            else:
                I.add("Initial.SurfaceShape.SmoothConvex", Order=order, CreateType=ct, Base=0.0,
                      Height=rng.uniform(0.2, 1.0), Center=[cx(), cy()], Width=w)
            shapes.append(kind)
    m["shapes"] = shapes

    # the substrate: flat, or a low feature away from any pinned edge
    I.set("Substrate", TotalHeight=H0)
    if not pinned and rng.random() < 0.5:
        I.add("Substrate.EllipseBump", Order=1, CreateType="add", Base=0.0,
              Height=rng.uniform(0.05, 0.4), Center=[cx(), cy()],
              Width=[lx * rng.uniform(0.3, 1.0), ly * rng.uniform(0.3, 1.0) if twoD else 1.0],
              SlopeWidth=rng.uniform(0.1, 0.3))
        m["substrate"] = True

    I.set("Fluid", Name="random")
    I.set("Fluid.Density", Value=rho)
    coef = rng.uniform(-0.6, 0.6) * sig
    I.set("Fluid.SurfaceTension", BaseValue=sig, BaseConc=0.5, Coef=coef)
    I.set("Fluid.Viscosity", BaseValue=mu, BaseConc=0.5, Exponent=rng.uniform(-3.0, 10.0))
    I.set("Fluid.Diffusivity", BaseValue=logu(rng, -11, -8.5), BaseConc=0.5, Exponent=rng.uniform(-7.0, 5.0))
    I.set("Fluid.MolarMass", Solvent=logu(rng, 1.5, 2.5), Resin=logu(rng, 2.5, 4.0))

    m["conv"] = rng.random() < 0.85
    I.set("Physics.ConcField", Model=m["model"], EnableConvection=m["conv"],
          EnableDiffusion=rng.random() < 0.85)
    I.set("Physics.LaplacePressure", Enable=rng.random() < 0.9 or disj)
    m["marangoni"] = rng.random() < 0.5
    I.set("Physics.SolutoCapillarity", Enable=m["marangoni"])
    I.set("Physics.VerticalGravity", Enable=rng.random() < 0.4, Acceleration=G * rng.choice([1.0, 1.0, 0.3, -0.3]))
    hg = rng.random() < 0.3
    ang = rng.uniform(0, 2 * math.pi)
    gm = logu(rng, -1.0, 1.0)
    I.set("Physics.HorizontalGravity", Enable=hg,
          Acceleration=[gm * math.cos(ang), gm * math.sin(ang) if twoD else 0.0])
    cf = rng.random() < 0.25
    I.set("Physics.Centrifugal", Enable=cf, RotationRate=logu(rng, 0.0, 2.0),
          RotationAxis=[cx(), cy()])
    ev = rng.random() < 0.35
    I.set("Physics.Evaporation", Enable=ev, Coef=logu(rng, -9.0, -6.5),
          ConcFunc=rng.choice(["Raoult", "PowerN"]), PowerN=rng.uniform(0.3, 2.0),
          SpatialFunc="uniform", EndConcValue=0.999)
    src = (not cz) and rng.random() < 0.12
    I.set("Physics.Source", Enable=src, Radius=lx * rng.uniform(0.05, 0.2),
          InletRate=H0 * logu(rng, -1.0, 1.0) / 1.0, Conc=rng.uniform(0.05, 0.95),
          TimeConstant=0.0, ArmLength=0.0, ArmCenter=[cx(), cy()], RotationRate=0.0,
          StartAngle=0.0, RampTime=0.0, PlateauTime=1e10, SubRegions=rng.randint(1, 4))
    if disj:
        n, mm = rng.choice([(3.0, 2.0), (4.0, 3.0), (9.0, 3.0)])
        I.set("Physics.DisjoiningPressure", Enable=True, ContactAngle=rng.uniform(5.0, 30.0),
              ExponentN=n, ExponentM=mm)
    m.update(evap=ev, source=src, hgrav=hg, centrifugal=cf)
    #  A rotating film is not periodic (its arm x - x0 jumps across a seam),
    #  and the solver refuses a Periodic edge with the centrifugal term
    #  (CheckParameters): such a pair is drawn as Symmetry walls. (No draw
    #  of the generator is spent here, so what every other seed draws does
    #  not depend on it.)
    if cf:
        for e in ("Xmin", "Xmax", "Ymin", "Ymax"):
            if bc[e] == "Periodic":
                bc[e] = "Symmetry"
                I.set("Boundary.%sBC" % e, Type="Symmetry")

    # --- the numerics
    L0 = math.sqrt(sig / (rho * G))
    T0 = mu * L0 ** 4 / (sig * H0 ** 3)
    dxn = min(lx / nx, ly / ny if twoD else 1e9) / L0
    steps = rng.randint(10, 20) if quick else rng.randint(20, 50)
    adaptive = rng.random() < 0.5
    #  The capillary time of a cell, T0 dx^4 (dimensionless dx), times a
    #  factor: an adaptive step finds its own size; a fixed one that is far
    #  too large only diverges at the first step.
    dt = T0 * dxn ** 4 * (logu(rng, 0.5, 3.0) if adaptive else logu(rng, -1.0, 1.5))
    m.update(dt=dt, steps=steps, adaptive=adaptive, T0=T0, U0=L0 / T0)
    conc = rng.choice(SCHEMES)
    m["conc_scheme"] = conc
    if cz:
        I.set("Numerics.Convection", Concentration=conc)
    else:
        m["height_scheme"] = rng.choice(SCHEMES)
        I.set("Numerics.Convection", Concentration=conc, Height=m["height_scheme"])
    I.set("Numerics", InitialDt=dt, EndTime=dt * steps * 1e3, MaxSteps=steps,
          SteadyStateCriterion=-1.0)
    I.set("Numerics.TimeStepControl", Enable=adaptive, MaxDtFactor=1e10, MinDt=1e-15,
          Criterion=rng.choice([1e-2, 1e-3]), KeepSteps=rng.randint(2, 10), ChangeFactor=0.5)
    I.set("Numerics.HeightEquation", NewtonMaxIterations=30, NewtonEpsilon=1e-11)
    for eq in ("Height", "Conc"):
        I.set("Numerics.%sEquation.LinearSolver" % eq, Backend="eigen", MaxIterations=5000, Epsilon=1e-13)
    I.set("Output", Dimensioned=True, Overwrite=True)
    fint = dt * rng.randint(3, 8)
    I.set("Output.FieldData", TimeInterval=fint, NeutralFile=True, VTKFile=True)
    #  History and CalcInfo at every step's worth of time; they also cap dt
    I.set("Output.GeneralHistory", TimeInterval=dt)
    I.set("Output.CalcInfo", TimeInterval=dt)
    I.set("Output.ProbeList", TimeInterval=dt)
    m["fint"] = fint
    return I, m


# ---------------------------------------------------------------------------
#  Transformations of a 2-D input
# ---------------------------------------------------------------------------
XY_PAIRS = {"Center", "Width", "SlopeWidth", "Wavelength", "Phase", "Acceleration", "RotationAxis", "ArmCenter"}


def transposed(I):
    """x and y exchanged in everything."""
    J = json.loads(json.dumps(I.t))
    g = J["Grid"]
    g["Xmin"], g["Ymin"] = g["Ymin"], g["Xmin"]
    g["Xmax"], g["Ymax"] = g["Ymax"], g["Xmax"]
    g["Nx"], g["Ny"] = g["Ny"], g["Nx"]
    swap = {"Xmin": "Ymin", "Ymin": "Xmin", "Xmax": "Ymax", "Ymax": "Xmax"}
    bcs = {e: J["Boundary.%sBC" % e]["Type"] for e in swap}
    for e in swap:
        J["Boundary.%sBC" % e]["Type"] = bcs[swap[e]]
    for name, body in J.items():
        for b in (body if isinstance(body, list) else [body]):
            for k, v in b.items():
                if k in XY_PAIRS and isinstance(v, list):
                    b[k] = [v[1], v[0]]
                if k == "Location":
                    b[k] = swap[v]
    K = Input()
    K.t = J
    return K


def mirrored_x(I):
    """x -> -x in everything."""
    J = json.loads(json.dumps(I.t))
    g = J["Grid"]
    g["Xmin"], g["Xmax"] = -g["Xmax"], -g["Xmin"]
    a, b = J["Boundary.XminBC"]["Type"], J["Boundary.XmaxBC"]["Type"]
    J["Boundary.XminBC"]["Type"], J["Boundary.XmaxBC"]["Type"] = b, a
    for name, body in J.items():
        for bd in (body if isinstance(body, list) else [body]):
            for k, v in bd.items():
                if k in ("Center", "Phase", "Acceleration", "RotationAxis", "ArmCenter") and isinstance(v, list):
                    bd[k] = [-v[0], v[1]]
                if k == "Location" and v in ("Xmin", "Xmax"):
                    bd[k] = "Xmax" if v == "Xmin" else "Xmin"
    K = Input()
    K.t = J
    return K


# ---------------------------------------------------------------------------
#  Running and reading
# ---------------------------------------------------------------------------
def run(exe, d, src, threads=1, timeout=300):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "in.toml"), "w", encoding="utf-8") as f:
        f.write(src)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_WAIT_POLICY="PASSIVE")
    t0 = time.time()
    try:
        r = subprocess.run([exe, "in.toml"], cwd=d, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return None, "timed out after %g s" % timeout, time.time() - t0
    return r.returncode, r.stdout + "\n" + r.stderr, time.time() - t0


def table(path):
    """A History/CalcInfo file as (column names, rows of strings)."""
    with open(path, encoding="utf-8") as f:
        head = [c.split("[")[0].strip() for c in f.readline().rstrip("\n").split("\t")]
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    return head, rows


def field(path):
    """A field*.dat as {column: [floats]}; rows x fastest."""
    with open(path, encoding="utf-8") as f:
        head = [c.split("[")[0] for c in f.readline().lstrip("#").split()]
        cols = {n: [] for n in head}
        for line in f:
            v = line.split()
            if len(v) == len(head):
                for n, x in zip(head, v):
                    cols[n].append(float(x))
    return cols


def nonfinite(path):
    with open(path, "rb") as f:
        return re.search(rb"(?i)\b-?(nan|inf)\b", f.read()) is not None


def column(head, rows, name):
    k = head.index(name)
    return [float(r[k]) for r in rows]


# ---------------------------------------------------------------------------
#  The properties of one run
# ---------------------------------------------------------------------------
def interior(cols, m):
    """Indices of the rows that are cell centers (not wall faces)."""
    nxp = m["nx"] + 2
    nyp = m["ny"] + 2 if m["twoD"] else 1
    n = len(cols["xx"])
    out = []
    for r in range(n):
        i = r % nxp
        j = (r // nxp) % nyp
        if i in (0, nxp - 1):
            continue
        if m["twoD"] and j in (0, nyp - 1):
            continue
        out.append(r)
    return out


def check_run(d, m, rc, log):
    """The list of violated properties (empty when all hold)."""
    bad = []
    rd = os.path.join(d, "run")
    if rc is None:
        return ["ends: " + log]
    if re.search(r"runtime error:|AddressSanitizer|terminate called|Segmentation", log):
        bad.append("ends: sanitizer report or abort: " + log[-600:])
    if rc not in (0, 1):
        return bad + ["ends: exit code %d: %s" % (rc, log[-600:])]
    intended = None
    if rc == 1:
        intended = next((s for s in INTENDED if s in log), None)
        if not intended:
            return bad + ["ends: exit 1 without a stop the solver means: %s" % log[-800:]]
    m["stop"] = intended or "ok"

    files = sorted(glob.glob(os.path.join(rd, "field*.dat")))
    for p in files + [os.path.join(rd, "History.dat"), os.path.join(rd, "CalcInfo.dat")]:
        if os.path.exists(p) and nonfinite(p) and not intended:
            bad.append("finite: NaN or inf in %s" % os.path.basename(p))
    #  Every row of a time-series file has the columns its header names,
    #  no label more or less
    for fn in ("History.dat", "CalcInfo.dat"):
        p = os.path.join(rd, fn)
        if os.path.exists(p):
            hd, rows = table(p)
            if any(len(r) != len(hd) for r in rows):
                bad.append("files: %s has rows of %s columns under a header of %d"
                           % (fn, sorted({len(r) for r in rows}), len(hd)))
    #  The clock is the sum of the steps solved
    ji = os.path.join(rd, "JobInfo.txt")
    if os.path.exists(ji):
        vals = {}
        for line in open(ji, encoding="utf-8"):
            k, _, v = line.partition(" :")
            if k.strip() in ("Sum of the steps solved", "Time at the start of the loop",
                             "Time at the end of the loop"):
                vals[k.strip()] = float(v)
        if len(vals) == 3:
            span = vals["Time at the end of the loop"] - vals["Time at the start of the loop"]
            if span > 0 and abs(span - vals["Sum of the steps solved"]) > 1e-12 * span:
                bad.append("clock: the time run %.17g is not the sum of the steps solved %.17g"
                           % (span, vals["Sum of the steps solved"]))
    if intended:
        return bad       # the run says it broke down; what it wrote last may be broken
    if not files:
        return bad + ["ends: no field file written"]

    # film and bounds, over every field written
    floor = 0.5 * m["hstar"] if m["disj"] else 0.0
    cmin, cmax, hmin_rel = 1.0, 0.0, float("inf")
    courant = 0.0
    hd, cr = table(os.path.join(rd, "CalcInfo.dat"))
    dtmax = max(column(hd, cr, "deltaT")) if cr else m["dt"]
    dx = min(m["lx"] / m["nx"], m["ly"] / m["ny"] if m["twoD"] else 1e30)
    for p in files:
        c = field(p)
        rows = interior(c, m)
        hmax = max(c["HH"][r] for r in rows)
        hmin = min(c["HH"][r] for r in rows)
        if hmax > 0:
            hmin_rel = min(hmin_rel, hmin / hmax)
        if not (hmin > floor * (m["H0"] * 1e6)):
            bad.append("film: min h = %.3g um in %s (floor %.3g um)"
                       % (hmin, os.path.basename(p), floor * m["H0"] * 1e6))
            break
        cmin = min(cmin, min(c["CC"][r] for r in rows))
        cmax = max(cmax, max(c["CC"][r] for r in rows))
        for u in ("Ux", "Uy"):
            if u in c:
                courant = max(courant, max(abs(v) for v in c[u]) * dtmax / dx)
    m["cmin"], m["cmax"], m["courant"], m["hmin_rel"] = cmin, cmax, courant, hmin_rel
    if m["conc_scheme"] in BOUNDED and courant <= 1.0 and not (cmin >= -1e-12 and cmax <= 1 + 1e-12):
        bad.append("bounds: c in [%.6g, %.6g] with %s at Courant %.2g"
                   % (cmin, cmax, m["conc_scheme"], courant))

    # conservation
    hd, hr = table(os.path.join(rd, "History.dat"))
    closed = all(m["bc"][e] in CLOSED for e in (("Xmin", "Xmax", "Ymin", "Ymax") if m["twoD"]
                                                 else ("Xmin", "Xmax")))
    if hr:
        #  Both models write the balance (ConcZ counting what
        #  evaporates and what the 1-D centrifugal term flings sideways)
        verr = max(abs(v) for v in column(hd, hr, "VolErr"))
        serr = max(abs(v) for v in column(hd, hr, "ResinVolErr"))
        m["volerr"], m["soluteerr"] = verr, serr
        #  In 1-D the centrifugal term flings the film into the y that is
        #  not solved for; VolErr and ResinVolErr count it as outflowed, so
        #  the film is closed all the same. (The upwinded terms give a wall
        #  face no mobility, UpwindKernelBase::WallFace: with the face value
        #  of the film taken from the ghosts, an in-plane body force goes
        #  through a Symmetry or a Pinning wall.) A Periodic edge with the
        #  centrifugal term is refused by the solver, and not drawn
        #  (generate).
        #  (A wall face carries no flux and no velocity of any term,
        #  CZ_HeightEquation_2D::WallFace; a body-force flux taken from the
        #  ghosts there gives VolErr 2e-5 between Symmetry edges, seed 1015.)
        m["closed"] = closed
        #  An open edge as well: what crosses a FreshInlet or an Outlet is
        #  counted from the face fluxes of each step (counted from the
        #  velocities and face values of the field after it, the solute
        #  balance was found off by up to 1.1e-2, seed 1033; corrected)
        if not m["evap"] and not verr <= VOL_TOL:
            bad.append("volume: |VolErr| = %.3g (%s edges)" % (verr, "closed" if closed else "open"))
        #  A film that dries out too: the evaporation is limited to the
        #  solvent there is (C_Evaporation::SolventLimit), so c reaches 1
        #  without passing it and the cap cuts no solute away (seed 1206: a
        #  pinned film 95 % evaporated; a cap that took the solute above
        #  c = 1 out of the film gives ResinVolErr -0.21).
        #  With a nozzle too: ResinVolErr counts the solute it supplies
        #  (SuppliedResin; without it, seed 1122 gives ResinVolErr 0.057).
        #  (Without the advection the flow thins and thickens the film under
        #  a solute it does not carry: not a balance, and not checked.)
        if m["model"] != "CZ" and m.get("conv", True) and not serr <= SOLUTE_TOL:
            bad.append("solute: |ResinVolErr| = %.3g (%s edges)" % (serr, "closed" if closed else "open"))
        #  The well-mixed concentration without evaporation: uniform it
        #  stays uniform, the face solute flux being c_f times the volume
        #  flux of the film-height equation
        if m["model"] != "CZ" and not m["evap"]:
            drift = max(abs(v - m["c0"]) for v in field(files[-1])["CC"])
            m["wm_drift"] = drift
            if not drift <= WM_UNIFORM_TOL:
                bad.append("uniform: the uniform well-mixed concentration drifted by %.3g" % drift)
        #  ConcZ: from the field files, which carry every digit (History.dat
        #  has seven). The concentration equation has to be in conservative
        #  form, and its surface condition has to keep the conormal part of
        #  the diffusive flux: otherwise a drying ConcZ film loses 3 to 25 %
        #  of the solute its evaporated fraction should have concentrated
        #  (seeds 2038, 2299), and a uniform c drifts by ~1e-6 (seed 2104).
        if m["model"] == "CZ":
            s0, s1, drift = cz_balance(files[0], files[-1], m)
            m["cz_solute"], m["cz_drift"] = abs(s1 / s0 - 1), drift
            #  The solute flung sideways in 1-D left the field files: then
            #  the balance of History.dat, which counts it (computed to
            #  every digit, written to seven)
            if m["centrifugal"] and not m["twoD"]:
                m["cz_solute"] = serr
            #  Also once the field has reached c = 1: c is kept at or below 1
            #  without taking solute away (CZ_HeightEquation_2D::
            #  RedistributeExcess; a cap of the field alone drops the solute
            #  above it, seed 1118).
            if closed and not m["source"] and not m["cz_solute"] <= CZ_SOLUTE_TOL:
                bad.append("solute: the ConcZ solute changed by %.3g with closed edges" % m["cz_solute"])
            if not m["evap"] and not drift <= CZ_UNIFORM_TOL:
                bad.append("uniform: the uniform ConcZ concentration drifted by %.3g" % drift)
    return bad


def cz_balance(first, last, m):
    """The ConcZ solute sum over the columns of h times the mean of c over
    the cells through the film, in the first and the last field file, and
    the largest |c - c0| of the last one (every row)."""
    out = []
    for p in (first, last):
        c = field(p)
        nxp = m["nx"] + 2
        nyp = m["ny"] + 2 if m["twoD"] else 1
        n2 = nxp * nyp
        nz = len(c["xx"]) // n2
        tot = 0.0
        for r in interior(c, m):
            if r >= n2:
                continue
            col = [c["CC"][k * n2 + r] for k in range(1, nz - 1)]
            tot += c["HH"][r] * sum(col) / len(col)
        out.append(tot)
    drift = max(abs(v - m["c0"]) for v in field(last)["CC"])
    return out[0], out[1], drift


def outputs(rd):
    """The output files that must not depend on the thread count."""
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(rd, "*"))
                  if os.path.basename(p) not in ("JobInfo.txt", "in.toml", "Warnings.txt", "FatalErrors.txt"))


def check_threads(exe, d, src, timeout):
    rc, log, _ = run(exe, os.path.join(d, "t2"), src, threads=2, timeout=timeout)
    if rc != 0:
        return ["threads: the 2-thread run ended rc=%s: %s" % (rc, log[-400:])]
    a, b = os.path.join(d, "run"), os.path.join(d, "t2", "run")
    fa, fb = outputs(a), outputs(b)
    if fa != fb:
        return ["threads: different files written: %s" % sorted(set(fa) ^ set(fb))[:5]]
    for f in fa:
        with open(os.path.join(a, f), "rb") as x, open(os.path.join(b, f), "rb") as y:
            if x.read() != y.read():
                return ["threads: %s differs between 1 and 2 threads" % f]
    return []


def check_symmetry(exe, d, I, m, which, timeout):
    import equivalence as eq
    J = transposed(I) if which == "transpose" else mirrored_x(I)
    rc, log, _ = run(exe, os.path.join(d, which), J.text(), timeout=timeout)
    if rc != 0:
        return ["%s: the transformed run ended rc=%s: %s" % (which, rc, log[-400:])]
    a0, a1 = eq.last_field(os.path.join(d, "run"))
    b0, b1 = eq.last_field(os.path.join(d, which, "run"))
    if os.path.basename(a1) != os.path.basename(b1):
        return ["%s: the runs wrote different numbers of fields" % which]
    A, A0, B = eq.read_field(a1), eq.read_field(a0), eq.read_field(b1)
    names = eq.CZ_COLS if m["model"] == "CZ" else eq.WM_COLS
    if which == "transpose":
        pmap, cols = (lambda i, j, k, dd: (j, i, k)), eq.transposed(names)
        xy = lambda x, y: (y, x)
    else:
        pmap, cols = (lambda i, j, k, dd: (dd[0] - 1 - i, j, k)), eq.mirrored_x(names)
        xy = lambda x, y: (-x, y)
    why = eq.coords_agree(A, B, pmap, xy)
    if why:
        return ["%s: %s" % (which, why)]
    res = compare(A, B, A0, pmap, cols, m["U0"])
    worst = max(v for v, _ in res.values())
    m[which] = worst
    if not worst <= 1.0:
        return ["%s: fields differ by %s" % (
            which, ", ".join("%s %.1e of its %s" % (k, v * t[0], t[1]) for k, (v, t) in res.items()))]
    return []


def compare(A, B, A0, point_map, col_map, u0):
    """For each column, max |a - T(b)| over the points as a fraction of what
    it may be: 1e-7 of how far the column moved over the run (as
    tools/equivalence.py measures it), but of at least 1e-2 of its size --
    the Newton iteration stops at 1e-6 of a step's first correction, so two
    orderings of a film that hardly moves differ by about 1e-10 of it. The
    velocities are third derivatives of the film and amplify that: they are
    held to 1e-4 of the largest of their change, the largest velocity and
    1e-2 of the velocity scale L0/T0 (a film at rest has velocities of
    rounding size only)."""
    ca, nx, ny, nz = A
    cb, mx, my, mz = B
    c0 = A0[0]
    vel = [n for n in ("Ux", "Uy", "Uz") if n in ca]
    umax = max([abs(v) for n in vel for v in ca[n]] + [0.0])
    res = {}
    for name, (bname, sign) in col_map.items():
        worst, change, size = 0.0, 0.0, 0.0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    r = (k * ny + j) * nx + i
                    i2, j2, k2 = point_map(i, j, k, (nx, ny, nz))
                    r2 = (k2 * my + j2) * mx + i2
                    a, b = ca[name][r], sign * cb[bname][r2]
                    worst = max(worst, abs(a - b))
                    change = max(change, abs(a - c0[name][r]))
                    size = max(size, abs(a), abs(b))
        if name in vel:
            scale, tol, what = max(change, umax, 1e-2 * u0), 1e-4, "scale"
        elif change > 1e-2 * size:
            scale, tol, what = change, EQUIV_TOL, "change"
        else:
            scale, tol, what = 1e-2 * size, EQUIV_TOL, "size/100"
        res[name] = ((worst / scale / tol) if scale > 0 else 0.0, (tol, what))
    return res


def check_restart(exe, d, I, m, timeout):
    import restart_roundtrip as rr
    path = os.path.join(d, "chain_case")
    with open(path + ".toml", "w", encoding="utf-8") as f:
        f.write(I.text())
    steps = m["steps"]
    #  The straight run ends on EndTime (MaxSteps of a part counts from that
    #  part's start), the time the run above reached in its steps (an
    #  adaptive step may have shrunk by orders); the first two parts take a
    #  third of the steps each.
    ent = rr.pvd_entries(os.path.join(d, "run"))
    tend = float(ent[-1][0]) if ent else m["dt"] * steps
    k = max(2, steps // 3)
    v = rr.chain_case(path, {}, num(tend), num(m["fint"]), num(m["dt"]), [k, k])
    problems, summary = rr.run_chain(exe, d, "chain", v, timeout)
    if problems and problems[0].startswith("straight run") and any(s in problems[0] for s in INTENDED):
        return []       # the run without MaxSteps went further, and stopped as the solver means
    return ["restart: " + p for p in problems]


def subsets(seed, m, quick):
    """Which of the expensive properties this seed checks."""
    r = random.Random(seed * 7919 + 1)
    out = []
    #  (ConcZ with Marangoni is drawn too: a transposed run that differs by
    #  1e-5..1e-4 is not the tolerance of the coupling passes but a drift of
    #  the uniform concentration.)
    sym = m["twoD"]
    if sym and r.random() < (1.0 if quick else 0.3):
        out.append(r.choice(["transpose", "mirror"]))
    if r.random() < (0.3 if quick else 0.15):
        out.append("threads")
    #  (A pinned edge with evaporation is drawn too. Seed 1035, a film that
    #  dries to c = 1 against a pinned wall, continues from its state file
    #  as the run straight through only with the cap at 1 on the
    #  concentration equation's own state and the evaporation limited to
    #  the solvent there is; a cap on the field alone puts c 3e-3 apart an
    #  output later. Seed 1206: the check of the pinned edges is of the
    #  initial condition, and on a restart of a film whose edge has evolved
    #  it is a warning.)
    if r.random() < (0.2 if quick else 0.12):
        out.append("restart")
    return out


def one(exe, work, seed, quick, keep, timeout, force=()):
    I, m = generate(seed, quick)
    d = os.path.join(work, "s%d" % seed)
    shutil.rmtree(d, ignore_errors=True)
    src = I.text()
    rc, log, sec = run(exe, d, src, timeout=timeout)
    m["sec"] = round(sec, 2)
    m["rc"] = rc
    bad = check_run(d, m, rc, log)
    extra = list(force) or subsets(seed, m, quick)
    #  The symmetry checks replace run/ by a fixed-dt run: they go last.
    extra.sort(key=lambda x: x in ("transpose", "mirror"))
    if not bad and m.get("stop") == "ok":
        for x in extra:
            if x == "threads":
                bad += check_threads(exe, d, src, timeout)
            elif x in ("transpose", "mirror"):
                fixed = Input()
                fixed.t = json.loads(json.dumps(I.t))
                fixed.t["Numerics.TimeStepControl"]["Enable"] = False
                rc2, log2, _ = run(exe, os.path.join(d, "sym_a"), fixed.text(), timeout=timeout)
                if rc2 != 0:
                    if not any(s in log2 for s in INTENDED):
                        bad.append("%s: fixed-dt run ended rc=%s: %s" % (x, rc2, log2[-300:]))
                    continue
                shutil.rmtree(os.path.join(d, "run_adaptive"), ignore_errors=True)
                os.rename(os.path.join(d, "run"), os.path.join(d, "run_adaptive"))
                shutil.copytree(os.path.join(d, "sym_a", "run"), os.path.join(d, "run"))
                bad += check_symmetry(exe, d, fixed, m, x, timeout)
            elif x == "restart":
                bad += check_restart(exe, d, I, m, timeout)
    m["checked"] = extra if m.get("stop") == "ok" else []
    if not bad and not keep:
        shutil.rmtree(d, ignore_errors=True)
    return bad, m


QUICK_SEEDS = list(range(1, 21))
#  Seeds of the full-size generator kept as regression cases, run with
#  --quick as well; each must end normally (not merely with a stop the
#  solver means). What each one exercises, and where it is handled:
#  1002: a failed solve of the concentration, whose NaN must not stay in the
#  guess of every retry (C_SemiImplicitEquation::Solve). 1171: inlets on two
#  adjacent edges, whose corner the Newton correction must leave where it is
#  (BoundaryCondition2D::ApplyCorners). 2283: a restart from a time the .pvd
#  rounds up (0.14155857 written as 0.141559), whose entry the continued
#  .pvd must keep (C_PVD_File::Resume). 1065: in-plane gravity toward a
#  pinned edge, which must not take liquid through it (VolErr 6.4e-5
#  otherwise; UpwindKernelBase::WallFace). 1122: VolErr with a nozzle (0.033
#  if it were not counted).
#  2038, 2299: a drying ConcZ film, which must keep its solute; 2104: a
#  uniform ConcZ concentration, which must not drift, and mirror
#  symmetrically (CZ_ConcEquation_2D, the conservative form). 1015: ConcZ
#  in-plane gravity and the centrifugal flux at a Symmetry wall (VolErr 2e-5
#  if they went through it; CZ_HeightEquation_2D::WallFace). 1118, 1206: a
#  film dried to c = 1, which must keep the solute above the cap (ConcZ
#  4e-4, well-mixed ResinVolErr -0.21 if it did not;
#  C_Evaporation::SolventLimit, the enrichment of the well-mixed
#  concentration equation taken with the film-height equation's
#  evaporation, CZ_HeightEquation_2D::RedistributeExcess).
#  1035 and 1206 with the restart chain: a film dried to c = 1 against a
#  pinned wall, and one whose film at the pinned wall grows past a tenth of
#  its thickest (subsets()).
#  1054, 2105: a ConcZ film that evaporates a sizeable part of itself in one
#  step, whose solute must change by no more than rounding (1e-16): the
#  shift that puts the linear solver's residual back into the solute has to
#  include the part the shift itself changes through the surface term, a
#  quarter of it, 1e-10 (CZ_ConcEquation_2D::SolveStage). 1282: a drying
#  film with Periodic edges, whose seam takes the ghosts of the change to be
#  the copies (the same place).
REGRESSION_SEEDS = [1002, 1065, 1122, 1171, 2283, 2038, 2299, 2104, 1015, 1118, 1206, 1035,
                    1054, 2105, 1282]
REGRESSION_CHECKS = {1035: ("restart",), 1206: ("restart",)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--seed", type=int, default=1, help="first seed")
    ap.add_argument("--count", type=int, default=50)
    ap.add_argument("--only", type=int, action="append", help="run only this seed (repeatable)")
    ap.add_argument("--check", action="append", default=[],
                    help="with --only: also check transpose, mirror, threads or restart")
    ap.add_argument("--quick", action="store_true", help="the fixed slice of the ctest entry")
    ap.add_argument("--keep", action="store_true", help="keep the run directories")
    ap.add_argument("--timeout", type=float, default=300)
    ap.add_argument("--log", help="append one JSON line per seed to this file")
    ap.add_argument("--print-input", type=int, help="print the input of this seed and stop")
    a = ap.parse_args()
    if a.print_input is not None:
        sys.stdout.write(generate(a.print_input, a.quick)[0].text())
        return 0
    exe = os.path.abspath(a.exe)
    os.makedirs(a.workdir, exist_ok=True)
    seeds = a.only or (QUICK_SEEDS + REGRESSION_SEEDS if a.quick else list(range(a.seed, a.seed + a.count)))
    nbad = 0
    stops = {}
    t0 = time.time()
    logf = open(a.log, "a", encoding="utf-8") if a.log else None
    for s in seeds:
        t1 = time.time()
        regression = a.quick and s in REGRESSION_SEEDS
        bad, m = one(exe, a.workdir, s, a.quick and not regression, a.keep, a.timeout,
                     a.check if a.only else REGRESSION_CHECKS.get(s, ()) if regression else ())
        m["sec"] = round(time.time() - t1, 1)
        if regression and not bad and m.get("stop") != "ok":
            bad = ["ends: a regression seed that must end normally stopped with: %s" % m.get("stop")]
        stops[m.get("stop", "violation")] = stops.get(m.get("stop", "violation"), 0) + 1
        desc = "%s %s %dx%dx%d %s%s%s" % (m["model"], "2-D" if m["twoD"] else "1-D", m["nx"], m["ny"],
                                          m["nz"], "/".join(m["bc"][e][:3] for e in ("Xmin", "Xmax", "Ymin", "Ymax")),
                                          " +" + "+".join(m["checked"]) if m.get("checked") else "",
                                          "" if m.get("stop") in (None, "ok") else " [" + m["stop"].rstrip(":") + "]")
        if logf:
            logf.write(json.dumps({"seed": s, "bad": bad, "meta": m}, default=str) + "\n")
            logf.flush()
        if bad:
            nbad += 1
            print("FAIL seed %d  %s  (%.1f s)" % (s, desc, m.get("sec", 0)))
            for b in bad:
                print("     " + b[:1000])
        else:
            print("ok   seed %d  %s  (%.1f s)" % (s, desc, m.get("sec", 0)))
        sys.stdout.flush()
    print("%d seed(s), %d with a violation; stops: %s; %.0f s"
          % (len(seeds), nbad, ", ".join("%s %d" % kv for kv in sorted(stops.items())), time.time() - t0))
    return 1 if nbad else 0


if __name__ == "__main__":
    sys.exit(main())
