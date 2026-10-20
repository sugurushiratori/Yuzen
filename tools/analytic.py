#!/usr/bin/env python3
"""Compare the solver with problems whose answer is known in closed form.

    python3 tools/analytic.py --exe build/bin/yuzen --workdir /tmp/an
    python3 tools/analytic.py --exe build/bin/yuzen --workdir /tmp/an --only leveling_1d

Each check sets up a problem that the lubrication equations answer in closed
form -- exactly, or for the linear ones to first order in an amplitude of
1e-3 of the film -- runs the solver on it and compares. Where the error
should fall at a known rate, it is the rate that is checked, from two runs:

    leveling_1d    a small sinusoid on a flat film under capillary pressure
                    decays as exp(-sigma h0^3 k^4 t / 3mu); the error in the
                    rate must fall at second order in dx
    leveling_2d    the same in plan view, cos(kx x) cos(ky y)
    leveling_dt    the same rate, second order in dt (Crank-Nicolson)
    leveling_cz    ConcZ gives the rate well-mixed gives
    eq_slip        the equilibrium droplet of validation.md section 5.3 with
                    Navier slip of 1 um (5 times the precursor film): it
                    stays where it is (slip changes the dynamics only)
    energy_*       Output.Energy (Energy.dat) on nonlinear runs: a 30 %
                    sinusoid leveling (periodic; with slip; between Symmetry
                    walls with vertical gravity; in 2-D) and a film
                    dewetting onto its precursor: the balance of the run,
                    F(T) - F(0) + int D dt (the steps' (F^{n+1} - F^n)/dt +
                    (D^n + D^{n+1})/2 summed), relative to int D dt, at
                    second order in dt, and F never rising from step to step
    leveling_slip_1d  leveling_1d with Navier slip (Physics.Slip, lam = H0/2):
                    the rate sigma (H0^3/3 + lam H0^2) k^4 / mu; second order
                    in dx, and _dt in dt; _err the error at 32 cells
    gravity_1d      with vertical gravity the rate gains rho g h0^3 k^2 / 3mu
    potential_1d    the PHI column of WritePotential on a sinusoid under
                    capillary and disjoining pressure and vertical gravity,
                    against sigma lap h + Pi(h) - rho g h (lap the second
                    difference), in Pa
    spinodal_1d_unstable  a film thicker than 1.5 h* under capillary and
                    disjoining pressure: a long sinusoid grows at
                    omega = -(h0^3/3mu) k^2 (sigma k^2 - Pi'(h0)) > 0 (linear
                    stability); second order in dx, and _dt in dt
    spinodal_1d_stable    the same below 1.5 h*: it decays at that omega < 0
    wave_1d         with in-plane gravity the sinusoid travels at
                    rho g h0^2 / mu; the upwinding (UTOPIA) damps it by an
                    amount that falls at third order in dx
    wave_1d_speed   ... and the speed is right
    wave_2d_speed   the same along Acceleration in 2-D
    wave_*_koren    the same with the default Koren-limited scheme, which
                    clips the extrema of the wave: a larger speed error,
                    falling at about third order
    wave_cz         ConcZ with in-plane gravity: the decay rate, second order
                    in dx (ConcZ does not upwind the term)
    wave_cz_speed   ... the speed is right
    spin_1d         a flat film on a rotating substrate thins as
                    h0 / sqrt(1 + 4 rho w^2 h0^2 t / 3mu) (a strip of it,
                    an Outlet at either end)
    evap_wm         a flat film dries, h' = -E (1 - chi(c)) with c h fixed,
                    against the ODE integrated here; second order in dt
    evap_cz         the same with ConcZ, the solute kept uniform through the
                    thickness by a fast diffusion
    slab_cz         ConcZ, slow drying: the profile through the thickness
                    against diffusion in a slab with a constant flux at the
                    top (Carslaw and Jaeger's series); second order in Nz
    evap_threshold  one step of a flat film losing E dt = 0.3, 0.45, 0.75 h0:
                    the cap h / (2 dt) acts from E dt > h / 2, at either end
                    of the Crank-Nicolson step, in the film-height and the
                    concentration equations alike (and in 2-D)
    evap_solvent    the same with c0 = 0.5: the evaporation is limited to the
                    solvent there is, E dt <= h (1 - c), at either end
    dried_wm/_cz    a film dried out to c = 1 by a rate that does not fall
                    with c: the solute is kept to rounding and c does not
                    pass 1 (the check found the cap at 1 cutting the
                    solute above it; corrected);
                    _adapt: the same with dt lengthened by the time-step
                    control (the steps' last solves took the next step's dt)
    negpow_rim      a flat film with no flow and SpatialFunc "negPow": the
                    thickness lost in one step is the profile, against its
                    definition across the blend at NegPowRmax; the slope is
                    carried into the blend and the rate is constant past it
    inlet_outlet    a bump carried out through an Outlet, the film fed
                    through a FreshInlet: VolErr and ResinVolErr at
                    rounding, what crossed the edges counted from the face
                    fluxes each step put through them (counted from the
                    velocities and face values of the field after each
                    step, they were off at the order of the scheme, 1.8 in
                    dx; found, corrected); the film uniform again once
                    the bump has left
    solute_balance  the same film drying, so that c varies along it and
                    diffuses through the edges: VolErr and ResinVolErr at
                    rounding with UTOPIA, Koren, the time step adapted, and
                    without the concentration; in 2-D uniform in y, and
                    open on all four sides with every kind of open corner
    wall_balance    a closed box, a strong in-plane gravity toward its
                    walls and a centrifugal force, every convection
                    scheme: nothing crosses a Symmetry wall, |VolErr| at
                    rounding (the body forces went through it, found by
                    tools/proptest.py)
    wall_balance_pinned  the same toward a Pinning wall, the film brought
                    down to it by an EdgeBump (and the 1-D volume took the
                    ghosts in: a pinned film seemed to lose 5e-5)
    wall_balance_slip  wall_balance and wall_balance_pinned (1-D, 2-D) with
                    Navier slip
    wall_balance_cz  a closed ConcZ box with every term of the model, the
                    body forces toward the walls: volume and solute to
                    rounding, 1-D and 2-D, with evaporation or without
    conc_overshoot  a sharp front of the well-mixed concentration, 0.99 next
                    to 0.3, carried by a uniform flow (UTOPIA and Central4,
                    1-D; UTOPIA, 2-D, along a diagonal): the schemes overshoot
                    past c = 1, and the solute above it is moved to the cells
                    around, so ResinVolErr stays at rounding and c <= 1 (the
                    cap at 1 cut it, 0.7 % of the solute over one crossing
                    of the domain; found, corrected)
    conc_overshoot_koren  the same with Koren past Courant 1 (1.5), where it
                    overshoots as well: the limiter is frozen at the iterate
                    as solved, not at the field capped at 1, or the flux form
                    is lost (+1 % of the solute; found, corrected); the step
                    completed with the limiter settled at its own answer and
                    the solute formed from one flux per face, to rounding
    conc_drying     the front of conc_overshoot, 0.9 next to 0.3, on a film
                    drying by Raoult's law as it is carried (UTOPIA and Koren,
                    1-D; UTOPIA, 2-D): evaporation takes solvent only, so
                    the solute is kept, |ResinVolErr| at rounding (the
                    concentration equation solved for c lost 1e-5 to 3e-4
                    of it; found, corrected: it is the balance of h c)
    conc_uniform    a uniform concentration in a closed box with a nozzle,
                    a strong in-plane gravity, the vertical one, the
                    capillary pressure (and the spin in 1-D, which flings
                    liquid sideways): it stays uniform to rounding, the
                    face solute flux being c_f times the very volume flux
                    of the film-height equation
    conc_undershoot a front 0.7 next to 0.005 (UTOPIA, KK, Central4, 1-D;
                    UTOPIA, 2-D; ConcZ, UTOPIA, 1-D): the schemes undershoot
                    below c = 0, and the deficit is filled from the cells
                    around, so c >= 0 and |ResinVolErr| stays at rounding
                    (left below 0, the undershoot passed the divergence
                    test of 0.1 and stopped the run; found, corrected)
    spin_balance    a 1-D closed box spinning (every edge Symmetry, or one
                    Pinning): the liquid the centrifugal term flings into
                    the y that is not solved for is counted in OutflowedVol,
                    its solute in OutflowedResin, so VolErr and ResinVolErr
                    are at rounding, well-mixed (fixed and adaptive dt) and
                    ConcZ (with evaporation, c stratified, too); found
                    uncounted by tools/proptest.py, corrected
    nozzle_balance  a nozzle in a closed box, adaptive dt: VolErr counts
                    what it supplies, to rounding; ResinVolErr the solute
                    that came with it (SuppliedResin: the liquid arrives at
                    the concentration it lands on), to rounding as well
                    (found uncounted, ResinVolErr 1.9 in nozzle_1d_fixed;
                    corrected). Also nozzle_1d_rotating with its arm turning,
                    the film thickened from 10 um to 1.1 mm over 3.3 um
                    cells: found VolErr 4e-5, the capillary term's A h
                    summing over the cells to eps |A| |h| (a 1-D Symmetry
                    ghost off its mirror by rounding added to it);
                    corrected (A h of the capillary and vertical-gravity
                    terms formed face by face, mirror ghosts exact)
    nozzle_supplied the same with the liquid supplied at its own
                    concentration (Physics.Source.ConcType = "supplied"):
                    SuppliedResin is the sum of dt S Conc the solute balance
                    took in, so ResinVolErr stays at rounding (1-D, 2-D, and
                    the thick film with its arm turning); so does VolErr,
                    the thick film's stagnated Newton steps being cut (found
                    taken with a warning, VolErr 2.2e-8; corrected)
    nozzle_mixing   liquid at Conc supplied onto a flat film that does not
                    flow, over the whole domain (h = h0 + S t, c = (h0 c0 +
                    S t Conc)/(h0 + S t)) and under a small nozzle cell by
                    cell (c = (h0 c0 + (h - h0) Conc)/h), 1-D and 2-D; with
                    "local", c stays c0. Exact per step (S constant), so
                    the error is rounding
    leveling_logout  the shipped case leveling_1d_logout, written at 141
                    times spaced logarithmically over 7 decades and landed
                    on (Output.FieldData.LogSchedule, ExactTimes): the times
                    to rounding, and the decay at each as at a fixed step
    leveling_adaptive  the capillary decay of leveling_1d with the time
                    step adapted (cut back from a step too long, then
                    lengthened every third step), well-mixed and ConcZ: the
                    amplitude at the times the outputs report follows the
                    decay of the discrete operator to the error the
                    tolerance allows, and the clock is the sum of the steps
                    solved (JobInfo.txt). The check found the clock adding
                    the lengthened step, which had not been solved: the decay
                    seemed 0.8 % slow (well-mixed; ConcZ moved the film by
                    the lengthened step as well, 8e-5); corrected. With the
                    truncation estimate differencing over the steps actually
                    taken the control takes 69 steps here (385 with the wrong
                    step), for 2.9e-4

The linear checks use the Fourier coefficient of the wave number that was
put in, taken over the cell centers of a periodic domain.

Only the Python standard library is used, so the test runs wherever ctest
does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import cmath
import glob
import math
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(os.path.dirname(HERE), "cases")
sys.path.insert(0, HERE)
from concz_vs_wm import edit_toml               # noqa: E402
from equivalence import (read_field, drop_section, add_section, arr,  # noqa: E402
                         fixed_steps, run)

RHO, MU, SIG, G = 1000.0, 0.1, 30.0e-3, 9.80665


def case(name):
    with open(os.path.join(CASES, name + ".toml")) as f:
        return f.read()


def fields(d, interval):
    """(time, field) for every field*.dat. The runs here have a fixed time
    step that divides the output interval, so file n is at n * interval."""
    out = []
    for p in sorted(glob.glob(os.path.join(d, "field*.dat"))):
        n = int(re.search(r"(\d+)\.dat$", p).group(1))
        out.append((n * interval, read_field(p)))
    return out


def quiet(src, twoD=False, periodic=True):
    """Every term off but the capillary pressure; the concentration does not
    enter anything; periodic edges."""
    s = edit_toml(src, {}, {
        "Physics.SolutoCapillarity": {"Enable": "false"},
        "Physics.VerticalGravity": {"Enable": "false"},
        "Physics.HorizontalGravity": {"Enable": "false"},
        "Physics.Centrifugal": {"Enable": "false"},
        "Physics.Evaporation": {"Enable": "false"},
        "Fluid.Viscosity": {"Exponent": "0.0"},
        "Fluid.SurfaceTension": {"Coef": "0.0"},
        "Initial.SurfaceShape": {"PrecursorRatio": "0.0"},
    })
    if periodic:
        kind = '"Periodic"'
        s = edit_toml(s, {}, {"Boundary.XminBC": {"Type": kind}, "Boundary.XmaxBC": {"Type": kind},
                              "Boundary.YminBC": {"Type": kind if twoD else '"Symmetry"'},
                              "Boundary.YmaxBC": {"Type": kind if twoD else '"Symmetry"'}})
    return s


def wave_film(src, h0_frac, eps, L, Ly=None):
    """A film h = T (h0_frac + eps cos(2 pi x / L) [cos(2 pi y / Ly)]) on a
    flat substrate."""
    for sec in ("Initial.SurfaceShape.EllipseBump", "Substrate.EllipseBump"):
        src = drop_section(src, sec)
    src = add_section(src, "Initial.SurfaceShape.PeriodicWave", {
        "Order": "1", "CreateType": '"stack"', "Base": repr(h0_frac), "Amplitude": repr(eps),
        "Phase": arr(0.0, 0.0), "Wavelength": arr(L, Ly if Ly else 1.0)}, "Initial.SurfaceShape")
    return src


def mode(f, nx, ny, kx, ky, name="HH"):
    """Complex Fourier coefficient of exp(i (kx x + ky y)) over the cell
    centers (coordinates in mm, kx and ky in 1/m)."""
    s, n = 0j, 0
    for j in (range(1, ny - 1) if ny > 1 else [0]):
        for i in range(1, nx - 1):
            r = j * nx + i
            x = f["xx"][r] * 1e-3
            y = f["yy"][r] * 1e-3 if "yy" in f else 0.0
            s += f[name][r] * cmath.exp(-1j * (kx * x + ky * y))
            n += 1
    return 2 * s / n


def linear_wave(tag, workdir, exe, dim, nx, dt, steps, extra, rate, speed=(0.0, 0.0),
                ky_factor=None, model="WM", L=0.25e-3, precursor=None):
    """Run a small sinusoid and compare its amplitude and phase with
    exp(-rate t) and a translation by speed * t (m/s along x, y). A negative
    rate is a growth. `precursor` sets Initial.SurfaceShape.PrecursorRatio
    (the film is then TotalHeight (Base (1 - p) + p), the sinusoid on it)."""
    src = case("concz_1d_reg" if model == "CZ" else "leveling_1d_reg")
    s = quiet(src, twoD=(dim == 2))
    if precursor is not None:
        s = edit_toml(s, {}, {"Initial.SurfaceShape": {"PrecursorRatio": repr(precursor)}})
    s = edit_toml(s, {}, {"Grid": {"Nx": str(nx), "Xmin": repr(-L / 2), "Xmax": repr(L / 2)}})
    Ly = None
    if dim == 2:
        Ly = L * ky_factor
        ny = max(4, int(round(nx * ky_factor)))
        s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": str(ny),
                                       "Ymin": repr(-Ly / 2), "Ymax": repr(Ly / 2)}})
    s = wave_film(s, 0.5, 1e-3, L, Ly)
    s = edit_toml(s, {}, extra)
    interval = dt * steps / 4
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, tag)
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), interval)
    kx = 2 * math.pi / L
    ky = 2 * math.pi / Ly if Ly else 0.0
    (t0, (f0, nx_, ny_, _)), (t1, (f1, _, _, _)) = fs[0], fs[-1]
    #  In 2-D the initial shape is cos kx x cos ky y, i.e. the four modes
    #  (+-kx, +-ky); the (kx, ky) one is followed.
    a0 = mode(f0, nx_, ny_, kx, ky)
    a1 = mode(f1, nx_, ny_, kx, ky)
    t = t1 - t0
    got_rate = -math.log(abs(a1) / abs(a0)) / t
    #  The phase moves by -(kx cx + ky cy) t for a translation at (cx, cy)
    dphase = cmath.phase(a1 / a0)
    want_phase = -(kx * speed[0] + ky * speed[1]) * t
    err_rate = abs(got_rate - rate) / abs(rate)
    err_phase = abs(dphase - want_phase) / (abs(want_phase) if want_phase else 1.0)
    return (err_rate, err_phase, got_rate, rate, dphase, want_phase), None


def k4(L, Ly=None):
    kx = 2 * math.pi / L
    ky = 2 * math.pi / Ly if Ly else 0.0
    return (kx * kx + ky * ky) ** 2, kx * kx + ky * ky


H0 = 5e-6          # the film: TotalHeight 10 um, Base 0.5


def check_leveling_1d(workdir, exe, nx=32, dt=2e-5, steps=200):
    kk4, kk2 = k4(0.25e-3)
    rate = SIG * H0 ** 3 * kk4 / (3 * MU)
    return linear_wave("leveling_1d_n%d" % nx, workdir, exe, 1, nx, dt, steps, {}, rate)


SLIP = 2.5e-6      # the slip length of leveling_slip_1d [m]: Lambda = 0.25, half the film


def check_leveling_slip_1d(workdir, exe, nx=32, dt=2e-5, steps=200, lam=SLIP):
    """leveling_1d with Navier slip: the rate sigma (H0^3/3 + lam H0^2) k^4 / mu."""
    kk4, kk2 = k4(0.25e-3)
    rate = SIG * (H0 ** 3 / 3 + lam * H0 ** 2) * kk4 / MU
    return linear_wave("leveling_slip_1d_n%d_dt%g" % (nx, dt), workdir, exe, 1, nx, dt, steps,
                       {"Physics.Slip": {"Enable": "true", "Length": repr(lam)}}, rate)


def check_eq_slip(workdir, exe, slip=1e-6, T=0.0102):
    """The equilibrium droplet of validation.md section 5.3 (cases/eq/eq_run.toml,
    restarted from eq_n2000.dat, the steady-state test off) advanced to
    T without and with Navier slip: slip changes how fast the film moves,
    not where it rests. Returned: the largest |h(T) - h(0)| / max h with
    slip, and the same without."""
    out = []
    for lam in (None, slip):
        s = case("eq/eq_run")
        s = edit_toml(s, {}, {"Numerics": {"EndTime": repr(T), "SteadyStateCriterion": "-1"},
                              "Output.FieldData": {"VTKFile": "false"}})
        if lam:
            s = edit_toml(s, {}, {"Physics.Slip": {"Enable": "true", "Length": repr(lam)}})
        d = os.path.join(workdir, "eq_slip_%s" % ("on" if lam else "off"))
        os.makedirs(d, exist_ok=True)
        shutil.copy(os.path.join(CASES, "eq", "eq_n2000.dat"), d)
        rc, log = run(exe, d, "a", s, 1800)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        fs = sorted(glob.glob(os.path.join(d, "a", "field*.dat")))
        f0, f1 = read_field(fs[0])[0], read_field(fs[-1])[0]
        out.append(max(abs(a - b) for a, b in zip(f0["HH"], f1["HH"])) / max(f0["HH"]))
    return (out[1], out[0]), None


def read_energy(d):
    """Energy.dat as a list of rows (t, F, D, int D dt, balance)."""
    rows = []
    with open(os.path.join(d, "Energy.dat")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    return rows


def energy_src(kind, nx, dt, steps):
    """The runs of the Output.Energy checks, nonlinear (a 30 % sinusoid):
    leveling   capillary pressure, periodic
    slip       the same with Navier slip of half the film
    box        capillary pressure and vertical gravity between Symmetry walls
    dewetting  capillary and disjoining pressure (10 deg) on a film 5.5 h*
               thick, spinodally unstable: a 1 mm sinusoid of 30 % grows,
               the film ruptures onto the precursor at about 10 s and the
               liquid gathers into a drop (40 s)
    saturate   the same on a film 2 h* thick (the unstable one of
               spinodal_1d), which does not rupture: the sinusoid grows
               into a steady film, where the dissipation goes to zero
    2d         capillary pressure and vertical gravity in 2-D, periodic"""
    twoD = kind == "2d"
    L = SPIN_L if kind in ("dewetting", "saturate") else 0.25e-3
    s = quiet(case("leveling_2d_reg" if twoD else "leveling_1d_reg"), twoD=twoD, periodic=(kind != "box"))
    if kind in ("dewetting", "saturate"):
        s = edit_toml(s, {}, {"Initial.SurfaceShape": {"PrecursorRatio": "0.1" if kind == "dewetting" else repr(1.0 / 3.0)},
                              "Physics.DisjoiningPressure": {"Enable": "true", "ContactAngle": repr(SPIN_THETA),
                                                             "ExponentN": "3.0", "ExponentM": "2.0"}})
    s = edit_toml(s, {}, {"Grid": {"Nx": str(nx), "Xmin": repr(-L / 2), "Xmax": repr(L / 2)}})
    Ly = None
    if twoD:
        Ly = 1.5 * L
        s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": str(int(round(1.5 * nx))),
                                       "Ymin": repr(-Ly / 2), "Ymax": repr(Ly / 2)}})
    s = wave_film(s, 0.5, 0.15, L, Ly)
    if kind in ("box", "2d"):
        s = edit_toml(s, {}, {"Physics.VerticalGravity": {"Enable": "true", "Acceleration": "5000.0"}})
    if kind == "slip":
        s = edit_toml(s, {}, {"Physics.Slip": {"Enable": "true", "Length": repr(SLIP)}})
    s = fixed_steps(s, dt, steps, dt * steps / 8)
    #  A row of Energy.dat at every step, so that a rise of F between two
    #  steps is seen
    return edit_toml(s, {}, {"Output.Energy": {"Enable": "true"},
                             "Output.GeneralHistory": {"TimeInterval": repr(dt)}})


def check_energy(workdir, exe, kind, nx, dt, steps):
    """One run: (largest rise of F between rows relative to F's whole fall,
    the balance of the run, |F(T) - F(0) + int_0^T D dt| / int_0^T D dt,
    the fall of F relative to its start, the largest |Balance| of a step
    relative to the largest D)."""
    d = os.path.join(workdir, "energy_%s_n%d_dt%g" % (kind, nx, dt))
    rc, log = run(exe, d, "a", energy_src(kind, nx, dt, steps), 1800)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    rows = read_energy(os.path.join(d, "a"))
    F = [r[1] for r in rows]
    fall = F[0] - min(F)
    rise = max(max(b - a for a, b in zip(F, F[1:])), 0.0) / fall
    #  A negative dissipation is a rise of F in the making
    if min(r[2] for r in rows) < 0:
        rise = max(rise, -min(r[2] for r in rows) / max(r[2] for r in rows))
    dmax = max(r[2] for r in rows)
    run_bal = abs(F[-1] - F[0] + rows[-1][3]) / rows[-1][3]
    step_bal = max(abs(r[4]) for r in rows[1:]) / dmax
    return (rise, run_bal, fall / abs(F[0]) if F[0] else fall, step_bal), None


_ENERGY = {}


def energy_order(kind, nx, dt, steps):
    """The balance (F^{n+1} - F^n)/dt + (D^n + D^{n+1})/2, relative to D,
    at dt and dt/2 over the same time: its order in dt; and F never rises.
    Made once for each set of arguments."""
    key = (kind, nx, dt, steps)
    if key not in _ENERGY:
        _ENERGY[key] = _energy_order(kind, nx, dt, steps)
    return _ENERGY[key]


def _energy_order(kind, nx, dt, steps):
    """The balance of the whole run at dt and dt/2: its order."""
    a, why = check_energy(W, X, kind, nx, dt, steps)
    if a is None:
        raise RuntimeError(why)
    b, why = check_energy(W, X, kind, nx, dt / 2, 2 * steps)
    if b is None:
        raise RuntimeError(why)
    return a, b, order(a[1], b[1])


def check_leveling_2d(workdir, exe, nx=24, dt=2e-5, steps=100):
    kk4, kk2 = k4(0.25e-3, 0.25e-3 * 1.5)
    rate = SIG * H0 ** 3 * kk4 / (3 * MU)
    return linear_wave("leveling_2d_n%d" % nx, workdir, exe, 2, nx, dt, steps, {}, rate,
                       ky_factor=1.5)


def check_gravity_1d(workdir, exe, nx=32, dt=2e-5, steps=200):
    kk4, kk2 = k4(0.25e-3)
    g = 5000.0
    rate = SIG * H0 ** 3 * kk4 / (3 * MU) + RHO * g * H0 ** 3 * kk2 / (3 * MU)
    return linear_wave("gravity_1d_n%d" % nx, workdir, exe, 1, nx, dt, steps,
                       {"Physics.VerticalGravity": {"Enable": "true", "Acceleration": repr(g)}}, rate)


#  The linear stability of a flat film under capillary and disjoining
#  pressure. With Pi(h) = B [(h*/h)^n - (h*/h)^m], B = (n-1)(m-1) sigma
#  theta_e^2 / (2 h* (n-m)) (the calibration of NonDimensionalNumbers), a
#  sinusoid of wavenumber k on a film of mean thickness h0 grows at
#
#      omega = -(h0^3 / 3 mu) k^2 (sigma k^2 - Pi'(h0)),
#      Pi'(h) = -B [n (h*/h)^n - m (h*/h)^m] / h,
#
#  (in the dimensionless form of the instructions, -H0^3 k^2 (k^2 + Bo -
#  Dp Pi'(H0)) / 3, Bo = 0 here), so the film is spinodally unstable at long
#  wavelengths where Pi'(h0) > 0, i.e. h0 > h* (n/m)^(1/(n-m)) = 1.5 h* for
#  (3, 2). The film of the check is TotalHeight (Base (1 - p) + p) with
#  Base 0.5 and p the PrecursorRatio: h*/h0 = 1/2 for p = 1/3 (unstable),
#  4/5 for p = 2/3 (stable at every wavelength).
SPIN_L = 1.0e-3      # the wavelength
SPIN_THETA = 10.0    # ContactAngle [deg]


def spinodal_rate(p, L=SPIN_L, theta=SPIN_THETA, n=3.0, m=2.0, total=10e-6):
    """The growth rate omega (> 0: the sinusoid grows)."""
    hstar = p * total
    h0 = total * (0.5 * (1.0 - p) + p)
    th = theta * math.pi / 180.0
    B = (n - 1) * (m - 1) * SIG * th * th / (2.0 * hstar * (n - m))
    r = hstar / h0
    dpi = -B * (n * r ** n - m * r ** m) / h0
    k2 = (2 * math.pi / L) ** 2
    return -(h0 ** 3 / (3 * MU)) * k2 * (SIG * k2 - dpi)


def check_spinodal(workdir, exe, p, nx=32, dt=1e-3, steps=4000):
    rate = -spinodal_rate(p)        # linear_wave takes a decay rate
    extra = {"Physics.DisjoiningPressure": {"Enable": "true", "ContactAngle": repr(SPIN_THETA),
                                            "ExponentN": "3.0", "ExponentM": "2.0"}}
    return linear_wave("spinodal_p%.3f_n%d_dt%g" % (p, nx, dt), workdir, exe, 1, nx, dt, steps,
                       extra, rate, L=SPIN_L, precursor=p)


def check_potential(workdir, exe, nx=32, p=1.0 / 3.0, g=500.0):
    """Output.FieldData.WritePotential: the PHI column of the initial field
    (a 1 mm sinusoid on a film under capillary and disjoining pressure and
    vertical gravity) against sigma lap(h) + Pi(h) - rho g h at the cell
    centers, lap the second difference: the relative error over the cells."""
    src = case("leveling_1d_reg")
    s = quiet(src)
    s = edit_toml(s, {}, {"Initial.SurfaceShape": {"PrecursorRatio": repr(p)},
                          "Grid": {"Nx": str(nx), "Xmin": repr(-SPIN_L / 2), "Xmax": repr(SPIN_L / 2)},
                          "Physics.DisjoiningPressure": {"Enable": "true", "ContactAngle": repr(SPIN_THETA),
                                                         "ExponentN": "3.0", "ExponentM": "2.0"},
                          "Physics.VerticalGravity": {"Enable": "true", "Acceleration": repr(g)},
                          "Output.FieldData": {"WritePotential": "true"}})
    s = wave_film(s, 0.5, 1e-2, SPIN_L)
    s = fixed_steps(s, 1e-6, 1, 1e-6)
    d = os.path.join(workdir, "potential_1d")
    rc, log = run(exe, d, "a", s, 300)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    path = os.path.join(d, "a", "field00000.dat")
    with open(path) as f:
        head = f.readline()
    if "PHI[Pa]" not in head:
        return None, "no PHI[Pa] column in " + path
    f0, n_, _, _ = read_field(path)
    total = 10e-6
    hstar = p * total
    th = SPIN_THETA * math.pi / 180.0
    B = 2.0 * SIG * th * th / (2.0 * hstar)     # (n-1)(m-1)/(n-m) = 2 for (3, 2)
    xs = [x * 1e-3 for x in f0["xx"]]
    hs = [h * 1e-6 for h in f0["HH"]]
    phi = f0["PHI"]
    err, scale = 0.0, 0.0
    for i in range(2, len(xs) - 2):
        dxl, dxr = xs[i] - xs[i - 1], xs[i + 1] - xs[i]
        lap = 2.0 * ((hs[i + 1] - hs[i]) / dxr - (hs[i] - hs[i - 1]) / dxl) / (dxl + dxr)
        r = hstar / hs[i]
        want = SIG * lap + B * (r ** 3 - r ** 2) - RHO * g * hs[i]
        err = max(err, abs(phi[i] - want))
        scale = max(scale, abs(want))
    return err / scale, None


def scheme(height):
    """Numerics.Convection for the upwinded fluxes of the film height, or
    nothing for the default"""
    return {"Numerics.Convection": {"Height": '"%s"' % height}} if height else {}


def check_wave_1d(workdir, exe, nx=32, dt=2e-5, steps=200, height=None):
    kk4, kk2 = k4(0.25e-3)
    g = 5.0e4
    rate = SIG * H0 ** 3 * kk4 / (3 * MU)
    c = RHO * g * H0 ** 2 / MU
    extra = {"Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(g, 0.0)}}
    extra.update(scheme(height))
    return linear_wave("wave_1d_%s_n%d" % (height or "default", nx), workdir, exe, 1, nx, dt, steps,
                       extra, rate, speed=(c, 0.0))


def check_wave_2d(workdir, exe, nx=24, dt=2e-5, steps=100, height=None):
    kk4, kk2 = k4(0.25e-3, 0.25e-3 * 1.5)
    gx, gy = 4.0e4, 3.0e4
    rate = SIG * H0 ** 3 * kk4 / (3 * MU)
    c = (RHO * gx * H0 ** 2 / MU, RHO * gy * H0 ** 2 / MU)
    extra = {"Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(gx, gy)}}
    extra.update(scheme(height))
    return linear_wave("wave_2d_%s_n%d" % (height or "default", nx), workdir, exe, 2, nx, dt, steps,
                       extra, rate, speed=c, ky_factor=1.5)


def check_leveling_cz(workdir, exe, nx=32, dt=2e-5, steps=200):
    kk4, kk2 = k4(0.25e-3)
    rate = SIG * H0 ** 3 * kk4 / (3 * MU)
    return linear_wave("leveling_cz_n%d" % nx, workdir, exe, 1, nx, dt, steps, {}, rate,
                       model="CZ")


def check_wave_cz(workdir, exe, nx=32, dt=2e-5, steps=200):
    kk4, kk2 = k4(0.25e-3)
    g = 5.0e4
    rate = SIG * H0 ** 3 * kk4 / (3 * MU)
    c = RHO * g * H0 ** 2 / MU
    return linear_wave("wave_cz_n%d" % nx, workdir, exe, 1, nx, dt, steps,
                       {"Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(g, 0.0)}},
                       rate, speed=(c, 0.0), model="CZ")


def flat(src, extra, model="WM", nx=10):
    """A flat film of h0 = 5 um, no precursor, periodic in x, the given terms
    on top of the capillary pressure (which does nothing on a flat film)."""
    s = quiet(src)
    s = wave_film(s, 0.5, 0.0, 1.0)
    s = edit_toml(s, {}, dict({"Grid": {"Nx": str(nx)},
                               "Physics.ConcField": {"Model": '"%s"' % model}}, **extra))
    return s


def check_spin_1d(workdir, exe, dt=1e-3, steps=300):
    """A strip of a flat film that goes on beyond the domain (an Outlet at
    either end: a rotating film is not periodic, and the solver refuses a
    Periodic edge with the centrifugal term)."""
    w = 2000.0
    s = flat(case("leveling_1d_reg"), {"Physics.Centrifugal": {
        "Enable": "true", "RotationRate": repr(w), "RotationAxis": arr(0.0, 0.0)},
        "Boundary.XminBC": {"Type": '"Outlet"'}, "Boundary.XmaxBC": {"Type": '"Outlet"'}})
    interval = dt * steps / 3
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, "spin_1d")
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    worst = 0.0
    for t, (f, nx, ny, nz) in fields(os.path.join(d, "a"), interval):
        want = H0 / math.sqrt(1.0 + 4.0 * RHO * w * w * H0 * H0 * t / (3.0 * MU))
        for i in range(1, nx - 1):
            worst = max(worst, abs(f["HH"][i] * 1e-6 - want) / (H0 - want if t else 1.0))
    return worst, None


def evap_ode(c0, mr, E0, T, n=200000):
    """h(t) of a drying flat film, h' = -E0 (1 - chi(c)), c h = c0 h0,
    by RK4 in n steps; returns (h/h0 at T, c at T)."""
    def rhs(h):
        c = c0 / h
        chi = (c / mr) / (c / mr + (1.0 - c))
        return -E0 * (1.0 - chi) / H0
    h, dt = 1.0, T / n
    for _ in range(n):
        k1 = rhs(h); k2 = rhs(h + 0.5 * dt * k1); k3 = rhs(h + 0.5 * dt * k2); k4_ = rhs(h + dt * k3)
        h += dt * (k1 + 2 * k2 + 2 * k3 + k4_) / 6.0
    return h, c0 / h


def check_evap(workdir, exe, model, dt=1e-3, steps=400):
    E0 = 5e-6
    extra = {"Physics.Evaporation": {"Enable": "true", "Coef": repr(E0),
                                     "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'},
             "Fluid.Diffusivity": {"BaseValue": "1e-6", "Exponent": "0.0"}}
    base = case("concz_1d_reg") if model == "CZ" else case("leveling_1d_reg")
    s = flat(base, extra, model)
    interval = dt * steps / 2
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, "evap_" + model)
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    mr = 2000.0 / 100.0
    t, (f, nx, ny, nz) = fields(os.path.join(d, "a"), interval)[-1]
    hr, cr = evap_ode(0.5, mr, E0, t)
    i = nx // 2
    h = f["HH"][i] * 1e-6 / H0
    c = (sum(f["CC"][k * nx + i] for k in range(1, nz - 1)) / (nz - 2)) if model == "CZ" else f["CC"][i]
    return max(abs(h - hr) / (1 - hr), abs(c - cr) / (cr - 0.5)), None


def check_evap_wm(workdir, exe, **kw):
    return check_evap(workdir, exe, "WM", **kw)


def check_evap_cz(workdir, exe, **kw):
    return check_evap(workdir, exe, "CZ", **kw)


def check_slab_cz(workdir, exe, nz=20, dt=1e-3, steps=200):
    """Slow drying of a flat ConcZ film: the solute that the receding
    surface leaves behind diffuses down into a slab with no flux at the
    bottom. To first order in E h / D the profile relative to its mean is

        (q h / D) [ (3 zeta^2 - 1)/6 - 2/pi^2 sum (-1)^n/n^2 exp(-n^2 pi^2 D t/h^2) cos(n pi zeta) ]

    with q = c E the flux (Carslaw and Jaeger, the slab with a constant
    flux at one face). E h / D is 5e-5 here, so the linearization costs
    less than that."""
    E0, D = 1e-9, 1e-10
    extra = {"Grid": {"Nz": str(nz)},
             "Physics.Evaporation": {"Enable": "true", "Coef": repr(E0),
                                     "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'},
             "Fluid.Diffusivity": {"BaseValue": repr(D), "Exponent": "0.0"}}
    s = flat(case("concz_1d_reg"), extra, "CZ")
    interval = dt * steps / 2
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, "slab_cz_nz%d" % nz)
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), interval)
    (t0, (f0, nx, ny, nz_)), (t1, (f1, _, _, _)) = fs[0], fs[-1]
    i = nx // 2
    h0, h1 = f0["HH"][i] * 1e-6, f1["HH"][i] * 1e-6
    col = [f1["CC"][k * nx + i] for k in range(nz_)]
    cells = col[1:-1]
    mean = sum(cells) / len(cells)
    q = mean * (h0 - h1) / t1              # the flux the surface actually took
    h = 0.5 * (h0 + h1)
    worst, amp = 0.0, q * h / D
    for k, c in enumerate(cells):
        z = (k + 0.5) / len(cells)
        ser = sum((-1) ** n / n ** 2 * math.exp(-n * n * math.pi ** 2 * D * t1 / (h * h))
                  * math.cos(n * math.pi * z) for n in range(1, 200))
        want = amp * ((3 * z * z - 1) / 6 - 2 / math.pi ** 2 * ser)
        worst = max(worst, abs((c - mean) - want))
    return worst / (amp / 3), None


def check_nozzle_arm(workdir, exe, dt=1e-3):
    """A nozzle on a rotating arm over a film with no flow at all: every point
    just collects what falls on it, so the centroid of what arrived between
    two outputs is the mean position of the nozzle over that interval.

    The arm (length 0.8 mm, pivot at the origin, 1 rad/s) ramps up over
    Ts = 0.3 s, turns steadily for Tm = 0.5 s and ramps down over another
    0.3 s. Its angle is the integral of the rate, 0.15 rad at the end of the
    ramp-up, and it stops at omega (Ts + Tm) = 0.8 rad. Returns the largest
    error of the centroid over the steady turn [0.3, 0.6] s and after the
    stop [1.2, 1.5] s, relative to the arm length."""
    L, w, Ts, Tm = 0.8e-3, 1.0, 0.3, 0.5
    s = flat(case("leveling_1d_reg"), {
        "Grid": {"Nx": "400", "Xmin": "-1.0e-3", "Xmax": "1.0e-3"},
        "Physics.LaplacePressure": {"Enable": "false"},
        "Physics.Source": {"Enable": "true", "Radius": "20e-6", "InletRate": "1e-5",
                           "TimeConstant": "0.0", "ArmLength": repr(L),
                           "ArmCenter": arr(0.0, 0.0), "RotationRate": repr(w),
                           "StartAngle": "0.0", "RampTime": repr(Ts),
                           "PlateauTime": repr(Tm), "SubRegions": "8"}}, nx=400)
    s = edit_toml(s, {}, {"Grid": {"Xmin": "-1.0e-3", "Xmax": "1.0e-3"}})
    interval, steps = 0.3, int(round(1.5 / dt))
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, "nozzle_arm")
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), interval)

    def centroid(k0, k1):
        (_, (f0, nx, _, _)), (_, (f1, _, _, _)) = fs[k0], fs[k1]
        num = den = 0.0
        for i in range(1, nx - 1):
            dh = f1["HH"][i] - f0["HH"][i]
            num += f0["xx"][i] * 1e-3 * dh
            den += dh
        return num / den

    steady = L * (math.sin(Ts / 2 + w * 0.3) - math.sin(Ts / 2)) / (w * 0.3)
    stop = L * math.cos(w * (Ts + Tm))
    return max(abs(centroid(1, 2) - steady), abs(centroid(4, 5) - stop)) / L, None


def dry_profile_reference(h0, c0, E0, mr, D0, De, c0D, T, n=600, steps=3000):
    """The concentration through a flat drying film, by a finite-volume
    solution of its own: in zeta = z/h the solute per cell is conserved,

        (h c)_t = d/dzeta [ h_t zeta c + (D(c)/h) c_zeta ],

    with no flux at the bottom or at the surface (the solvent leaves, the
    solute stays: D c_z = c E at the surface), h_t = -E0 (1 - chi(c_s)), and
    D = D0 exp(De (c - c0D)). Implicit midpoint in time, each step solved by
    iterating the tridiagonal system; n cells. Returns (h(T), c at the cell
    centers)."""
    dz = 1.0 / n
    c = [c0] * n
    h = h0
    dt = T / steps

    def evap(cs):
        chi = (cs / mr) / (cs / mr + (1.0 - cs))
        return E0 * (1.0 - chi)

    def surface(cc, hh):
        #  c_s from the top cell and the flux condition (D/h) c_zeta = c_s E
        cs = cc[-1]
        for _ in range(3):
            D = D0 * math.exp(De * (cs - c0D))
            cs = cc[-1] / (1.0 - 0.5 * dz * hh * evap(cs) / D)
        return cs

    for _ in range(steps):
        cmid, hmid = c[:], h
        for _it in range(3):
            cs = surface(cmid, hmid)
            ht = -evap(cs)
            hnew = h + dt * ht
            hm = 0.5 * (h + hnew)
            #  (hnew c_new - h c_old)/dt = 0.5 (L c_new + L c_old), L frozen at the midpoint
            a_, b_, cc_, r_ = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
            coef = []
            for j in range(n - 1):
                zf = (j + 1) * dz
                Df = D0 * math.exp(De * (0.5 * (cmid[j] + cmid[j + 1]) - c0D))
                coef.append((ht * zf * 0.5, Df / (hm * dz)))      # F = u (c_j + c_j+1) + k (c_j+1 - c_j)
            def L(v, j):
                out = 0.0
                if j < n - 1:
                    u, k = coef[j]; out += u * (v[j] + v[j + 1]) + k * (v[j + 1] - v[j])
                if j > 0:
                    u, k = coef[j - 1]; out -= u * (v[j - 1] + v[j]) + k * (v[j] - v[j - 1])
                return out / dz
            for j in range(n):
                b_[j] = hnew / dt
                r_[j] = h * c[j] / dt + 0.5 * L(c, j)
                if j < n - 1:
                    u, k = coef[j]; b_[j] -= 0.5 * (u - k) / dz; cc_[j] = -0.5 * (u + k) / dz
                if j > 0:
                    u, k = coef[j - 1]; b_[j] -= 0.5 * (-u - k) / dz; a_[j] = -0.5 * (-u + k) / dz
            #  Thomas
            cp, dp = [0.0] * n, [0.0] * n
            cp[0] = cc_[0] / b_[0]; dp[0] = r_[0] / b_[0]
            for j in range(1, n):
                m = b_[j] - a_[j] * cp[j - 1]
                cp[j] = cc_[j] / m if j < n - 1 else 0.0
                dp[j] = (r_[j] - a_[j] * dp[j - 1]) / m
            cn = [0.0] * n
            cn[-1] = dp[-1]
            for j in range(n - 2, -1, -1):
                cn[j] = dp[j] - cp[j] * cn[j + 1]
            cmid = [0.5 * (x + y) for x, y in zip(c, cn)]
            hmid = hm
        c, h = cn, hnew
    return h, c


def check_dried_out(workdir, exe, model, twoD=False, a=0.6, dt=0.01, steps=40, adaptive=False):
    """A flat film dried out to c = 1: a rate that does not fall with c
    (PowerN with exponent 0, which never stops by itself) taking E dt = a h0
    a step, c0 = 0.5, 40 steps. The evaporation is limited to the solvent
    there is (C_Evaporation::SolventLimit), so the film approaches c = 1 --
    a third of the solvent left after each step once the limit acts -- and
    reaches it to rounding without passing it. Returns (the change of the
    solute sum h mean(c) between the first and the last field, the largest
    c - 1, 1 - the smallest c of the last field). A cap of the field at
    c = 1 that took the solute above it away would show as a loss here.

    adaptive: the time-step control on, so that dt is lengthened as the film
    dries (the run must say "dt+"). TimeStepAdjustment lengthens DT().now
    before the step is finished, and the last solve of the well-mixed
    concentration once took that, the next step's dt, with the evaporation
    of this one: the solute grew at every lengthening (ResinVolErr 1e-5 in
    wet_coupled)."""
    extra = {"Physics.Evaporation": {"Enable": "true", "Coef": repr(a * H0 / dt),
                                     "ConcFunc": '"PowerN"', "PowerN": "0.0", "ConcR": "0.0",
                                     "SpatialFunc": '"uniform"',
                                     #  (no dry-out stop: c does not pass 1)
                                     "EndConcType": '"min"', "EndConcValue": "1.0"},
             "Fluid.Diffusivity": {"BaseValue": "1e-9", "Exponent": "0.0"}}
    if twoD:
        extra["Grid"] = {"Nx": "6", "Ny": "6"}
    if model == "CZ":
        base = case("concz_2d" if twoD else "concz_1d_reg")
        extra.setdefault("Grid", {})["Nz"] = "8"
    else:
        base = case("leveling_2d_reg" if twoD else "leveling_1d_reg")
    s = flat(base, extra, model, nx=6 if twoD else 10)
    s = fixed_steps(s, dt, steps, dt * steps)
    if adaptive:
        s = edit_toml(s, {}, {"Numerics.TimeStepControl": {"Enable": "true", "KeepSteps": "2",
                                                            "ChangeFactor": "0.75", "Criterion": "0.5",
                                                            "MaxDtFactor": "1e10", "MinDt": "1e-15"}})
    d = os.path.join(workdir, "dried_%s_%s%s" % (model, "2d" if twoD else "1d", "_adapt" if adaptive else ""))
    rc, log = run(exe, d, "a", s, 600)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    if adaptive and "dt+" not in log:
        return None, "the time step was never lengthened: " + log[-200:]
    fs = fields(os.path.join(d, "a"), dt * steps)
    sol = []
    for _, (f, nx, ny, nz) in (fs[0], fs[-1]):
        cmax, cmin = 0.0, 2.0       # (of the last field, in the end)
        cols = [(i, j) for j in (range(1, ny - 1) if ny > 1 else [0]) for i in range(1, nx - 1)]
        n2 = nx * ny
        tot = 0.0
        for i, j in cols:
            r = j * nx + i
            cs = [f["CC"][k * n2 + r] for k in range(1, nz - 1)] if model == "CZ" else [f["CC"][r]]
            tot += f["HH"][r] * sum(cs) / len(cs)
            cmax, cmin = max(cmax, max(cs)), min(cmin, min(cs))
        sol.append(tot)
    return (abs(sol[1] / sol[0] - 1.0), max(cmax - 1.0, 0.0), 1.0 - cmin), None


def check_dry_cz(workdir, exe, nz=16, dt=1e-4, T=0.05):
    """A flat ConcZ film drying at a Peclet number E h/D of order 1, so that
    the solute piles up under the surface -- the regime slab_cz, at 5e-5,
    stays out of -- with a diffusivity that depends on the concentration.
    The profile the vertical advection of the zeta mapping and the diffusion
    shape together is compared with dry_profile_reference, a solution of the
    same equations by other means on 600 cells. Returns the largest error
    over the cells relative to the rise of the surface concentration."""
    E0, D0, De = 2e-5, 1e-10, -3.0
    extra = {"Grid": {"Nz": str(nz)},
             "Physics.Evaporation": {"Enable": "true", "Coef": repr(E0),
                                     "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'},
             "Fluid.Diffusivity": {"BaseValue": repr(D0), "Exponent": repr(De), "BaseConc": "0.5"}}
    s = flat(case("concz_1d_reg"), extra, "CZ")
    steps = int(round(T / dt))
    s = fixed_steps(s, dt, steps, T)
    d = os.path.join(workdir, "dry_cz_nz%d" % nz)
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), T)
    (t0, (f0, nx, ny, nz_)), (t1, (f1, _, _, _)) = fs[0], fs[-1]
    i = nx // 2
    h0 = f0["HH"][i] * 1e-6
    cells = [f1["CC"][k * nx + i] for k in range(nz_)][1:-1]
    href, cref = dry_profile_reference(h0, 0.5, E0, 2000.0 / 100.0, D0, De, 0.5, t1)
    n = len(cref)

    def at(z):
        x = z * n - 0.5
        j = min(max(int(math.floor(x)), 0), n - 2)
        w = x - j
        return (1 - w) * cref[j] + w * cref[j + 1]
    rise = max(cref) - 0.5
    worst = max(abs(c - at((k + 0.5) / len(cells))) for k, c in enumerate(cells))
    herr = abs(f1["HH"][i] * 1e-6 - href) / (h0 - href)
    return max(worst / rise, herr), None


def check_evap_cap(workdir, exe):
    """A flat film drying so fast that one step would take more than half of
    it, and more than the solvent it holds (c0 = 0.5): the evaporation is
    capped at h / (2 dt) and at the solvent h (1 - c) / dt, at the start of
    the step and at its end alike. The start takes 0.5 h0 (both caps), the
    end h1 - c0 h0 (the solvent), so h1 = h0 - (0.5 h0 + h1 - 0.5 h0) / 2,
    h1 = 2/3 h0 (evap_cap_step) -- and the solute h c must not change. (The
    cap of h / 2 alone gives 0.6 h0; a concentration equation that enriched
    the solute at the uncapped rate would make 30 % more solute in this one
    step.) Returns the larger of the two errors."""
    E0, dt = 1e-3, 0.01
    s = flat(case("leveling_1d_reg"), {
        "Physics.Evaporation": {"Enable": "true", "Coef": repr(E0),
                                "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'},
        "Fluid.Diffusivity": {"BaseValue": "1e-6", "Exponent": "0.0"}}, "WM")
    s = fixed_steps(s, dt, 1, dt)
    d = os.path.join(workdir, "evap_cap")
    rc, log = run(exe, d, "a", s, 300)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), dt)
    (_, (f0, nx, _, _)), (_, (f1, _, _, _)) = fs[0], fs[1]
    i = nx // 2
    h0, h1 = f0["HH"][i], f1["HH"][i]
    wh, _ = evap_cap_step(1e3, 0.5)
    #  h1 to where the Newton iteration stops, the solute to rounding (as
    #  in check_evap_threshold)
    return max(abs(h1 / h0 - wh), 1e6 * abs(h1 * f1["CC"][i] / (h0 * f0["CC"][i]) - 1.0)), None


def evap_cap_step(a, c0):
    """One Crank-Nicolson step of a flat film losing E dt = a h0, the rate
    independent of c, capped at h / (2 dt) and at the solvent h (1 - c) / dt
    at the start of the step (h0, c0) and at its end (h1, c1) alike. Returns
    (h1 / h0, c1 / c0): the film height from h1 = h0 - dt (e(h0) + e(h1)) / 2
    (solved by fixed point; the map contracts by 1/2), and c1 / c0 = h0 / h1:
    the concentration equation enriches the solute with the evaporation the
    film-height equation took out, over the film it ends on, so the solute h c
    of a film that only dries is kept exactly."""
    e0 = min(a, 0.5, 1.0 - c0)
    h1 = 1.0
    for _ in range(200):
        h1 = 1.0 - 0.5 * e0 - 0.5 * min(a, 0.5 * h1, max(h1 - c0, 0.0))
    return h1, 1.0 / h1


def check_evap_threshold(workdir, exe, a, twoD=False, c0=0.1):
    """The evaporation caps' thresholds: a flat film losing E dt = a h0 in
    one step, the rate independent of c (PowerN with exponent 0), against
    evap_cap_step(a, c0). With c0 = 0.1 the solvent limit stays out of it and
    the cap of h / 2 is tested: a = 0.3 is never capped (h1 = 0.7 h0); 0.45 is
    below the threshold h / 2 at the start of the step and above it at the
    end (0.62 h0 -- a threshold of h would give 0.55); 0.75, between h / 2 and
    h, is above it at both (0.6 h0 -- a threshold of h would give 0.5). With
    c0 = 0.5 the solvent is the tighter limit: 0.3 passes it at the end of
    the step (h1 = 0.7333 h0, against 0.7 without it), 0.75 at both ends
    (2/3 h0). The concentration after the step is h0 / h1 times c0 in every
    case: the solute h c is kept exactly, because the concentration
    equation takes what the film-height equation took out (with the rate
    itself, implicitly, it is kept only to the order of the scheme below the
    cap, 1.0245 at a = 0.3). Returns the larger of the two errors."""
    dt = 0.01
    extra = {"Physics.Evaporation": {"Enable": "true", "Coef": repr(a * H0 / dt),
                                     "ConcFunc": '"PowerN"', "PowerN": "0.0",
                                     "SpatialFunc": '"uniform"'},
             "Fluid.Diffusivity": {"BaseValue": "1e-6", "Exponent": "0.0"},
             "Initial.Concentration": {"Value": repr(c0)}}
    if twoD:
        extra["Grid"] = {"Nx": "6", "Ny": "6"}
    s = flat(case("leveling_2d_reg" if twoD else "leveling_1d_reg"), extra, "WM")
    s = fixed_steps(s, dt, 1, dt)
    d = os.path.join(workdir, "evap_threshold_%s_%g_%g" % ("2d" if twoD else "1d", a, c0))
    rc, log = run(exe, d, "a", s, 300)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), dt)
    (_, (f0, nx, ny, _)), (_, (f1, _, _, _)) = fs[0], fs[1]
    i = (ny // 2) * nx + nx // 2 if twoD else nx // 2
    h0, h1 = f0["HH"][i], f1["HH"][i]
    wh, _ = evap_cap_step(a * H0 / (h0 * 1e-6), c0)
    #  h1 to where the Newton iteration stops (NewtonEpsilonRel; the caps act
    #  on its guess); the solute h c to rounding, weighted so that the
    #  one bound of the plan (1e-6) holds it to 1e-12
    return max(abs(h1 / h0 - wh), 1e6 * abs(h1 * f1["CC"][i] / (h0 * f0["CC"][i]) - 1.0)), None


def negpow_profile(r, rmax, n, dr):
    """f_space of SpatialFunc = "negPow" at the radius r, evaluated here from
    its definition: the negative power k (1 - r / rmax')^(-n),
    k = (1 - n)(2 - n) / 2, rmax' = rmax (1 + 1e-3), at most 10, up to
    rmax - eps (eps = 5 cells); from there a quarter cosine that leaves it
    with its value and slope -- the slope a one-cell backward difference --
    and rises to a maximum at rmax; constant beyond."""
    def negpow(x):
        u = 1.0 - x / (rmax * (1.0 + 1e-3))
        if u < 1e-12:
            u = 1e-10
        return min(0.5 * (1.0 - n) * (2.0 - n) * u ** (-n), 10.0)
    eps = 5.0 * dr
    if r < rmax - eps:
        return negpow(r)
    fre = negpow(rmax - eps)
    slope = (fre - negpow(rmax - eps - dr)) / dr
    return fre + 2.0 * eps / math.pi * slope * math.cos(0.5 * math.pi * (min(r, rmax) - rmax) / eps)


def check_negpow_rim(workdir, exe, nx=100, expo=0.25):
    """The "negPow" evaporation profile across its rim: a flat film on
    [-0.5, 0.5] mm with no flow at all (no capillary pressure), the
    rate independent of c and h, and NegPowRmax = 0.3 mm, so that the cells
    run from the negative-power curve through the blend (0.25 .. 0.3 mm) to
    3.9 blend widths past the rim. One Crank-Nicolson step then removes
    E0 f_space(|x|) dt from every cell exactly. Returns (the largest
    difference from negpow_profile, relative to its value at the rim;
    the slope jump where the blend takes over, relative to the slope;
    the spread of f past rmax, relative to its value there)."""
    dt, rmax, L = 1e-3, 0.3e-3, 1.0e-3
    E0 = 0.04 * H0 / dt
    s = flat(case("leveling_1d_reg"), {
        "Grid": {"Nx": str(nx), "Xmin": repr(-L / 2), "Xmax": repr(L / 2)},
        "Physics.LaplacePressure": {"Enable": "false"},
        "Physics.Evaporation": {"Enable": "true", "Coef": repr(E0),
                                "ConcFunc": '"PowerN"', "PowerN": "0.0",
                                "SpatialFunc": '"negPow"', "NegPowExponent": repr(expo),
                                "NegPowRmax": repr(rmax)}}, "WM", nx=nx)
    s = fixed_steps(s, dt, 1, dt)
    d = os.path.join(workdir, "negpow_rim")
    rc, log = run(exe, d, "a", s, 300)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), dt)
    (_, (f0, n, _, _)), (_, (f1, _, _, _)) = fs[0], fs[1]
    dx = L / nx
    xs = [-L / 2 + (i - 0.5) * dx for i in range(1, n - 1)]     # the cells of rows 1 .. n-2
    got = [(f0["HH"][i] - f1["HH"][i]) * 1e-6 / (E0 * dt) for i in range(1, n - 1)]
    want = [negpow_profile(abs(x), rmax, expo, dx) for x in xs]
    if any(abs(x - 1e-3 * f0["xx"][i]) > 1e-3 * dx for i, x in zip(range(1, n - 1), xs)):
        return None, "the cell centers are not where this check puts them"
    scale = negpow_profile(rmax, rmax, expo, dx)
    err = max(abs(g - w) for g, w in zip(got, want)) / scale
    #  x > 0: the last two cells on the negative-power curve and the first two of the
    #  blend. Their slopes differ by the curvature, O(dx / (rmax - r)); a
    #  blend that did not take over the curve's slope would be off by O(1).
    r = [(x, g) for x, g in zip(xs, got) if x > 0]
    k = max(i for i, (x, _) in enumerate(r) if x < rmax - 5 * dx)
    sl = (r[k][1] - r[k - 1][1]) / dx
    sr = (r[k + 2][1] - r[k + 1][1]) / dx
    beyond = [g for x, g in r if x >= rmax]
    return (err, abs(sr - sl) / abs(sl), (max(beyond) - min(beyond)) / scale), None


def inlet_outlet_case(nx, dt, T, twoD=False, corners=False):
    """The film of check_inlet_outlet: a bump of 1.5 um on 5 um, fed through
    a "FreshInlet" at x = 0 and drained through an "Outlet" at x = 1 mm by an
    in-plane gravity of 1e4 m/s^2. twoD: the same, uniform in y (4 cells, y
    periodic). corners: 2-D, nx by nx cells over 1 mm by 1 mm, a round bump
    and the gravity along the diagonal, fed through a FreshInlet at Xmin and
    at Ymin and drained through an Outlet at Xmax and at Ymax, so that every
    kind of open corner (inlet-inlet, inlet-outlet, outlet-outlet) is
    crossed."""
    s = quiet(case("leveling_1d_reg"), periodic=False)
    for sec in ("Initial.SurfaceShape.EllipseBump", "Substrate.EllipseBump"):
        s = drop_section(s, sec)
    s = add_section(s, "Initial.SurfaceShape.EllipseBump", {
        "Order": "1", "CreateType": '"stack"', "Base": "0.5", "Height": "0.15",
        "Center": arr(0.3e-3, 0.3e-3 if corners else 0.0),
        "Width": arr(0.3e-3, 0.3e-3 if corners else 1e30 if twoD else 1.0),
        "SlopeWidth": "0.3"}, "Initial.SurfaceShape")
    g = arr(0.7e4, 0.7e4) if corners else arr(1.0e4, 0.0)
    s = edit_toml(s, {}, {"Grid": {"Nx": str(nx), "Xmin": "0.0", "Xmax": "1.0e-3"},
                          "Boundary.XminBC": {"Type": '"FreshInlet"'},
                          "Boundary.XmaxBC": {"Type": '"Outlet"'},
                          "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": g}})
    if corners:
        s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": str(nx),
                                       "Ymin": "0.0", "Ymax": "1.0e-3"},
                              "Boundary.YminBC": {"Type": '"FreshInlet"'},
                              "Boundary.YmaxBC": {"Type": '"Outlet"'}})
    elif twoD:
        s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": "4",
                                       "Ymin": "-0.1e-3", "Ymax": "0.1e-3"},
                              "Boundary.YminBC": {"Type": '"Periodic"'},
                              "Boundary.YmaxBC": {"Type": '"Periodic"'}})
    return fixed_steps(s, dt, int(round(T / dt)), T / 4)


def max_balance(d):
    """The largest |VolErr| and |ResinVolErr| over the rows of History.dat."""
    return max_volerr(d), max_volerr(d, "ResinVolErr")


def check_solute_balance(workdir, exe, variant="UTOPIA", nx=100, dt=2.5e-4, T=0.4, twoD=False,
                         corners=False):
    """The film of check_inlet_outlet, drying (E = 1e-7 m/s, Raoult) while
    it flows, so that the concentration rises downstream of the fresh inlet
    and the solute diffuses through the edges as well as flowing through
    them. variant: the concentration scheme (UTOPIA, Koren: a limited scheme,
    the step completed with its limiter settled), "adaptive" (UTOPIA, the
    time step cut back and lengthened by the control: only the committed
    steps count), or "none" (no concentration: the volume alone). Returns
    the largest |VolErr| and |ResinVolErr| (0 for "none") of History.dat:
    what crossed the open edges is counted from the face fluxes each step
    put through them, so the two close to rounding. (They were counted from
    the velocities and the face values of the field after each step,
    integrated by the trapezoidal rule: an error of the order of the scheme,
    VolErr -5.7e-7 and ResinVolErr -6.1e-7 here at 100 cells, falling at
    order 1.8 in dx; found, corrected.)"""
    s = edit_toml(inlet_outlet_case(nx, dt, T, twoD, corners), {}, {
        "Physics.Evaporation": {"Enable": "true", "Coef": "1e-7",
                                "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'}})
    if variant == "none":
        s = edit_toml(s, {}, {"Physics.ConcField": {"Model": '"None"'}})
    elif variant == "adaptive":
        s = edit_toml(s, {}, {"Numerics": {"InitialDt": repr(8 * dt)},
                              "Numerics.TimeStepControl": {"Enable": "true", "KeepSteps": "3",
                                                           "ChangeFactor": "0.75", "Criterion": "1e-4",
                                                           "MaxDtFactor": "1e10", "MinDt": "1e-15"}})
    else:
        s = edit_toml(s, {}, {"Numerics.Convection": {"Concentration": '"%s"' % variant}})
    d = os.path.join(workdir, "solute_balance_%s_%s_n%d"
                     % ("corners" if corners else "2d" if twoD else "1d", variant, nx))
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    if variant == "adaptive" and ("dt+" not in log or "dt-" not in log):
        return None, "the time step was not both cut and lengthened: " + log[-200:]
    v, r = max_balance(os.path.join(d, "a"))
    return (v, 0.0 if variant == "none" else r), None


def check_inlet_outlet(workdir, exe, nx=100, dt=2.5e-4, T=0.6):
    """A film on an incline, fed through a "FreshInlet" at x = 0 and drained
    through an "Outlet" at x = 1 mm, carries a bump of 1.5 um out of the
    domain (in-plane gravity 1e4 m/s^2, the bump leaves after about 0.3 s).
    Returns (the largest |VolErr| and |ResinVolErr| of History.dat, which
    count what crossed the edges from the face fluxes of each step; how far
    from uniform the film is at T, relative to the bump; how far the
    concentration moved from its uniform 0.5)."""
    s = inlet_outlet_case(nx, dt, T)
    d = os.path.join(workdir, "inlet_outlet_n%d" % nx)
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    volerr, resinerr = max_balance(os.path.join(d, "a"))
    t, (f, nx_, ny, nz) = fields(os.path.join(d, "a"), T / 4)[-1]
    hh = f["HH"][1:-1]
    flat = max(abs(h - 5.0) for h in hh) / 1.5
    conc = max(abs(c - 0.5) for c in f["CC"])
    return (volerr, resinerr, flat, conc), None


def open_balance(runs):
    """The largest |VolErr| or |ResinVolErr| of check_solute_balance over
    the runs (keyword arguments)."""
    worst = 0.0
    for kw in runs:
        v, err = check_solute_balance(W, X, **kw)
        if v is None:
            raise RuntimeError(err)
        worst = max(worst, max(v))
    return worst


def closed_box(twoD, pinned=None, scheme="Koren", g=300.0, spin=200.0, nozzle=False):
    """A closed box: every edge Symmetry, or one Pinning edge (pinned: its
    name) with the film brought down to it by an EdgeBump. An in-plane
    gravity of g m/s^2 toward Xmin (and Ymin in 2-D), the vertical gravity,
    in 2-D a rotation about a point off the center; the capillary pressure.
    nozzle: a nozzle at the center instead of the body forces."""
    s = quiet(case("leveling_2d_reg" if twoD else "leveling_1d_reg"), twoD, periodic=False)
    if pinned:
        for sec in ("Initial.SurfaceShape.EllipseBump", "Substrate.EllipseBump"):
            s = drop_section(s, sec)
        L = 1.0e-3
        s = add_section(s, "Initial.SurfaceShape.EdgeBump", {
            "Order": "1", "CreateType": '"add"', "Base": "0.5", "Height": "1.0",
            "Location": '"%s"' % pinned, "Alpha": "1.0", "InnerLength": repr(0.05 * L),
            "OuterLength": repr(0.15 * L)}, "Initial.SurfaceShape")
    else:
        #  A bump against the walls at Xmin and Ymin: on a film the same at
        #  two opposite walls what a leak takes out at one it puts in at the
        #  other
        s = drop_section(s, "Initial.SurfaceShape.EllipseBump")
        s = add_section(s, "Initial.SurfaceShape.EllipseBump", {
            "Order": "1", "CreateType": '"add"', "Base": "0.0", "Height": "0.7",
            "Center": arr(-0.35e-3, -0.3e-3 if twoD else 0.0),
            "Width": arr(0.7e-3, 0.6e-3 if twoD else 1.0),
            "SlopeWidth": "0.1e-3"}, "Initial.SurfaceShape")
    #  A precursor, so that a Symmetry box has a film everywhere; a thin one
    #  for a pinned edge, where the film has to come down to zero
    s = edit_toml(s, {}, {"Initial.SurfaceShape": {"PrecursorRatio": "0.01" if pinned else "0.3"}})
    edges = {e: '"Symmetry"' for e in ("Xmin", "Xmax", "Ymin", "Ymax")}
    if pinned:
        edges[pinned] = '"Pinning"'
    kv = {"Grid": {"Nx": "24" if twoD else "48"},
          "Physics.VerticalGravity": {"Enable": "true"},
          "Numerics.Convection": {"Height": '"%s"' % scheme, "Concentration": '"%s"' % scheme}}
    if twoD:
        kv["Grid"]["Ny"] = "20"
    for e, t in edges.items():
        kv["Boundary.%sBC" % e] = {"Type": t}
    if nozzle:
        kv["Physics.Source"] = {"Enable": "true", "Radius": "0.1e-3", "InletRate": "1e-5",
                                "TimeConstant": "0.0", "ArmLength": "0.0",
                                "ArmCenter": arr(0.1e-3, 0.0), "RotationRate": "0.0",
                                "StartAngle": "0.0", "RampTime": "0.0", "PlateauTime": "1e10",
                                "SubRegions": "4"}
    else:
        kv["Physics.HorizontalGravity"] = {"Enable": "true",
                                           "Acceleration": arr(-g, -0.6 * g if twoD else 0.0)}
        #  In 1-D too: the liquid it flings sideways, into the y that is
        #  not solved for, is counted in OutflowedVol (and its solute in
        #  OutflowedResin), so the box balances to rounding all the same
        if spin:
            kv["Physics.Centrifugal"] = {"Enable": "true", "RotationRate": repr(spin),
                                         "RotationAxis": arr(-0.2e-3, 0.1e-3)}
    return edit_toml(s, {}, kv)


def check_nozzle_area(workdir, exe, dt=1e-3, steps=5):
    """A fixed 2-D nozzle of radius R = 0.12 mm in a closed box, its area
    found from 20 x 20 points a cell, on cells of 25 x 31 um, of 17 x 188 um
    (11 times longer in y) and of 250 x 12.5 um: what it supplied over the
    run must be InletRate pi R^2 t, whatever the shape of the cells. Returns
    the largest |supplied / (InletRate pi R^2 t) - 1|. (The cells the nozzle
    covers were found by a test on the distance of their centers from it,
    R + 2 min(dx, dy), which cells longer than 3.9 times their width
    overreach: 19 % of the nozzle was left out on the second grid.)"""
    R, rate = 0.12e-3, 1e-5
    worst = 0.0
    for nx, ny in ((40, 24), (60, 4), (4, 60)):
        s = closed_box(True, nozzle=True)
        s = edit_toml(s, {}, {"Grid": {"Nx": str(nx), "Ny": str(ny), "Xmin": "-0.5e-3", "Xmax": "0.5e-3",
                                       "Ymin": "-0.3e-3", "Ymax": "0.45e-3"},
                              "Physics.Source": {"Radius": repr(R), "InletRate": repr(rate),
                                                 "ArmCenter": arr(0.1e-3, -0.05e-3), "SubRegions": "20"}})
        s = fixed_steps(s, dt, steps, dt * steps)
        d = os.path.join(workdir, "nozzle_area_%dx%d" % (nx, ny))
        rc, log = run(exe, d, "a", s, 900)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        hist = [l.split("\t") for l in open(os.path.join(d, "a", "History.dat"))]
        cols = [c.split("[")[0].strip() for c in hist[0]]
        t, sup = float(hist[-1][cols.index("Time")]), float(hist[-1][cols.index("SuppliedVol")])
        worst = max(worst, abs(sup / (rate * math.pi * R * R * t) - 1.0))
    return worst, None


def max_volerr(d, col="VolErr"):
    hist = [l.split("\t") for l in open(os.path.join(d, "History.dat"))]
    k = [c.split("[")[0].strip() for c in hist[0]].index(col)
    return max(abs(float(r[k])) for r in hist[1:])


def cz_box(s):
    """A case made ConcZ: 6 layers, and no convection scheme for the film
    height (the film-height equation of ConcZ takes none)."""
    s = edit_toml(s, {}, {"Physics.ConcField": {"Model": '"CZ"'}, "Grid": {"Nz": "6"}})
    out, cur = [], None
    for l in s.split("\n"):
        st = l.strip()
        if st.startswith("[") and st.endswith("]"):
            cur = st[1:-1]
        if cur == "Numerics.Convection" and st.startswith("Height"):
            continue
        out.append(l)
    return "\n".join(out)


def check_spin_balance(workdir, exe, dt=2e-3, steps=50):
    """A 1-D closed box spinning about an axis off its center (closed_box,
    every edge Symmetry, or pinned at Xmin): the liquid the centrifugal term
    flings sideways, into the y that is not solved for, is counted in
    OutflowedVol and its solute in OutflowedResin, so |VolErr| and
    |ResinVolErr| are at rounding -- well-mixed with a fixed and with an
    adaptive dt, ConcZ without and with evaporation (which stratifies c).
    Without the count VolErr is the volume flung off: -0.12 of the volume
    over the 0.1 s of the well-mixed box at 3000 rad/s."""
    worst = 0.0
    runs = []
    for pinned in (None, "Xmin"):
        runs.append(("wm_%s" % (pinned or "sym"), fixed_steps(closed_box(False, pinned, spin=3000.0),
                                                              dt, steps, dt * steps / 2), True))
    a = fixed_steps(closed_box(False, spin=3000.0), dt, steps, dt * steps / 2)
    a = edit_toml(a, {}, {"Numerics.TimeStepControl": {"Enable": "true", "KeepSteps": "3",
                                                        "ChangeFactor": "0.75", "Criterion": "1e-2",
                                                        "MaxDtFactor": "1e10", "MinDt": "1e-15"}})
    runs.append(("wm_adapt", a, True))
    for evap in (False, True):
        c = cz_box(closed_box(False, spin=3000.0))
        c = edit_toml(c, {}, {"Physics.Evaporation": {
            "Enable": "true" if evap else "false", "Coef": "2e-5",
            "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'}})
        runs.append(("cz_%s" % ("evap" if evap else "closed"), fixed_steps(c, 1e-3, 40, 2e-2), True))
    for tag, src, solute in runs:
        d = os.path.join(workdir, "spin_balance_%s" % tag)
        rc, log = run(exe, d, "a", src, 900)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        if max_volerr(os.path.join(d, "a"), "OutflowedVol") <= 0:
            return None, "%s: nothing was flung off" % tag
        worst = max(worst, max_volerr(os.path.join(d, "a")))
        if solute:
            worst = max(worst, max_volerr(os.path.join(d, "a"), "ResinVolErr"))
    return worst, None


def check_wall_balance(workdir, exe, twoD, pinned=None, schemes=("Koren", "KK", "UTOPIA", "QUICK"),
                       dt=2e-3, steps=50, slip=None):
    """closed_box() for each scheme, a fixed dt; the largest |VolErr|.
    slip: a Navier slip length [m] (Physics.Slip)."""
    worst = 0.0
    for sc in schemes:
        s = fixed_steps(closed_box(twoD, pinned, sc), dt, steps, dt * steps / 2)
        tag = ""
        if slip:
            s = edit_toml(s, {}, {"Physics.Slip": {"Enable": "true", "Length": repr(slip)}})
            tag = "_slip"
        d = os.path.join(workdir, "wall_%s_%s_%s%s" % ("2d" if twoD else "1d", pinned or "sym", sc, tag))
        rc, log = run(exe, d, "a", s, 900)
        if rc != 0:
            return None, "run failed: " + log[-400:]
        worst = max(worst, max_volerr(os.path.join(d, "a")))
    return worst, None


def cz_sums(f):
    """The volume sum h and the solute sum h mean(c) over the cell columns of
    a ConcZ field (the wall faces left out)."""
    nx, ny, nz = f[1], f[2], f[3]
    f = f[0]
    n2 = nx * ny
    vol = sol = 0.0
    for j in (range(1, ny - 1) if ny > 1 else [0]):
        for i in range(1, nx - 1):
            r = j * nx + i
            cs = [f["CC"][k * n2 + r] for k in range(1, nz - 1)]
            vol += f["HH"][r]
            sol += f["HH"][r] * sum(cs) / len(cs)
    return vol, sol


def check_wall_balance_cz(workdir, exe, twoD, evap, scheme="Koren", dt=1e-3, steps=40):
    """A closed ConcZ box (every edge Symmetry), every term of the model on:
    the capillary pressure over a relief, the Marangoni stress, the vertical
    gravity, an in-plane gravity toward the walls at Xmin (and Ymin), in 2-D
    a rotation about a point off the center, a viscosity and a diffusivity
    that depend on c; with evaporation (which makes c, and with it the
    Marangoni stress, vary) or without. Nothing may cross a wall: returns
    (the change of the volume, without evaporation, else 0; the change of
    the solute), from the first and the last field file, which carry every
    digit. (ConcZ took the in-plane body forces through a Symmetry wall from
    the ghosts: tools/proptest.py seed 1015.)"""
    s = closed_box(twoD, scheme=scheme)
    for sec in ("Substrate.EllipseBump",):
        s = drop_section(s, sec)
    s = add_section(s, "Substrate.EllipseBump", {
        "CreateType": '"add"', "Base": "0.0", "Height": "0.3",
        "Center": arr(0.1e-3, 0.05e-3 if twoD else 0.0),
        "Width": arr(0.5e-3, 0.4e-3 if twoD else 1.0), "SlopeWidth": "0.5"}, "Substrate")
    kv = {"Physics.ConcField": {"Model": '"CZ"'},
          "Grid": {"Nz": "6"},
          "Physics.SolutoCapillarity": {"Enable": "true"},
          "Fluid.SurfaceTension": {"Coef": "10e-3", "BaseConc": "0.5"},
          "Fluid.Viscosity": {"Exponent": "3.0", "BaseConc": "0.5"},
          "Fluid.Diffusivity": {"BaseValue": "1e-10", "Exponent": "-2.0", "BaseConc": "0.5"},
          "Physics.Evaporation": {"Enable": "true" if evap else "false", "Coef": "2e-5",
                                  "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'}}
    if not twoD:
        #  It flings the 1-D film into y, which the field files do not see
        #  (check_spin_balance counts it, from History.dat)
        kv["Physics.Centrifugal"] = {"Enable": "false"}
    s = edit_toml(s, {}, kv)
    #  The film-height equation of ConcZ takes no convection scheme
    out, cur = [], None
    for l in s.split("\n"):
        st = l.strip()
        if st.startswith("[") and st.endswith("]"):
            cur = st[1:-1]
        if cur == "Numerics.Convection" and st.startswith("Height"):
            continue
        out.append(l)
    s = "\n".join(out)
    s = fixed_steps(s, dt, steps, dt * steps)
    d = os.path.join(workdir, "wall_cz_%s_%s_%s" % ("2d" if twoD else "1d", "evap" if evap else "closed", scheme))
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), dt * steps)
    (v0, s0), (v1, s1) = cz_sums(fs[0][1]), cz_sums(fs[-1][1])
    return (0.0 if evap else abs(v1 / v0 - 1.0), abs(s1 / s0 - 1.0)), None


def overshoot_film(twoD, nx, dt, steps, scheme, evap=0.0):
    """A flat well-mixed film of h0 = 5 um, periodic, carried by in-plane
    gravity at the uniform velocity U = 1 cm/s along x (in 2-D at
    (U, 0.6 U)) over L = 1 mm; nothing but the flow moves the solute (D =
    1e-13 m^2/s, no Marangoni), so the solute is kept to rounding by the
    advection itself. evap > 0: the film dries by Raoult's law at that
    coefficient (m/s), which takes solvent only."""
    L, U, mu, rho = 1.0e-3, 1.0e-2, 0.1, 1000.0
    g = 3.0 * mu * U / (rho * H0 ** 2)
    s = quiet(case("leveling_1d_reg"), twoD)
    s = wave_film(s, 0.5, 0.0, L, L if twoD else None)
    kv = {"Grid": {"Nx": str(nx), "Xmin": repr(-L / 2), "Xmax": repr(L / 2)},
          "Physics.ConcField": {"Model": '"WM"', "EnableConvection": "true", "EnableDiffusion": "true"},
          "Physics.LaplacePressure": {"Enable": "true"},
          "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(g, 0.6 * g if twoD else 0.0)},
          "Fluid.Density": {"Value": repr(rho)},
          "Fluid.Viscosity": {"BaseValue": repr(mu), "Exponent": "0.0"},
          "Fluid.Diffusivity": {"BaseValue": "1.0e-13", "Exponent": "0.0"},
          "Numerics.Convection": {"Concentration": '"%s"' % scheme}}
    if evap > 0.0:
        kv["Physics.Evaporation"] = {"Enable": "true", "Coef": repr(evap),
                                     "ConcFunc": '"Raoult"', "SpatialFunc": '"uniform"'}
    if twoD:
        kv["Grid"].update({"Dimension": '"xyh"', "Ny": str(nx), "Ymin": repr(-L / 2), "Ymax": repr(L / 2)})
    s = edit_toml(s, {}, kv)
    return fixed_steps(s, dt, steps, dt * steps), L, U


def check_conc_overshoot(workdir, exe, scheme, courant, twoD=False, nx=50, passage=0.3,
                         evap=0.0, front=(0.99, 0.3)):
    """overshoot_film() from a sharp front, c = 0.99 on half of the domain
    and 0.3 on the other (in 2-D along a diagonal, so that both directions
    carry it), for a fraction `passage` of the time the flow takes to cross
    the domain, at the Courant number u dt/dx = courant. An unbounded scheme
    overshoots past c = 1 at the front; the excess must be moved to the cells
    around, not cut. Returns (the largest |ResinVolErr| of History.dat, the
    solute moved from above c = 1 as JobInfo.txt reports it, the largest c
    of the run - 1). The cap that cut it lost 0.7 % of the solute over one
    crossing with UTOPIA on 100 cells (found; corrected), and the limited
    schemes past Courant 1 a further 1 % through a limiter frozen at the
    capped field (found; corrected).

    evap > 0 dries the film as it goes (overshoot_film); front gives the
    two concentrations. The fourth value returned is the solute moved into
    cells below c = 0, the fifth the smallest c of the run."""
    tag = "overshoot_%s_%s_C%g_n%d_e%g_f%g_%g" % ("2d" if twoD else "1d", scheme, courant, nx,
                                                 evap, front[0], front[1])
    d = os.path.join(workdir, tag)
    src, L, U = overshoot_film(twoD, nx, 1.0, 1, scheme, evap)
    dt = courant * (L / nx) / U
    steps = max(1, int(round(passage * L / U / dt)))
    #  The film, from a run of one step, its concentration replaced
    rc, log = run(exe, os.path.join(d, "init"), "a", overshoot_film(twoD, nx, dt, 1, scheme, evap)[0], 300)
    if rc != 0:
        return None, "initial run failed: " + log[-400:]
    with open(os.path.join(d, "init", "a", "field00000.dat")) as f:
        lines = f.read().split("\n")
    names = [c.split("[")[0] for c in lines[0].lstrip("#").split()]
    ix, ic = names.index("xx"), names.index("CC")
    iy = names.index("yy") if "yy" in names else None
    out = [lines[0]]
    for ln in lines[1:]:
        v = ln.split()
        if len(v) == len(names):
            x = float(v[ix]) * 1e-3
            y = float(v[iy]) * 1e-3 if iy is not None else 0.0
            v[ic] = "%.17e" % (front[0] if math.cos(2 * math.pi * (x + y) / L) > 0 else front[1])
            ln = "\t".join(v)
        out.append(ln)
    with open(os.path.join(d, "start.dat"), "w") as f:
        f.write("\n".join(out))
    s = overshoot_film(twoD, nx, dt, steps, scheme, evap)[0]
    s = edit_toml(s, {}, {"Restart": {"Enable": "true", "ASCII": "true", "FileName": '"../start.dat"',
                                      "RestartTime": "0", "FileIndex": "0"}})
    rc, log = run(exe, os.path.join(d, "run"), "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    r = os.path.join(d, "run", "a")
    #  (ConcZ moves solute above c = 1 within the column and does not
    #  report it)
    moved = jobinfo_value(os.path.join(r, "JobInfo.txt"), "Solute moved from c > 1 [-]", 0.0)
    under = jobinfo_value(os.path.join(r, "JobInfo.txt"), "Solute moved to c < 0 [-]")
    return (max_volerr(r, "ResinVolErr"), moved, max_volerr(r, "maxC") - 1.0, under,
            min_column(r, "minC")), None


def min_column(d, col):
    hist = [l.split("\t") for l in open(os.path.join(d, "History.dat"))]
    k = [c.split("[")[0].strip() for c in hist[0]].index(col)
    return min(float(r[k]) for r in hist[1:])


def check_conc_uniform(workdir, exe, twoD, c0=0.37, dt=1e-3, steps=30):
    """A uniform c0 in closed_box() with a nozzle, the in-plane gravity
    (300 m/s^2 toward Xmin, and Ymin in 2-D) and, in 1-D, the spin as well.
    Returns the largest |c - c0| over the real cells of the fields written."""
    s = closed_box(twoD, nozzle=True)
    kv = {"Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(-300.0, -180.0 if twoD else 0.0)},
          "Physics.ConcField": {"Model": '"WM"', "EnableConvection": "true", "EnableDiffusion": "true"},
          "Numerics.Convection": {"Concentration": '"UTOPIA"'}}
    if not twoD:
        kv["Physics.Centrifugal"] = {"Enable": "true", "RotationRate": "200.0",
                                     "RotationAxis": arr(-0.2e-3, 0.1e-3)}
    s = edit_toml(s, {}, kv)
    s = set_initial_conc(s, c0)
    interval = dt * steps / 3
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, "conc_uniform_%s" % ("2d" if twoD else "1d"))
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    worst = 0.0
    fs = fields(os.path.join(d, "a"), interval)
    if len(fs) < 3:
        return None, "too few fields written"
    for t, (f, nx, ny, nz) in fs:
        worst = max(worst, max(abs(c - c0) for c in f["CC"]))
    return worst, None


def set_initial_conc(src, c0):
    """The initial concentration made uniform at c0."""
    return edit_toml(src, {}, {"Initial.Concentration": {"Type": '"uniform"', "Value": repr(c0)}})


def conc_undershoot(runs, cz=False):
    """check_conc_overshoot() from the front 0.7 | 0.005, which the unbounded
    schemes carry below c = 0: the largest |ResinVolErr|, or 1 (a failure)
    where nothing was moved into cells below 0 or c was left below 0."""
    worst = 0.0
    for sc, co, t in runs:
        if cz:
            v, err = check_conc_overshoot_cz(W, X, sc, co, front=(0.7, 0.005))
        else:
            v, err = check_conc_overshoot(W, X, sc, co, t, nx=24 if t else 100, front=(0.7, 0.005))
        if v is None:
            raise RuntimeError(err)
        if v[3] < 1e-4 or v[4] < 0.0:
            return 1.0
        worst = max(worst, v[0])
    return worst


def check_conc_overshoot_cz(workdir, exe, scheme, courant, front):
    """check_conc_overshoot() with the film made ConcZ (cz_box), 1-D."""
    global overshoot_film
    plain = overshoot_film
    def cz(*a, **k):
        src, L, U = plain(*a, **k)
        return cz_box(src), L, U
    overshoot_film = cz
    try:
        return check_conc_overshoot(os.path.join(workdir, "cz"), exe, scheme, courant, False, nx=100, front=front)
    finally:
        overshoot_film = plain


def conc_drying(runs):
    """check_conc_overshoot() from the front 0.9 | 0.3 on a film drying by
    Raoult's law (2e-5 m/s): the largest |ResinVolErr| of the runs."""
    worst = 0.0
    for sc, t in runs:
        v, err = check_conc_overshoot(W, X, sc, 0.4, t, nx=24 if t else 50, evap=2e-5, front=(0.9, 0.3))
        if v is None:
            raise RuntimeError(err)
        worst = max(worst, v[0])
    return worst


def conc_overshoot(runs):
    """The largest |ResinVolErr| of the runs (scheme, Courant number, 2-D),
    or 1 (a failure) where nothing was moved (the overshoot not provoked) or
    c was left above 1."""
    worst = 0.0
    for sc, co, t in runs:
        v, err = check_conc_overshoot(W, X, sc, co, t)
        if v is None:
            raise RuntimeError(err)
        if v[1] < 1e-4 or v[2] > 0.0:
            return 1.0
        worst = max(worst, v[0])
    return worst


def check_nozzle_balance(workdir, exe, twoD=False, dt=1e-3, steps=60):
    """A nozzle in a closed box, the time step adapted (it is lengthened
    after accepted steps, which is where DT().now is no longer the step's
    dt): the largest |VolErr| and |ResinVolErr|."""
    s = fixed_steps(closed_box(twoD, nozzle=True), dt, steps, dt * steps / 2)
    s = edit_toml(s, {}, {"Numerics.TimeStepControl": {"Enable": "true", "KeepSteps": "3",
                                                        "ChangeFactor": "0.75", "Criterion": "1e-2",
                                                        "MaxDtFactor": "1e10", "MinDt": "1e-15"}})
    d = os.path.join(workdir, "nozzle_balance_%s" % ("2d" if twoD else "1d"))
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    if "dt+" not in log:
        return None, "the time step was never lengthened: " + log[-200:]
    res = os.path.join(d, "a")
    return (max_volerr(res), max_volerr(res, "ResinVolErr")), None


def check_nozzle_thick(workdir, exe):
    """nozzle_1d_rotating with its arm turning (1 rad/s, so the nozzle
    moves half the domain and partly out of it): the film thickens from
    10 um to 1.1 mm over cells 3.3 um wide, the substrate spinning (the
    centrifugal term flings liquid into the y that is not solved for), the
    time step adapted. The largest |VolErr| and |ResinVolErr|. With A h
    formed as the product of the entries of A with h, the capillary term
    added eps |A| |h| per cell and step: VolErr 4e-5 at the end (of the
    initial volume; 2e-7 of the film); formed face by face it is at rounding."""
    s = case("nozzle_1d_rotating")
    s = edit_toml(s, {}, {"Physics.Source": {"RotationRate": "1.0"},
                          "Output.FieldData": {"TimeInterval": "1.0", "VTKFile": "false"}})
    d = os.path.join(workdir, "nozzle_thick")
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    res = os.path.join(d, "a")
    return (max_volerr(res), max_volerr(res, "ResinVolErr")), None


def check_nozzle_supplied(workdir, exe, twoD=False, cs=0.9):
    """check_nozzle_balance() (twoD = "thick": check_nozzle_thick()) with
    the liquid supplied at Physics.Source.Conc = cs (ConcType = "supplied")
    instead of at the concentration it lands on (the film starts at 0.5):
    the largest |VolErr| and |ResinVolErr|, SuppliedResin being the sum of
    dt S cs the solute balance took in, and how far c moved from 0.5 (the
    largest |c - 0.5| over the History rows). The plan checks |VolErr| of
    the closed boxes and the thick film alike. At cs = 0.1 the viscosity
    under the turning nozzle of the thick film falls e^4-fold (Exponent 10)
    and the Newton correction of some steps stagnates far from convergence
    (8e-4 to 1.9 times the first): such a step fails and is cut
    (C_NewtonKantorovich::NKLoop). The run also returns how many steps
    failed so (JobInfo.txt), which must be one or more on the thick film:
    found taken with a warning, the film keeping the residual as VolErr
    2.2e-8; corrected, VolErr 2e-13."""
    if twoD == "thick":
        s = case("nozzle_1d_rotating")
        s = edit_toml(s, {}, {"Physics.Source": {"RotationRate": "1.0", "ConcType": '"supplied"',
                                                 "Conc": repr(cs)},
                              "Output.FieldData": {"TimeInterval": "1.0", "VTKFile": "false"}})
        d = os.path.join(workdir, "nozzle_supplied_thick")
    else:
        s = fixed_steps(closed_box(twoD, nozzle=True), 1e-3, 60, 0.03)
        s = edit_toml(s, {}, {"Numerics.TimeStepControl": {"Enable": "true", "KeepSteps": "3",
                                                            "ChangeFactor": "0.75", "Criterion": "1e-2",
                                                            "MaxDtFactor": "1e10", "MinDt": "1e-15"},
                              "Physics.Source": {"ConcType": '"supplied"', "Conc": repr(cs)}})
        d = os.path.join(workdir, "nozzle_supplied_%s" % ("2d" if twoD else "1d"))
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    res = os.path.join(d, "a")
    hist = [l.split("\t") for l in open(os.path.join(res, "History.dat"))]
    cols = [c.split("[")[0].strip() for c in hist[0]]
    moved = max(abs(float(r[cols.index(k)]) - 0.5) for r in hist[1:] for k in ("minC", "maxC"))
    cuts = 0
    for l in open(os.path.join(res, "JobInfo.txt")):
        if l.startswith("Steps failed: the Newton correction stagnated :"):
            cuts = int(l.rsplit(":", 1)[1])
    return (max_volerr(res), max_volerr(res, "ResinVolErr"), moved, cuts), None


_SUPPLIED = {}


def nozzle_supplied(twoD):
    """check_nozzle_supplied, run once per film; a failure (1) where the
    supplied concentration did not move c (the mode not in effect), or
    where no step of the thick film failed on a stagnated correction (the
    path this film is there to exercise)."""
    if twoD not in _SUPPLIED:
        v, why = check_nozzle_supplied(W, X, twoD, cs=0.1 if twoD == "thick" else 0.9)
        if why:
            raise RuntimeError(why)
        _SUPPLIED[twoD] = v
    v = _SUPPLIED[twoD]
    if v[2] <= 1e-3 or (twoD == "thick" and v[3] < 1):
        return (1.0, 1.0, v[2], v[3])
    return v


def check_nozzle_mixing(workdir, exe, twoD, footprint, conc_type="supplied", h0=5e-6, c0=0.5,
                        cs=0.1, rate=1e-5, dt=1e-3, steps=300):
    """Liquid supplied onto a flat film that does not flow (no capillary
    pressure or any other term; no evaporation), at Physics.Source.Conc =
    cs (ConcType = "supplied", or "local"). footprint "all": a nozzle far
    wider than the domain, so S = InletRate in every cell and
        h(t) = h0 + S t,   c(t) = (h0 c0 + S t cs) / (h0 + S t);
    "part": a nozzle of radius 0.1 mm over 40 cells (2-D: 20 x 16), some
    covered in part, diffusion off, where every cell keeps its own balance,
        c = (h0 c0 + (h - h0) cs) / h
    with its own h, and h = h0 outside the footprint. Each step adds dt S to
    h and dt S cs to h c exactly (S is constant; no time-discretization
    error), so the error is rounding, to the digits of the field files.
    Returns the largest |c - exact| (and |h - exact| / h for "all") over
    the cells of the fields written, and the largest |c - c0| (how far the
    supply moved c)."""
    L = 1.0e-3
    s = flat(case("leveling_2d_reg" if twoD else "leveling_1d_reg"), {
        "Physics.LaplacePressure": {"Enable": "false"},
        "Physics.ConcField": {"Model": '"WM"', "EnableConvection": "true",
                              "EnableDiffusion": "true" if footprint == "all" else "false"},
        "Physics.Source": {"Enable": "true", "Radius": "20e-3" if footprint == "all" else "0.1e-3",
                           "InletRate": repr(rate), "ConcType": '"%s"' % conc_type, "Conc": repr(cs),
                           "TimeConstant": "0.0", "ArmLength": "0.0",
                           "ArmCenter": arr(0.05e-3, 0.03e-3 if twoD else 0.0), "RotationRate": "0.0",
                           "StartAngle": "0.0", "RampTime": "0.0", "PlateauTime": "1e10",
                           "SubRegions": "4"}}, nx=20 if twoD else 40)
    kv = {"Grid": {"Xmin": repr(-L / 2), "Xmax": repr(L / 2)}}
    if twoD:
        kv["Grid"].update({"Ny": "16", "Ymin": repr(-L / 2), "Ymax": repr(L / 2)})
    s = edit_toml(s, {}, kv)
    s = set_initial_conc(s, c0)
    interval = dt * steps / 3
    s = fixed_steps(s, dt, steps, interval)
    d = os.path.join(workdir, "nozzle_mixing_%s_%s_%s" % ("2d" if twoD else "1d", footprint, conc_type))
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    fs = fields(os.path.join(d, "a"), interval)
    if len(fs) < 4:
        return None, "too few fields written"
    want_cs = cs if conc_type == "supplied" else c0
    werr = moved = 0.0
    for t, (f, nx, ny, nz) in fs:
        for k in range(len(f["HH"])):
            i, j = k % nx, k // nx
            if i in (0, nx - 1) or (ny > 1 and j in (0, ny - 1)):
                continue
            h, c = f["HH"][k] * 1e-6, f["CC"][k]
            if footprint == "all":
                hx = h0 + rate * t
                werr = max(werr, abs(h - hx) / hx)
            else:
                hx = h
            cx = (h0 * c0 + (hx - h0) * want_cs) / hx
            werr = max(werr, abs(c - cx))
            moved = max(moved, abs(c - c0))
    return (werr, moved), None


def nozzle_mixing():
    """check_nozzle_mixing over 1-D and 2-D, both footprints, "supplied";
    and "local" (c stays c0). The largest error, or 1 (a failure) where the
    supplied liquid did not move c, or the local one did."""
    worst = 0.0
    for twoD in (False, True):
        for fp in ("all", "part"):
            v, why = check_nozzle_mixing(W, X, twoD, fp)
            if why:
                raise RuntimeError(why)
            if v[1] < 0.1:
                return 1.0
            worst = max(worst, v[0])
        v, why = check_nozzle_mixing(W, X, twoD, "part", conc_type="local")
        if why:
            raise RuntimeError(why)
        if v[1] > 0.0:
            return 1.0
        worst = max(worst, v[0])
    return worst


_NOZZLE = {}


def nozzle_balance(twoD):
    """check_nozzle_balance, run once per dimension for the two entries
    (twoD = "thick": check_nozzle_thick)."""
    if twoD not in _NOZZLE:
        v, why = check_nozzle_thick(W, X) if twoD == "thick" else check_nozzle_balance(W, X, twoD)
        if why:
            raise RuntimeError(why)
        _NOZZLE[twoD] = v
    return _NOZZLE[twoD]


def jobinfo_value(path, title, default=None):
    for line in open(path):
        if line.startswith(title + " :"):
            return float(line.split(":", 1)[1])
    if default is not None:
        return default
    raise ValueError("%s: no '%s'" % (path, title))


def check_leveling_adaptive(workdir, exe, model="WM", nx=32, criterion="1e-5"):
    """The sinusoid of leveling_1d over three e-folds, the time step adapted:
    started at a tenth of the output interval (too long: cut back), then
    lengthened every third accepted step. Returns (the largest relative time
    shift |t + ln(a/a0)/rate| / t over the outputs, |time run - sum of the
    steps solved| / time run). The rate is that of the discrete operator,
    (2 sin(k dx/2)/dx)^4 for k^4, so what is left is the error of the time
    integration."""
    L = 0.25e-3
    s = quiet(case("concz_1d_reg" if model == "CZ" else "leveling_1d_reg"))
    s = edit_toml(s, {}, {"Grid": {"Nx": str(nx), "Xmin": repr(-L / 2), "Xmax": repr(L / 2)}})
    s = wave_film(s, 0.5, 1e-3, L)
    dx = L / nx
    rate = SIG * H0 ** 3 * (2 * math.sin(math.pi / nx) / dx) ** 4 / (3 * MU)
    T = 3.0 / rate
    interval = T / 6
    s = fixed_steps(s, interval / 10, 10 ** 9, interval)
    s = edit_toml(s, {}, {"Numerics": {"EndTime": repr(T)},
                          "Numerics.TimeStepControl": {
                              "Enable": "true", "KeepSteps": "3", "ChangeFactor": "0.75",
                              "Criterion": criterion, "MaxDtFactor": "1e10", "MinDt": "1e-15"}})
    d = os.path.join(workdir, "leveling_adaptive_%s" % model)
    rc, log = run(exe, d, "a", s, 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    if "dt+" not in log or "dt-" not in log:
        return None, "the time step was not both cut back and lengthened"
    res = os.path.join(d, "a")
    ji = os.path.join(res, "JobInfo.txt")
    try:
        run_time = (jobinfo_value(ji, "Time at the end of the loop")
                    - jobinfo_value(ji, "Time at the start of the loop"))
        clock = abs(jobinfo_value(ji, "Sum of the steps solved") - run_time) / run_time
    except ValueError:
        clock = float("inf")        # a build that does not report the sum
    #  The field files and the History.dat rows are written at the same steps
    #  (the same interval); the time is taken from History.dat
    times = []
    for line in open(os.path.join(res, "History.dat")):
        if line.strip() and not line.startswith(("#", "Time")):
            times.append(float(line.split()[0]))
    fs = sorted(glob.glob(os.path.join(res, "field*.dat")))
    if len(fs) != len(times) or len(fs) < 6:
        return None, "%d field files and %d History rows" % (len(fs), len(times))
    kx = 2 * math.pi / L
    f0, nx0, ny0, _ = read_field(fs[0])
    a0 = abs(mode(f0, nx0, ny0, kx, 0.0))
    shift = 0.0
    for p, t in zip(fs[1:], times[1:]):
        f, nx_, ny_, _ = read_field(p)
        shift = max(shift, abs(t + math.log(abs(mode(f, nx_, ny_, kx, 0.0)) / a0) / rate) / t)
    return (shift, clock), None


_ADAPTIVE = {}


def leveling_adaptive(model="WM"):
    """check_leveling_adaptive, run once per model for the two entries."""
    if model not in _ADAPTIVE:
        v, why = check_leveling_adaptive(W, X, model)
        if why:
            raise RuntimeError(why)
        _ADAPTIVE[model] = v
    return _ADAPTIVE[model]


def state_time(path):
    """The time of a field*.state, dimensionless, to every digit."""
    with open(path, "rb") as f:
        m = re.search(rb"\ntime ([0-9.eE+-]+)", f.read())
    return float(m.group(1))


def check_leveling_logout(workdir, exe):
    """The shipped case leveling_1d_logout: a sinusoid of 0.375 mm on a
    periodic 5 um film, written at 141 times from 1e-7 to 1 s, 20 a decade,
    each landed on (ExactTimes). Returns the largest relative difference of
    an output time from its listed time, the largest relative error of the
    decay rate measured at any output against sigma h0^3 k^4 / 3 mu, and the
    steps taken."""
    d = os.path.join(workdir, "leveling_logout")
    rc, log = run(exe, d, "leveling_1d_logout", case("leveling_1d_logout"), 900)
    if rc != 0:
        return None, "run failed: " + log[-400:]
    job = os.path.join(d, "leveling_1d_logout")
    L0 = math.sqrt(SIG / (RHO * G))
    T0 = MU * L0 ** 4 / (SIG * 10e-6 ** 3)
    want = [1e-7 * 10 ** (i / 20.0) for i in range(141)]
    paths = sorted(glob.glob(os.path.join(job, "field*.dat")))
    if len(paths) != 143:
        return None, "%d field files, not 143 (the start, 141 listed, the end)" % len(paths)
    k = 2 * math.pi / 0.375e-3
    rate = SIG * H0 ** 3 * k ** 4 / (3 * MU)
    f0, nx, ny, _ = read_field(paths[0])
    a0 = abs(mode(f0, nx, ny, k, 0.0))
    terr, rerr = 0.0, 0.0
    for i, p in enumerate(paths[1:142]):
        t = state_time(p[:-4] + ".state") * T0
        terr = max(terr, abs(t - want[i]) / want[i])
        f, nx, ny, _ = read_field(p)
        got = -math.log(abs(mode(f, nx, ny, k, 0.0)) / a0) / t
        rerr = max(rerr, abs(got - rate) / rate)
    steps = int(jobinfo_value(os.path.join(job, "JobInfo.txt"), "Time steps", 0))
    return (terr, rerr, steps), None


_LOGOUT = {}


def leveling_logout():
    if not _LOGOUT:
        v, why = check_leveling_logout(W, X)
        if why:
            raise RuntimeError(why)
        _LOGOUT["v"] = v
        print("     leveling_1d_logout: %d steps for 141 outputs from 1e-7 to 1 s" % v[2])
    return _LOGOUT["v"]


def order(e_coarse, e_fine, ratio=2.0):
    return math.log(e_coarse / e_fine) / math.log(ratio)


def plan():
    """(name, what, function returning the measured value, lower, upper)."""
    def lv1():
        a, _ = check_leveling_1d(W, X, nx=16, dt=2e-6, steps=2000)
        b, _ = check_leveling_1d(W, X, nx=32, dt=2e-6, steps=2000)
        return order(a[0], b[0])
    def lv2():
        a, _ = check_leveling_2d(W, X, nx=12, dt=2e-6, steps=1000)
        b, _ = check_leveling_2d(W, X, nx=24, dt=2e-6, steps=1000)
        return order(a[0], b[0])
    def lvt():
        T, r = 0.16, {}
        for dt in (8e-3, 4e-3, 2.5e-4):
            v, _ = check_leveling_1d(W, X, nx=32, dt=dt, steps=int(round(T / dt)))
            r[dt] = v[2]
        return order(abs(r[8e-3] - r[2.5e-4]), abs(r[4e-3] - r[2.5e-4]))
    def lvs():
        a, _ = check_leveling_slip_1d(W, X, nx=16, dt=2e-6, steps=2000)
        b, _ = check_leveling_slip_1d(W, X, nx=32, dt=2e-6, steps=2000)
        return order(a[0], b[0])
    def lvst():
        T, r = 0.064, {}
        for dt in (3.2e-3, 1.6e-3, 1e-4):
            v, _ = check_leveling_slip_1d(W, X, nx=32, dt=dt, steps=int(round(T / dt)))
            r[dt] = v[2]
        return order(abs(r[3.2e-3] - r[1e-4]), abs(r[1.6e-3] - r[1e-4]))
    def spin(p):
        return lambda: order(check_spinodal(W, X, p, nx=16)[0][0], check_spinodal(W, X, p, nx=32)[0][0])
    def spin_dt(p):
        def f():
            T, r = 4.0, {}
            for dt in (0.4, 0.2, 0.0125):
                r[dt] = check_spinodal(W, X, p, nx=32, dt=dt, steps=int(round(T / dt)))[0][2]
            return order(abs(r[0.4] - r[0.0125]), abs(r[0.2] - r[0.0125]))
        return f
    def lvcz():
        a, _ = check_leveling_1d(W, X, nx=32)
        b, _ = check_leveling_cz(W, X, nx=32)
        return abs(a[2] - b[2]) / a[2]
    #  The linear third-order scheme, whose damping and dispersion are
    #  known, and the default (Koren), which clips the extrema of the wave
    def wv1(height=None):
        def f():
            a, _ = check_wave_1d(W, X, nx=32, height=height)
            b, _ = check_wave_1d(W, X, nx=64, height=height)
            return order(a[0], b[0])
        return f
    def wvph(height=None, nx=32):
        return lambda: check_wave_1d(W, X, nx=nx, height=height)[0][1]
    def wvpho(height=None):
        return lambda: order(check_wave_1d(W, X, nx=32, height=height)[0][1],
                             check_wave_1d(W, X, nx=64, height=height)[0][1])
    def wv2(height=None):
        return lambda: check_wave_2d(W, X, nx=24, height=height)[0][1]
    def wvcz():
        a, _ = check_wave_cz(W, X, nx=32)
        b, _ = check_wave_cz(W, X, nx=64)
        return order(a[0], b[0])
    def wvczph():
        a, _ = check_wave_cz(W, X, nx=32)
        return a[1]
    def evo(fn):
        def f():
            T = 0.4
            a, _ = fn(W, X, dt=1e-2, steps=int(round(T / 1e-2)))
            b, _ = fn(W, X, dt=5e-3, steps=int(round(T / 5e-3)))
            return order(a, b)
        return f
    return [
        ("leveling_1d", "capillary decay, order in dx (16, 32 cells)", lv1, 1.9, 2.1),
        ("leveling_2d", "capillary decay in 2-D, order in dx (12, 24)", lv2, 1.9, 2.1),
        ("leveling_dt", "capillary decay, order in dt (Crank-Nicolson)", lvt, 1.9, 2.1),
        ("leveling_cz", "ConcZ decay rate against well-mixed", lvcz, 0.0, 1e-8),
        ("leveling_slip_1d", "Navier slip (Lambda 0.25): decay (H0^3/3 + lam H0^2) k^4, order in dx (16, 32)",
         lvs, 1.9, 2.1),
        ("leveling_slip_dt", "the same, order in dt", lvst, 1.9, 2.1),
        ("energy_leveling", "Output.Energy, 30 % sinusoid, capillary, periodic: |F(T)-F(0)+int D| / int D, order in dt",
         lambda: energy_order("leveling", 32, 2e-3, 40)[2], 1.9, 2.1),
        ("energy_leveling_mono", "... F never rises (largest rise / fall of F, both dt)",
         lambda: max(energy_order("leveling", 32, 2e-3, 40)[0][0], energy_order("leveling", 32, 2e-3, 40)[1][0]), 0.0, 1e-12),
        ("energy_slip", "the same with Navier slip: order in dt", lambda: energy_order("slip", 32, 1e-3, 40)[2], 1.9, 2.1),
        ("energy_box", "capillary + vertical gravity between Symmetry walls: order in dt",
         lambda: energy_order("box", 32, 1e-3, 40)[2], 1.9, 2.1),
        ("energy_2d", "the same in 2-D, periodic: order in dt", lambda: energy_order("2d", 16, 1e-3, 40)[2], 1.9, 2.1),
        ("energy_dewetting", "dewetting (disjoining pressure, h0 = 5.5 h*) through rupture: order in dt",
         lambda: energy_order("dewetting", 128, 0.125, 320)[2], 1.9, 2.1),
        ("energy_dewetting_mono", "... F never rises (both dt)",
         lambda: max(energy_order("dewetting", 128, 0.125, 320)[0][0], energy_order("dewetting", 128, 0.125, 320)[1][0]),
         0.0, 1e-12),
        ("energy_saturate_mono", "capillary + disjoining, h0 = 2 h*, grown to a steady film over 40 s: F never rises, D >= 0",
         lambda: check_energy(W, X, "saturate", 64, 1.0, 40)[0][0], 0.0, 1e-12),
        ("eq_slip", "the equilibrium droplet (section 5.3) with slip 1 um, to 0.01 s: max|dh| / max h",
         lambda: check_eq_slip(W, X)[0][0], 0.0, 1e-4),
        ("leveling_slip_1d_err", "... its error at 32 cells (no slip: the rate is 2.5 times smaller)",
         lambda: check_leveling_slip_1d(W, X, nx=32, dt=2e-6, steps=2000)[0][0], 0.0, 2e-2),
        ("gravity_1d", "decay with vertical gravity, error at 32 cells",
         lambda: check_gravity_1d(W, X)[0][0], 0.0, 1e-2),
        ("potential_1d", "WritePotential: PHI against sigma lap h + Pi(h) - rho g h, initial field",
         lambda: check_potential(W, X)[0], 0.0, 1e-12),
        ("spinodal_1d_unstable", "disjoining pressure, h0 = 2 h*: growth rate, order in dx (16, 32)",
         spin(1.0 / 3.0), 1.9, 2.1),
        ("spinodal_1d_unstable_dt", "the same, order in dt", spin_dt(1.0 / 3.0), 1.9, 2.1),
        ("spinodal_1d_stable", "disjoining pressure, h0 = 1.25 h*: decay rate, order in dx (16, 32)",
         spin(2.0 / 3.0), 1.9, 2.1),
        ("spinodal_1d_stable_dt", "the same, order in dt", spin_dt(2.0 / 3.0), 1.9, 2.1),
        ("wave_1d", "UTOPIA: damping of the upwinded in-plane gravity, order in dx (32, 64)",
         wv1("UTOPIA"), 2.8, 3.3),
        ("wave_1d_speed", "UTOPIA: wave speed rho g h0^2/mu, error at 32 cells",
         wvph("UTOPIA"), 0.0, 1e-4),
        ("wave_2d_speed", "UTOPIA: wave speed along Acceleration in 2-D, error at 24 cells",
         wv2("UTOPIA"), 0.0, 2e-4),
        ("wave_1d_koren", "Koren (default): damping, order in dx (32, 64)", wv1(), 2.8, 3.6),
        ("wave_1d_koren_speed", "Koren: wave speed, error at 32 cells", wvph(), 0.0, 1e-3),
        ("wave_1d_koren_speed_order", "Koren: wave speed error, order in dx (32, 64)",
         wvpho(), 2.5, 3.5),
        ("wave_2d_koren_speed", "Koren: wave speed in 2-D, error at 24 cells", wv2(), 0.0, 2e-3),
        ("wave_cz", "ConcZ decay rate with in-plane gravity, order in dx (32, 64)",
         wvcz, 1.9, 2.3),
        ("wave_cz_speed", "ConcZ wave speed, error at 32 cells", wvczph, 0.0, 1e-2),
        ("spin_1d", "spin-off thinning h0/sqrt(1 + 4 rho w^2 h0^2 t/3mu)",
         lambda: check_spin_1d(W, X)[0], 0.0, 1e-5),
        ("evap_wm", "drying flat film against the ODE, order in dt (10, 5 ms)",
         evo(check_evap_wm), 1.9, 2.3),
        ("evap_cz", "the same with ConcZ, order in dt",
         evo(check_evap_cz), 1.9, 2.3),
        ("dry_cz", "ConcZ drying at Pe ~ 1, D(c): profile, order in Nz (16, 32)",
         lambda: order(check_dry_cz(W, X, nz=16)[0], check_dry_cz(W, X, nz=32)[0]), 1.8, 2.3),
        ("evap_cap", "evaporation capped at h/2 and at the solvent: h1 = 2/3 h0, solute kept",
         lambda: check_evap_cap(W, X)[0], 0.0, 1e-6),
        ("evap_threshold", "the cap acts from E dt > h/2, on h and on c: E dt/h0 = 0.3, 0.45, 0.75",
         lambda: max(check_evap_threshold(W, X, a)[0] for a in (0.3, 0.45, 0.75)), 0.0, 1e-6),
        ("evap_threshold_2d", "the same in 2-D (WM2D concentration equation)",
         lambda: max(check_evap_threshold(W, X, a, twoD=True)[0] for a in (0.45, 0.75)),
         0.0, 1e-6),
        ("evap_solvent", "the solvent limit E dt <= h (1 - c), c0 = 0.5: E dt/h0 = 0.3, 0.75, 1-D and 2-D",
         lambda: max(check_evap_threshold(W, X, a, twoD=t, c0=0.5)[0]
                     for a in (0.3, 0.75) for t in (False, True)), 0.0, 1e-6),
        ("dried_wm", "a film dried out to c = 1 (rate independent of c): solute kept, c <= 1, 1-D and 2-D",
         lambda: max(max(v[0], v[1] * 10.0, v[2] * 1e-3) for v in
                     (check_dried_out(W, X, "WM", t)[0] for t in (False, True))), 0.0, 1e-11),
        ("dried_wm_adapt", "... with dt lengthened by the time-step control (E dt/h0 = 0.05 at first)",
         lambda: max(max(v[0], v[1] * 10.0, v[2] * 1e-3) for v in
                     (check_dried_out(W, X, "WM", t, a=0.05, adaptive=True)[0] for t in (False, True))),
         0.0, 1e-11),
        ("dried_cz", "the same with ConcZ",
         lambda: max(max(v[0], v[1] * 10.0, v[2] * 1e-3) for v in
                     (check_dried_out(W, X, "CZ", t)[0] for t in (False, True))), 0.0, 1e-11),
        ("dried_cz_adapt", "... ConcZ with dt lengthened",
         lambda: max(max(v[0], v[1] * 10.0, v[2] * 1e-3) for v in
                     (check_dried_out(W, X, "CZ", t, a=0.05, adaptive=True)[0] for t in (False, True))),
         0.0, 1e-11),
        ("negpow_rim", "negPow evaporation across NegPowRmax, against its definition",
         lambda: check_negpow_rim(W, X)[0][0], 0.0, 1e-9),
        ("negpow_rim_slope", "... the slope carried into the blend (jump / slope)",
         lambda: check_negpow_rim(W, X)[0][1], 0.0, 0.2),
        ("negpow_rim_flat", "... and constant past NegPowRmax (spread / value)",
         lambda: check_negpow_rim(W, X)[0][2], 0.0, 1e-9),
        ("inlet_outlet", "FreshInlet/Outlet, a bump carried out: max |VolErr|, |ResinVolErr|",
         lambda: max(check_inlet_outlet(W, X, nx=100)[0][:2]), 0.0, 1e-12),
        ("inlet_outlet_left", "the same: film uniform again, c untouched, after the bump left",
         lambda: max(check_inlet_outlet(W, X, nx=100)[0][2:]), 0.0, 1e-3),
        ("solute_balance", "... drying: max |VolErr|, |ResinVolErr| (UTOPIA, Koren, adaptive dt, no c)",
         lambda: open_balance([dict(variant=v) for v in ("UTOPIA", "Koren", "adaptive", "none")]),
         0.0, 1e-12),
        ("solute_balance_2d", "the same in 2-D: uniform in y; open on all four sides (corners)",
         lambda: open_balance([dict(twoD=True), dict(variant="Koren", twoD=True),
                               dict(nx=24, dt=5e-4, corners=True),
                               dict(variant="none", nx=24, dt=5e-4, corners=True)]), 0.0, 1e-12),
        ("wall_balance", "closed box, gravity + spin toward Symmetry walls: |VolErr|, 4 schemes, 1-D",
         lambda: check_wall_balance(W, X, False)[0], 0.0, 1e-12),
        ("wall_balance_2d", "the same in 2-D, every edge and corner (with the centrifugal force)",
         lambda: check_wall_balance(W, X, True)[0], 0.0, 1e-12),
        ("wall_balance_pinned", "the same toward a Pinning wall (EdgeBump), 1-D: |VolErr|",
         lambda: check_wall_balance(W, X, False, "Xmin")[0], 0.0, 1e-12),
        ("wall_balance_pinned_2d", "... and in 2-D, pinned at Ymin",
         lambda: check_wall_balance(W, X, True, "Ymin", schemes=("Koren", "KK"))[0], 0.0, 1e-12),
        ("wall_balance_slip", "the closed boxes with Navier slip (lam = 3 um): |VolErr|, 1-D, 2-D, pinned",
         lambda: max(check_wall_balance(W, X, False, slip=3e-6)[0], check_wall_balance(W, X, True, slip=3e-6)[0],
                     check_wall_balance(W, X, False, "Xmin", slip=3e-6)[0],
                     check_wall_balance(W, X, True, "Ymin", schemes=("Koren", "KK"), slip=3e-6)[0]), 0.0, 1e-12),
        ("wall_balance_cz", "ConcZ closed box, every term (body forces toward the walls): volume, 1-D, 2-D",
         lambda: max(check_wall_balance_cz(W, X, t, False, scheme=sc)[0][0]
                     for t, sc in ((False, "Koren"), (False, "QUICK"), (True, "Koren"))), 0.0, 1e-12),
        ("wall_balance_cz_solute", "... its solute, without and with evaporation",
         lambda: max(check_wall_balance_cz(W, X, t, ev, scheme=sc)[0][1]
                     for t, sc in ((False, "Koren"), (False, "QUICK"), (True, "Koren"))
                     for ev in (False, True)), 0.0, 1e-11),
        ("conc_overshoot", "a front c = 0.99 | 0.3 carried by the flow (UTOPIA, Central4; UTOPIA 2-D): c > 1 moved, |ResinVolErr|",
         lambda: conc_overshoot([("UTOPIA", 0.4, False), ("Central4", 0.4, False), ("UTOPIA", 0.4, True)]),
         0.0, 1e-12),
        ("conc_overshoot_koren", "... Koren at Courant 1.5 (limiter frozen at the iterate as solved): |ResinVolErr|",
         lambda: conc_overshoot([("Koren", 1.5, False)]), 0.0, 1e-12),
        ("conc_drying", "the front 0.9 | 0.3 on a film drying as it is carried (UTOPIA, Koren; UTOPIA 2-D): |ResinVolErr|",
         lambda: conc_drying([("UTOPIA", False), ("Koren", False), ("UTOPIA", True)]), 0.0, 1e-12),
        ("conc_uniform", "a uniform c in a closed box, nozzle, strong gravity (and spin, 1-D): max |c - c0|, 1-D and 2-D",
         lambda: max(check_conc_uniform(W, X, t)[0] for t in (False, True)), 0.0, 1e-12),
        ("conc_undershoot", "a front 0.7 | 0.005 (UTOPIA, KK, Central4; UTOPIA 2-D): c < 0 filled, |ResinVolErr|",
         lambda: conc_undershoot([("UTOPIA", 0.4, False), ("KK", 0.4, False), ("Central4", 0.4, False),
                                  ("UTOPIA", 0.4, True)]), 0.0, 1e-12),
        ("conc_undershoot_cz", "... ConcZ (UTOPIA, 1-D): c < 0 filled, |ResinVolErr|",
         lambda: conc_undershoot([("UTOPIA", 0.4, False)], cz=True), 0.0, 1e-12),
        ("spin_balance", "1-D closed box spinning: |VolErr|, |ResinVolErr| count the outflow into y",
         lambda: check_spin_balance(W, X)[0], 0.0, 1e-12),
        ("nozzle_balance", "a nozzle in a closed box, adaptive dt, and on a turning arm thickening the film 100x: |VolErr|",
         lambda: max(nozzle_balance(False)[0], nozzle_balance(True)[0], nozzle_balance("thick")[0]), 0.0, 1e-12),
        ("nozzle_solute", "... |ResinVolErr| counts the solute that came with it",
         lambda: max(nozzle_balance(False)[1], nozzle_balance(True)[1], nozzle_balance("thick")[1]), 0.0, 1e-12),
        ("nozzle_supplied", "the same, the liquid supplied at Source.Conc (ConcType \"supplied\"; stagnated steps cut): |VolErr|, |ResinVolErr|",
         lambda: max(max(nozzle_supplied(t)[0] for t in (False, True, "thick")),
                     max(nozzle_supplied(t)[1] for t in (False, True, "thick"))), 0.0, 1e-12),
        ("nozzle_mixing", "liquid at Source.Conc on a film that does not flow: c = (h0 c0 + S t cs) / (h0 + S t), 1-D, 2-D",
         lambda: nozzle_mixing(), 0.0, 1e-12),
        ("leveling_logout_times", "LogSchedule with ExactTimes: |output time - listed| / listed",
         lambda: leveling_logout()[0], 0.0, 1e-12),
        ("leveling_logout", "... decay rate at each of the 141 outputs, 64 cells (as at a fixed dt)",
         lambda: leveling_logout()[1], 0.0, 2e-3),
        ("leveling_adaptive", "adaptive dt: decay at the reported times, |time shift| / t",
         lambda: leveling_adaptive()[0], 0.0, 5e-4),
        ("leveling_adaptive_clock", "... |time run - sum of the steps solved| / time run",
         lambda: leveling_adaptive()[1], 0.0, 1e-12),
        ("leveling_adaptive_cz", "the same with ConcZ: |time shift| / t",
         lambda: leveling_adaptive("CZ")[0], 0.0, 5e-4),
        ("leveling_adaptive_cz_clock", "... clock = sum of the steps solved",
         lambda: leveling_adaptive("CZ")[1], 0.0, 1e-12),
        ("nozzle_arm", "rotating nozzle: where the liquid lands (steady turn, stop)",
         lambda: check_nozzle_arm(W, X)[0], 0.0, 2e-3),
        ("nozzle_area", "2-D nozzle on square and on long thin cells: |supplied / (InletRate pi R^2 t) - 1|",
         lambda: check_nozzle_area(W, X)[0], 0.0, 1e-2),
        ("slab_cz", "ConcZ profile against the slab series, Nz = 20",
         lambda: check_slab_cz(W, X, nz=20)[0], 0.0, 1e-3),
    ]


W = X = None


def main():
    global W, X
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", action="append", help="run only this check (repeatable)")
    args = ap.parse_args()
    X = os.path.abspath(args.exe)
    W = args.workdir
    shutil.rmtree(W, ignore_errors=True)
    bad = 0
    for name, what, fn, lo, hi in plan():
        if args.only and name not in args.only:
            continue
        try:
            v = fn()
        except Exception as e:          # a failed run surfaces as a None result
            print("FAIL %-14s %s: %s" % (name, what, e))
            bad += 1
            continue
        ok = lo <= v <= hi
        bad += not ok
        print("%s %-14s %-62s %.3g  (allowed %g .. %g)" % ("ok  " if ok else "FAIL", name, what, v, lo, hi))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
