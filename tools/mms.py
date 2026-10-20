#!/usr/bin/env python3
"""The method of manufactured solutions for the film-height equation.

    python3 tools/mms.py --exe build/bin/yuzen --workdir /tmp/mms
    python3 tools/mms.py --exe build/bin/yuzen --workdir /tmp/mms --only wm_1d
    python3 tools/mms.py --exe build/bin/yuzen --workdir /tmp/mms --scheme Koren --each

--scheme sets [Numerics.Convection] Concentration and (well-mixed only)
Height for every run; --conc and --height set one of them. The thresholds
are the same whatever the scheme (--lo).

With [Numerics.ManufacturedSolution] enabled the solver adds to the
film-height equation the source (model/ManufacturedSource.cpp) that makes

    h = Mean + Amplitude cos(kx x - Omega t) cos(ky y)

an exact solution. The amplitude here is 30 % of the mean, so the mobility
h^3 varies by a factor of 6 over a wavelength and every nonlinear part of
the discretization is exercised -- which the linear checks of
tools/analytic.py, at an amplitude of 1e-3, cannot do. The capillary
pressure, vertical gravity and in-plane gravity are all on. The error is
the largest difference between the computed film and the manufactured one
at the end of the run, over the cell centers, relative to the amplitude.

For each model and dimension the error is measured on two grids at a time
step small enough not to matter, and at two time steps on a grid fine
enough not to matter, and it is the orders of convergence that are
checked: second in dx and second in dt.

    wm_1d    well-mixed, 1-D
    wm_2d    well-mixed, 2-D
    cz_1d    ConcZ, 1-D
    cz_2d    ConcZ, 2-D
    wm_all_*, cz_all_*   every term of the model (see full_terms, cz_terms)
    cz_conc_*            ConcZ with a concentration that varies in the
                         plane, mu(c), sigma(c), D(c) and Marangoni
                         (cz_conc_terms)
    cz_strat_*           the same with a concentration that varies through
                         the film as well (cz_strat_terms)
    cz_strat_1d_koren    cz_strat_1d with the limited scheme Koren for the
                         concentration (the default is UTOPIA)
    cz_strat_1d_central4 cz_strat_1d with the central scheme Central4
    wm_slip_1d, wm_all_slip_1d   wm_1d and wm_all_1d with Navier slip
                         (Physics.Slip, lam = 2.5 um, half the mean film):
                         mobility h^3/3 + lam h^2, h^2/2 + lam h for the
                         Marangoni flux (SLIP)

Only the Python standard library is used, so the test runs wherever ctest
does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import math
import os
import shutil
import sys
try:
    import tomllib                              # Python 3.11
except ImportError:                             # pragma: no cover
    tomllib = None

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from concz_vs_wm import edit_toml               # noqa: E402
from equivalence import arr, fixed_steps, run   # noqa: E402
from analytic import case, quiet, wave_film, fields  # noqa: E402

MEAN, AMP = 5.0e-6, 1.5e-6        # m; TotalHeight 10 um, Base 0.5, Amplitude 0.15
LX = 0.25e-3                      # wavelength in x, m
LY = 0.375e-3                     # and in y (2-D)
OMEGA = 60.0                      # rad/s
T_END = 0.02                      # s


def exact(x, y, t, twoD):
    ky = 2 * math.pi / LY if twoD else 0.0
    return MEAN + AMP * math.cos(2 * math.pi / LX * x - OMEGA * t) * math.cos(ky * y)


#  Every term of the well-mixed model at once, each tuned so that it moves
#  h_t or c_t by about as much as the others (5e-5 m/s and 3 /s; the film
#  itself moves by A Omega = 9e-5 m/s): the capillary pressure over a
#  substrate with a relief of 2 um, both gravities, the centrifugal term
#  about an axis off the domain's center, the disjoining pressure with the
#  film at 1.75 to 3.25 times the precursor, the Marangoni flux, a viscosity
#  that changes by a factor e^0.6 and a diffusivity by e^-0.4 over the
#  concentration's range, the convection and diffusion of the
#  concentration, and the evaporation (with the precursor correction), which
#  both thins the film and enriches the solute.
#  Navier slip of half the mean film, so that lam h^2 is as large as h^3/3
#  and the face mobility of the capillary term is far from h^3
SLIP = {"Physics.Slip": {"Enable": "true", "Length": "2.5e-6"}}


def full_terms(twoD, law):
    t = {
        "Initial.SurfaceShape": {"PrecursorRatio": "0.2"},
        "Fluid.Viscosity": {"Exponent": "3.0", "BaseConc": "0.5"},
        "Fluid.SurfaceTension": {"Coef": "0.004", "BaseConc": "0.5"},
        "Fluid.Diffusivity": {"BaseValue": "3.0e-8", "Exponent": "-2.0", "BaseConc": "0.5"},
        "Physics.ConcField": {"EnableConvection": "true", "EnableDiffusion": "true"},
        "Physics.VerticalGravity": {"Enable": "true", "Acceleration": "3.0e4"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(6.0e3, 4.0e3 if twoD else 0.0)},
        "Physics.Centrifugal": {"Enable": "true", "RotationRate": repr(math.sqrt(1.2e7)),
                                "RotationAxis": arr(-0.3e-3, 0.1e-3)},
        "Physics.SolutoCapillarity": {"Enable": "true"},
        "Physics.DisjoiningPressure": {"Enable": "true", "ContactAngle": "25.0",
                                       "ExponentN": "3.0", "ExponentM": "2.0"},
        "Physics.Evaporation": {"Enable": "true", "NegPowExponent": "0.25", "NegPowRmax": "1.0e-3",
                                "PowerN": "2.0", "ConcR": "0.1"},
        "Numerics.ManufacturedSolution": {"ConcMean": "0.5", "ConcAmplitude": "0.1",
                                          "ConcOmega": "40.0", "ConcPhase": "0.7",
                                          "ConcWavelength": arr(LX, LY),
                                          "SubstrateAmplitude": "2.0e-6",
                                          "SubstrateWavelength": arr(LX, LY)},
    }
    if law == "PowerN":
        t["Physics.Evaporation"].update({"Coef": "3.0e-3", "ConcFunc": '"PowerN"', "SpatialFunc": '"uniform"'})
    else:
        t["Physics.Evaporation"].update({"Coef": "6.0e-4", "ConcFunc": '"Raoult"', "SpatialFunc": '"negPow"'})
    if twoD:
        #  In 2-D a tenth of the diffusivity and half the relief. The
        #  Kawamura-Kuwahara upwinding of the convection of c (and of the
        #  body forces and the Marangoni flux in the film-height equation)
        #  adds a third-order dissipation, -|u| dx^3/4 c_xxxx, which is of the
        #  opposite sign to the second-order error of the diffusion
        #  (+D dx^2/12 c_xxxx); the two cancel at a cell Peclet number u dx/D
        #  near 1/3, and around that point the order of c is below 2 on every
        #  grid 2-D can afford: 1.66, 1.58, 1.81 on 16 -> 32 -> 64 -> 128
        #  with the diffusivity of 1-D, 1.99, 1.98, 1.95 with the upwind
        #  part of the schemes switched off (Central2). With D = 3e-9 the
        #  cell Peclet number on 24 to 64 cells is well above that point,
        #  and the order of c comes down toward 2 from above (2.6, 2.6,
        #  2.5, 2.3 on 24 -> 32 -> 48 -> 64 -> 96, 2.0 on 96 -> 128, where
        #  the cancellation is approached again); D k^2 is still 7 % of the
        #  rate of the manufactured c. In the film-height equation the
        #  upwinding has the same effect against the second-order error of
        #  the capillary term, and over the relief of 2 um it kept h below 2
        #  as well (1.96, 1.89 on 32 -> 64 -> 128); over 1 um h converges at
        #  2.12, 2.11, 2.08, 2.06, 1.97 on 24 -> 32 -> 48 -> 64 -> 96 -> 128
        #  (KK for the concentration, Koren for the film height).
        t["Fluid.Diffusivity"]["BaseValue"] = "3.0e-9"
        t["Numerics.ManufacturedSolution"]["SubstrateAmplitude"] = "1.0e-6"
    return t


#  ConcZ has neither the disjoining pressure nor a concentration that varies
#  in the plane here (cz_conc_terms adds that); everything else that acts on
#  the film does.
def cz_terms(twoD):
    return {
        "Physics.VerticalGravity": {"Enable": "true", "Acceleration": "3.0e4"},
        "Physics.HorizontalGravity": {"Enable": "true", "Acceleration": arr(6.0e3, 4.0e3 if twoD else 0.0)},
        "Physics.Centrifugal": {"Enable": "true", "RotationRate": repr(math.sqrt(1.2e7)),
                                "RotationAxis": arr(-0.3e-3, 0.1e-3)},
        "Numerics.ManufacturedSolution": {"SubstrateAmplitude": "2.0e-6",
                                          "SubstrateWavelength": arr(LX, LY)},
    }


#  ConcZ with a concentration that varies strongly in the plane -- from 0.3
#  to 0.7, the same at every height in the film -- and everything that it
#  drives: the viscosity changes by a factor e^1.2 and the diffusivity by
#  e^-0.8 over that range, the surface tension enters the capillary
#  mobility, and its gradient drives the Marangoni flux, all on top of the
#  terms of cz_terms. The concentration has its own source, which depends on
#  the height zeta in the film through the velocity profile, and its error
#  is measured at every zeta. Evaporation (whose surface flux would put a
#  profile through the thickness) stays out.
def cz_conc_terms(twoD):
    t = cz_terms(twoD)
    t.update({
        "Fluid.Viscosity": {"Exponent": "3.0", "BaseConc": "0.5"},
        "Fluid.SurfaceTension": {"Coef": "0.004", "BaseConc": "0.5"},
        "Fluid.Diffusivity": {"BaseValue": "3.0e-8", "Exponent": "-2.0", "BaseConc": "0.5"},
        "Physics.SolutoCapillarity": {"Enable": "true"},
    })
    t["Numerics.ManufacturedSolution"].update({"ConcMean": "0.5", "ConcAmplitude": "0.2",
                                               "ConcOmega": "40.0", "ConcPhase": "0.7",
                                               "ConcWavelength": arr(LX, LY)})
    return t


#  ConcZ with a concentration that varies through the film as well: on top of
#  the c(x, y, t) of cz_conc_terms, 0.1 cos(kx x - 30 t) cos(ky y) cos(pi zeta),
#  which has dc/dzeta = 0 at the substrate and at the surface -- the
#  conditions ConcZ imposes there without evaporation. The surface value then
#  differs from the depth mean by 0.1 (and its gradient from the mean's), and
#  the viscosity changes by e^0.6 through the thickness, so the surface
#  tension, the Marangoni stress and the depth integral of 1/mu that set the
#  film's flux are each checked where the surface value, the mean and the
#  profile differ; the vertical advection of the mapping, its metric terms and
#  the diffusion through the film enter the concentration equation. The
#  diffusivity is 20 times lower than in cz_conc_terms, so that the
#  diffusion through the film (D pi^2/h^2, 6e2 /s) stays within an order of
#  the other terms of c_t rather than a thousand times larger -- and no
#  lower, because the state must be one the model keeps. Where the surface
#  tension rises with c and c rises toward the surface (sigma_c dc/dzeta >
#  0, over half of each wavelength here), the Marangoni flow and the
#  viscosity's dependence on c make short waves grow, and only the diffusion
#  cuts them off: tools/cz_stability.py, a frozen-coefficient analysis of the
#  linearized equations about this state, gives growth of up to 7e2 /s at
#  3e-10 m^2/s (waves of 1 to 3 cells of 128 in 1-D), which an
#  upwind-biased scheme damped on these grids and a central one did not.
#  At 1.5e-9 every wave from a quarter of the state's wavelength down to
#  two cells of a grid twice as fine as the finest here decays, in 1-D
#  (256 cells) and in 2-D (64, in every direction), at 3.5 /s or faster;
#  2-D needs at least about 1.25e-9. The cells through the
#  film are refined with the grid (_Strat.Nz: nx -> Nz); the dt runs share
#  one, whose error cancels between them. Each of h, c and u must converge on
#  its own in dt as well as in dx.
def cz_strat_terms(twoD, nz):
    t = cz_conc_terms(twoD)
    t["Fluid.Diffusivity"] = {"BaseValue": "1.5e-9", "Exponent": "-2.0", "BaseConc": "0.5"}
    t["Numerics.ManufacturedSolution"].update({"ConcStratAmplitude": "0.1", "ConcStratOmega": "30.0",
                                               "ConcStratWavelength": arr(LX, LY)})
    t["_Strat"] = {"Nz": nz}
    return t


def exact_c(x, y, t, twoD, ms, zeta=0.0):
    kx = 2 * math.pi / LX
    ky = 2 * math.pi / LY if twoD else 0.0
    c = (float(ms["ConcMean"]) + float(ms["ConcAmplitude"])
         * math.cos(kx * x - float(ms["ConcOmega"]) * t + float(ms["ConcPhase"])) * math.cos(ky * y))
    if "ConcStratAmplitude" in ms:
        #  (ConcStratWavelength is (LX, LY) in cz_strat_terms)
        c += (float(ms["ConcStratAmplitude"]) * math.cos(kx * x - float(ms["ConcStratOmega"]) * t)
              * math.cos(ky * y) * math.cos(math.pi * zeta))
    return c


def exact_u(x, y, zeta, t, twoD, inp):
    """The in-plane velocity of ConcZ at the height zeta in [0, 1] of the
    film, for the manufactured film, substrate and concentration, from the
    parameters of the run's input: the lubrication profile with the
    viscosity of the local concentration,

        u = int_0^zeta [h^2 (1 - s) b + h tau] / mu(c(s)) ds,
        b = sigma(c_s) grad lap(h + s) - rho g grad(h + s) + rho g_par
            + rho Omega^2 (r - axis),       tau = ms sigma_c grad c_s,

    c_s the surface value c(zeta = 1), whose depth integral is the flux of
    the film-height equation. With c the same at every zeta it is
    [h^2 (zeta - zeta^2/2) b + h zeta tau] / mu(c). Written out by hand
    here, independently of the solver and of mms_generate.py; the integral
    by Simpson's rule on 64 intervals."""
    def on(sec):
        return sec.get("Enable", True)
    fl, ph, m = inp["Fluid"], inp["Physics"], inp["Numerics"]["ManufacturedSolution"]
    tp = 2 * math.pi
    kx, ky = tp / m["Wavelength"][0], (tp / m["Wavelength"][1] if twoD else 0.0)
    ckx, cky = tp / m["ConcWavelength"][0], (tp / m["ConcWavelength"][1] if twoD else 0.0)
    skx, sky = tp / m["SubstrateWavelength"][0], (tp / m["SubstrateWavelength"][1] if twoD else 0.0)
    A, SA, CA = m["Amplitude"], m.get("SubstrateAmplitude", 0.0), m["ConcAmplitude"]
    CB = m.get("ConcStratAmplitude", 0.0)
    bl = m.get("ConcStratWavelength", [1.0, 1.0])
    bkx, bky = tp / bl[0], (tp / bl[1] if twoD else 0.0)
    th = kx * x - m["Omega"] * t
    h = m["Mean"] + A * math.cos(th) * math.cos(ky * y)
    hx, hy = -A * kx * math.sin(th) * math.cos(ky * y), -A * ky * math.cos(th) * math.sin(ky * y)
    sx = -SA * skx * math.sin(skx * x) * math.cos(sky * y)
    sy = -SA * sky * math.cos(skx * x) * math.sin(sky * y)
    tc = ckx * x - m["ConcOmega"] * t + m["ConcPhase"]
    tb = bkx * x - m.get("ConcStratOmega", 0.0) * t
    c0 = m["ConcMean"] + CA * math.cos(tc) * math.cos(cky * y)       # the part the same at every zeta
    cb = CB * math.cos(tb) * math.cos(bky * y)                        # times cos(pi zeta)
    #  The surface value and its gradient (cos(pi) = -1)
    c = c0 - cb
    cx = -CA * ckx * math.sin(tc) * math.cos(cky * y) + CB * bkx * math.sin(tb) * math.cos(bky * y)
    cy = -CA * cky * math.cos(tc) * math.sin(cky * y) + CB * bky * math.cos(tb) * math.sin(bky * y)
    ms = 1.0 if on(ph["SolutoCapillarity"]) else 0.0
    st, vi = fl["SurfaceTension"], fl["Viscosity"]
    sig = st["BaseValue"] + ms * st["Coef"] * (c - st["BaseConc"])

    def rmu(z):
        return math.exp(-vi["Exponent"] * (c0 + cb * math.cos(math.pi * z) - vi["BaseConc"])) / vi["BaseValue"]
    rho = fl["Density"]["Value"]
    cap = 1.0 if on(ph.get("LaplacePressure", {})) else 0.0
    g = ph["VerticalGravity"]["Acceleration"] if on(ph["VerticalGravity"]) else 0.0
    gv = ph["HorizontalGravity"]["Acceleration"] if on(ph["HorizontalGravity"]) else [0.0, 0.0]
    ce = ph["Centrifugal"]
    om2 = ce["RotationRate"] ** 2 if on(ce) else 0.0
    ax = (list(ce.get("RotationAxis", [])) + [0.0, 0.0])[:2]
    k2, sk2 = kx * kx + ky * ky, skx * skx + sky * sky
    bx = (cap * sig * (-k2 * hx - sk2 * sx) - rho * g * (hx + sx) + rho * gv[0]
          + rho * om2 * (x - ax[0]))
    by = (cap * sig * (-k2 * hy - sk2 * sy) - rho * g * (hy + sy) + rho * (gv[1] if twoD else 0.0)
          + rho * om2 * (y - ax[1]))
    #  int_0^zeta (1 - s)/mu ds and int_0^zeta 1/mu ds
    n = 64
    i1 = i0 = 0.0
    for j in range(n + 1):
        z = zeta * j / n
        w = (1 if j in (0, n) else (4 if j % 2 else 2)) * zeta / (3 * n)
        r = rmu(z)
        i1 += w * (1 - z) * r
        i0 += w * r
    prof, shear = h * h * i1, ms * st["Coef"] * h * i0
    return prof * bx + shear * cx, prof * by + shear * cy

#  [Numerics.Convection] for every run (main() fills it from the options)
CONVECTION = {}


def run_mms(workdir, exe, model, twoD, nx, dt, extra=None, tag="", after=None):
    src = case("concz_1d_reg" if model == "CZ" else "leveling_1d_reg")
    s = quiet(src, twoD=twoD)
    s = edit_toml(s, {}, {"Grid": {"Nx": str(nx), "Xmin": repr(-LX / 2), "Xmax": repr(LX / 2)},
                          "Physics.ConcField": {"Model": '"%s"' % model},
                          "Physics.VerticalGravity": {"Enable": "true", "Acceleration": "2000.0"},
                          "Physics.HorizontalGravity": {"Enable": "true",
                                                        "Acceleration": arr(3.0e4, 2.0e4 if twoD else 0.0)}})
    extra = dict(extra or {})
    strat = extra.pop("_Strat", None)
    #  The convection schemes of --scheme / --conc / --height (CZ has no
    #  upwinded term in its film-height equation, and refuses Height)
    cv = {}
    if CONVECTION.get("Concentration"):
        cv["Concentration"] = '"%s"' % CONVECTION["Concentration"]
    if CONVECTION.get("Height") and model != "CZ":
        cv["Height"] = '"%s"' % CONVECTION["Height"]
    if cv:
        extra["Numerics.Convection"] = cv
    if model == "CZ":
        #  Where the concentration is the same at every height, a few cells
        #  through the thickness are as good as many; where it is not, their
        #  number goes with the grid (cz_strat_terms)
        s = edit_toml(s, {}, {"Grid": {"Nz": str(strat["Nz"][nx]) if strat else "4"}})
    ny = 1
    if twoD:
        ny = int(round(nx * LY / LX))
        s = edit_toml(s, {}, {"Grid": {"Dimension": '"xyh"', "Ny": str(ny),
                                       "Ymin": repr(-LY / 2), "Ymax": repr(LY / 2)}})
    s = wave_film(s, MEAN / 10e-6, AMP / 10e-6, LX, LY if twoD else None)
    s = edit_toml(s, {}, {"Numerics.ManufacturedSolution": {
        "Enable": "true", "Mean": repr(MEAN), "Amplitude": repr(AMP),
        "Omega": repr(OMEGA), "Wavelength": arr(LX, LY)}})
    if extra:
        s = edit_toml(s, {}, extra)
    steps = int(round(T_END / dt))
    s = fixed_steps(s, dt, steps, T_END)
    if after:
        #  Keys set over those of fixed_steps (tools/features.py coupling_cz
        #  runs the state with the time-step control on)
        s = edit_toml(s, {}, after)
    d = os.path.join(workdir, "%s%s_%s_n%d_dt%g" % (model, tag, "2d" if twoD else "1d", nx, dt))
    rc, log = run(exe, d, "a", s, 1800)
    if rc != 0:
        raise RuntimeError("run failed:\n" + log[-600:])
    t, (f, nx_, ny_, nz_) = fields(os.path.join(d, "a"), T_END)[-1]
    ms = (extra or {}).get("Numerics.ManufacturedSolution", {})
    ca = float(ms.get("ConcAmplitude", "0"))
    worst, hs = 0.0, []
    parts = {"h": 0.0}
    #  Each quantity on its own, for the dt orders of cz_strat_*: h, c and
    #  u (less the exact profile) as computed; the grid error is the same in
    #  all the runs and drops out of their differences.
    each = {"h": [], "c": [], "u": []}
    #  In ConcZ with a manufactured concentration the velocity profile that
    #  carries it is checked as well, at the cell centers of every real level
    #  (1 .. nz_-2; 0 and nz_-1 are the substrate and the surface), at the
    #  time of the file. (The field holds the velocity of that time, not of
    #  the last step's midpoint half a step earlier; compared there, the grid
    #  error of the velocity would drift between runs with different steps,
    #  and its dt order could not be measured.)
    vel = []
    if model == "CZ" and ca > 0 and tomllib is not None:
        with open(os.path.join(d, "a.toml"), "rb") as fh:
            inp = tomllib.load(fh)
        tv, nzr = t, nz_ - 2
    for j in (range(1, ny_ - 1) if ny_ > 1 else [0]):
        for i in range(1, nx_ - 1):
            r = j * nx_ + i
            x = f["xx"][r] * 1e-3
            y = f["yy"][r] * 1e-3 if "yy" in f else 0.0
            h = f["HH"][r] * 1e-6
            hs.append(h)
            each["h"].append(h / AMP)
            parts["h"] = max(parts["h"], abs(h - exact(x, y, t, twoD)) / AMP)
            if ca > 0:
                #  The concentration, measured against its own amplitude, at
                #  every height in the film (ConcZ writes a row per level)
                for k in range(nz_):
                    #  (level k at zeta = (k - 1/2)/Nz; the ghosts 0 and
                    #  nz_ - 1 mirror the levels next to them, as the
                    #  manufactured c does about zeta = 0 and 1)
                    ce = exact_c(x, y, t, twoD, ms, (k - 0.5) / (nz_ - 2))
                    c = f["CC"][k * nx_ * ny_ + r]
                    hs.append(c * AMP / ca)
                    each["c"].append(c / ca)
                    parts["c"] = max(parts.get("c", 0.0), abs(c - ce) / ca)
            if model == "CZ" and ca > 0 and tomllib is not None:
                for k in range(1, nz_ - 1):
                    q = k * nx_ * ny_ + r
                    ue = exact_u(x, y, (k - 0.5) / nzr, tv, twoD, inp)
                    vel.append(((k, 0), f["Ux"][q], ue[0]))
                    if twoD:
                        vel.append(((k, 1), f["Uy"][q], ue[1]))
    each["u"] = [(key, a - e, e) for key, a, e in vel]
    if vel:
        #  Each component at each level against its own largest value: the
        #  velocity near the substrate is a fraction of that at the surface,
        #  and an error there would be lost against the surface's scale.
        scale, err = {}, {}
        for key, a, e in vel:
            scale[key] = max(scale.get(key, 0.0), abs(e))
            err[key] = max(err.get(key, 0.0), abs(a - e))
        parts["u"] = max(err[key] / scale[key] for key in err)
    worst = max(parts.values())
    return worst, hs, parts, each


#  (name, model, 2-D, grids for the dx order (at dt_fine), dt_fine,
#   grid for the dt order, time steps for the dt order)
PLAN = [("wm_1d", "WM", False, (32, 64), 2.5e-5, 128, (1e-3, 5e-4), None),
        #  (wm_2d: h on 32 cells at 2.5e-4 s is 0.4 % off that at 1e-4 s)
        ("wm_2d", "WM", True, (16, 32), 2.5e-4, 32, (2e-3, 1e-3), None),
        ("cz_1d", "CZ", False, (32, 64), 2.5e-5, 128, (1e-3, 5e-4), None),
        ("cz_2d", "CZ", True, (16, 32), 1e-4, 32, (2e-3, 1e-3), None),
        ("wm_all_1d", "WM", False, (64, 128), 2.5e-5, 128, (1e-3, 5e-4), full_terms(False, "PowerN")),
        #  wm_all_2d on 24 and 48 cells, at a step whose own error is a
        #  hundredth of theirs (1 % of c on 48 cells at 2.5e-4 s against
        #  1e-4 s)
        ("wm_all_2d", "WM", True, (24, 48), 2.5e-4, 32, (2e-3, 1e-3), full_terms(True, "Raoult")),
        ("wm_slip_1d", "WM", False, (32, 64), 2.5e-5, 128, (1e-3, 5e-4), dict(SLIP)),
        ("wm_all_slip_1d", "WM", False, (64, 128), 2.5e-5, 128, (1e-3, 5e-4),
         dict(full_terms(False, "PowerN"), **SLIP)),
        ("cz_all_1d", "CZ", False, (32, 64), 2.5e-5, 128, (1e-3, 5e-4), cz_terms(False)),
        ("cz_all_2d", "CZ", True, (16, 32), 1e-4, 32, (2e-3, 1e-3), cz_terms(True)),
        ("cz_conc_1d", "CZ", False, (32, 64), 2.5e-5, 128, (1e-3, 5e-4), cz_conc_terms(False)),
        #  (cz_conc_2d: the dx runs at the dt order's reference step, which
        #  with the coupling of the film and the concentration implicit has
        #  an error of 1e-4 of theirs, so that the finer is made once)
        ("cz_conc_2d", "CZ", True, (16, 32), 2.5e-4, 32, (2e-3, 1e-3), cz_conc_terms(True)),
        #  cz_strat_*: the dt order in 1-D on the grid of cz_conc_1d. (The
        #  Marangoni velocity of the solute has to carry the flux the film
        #  moves with; one that does not limits the step as dx^2 -- at 128
        #  cells in 1-D 5e-4 s diverges.) The dx runs take a step whose own
        #  error is a hundredth of theirs. In 2-D the diffusion through the
        #  film, 2.2e3 /s where the film is thinnest and c lowest, makes 2e-3 s
        #  a step of 4.4 relaxation times: u converged at 1.7 on 2e-3 -> 1e-3,
        #  2.3 on 1e-3 -> 5e-4. The dt order is therefore on 1e-3 -> 5e-4,
        #  and on 16 cells, where the reference at 1.25e-4 costs an eighth of
        #  what it would on 32. Nz is 6 and 12: with 4 and 8 the error of c
        #  through the film held its dx order at 1.83 (1.88 with 6 and 12,
        #  1.90 with 8 and 16).
        ("cz_strat_1d", "CZ", False, (32, 64), 1e-4, 128, (1e-3, 5e-4),
         cz_strat_terms(False, {32: 8, 64: 16, 128: 16})),
        ("cz_strat_2d", "CZ", True, (16, 32), 2.5e-4, 16, (1e-3, 5e-4),
         cz_strat_terms(True, {16: 6, 32: 12})),
        #  cz_strat_1d with a limited scheme for the solute, whose stencil
        #  through the thickness (the far-upwind cell of a face, read from
        #  the values themselves) none of the other cases reaches: the
        #  advection tests carry the solute in the plane only. Its dt order
        #  is taken on 2.5e-4 -> 1.25e-4. The conservative form solves the
        #  step until the limiter is that of the midpoint of the step itself
        #  (CZ_ConcEquation_2D.cpp, SettleLimiter): the implicit midpoint
        #  rule of a right-hand side with kinks, whose asymptotic range
        #  starts later. On 1e-3 -> 5e-4 (against 1.25e-4) c, h and u
        #  converge at 2.05, 2.22 and 3.26, u not yet in its asymptotic
        #  range; on 2.5e-4 -> 1.25e-4 at 2.06, 2.14 and 2.79.
        ("cz_strat_1d_koren", "CZ", False, (32, 64), 1e-4, 128, (2.5e-4, 1.25e-4),
         dict(cz_strat_terms(False, {32: 8, 64: 16, 128: 16}),
              **{"Numerics.Convection": {"Concentration": '"Koren"'}})),
        #  cz_strat_1d with the central scheme Central4 for the solute,
        #  which dissipates nothing in any direction: on a state the model
        #  does not keep it follows the growth of the shortest waves
        #  (cz_strat_terms; at 3e-10 m^2/s c converged at 0.48 in dt).
        ("cz_strat_1d_central4", "CZ", False, (32, 64), 1e-4, 128, (1e-3, 5e-4),
         dict(cz_strat_terms(False, {32: 8, 64: 16, 128: 16}),
              **{"Numerics.Convection": {"Concentration": '"Central4"'}}))]


def dt_each(runs, ref):
    """The error of each of h, c and u in each run of the dt order against
    the reference run: h and c relative to their amplitudes, u to the largest
    exact value of each component at each level, as in the dx order."""
    out = []
    for r in runs:
        e = {}
        for k in ("h", "c"):
            if ref[k]:
                e[k] = max(abs(a - b) for a, b in zip(r[k], ref[k]))
        if ref["u"]:
            scale, err = {}, {}
            for (key, d, ex), (_, d0, _) in zip(r["u"], ref["u"]):
                scale[key] = max(scale.get(key, 0.0), abs(ex))
                err[key] = max(err.get(key, 0.0), abs(d - d0))
            e["u"] = max(err[key] / scale[key] for key in err)
        out.append(e)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", action="append")
    ap.add_argument("--lo", type=float, default=1.8, help="lowest order accepted")
    ap.add_argument("--scheme", help="Numerics.Convection: this scheme for the concentration"
                    " and (well-mixed) the film height, e.g. Koren")
    ap.add_argument("--conc", help="Numerics.Convection.Concentration alone")
    ap.add_argument("--height", help="Numerics.Convection.Height alone (well-mixed cases)")
    ap.add_argument("--each", action="store_true",
                    help="show the dt order of each of h, c and u in every case (checked only"
                         " where it always is, cz_strat_*)")
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
    for name, model, twoD, grids, dt_fine, n_fine, dts, extra in PLAN:
        if args.only and name not in args.only:
            continue
        tag = "_" + name
        try:
            #  A run the dx and the dt orders both ask for (the same grid and
            #  step) is made once.
            done = {}

            def once(n, dt):
                if (n, dt) not in done:
                    done[(n, dt)] = run_mms(args.workdir, exe, model, twoD, n, dt, extra, tag)
                return done[(n, dt)]
            rx = [once(n, dt_fine) for n in grids]
            ex = [r[0] for r in rx]
            #  The dt errors are measured against a run at a quarter of the
            #  smaller step on the same grid, so that the grid error, the same
            #  in all of them, drops out.
            rt = [once(n_fine, dt) for dt in dts]
            rr = once(n_fine, dts[1] / 4)
            ht, ref = [r[1] for r in rt], rr[1]
        except RuntimeError as e:
            print("FAIL %-6s %s" % (name, e))
            bad += 1
            continue
        px = math.log(ex[0] / ex[1]) / math.log(grids[1] / grids[0])
        et = [max(abs(a - b) for a, b in zip(h, ref)) / AMP for h in ht]
        pt = math.log(et[0] / et[1]) / math.log(dts[0] / dts[1]) if et[1] > 0 else float("inf")
        ok = px >= args.lo and pt >= args.lo
        #  With more than one quantity measured, the order of each is shown.
        #  In the cases that check the ConcZ velocity each has to converge on
        #  its own: a first-order error in one of them can hide for a while
        #  behind a larger second-order one in another (the combined order
        #  of cz_conc_2d with the bottom viscosity of the y profile taken from
        #  the wrong neighbor is 1.56, its velocity alone 1.28). So do h
        #  and c in wm_all_1d: with the periodic seam face of the
        #  concentration equation given the velocity of the far face (which
        #  for the centrifugal force, whose arm x - x0 is not periodic, is
        #  not the same) c alone converged at 1.24 on 32 -> 64 and 0.92 on
        #  64 -> 128 while the combined order passed. wm_all_1d is on 64 and
        #  128 cells: on 32
        #  the third-order dissipation of the Kawamura-Kuwahara convection
        #  of c, of the other sign, still takes a quarter off the error, and
        #  c converges at 1.67 (1.85 on 64 -> 128, 1.92 on 128 -> 256).
        #  wm_all_2d cannot afford grids that fine, and with the parameters
        #  of 1-D c alone converged at 1.41 on 16 -> 32 and 1.73 on 32 -> 64
        #  (h 1.79); with its own diffusivity and relief (full_terms) each
        #  of h and c is required to converge on 24 -> 48.
        each = ""
        if len(rx[0][2]) > 1:
            po = {k: math.log(rx[0][2][k] / rx[1][2][k]) / math.log(grids[1] / grids[0]) for k in rx[0][2]}
            if "u" in po or name in ("wm_all_1d", "wm_all_2d", "wm_all_slip_1d"):
                ok = ok and min(po.values()) >= args.lo
            each = "  [%s]" % ", ".join("%s %.2e -> %.2e: %.2f" % (k, rx[0][2][k], rx[1][2][k], po[k])
                                        for k in sorted(po))
        if (extra and "_Strat" in extra) or args.each:
            #  cz_strat_*: each of h, c and u on its own in dt as well
            de = dt_each([r[3] for r in rt], rr[3])
            pd = {k: (math.log(de[0][k] / de[1][k]) / math.log(dts[0] / dts[1])
                      if de[0][k] > 0 and de[1][k] > 0 else float("inf")) for k in de[0]}
            if extra and "_Strat" in extra:
                ok = ok and min(pd.values()) >= args.lo
            each += "  dt [%s]" % ", ".join("%s %.2e -> %.2e: %.2f" % (k, de[0][k], de[1][k], pd[k])
                                            for k in sorted(pd))
        bad += not ok
        print("%s %-9s error %.2e -> %.2e (order %.2f in dx), dt error %.2e -> %.2e (order %.2f in dt)%s"
              % ("ok  " if ok else "FAIL", name, ex[0], ex[1], px, et[0], et[1], pt, each))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
