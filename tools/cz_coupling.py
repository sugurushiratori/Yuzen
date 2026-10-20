#!/usr/bin/env python3
"""The gain of the coupling passes of the ConcZ model, per Fourier mode, on a
flat film with a concentration stratified through it.

    python3 tools/cz_coupling.py --dt 5e-3 --cells 480 --length 1e-3 --dc 0.2
    python3 tools/cz_coupling.py --dt 1e-3 --sigma-c -10e-3 --nz 16

With the Marangoni term on, every step of ConcZ solves the film and the
concentration again, one stage each, at the midpoint concentration c_m the
last pass found (CZ_HeightEquation_2D::AfterTimeStepTrial). Linearized about
a flat film of thickness h at rest, with c0(zeta) through it, a perturbation
exp(i k x) of the film eta and of the concentration gamma(zeta) is moved by
the plain pass (no relaxation) as

    gamma_m  ->  b + v gamma_m,s ,     lambda = v_s,

where gamma_m,s is its surface value: the film sees c_m only through the
Marangoni stress (and the surface tension) of the surface, and the
concentration only through the zeta velocity that the Marangoni flow and the
film's response to it make, which carries the stratification c0_zeta. The
pass is therefore rank one per mode, and lambda is its only non-zero
eigenvalue: the plain pass converges iff |lambda| < 1 for every mode. Its
fixed point is the implicit midpoint (Crank--Nicolson) step of the linear
system, and

    1 - lambda = det(I - dt/2 A) / [(1 + a) det(B)],

with A the linearized right-hand side, 1 + a the film's own
Crank--Nicolson factor (a = dt sigma h^3 k^4 / (6 mu)) and B that of the
diffusion through the film, both positive. So lambda > 1 only for a mode
whose linearized growth rate exceeds 2/dt, where the step itself turns a
growing perturbation into one that changes sign.

The differences are those of the solver: k^2 is the film's two-point symbol
(2 - 2 cos k dx)/dx^2, and the velocity's surface gradient, a fourth-order
face difference, is m_k = (27 sin(k dx/2) - sin(3 k dx/2)) / (24 sin(k dx/2))
times it (7/6 on the shortest wave). A viscosity that depends on c and the
flows of a film that is not at rest are left out. The ratio of successive
plain passes of a run (YUZEN_DEBUG_COUPLING) approaches the largest |lambda|
of the waves that hold the residual (ch-verification, "The coupling passes
of CZ").

Printed for a selection of the wavelengths the grid holds (2 cells to the
domain): lambda, the largest real growth rate of A and 1 - det(...) above as
a check; then the largest and smallest lambda. The cross-film profile is
linear (--dc is c at the surface less c at the substrate), with --evap-pe the
steady boundary layer exp(Pe zeta) of an evaporating film, or with --profile
the values of a column of a field file. Needs NumPy; not part of ctest.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import math
import sys

import numpy as np


def gain(k2, a, m=1.0):
    """lambda, the largest real growth rate and 1 - det(...) for the
    discrete wavenumber squared k2 (1/m^2) of the film's two-point
    differences; m is the velocity's surface gradient over the film's."""
    n, h, dt = a.nz, a.h, a.dt
    z = (np.arange(n) + 0.5) / n
    dz = 1.0 / n
    # c0_zeta at the cell centers, per unit zeta
    if a.profile:
        c = np.array([float(v) for v in a.profile.split(",")])
        n = len(c)
        z = (np.arange(n) + 0.5) / n
        dz = 1.0 / n
        cz = np.gradient(c, dz)
    elif a.evap_pe:
        pe = a.evap_pe
        cz = a.dc * pe * np.exp(pe * (z - 1.0)) / (1.0 - math.exp(-pe))
    else:
        cz = np.full(n, a.dc)
    lap = (np.diag(np.full(n - 1, 1.0), -1) + np.diag(np.full(n - 1, 1.0), 1)
           - np.diag(np.r_[1.0, np.full(n - 2, 2.0), 1.0])) / (dz * dz)   # no flux at both ends
    es = np.zeros(n)
    es[-1] = 1.0                              # the surface value: that of the top cell
    k4 = k2 * k2
    mu, s0, sc, dif = a.mu, a.sigma, a.sigma_c, a.diffusivity
    # zeta velocity of the perturbation (1/s): eta part and surface-c part
    om_eta = s0 * k4 * h * h / mu * z * (1 - z) * (2 - z) / 6.0
    om_c = -m * sc * k2 * h / (2 * mu) * z * (1 - z)
    lc = dif / (h * h) * lap - dif * k2 * np.eye(n)        # diffusion of gamma
    # film: eta_t = -fe eta + fc gamma_s
    fe = s0 * h ** 3 * k4 / (3 * mu)
    fc = sc * k2 * h * h / (2 * mu)
    # the pass: the film by Crank--Nicolson with the Marangoni stress of c_m,
    # the velocity from the film halfway, the concentration by Crank--Nicolson
    # with that velocity; the midpoint of the concentration it reaches
    aa = 0.5 * dt * fe
    deta = 0.5 * dt * fc / (1 + aa)          # d eta_mid / d gamma_m,s
    dom = om_eta * deta + om_c                # d Omega / d gamma_m,s
    bm = np.eye(n) - 0.5 * dt * lc
    v = -0.5 * dt * np.linalg.solve(bm, cz * dom)
    lam = v @ es
    # the linear system itself, y = (eta, gamma): its growth rates
    am = np.zeros((n + 1, n + 1))
    am[0, 0] = -fe
    am[0, 1:] = fc * es
    am[1:, 0] = -cz * om_eta
    am[1:, 1:] = lc - np.outer(cz * om_c, es)
    ev = np.linalg.eigvals(am)
    grow = max(ev.real)
    # 1 - lambda from the determinants, as a check of the two above
    det = np.linalg.det(np.eye(n + 1) - 0.5 * dt * am) / ((1 + aa) * np.linalg.det(bm))
    return lam, grow, 1.0 - det


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dt", type=float, required=True, help="time step, s")
    ap.add_argument("--h", type=float, default=10e-6, help="film thickness, m")
    ap.add_argument("--length", type=float, default=1e-3, help="domain length, m")
    ap.add_argument("--cells", type=int, default=128, help="cells over the domain")
    ap.add_argument("--nz", type=int, default=11, help="cells through the film (Grid.Nz)")
    ap.add_argument("--mu", type=float, default=0.1, help="viscosity, Pa s")
    ap.add_argument("--sigma", type=float, default=30e-3, help="surface tension, N/m")
    ap.add_argument("--sigma-c", type=float, default=10e-3, help="Fluid.SurfaceTension.Coef, N/m")
    ap.add_argument("--diffusivity", type=float, default=1e-9, help="m^2/s")
    ap.add_argument("--dc", type=float, default=0.1, help="c at the surface less c at the substrate")
    ap.add_argument("--evap-pe", type=float, default=0.0,
                    help="a profile exp(Pe zeta) of that size instead of a linear one")
    ap.add_argument("--profile", help="c at the cell centers through the film, substrate to surface, "
                    "comma-separated (instead of --dc; sets the number of cells)")
    ap.add_argument("--quiet", action="store_true", help="the extremes only")
    a = ap.parse_args()
    dx = a.length / a.cells
    rows = []
    for m in range(1, a.cells // 2 + 1):
        k = 2 * math.pi * m / a.length
        k2 = (2 - 2 * math.cos(k * dx)) / (dx * dx)
        th = k * dx
        m4 = (27 * math.sin(th / 2) - math.sin(1.5 * th)) / (24 * math.sin(th / 2))
        lam, grow, lam_det = gain(k2, a, m4)
        rows.append((a.cells / m, lam, grow, lam_det))
    if not a.quiet:
        print("%12s %12s %14s %12s" % ("cells/wave", "lambda", "growth 1/s", "check"))
        for r in rows[::max(1, len(rows) // 16)] + [rows[-1]]:
            print("%12.2f %12.4g %14.4g %12.4g" % r)
    lmax = max(rows, key=lambda r: r[1])
    lmin = min(rows, key=lambda r: r[1])
    gmax = max(r[2] for r in rows)
    print("lambda in [%.4g, %.4g] (at %.3g and %.3g cells a wave); largest growth %.4g /s, dt/2 times it %.4g"
          % (lmin[1], lmax[1], lmin[0], lmax[0], gmax, 0.5 * a.dt * gmax))
    return 0


if __name__ == "__main__":
    sys.exit(main())
