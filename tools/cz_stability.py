#!/usr/bin/env python3
"""Short-wave stability of the ConcZ model about the stratified manufactured
state of tools/mms.py cz_strat_1d and cz_strat_2d.

    python3 tools/cz_stability.py --check               # cz_strat_1d
    python3 tools/cz_stability.py --twod --check        # cz_strat_2d
    python3 tools/cz_stability.py --diffusivity 3e-10   # the state at a lower D

A frozen-coefficient (WKB) analysis: at each x (and y) and t of the run the
state is taken as uniform in the plane, and a perturbation
exp(i k.x + lam t) of the film height and of the concentration at every
level through the film is linearized in the lubrication velocity profile
(viscosity mu(c), the Marangoni stress of the surface value, the capillary
pressure and the gravities of the case), the vertical velocity of the zeta
grid that carries the stratification, and the diffusion; finite differences
through the film. In 2-D the wave vector k takes every direction (12 of
them over a half turn). Printed is the largest Re(lam) over x, y, t and the
direction for each |k|, in 1/s: positive means that waves of that length
grow in the model itself, whatever the discretization (ch-verification,
"Choosing the convection scheme"). The wavenumbers are those of the grids
of the case and of one twice as fine, from the longest wave the analysis
means anything for (a quarter of the state's wavelength; longer waves are
the scale of the state itself, where the coefficients cannot be frozen; there
the rates are up to 5 /s for any diffusivity up to 3e-9) to the shortest one the finest grid
holds (two cells). With --check the exit status is 1 if any of them grows.

At 1.5e-9 m^2/s, the diffusivity of mms.py, the largest rate is -7.8 /s in
1-D and -3.5 /s in 2-D; at 3e-10, +710 /s at waves of 2 to 6 cells of 256
(e^14 over the run), and 2-D needs at least about 1.25e-9.

Needs NumPy; not part of ctest.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import math
import sys

import numpy as np

LX, LY = 0.25e-3, 0.375e-3       # wavelengths of the manufactured state, m
KX, KY = 2 * math.pi / LX, 2 * math.pi / LY
MEAN, AMP, OMEGA = 5e-6, 1.5e-6, 60.0
RHO, GZ = 1000.0, 3.0e4
SIG0 = 30e-3
MU0, D_EXP = 0.1, -2.0


class State:
    """The manufactured state of cz_strat_1d (twod False) or cz_strat_2d."""
    def __init__(self, twod):
        self.ky = KY if twod else 0.0
        self.g = np.array([6.0e3, 4.0e3 if twod else 0.0])

    def h(self, x, y, t):
        return MEAN + AMP * math.cos(KX * x - OMEGA * t) * math.cos(self.ky * y)

    def c(self, x, y, zeta, t):
        cy = math.cos(self.ky * y)
        return (0.5 + 0.2 * math.cos(KX * x - 40 * t + 0.7) * cy
                + 0.1 * math.cos(KX * x - 30 * t) * cy * math.cos(math.pi * zeta))

    def grad(self, f, x, y, eps=1e-9):
        return np.array([(f(x + eps, y) - f(x - eps, y)) / (2 * eps),
                         (f(x, y + eps) - f(x, y - eps)) / (2 * eps)])


def growth(st, x, y, t, ks, dirs, a, n=32):
    """The largest Re(lam) at each |k| over the directions, about the state
    at (x, y, t)."""
    z = (np.arange(n) + 0.5) / n
    dz = 1.0 / n
    h = st.h(x, y, t)
    hg = st.grad(lambda p, q: st.h(p, q, t), x, y)
    c0 = np.array([st.c(x, y, s, t) for s in z])
    c0g = np.array([st.grad(lambda p, q: st.c(p, q, s, t), x, y) for s in z])      # (n, 2)
    c0z = np.array([(st.c(x, y, s + 1e-6, t) - st.c(x, y, s - 1e-6, t)) / 2e-6 for s in z])
    csg = st.grad(lambda p, q: st.c(p, q, 1.0, t), x, y)
    mu = MU0 * np.exp(a.mu_exponent * (c0 - 0.5))
    dif = a.diffusivity * np.exp(D_EXP * (c0 - 0.5))
    sig = SIG0 + a.sigma_c * (st.c(x, y, 1.0, t) - 0.5)
    g0 = RHO * st.g - RHO * GZ * hg                # -grad p of the state
    zc = z * h
    du0 = (np.outer(h - zc, g0) + a.sigma_c * csg) / mu[:, None]    # du/dz, (n, 2)
    u0 = np.cumsum(du0 * dz * h, axis=0) - 0.5 * du0 * dz * h
    us = u0[-1] + 0.5 * du0[-1] * dz * h
    integ = np.tril(np.ones((n, n))) * dz - np.eye(n) * 0.5 * dz   # int_0^zeta
    lap = (np.diag(np.full(n - 1, 1.0), -1) + np.diag(np.full(n - 1, 1.0), 1)
           - np.diag(np.r_[1.0, np.full(n - 2, 2.0), 1.0])) / (dz * dz)   # no flux at both ends
    ws = np.zeros(n); ws[-1] = 1.5; ws[-2] = -0.5                  # surface value
    m = n + 1
    bm = np.eye(m, dtype=complex)
    bm[1:, 0] = -z * c0z / h
    out = []
    for k in ks:
        best = -math.inf
        for e in dirs:
            kv = k * e
            # du(zeta) = uh dh + uc dc, each of the two components
            dgdh = -1j * kv * (sig * k * k + RHO * GZ)
            uh = [h * (integ @ ((dgdh[d] * (h - zc) + g0[d]) / mu)) for d in range(2)]
            uc = [h * (integ @ (np.outer(1.0 / mu, a.sigma_c * 1j * kv[d] * ws)
                                - np.diag(du0[:, d] * a.mu_exponent))) for d in range(2)]
            qh = sum(kv[d] * (h * np.sum(uh[d]) * dz + us[d]) for d in range(2))
            qc = sum(kv[d] * h * np.sum(uc[d], axis=0) * dz for d in range(2))
            # zeta flux: dW = -i k.[h int du + dh int u0] - zeta lam dh
            wh = -1j * sum(kv[d] * (h * (integ @ uh[d]) + integ @ u0[:, d]) for d in range(2))
            wc = -1j * h * sum(kv[d] * (integ @ uc[d]) for d in range(2))
            ku0 = u0 @ kv
            am = np.zeros((m, m), complex)
            am[0, 0] = -1j * qh
            am[0, 1:] = -1j * qc
            am[1:, 1:] = (-sum(uc[d] * c0g[:, d][:, None] for d in range(2)) - wc * (c0z / h)[:, None]
                          + (dif / (h * h))[:, None] * lap
                          + np.diag(-1j * ku0 - dif * k * k))
            am[1:, 0] = -sum(uh[d] * c0g[:, d] for d in range(2)) - wh * c0z / h
            lam = np.linalg.eigvals(np.linalg.solve(bm, am))
            best = max(best, lam.real.max())
        out.append(best)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--twod", action="store_true", help="cz_strat_2d (the default is cz_strat_1d)")
    ap.add_argument("--diffusivity", type=float, default=1.5e-9, help="Fluid.Diffusivity.BaseValue, m^2/s (that of mms.py)")
    ap.add_argument("--sigma-c", type=float, default=0.004, help="Fluid.SurfaceTension.Coef, N/m")
    ap.add_argument("--mu-exponent", type=float, default=3.0, help="Fluid.Viscosity.Exponent")
    ap.add_argument("--cells", type=int, help="the finest grid, cells per wavelength in x"
                    " (the default is twice the finest of mms.py: 256 in 1-D, 64 in 2-D)")
    ap.add_argument("--nk", type=int, default=16, help="wavenumbers, spaced geometrically")
    ap.add_argument("--check", action="store_true", help="exit 1 if any wavenumber grows")
    a = ap.parse_args()
    st = State(a.twod)
    cells = a.cells or (64 if a.twod else 256)
    kmin = 4 * KX                                   # a quarter of the state's wavelength
    kmax = math.pi / (LX / cells)                   # two cells of the finest grid
    ks = list(np.geomspace(kmin, kmax, a.nk))
    if a.twod:
        dirs = [np.array([math.cos(p), math.sin(p)]) for p in np.linspace(0, math.pi, 12, endpoint=False)]
        pts = [(-LX / 2 + i * LX / 8, -LY / 2 + j * LY / 8) for i in range(8) for j in range(8)]
        ts = np.linspace(0, 0.02, 5)
    else:
        dirs = [np.array([1.0, 0.0])]
        pts = [(-LX / 2 + i * LX / 16, 0.0) for i in range(16)]
        ts = np.linspace(0, 0.02, 9)
    best = [-math.inf] * len(ks)
    for t in ts:
        for x, y in pts:
            r = growth(st, x, y, t, ks, dirs, a)
            best = [max(p, q) for p, q in zip(best, r)]
    print("cz_strat_%s, D = %g m^2/s, sigma_c = %g N/m, mu exponent %g"
          % ("2d" if a.twod else "1d", a.diffusivity, a.sigma_c, a.mu_exponent))
    print("%12s %14s %16s" % ("k [1/m]", "cells of %d" % cells, "max Re(lam) 1/s"))
    for k, b in zip(ks, best):
        print("%12.3e %14.2f %16.1f" % (k, 2 * math.pi / k / (LX / cells), b))
    worst = max(best)
    print("largest growth %.1f /s (e^%.2f over the run of 0.02 s)" % (worst, worst * 0.02))
    return 1 if a.check and worst > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
