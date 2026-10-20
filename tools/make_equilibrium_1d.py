#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For validating the disjoining pressure: construct a 1D equilibrium droplet
shape and write it to an ASCII restart file.

--------------------------------------------------------------------
Equilibrium condition
--------------------------------------------------------------------
The film-thickness equation is

    dH/dT = -d/dX [ H^3/(3 mu) d/dX ( d2H/dX2 + Dp * Pi(H) ) ]

so zero flux = equilibrium means

    d2H/dX2 + Dp * Pi(H) = K        (K is a constant = the pressure)

    Pi(H) = (H*/H)^n - (H*/H)^m

Multiplying by H' and integrating once gives

    (1/2) H'^2 = K (H - Hf) - Dp [ w(H) - w(Hf) ],   w(H) = int_{H*}^{H} Pi ds

On the flat film Hf surrounding the droplet H' = H'' = 0, hence
K = Dp * Pi(Hf). At the apex Hmax we also have H' = 0, so

    K (Hmax - Hf) - Dp [ w(Hmax) - w(Hf) ] = 0

This is the relation linking Hf and Hmax. Given Hmax, Hf is found by
bisection.

--------------------------------------------------------------------
Contact angle
--------------------------------------------------------------------
In the limit Hmax >> H* the above gives (H')^2 -> As * theta_e^2, so the
physical slope eps*|H'| matches theta_e. The expression for Dp

    Dp = As (n-1)(m-1) / (2 H* (n-m)) theta_e^2

has been chosen so that this holds. This script actually builds the
equilibrium shape and reports the measured maximum slope, so the agreement
can be checked numerically.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""

import argparse
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp


# ------------------------------------------------------------------
#  Disjoining pressure and its integral
# ------------------------------------------------------------------
def make_pi(hstar, n, m):
    def Pi(H):
        r = hstar / H
        return r**n - r**m

    def dPi(H):
        r = hstar / H
        return -(n * r**n - m * r**m) / H

    def w(H):
        #  int_{H*}^{H} Pi ds  (n != 1, m != 1)
        a = hstar**n / (1.0 - n) * (H ** (1.0 - n) - hstar ** (1.0 - n))
        b = hstar**m / (1.0 - m) * (H ** (1.0 - m) - hstar ** (1.0 - m))
        return a - b

    return Pi, dPi, w


# ------------------------------------------------------------------
#  Determine the flat film thickness Hf from Hmax
# ------------------------------------------------------------------
def solve_hf(Hmax, hstar, n, m, Dp):
    Pi, dPi, w = make_pi(hstar, n, m)

    def G(Hf):
        K = Dp * Pi(Hf)
        return K * (Hmax - Hf) - Dp * (w(Hmax) - w(Hf))

    lo = hstar * (1.0 + 1e-12)
    #  Pi attains its minimum at H = H*(n/m)^(1/(n-m)); K decreases
    #  monotonically up to that point.
    hi = hstar * (n / m) ** (1.0 / (n - m))
    if G(lo) * G(hi) > 0:
        raise RuntimeError(
            "cannot bracket Hf (Hmax = %g); try a smaller Hmax." % Hmax
        )
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if G(lo) * G(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------------
#  Integration of the shape
#
#  The second-order ODE is solved inwards starting from the foot
#  (H = Hf + eta). At the foot H'' = -Dp*Pi'(Hf) * eta = lambda^2 * eta,
#  so the growing solution eta ~ exp(lambda X) is used as the initial
#  condition. Integrating in this direction is numerically stable, since
#  errors along the decaying solution die out on their own.
# ------------------------------------------------------------------
def integrate_profile(Hmax, hstar, n, m, Dp, eta0=1e-10, rtol=1e-12, atol=1e-16):
    Pi, dPi, w = make_pi(hstar, n, m)
    Hf = solve_hf(Hmax, hstar, n, m, Dp)
    K = Dp * Pi(Hf)

    lam2 = -Dp * dPi(Hf)
    if lam2 <= 0:
        raise RuntimeError("the foot is unstable (lambda^2 <= 0)")
    lam = math.sqrt(lam2)

    def rhs(X, y):
        return [y[1], K - Dp * Pi(y[0])]

    def apex(X, y):
        return y[1]

    apex.terminal = True
    apex.direction = -1

    y0 = [Hf + eta0, lam * eta0]
    #  Length needed for the exponential growth at the foot, plus the length
    #  of the outer parabolic cap. If it is not enough, double it and solve
    #  again.
    Xend = math.log(1.0 / eta0) / lam + 4.0 * math.sqrt(2.0 * Hmax / max(-K, 1e-300))
    sol = None
    for _ in range(12):
        sol = solve_ivp(
            rhs, (0.0, Xend), y0, events=apex, rtol=rtol, atol=atol,
            dense_output=True, max_step=0.2 / lam,
        )
        if sol.t_events[0].size:
            break
        Xend *= 2.0
    if sol is None or not sol.t_events[0].size:
        raise RuntimeError("the apex was not reached")

    Xapex = sol.t_events[0][0]
    Happex = sol.y_events[0][0][0]

    return sol, Xapex, Happex, Hf, K, lam


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="ASCII restart file to write")
    #  Grid (match [Grid] in the input file)
    ap.add_argument("--xmin", type=float, default=-0.3e-3)
    ap.add_argument("--xmax", type=float, default=+0.3e-3)
    ap.add_argument("--nx", type=int, default=600)
    #  Scales (match [Fluid] in the input file)
    ap.add_argument("--sigma", type=float, default=30.0e-3)
    ap.add_argument("--rho", type=float, default=1000.0)
    ap.add_argument("--h0", type=float, default=10e-6, help="Initial.SurfaceShape.TotalHeight")
    ap.add_argument("--length-scale", type=float, default=None,
                    help="LengthScale used when LengthScaleType='custom'")
    #  Disjoining pressure (match [Physics.DisjoiningPressure] in the input file)
    ap.add_argument("--hstar", type=float, default=0.2, help="PrecursorRatio")
    ap.add_argument("--theta", type=float, default=10.0, help="ContactAngle [deg]")
    ap.add_argument("--n", type=float, default=3.0)
    ap.add_argument("--m", type=float, default=2.0)
    #  Droplet
    ap.add_argument("--hmax", type=float, default=1.0, help="apex film thickness (dimensionless)")
    ap.add_argument("--conc", type=float, default=0.5)
    ap.add_argument("--sweep", action="store_true",
                    help="only print a table of the effective contact angle over a range of Hmax/H* (writes no file)")
    args = ap.parse_args()

    #  ---- Scales ----
    L0 = args.length_scale
    if L0 is None:
        L0 = math.sqrt(args.sigma / (args.rho * 9.80665))
    h0 = args.h0
    eps = h0 / L0
    As = eps ** (-2.0)

    hstar, n, m = args.hstar, args.n, args.m
    th = args.theta * math.pi / 180.0
    Dp = As * (n - 1.0) * (m - 1.0) / (2.0 * hstar * (n - m)) * th * th

    if args.sweep:
        Pi, dPi, w = make_pi(hstar, n, m)
        print("theta_e = %g deg,  n = %g, m = %g" % (args.theta, n, m))
        print("  Hmax/H*   Hf/H*    apparent angle[deg]  max slope[deg]  foot decay length/h* (physical)")
        for ratio in (2, 5, 10, 20, 50, 100, 200, 500, 1000):
            Hmax = hstar * ratio
            try:
                Hf = solve_hf(Hmax, hstar, n, m, Dp)
            except RuntimeError:
                continue
            K = Dp * Pi(Hf)
            Hs = np.linspace(Hf, Hmax, 200001)
            E = (K * (Hs - Hf) - Dp * (w(Hs) - w(Hf))).max()
            smax = eps * math.sqrt(max(2.0 * E, 0.0)) * 180.0 / math.pi
            cap = eps * math.sqrt(max(-2.0 * K * Hmax, 0.0)) * 180.0 / math.pi
            lam = math.sqrt(-Dp * dPi(Hf))
            print("  %7d  %7.4f  %13.4f  %13.4f  %12.2f"
                  % (ratio, Hf / hstar, cap, smax, (1.0 / lam) / (hstar * eps)))
        print("\napparent angle = slope obtained by extrapolating the outer parabola to H=0 "
              "(matches theta_e*(1-H*/2Hmax)).")
        print("max slope      = actual maximum slope of the shape. For small Hmax/H* it is "
              "considerably smaller than theta_e.")
        return

    #  ---- Shape ----
    sol, Xapex, Happex, Hf, K, lam = integrate_profile(args.hmax, hstar, n, m, Dp)

    #  ---- Grid ----
    #  Match the ordering written by OutputNeutral. C_Grid1D has two ghost
    #  layers on each side, and output/restart files hold nx()-2 = Nx+2 rows,
    #  ordered as
    #      Xmin, Nx cell centers, Xmax
    #  (the two end points are the boundaries themselves).
    nx = args.nx
    dx_dim = (args.xmax - args.xmin) / nx
    x_dim = np.empty(nx + 2)
    x_dim[0] = args.xmin
    x_dim[1:-1] = args.xmin + (np.arange(nx) + 0.5) * dx_dim
    x_dim[-1] = args.xmax
    nrow = nx + 2
    X = x_dim / L0                      # dimensionless
    r = np.abs(X)                       # symmetric shape with the apex at X=0

    #  Distance from the apex to the foot
    Xtail = Xapex                       # X=0 of sol is the foot, Xapex the apex
    H = np.full(nrow, Hf)
    inside = r < Xtail
    H[inside] = sol.sol(Xapex - r[inside])[0]

    if r.max() < Xtail:
        print("!! the droplet extends beyond the computational domain; reduce --hmax or widen the domain.",
              file=sys.stderr)

    #  ---- Diagnostics ----
    #  Maximum slope of the continuous solution (physical)
    Xs = np.linspace(0.0, Xapex, 20001)
    Hp = sol.sol(Xs)[1]
    slope_cont = eps * np.abs(Hp).max()
    #  Maximum slope after discretization (physical)
    slope_disc = eps * np.abs(np.diff(H) / (dx_dim / L0)).max()

    #  In the outer region (where the disjoining pressure is negligible) the
    #  shape is a parabola with H'' = K. The slope obtained by extrapolating
    #  it to H=0 corresponds to the "apparent contact angle" measured in
    #  experiments:
    #      theta_cap = eps sqrt(-2 K Hmax)
    #  For Hmax >> H*, theta_cap -> theta_e (1 - H*/(2 Hmax)).
    cap_cont = eps * math.sqrt(max(-2.0 * K * Happex, 0.0))
    cap_pred = th * (1.0 - hstar / (2.0 * Happex))

    vol = np.sum(H - Hf) * dx_dim / L0

    print("######## equilibrium droplet ########")
    print(" L0 [m]          = %.6e" % L0)
    print(" eps = h0/L0     = %.6e" % eps)
    print(" As              = %.6e" % As)
    print(" Dp              = %.6e" % Dp)
    print(" H*              = %.6f   (%.4f um)" % (hstar, hstar * h0 * 1e6))
    print(" Hf (flat film)  = %.9f   (%.6f um)" % (Hf, Hf * h0 * 1e6))
    print(" K  (pressure)   = %.6e" % K)
    print(" Hmax (apex)     = %.6f   (%.4f um)" % (Happex, Happex * h0 * 1e6))
    print(" droplet half-width = %.6e  (%.3f um)" % (Xapex, Xapex * L0 * 1e6))
    print(" excess volume      = %.6e" % vol)
    print(" grid points / half-width = %.1f" % (Xapex * L0 / dx_dim))
    print(" foot decay length        = %.3f um" % (L0 / lam * 1e6))
    print("-------- contact angle --------")
    d = 180.0 / math.pi
    print(" prescribed theta_e             = %.6f deg" % args.theta)
    print(" apparent angle (outer parabola)= %.6f deg" % (cap_cont * d))
    print("   asymptotic th(1-H*/2Hmax)    = %.6f deg" % (cap_pred * d))
    print(" max slope, continuous solution = %.6f deg" % (slope_cont * d))
    print(" max slope on the grid          = %.6f deg" % (slope_disc * d))
    print("####################################")

    #  ---- Output (same columns as OutputNeutral: xx HS Ux CC HH SS) ----
    #  With Output.Dimensioned = true, film thickness is in um and x in mm.
    with open(args.out, "w") as f:
        f.write("xx \t HS \t Ux \t CC \t HH \t SS\n")
        for i in range(nrow):
            hh = H[i] * h0 * 1e6
            #  The second difference is divided by dX^2, so dropping digits
            #  corrupts the pressure. Write the full double precision.
            f.write("%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n"
                    % (x_dim[i] * 1e3, hh, 0.0, args.conc, hh, 0.0))
    print("wrote %s (%d rows)" % (args.out, nrow))


if __name__ == "__main__":
    main()
