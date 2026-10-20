#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Measure the contact line position, the apparent contact angle and the
Cox-Voinov (Tanner) relation from a spreading calculation.

--------------------------------------------------------------------
What is measured
--------------------------------------------------------------------
The Neutral output (field*****.dat) is read in order, and at each time

  - a parabola is least-squares fitted to the outer part of the droplet
    (the region where disjoining pressure is not active)
  - the base radius R and the slope theta obtained by extrapolating that
    parabola to h = 0
      h = a + b x^2   ->   R = sqrt(-a/b),  theta = 2 sqrt(-a b)
    theta corresponds to the apparent contact angle measured in
    experiments.
  - U = dR/dt and Ca = mu U / sigma

are reported. Finally a straight line is fitted to Cox-Voinov

    theta^3 - theta_e^3 = 9 ln(L/l) * Ca

taking theta_e to be theta at the final time (the equilibrium value).
This relation holding is the substance of the so-called Tanner law
(R ~ t^(1/7)).

--------------------------------------------------------------------
Usage
--------------------------------------------------------------------
    python3 tools/measure_spreading.py runs/spread/spread --dt-out 1e-3

--dt-out must match Output.FieldData.TimeInterval in the input file.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""

import argparse
import glob
import math
import os

import numpy as np


def fit_cap(path, hcut_um):
    """Return (R [m], theta [rad], Hmax [m]) from a single field file."""
    d = np.loadtxt(path, skiprows=1)
    x = d[:, 0] * 1e-3          # mm -> m
    h = d[:, 4] * 1e-6          # um -> m
    sel = h > hcut_um * 1e-6
    if sel.sum() < 8:
        return None
    c = np.polyfit(x[sel], h[sel], 2)
    b, a = c[0], c[2]           # h = b x^2 + c1 x + a
    if not (b < 0 and a > 0):
        return None
    return math.sqrt(-a / b), 2.0 * math.sqrt(-a * b), h.max()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resultdir", help="directory containing field*****.dat")
    ap.add_argument("--dt-out", type=float, default=1e-3,
                    help="Output.FieldData.TimeInterval [s]")
    ap.add_argument("--hcut", type=float, default=1.0,
                    help="lower film thickness bound for the parabola fit [um]")
    ap.add_argument("--mu", type=float, default=0.1, help="viscosity [Pa s]")
    ap.add_argument("--sigma", type=float, default=30e-3,
                    help="surface tension [N/m]")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--tmax", type=float, default=1e30)
    ap.add_argument("--csv", default=None, help="write out the time series")
    args = ap.parse_args()

    fs = sorted(glob.glob(os.path.join(args.resultdir, "field*.dat")))
    if not fs:
        raise SystemExit("no field*****.dat found: " + args.resultdir)

    T, R, TH, HM = [], [], [], []
    for k, f in enumerate(fs):
        if k % args.stride:
            continue
        r = fit_cap(f, args.hcut)
        if r is None:
            continue
        T.append(k * args.dt_out)
        R.append(r[0]); TH.append(r[1]); HM.append(r[2])

    T = np.array(T); R = np.array(R); TH = np.array(TH); HM = np.array(HM)
    if len(T) < 5:
        raise SystemExit("not enough points")

    U = np.gradient(R, T)
    Ca = args.mu * U / args.sigma
    the = TH[-1]                      # equilibrium value
    deg = 180.0 / math.pi

    print("######## spreading ########")
    print(" files = %d,  t = %.4g .. %.4g s" % (len(T), T[0], T[-1]))
    print(" theta: %.4f -> %.4f deg   (final time taken as equilibrium)"
          % (TH[0] * deg, the * deg))
    print(" R    : %.3f -> %.3f um" % (R[0] * 1e6, R[-1] * 1e6))
    print(" volume proxy theta*R^2 : %.5g -> %.5g (constant means volume conserved)"
          % (TH[0] * R[0] ** 2, the * R[-1] ** 2))
    print()
    print("   t[s]        R[um]     theta[deg]    U[um/s]        Ca")
    n = max(1, len(T) // 15)
    for i in range(0, len(T), n):
        print("  %.5f  %10.3f  %11.4f  %11.3f  %11.3e"
              % (T[i], R[i] * 1e6, TH[i] * deg, U[i] * 1e6, Ca[i]))

    #  Cox-Voinov
    m = (U > 0) & (T > T[1]) & (T >= args.tmin) & (T <= args.tmax)
    #  Exclude the interval close enough to equilibrium, where
    #  theta^3 - theta_e^3 is buried in round-off error
    m &= (TH - the) > 0.02 * (TH[0] - the)
    y = TH ** 3 - the ** 3
    if m.sum() >= 5:
        A = np.polyfit(Ca[m], y[m], 1)
        rr = np.corrcoef(Ca[m], y[m])[0, 1]
        print()
        print("-------- Cox-Voinov --------")
        print(" theta^3 - theta_e^3 = %.4f * Ca %+.3e" % (A[0], A[1]))
        print("   => 9 ln(L/l) = %.3f  =>  ln(L/l) = %.3f  (L/l = %.1f)"
              % (A[0], A[0] / 9.0, math.exp(A[0] / 9.0)))
        print(" correlation coefficient = %.5f  (%d points, Ca = %.2e .. %.2e)"
              % (rr, m.sum(), Ca[m].min(), Ca[m].max()))
        print(" The closer the intercept is to 0, the better the choice of theta_e.")
    else:
        print("\nnot enough points to fit Cox-Voinov")

    if args.csv:
        np.savetxt(args.csv,
                   np.c_[T, R, TH, U, Ca, HM],
                   header="t[s] R[m] theta[rad] U[m/s] Ca Hmax[m]")
        print("\nwrote " + args.csv)


if __name__ == "__main__":
    main()
