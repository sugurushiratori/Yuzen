#!/usr/bin/env python3
"""Generate yuzen/src/model/ManufacturedSourceTerms.inc with SymPy.

    python3 tools/mms_generate.py            # rewrite the file
    python3 tools/mms_generate.py --check    # fail if it is out of date

The method of manufactured solutions needs h_t + div Q and the corresponding
residual of the concentration equation for a film, a concentration and a
substrate given in closed form. With every term of the well-mixed model in
play -- capillary pressure over a substrate, both gravities, the centrifugal
term, the disjoining pressure, the Marangoni flux and properties that depend
on the concentration -- the fourth derivatives are too long to write out by
hand without a mistake, and a mistake in the source would look exactly like
a mistake in the solver. So they are derived here, symbolically, from the
equations as the manual states them (theory.md, section 1), and written out
as C++. SymPy is needed only to regenerate the file; the build and the tests
use the generated file as it stands.

Everything is dimensional (SI). The equations:

    h_t + div Q = S_h - E
    Q = h^3/(3 mu) [ cap sigma grad lap(h + s) - rho g grad(h + s)
                     + P grad Pi(h) + rho g_par + rho Om^2 (r - axis) ]
        + ms h^2/(2 mu) sigma_c grad c

(and with a Navier slip length lam, MmsHeightResidualSlip and
MmsConcResidualSlip: h^3/3 + lam h^2 for h^3/3 and h^2/2 + lam h for h^2/2)
    c_t + conv (Q/h) . grad c - diff (1/h) div(D h grad c) = S_c + E c / h

    sigma = sig0 + ms sigma_c (c - c0s),   mu = mu0 exp(mue (c - c0m)),
    D = D0 exp(De (c - c0D)),              Pi = (hs/h)^pn - (hs/h)^pm

and, for ConcZ, whose concentration c(x, y, zeta, t) lives through the
film (zeta = (z - s)/h), see build_cz(). There the manufactured c also varies
through the thickness, the viscosity mu(c(zeta)) with it, and neither the
flux nor the velocity profile has a closed form: what is generated is the
integrand of each depth integral, closed-form in (x, y, s, t), and
ManufacturedSource.cpp integrates them by Gauss-Legendre quadrature.

The generated functions return h_t + div Q and the left-hand side of the
concentration equation without its evaporation term; the evaporation is
added in ManufacturedSource.cpp, where it needs no derivatives. In 1-D the
same expressions are used with every y wave number set to 0 and y = 0: the
y-derivative of the centrifugal flux then leaves rho Om^2 h^3/(3 mu), which
is the sideways outflow the 1-D model accounts for.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "yuzen", "src", "model", "ManufacturedSourceTerms.inc")

PARAMS = ("M A w kx ky CM CA cw ckx cky cph CB bw bkx bky SB SA skx sky "
          "rho sig0 sigc c0s mu0 mue c0m D0 De c0D g gx gy Om2 ax ay "
          "P hs pn pm cap ms conv diff lam").split()


def build(slip=False):
    """The well-mixed model; slip: with the Navier slip length lam [m] at
    the substrate (mobility h^3/3 + lam h^2, shear mobility h^2/2 + lam h).
    The no-slip expressions are generated without lam, so that they stay
    what they were."""
    x, y, t = sp.symbols("x y t", real=True)
    p = {n: sp.Symbol("p." + n, real=True) for n in PARAMS}
    h = p["M"] + p["A"] * sp.cos(p["kx"] * x - p["w"] * t) * sp.cos(p["ky"] * y)
    c = p["CM"] + p["CA"] * sp.cos(p["ckx"] * x - p["cw"] * t + p["cph"]) * sp.cos(p["cky"] * y)
    s = p["SB"] + p["SA"] * sp.cos(p["skx"] * x) * sp.cos(p["sky"] * y)

    sig = p["sig0"] + p["ms"] * p["sigc"] * (c - p["c0s"])
    mu = p["mu0"] * sp.exp(p["mue"] * (c - p["c0m"]))
    D = p["D0"] * sp.exp(p["De"] * (c - p["c0D"]))
    Pi = (p["hs"] / h) ** p["pn"] - (p["hs"] / h) ** p["pm"]

    f = h + s
    lap = sp.diff(f, x, 2) + sp.diff(f, y, 2)
    mob = (h ** 3 / 3 + p["lam"] * h ** 2) / mu if slip else h ** 3 / (3 * mu)
    bx = (p["cap"] * sig * sp.diff(lap, x) - p["rho"] * p["g"] * sp.diff(f, x)
          + p["P"] * sp.diff(Pi, x) + p["rho"] * p["gx"] + p["rho"] * p["Om2"] * (x - p["ax"]))
    by = (p["cap"] * sig * sp.diff(lap, y) - p["rho"] * p["g"] * sp.diff(f, y)
          + p["P"] * sp.diff(Pi, y) + p["rho"] * p["gy"] + p["rho"] * p["Om2"] * (y - p["ay"]))
    mar = (p["ms"] * (h ** 2 / 2 + p["lam"] * h) / mu * p["sigc"] if slip
           else p["ms"] * h ** 2 / (2 * mu) * p["sigc"])
    Qx = mob * bx + mar * sp.diff(c, x)
    Qy = mob * by + mar * sp.diff(c, y)

    Sh = sp.diff(h, t) + sp.diff(Qx, x) + sp.diff(Qy, y)
    #  In 1-D the flux in y is not part of the velocity (it is the sideways
    #  outflow); with ky = cky = 0, grad c has no y component, so u_y . c_y
    #  vanishes by itself.
    ux, uy = Qx / h, Qy / h
    Sc = (sp.diff(c, t)
          + p["conv"] * (ux * sp.diff(c, x) + uy * sp.diff(c, y))
          - p["diff"] / h * (sp.diff(D * h * sp.diff(c, x), x) + sp.diff(D * h * sp.diff(c, y), y)))

    return x, y, t, Sh, Sc



def build_cz():
    """ConcZ. The concentration lives on (x, y, zeta), zeta = (z - s)/h in
    [0, 1], and the manufactured one varies through the thickness:

        c = CM + CA cos(ckx x - cw t + cph) cos(cky y)
               + CB cos(bkx x - bw t) cos(bky y) cos(pi zeta),

    with dc/dzeta = 0 at the substrate and at the surface -- the conditions
    ConcZ imposes there without evaporation. The surface value c_s = c(zeta
    = 1) sets the surface tension and drives the Marangoni stress; the
    viscosity is mu(c) at each height. The lubrication velocity is then

        u(zeta) = int_0^zeta K(s) ds,  K(s) = [h^2 (1 - s) b + h tau] / mu(c(s)),
        b = cap sigma(c_s) grad lap(h + s) - rho g grad(h + s) + rho g_par
            + rho Om^2 (r - axis),       tau = ms sigma_c grad c_s,

    and the film-height equation carries its depth integral,

        Q = h int_0^1 u dzeta = int_0^1 (1 - s) h K(s) ds,
        h_t + div Q = h_t + int_0^1 (1 - s) Dv(s) ds,  Dv = div(h K)

    (the divergence at fixed s). The concentration equation, as
    CZ_ConcEquation_2D writes it on the moving zeta grid, is

        c_t + u c_x + v c_y + Omega c_zeta - L c = S_c(zeta),

    all derivatives at fixed zeta, with the zeta velocity of the mapping

        Omega = -[div(h int_0^zeta u) + zeta h_t] / h
              = -[int_0^zeta (zeta - s) Dv(s) ds + zeta h_t] / h,

    and the diffusion at fixed z mapped onto zeta,

        L c = (d_x - A d_zeta)[D (c_x - A c_zeta)] + (the same in y)
              + d_zeta(D c_zeta) / h^2,     A = (zeta h_x + s_x) / h.

    Returned: the integrands Kx, Ky and Dv at s, h_t, and the closed-form
    part of the concentration residual, c_t - L c - zeta h_t c_zeta / h,
    with c_x, c_y and c_zeta for the rest."""
    x, y, t = sp.symbols("x y t", real=True)
    zeta, sv = sp.symbols("zeta s", real=True)
    p = {n: sp.Symbol("p." + n, real=True) for n in PARAMS}
    h = p["M"] + p["A"] * sp.cos(p["kx"] * x - p["w"] * t) * sp.cos(p["ky"] * y)
    sub = p["SB"] + p["SA"] * sp.cos(p["skx"] * x) * sp.cos(p["sky"] * y)

    def conc(z):
        return (p["CM"] + p["CA"] * sp.cos(p["ckx"] * x - p["cw"] * t + p["cph"]) * sp.cos(p["cky"] * y)
                + p["CB"] * sp.cos(p["bkx"] * x - p["bw"] * t) * sp.cos(p["bky"] * y) * sp.cos(sp.pi * z))
    c, cs = conc(zeta), conc(1)
    sig = p["sig0"] + p["ms"] * p["sigc"] * (cs - p["c0s"])
    m = sp.exp(-p["mue"] * (conc(sv) - p["c0m"])) / p["mu0"]
    D = p["D0"] * sp.exp(p["De"] * (c - p["c0D"]))

    f = h + sub
    lap = sp.diff(f, x, 2) + sp.diff(f, y, 2)
    bx = (p["cap"] * sig * sp.diff(lap, x) - p["rho"] * p["g"] * sp.diff(f, x)
          + p["rho"] * p["gx"] + p["rho"] * p["Om2"] * (x - p["ax"]))
    by = (p["cap"] * sig * sp.diff(lap, y) - p["rho"] * p["g"] * sp.diff(f, y)
          + p["rho"] * p["gy"] + p["rho"] * p["Om2"] * (y - p["ay"]))
    Kx = (h ** 2 * (1 - sv) * bx + h * p["ms"] * p["sigc"] * sp.diff(cs, x)) * m
    Ky = (h ** 2 * (1 - sv) * by + h * p["ms"] * p["sigc"] * sp.diff(cs, y)) * m
    #  In 1-D (ky = cky = bky = 0, y = 0) d/dy of h Ky leaves the sideways
    #  centrifugal outflow rho Om^2 h^3 (1 - s) m, which ConcZ counts in the
    #  film height and in its vertical velocity.
    Dv = sp.diff(h * Kx, x) + sp.diff(h * Ky, y)
    Ht = sp.diff(h, t)

    A_x = (zeta * sp.diff(h, x) + sp.diff(sub, x)) / h
    A_y = (zeta * sp.diff(h, y) + sp.diff(sub, y)) / h
    cz = sp.diff(c, zeta)
    Fx = D * (sp.diff(c, x) - A_x * cz)
    Fy = D * (sp.diff(c, y) - A_y * cz)
    Lc = (sp.diff(Fx, x) - A_x * sp.diff(Fx, zeta) + sp.diff(Fy, y) - A_y * sp.diff(Fy, zeta)
          + sp.diff(D * cz, zeta) / h ** 2)
    Loc = sp.diff(c, t) - Lc - zeta * Ht / h * cz
    return {"Kx": Kx, "Ky": Ky, "Dv": Dv, "Ht": Ht, "Loc": Loc,
            "cx": sp.diff(c, x), "cy": sp.diff(c, y), "cz": cz}

def ccode(name, expr, args="double x, double y, double t"):
    reps, (red,) = sp.cse(expr, symbols=sp.numbered_symbols("q"), optimizations="basic")
    lines = ["inline double %s(const MmsParams &p, %s)" % (name, args), "{"]
    for sym, e in reps:
        lines.append("\tconst double %s = %s;" % (sym, sp.ccode(e)))
    lines.append("\treturn %s;" % sp.ccode(red))
    lines.append("}")
    return "\n".join(lines)


def text():
    x, y, t, Sh, Sc = build()
    _, _, _, ShS, ScS = build(slip=True)
    z = build_cz()
    # The head of the generated file, license header included. The SPDX tag
    # below is text written to that file, not a tag of this one, so `reuse
    # lint` is told to skip it: REUSE-IgnoreStart
    head = ("//  Generated by tools/mms_generate.py -- do not edit by hand.\n"
            "//  h_t + div Q and the concentration residual without evaporation of\n"
            "//  the well-mixed model; for ConcZ the integrands of the depth\n"
            "//  integrals and the closed-form part of the concentration residual\n"
            "//  (see that script and ManufacturedSource.h).\n"
            "//\n"
            "//  Copyright 2017-2026 Suguru Shiratori, Tokyo City University\n"
            "//  SPDX-License-Identifier: Apache-2.0\n")
    # REUSE-IgnoreEnd
    fs = "double x, double y, double s, double t"
    zs = "double x, double y, double zeta, double t"
    parts = [ccode("MmsHeightResidual", Sh), ccode("MmsConcResidual", Sc),
             ccode("MmsCZKx", z["Kx"], fs), ccode("MmsCZKy", z["Ky"], fs), ccode("MmsCZDv", z["Dv"], fs),
             ccode("MmsCZHt", z["Ht"]), ccode("MmsCZLocal", z["Loc"], zs),
             ccode("MmsCZcx", z["cx"], zs), ccode("MmsCZcy", z["cy"], zs), ccode("MmsCZcz", z["cz"], zs),
             ccode("MmsHeightResidualSlip", ShS), ccode("MmsConcResidualSlip", ScS)]
    return head + "\n" + "\n\n".join(parts) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    new = text()
    if args.check:
        old = open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        if old != new:
            print("%s is out of date: run tools/mms_generate.py" % OUT)
            return 1
        print("up to date")
        return 0
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(new)
    print("wrote %s (%d lines)" % (OUT, new.count("\n")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
