#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Draw the figures for the shape-function reference of the Yuzen manual.

The formulas here are transcriptions of
yuzen/src/commonFunctions/ShapeFunctions.cpp and the combination rule of
GenerateShape.cpp.  They are kept literal on purpose: the point of these
figures is to show what the solver actually builds, so the definitions must
track the source, not an idealization of it.

    python3 docs/manual/figures/make_shapes.py

writes shapes-1d.pdf, shapes-createtype.pdf and film-geometry.pdf next to
itself.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))

# --- the Yuzen palette -----------------------------------------------------
INK      = "#22201D"
CELADON  = "#5E8C78"
SHADOW   = "#48705E"
AMBER    = "#C08A3E"
LAPIS    = "#2E4A7D"
PAPER    = "#F4F1EA"
GRAY     = "#9A968E"

plt.rcParams.update({
    "font.family":      "serif",
    "font.serif":       ["TeX Gyre Pagella", "Palatino", "DejaVu Serif"],
    "font.size":        8.0,
    "axes.titlesize":   8.5,
    "axes.labelsize":   8.0,
    "axes.edgecolor":   INK,
    "axes.linewidth":   0.6,
    "axes.labelcolor":  INK,
    "text.color":       INK,
    "xtick.color":      INK,
    "ytick.color":      INK,
    "xtick.labelsize":  7.0,
    "ytick.labelsize":  7.0,
    "xtick.direction":  "out",
    "ytick.direction":  "out",
    "legend.fontsize":  7.0,
    "legend.frameon":   False,
    "figure.facecolor": "white",
    "savefig.facecolor":"white",
    "savefig.bbox":     "tight",
    "savefig.pad_inches": 0.02,
    # text as outlines, so the SVG needs no font and converts faithfully
    "svg.fonttype":     "path",
})


def save(fig, stem):
    """Write both forms: PDF for the LaTeX build, SVG for the HTML build."""
    out = os.path.join(HERE, stem + ".pdf")
    fig.savefig(out)
    fig.savefig(os.path.join(HERE, stem + ".svg"))
    plt.close(fig)
    return out


# --- the shape functions, transcribed from ShapeFunctions.cpp --------------
def sigmoid(x, x0, eps):
    return 1.0 / (1.0 + np.exp(-(x - x0) / eps))


def rect_profile_base(p, center, width, slope):
    """RectProfile_base(): note the 0.1 factor applied to the slope."""
    return (sigmoid(p, center - 0.5 * width, slope * 0.1)
            * (1.0 - sigmoid(p, center + 0.5 * width, slope * 0.1)))


def rect_bump(x, base=0.0, amp=1.0, cx=0.0, wx=1.0, sx=0.1):
    return amp * rect_profile_base(x, cx, wx, sx) + base


def ellipse_bump(x, base=0.0, amp=1.0, cx=0.0, wx=0.5, slope=0.1):
    r = np.abs((x - cx) / wx)
    return amp * rect_profile_base(r, 0.0, 1.0, slope) + base


def flat_plane(x, amp=0.0, cx=0.0, angle=0.5):
    return (x - cx) * np.tan(angle) + amp


def mexican_hat(x, base=0.0, amp=1.0, cx=0.0, L=0.5, beta=1.0):
    u = ((x - cx) / L) ** 2
    return amp * (1.0 - u) * np.exp(-beta * u) + base


def smooth_convex(x, base=0.0, amp=1.0, cx=0.0, L=0.6):
    u = ((x - cx) / L) ** 2
    out = np.where(np.sqrt(u) > 1.0, base, amp * (1.0 - u) ** 2 + base)
    return out


def periodic_wave(x, base=0.0, amp=1.0, cx=0.0, L=1.0):
    return amp * np.cos(2.0 * np.pi * (x - cx) / L) + base


def edge_bump(x, base=0.0, amp=1.0, alpha=1.0, Lin=0.3, Lout=0.6, rs=None,
              location="Xmax", xmin=-1.0, xmax=1.0):
    # rs is not an input key: GenerateShape.cpp derives it from the domain edge,
    # rs = Xmax - Lout for an Xmax bump and Xmin + Lout for an Xmin one.
    if rs is None:
        rs = (xmax - Lout) if location == "Xmax" else (xmin + Lout)
    aa = amp - base
    eps = 0.1 * Lin
    pp = (x - rs) if location == "Xmax" else -(x - rs)
    ff = aa * (1.0 - (pp / Lin) ** 2) * np.exp(-1.0 / alpha * (pp / Lin) ** 2) + base
    gg = amp * (1.0 - (pp / Lout) ** 2)      # ShapeFunctions.cpp; holds for amp == base
    sig = sigmoid(pp, 0.0, eps)
    return (1.0 - sig) * ff + sig * gg


# --- figure 1: the seven shapes in 1-D -------------------------------------
def figure_shapes_1d():
    x = np.linspace(-1.0, 1.0, 2000)

    panels = [
        ("RectBump", rect_bump(x),
         r"$\mathrm{Width}=1.0$, $\mathrm{SlopeWidth}=0.1$"),
        ("EllipseBump", ellipse_bump(x),
         r"$\mathrm{Width}=0.5$, $\mathrm{SlopeWidth}=0.1$"),
        ("SmoothConvex", smooth_convex(x),
         r"$L=0.6$;  exactly $\mathrm{Base}$ outside"),
        ("MexicanHat", mexican_hat(x),
         r"$L=0.5$, $\beta=1$;  goes below $\mathrm{Base}$"),
        ("PeriodicWave", periodic_wave(x),
         r"$\mathrm{Wavelength}=1.0$"),
        ("FlatPlane", flat_plane(x),
         r"$\mathrm{Angle}=0.5$ rad;  $\mathrm{Base}$ is unused"),
        ("EdgeBump", edge_bump(x, base=1.0, amp=1.0),
         r"$\mathrm{Base}=\mathrm{Height}=1$ (dashed: $\mathrm{Height}=1.3$), pinned at $x_{\max}$"),
    ]

    fig = plt.figure(figsize=(6.6, 5.9))
    gs = GridSpec(4, 2, figure=fig, hspace=0.95, wspace=0.22)

    for i, (name, y, note) in enumerate(panels):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        ax.axhline(0.0, color=GRAY, lw=0.5, ls=(0, (3, 3)), zorder=1)
        if name == "EdgeBump":
            ax.plot(x, edge_bump(x, base=1.0, amp=1.3), color=GRAY, lw=0.9,
                    ls=(0, (4, 2)), zorder=2)
        ax.plot(x, y, color=CELADON, lw=1.7, zorder=3, solid_capstyle="round")
        ax.fill_between(x, 0.0, y, color=CELADON, alpha=0.13, zorder=2)
        ax.set_title(name, loc="left", color=SHADOW, fontweight="bold", pad=15)
        ax.text(0.0, 1.015, note, transform=ax.transAxes, fontsize=6.4,
                color=INK, alpha=0.75, va="bottom")
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-1.25, 1.45)
        ax.set_yticks([-1, 0, 1])
        ax.set_xticks([-1, 0, 1])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    # legend panel in the last cell
    ax = fig.add_subplot(gs[3, 1])
    ax.axis("off")
    ax.text(0.0, 0.99,
            "Five of the seven are\n"
            r"   $\mathrm{Amp}\times f(x) + \mathrm{Base}$." "\n"
            "FlatPlane adds Amp to a tilted\n"
            "plane and ignores Base;\n"
            "EdgeBump holds Base inside\n"
            "and falls from Height to zero\n"
            "at the pinned edge.\n\n"
            r"Drawn at $\mathrm{Amp}=1$, $\mathrm{Base}=0$," "\n"
            r"$\mathrm{Center}=0$, EdgeBump at" "\n"
            r"$L_\mathrm{in}=0.3$, $L_\mathrm{out}=0.6$."
            "\n\n"
            "The axis is arbitrary here. In an\n"
            "input file the horizontal\n"
            "parameters are in meters, in the\n"
            "Grid.Xmin to Grid.Xmax\n"
            "coordinates, except\n"
            "EllipseBump.SlopeWidth, which\n"
            "is in normalized radius.",
            transform=ax.transAxes, va="top", fontsize=6.4, color=INK,
            linespacing=1.45)

    return save(fig, "shapes-1d")


# --- figure 2: how CreateType combines shapes ------------------------------
def figure_createtype():
    """The accumulator starts at 0 and each shape folds into it in Order,
    by Max (add), Min (remove) or += (stack).  See GenerateShape.cpp."""
    x = np.linspace(-1.0, 1.0, 3000)
    a = smooth_convex(x, amp=1.0, cx=-0.22, L=0.60)
    b = smooth_convex(x, amp=0.75, cx=0.22, L=0.50)

    combos = [
        ("add", np.maximum(np.maximum(0.0, a), b), CELADON,
         "the pointwise maximum of the two:\ntheir union"),
        ("stack", 0.0 + a + b, CELADON,
         "the sum: where they overlap,\nthe heights add"),
        ("remove", np.minimum(np.maximum(0.0, a), b), AMBER,
         "the pointwise minimum: their\nintersection, which is how you carve"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(6.6, 2.15), sharey=True)
    for ax, (name, y, color, note) in zip(axes, combos):
        ax.plot(x, a, color=GRAY, lw=0.9, ls=(0, (4, 2)), zorder=2)
        ax.plot(x, b, color=GRAY, lw=0.9, ls=(0, (4, 2)), zorder=2)
        ax.fill_between(x, 0.0, y, color=color, alpha=0.18, lw=0, zorder=3)
        ax.plot(x, y, color=color, lw=2.0, zorder=4, solid_capstyle="round")
        ax.axhline(0.0, color=GRAY, lw=0.5, ls=(0, (3, 3)), zorder=1)
        ax.set_title('CreateType = "%s"' % name, loc="left",
                     color=SHADOW, fontweight="bold", pad=4)
        ax.text(0.0, -0.20, note, transform=ax.transAxes, fontsize=6.6,
                color=INK, alpha=0.85, va="top", linespacing=1.5)
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-0.08, 1.85)
        ax.set_xticks([-1, 0, 1])
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_yticks([0, 1])
    fig.subplots_adjust(bottom=0.30, wspace=0.18)

    return save(fig, "shapes-createtype")


# --- figure 3: the film, the substrate and the precursor -------------------
def figure_film_geometry():
    x = np.linspace(-1.0, 1.0, 4000)

    # substrate: one broad low feature
    s = 0.20 * smooth_convex(x, amp=1.0, cx=0.35, L=0.55)
    # film: a droplet sitting on a thin precursor film
    hstar = 0.07
    h = hstar + 0.85 * smooth_convex(x, amp=1.0, cx=-0.30, L=0.45)

    fig, ax = plt.subplots(figsize=(6.6, 2.6))

    ax.fill_between(x, 0.0, s, color=INK, alpha=0.13, lw=0, zorder=2)
    ax.plot(x, s, color=INK, lw=1.4, zorder=4)
    ax.fill_between(x, s, s + h, color=CELADON, alpha=0.32, lw=0, zorder=3)
    ax.plot(x, s + h, color=SHADOW, lw=1.8, zorder=5)

    i = lambda v: int(np.argmin(np.abs(x - v)))

    ax.annotate("free surface  $h+s$   (column HS)",
                xy=(-0.30, (s + h)[i(-0.30)]),
                xytext=(-0.97, 1.30), fontsize=7.4, color=SHADOW,
                arrowprops=dict(arrowstyle="-", color=SHADOW, lw=0.7,
                                shrinkA=0, shrinkB=3))
    ax.annotate("substrate  $s$   (column SS)",
                xy=(0.45, s[i(0.45)]),
                xytext=(0.20, 0.92), fontsize=7.4, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.7,
                                shrinkA=0, shrinkB=3))

    # film thickness
    x0 = -0.30
    ax.annotate("", xy=(x0, s[i(x0)]), xytext=(x0, (s + h)[i(x0)]),
                arrowprops=dict(arrowstyle="<->", color=AMBER, lw=1.2))
    ax.text(x0 - 0.045, s[i(x0)] + 0.5 * h[i(x0)],
            "film thickness $h$\n(column HH)",
            fontsize=7.4, color=AMBER, va="center", ha="right",
            linespacing=1.4)

    # precursor film, drawn on a magnified inset-style callout
    xp = 0.80
    ax.annotate("", xy=(xp, s[i(xp)]), xytext=(xp, (s + h)[i(xp)]),
                arrowprops=dict(arrowstyle="<->", color=LAPIS, lw=1.2))
    ax.plot([xp, xp], [(s + h)[i(xp)], 0.52], color=LAPIS, lw=0.6,
            ls=(0, (2, 2)))
    ax.text(xp, 0.55, "precursor film $h_*$\n(never zero)", fontsize=7.4,
            color=LAPIS, ha="center", va="bottom", linespacing=1.4)

    # contact line
    xc = 0.145
    ax.plot([xc], [(s + h)[i(xc)]], marker="o", ms=3.4, color=AMBER, zorder=6)
    ax.annotate("contact line", xy=(xc, (s + h)[i(xc)]),
                xytext=(xc - 0.06, 0.62), fontsize=7.4, color=AMBER,
                arrowprops=dict(arrowstyle="-", color=AMBER, lw=0.7,
                                shrinkA=0, shrinkB=3))

    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(0.0, 1.55)
    ax.set_xlabel("in-plane coordinate")
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    return save(fig, "film-geometry")


if __name__ == "__main__":
    for f in (figure_shapes_1d(), figure_createtype(), figure_film_geometry()):
        print("wrote", f)
