#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Draw the example-result figures of the Yuzen manual and the README.

Unlike docs/manual/figures/make_shapes.py, which draws the shape functions
from their formulas, everything here is drawn from output the solver actually
produced. Nothing is idealized or redrawn by hand: if a figure looks wrong,
the run behind it is wrong.

    python3 docs/figures/make_examples.py            # all of them
    python3 docs/figures/make_examples.py --run      # run the cases first

The runs take a few minutes altogether. Without --run the script expects the
result directories to be under the working directory already, which is what
the release workflow does.

Writes, next to itself:

    leveling-1d.{png,pdf,svg}   the 1-D film height over a substrate step
    leveling-1d.gif             the same as an animation
    leveling-2d.{png,pdf,svg}   the 2-D free surface at three times
    leveling-2d.gif             the same as an animation
    wetting.{png,pdf,svg}        a droplet relaxing to its contact angle
    concz.{png,pdf,svg}          concentration through the film thickness

The four forms have one consumer each: .png for the repository README, .pdf
for the LaTeX manual, .svg for the HTML manual, which inlines it, and .gif for
the two cases worth animating.

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

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CASES = os.path.join(ROOT, "cases")

# --- the Yuzen palette, as in docs/manual/figures/make_shapes.py -----------
INK = "#22201D"
CELADON = "#5E8C78"
SHADOW = "#3C5F4E"
AMBER = "#C08A3E"
LAPIS = "#2E4A7D"
PAPER = "#F4F1EA"
GRAY = "#9A968E"

#  A glaze-like ramp for the surface maps: from the paper of an unglazed body,
#  through celadon, to the dark of a thick pour. It runs light to dark so that
#  it survives being printed in gray.
GLAZE = LinearSegmentedColormap.from_list(
    "glaze", [PAPER, "#BBD0C4", CELADON, SHADOW, "#233A30"])

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["TeX Gyre Pagella", "Palatino", "DejaVu Serif"],
    "font.size": 8.0,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8.0,
    "axes.edgecolor": INK,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.fontsize": 7.0,
    "legend.frameon": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "svg.fonttype": "path",
    #  Without a fixed salt the element ids in an SVG are randomized per
    #  session, so two runs of this script over the same results would differ.
    "svg.hashsalt": "yuzen",
})


# ---------------------------------------------------------------------------
#  reading what the solver wrote
# ---------------------------------------------------------------------------

def fields(resultdir):
    """Every field*.dat of a result directory, in time order."""
    fs = sorted(glob.glob(os.path.join(resultdir, "field*.dat")))
    if not fs:
        sys.exit("make_examples.py: no field*.dat under %s.\n"
                 "  Run the cases first: %s --run"
                 % (resultdir, os.path.basename(__file__)))
    return fs


def columns(path):
    """Column names of a field file, with the unit stripped off each."""
    head = open(path, encoding="utf-8").readline().lstrip("#").strip()
    return [c.split("[")[0].strip() for c in head.split("\t") if c.strip()]


def read(path):
    """A field file as {name: column}, keyed on the bare column name."""
    d = np.loadtxt(path, skiprows=1)
    return dict(zip(columns(path), d.T))


def grid2d(f, name):
    """One column of a 2-D field file as (x, y, Z) ready for pcolormesh."""
    xs, ys = np.unique(f["xx"]), np.unique(f["yy"])
    z = f[name].reshape(len(ys), len(xs)) if f["xx"][1] != f["xx"][0] \
        else f[name].reshape(len(xs), len(ys)).T
    return xs, ys, z


def times(resultdir):
    """The time of each field file, read from the ParaView collection.

    The `.pvd` carries one `timestep` per written field, which is the only
    place the time of a given `field*.dat` is recorded; `History.dat` is
    written on its own interval and does not line up with it.
    """
    pvd = glob.glob(os.path.join(resultdir, "*.pvd"))
    if not pvd:
        return None
    import re
    text = open(pvd[0], encoding="utf-8").read()
    return np.array([float(t) for t in
                     re.findall(r'timestep="([^"]+)"', text)])


def label_time(t):
    """A time in seconds, written the way the run log writes it."""
    return "$t = %s$ s" % ("0" if t == 0 else ("%.3g" % t))


# ---------------------------------------------------------------------------
#  running the cases
# ---------------------------------------------------------------------------

#  The cases behind the figures, and the field-output interval each one needs.
#  The shipped cases write a field only every so often, which is right for them
#  and too coarse for an animation; None means the case's own interval is left
#  alone. Nothing else about any case is touched.
PLAN = [("allphys_2d_reg", "6e-3"),
        ("leveling_1d_reg", "4e-3"),
        ("concz_2d", None),
        ("wet_static", None)]


def run_cases(workdir, plan=None):
    """Run the cases the figures are drawn from, into workdir.

    Each case is run with `workdir` as the working directory, so the solver
    puts its results in `workdir/<JobName>` -- and every case sets `JobName`
    to its own file name, which is where the figures look for them.
    """
    exe = os.path.join(ROOT, "build", "bin", "yuzen")
    if not os.path.isfile(exe):
        exe = os.path.join(ROOT, "build", "bin", "Release", "yuzen.exe")
    if not os.path.isfile(exe):
        sys.exit("make_examples.py: build/bin/yuzen not found; build it first.")

    for case, interval in (plan or PLAN):
        src = open(os.path.join(CASES, case + ".toml"), encoding="utf-8").read()
        if interval:
            src, n = re.subn(r"(\[Output\.FieldData\]\n)TimeInterval = [^\n]*",
                             lambda m: m.group(1) + "TimeInterval = " + interval,
                             src)
            if n != 1:
                sys.exit("make_examples.py: could not set the field interval of"
                         " %s (matched %d times)." % (case, n))
        open(os.path.join(workdir, case + ".toml"), "w",
             encoding="utf-8").write(src)
        #  a case that starts from a restart file needs the file beside it
        for dat in set(re.findall(r'"([^"]+\.dat)"', src)):
            if os.path.exists(os.path.join(CASES, dat)):
                shutil.copy(os.path.join(CASES, dat), workdir)
        print("  running %s ..." % case, flush=True)
        r = subprocess.run([exe, case + ".toml"], cwd=workdir,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            sys.exit("make_examples.py: %s failed:\n%s" % (case, r.stdout[-2000:]))
    return workdir


# ---------------------------------------------------------------------------
#  output
# ---------------------------------------------------------------------------

#  PDF and SVG carry a creation date unless one is suppressed, which would
#  make every rerun of this script a change to every figure in the repository
#  even when the runs behind them are identical. PNG carries none.
NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}


def save(fig, stem, png_dpi=200):
    """PNG for the README and the HTML, PDF for LaTeX, SVG for the record."""
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(HERE, "%s.%s" % (stem, ext)),
                    metadata=NO_DATE[ext])
    fig.savefig(os.path.join(HERE, stem + ".png"), dpi=png_dpi)
    plt.close(fig)
    print("  wrote %s.{png,pdf,svg}" % stem)


def save_gif(frames, stem, fps=12):
    """Write a GIF from a list of RGB arrays, via Pillow."""
    from PIL import Image
    imgs = [Image.fromarray(f).convert("P", palette=Image.ADAPTIVE, colors=128)
            for f in frames]
    out = os.path.join(HERE, stem + ".gif")
    imgs[0].save(out, save_all=True, append_images=imgs[1:],
                 duration=int(1000 / fps), loop=0, optimize=True)
    kb = os.path.getsize(out) / 1024.0
    print("  wrote %s.gif (%d frames, %.0f kB)" % (stem, len(imgs), kb))


def render(fig):
    """A figure as an RGB array, for assembling into a GIF."""
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    return buf.reshape(h, w, 4)[:, :, :3].copy()


# ---------------------------------------------------------------------------
#  the figures
# ---------------------------------------------------------------------------

def fig_leveling_2d(work):
    """The 2-D free surface at three times of the run."""
    d = os.path.join(work, "allphys_2d_reg")
    fs = fields(d)
    t = times(d)
    pick = [0, len(fs) // 4, len(fs) - 1]
    frames = [read(fs[i]) for i in pick]

    fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.3))
    lo = min(f["HS"].min() for f in frames)
    hi = max(f["HS"].max() for f in frames)
    for ax, i, f in zip(axes, pick, frames):
        xs, ys, z = grid2d(f, "HS")
        m = ax.pcolormesh(xs, ys, z, cmap=GLAZE, vmin=lo, vmax=hi,
                          shading="gouraud", rasterized=True)
        ax.contour(xs, ys, z, levels=8, colors=INK, linewidths=0.3, alpha=0.5)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$ [mm]")
        ax.set_title(label_time(t[i]) if t is not None else "frame %d" % i)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axes[0].set_ylabel("$y$ [mm]")
    for ax in axes[1:]:
        ax.set_yticklabels([])
    cb = fig.colorbar(m, ax=axes, fraction=0.030, pad=0.02)
    cb.set_label(r"free surface $h+s$  [$\mu$m]")
    cb.outline.set_linewidth(0.6)
    save(fig, "leveling-2d")


def gif_leveling_2d(work):
    """The 2-D free surface leveling, as an animation."""
    d = os.path.join(work, "allphys_2d_reg")
    fs = fields(d)
    if len(fs) < 4:
        print("  (skipping leveling-2d.gif: only %d frames)" % len(fs))
        return
    t = times(d)
    all_hs = [read(f)["HS"] for f in fs]
    lo = min(h.min() for h in all_hs)
    hi = max(h.max() for h in all_hs)

    frames = []
    for i, path in enumerate(fs):
        f = read(path)
        fig, ax = plt.subplots(figsize=(3.0, 2.6), dpi=110)
        xs, ys, z = grid2d(f, "HS")
        m = ax.pcolormesh(xs, ys, z, cmap=GLAZE, vmin=lo, vmax=hi,
                          shading="gouraud")
        ax.contour(xs, ys, z, levels=8, colors=INK, linewidths=0.3, alpha=0.5)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$ [mm]")
        ax.set_ylabel("$y$ [mm]")
        ax.set_title(label_time(t[i]) if t is not None
                     else "frame %d/%d" % (i + 1, len(fs)))
        cb = fig.colorbar(m, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label(r"$h+s$ [$\mu$m]")
        fig.tight_layout()
        frames.append(render(fig))
        plt.close(fig)
    save_gif(frames, "leveling-2d")


def fig_leveling_1d(work):
    """The 1-D film height at a few times, one curve per output."""
    d = os.path.join(work, "leveling_1d_reg")
    fs = fields(d)
    t = times(d)

    fig, ax = plt.subplots(figsize=(5.4, 2.6))
    pick = np.linspace(0, len(fs) - 1, min(6, len(fs))).astype(int)
    for k, i in enumerate(pick):
        f = read(fs[i])
        c = GLAZE(0.15 + 0.8 * k / max(1, len(pick) - 1))
        lab = None
        if t is not None and k in (0, len(pick) - 1):
            lab = label_time(t[i])
        ax.plot(f["xx"], f["HS"], color=c, lw=1.2, label=lab)
    f0 = read(fs[0])
    ax.fill_between(f0["xx"], 0, f0["SS"], color=GRAY, alpha=0.35, lw=0)
    ax.plot(f0["xx"], f0["SS"], color=GRAY, lw=0.8)
    ax.set_xlabel("$x$ [mm]")
    ax.set_ylabel(r"height  [$\mu$m]")
    ax.set_title("capillary leveling over a substrate step")
    ax.legend(loc="upper left")
    save(fig, "leveling-1d")


def gif_leveling_1d(work):
    """The 1-D film height, as an animation."""
    d = os.path.join(work, "leveling_1d_reg")
    fs = fields(d)
    if len(fs) < 4:
        print("  (skipping leveling-1d.gif: only %d frames)" % len(fs))
        return
    t = times(d)
    hi = max(read(f)["HS"].max() for f in fs) * 1.05

    frames = []
    f0 = read(fs[0])
    for i, path in enumerate(fs):
        f = read(path)
        fig, ax = plt.subplots(figsize=(3.6, 2.2), dpi=110)
        ax.fill_between(f0["xx"], 0, f0["SS"], color=GRAY, alpha=0.35, lw=0)
        ax.fill_between(f["xx"], f["SS"], f["HS"], color=CELADON, alpha=0.75,
                        lw=0)
        ax.plot(f["xx"], f["HS"], color=SHADOW, lw=1.2)
        ax.set_xlim(f["xx"].min(), f["xx"].max())
        ax.set_ylim(0, hi)
        ax.set_xlabel("$x$ [mm]")
        ax.set_ylabel(r"height [$\mu$m]")
        ax.set_title(label_time(t[i]) if t is not None
                     else "frame %d/%d" % (i + 1, len(fs)))
        fig.tight_layout()
        frames.append(render(fig))
        plt.close(fig)
    save_gif(frames, "leveling-1d")


def concz_slice(f):
    """A CZ field as (x, z, c) on a plane of constant y.

    `zz` is the height of the point, so the mesh comes straight out of the
    file; it follows the free surface, which is why z is two-dimensional and
    not a single axis. The slice is taken at the y nearest the middle of the
    domain.
    """
    nx = len(np.unique(f["xx"]))
    ny = len(np.unique(f["yy"])) if "yy" in f else 1
    nz = len(f["CC"]) // (nx * ny)
    shape = (nz, ny, nx)
    xs = f["xx"].reshape(shape)[0, 0, :]
    j = int(np.argmin(np.abs(f["yy"].reshape(shape)[0, :, 0]))) if ny > 1 else 0
    z = f["zz"].reshape(shape)[:, j, :]
    cc = f["CC"].reshape(shape)[:, j, :]
    hs = f["HS"].reshape(shape)[0, j, :]
    ss = f["SS"].reshape(shape)[0, j, :]
    return xs, z, cc, hs, ss


def fig_concz(work):
    """Concentration through the film thickness: what the CZ model is for."""
    d = os.path.join(work, "concz_2d")
    fs = fields(d)
    xs, z, cc, hs, ss = concz_slice(read(fs[-1]))
    x2 = np.repeat(xs[None, :], z.shape[0], axis=0)
    #  The depth average over the film, which is the single value a well-mixed
    #  run would carry in place of the whole profile.
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    bar = trapz(cc, z, axis=0) / np.maximum(hs - ss, 1e-300)
    #  The differences that matter are in the seventh decimal, so the value
    #  the film started at is taken out and named, rather than left for the
    #  axis to carry in a corner label.
    c0 = float(np.max(read(fs[0])["CC"]))
    cc = (cc - c0) * 1e6
    bar = (bar - c0) * 1e6

    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.3),
                             gridspec_kw={"width_ratios": [1.5, 1.3, 0.8],
                                          "wspace": 0.80})

    #  (a) where the film is, and where the solute sits in it. The layer is
    #  confined to where the film thins; over the flat film there is nothing
    #  to see, which is what panel (b) makes explicit.
    ix = int(np.argmax(cc.max(axis=0) - cc.min(axis=0)))
    m = axes[0].pcolormesh(x2, z, cc, cmap=GLAZE, shading="gouraud",
                           rasterized=True)
    axes[0].plot(xs, hs, color=INK, lw=0.9)
    axes[0].fill_between(xs, 0, ss, color=GRAY, alpha=0.45, lw=0)
    axes[0].plot(xs, ss, color=GRAY, lw=0.8)
    axes[0].set_xlim(xs.min(), xs.max())
    axes[0].set_ylim(0, hs.max() * 1.08)
    axes[0].set_xlabel("$x$ [mm]")
    axes[0].set_ylabel(r"$z$ [$\mu$m]")
    axes[0].set_title(r"$10^{6}\,(c - c_0)$,  $c_0 = %.6f$" % c0)
    axes[0].xaxis.set_major_locator(MaxNLocator(nbins=4))
    cb = fig.colorbar(m, ax=axes[0], fraction=0.046, pad=0.025)
    #  The values sit on a large common offset, and matplotlib would write it
    #  as a separate label above the bar, which lands on top of the panel.
    #  Plain six-digit ticks are wider but they stay where they are put.
    cb.locator = MaxNLocator(nbins=4)
    cb.update_ticks()
    cb.ax.tick_params(labelsize=6.0, pad=1.5)
    cb.outline.set_linewidth(0.6)

    #  (b) the point of the model: the free-surface value is what sets the
    #  surface tension, and it is not the depth average a well-mixed run
    #  would use in its place.
    axes[1].plot(xs, cc[-1], color=SHADOW, lw=1.4, label="free surface")
    axes[1].plot(xs, bar, color=AMBER, lw=1.2, ls="--", label="depth average")
    axes[1].plot(xs, cc[0], color=GRAY, lw=1.0, label="substrate")
    axes[1].set_xlim(xs.min(), xs.max())
    axes[1].set_xlabel("$x$ [mm]")
    axes[1].set_ylabel(r"$10^{6}\,(c - c_0)$")
    axes[1].xaxis.set_major_locator(MaxNLocator(nbins=4))
    axes[1].set_title("surface value vs. depth average")
    axes[1].set_ylim(cc.min(), cc.min() + 1.45 * (cc[-1].max() - cc.min()))
    axes[1].legend(loc="upper center", ncol=1, fontsize=6.0)
    axes[1].yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))

    #  (c) a profile through the thickness where the stratification is
    #  strongest, drawn on the grid the model actually uses.
    axes[2].plot(cc[:, ix], z[:, ix], color=SHADOW, lw=1.4, marker="o",
                 ms=2.0, mfc=PAPER, mew=0.6)
    axes[2].axvline(bar[ix], color=AMBER, lw=1.0, ls="--")
    axes[2].set_xlabel(r"$10^{6}\,(c - c_0)$")
    axes[2].set_ylabel(r"$z$  [$\mu$m]")
    axes[2].set_title("profile at $x = %.2f$ mm" % xs[ix])
    axes[2].xaxis.set_major_locator(MaxNLocator(nbins=3, steps=[1, 2, 5, 10]))
    axes[2].tick_params(axis="x", labelsize=6.0)
    save(fig, "concz")


def fig_wetting(work):
    """An equilibrium droplet: the contact line the disjoining term makes."""
    d = os.path.join(work, "wet_static")
    fs = fields(d)
    first, last = read(fs[0]), read(fs[-1])

    #  The precursor thickness h* the input asks for (TotalHeight times
    #  PrecursorRatio). With h*/h0 = 0.01 it is invisible at the scale of the
    #  drop, so an inset shows the contact line and the precursor film.
    src = open(os.path.join(work, "wet_static.toml"), encoding="utf-8").read()
    shape = src.split("[Initial.SurfaceShape]", 1)[1]
    num = lambda k: float(re.search(r"^%s\s*=\s*([-0-9.eE+]+)" % k, shape, re.M).group(1))
    hstar = num("TotalHeight") * num("PrecursorRatio") * 1e6       # [um]

    fig, ax = plt.subplots(figsize=(4.8, 2.4))
    ax.plot(first["xx"], first["HS"], color=GRAY, lw=1.0, ls="--",
            label="initial")
    ax.fill_between(last["xx"], 0, last["HS"], color=CELADON, alpha=0.7, lw=0)
    ax.plot(last["xx"], last["HS"], color=SHADOW, lw=1.4, label="relaxed")
    ax.axhline(hstar, color=AMBER, lw=0.9, ls=":")
    ins = ax.inset_axes([0.06, 0.42, 0.34, 0.5])
    ins.plot(first["xx"], first["HS"], color=GRAY, lw=0.9, ls="--")
    ins.fill_between(last["xx"], 0, last["HS"], color=CELADON, alpha=0.7, lw=0)
    ins.plot(last["xx"], last["HS"], color=SHADOW, lw=1.2)
    ins.axhline(hstar, color=AMBER, lw=0.9, ls=":")
    ins.annotate(r"$h_*$", xy=(0.10, hstar), xytext=(-2, 2), textcoords="offset points",
                 ha="right", va="bottom", color=AMBER, fontsize=6.5)
    ins.set_xlim(0.05, 0.10)
    ins.set_ylim(0, 1.5)
    ins.tick_params(labelsize=5.5, length=2)
    ins.set_title("contact line", fontsize=6.5, pad=2)
    ax.set_xlim(last["xx"].min(), last["xx"].max())
    ax.set_ylim(0, first["HS"].max() * 1.12)
    ax.set_xlabel("$x$ [mm]")
    ax.set_ylabel(r"height  [$\mu$m]")
    ax.set_title("a droplet relaxing to its equilibrium contact angle")
    ax.legend(loc="upper right")
    save(fig, "wetting")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="store_true",
                    help="run the cases first (a few minutes)")
    ap.add_argument("--workdir", default=os.path.join(ROOT, "runs", "figures"),
                    help="where the result directories are (or go)")
    args = ap.parse_args()

    work = args.workdir
    if args.run:
        os.makedirs(work, exist_ok=True)
        print("running the cases into %s" % work)
        run_cases(work)

    print("drawing:")
    fig_leveling_2d(work)
    gif_leveling_2d(work)
    fig_leveling_1d(work)
    gif_leveling_1d(work)
    fig_concz(work)
    fig_wetting(work)
    return 0


if __name__ == "__main__":
    sys.exit(main())
