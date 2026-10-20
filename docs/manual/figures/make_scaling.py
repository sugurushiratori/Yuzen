#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Draw the measured thread scaling of the Yuzen manual (Parallel execution).

Reads the TSV written by tools/scaling/scaling.sh and plots the speed-up
S = T(1)/T(p) of the time-step loop against the thread count p, for each
case size and backend, the fastest of the repeats of each point:

    python3 docs/manual/figures/make_scaling.py [tools/scaling/results/scaling-mc076.tsv]

writes parallel-scaling.pdf and parallel-scaling.svg next to itself. Run it
again after a new measurement; tab:scaling is printed by
tools/scaling/summarize.py from the same file.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
DEFAULT = os.path.join(ROOT, "tools", "scaling", "results", "scaling-mc076.tsv")

INK = "#22201D"
GRAY = "#9A968E"
#  The two backends: a blue and an amber of the logo palette, stepped so that
#  both clear 3:1 on white and stay apart under color-vision deficiency;
#  eigen is dashed as well, so color is never the only cue.
BACKEND = {"amgcl": ("#3D66B0", "-"), "eigen": ("#B27A28", (0, (4, 2)))}
MARKERS = ["o", "s", "^", "D"]

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["TeX Gyre Pagella", "DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "font.size": 8, "axes.titlesize": 8.5,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "yuzen",
})


def read(path):
    host = ""
    rows = []
    with open(path, encoding="utf-8") as f:
        lines = []
        for line in f:
            if line.startswith("#"):
                if line.startswith("# cpu"):
                    host = line.split("\t", 1)[1].strip()
                continue
            lines.append(line)
    for r in csv.DictReader(lines, delimiter="\t"):
        if r.get("status", "ok") != "ok":
            continue
        rows.append(r)
    best = {}
    for r in rows:
        k = (r["case"], r["backend"], int(r["threads"]))
        t = float(r["wall_s"])
        best[k] = min(best.get(k, t), t)
    return host, best


def label(case):
    n = case.split("_n")[1].lstrip("0")
    return r"$N_x = %s$" % n if case.startswith("bench1d") else r"$%s^2$" % n


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    host, best = read(path)
    cases = sorted({k[0] for k in best}, key=lambda c: (c[:7], int(c.split("_n")[1])))
    threads = sorted({k[2] for k in best})
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for ax, dim, title in ((axes[0], "bench1d", "1-D films"), (axes[1], "bench2d", "2-D films")):
        ax.plot(threads, threads, color=GRAY, lw=0.8, ls=(0, (1, 2)), zorder=1)
        ax.text(threads[-1] / 1.25, threads[-1] / 1.05, "ideal", color=GRAY, fontsize=6.5,
                ha="right", va="top", rotation=0)
        sizes = [c for c in cases if c.startswith(dim)]
        for i, c in enumerate(sizes):
            for be, (color, ls) in BACKEND.items():
                ps = [p for p in threads if (c, be, p) in best]
                if 1 not in ps:
                    continue
                s = [best[(c, be, 1)] / best[(c, be, p)] for p in ps]
                ax.plot(ps, s, color=color, ls=ls, lw=1.3, marker=MARKERS[i % 4], ms=4.2,
                        mec="white", mew=0.6, zorder=3)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.set_xticks(threads)
        ax.set_xticklabels([str(p) for p in threads])
        ax.set_yticks([0.25, 0.5, 1, 2, 4, 8, 16, 32])
        ax.set_yticklabels(["1/4", "1/2", "1", "2", "4", "8", "16", "32"])
        ax.set_ylim(0.15, 140)
        ax.minorticks_off()
        ax.axhline(1.0, color=GRAY, lw=0.5, zorder=0)
        ax.grid(axis="y", color="#E6E3DC", lw=0.5, zorder=0)
        ax.set_title(title, loc="left", color=INK, fontweight="bold")
        ax.set_xlabel("threads $p$")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        # legend: sizes by marker (in ink), backends by color and dash
        handles = [plt.Line2D([], [], color=INK, lw=0, marker=MARKERS[i % 4], ms=4.2, label=label(c))
                   for i, c in enumerate(sizes)]
        handles += [plt.Line2D([], [], color=col, ls=ls, lw=1.3, label='"%s"' % be)
                    for be, (col, ls) in BACKEND.items()]
        ax.legend(handles=handles, loc="upper left", frameon=False, ncol=2, handlelength=2.4,
                  columnspacing=1.0, labelcolor=INK)
    axes[0].set_ylabel(r"speed-up $S = T(1)/T(p)$")
    if host:
        fig.text(0.99, 0.01, host, ha="right", va="bottom", fontsize=6, color=GRAY)
    fig.tight_layout()
    out = os.path.join(HERE, "parallel-scaling.pdf")
    fig.savefig(out)
    fig.savefig(os.path.join(HERE, "parallel-scaling.svg"))
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
