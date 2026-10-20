#!/usr/bin/env python3
"""Write the seed corpora and the dictionary of the fuzz harnesses.

    python3 tools/fuzz_corpus.py --out DIR [--exe build/bin/yuzen]

The harnesses are in yuzen/tests/fuzz (built with -DYUZEN_BUILD_FUZZERS=ON).
This writes

  DIR/toml/     every case of cases/ (and cases/eq), as they are
  DIR/inp/      the same cases written in the block-structured .inp syntax, and the
                .inp samples kept in the repository's git history (when git
                is there)
  DIR/restart/  with --exe: field*.dat, .vts and field*.state written by
                short runs on the grids fuzz_restart builds, each prefixed
                with the byte that selects its reader in fuzz_restart.cpp
  DIR/input.dict    block paths, key names and accepted values (from the
                    parser, through tools/doc_keys.py), in libFuzzer's format
  DIR/restart.dict  the tags and keys of the restart files

Only the Python standard library is used.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

try:
    import tomllib
except ImportError:          # Python < 3.11
    tomllib = None

#  The grids of fuzz_restart.cpp: (selector of the .dat, of the .vts, model,
#  Dimension, Nx, Ny, Nz). Keep in step with F() there.
RESTART_GRIDS = [
    (0, 2, "WM", "xh", 8, 1, 1),
    (1, 3, "WM", "xyh", 6, 5, 1),
    (6, 4, "CZ", "xh", 6, 1, 4),
    (7, 5, "CZ", "xyh", 5, 4, 3),
]
STATE_SELECTOR = 8

RESTART_TEMPLATE = """SolverName = "Yuzen-1.0"
JobName = "r"
[Grid]
Dimension = "{dim}"
Xmin = 0
Xmax = 1e-3
Ymin = 0
Ymax = 1e-3
Nx = {nx}
Ny = {ny}
Nz = {nz}
[Initial.Concentration]
Type = "uniform"
Value = 0.5
[Initial.SurfaceShape]
TotalHeight = 10e-6
PrecursorRatio = 0.1
[Initial.SurfaceShape.EllipseBump]
Height = 0.5
Center = [0.5e-3, {cy}]
Width = [0.4e-3, 0.4e-3]
SlopeWidth = 0.1
[Boundary.XminBC]
[Boundary.XmaxBC]
[Fluid]
Name = "f"
[Fluid.Density]
Value = 1000.0
[Fluid.SurfaceTension]
BaseValue = 30.0e-3
[Fluid.Viscosity]
BaseValue = 100.0e-3
[Fluid.Diffusivity]
BaseValue = 1.0e-9
[Fluid.MolarMass]
[Numerics.ConcEquation.LinearSolver]
[Numerics.HeightEquation.LinearSolver]
[Numerics.TimeStepControl]
[Output]
Overwrite = true
[Output.CalcInfo]
[Output.GeneralHistory]
[Physics.ConcField]
Model = "{model}"
[Numerics]
InitialDt = 1e-4
EndTime = 3e-4
[Output.FieldData]
TimeInterval = 1e-4
NeutralFile = true
VTKFile = true
"""


def inp_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return '"%s"' % v
    if isinstance(v, list):
        return "(" + ", ".join(inp_value(x) for x in v) + ")"
    return repr(v)


def to_inp(doc, indent=""):
    """A parsed TOML document in the .inp syntax ConfigBlock reads."""
    out = []
    for k, v in doc.items():
        if isinstance(v, dict):
            out.append("%sblock %s" % (indent, k))
            out.extend(to_inp(v, indent + "  "))
            out.append("%send" % indent)
        elif isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
            for x in v:
                out.append("%sblock %s" % (indent, k))
                out.extend(to_inp(x, indent + "  "))
                out.append("%send" % indent)
        else:
            out.append("%s%s = %s" % (indent, k, inp_value(v)))
    return out


def case_files():
    return sorted(glob.glob(os.path.join(ROOT, "cases", "*.toml"))
                  + glob.glob(os.path.join(ROOT, "cases", "*", "*.toml")))


def write_inputs(out):
    tdir, idir = os.path.join(out, "toml"), os.path.join(out, "inp")
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(idir, exist_ok=True)
    n = 0
    for f in case_files():
        base = os.path.splitext(os.path.basename(f))[0]
        shutil.copy(f, os.path.join(tdir, base + ".toml"))
        n += 1
        if tomllib:
            with open(f, "rb") as fh:
                doc = tomllib.load(fh)
            with open(os.path.join(idir, base + ".inp"), "w") as fh:
                fh.write("\n".join(to_inp(doc)) + "\n")
    #  The .inp samples removed in b94a8a0
    try:
        names = subprocess.run(["git", "-C", ROOT, "ls-tree", "-r", "--name-only", "b94a8a0^"],
                               capture_output=True, text=True, check=True).stdout.split()
        for name in names:
            if name.endswith(".inp"):
                body = subprocess.run(["git", "-C", ROOT, "show", "b94a8a0^:" + name],
                                      capture_output=True, check=True).stdout
                with open(os.path.join(idir, "hist_" + name.replace("/", "_")), "wb") as fh:
                    fh.write(body)
    except (OSError, subprocess.CalledProcessError):
        pass
    return n


def write_restart(out, exe):
    rdir = os.path.join(out, "restart")
    work = os.path.join(out, "restart_runs")
    os.makedirs(rdir, exist_ok=True)
    os.makedirs(work, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    for sel_dat, sel_vts, model, dim, nx, ny, nz in RESTART_GRIDS:
        name = "%s_%s" % (model, dim)
        d = os.path.join(work, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "r.toml"), "w") as fh:
            #  A 1-D film is the line y = 0 of the shapes: the bump is centered on it.
            cy = "0.5e-3" if dim == "xyh" else "0.0"
            fh.write(RESTART_TEMPLATE.format(model=model, dim=dim, nx=nx, ny=ny, nz=nz, cy=cy))
        p = subprocess.run([exe, "r.toml"], cwd=d, env=env, capture_output=True, text=True)
        if p.returncode != 0:
            sys.exit("fuzz_corpus: the restart run %s failed:\n%s%s" % (name, p.stdout[-2000:], p.stderr))
        for k in (0, 2):
            stem = os.path.join(d, "r", "field%05d" % k)
            for ext, sel in ((".dat", sel_dat), (".vts", sel_vts), (".state", STATE_SELECTOR)):
                with open(stem + ext, "rb") as fh:
                    body = fh.read()
                with open(os.path.join(rdir, "%s_%d%s" % (name, k, ext)), "wb") as fh:
                    fh.write(bytes([sel]) + body)
    shutil.rmtree(work, ignore_errors=True)


def dict_entry(s):
    return '"' + "".join(c if (32 <= ord(c) < 127 and c not in '"\\') else "\\x%02x" % ord(c)
                         for c in s) + '"'


def write_dicts(out):
    import doc_keys
    blocks, keys, _ = doc_keys.code_model()
    words = set()
    for b in blocks:
        if b:
            words.add("[%s]" % b)
            words.add("[[%s]]" % b)
            words.add("block %s" % b.split(".")[-1])
            words.add(b.split(".")[-1])
    for k in keys.values():
        words.add(k.name)
        words.add(k.name + " = ")
        for v in (k.enum or []):
            words.add('"%s"' % v)
    words.update(["end", "true", "false", "nan", "inf", "-inf", "1e308", "-1e308", "0",
                  "-1", "1e-308", "4e-324", "2147483648", "9223372036854775808", "[]", "( )",
                  "{}", "=", "\"\"", "'''", '"""', "0x7fffffff", "1_000", "+inf"])
    with open(os.path.join(out, "input.dict"), "w") as fh:
        for w in sorted(words):
            fh.write(dict_entry(w) + "\n")
    rwords = ["# Yuzen restart state 1", "binary", "vector ", "text ", "model",
              "Well-Mixed-1D", "Well-Mixed-2D", "ConcZ-2D", "ConcZ-3D", "time ", "step ",
              "dt ", "out.Field", "out.Probe", "out.History", "out.Calc", "steadyCount",
              "adjustCount.height", "height.old", "height.correction",
              "conc.previous", "conc.previousDt", "field.h", "field.s", "field.c",
              "volume0", "evaporated", "outflowed",
              "<YuzenData", "/>", "Nx=\"", "Ny=\"", "Nz=\"", "Model=\"ConcZ\"",
              "<AppendedData encoding=\"raw\">", "   _", "Name=\"FilmThickness\"",
              "\xff\xff\xff\xff", "\x00\x00\x00\x00", "nan", "inf", "-1", "1e308", "\t", "#"]
    with open(os.path.join(out, "restart.dict"), "w") as fh:
        for w in rwords:
            fh.write(dict_entry(w) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--exe", help="yuzen binary, for the restart corpus")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    n = write_inputs(a.out)
    write_dicts(a.out)
    if a.exe:
        write_restart(a.out, os.path.abspath(a.exe))
    print("fuzz_corpus: %d cases, corpora and dictionaries in %s" % (n, a.out))


if __name__ == "__main__":
    main()
