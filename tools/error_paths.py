#!/usr/bin/env python3
"""Check what the solver does when the user, the disk or the job scheduler
gets in the way.

    python3 tools/error_paths.py --exe build/bin/yuzen --workdir /tmp/ep

Every case here ran, before it was looked at, to a wrong end: a result
directory that could not be written to ran to the end with exit code 0 and
no output at all; a full disk left empty field files and empty .state files
behind a run that ended with exit code 0; a JobName with a directory in it
lost its .pvd; a restart file that was not there left a result directory set
up for a run that never started; running out of memory dumped core; and a
ConcZ .pvd was never closed.

  cli         no argument, an unknown option, two input files: exit code 1
              and the usage; --help and --version: exit code 0
  input       an input file that is missing, a directory, empty or binary
              garbage (and, run as a user other than root, unreadable):
              refused with the path, with and without --check-input, and no
              result directory made
  resultdir   a result directory that exists (Output.Overwrite = false),
              a file in its place, a parent that is a file, and (as a user
              other than root) one that cannot be written to: refused before
              anything is computed, with the path, JobName and the reason;
              a JobName with missing parents is made, its .pvd inside it
  restart     Restart.FileName missing: refused before the result directory
              is made, noted by --check-input; a truncated .vts or .dat
              refused; a truncated .state not continued from, with a warning
  diskfull    a file size limit (RLIMIT_FSIZE, as ulimit -f, with SIGXFSZ
              ignored by the solver) that the first .vts exceeds, and a
              History.dat on /dev/full: exit code 1, the file and the reason
              on standard error and in FatalErrors.txt, no field file or
              .state cut short, "abort" in JobInfo.txt
  sigint      Ctrl-C (SIGINT, a well-mixed run) and SIGTERM (a ConcZ run)
              in the middle of a run: the step finishes, the state reached
              is written, exit code 128 + the signal, "interrupted" in
              JobInfo.txt, a .pvd that is valid XML; the run continued from
              its last field file and .state is, to the bit, the run made
              straight through (as tools/restart_roundtrip.py --chain)
  sigkill     the run killed at once (SIGKILL) at a random moment, possibly
              in the middle of writing: no field file, .state or .pvd is
              left cut short, and the run continued from the last field
              file with its .state is again, to the bit, the straight run
  paths       a JobName with spaces and Japanese characters, a long path, an
              input file given by a relative path from another directory and
              through a symbolic link, and a Restart.FileName with spaces and
              Japanese characters

The cases that need signals or resource limits are skipped where the
platform has none (Windows). Only the Python standard library is used, so
the test runs wherever ctest does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import restart_roundtrip as rr  # noqa: E402

CASES = os.path.join(os.path.dirname(HERE), "cases")
POSIX = os.name == "posix"
ROOT = POSIX and os.geteuid() == 0

#  A small 1-D case: well under a second a run
SMALL = {"Grid": {"Nx": 60}, "Output": {"Overwrite": "false"},
         "Numerics": {"EndTime": "0.02"}}


def small(job, extra=None):
    sect = {k: dict(v) for k, v in SMALL.items()}
    sect.setdefault("", {})["JobName"] = '"%s"' % job
    for k, v in (extra or {}).items():
        sect.setdefault(k, {}).update(v)
    return rr.toml_edit(open(os.path.join(CASES, "smoke_1d.toml"), encoding="utf-8").read(), sect)


def write(path, text, mode="w"):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, mode, **({} if "b" in mode else {"encoding": "utf-8"})) as f:
        f.write(text)


def run(exe, args, cwd, timeout=300, preexec=None, env=None):
    e = dict(os.environ)
    e.setdefault("OMP_NUM_THREADS", "1")
    e.setdefault("OMP_WAIT_POLICY", "PASSIVE")
    if env:
        e.update(env)
    try:
        r = subprocess.run([exe] + args, cwd=cwd, capture_output=True, timeout=timeout,
                           preexec_fn=preexec, env=e)
    except subprocess.TimeoutExpired:
        return -1, "", "timed out"
    return (r.returncode, r.stdout.decode("utf-8", "replace"),
            r.stderr.decode("utf-8", "replace"))


class Checks:
    def __init__(self, name):
        self.name, self.problems, self.notes = name, [], []

    def expect(self, cond, what):
        if not cond:
            self.problems.append(what)
        return cond


def short(s, n=300):
    s = s.strip().replace("\n", " | ")
    return s if len(s) <= n else s[:n] + "..."


#--------------------------------------------------------------------------
def case_cli(exe, wd):
    c = Checks("cli")
    write(os.path.join(wd, "a.toml"), small("a"))
    for args, what in (([], "no argument"), (["--frobnicate", "a.toml"], "an unknown option"),
                       (["a.toml", "a.toml"], "two input files")):
        rc, out, err = run(exe, args, wd)
        c.expect(rc == 1 and "usage:" in err, "%s: rc=%d, %s" % (what, rc, short(err)))
    c.expect(not os.path.exists(os.path.join(wd, "a")), "a result directory made by a command-line error")
    rc, out, err = run(exe, ["--help"], wd)
    c.expect(rc == 0 and "usage:" in out and "--check-input" in out, "--help: rc=%d %s" % (rc, short(out)))
    rc, out, err = run(exe, ["--version"], wd)
    c.expect(rc == 0 and re.match(r"Yuzen-\d+\.\d+\.\d+\s*$", out), "--version: rc=%d %s" % (rc, short(out)))
    return c


def case_input(exe, wd):
    c = Checks("input")
    os.makedirs(os.path.join(wd, "adir.toml"), exist_ok=True)
    write(os.path.join(wd, "empty.toml"), "")
    write(os.path.join(wd, "garbage.toml"), bytes(range(256)) * 8, "wb")
    write(os.path.join(wd, "unreadable.toml"), small("unreadable"))
    files = [("missing.toml", "does not exist"), ("adir.toml", "is a directory"),
             ("empty.toml", "SolverName"), ("garbage.toml", "garbage.toml")]
    if POSIX and not ROOT:
        os.chmod(os.path.join(wd, "unreadable.toml"), 0)
        files.append(("unreadable.toml", "cannot be read"))
    else:
        c.notes.append("unreadable input not tried (run as root or not POSIX)")
    for f, msg in files:
        for args in ([f], ["--check-input", f]):
            rc, out, err = run(exe, args, wd)
            c.expect(rc == 1 and msg in err, "%s: rc=%d, %s" % (" ".join(args), rc, short(err)))
    c.expect(not [d for d in os.listdir(wd) if d not in [f for f, _ in files] + ["unreadable.toml"]],
             "files left by refused inputs: %s" % os.listdir(wd))
    return c


def case_resultdir(exe, wd):
    c = Checks("resultdir")
    write(os.path.join(wd, "a.toml"), small("res"))
    os.makedirs(os.path.join(wd, "res"))
    write(os.path.join(wd, "res", "keep.txt"), "x")
    rc, out, err = run(exe, ["a.toml"], wd)
    c.expect(rc == 1 and '"res"' in err and "Output.Overwrite" in err,
             "existing directory: rc=%d %s" % (rc, short(err)))
    c.expect(os.listdir(os.path.join(wd, "res")) == ["keep.txt"], "existing directory written into")

    write(os.path.join(wd, "b.toml"), small("afile"))
    write(os.path.join(wd, "afile"), "x")
    rc, out, err = run(exe, ["b.toml"], wd)
    c.expect(rc == 1 and '"afile"' in err and "not a directory" in err,
             "a file in its place: rc=%d %s" % (rc, short(err)))

    write(os.path.join(wd, "d.toml"), small("afile/sub"))
    rc, out, err = run(exe, ["d.toml"], wd)
    c.expect(rc == 1 and '"afile/sub"' in err and "cannot be created" in err,
             "a parent that is a file: rc=%d %s" % (rc, short(err)))

    if POSIX and not ROOT:
        os.makedirs(os.path.join(wd, "ro"))
        os.chmod(os.path.join(wd, "ro"), 0o555)
        write(os.path.join(wd, "e.toml"), small("ro", {"Output": {"Overwrite": "true"}}))
        rc, out, err = run(exe, ["e.toml"], wd)
        c.expect(rc == 1 and "cannot be written to" in err, "read-only directory: rc=%d %s" % (rc, short(err)))
        write(os.path.join(wd, "f.toml"), small("ro/sub"))
        rc, out, err = run(exe, ["f.toml"], wd)
        c.expect(rc == 1 and "cannot be created" in err, "read-only parent: rc=%d %s" % (rc, short(err)))
        os.chmod(os.path.join(wd, "ro"), 0o755)
    else:
        c.notes.append("read-only directories not tried (run as root or not POSIX)")

    #  missing parents are made, and the .pvd goes inside, named after the
    #  last part of JobName
    write(os.path.join(wd, "g.toml"), small("p/q/r"))
    rc, out, err = run(exe, ["g.toml"], wd)
    c.expect(rc == 0 and os.path.isfile(os.path.join(wd, "p", "q", "r", "r.pvd")),
             "JobName p/q/r: rc=%d, files %s" % (rc, sorted(os.listdir(os.path.join(wd, "p", "q", "r")))
                                                  if os.path.isdir(os.path.join(wd, "p", "q", "r")) else "none"))
    return c


def case_restart(exe, wd):
    c = Checks("restart")
    #  a source run to restart from
    write(os.path.join(wd, "src.toml"), small("src"))
    rc, out, err = run(exe, ["src.toml"], wd)
    if not c.expect(rc == 0, "source run: rc=%d %s" % (rc, short(err))):
        return c
    ent = rr.pvd_entries(os.path.join(wd, "src"))
    t, f = ent[-1]
    idx = int(re.search(r"(\d+)\.vts$", f).group(1))

    def restart_input(job, fname, ascii_=False):
        return small(job, {"Output": {"Overwrite": "true"},
                           "Restart": {"Enable": "true", "ASCII": "true" if ascii_ else "false",
                                       "FileName": '"%s"' % fname, "RestartTime": t,
                                       "FileIndex": str(idx + 1)}})

    write(os.path.join(wd, "m.toml"), restart_input("m", "nothere/field00003.vts"))
    rc, out, err = run(exe, ["m.toml"], wd)
    c.expect(rc == 1 and "Restart.FileName" in err and "does not exist" in err,
             "missing restart file: rc=%d %s" % (rc, short(err)))
    c.expect(not os.path.exists(os.path.join(wd, "m")), "a result directory made for a missing restart file")
    rc, out, err = run(exe, ["--check-input", "m.toml"], wd)
    c.expect(rc == 0 and "Note: the restart file" in out, "--check-input, missing restart file: rc=%d %s"
             % (rc, short(out[-300:])))

    #  truncated field files are refused, by both readers
    base = os.path.join(wd, "src", f[:-4])
    for ext, ascii_ in ((".vts", False), (".dat", True)):
        os.makedirs(os.path.join(wd, "cut"), exist_ok=True)
        data = open(base + ext, "rb").read()
        write(os.path.join(wd, "cut", "field" + ext), data[:len(data) * 2 // 3], "wb")
        shutil.copy(base + ".state", os.path.join(wd, "cut", "field.state"))
        write(os.path.join(wd, "t.toml"), restart_input("t" + ext[1:], "cut/field" + ext, ascii_))
        rc, out, err = run(exe, ["t.toml"], wd)
        c.expect(rc == 1 and "Restart" in err and "cut/field" + ext in err,
                 "truncated %s: rc=%d %s" % (ext, rc, short(err)))

    #  a truncated .state beside a whole field file: the field alone, with a warning
    os.makedirs(os.path.join(wd, "cut2"), exist_ok=True)
    shutil.copy(base + ".vts", os.path.join(wd, "cut2", "field.vts"))
    st = open(base + ".state", "rb").read()
    write(os.path.join(wd, "cut2", "field.state"), st[:len(st) - 100], "wb")
    write(os.path.join(wd, "u.toml"), restart_input("u", "cut2/field.vts"))
    rc, out, err = run(exe, ["u.toml"], wd)
    warn = open(os.path.join(wd, "u", "Warnings.txt"), encoding="utf-8").read() if os.path.isdir(os.path.join(wd, "u")) else ""
    c.expect(rc == 0 and "continuing exactly" not in out and "is damaged" in out and "is damaged" in warn,
             "truncated .state: rc=%d %s / %s" % (rc, short(out[-400:]), short(warn)))
    return c


def no_cut_files(c, d, what):
    """Every field file, .state and .pvd in d is whole (the readers and
    the format say so); a stray .tmp is allowed."""
    for p in glob.glob(os.path.join(d, "field*.dat")):
        rows = [ln for ln in open(p, encoding="utf-8").read().split("\n")[1:] if ln.strip()]
        c.expect(rows and all(len(r.split()) == len(rows[0].split()) for r in rows)
                 and open(p, encoding="utf-8").read().endswith("\n"), "%s: %s cut short" % (what, os.path.basename(p)))
    for p in glob.glob(os.path.join(d, "field*.vts")):
        c.expect(open(p, "rb").read().rstrip().endswith(b"</VTKFile>"),
                 "%s: %s cut short" % (what, os.path.basename(p)))
    for p in glob.glob(os.path.join(d, "field*.state")):
        c.expect(os.path.getsize(p) > 100 and open(p, "rb").read(21) == b"# Yuzen restart state",
                 "%s: %s cut short" % (what, os.path.basename(p)))
    for p in glob.glob(os.path.join(d, "*.pvd")):
        try:
            ET.parse(p)
        except ET.ParseError as e:
            c.expect(False, "%s: %s is not valid XML (%s)" % (what, os.path.basename(p), e))


def case_diskfull(exe, wd):
    c = Checks("diskfull")
    try:
        import resource
    except ImportError:
        c.notes.append("skipped: no resource limits on this platform")
        return c

    #  The first .vts is about 24 kB, the .dat 9 kB: a 16 kB limit on the
    #  size of a file lets the input copy, JobInfo.txt and the .dat through
    #  and stops the .vts.
    def limit():
        resource.setrlimit(resource.RLIMIT_FSIZE, (16384, 16384))

    write(os.path.join(wd, "a.toml"), small("full"))
    rc, out, err = run(exe, ["a.toml"], wd, preexec=limit)
    d = os.path.join(wd, "full")
    fatal = open(os.path.join(d, "FatalErrors.txt"), encoding="utf-8").read() if os.path.isdir(d) else ""
    c.expect(rc == 1 and "Could not write" in err and "field00000.vts" in err and "field00000.vts" in fatal,
             "file size limit: rc=%d %s" % (rc, short(err)))
    left = sorted(os.listdir(d)) if os.path.isdir(d) else []
    c.expect(not [p for p in left if re.match(r"field.*\.(vts|state)$", p)],
             "file size limit: a field file or .state left: %s" % left)
    c.expect(not [p for p in left if p.endswith(".tmp")], "file size limit: .tmp left: %s" % left)
    no_cut_files(c, d, "file size limit")
    if os.path.isfile(os.path.join(d, "JobInfo.txt")):
        c.expect("abort" in open(os.path.join(d, "JobInfo.txt"), encoding="utf-8").read(), "file size limit: JobInfo.txt not 'abort'")

    if os.path.exists("/dev/full"):
        d = os.path.join(wd, "devfull")
        os.makedirs(d)
        os.symlink("/dev/full", os.path.join(d, "History.dat"))
        write(os.path.join(wd, "b.toml"), small("devfull", {"Output": {"Overwrite": "true"}}))
        rc, out, err = run(exe, ["b.toml"], wd)
        c.expect(rc == 1 and "History.dat" in err and "No space left" in err,
                 "History.dat on /dev/full: rc=%d %s" % (rc, short(err)))
    return c


#--------------------------------------------------------------------------
#  Interruption. The cases are those of restart_roundtrip.py --chain: wm1d
#  (adaptive dt, evaporation, every body force, Marangoni, a probe) and, for
#  SIGINT, cz1d (the ConcZ model).
def chain_input(job, restart=None, chain="wm1d"):
    CHAIN = rr.CHAIN[chain]
    s = {k: dict(x) for k, x in CHAIN["sect"].items()}
    s.setdefault("", {})["JobName"] = '"%s"' % job
    if restart:
        s["Restart"] = restart
    base = open(os.path.join(CASES, CHAIN["case"] + ".toml"), encoding="utf-8").read()
    return rr.toml_edit(base, s, rr.PROBE + CHAIN.get("append", ""))


def straight_run(exe, wd, chain="wm1d"):
    d = os.path.join(wd, "U_" + chain)
    write(os.path.join(d, "in.toml"), chain_input("run", chain=chain))
    rc, out, err = run(exe, ["in.toml"], d)
    return rc, os.path.join(d, "run"), err


def stop_midway(exe, cwd, sig, wait_for, delay=0.0):
    """Start a run, wait until the file wait_for appears, wait delay
    seconds more, send sig; return (rc, stdout, stderr, was it still
    running)."""
    e = dict(os.environ, OMP_NUM_THREADS="1", OMP_WAIT_POLICY="PASSIVE")
    p = subprocess.Popen([exe, "in.toml"], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=e)
    t0 = time.time()
    while not os.path.exists(wait_for) and p.poll() is None and time.time() - t0 < 300:
        time.sleep(0.005)
    if delay:
        time.sleep(delay)
    running = p.poll() is None
    if running:
        p.send_signal(sig)
    out, err = p.communicate(timeout=300)
    return p.returncode, out.decode("utf-8", "replace"), err.decode("utf-8", "replace"), running


def continue_and_compare(c, exe, cwd, d, ud, what, chain="wm1d"):
    """Restart from the last field file of d with its .state, into d, and
    compare everything with the straight run ud."""
    ent = [(t, f) for t, f in rr.pvd_entries(d) if os.path.exists(os.path.join(d, f[:-4] + ".state"))]
    if not c.expect(ent, "%s: no field file with its .state to restart from" % what):
        return
    t, f = ent[-1]
    stop_times = {t}
    idx = int(re.search(r"(\d+)\.vts$", f).group(1))
    restart = {"Enable": "true", "ASCII": "false", "FileName": '"run/%s"' % f,
               "RestartTime": t, "FileIndex": str(idx + 1)}
    write(os.path.join(cwd, "in2.toml"), chain_input("run", restart, chain))
    rc, out, err = run(exe, ["in2.toml"], cwd)
    if not c.expect(rc == 0 and "continuing exactly" in out and "taken to the bit" in out,
                    "%s: the continued run: rc=%d %s" % (what, rc, short(err or out[-300:]))):
        return
    no_cut_files(c, d, what + ", continued")

    def is_stop(tk):
        return any(abs(float(tk) - float(e)) <= 5e-6 * abs(float(e)) for e in stop_times)

    ue = rr.pvd_entries(ud)
    ce = dict(rr.pvd_entries(d))
    missing = [t for t, _ in ue if t not in ce]
    surplus = [t for t in ce if t not in dict(ue) and not is_stop(t)]
    c.expect(not missing, "%s: field times missing: %s" % (what, missing[:4]))
    c.expect(not surplus, "%s: field times the straight run has not: %s" % (what, surplus[:4]))
    for t, f in ue:
        if t in ce:
            a = open(os.path.join(ud, f[:-4] + ".dat"), encoding="utf-8").read()
            b = open(os.path.join(d, ce[t][:-4] + ".dat"), encoding="utf-8").read()
            c.expect(a == b, "%s: t = %s: field*.dat differs" % (what, t))
            c.expect(rr.vts_payload(os.path.join(ud, f)) == rr.vts_payload(os.path.join(d, ce[t])),
                     "%s: t = %s: .vts arrays differ" % (what, t))
    for fn in ("History.dat", "CalcInfo.dat", "Probe00.dat"):
        ur, _ = rr.timeline_rows(os.path.join(ud, fn))
        cr, corder = rr.timeline_rows(os.path.join(d, fn))
        c.expect(len(corder) == len(set(corder)), "%s: %s: a time written twice" % (what, fn))
        diff = [tk for tk in ur if tk in cr and cr[tk] != ur[tk]]
        miss = [tk for tk in ur if tk not in cr]
        surp = [tk for tk in cr if tk not in ur and not is_stop(tk)]
        c.expect(not diff and not miss and not surp, "%s: %s: %d rows differ, %d missing, %d extra"
                 % (what, fn, len(diff), len(miss), len(surp)))


def case_sigint(exe, wd, ud):
    c = Checks("sigint")
    for chain, sig, name in (("wm1d", signal.SIGINT, "SIGINT"), ("cz1d", signal.SIGTERM, "SIGTERM")):
        u = ud if chain == "wm1d" else None
        if u is None:
            rc, u, err = straight_run(exe, wd, chain)
            if not c.expect(rc == 0, "%s: straight run: rc=%d %s" % (chain, rc, short(err))):
                continue
        cwd = os.path.join(wd, "C_" + chain)
        write(os.path.join(cwd, "in.toml"), chain_input("run", chain=chain))
        d = os.path.join(cwd, "run")
        rc, out, err, running = stop_midway(exe, cwd, sig, os.path.join(d, "field00001.state"))
        if not running:
            c.notes.append("%s: the run ended before the signal was sent" % chain)
            continue
        c.expect(rc == 128 + sig and "Stopped by " + name in out,
                 "%s: rc=%d, %s" % (chain, rc, short(out[-300:] + err)))
        info = open(os.path.join(d, "JobInfo.txt"), encoding="utf-8").read()
        c.expect("interrupted (%s)" % name in info, "%s: JobInfo.txt does not say interrupted" % chain)
        no_cut_files(c, d, name + " " + chain)
        c.expect(not glob.glob(os.path.join(d, "*.tmp")), "%s %s: .tmp files left" % (name, chain))
        #  the state reached is written: the last entry of the .pvd has its .state
        ent = rr.pvd_entries(d)
        c.expect(ent and os.path.exists(os.path.join(d, ent[-1][1][:-4] + ".state")),
                 "%s %s: the state reached is not written" % (name, chain))
        continue_and_compare(c, exe, cwd, d, u, name + " " + chain, chain)
    return c


def case_sigkill(exe, wd, ud):
    c = Checks("sigkill")
    if not hasattr(signal, "SIGKILL"):
        c.notes.append("skipped: no SIGKILL on this platform")
        return c
    #  at a moment that varies from one try to the next: right after a
    #  field file appears, when the next outputs are being written
    for k, delay in enumerate((0.0, 0.002, 0.01)):
        cwd = os.path.join(wd, "K%d" % k)
        write(os.path.join(cwd, "in.toml"), chain_input("run"))
        d = os.path.join(cwd, "run")
        rc, out, err, running = stop_midway(exe, cwd, signal.SIGKILL, os.path.join(d, "field00001.vts"), delay)
        if not running:
            c.notes.append("try %d: the run ended before the signal" % k)
            continue
        c.expect(rc == -signal.SIGKILL, "try %d: rc=%d" % (k, rc))
        no_cut_files(c, d, "SIGKILL try %d" % k)
        continue_and_compare(c, exe, cwd, d, ud, "SIGKILL try %d" % k)
    return c


def case_paths(exe, wd):
    c = Checks("paths")
    #  an input file in a directory with spaces and Japanese characters, run
    #  from a sibling directory through a relative path, writing into a
    #  JobName with both, and missing parents
    indir = os.path.join(wd, "入力 ファイル")
    job = "結果 dir/ラン 1"
    write(os.path.join(indir, "ケース 1.toml"), small(job))
    os.makedirs(os.path.join(wd, "work"))
    rc, out, err = run(exe, [os.path.join("..", "入力 ファイル", "ケース 1.toml")], os.path.join(wd, "work"))
    d = os.path.join(wd, "work", "結果 dir", "ラン 1")
    c.expect(rc == 0 and os.path.isfile(os.path.join(d, "ラン 1.pvd"))
             and os.path.isfile(os.path.join(d, "ケース 1.toml")) and os.path.isfile(os.path.join(d, "field00004.state")),
             "non-ASCII paths: rc=%d %s %s" % (rc, short(err), os.listdir(d) if os.path.isdir(d) else ""))

    #  a restart from there, by a FileName with spaces and Japanese characters
    ent = rr.pvd_entries(d)
    if c.expect(ent, "non-ASCII paths: no .pvd entries"):
        t, f = ent[-2]
        idx = int(re.search(r"(\d+)\.vts$", f).group(1))
        write(os.path.join(wd, "work", "再開.toml"), small("再開 run", {"Restart": {
            "Enable": "true", "ASCII": "true", "FileName": '"%s/%s"' % (job, f[:-4] + ".dat"),
            "RestartTime": t, "FileIndex": str(idx + 1)}}))
        rc, out, err = run(exe, ["再開.toml"], os.path.join(wd, "work"))
        c.expect(rc == 0 and "continuing exactly" in out, "non-ASCII Restart.FileName: rc=%d %s"
                 % (rc, short(err or out[-300:])))

    #  a long path: the result directory nested past 300 characters (not on
    #  Windows, where paths past 260 characters need a system setting)
    if POSIX:
        long_job = "/".join(["d%02d_" % k + "x" * 40 for k in range(8)])
        write(os.path.join(wd, "long.toml"), small(long_job))
        rc, out, err = run(exe, [os.path.abspath(os.path.join(wd, "long.toml"))], wd)
        c.expect(rc == 0 and os.path.isfile(os.path.join(wd, long_job, "field00004.vts")),
                 "long path (%d characters): rc=%d %s"
                 % (len(os.path.abspath(os.path.join(wd, long_job))), rc, short(err)))
    else:
        c.notes.append("long path not tried (Windows)")

    #  an input file through a symbolic link, into a result directory that is one
    if POSIX:
        os.makedirs(os.path.join(wd, "real_out"))
        os.symlink("real_out", os.path.join(wd, "linked_out"))
        write(os.path.join(wd, "real.toml"), small("linked_out", {"Output": {"Overwrite": "true"}}))
        os.symlink("real.toml", os.path.join(wd, "link.toml"))
        rc, out, err = run(exe, ["link.toml"], wd)
        c.expect(rc == 0 and os.path.isfile(os.path.join(wd, "real_out", "linked_out.pvd")),
                 "symbolic links: rc=%d %s" % (rc, short(err)))
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--only", nargs="*", help="run only these cases")
    args = ap.parse_args()
    exe = os.path.abspath(args.exe)
    wd = os.path.abspath(args.workdir)
    if os.path.isdir(wd):
        for root, dirs, _ in os.walk(wd):
            for x in dirs:
                try:
                    os.chmod(os.path.join(root, x), 0o755)
                except OSError:
                    pass
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)

    cases = [("cli", case_cli), ("input", case_input), ("resultdir", case_resultdir),
             ("restart", case_restart), ("diskfull", case_diskfull), ("paths", case_paths)]
    signals = [("sigint", case_sigint), ("sigkill", case_sigkill)]
    bad = 0
    ud = None
    for name, fn in cases + signals:
        if args.only and name not in args.only:
            continue
        cwd = os.path.join(wd, name)
        os.makedirs(cwd)
        if fn in (case_sigint, case_sigkill):
            if not POSIX:
                print("skip %-10s no POSIX signals" % name)
                continue
            if ud is None:
                rc, ud, err = straight_run(exe, wd)
                if rc != 0:
                    print("FAIL %-10s straight run: rc=%d %s" % (name, rc, short(err)))
                    bad += 1
                    break
            c = fn(exe, cwd, ud)
        else:
            c = fn(exe, cwd)
        if c.problems:
            bad += 1
            print("FAIL %-10s %s" % (name, "\n                ".join(c.problems)))
        else:
            print("ok   %-10s %s" % (name, "; ".join(c.notes)))
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
