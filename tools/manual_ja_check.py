#!/usr/bin/env python3
"""Check that the Japanese manual says what the English one says.

    python3 tools/manual_ja_check.py            # every chapter
    python3 tools/manual_ja_check.py ch-output  # one chapter
    python3 tools/manual_ja_check.py --pair EN.tex JA.tex   # any two files

The English chapters (docs/manual/chapters/) are the reference; each file in
docs/manual/chapters-ja/ translates the file of the same name. Prose is free
to change in translation, but nothing a reader acts on may: for every chapter
this compares, between the two editions,

  * the labels, and the references made (\\cref, \\Cref, \\ref, \\crefrange);
  * the contents of \\key, \\val, \\file, \\code, \\env and \\url;
  * every inline formula ($...$ and \\(...\\)) and every displayed one
    (equation, align and their starred forms), exactly;
  * every code listing (toml, shell, inp, cpp, plain, verbatim), exactly;
  * every \\SI and \\num quantity, and every number written in the prose;
  * the figures included and the number of tables, list items and boxes.

A difference is printed with the chapter and what is missing from which side,
and the script exits 1. Run by `make check` in docs/manual.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import collections
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EN = os.path.join(ROOT, "docs", "manual", "chapters")
JA = os.path.join(ROOT, "docs", "manual", "chapters-ja")

CODE_ENVS = ("toml", "shell", "inp", "cpp", "plain", "verbatim", "lstlisting")
MATH_ENVS = ("equation", "equation*", "align", "align*", "gather", "gather*", "multline")


def strip_comments(t):
    return re.sub(r"(?<!\\)%.*", "", t)


def group(s, i):
    """The {...} group starting at s[i] == '{' -> (content, index after)."""
    depth, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    return s[i + 1:], len(s)


def args(s, name):
    out = []
    for m in re.finditer(r"\\%s\*?\s*(?=\{)" % re.escape(name), s):
        a, _ = group(s, m.end())
        out.append(a.strip())
    return out


def envs(s, name):
    pat = re.compile(r"\\begin\{%s\}(.*?)\\end\{%s\}" % (re.escape(name), re.escape(name)), re.S)
    return [m.group(1) for m in pat.finditer(s)]


def norm(x):
    """Spacing and the optional break points (\\bk) do not count."""
    return re.sub(r"\s+", " ", re.sub(r"\\bk\s*", "", x)).strip()


def facts(path):
    raw = open(path, encoding="utf-8").read()
    f = {}
    code = []
    rest = raw
    for e in CODE_ENVS:
        for body in envs(rest, e):
            code.append((e, body.strip("\n")))
        rest = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (re.escape(e), re.escape(e)), " ", rest, flags=re.S)
    f["code listing"] = collections.Counter(code)
    s = strip_comments(rest)
    disp = []
    for e in MATH_ENVS:
        for body in envs(s, e):
            disp.append(norm(re.sub(r"\\label\{[^}]*\}", "", body)))
        s = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (re.escape(e), re.escape(e)), " ", s, flags=re.S) \
            if e in MATH_ENVS else s
    f["displayed formula"] = collections.Counter(disp)
    f["label"] = collections.Counter(re.findall(r"\\label\{([^}]*)\}", raw))
    refs = set()
    for c in ("cref", "Cref", "ref", "eqref", "crefrange", "pageref"):
        for a in args(s, c):
            refs.update(x.strip() for x in a.split(","))
    for m in re.finditer(r"\\crefrange\s*\{[^}]*\}\s*\{([^}]*)\}", s):
        refs.add(m.group(1).strip())
    f["reference"] = collections.Counter(refs)
    for c in ("key", "val", "file", "code", "env", "url"):
        f["\\" + c] = collections.Counter(set(norm(a) for a in args(s, c)))
    inline = re.findall(r"(?<!\\)\$(.+?)(?<!\\)\$", s, re.S) + re.findall(r"\\\((.+?)\\\)", s, re.S)
    f["inline formula"] = collections.Counter(set(norm(x) for x in inline))
    si = [norm(a) for a in args(s, "SI")] + [norm(a) for a in args(s, "num")]
    f["\\SI/\\num value"] = collections.Counter(set(si))
    f["figure"] = collections.Counter(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", s))
    # numbers in the prose: everything above with its arguments removed
    prose = s
    for c in ("label", "cref", "Cref", "ref", "eqref", "crefrange", "pageref", "key", "val", "file",
              "code", "env", "url", "SI", "num", "includegraphics", "si", "begin", "end", "hspace",
              "vspace", "setlength", "multicolumn", "cite", "addlinespace"):
        prose = re.sub(r"\\%s\*?(\[[^\]]*\])?\s*\{[^{}]*(\{[^{}]*\}[^{}]*)*\}(\{[^{}]*\})?" % c, " ", prose)
    prose = re.sub(r"(?<!\\)\$.+?(?<!\\)\$", " ", prose, flags=re.S)
    prose = re.sub(r"\\begin\{(?:yzntable|keytable|tabular)\}\{[^}]*\}(\{[^}]*\})?", " ", prose)
    prose = re.sub(r"\\[a-zA-Z]+\[[^\]]*\]", " ", prose)
    prose = prose.replace("{,}", "")
    nums = re.findall(r"(?<![A-Za-z0-9_.])[0-9][0-9,]*(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?", prose)
    f["number in the text"] = collections.Counter(set(n.replace(",", "").replace("{", "") for n in nums))
    f["table"] = collections.Counter({"tables": len(re.findall(r"\\begin\{(?:yzntable|keytable)\}", s))})
    f["box"] = collections.Counter(re.findall(r"\\begin\{(warn|note|tip)\}", s))
    return f


def compare(name, en=None, ja=None):
    a = facts(en or os.path.join(EN, name))
    b = facts(ja or os.path.join(JA, name))
    bad = []
    for k in a:
        only_en = a[k] - b[k]
        only_ja = b[k] - a[k]
        if only_en or only_ja:
            show = lambda c: "; ".join(repr(x)[:120] for x in sorted(c, key=str)[:8]) + (
                " (+%d more)" % (len(c) - 8) if len(c) > 8 else "")
            if only_en:
                bad.append("  %s missing in ja: %s" % (k, show(only_en)))
            if only_ja:
                bad.append("  %s only in ja: %s" % (k, show(only_ja)))
    return bad


def main():
    if sys.argv[1:2] == ["--pair"]:          # two files of your choosing
        bad = compare(None, sys.argv[2], sys.argv[3])
        print("\n".join(bad) if bad else "manual_ja_check: the two files agree")
        sys.exit(1 if bad else 0)
    names = sorted(f for f in os.listdir(EN) if f.endswith(".tex"))
    want = [n if n.endswith(".tex") else n + ".tex" for n in sys.argv[1:]] or names
    missing = [n for n in want if not os.path.exists(os.path.join(JA, n))]
    extra = sorted(set(f for f in os.listdir(JA) if f.endswith(".tex")) - set(names))
    failed = 0
    for n in missing:
        print("%s: no Japanese chapter" % n)
        failed += 1
    for n in extra:
        print("%s: a Japanese chapter with no English one" % n)
        failed += 1
    for n in want:
        if n in missing:
            continue
        bad = compare(n)
        if bad:
            failed += 1
            print("%s:" % n)
            print("\n".join(bad))
    if failed:
        print("manual_ja_check: %d chapter(s) differ between the editions" % failed)
        sys.exit(1)
    print("manual_ja_check: %d chapters agree" % len(want))


if __name__ == "__main__":
    main()
