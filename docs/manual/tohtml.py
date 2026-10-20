#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the HTML edition of the Yuzen manual from the same LaTeX sources.

    python3 docs/manual/tohtml.py            # writes yuzen-manual.html
    python3 docs/manual/tohtml.py --lang ja  # writes yuzen-manual-ja.html
    make html

The point of this script is that there is exactly one source of truth. The
PDF and the HTML are both built from chapters/*.tex, so they cannot drift.
It understands the markup that STYLE.md defines and nothing else: if a chapter
uses a construct not listed here, the build fails loudly rather than dropping
it silently.

Mathematics is pre-rendered to SVG with MathJax at build time, so the finished
page needs no network and no JavaScript to display an equation correctly. The
JavaScript that is there does one thing: full-text search.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))

#  The edition: English (chapters/) or Japanese (chapters-ja/). Everything the
#  converter itself writes into the page comes from WORDS.
LANG = "ja" if "--lang" in sys.argv and sys.argv[sys.argv.index("--lang") + 1:][:1] == ["ja"] else "en"
MAIN_TEX = "yuzen-manual-ja.tex" if LANG == "ja" else "yuzen-manual.tex"
OUT_HTML = "yuzen-manual-ja.html" if LANG == "ja" else "yuzen-manual.html"
CH = os.path.join(HERE, "chapters-ja" if LANG == "ja" else "chapters")
WORDS = {
    "en": dict(lang="en", title="Yuzen User Manual", sub="User Manual &middot; Version 1.0",
               desc="Yuzen: a thin liquid film flow solver based on the lubrication approximation. "
                    "Full manual: input reference, theory, numerical method and verification.",
               search="Search the manual…  (press /)", hint="Full text, offline. Esc clears.",
               none="Nothing found.", dark="Dark", light="Light", prev="Previous", next="Next",
               chapter="Chapter", appendix="Appendix", table="Table", figure="Figure",
               warn="Warning", note="Note", tip="Rule of thumb", required="required",
               khead=["Key", "Type", "Default", "Meaning, units, range"],
               other="日本語版", other_href="yuzen-manual-ja.html"),
    "ja": dict(lang="ja", title="Yuzen ユーザーマニュアル", sub="ユーザーマニュアル &middot; Version 1.0",
               desc="Yuzen: 潤滑近似に基づく液膜流れソルバー。入力リファレンス、理論、数値解法、検証を含む完全版マニュアル。",
               search="マニュアルを検索…  (/ キー)", hint="全文検索、オフラインで動作。Esc で消去。",
               none="見つかりませんでした。", dark="ダーク", light="ライト", prev="前へ", next="次へ",
               chapter="第%s章", appendix="付録%s", table="表", figure="図",
               warn="警告", note="注", tip="目安", required="必須",
               khead=["キー", "型", "既定値", "意味・単位・範囲"],
               other="English", other_href="yuzen-manual.html"),
}[LANG]
#  Figures come from two places: the ones drawn from formulas live beside the
#  manual, and the ones drawn from solver output live in docs/figures/ so that
#  the repository README can show them too. The key is the path exactly as
#  \includegraphics writes it.
FIGDIRS = [("figures/", os.path.join(HERE, "figures")),
           ("../figures/", os.path.join(HERE, os.pardir, "figures"))]
#  Scratch space for the build. It sits beside this script rather than in a
#  fixed absolute path, so the script works wherever the repository is checked
#  out. Everything in it is a by-product; .gitignore covers it.
BUILD = os.path.join(HERE, ".build")


#  Files the renderer loads. Checking these is enough: a damaged install
#  never damages only the files nobody reads.
MATHJAX_FILES = ("js/mathjax.js", "js/util/Options.js", "js/input/tex.js",
                 "js/output/svg.js", "js/adaptors/liteAdaptor.js",
                 "js/handlers/html.js")


def check_mathjax(pkg, node=None):
    """Return a description of what is wrong with a mathjax-full install.

    Returns None when it looks intact. The manifest is parsed here in Python
    rather than through node, because a damaged package.json breaks node's own
    module resolver before any of our code could report on it. The JavaScript
    is then handed to node to parse, which is the only reliable test: counting
    brackets calls valid minified files broken, since a brace inside a string
    or a regular expression counts all the same.
    """
    manifest = os.path.join(pkg, "package.json")
    if not os.path.isfile(manifest):
        return "%s is missing" % manifest
    try:
        json.load(open(manifest, encoding="utf-8"))
    except Exception as e:
        return "%s is not valid JSON (%s)" % (manifest, e)

    missing = []
    for rel in MATHJAX_FILES:
        f = os.path.join(pkg, rel)
        if not os.path.isfile(f):
            missing.append("%s is missing" % f)
        elif os.path.getsize(f) == 0:
            missing.append("%s is empty" % f)
    if missing:
        return missing[0]
    if node is None:
        return None

    #  Parse without executing: a syntax error means the file was written
    #  incompletely, which is how a failed install shows itself.
    probe = (
        "const fs=require('fs'),vm=require('vm'),path=require('path');"
        "for (const rel of process.argv.slice(2)) {"
        "  const f=path.join(%s,rel);"
        "  try { new vm.Script(fs.readFileSync(f,'utf8'), {filename:f}); }"
        "  catch (e) { console.log(f+' :: '+e.message); process.exit(3); } }"
        % json.dumps(pkg))
    r = subprocess.run([node, "-e", probe, "--"] + list(MATHJAX_FILES),
                       capture_output=True, text=True)
    if r.returncode == 3:
        f, _, why = r.stdout.strip().partition(" :: ")
        return "%s does not parse as JavaScript (%s)" % (f, why)
    return None


def report_install():
    """Say which mathjax-full install would be used, and whether it is intact.

    Run as:  python3 tohtml.py --check-install   (or: make doctor)
    """
    node = shutil.which("node")
    print("node:        %s" % (node or "NOT FOUND"))
    seen = []
    for d in (HERE, os.path.dirname(HERE), os.path.dirname(os.path.dirname(HERE)),
              os.getcwd()):
        cand = os.path.join(d, "node_modules")
        if cand not in seen:
            seen.append(cand)
    for cand in os.environ.get("NODE_PATH", "").split(os.pathsep):
        if cand and cand not in seen:
            seen.append(cand)
    found = False
    for cand in seen:
        pkg = os.path.join(cand, "mathjax-full")
        if not os.path.isdir(pkg):
            continue
        found = True
        bad = check_mathjax(pkg, node)
        print("mathjax-full: %s\n  %s" % (pkg, bad or "intact"))
    if not found:
        print("mathjax-full: not installed in any of:")
        for cand in seen:
            print("  %s" % cand)
    return 0


def node_environment():
    """Locate node and mathjax-full, or explain what is missing and stop.

    mathjax-full is resolved through NODE_PATH rather than the working
    directory, so it can be installed here, at the top of the repository, or
    globally, and this script finds it either way.
    """
    node = shutil.which("node")
    if node is None:
        sys.exit(
            "tohtml.py: node was not found.\n"
            "  The HTML edition pre-renders every formula to SVG with MathJax,\n"
            "  which runs under Node.\n"
            "    Debian/Ubuntu:  sudo apt install nodejs npm\n"
            "    macOS:          brew install node\n"
            "  then, in docs/manual:  npm install mathjax-full\n"
            "  The built PDF and HTML are attached to every release, so there\n"
            "  is no need to install this merely to read the manual.")

    #  Candidate node_modules directories, nearest first. An install made on
    #  the Linux filesystem and pointed at with NODE_PATH is checked too, since
    #  that is the way out when the repository sits on a Windows drive.
    roots, seen = [], []
    for d in (HERE, os.path.dirname(HERE), os.path.dirname(os.path.dirname(HERE)),
              os.getcwd()):
        cand = os.path.join(d, "node_modules")
        if cand not in seen:
            seen.append(cand)
    for cand in os.environ.get("NODE_PATH", "").split(os.pathsep):
        if cand and cand not in seen:
            seen.append(cand)
    for cand in seen:
        if os.path.isdir(os.path.join(cand, "mathjax-full")):
            roots.append(cand)

    if not roots:
        sys.exit(
            "tohtml.py: node is installed but the mathjax-full package is not.\n"
            "  In docs/manual:  npm install mathjax-full\n"
            "  (a global install works too). The built PDF and HTML are\n"
            "  attached to every release if you only want to read the manual.")

    #  Use the first install that is actually intact. A damaged one is worse
    #  than a missing one: node fails deep inside MathJax, or even inside its
    #  own module resolver while reading package.json, and the message reads
    #  like a bug in MathJax or in this script rather than a broken file.
    #  npm onto a Windows drive from WSL (/mnt/c/...) is how this happens, and
    #  it can recur however many times the install is repeated.
    good, broken = None, []
    for cand in roots:
        bad = check_mathjax(os.path.join(cand, "mathjax-full"), node)
        if bad is None:
            good = cand
            break
        broken.append(bad)
    if good is None:
        sys.exit(
            "tohtml.py: the installed mathjax-full package is damaged --\n"
            "  %s\n"
            "  npm writes files incompletely when it installs onto a Windows\n"
            "  drive from WSL, which is where this package sits. Reinstalling\n"
            "  in the same place often fails the same way.\n"
            "\n"
            "  Install it on the Linux filesystem instead and point node at it:\n"
            "    (cd ~ && npm install mathjax-full)\n"
            "    NODE_PATH=~/node_modules make html\n"
            "\n"
            "  Or do not build it at all: the PDF and the HTML edition are\n"
            "  attached to every release, and CI rebuilds both on a tag."
            % "\n  ".join(broken))

    env = dict(os.environ)
    env["NODE_PATH"] = os.pathsep.join(
        [good] + ([env["NODE_PATH"]] if env.get("NODE_PATH") else []))
    return node, env, os.path.join(good, "mathjax-full")

# ---------------------------------------------------------------------------
#  1. read the chapter order out of the main file, so the two never disagree
# ---------------------------------------------------------------------------

def read_structure():
    main = open(os.path.join(HERE, MAIN_TEX), encoding="utf-8").read()
    body = main.split(r"\begin{document}", 1)[1]
    items, part, appendix = [], None, False
    for m in re.finditer(r"\\(part|include|appendix|backmatter|mainmatter|frontmatter)"
                         r"(?:\{([^}]*)\})?", body):
        cmd, arg = m.group(1), m.group(2)
        if cmd == "part":
            part = arg
        elif cmd == "appendix":
            appendix = True
        elif cmd == "include":
            name = arg.split("/")[-1]
            if name in ("titlepage",):
                continue
            items.append({"file": name, "part": part, "appendix": appendix})
    return items


# ---------------------------------------------------------------------------
#  2. maths: collect every snippet, render once with MathJax, cache by text
# ---------------------------------------------------------------------------

MACROS = r"""
\newcommand{\As}{\mathrm{As}}
\newcommand{\Bo}{\mathrm{Bo}}
\newcommand{\Dpn}{\mathrm{Dp}}
\newcommand{\Ma}{\mathrm{Ma}}
\newcommand{\Spn}{\mathrm{Sp}}
\newcommand{\Ev}{\mathrm{Ev}}
\newcommand{\Pen}{\mathrm{Pe}}
\newcommand{\Df}{\mathrm{Df}}
\newcommand{\Mr}{\mathrm{Mr}}
\newcommand{\Ca}{\mathrm{Ca}}
\newcommand{\Pgx}{\mathrm{Pgx}}
\newcommand{\Pgy}{\mathrm{Pgy}}
\newcommand{\Pvec}{\mathbf{Pg}}
\newcommand{\Qvec}{\mathbf{Q}}
\newcommand{\uvec}{\mathbf{u}}
\newcommand{\er}{\mathbf{e}_r}
\newcommand{\hstar}{h_{*}}
\newcommand{\Hstar}{H_{*}}
\newcommand{\Hmax}{H_{\mathrm{max}}}
\newcommand{\Hf}{H_{\mathrm{f}}}
\newcommand{\thetae}{\theta_{e}}
\newcommand{\thetacap}{\theta_{\mathrm{cap}}}
\newcommand{\Pihat}{\hat{\Pi}}
\newcommand{\epsmach}{\varepsilon_{\mathrm{mach}}}
\newcommand{\dtexp}{\Delta t_{\mathrm{explicit}}}
\newcommand{\Rey}{\mathrm{Re}}
"""

_math_jobs = []          # (display:bool, tex:str)
_math_index = {}         # (display, tex) -> slot


def math_slot(tex, display):
    key = (display, tex)
    if key not in _math_index:
        _math_index[key] = len(_math_jobs)
        _math_jobs.append({"tex": tex, "display": display})
    return _math_index[key]


def render_all_math():
    """One node process for the whole document; ~400 snippets in a second."""
    if not _math_jobs:
        return []
    script = r"""
const fs = require('fs');
const path = require('path');
//  The package directory is passed in rather than resolved by name: node
//  searches the nearest node_modules first, so a damaged install beside the
//  repository would win over the good one we actually checked.
const MJ = process.argv[5];
const mj = p => require(path.join(MJ, p));
const {mathjax} = mj('js/mathjax.js');
const {TeX} = mj('js/input/tex.js');
const {SVG} = mj('js/output/svg.js');
const {liteAdaptor} = mj('js/adaptors/liteAdaptor.js');
const {RegisterHTMLHandler} = mj('js/handlers/html.js');
mj('js/input/tex/base/BaseConfiguration.js');
mj('js/input/tex/ams/AmsConfiguration.js');
mj('js/input/tex/newcommand/NewcommandConfiguration.js');
mj('js/input/tex/configmacros/ConfigMacrosConfiguration.js');
const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const jobs  = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const macros = fs.readFileSync(process.argv[3], 'utf8');
const tex = new TeX({packages: ['base', 'ams', 'newcommand', 'configmacros']});
const svg = new SVG({fontCache: 'local'});
const doc = mathjax.document('', {InputJax: tex, OutputJax: svg});
doc.convert(macros, {display: true});          // define the macros once
const out = jobs.map(j => {
  try {
    return {ok: true, svg: adaptor.outerHTML(
      doc.convert(j.tex, {display: j.display, em: 16, ex: 8, containerWidth: 800}))};
  } catch (e) {
    return {ok: false, err: String(e && e.message || e)};
  }
});
fs.writeFileSync(process.argv[4], JSON.stringify(out));
"""
    node, env, pkg = node_environment()
    tmp = os.path.join(BUILD, "mathbuild")
    os.makedirs(tmp, exist_ok=True)
    for name, data in (("jobs.json", json.dumps(_math_jobs)),
                       ("macros.tex", MACROS),
                       ("run.js", script)):
        open(os.path.join(tmp, name), "w", encoding="utf-8").write(data)
    r = subprocess.run([node, os.path.join(tmp, "run.js"),
                        os.path.join(tmp, "jobs.json"),
                        os.path.join(tmp, "macros.tex"),
                        os.path.join(tmp, "out.json"),
                        pkg],
                       cwd=tmp, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("MathJax failed:\n" + r.stderr[-3000:])
    res = json.load(open(os.path.join(tmp, "out.json")))
    bad = [(_math_jobs[i]["tex"], x.get("err")) for i, x in enumerate(res) if not x["ok"]]
    # MathJax does not throw on an unknown macro: it renders an error node in
    # place of the formula. That is a silent mistranslation, so treat it as a
    # failure too.
    bad += [(_math_jobs[i]["tex"],
             re.search(r'data-mjx-error="([^"]*)"', x["svg"]).group(1))
            for i, x in enumerate(res)
            if x["ok"] and "data-mjx-error" in x["svg"]]
    if bad:
        for t, e in bad[:20]:
            print("  MATH FAILED:", t[:90], "->", e, file=sys.stderr)
        sys.exit("%d maths snippets failed to render" % len(bad))
    return [x["svg"] for x in res]


# ---------------------------------------------------------------------------
#  3. inline markup
# ---------------------------------------------------------------------------

SEMANTIC = {"key": "k", "val": "v", "file": "f", "code": "c", "env": "e"}

SIUNITS = {
    r"\micro\metre": "µm", r"\milli\metre": "mm", r"\metre": "m",
    r"\milli\metre\per\second": "mm/s",
    # standalone \per: siunitx (per-mode=symbol) supplies the numerator 1
    r"\milli\second": "ms", r"\micro\second": "µs", r"\nano\second": "ns", r"\cubic\metre": "m³",
    r"\second": "s", r"\per\second\squared": "1/s²", r"\per\second": "1/s",
    r"\metre\per\second\squared": "m/s²", r"\metre\per\second": "m/s",
    r"\pascal\second": "Pa·s", r"\percent": "%", r"\giga\hertz": "GHz",
    r"\giga\byte": "GB", r"\mega\byte": "MB", r"\kilo\byte": "kB",
    r"\squared": "²", r"\kilo\gram": "kg", r"\newton\per\metre": "N/m",
    r"\kilogram\per\metre\cubed": "kg/m³", r"\gram\per\mole": "g/mol",
    r"\metre\squared\per\second": "m²/s", r"\metre\cubed\per\metre": "m³/m", r"\metre\cubed": "m³",
    r"\radian\per\second": "rad/s", r"\radian": "rad",
    r"\kilogram": "kg", r"\gram": "g", r"\mole": "mol", r"\cubed": "³",
    r"\hour": "h",
}

#  Text-mode accents, as combining characters: \'e is e + U+0301.
ACCENTS = {"'": "́", "`": "̀", "^": "̂", '"': "̈",
           "~": "̃", "=": "̄"}

SYMBOLS = [
    (r"\textcopyright{}", "©"), (r"\textcopyright", "©"),
    (r"\LaTeX", "LaTeX"), (r"\TeX", "TeX"),
    (r"\dots\ ", "… "), (r"\dots", "…"), (r"\ldots", "…"),
    (r"\textendash", "–"), (r"\textemdash", "—"),
    (r"\,", " "), (r"\;", " "), (r"\ ", " "), (r"\%", "%"),
    (r"\&", "&"), (r"\#", "#"), (r"\_", "_"), (r"\$", "$"),
    (r"\{", "{"), (r"\}", "}"),
    (r"\textasciitilde", "~"), (r"\textbackslash", "\\"),
    (r"\emph", ""),  # handled before this, defensive
]


UNKNOWN_UNITS = set()


def si_unit(u):
    """\\si / \\SI unit argument -> the text siunitx would set."""
    for k, v in sorted(SIUNITS.items(), key=lambda kv: -len(kv[0])):
        u = u.replace(k, v)
    # an unknown unit macro must not be silently flattened into prose
    # ("kilogrampermcubed"); report it like any other unhandled command.
    UNKNOWN_UNITS.update(re.findall(r"\\[a-zA-Z]+", u))
    return u.replace("\\", "")


def math_text(t):
    """\\SI, \\si, \\ang and the semantic code macros also occur inside maths.
    MathJax knows none of them, and it does not raise on an unknown control
    sequence -- it renders an error box -- so expand them to plain TeX here."""
    while True:
        mm = re.search(r"\\(SI|si|ang)\{", t)
        if not mm:
            break
        a, j = take_group(t, mm.end() - 1)
        if mm.group(1) == "SI":
            b, j = take_group(t, j)
            rep = r"%s\,\text{%s}" % (a, si_unit(b)) if b.strip() else a
        elif mm.group(1) == "si":
            rep = r"\text{%s}" % si_unit(a)
        else:
            rep = r"%s^{\circ}" % a
        t = t[:mm.start()] + rep + t[j:]
    for name in SEMANTIC:
        t = replace_cmd(t, name, lambda a: r"\texttt{%s}" % a)
    return t


def split_cells(line):
    """Split one table row at its top-level & into cells.

    An escaped character (\\&, \\{, \\}) is text, not a separator or a brace,
    and a control space (\\ ) that ends a cell is kept as a space, so that the
    strip() applied to each cell does not leave a lone backslash behind.
    """
    cells, depth, cur = [], 0, ""
    for chunk in re.split(r"(\\.|&|\{|\})", line, flags=re.S):
        if chunk == "{":
            depth += 1
        elif chunk == "}":
            depth -= 1
        elif chunk == "&" and depth == 0:
            cells.append(cur)
            cur = ""
            continue
        cur += chunk
    cells.append(cur)
    return [re.sub(r"\\\s*$", " ", c.strip()) for c in cells]


def take_group(s, i):
    """s[i] == '{'; return (content, index just past the matching '}')."""
    assert s[i] == "{", s[max(0, i-30):i+30]
    d, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        if c == "{":
            d += 1
        elif c == "}":
            d -= 1
            if d == 0:
                return s[i+1:j], j+1
        j += 1
    raise ValueError("unbalanced brace")


def cmd_with_arg(s, name):
    """Yield (start, end, arg) for every \\name{...} in s, outermost first."""
    out, i = [], 0
    pat = "\\" + name
    while True:
        k = s.find(pat, i)
        if k < 0:
            break
        j = k + len(pat)
        # not a prefix of a longer command name
        if j < len(s) and (s[j].isalpha()):
            i = k + 1
            continue
        while j < len(s) and s[j] in " \n":
            j += 1
        if j >= len(s) or s[j] != "{":
            i = k + 1
            continue
        arg, end = take_group(s, j)
        out.append((k, end, arg))
        i = end
    return out


def replace_cmd(s, name, fn):
    for k, end, arg in reversed(cmd_with_arg(s, name)):
        s = s[:k] + fn(arg) + s[end:]
    return s


# cleveref is loaded with [capitalise], so every reference is capitalized in
# the PDF, a multi-reference is merged under one plural word, and an equation
# number is parenthesized.
CREF_WORD = {"ch": ("Chapter", "Chapters"), "sec": ("Section", "Sections"),
             "tab": ("Table", "Tables"), "fig": ("Figure", "Figures"),
             "eq": ("Eq.", "Eqs.")}


def ref_word(kind, plural, appendix=False):
    one, many = CREF_WORD.get(kind, ("Section", "Sections"))
    if kind == "ch" and appendix:
        one, many = "Appendix", "Appendices"
    return many if plural else one


def ja_ref(kind, n):
    """第3章, 付録A, 3.2節, 表3.1, 図3.1, 式(3.1): the words of the Japanese PDF."""
    if kind == "ch":
        return ("付録%s" if n[:1].isalpha() else "第%s章") % n
    return {"tab": "表%s", "fig": "図%s", "eq": "式(%s)"}.get(kind, "%s節") % n


def num_key(n):
    """Sort key matching cleveref's ordering: 1 < 1.2 < 2 < A < A.1."""
    out = []
    for part in str(n).split("."):
        out.append((1, part) if part[:1].isalpha() else (0, int(part or 0)))
    return out


def delatex(a):
    """LaTeX escapes inside a code span are literal characters in HTML."""
    # \bk is a discretionary break: like any control word it also swallows the
    # whitespace that follows it, so "a/\bk b" is "a/b", not "a/ b".
    a = re.sub(r"\\bk\s*", "", a)
    for k, v in ((r"\_", "_"), (r"\&", "&"), (r"\#", "#"), (r"\%", "%"),
                 (r"\$", "$"), (r"\{", "{"), (r"\}", "}"),
                 (r"\textbackslash", "\\"), (r"\dots", "…"),
                 (r"\,", "\u202f"), (r"\ ", " ")):
        a = a.replace(k, v)
    return a


def inline(s, ctx):
    """LaTeX inline markup -> HTML. ctx carries the label->number maps."""
    # Footnotes come out first, as raw LaTeX: their bodies are rendered by
    # their own top-level inline() call, so they must not share this one's
    # sentinel store.
    for k, end, arg in reversed(cmd_with_arg(s, "footnote")):
        ctx["footnotes"].append(arg)
        s = s[:k] + "\x03%d\x03" % len(ctx["footnotes"]) + s[end:]

    # maths: protect it behind sentinels
    store = []

    def stash(payload):
        store.append(payload)
        return "\x00%d\x00" % (len(store) - 1)

    # \[ ... \] and $ ... $
    s = re.sub(r"\\\[(.+?)\\\]",
               lambda m: stash(("math", m.group(1).strip(), True)), s, flags=re.S)

    def dollar(m):
        return stash(("math", m.group(1), False))
    s = re.sub(r"(?<!\\)\$(.+?)(?<!\\)\$", dollar, s, flags=re.S)

    # \SI{n}{unit} and \si{unit} become text
    for name, wrap in (("SI", True), ("si", False)):
        while True:
            k = s.find("\\" + name + "{")
            if k < 0:
                break
            a, j = take_group(s, k + len(name) + 1)
            if wrap:
                b, j = take_group(s, j)
                s = s[:k] + stash(("txt", a + " " + si_unit(b))) + s[j:]
            else:
                s = s[:k] + stash(("txt", si_unit(a))) + s[j:]
    # \ang{10}
    while True:
        k = s.find("\\ang{")
        if k < 0:
            break
        a, j = take_group(s, k + 4)
        s = s[:k] + stash(("txt", a + "°")) + s[j:]

    # \texorpdfstring{tex}{pdf} -> the TeX form
    for k, end, arg in reversed(cmd_with_arg(s, "texorpdfstring")):
        j = end
        while j < len(s) and s[j] in " \n":
            j += 1
        if j < len(s) and s[j] == "{":
            _, j = take_group(s, j)
        s = s[:k] + arg + s[j:]
    s = replace_cmd(s, "num", lambda a: a)
    for drop in ("\\noindent", "\\hfill", "\\bigskip", "\\medskip", "\\smallskip",
                 "\\bfseries", "\\sffamily", "\\normalfont", "\\centering",
                 "\\quad", "\\qquad", "\\linebreak", "\\allowbreak",
                 "\\cleardoublepage", "\\clearpage", "\\par "):
        s = s.replace(drop, " " if drop in ("\\quad", "\\qquad") else "")

    # semantic code spans
    for name, cls in SEMANTIC.items():
        s = replace_cmd(s, name, lambda a, c=cls: stash(("span", c, a)))
    s = replace_cmd(s, "reqd", lambda a: a)  # never takes an argument; safety
    s = s.replace("\\reqd", stash(("span", "req", WORDS["required"])))
    s = re.sub(r"\\bk\s*", "", s)
    s = s.replace("\\dash", "–")

    # emphasis
    s = replace_cmd(s, "textbf", lambda a: "\x01b\x01" + a + "\x01/b\x01")
    s = replace_cmd(s, "emph", lambda a: "\x01i\x01" + a + "\x01/i\x01")
    s = replace_cmd(s, "textit", lambda a: "\x01i\x01" + a + "\x01/i\x01")
    s = replace_cmd(s, "textsf", lambda a: a)
    s = replace_cmd(s, "texttt", lambda a: stash(("span", "c", a)))
    s = replace_cmd(s, "textcolor", lambda a: a)
    s = replace_cmd(s, "mbox", lambda a: a)
    s = replace_cmd(s, "url", lambda a: stash(("url", a)))
    s = replace_cmd(s, "href", lambda a: a)

    # cross references
    def cref(arg, cap):
        parts = [p.strip() for p in arg.split(",") if p.strip()]
        seen = []
        for p in parts:
            num = ctx["nums"].get(p)
            if num is None:
                ctx["missing"].add(p)
                num = "?"
            seen.append((p.split(":")[0], p, num))
        kind = seen[0][0]
        try:
            seen.sort(key=lambda x: num_key(x[2]))
        except (ValueError, IndexError):
            pass
        appendix = all(x[2][:1].isalpha() for x in seen)
        # cleveref is loaded with [capitalise], so \cref and \Cref both
        # capitalize; the distinction does not survive into the output.
        if LANG == "ja":
            links = [stash(("ref", lbl, "", ja_ref(k, n))) for k, lbl, n in seen]
            return "と".join(links) if len(links) <= 2 else "、".join(links[:-1]) + "と" + links[-1]
        word = ref_word(kind, len(seen) > 1, appendix)
        links = [stash(("ref", lbl, "", ("(%s)" % n) if k == "eq" else n))
                 for k, lbl, n in seen]
        if len(links) == 1:
            body = links[0]
        elif len(links) == 2:
            body = links[0] + " and " + links[1]
        else:
            body = ", ".join(links[:-1]) + " and " + links[-1]
        return word + "\u00a0" + body

    s = replace_cmd(s, "Cref", lambda a: cref(a, True))
    s = replace_cmd(s, "cref", lambda a: cref(a, False))
    # \crefrange{first}{last} -> "chapters 4 to 9"
    for k, end, a in reversed(cmd_with_arg(s, "crefrange")):
        j = end
        while j < len(s) and s[j] in " \n":
            j += 1
        b, j = take_group(s, j) if (j < len(s) and s[j] == "{") else ("", end)
        kind = a.split(":")[0]
        n1 = ctx["nums"].get(a, "?")
        n2 = ctx["nums"].get(b, "?")
        for lbl, n in ((a, n1), (b, n2)):
            if n == "?":
                ctx["missing"].add(lbl)
        appendix = n1[:1].isalpha() and n2[:1].isalpha()
        if LANG == "ja":
            s = s[:k] + "%s〜%s" % (stash(("ref", a, "", ja_ref(kind, n1))),
                                    stash(("ref", b, "", ja_ref(kind, n2)))) + s[j:]
            continue
        word = ref_word(kind, True, appendix)
        s = s[:k] + "%s\u00a0%s to %s" % (
            word,
            stash(("ref", a, "", n1)),
            stash(("ref", b, "", n2))) + s[j:]
    s = replace_cmd(s, "label", lambda a: "")

    # accented letters (P\'eclet, y\=u, \'{e}) -> precomposed characters
    s = re.sub(r"\\([`'^\"~=])\{?([A-Za-z])\}?",
               lambda m: unicodedata.normalize("NFC", m.group(2) + ACCENTS[m.group(1)]), s)
    # a forced line break (\\ or \\[2pt]) in running text
    s = re.sub(r"\\\\(?:\[[^\]]*\])?\s*", "\x01br\x01", s)

    # leftovers that are pure text
    for a, b in SYMBOLS:
        s = s.replace(a, b)
    s = s.replace("---", "—").replace("--", "–")
    s = s.replace("``", "“").replace("''", "”").replace("~", " ")

    # anything still backslashed is a bug in this converter, not in the source
    leftover = re.findall(r"\\[a-zA-Z]+|\\[`'^\"=]", s)
    if leftover:
        ctx["unhandled"].update(leftover)

    s = html.escape(s)
    s = s.replace("\x01b\x01", "<strong>").replace("\x01/b\x01", "</strong>")
    s = s.replace("\x01i\x01", "<em>").replace("\x01/i\x01", "</em>")
    s = s.replace("\x01br\x01", "<br>")

    # put the stashed pieces back
    def unstash(m):
        item = store[int(m.group(1))]
        kind = item[0]
        if kind == "math":
            slot = math_slot(math_text(item[1]), item[2])
            return '\x02M%d\x02' % slot
        if kind == "txt":
            return html.escape(item[1])
        if kind == "span":
            return '<code class="%s">%s</code>' % (item[1],
                                                   html.escape(delatex(item[2])))
        if kind == "url":
            u = item[1]
            return '<a href="%s" class="ext">%s</a>' % (html.escape(u), html.escape(u))
        if kind == "ref":
            txt = ("%s\u00a0%s" % (item[2], item[3])) if item[2] else str(item[3])
            return '<a href="#%s" class="xref">%s</a>' % (
                html.escape(item[1]), txt)
        raise AssertionError(kind)
    while "\x00" in s:
        s = re.sub(r"\x00(\d+)\x00", unstash, s)
    s = re.sub(r"\x03(\d+)\x03",
               lambda m: '<sup class="fn">%s</sup>' % m.group(1), s)
    return s


# ---------------------------------------------------------------------------
#  4. block structure
# ---------------------------------------------------------------------------

ENV_BOX = {"warn": ("warn", WORDS["warn"]),
           "note": ("note", WORDS["note"]),
           "tip":  ("tip",  WORDS["tip"])}
ENV_CODE = {"shell": "shell", "toml": "toml", "inp": "inp",
            "cpp": "cpp", "plain": "plain"}


def split_env(text, name):
    """Return (before, optarg, body, after) for the first \\begin{name}."""
    b = "\\begin{%s}" % name
    e = "\\end{%s}" % name
    k = text.find(b)
    if k < 0:
        return None
    j = k + len(b)
    opt = None
    if j < len(text) and text[j] == "[":
        d, m = 0, j
        while m < len(text):
            if text[m] == "[":
                d += 1
            elif text[m] == "]":
                d -= 1
                if d == 0:
                    break
            m += 1
        opt = text[j+1:m]
        j = m + 1
    depth, m = 1, j
    while depth:
        nb, ne = text.find(b, m), text.find(e, m)
        if ne < 0:
            raise ValueError("unclosed " + name)
        if 0 <= nb < ne:
            depth += 1
            m = nb + len(b)
        else:
            depth -= 1
            m = ne + len(e)
    return text[:k], opt, text[j:ne], text[m:]


def parse_table(body, kind, ctx):
    """keytable / yzntable -> (caption_html, label, rows)."""
    body = body.strip()
    cap, j = take_group(body, body.index("{"))
    rest = body[j:]
    if kind == "yzntable":            # second argument is the column spec
        while rest and rest[0] in " \n":
            rest = rest[1:]
        _, j2 = take_group(rest, 0)
        rest = rest[j2:]
    label = None
    m = re.search(r"\\label\{([^}]*)\}", cap)
    if m:
        label = m.group(1)
        cap = cap[:m.start()] + cap[m.end():]
    rows = []
    #  Rules and vertical-spacing commands carry no meaning in HTML, so drop them.
    #  (\addlinespace is used to separate groups part-way through a table.)
    rest = (rest.replace("\\midrule", "").replace("\\toprule", "")
                .replace("\\bottomrule", "").replace("\\addlinespace", ""))
    rest = rest.replace("\\bottomrule", "").replace("\\endhead", "")
    for line in re.split(r"\\\\", rest):
        line = line.strip()
        if not line:
            continue
        cells = split_cells(line)
        row = []
        for c in cells:
            c = c.strip()
            #  \multicolumn{n}{spec}{text}: one cell spanning n columns,
            #  carried as (n, html) so that the renderer can give it colspan.
            if c.startswith("\\multicolumn{"):
                n, j = take_group(c, len("\\multicolumn"))
                _, j = take_group(c, j)
                txt, _ = take_group(c, j)
                row.append((int(n), inline(txt.strip(), ctx)))
            else:
                row.append(inline(c, ctx))
        rows.append(row)
    if kind == "keytable":
        header = list(WORDS["khead"])
    else:
        header = rows[0] if rows else []
        rows = rows[1:]
    return inline(cap, ctx), label, header, rows


def split_items(body):
    """Split on \\item, but only at the top level: an \\item belonging to a
    nested list must stay inside its own environment."""
    parts, depth, last, i = [], 0, 0, 0
    tok = re.compile(r"\\begin\{|\\end\{|\\item\b")
    while True:
        m = tok.search(body, i)
        if not m:
            break
        t = m.group(0)
        if t == "\\begin{":
            depth += 1
        elif t == "\\end{":
            depth -= 1
        elif depth == 0:
            parts.append(body[last:m.start()])
            last = m.end()
        i = m.end()
    parts.append(body[last:])
    return parts


def parse_list(body, tag, ctx):
    items = []
    # \item[label] text  |  \item text
    parts = split_items(body)
    for p in parts[1:]:
        p = p.strip()
        term = None
        if p.startswith("["):
            d, m = 0, 0
            while m < len(p):
                if p[m] == "[":
                    d += 1
                elif p[m] == "]":
                    d -= 1
                    if d == 0:
                        break
                m += 1
            term = p[1:m]
            p = p[m+1:]
        items.append((inline(term, ctx) if term else None, blocks(p, ctx)))
    return items


def blocks(text, ctx):
    """Turn a run of LaTeX into a list of HTML strings."""
    out = []
    while True:
        # find the earliest environment of any kind
        best = None
        for name in list(ENV_BOX) + list(ENV_CODE) + \
                ["keytable", "yzntable", "itemize", "enumerate", "description",
                 "equation", "align", "figure", "quote", "equation*", "align*",
                 "center", "tabular"]:
            k = text.find("\\begin{%s}" % name)
            if k >= 0 and (best is None or k < best[0]):
                best = (k, name)
        if best is None:
            out += paragraphs(text, ctx)
            return out
        k, name = best
        before, opt, body, after = split_env(text, name)
        out += paragraphs(before, ctx)

        if name in ENV_BOX:
            cls, default = ENV_BOX[name]
            title = inline(opt, ctx) if opt else default
            out.append('<div class="box %s"><div class="box-t">%s</div>%s</div>'
                       % (cls, title, "".join(blocks(body, ctx))))
        elif name in ENV_CODE:
            out.append('<pre class="code %s"><code>%s</code></pre>'
                       % (ENV_CODE[name], html.escape(body.strip("\n"))))
        elif name in ("keytable", "yzntable"):
            cap, label, header, rows = parse_table(body, name, ctx)
            n = ctx["nums"].get(label, "")
            idattr = ' id="%s"' % label if label else ""
            th = "".join("<th>%s</th>" % c for c in header)
            tb = "".join("<tr>%s</tr>" % "".join(
                             '<td colspan="%d">%s</td>' % c if isinstance(c, tuple)
                             else "<td>%s</td>" % c for c in r)
                         for r in rows)
            out.append('<figure class="tbl"%s><figcaption><b>%s</b> %s'
                       '</figcaption><div class="scroll"><table>'
                       '<thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>'
                       '</figure>' % (idattr, WORDS["table"] + ("" if LANG == "ja" else " ") + n, cap, th, tb))
        elif name in ("itemize", "enumerate"):
            tag = "ul" if name == "itemize" else "ol"
            items = parse_list(body, tag, ctx)
            lis = "".join("<li>%s%s</li>" % (
                ("<b>%s</b> " % t) if t else "", "".join(b)) for t, b in items)
            out.append("<%s>%s</%s>" % (tag, lis, tag))
        elif name == "description":
            items = parse_list(body, "dl", ctx)
            dd = "".join("<dt>%s</dt><dd>%s</dd>" % (t or "", "".join(b))
                         for t, b in items)
            out.append("<dl>%s</dl>" % dd)
        elif name in ("equation", "align", "equation*", "align*"):
            labs = re.findall(r"\\label\{([^}]*)\}", body)
            lab = labs[0] if labs else None
            body = re.sub(r"\\label\{[^}]*\}", "", body)
            tex = body.strip()
            if name.startswith("align"):
                tex = "\\begin{aligned}%s\\end{aligned}" % tex
            slot = math_slot(math_text(tex), True)
            num = ctx["nums"].get(lab)
            anchors = "".join('<span id="%s"></span>' % L for L in labs)
            tagspan = ('<span class="eqno">(%s)</span>' % num) if num else ""
            out.append('<div class="eq">%s\x02M%d\x02%s</div>'
                       % (anchors, slot, tagspan))
        elif name == "figure":
            m = re.search(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", body)
            src = m.group(1) if m else ""
            capm = cmd_with_arg(body, "caption")
            cap, lab = "", None
            if capm:
                raw = capm[0][2]
                lm = re.search(r"\\label\{([^}]*)\}", raw)
                if lm:
                    lab = lm.group(1)
                    raw = raw[:lm.start()] + raw[lm.end():]
                cap = inline(raw, ctx)
            n = ctx["nums"].get(lab, "")
            idattr = ' id="%s"' % lab if lab else ""
            out.append('<figure class="fig"%s>%s<figcaption><b>%s</b> %s'
                       '</figcaption></figure>'
                       % (idattr, figure_svg(ctx, src),
                          WORDS["figure"] + ("" if LANG == "ja" else " ") + n, cap))
        elif name == "center":
            out += blocks(body, ctx)
        elif name == "tabular":
            b = body.strip()
            _, j2 = take_group(b, b.index("{"))
            rows = []
            for line in re.split(r"\\\\", b[j2:]):
                line = line.strip()
                if not line or line.startswith("\\"):
                    continue
                rows.append([inline(c.strip(), ctx) for c in split_cells(line)])
            tb = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % c for c in r)
                         for r in rows)
            out.append('<div class="scroll plainrows"><table><tbody>%s</tbody>'
                       '</table></div>' % tb)
        elif name == "quote":
            out.append('<blockquote>%s</blockquote>'
                       % "".join(blocks(body, ctx)))
        text = after


PARBREAK = re.compile(r"\n\s*\n")


def paragraphs(text, ctx):
    out = []
    # headings first
    pieces = re.split(r"(\\(?:section|subsection|subsubsection|paragraph)\*?\s*\{)",
                      text)
    buf = pieces[0]
    i = 1
    chunks = [("", buf)]
    while i < len(pieces):
        kind = re.search(r"\\(\w+)", pieces[i]).group(1)
        arg, j = take_group("{" + pieces[i+1], 0)
        chunks.append((kind, "{%s}%s" % (arg, pieces[i+1][j-1:])))
        i += 2
    for kind, chunk in chunks:
        if kind:
            starred = chunk.startswith("*")
            if starred:
                chunk = chunk[1:]
            arg, j = take_group(chunk, 0)
            rest = chunk[j:]
            lab = None
            lm = re.search(r"\\label\{([^}]*)\}", rest)
            if lm and lm.start() < 40:
                lab = lm.group(1)
                rest = rest[:lm.start()] + rest[lm.end():]
            if kind == "paragraph":
                out.append('<p class="runin"><b>%s</b> %s</p>'
                           % (inline(arg, ctx), ""))
                # the run-in head introduces the paragraph that follows it
                out[-1] = '<p class="runin"><b>%s</b> ' % inline(arg, ctx)
                body = PARBREAK.split(rest.strip(), 1)
                out[-1] += inline(body[0].strip(), ctx) + "</p>"
                rest = body[1] if len(body) > 1 else ""
            else:
                lvl = {"section": 2, "subsection": 3, "subsubsection": 4}[kind]
                c = ctx["counter"]
                if starred:
                    num = ""
                else:
                    # the same counter LaTeX keeps, so an unlabeled heading
                    # is numbered too
                    if lvl == 2:
                        c[0] += 1; c[1] = 0; c[2] = 0
                        num = "%s.%d" % (ctx["chnum"], c[0])
                    elif lvl == 3:
                        c[1] += 1; c[2] = 0
                        num = "%s.%d.%d" % (ctx["chnum"], c[0], c[1])
                    else:
                        c[2] += 1
                        num = "%s.%d.%d.%d" % (ctx["chnum"], c[0], c[1], c[2])
                if not lab:
                    ctx["auto"] += 1
                    lab = "%s-h%d" % (ctx["chlabel"], ctx["auto"])
                idattr = ' id="%s"' % lab
                ctx["outline"].append((lvl, num, inline(arg, ctx), lab))
                out.append('<span data-sec="%s"></span>'
                           '<h%d%s><span class="n">%s</span> %s</h%d>'
                           % (lab, lvl, idattr, num, inline(arg, ctx), lvl))
            chunk = rest
        for para in PARBREAK.split(chunk):
            para = para.strip()
            if not para:
                continue
            out.append("<p>%s</p>" % inline(para, ctx))
    return out


# ---------------------------------------------------------------------------
#  5. numbering: walk the sources once to assign chapter/section/table numbers
# ---------------------------------------------------------------------------

def strip_comments(text):
    """Remove LaTeX comments, honoring the rule that a trailing % also eats
    the following newline and leading whitespace. Code environments are left
    alone: a # or % inside a listing is data, not a comment."""
    holes = []

    def hide(m):
        holes.append(m.group(0))
        return "\x04%d\x04" % (len(holes) - 1)
    for name in ENV_CODE:
        text = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (name, name),
                      hide, text, flags=re.S)
    text = re.sub(r"(?<!\\)%.*?\n\s*", "", text)
    text = re.sub(r"\x04(\d+)\x04", lambda m: holes[int(m.group(1))], text)
    return text


def strip_verbatim(text):
    """Blank out code environments so their contents cannot be mistaken for
    markup during the numbering pass."""
    for name in ENV_CODE:
        text = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (name, name),
                      "", text, flags=re.S)
    return text


def number_everything(items):
    nums, titles = {}, {}
    chn, appn = 0, 0
    for it in items:
        src = open(os.path.join(CH, it["file"] + ".tex"), encoding="utf-8").read()
        src = strip_verbatim(src)
        starred = bool(re.search(r"\\chapter\*", src))
        m = re.search(r"\\chapter\*?\s*\{", src)
        if not m:
            it["skip"] = True
            continue
        title, j = take_group(src, m.end() - 1)
        it["title"] = title
        lm = re.search(r"\\label\{([^}]*)\}", src[j:j+200])
        it["label"] = lm.group(1) if lm else it["file"]
        if starred:
            it["num"] = ""
        elif it["appendix"]:
            it["num"] = chr(ord("A") + appn)
            appn += 1
        else:
            chn += 1
            it["num"] = str(chn)
        nums[it["label"]] = it["num"]
        titles[it["label"]] = title
        # sections
        s, ss, sss = 0, 0, 0
        tabn, fign, eqn = 0, 0, 0
        for mm in re.finditer(
                r"\\(section|subsection|subsubsection)\s*\{|"
                r"\\begin\{(keytable|yzntable|figure|equation|align)\}", src):
            if mm.group(1) == "section":
                s += 1; ss = 0; sss = 0
                cur = "%s.%d" % (it["num"], s)
            elif mm.group(1) == "subsection":
                ss += 1; sss = 0
                cur = "%s.%d.%d" % (it["num"], s, ss)
            elif mm.group(1) == "subsubsection":
                sss += 1
                cur = "%s.%d.%d.%d" % (it["num"], s, ss, sss)
            else:
                env = mm.group(2)
                endm = src.find("\\end{%s}" % env, mm.end())
                tail = src[mm.end():endm if endm > 0 else mm.end() + 4000]
                if env in ("keytable", "yzntable"):
                    tail = tail[:tail.find("\\\\") if "\\\\" in tail else len(tail)]
                lb = re.search(r"\\label\{([^}]*)\}", tail)
                if env in ("keytable", "yzntable"):
                    tabn += 1
                    cur = "%s.%d" % (it["num"], tabn)
                elif env == "figure":
                    fign += 1
                    cur = "%s.%d" % (it["num"], fign)
                else:
                    # an align may carry one label per line
                    labs = [x for x in re.findall(r"\\label\{([^}]*)\}", tail)]
                    if len(labs) > 1:
                        for L in labs:
                            eqn += 1
                            nums[L] = "%s.%d" % (it["num"], eqn)
                        continue
                    eqn += 1
                    cur = "%s.%d" % (it["num"], eqn)
                if lb:
                    nums[lb.group(1)] = cur
                continue
            tail = src[mm.end():mm.end() + 200]
            lb = re.search(r"\\label\{([^}]*)\}", tail)
            if lb:
                nums[lb.group(1)] = cur
                titles[lb.group(1)] = take_group("{" + src[mm.end():], 0)[0]
    return nums, titles


# ---------------------------------------------------------------------------
#  6. figures: SVG, inlined
# ---------------------------------------------------------------------------

def figure_svg(ctx, src):
    """The inlined SVG of one \\includegraphics path, or a clear failure."""
    if src not in ctx["figs"]:
        sys.exit("tohtml.py: no figure for \\includegraphics{%s}.\n"
                 "  Looked in: %s" % (src, ", ".join(d for _, d in FIGDIRS)))
    return ctx["figs"][src]


def load_figures():
    """Inline SVG for every figure the chapters can reference.

    A matching .svg beside the .pdf is preferred and is what the figure
    scripts write; the logo comes from ../logo/. Converting the PDF is the
    last resort, because pdftocairo turns text into glyph references that do
    not survive without the original font.
    """
    out = {}
    def clean(body, stem, cls="figsvg"):
        body = re.sub(r"<\?xml[^>]*\?>", "", body)
        body = re.sub(r"<!DOCTYPE[^>]*>", "", body, flags=re.S)
        body = re.sub(r"<!--.*?-->", "", body, flags=re.S)
        body = re.sub(r'(<svg[^>]*?)\s(?:width|height)="[^"]*"', r"\1", body)
        body = re.sub(r'(<svg[^>]*?)\s(?:width|height)="[^"]*"', r"\1", body)
        #  Every matplotlib figure names its groups and glyphs the same way
        #  (figure_1, DejaVuSans-30, ...). Inlined side by side, those ids
        #  would repeat, so prefix the ids of each figure, and the references
        #  to them, with the figure's own name.
        ids = set(re.findall(r'\bid="([^"]+)"', body))
        if ids:
            pre = re.sub(r"[^A-Za-z0-9_-]", "-", stem) + "--"
            body = re.sub(r'\bid="([^"]+)"', lambda m: 'id="%s%s"' % (pre, m.group(1)), body)
            body = re.sub(r'((?:xlink:)?href="#|url\(#)([^")]+)',
                          lambda m: m.group(1) + (pre + m.group(2) if m.group(2) in ids
                                                  else m.group(2)), body)
        return body.replace("<svg ", '<svg class="%s" ' % cls, 1).strip()

    logo = os.path.join(HERE, "..", "logo")
    for prefix, fdir in FIGDIRS:
        if not os.path.isdir(fdir):
            continue
        for f in sorted(os.listdir(fdir)):
            if not f.endswith(".pdf"):
                continue
            stem = f[:-4]
            svg = os.path.join(fdir, stem + ".svg")
            alt = os.path.join(logo, stem + ".svg")
            if os.path.exists(svg):
                body = open(svg, encoding="utf-8").read()
            elif os.path.exists(alt):
                body = open(alt, encoding="utf-8").read()
            else:
                if shutil.which("pdftocairo") is None:
                    sys.exit(
                        "tohtml.py: %s has no .svg beside it and pdftocairo was\n"
                        "  not found to convert it (Debian/Ubuntu: sudo apt install\n"
                        "  poppler-utils). Adding %s%s.svg also fixes this,\n"
                        "  and gives a better result." % (f, prefix, stem))
                tmp = os.path.join(BUILD, "figsvg", "%s.svg" % stem)
                os.makedirs(os.path.dirname(tmp), exist_ok=True)
                subprocess.run(["pdftocairo", "-svg", os.path.join(fdir, f), tmp],
                               check=True, capture_output=True)
                body = open(tmp, encoding="utf-8").read()
                print("  note: %s was converted from PDF; add an .svg for fidelity"
                      % f, file=sys.stderr)
            out[prefix + f] = clean(body, stem)
    return out


# ---------------------------------------------------------------------------
#  7. the page
# ---------------------------------------------------------------------------

CSS = r"""
:root{
  --ink:#22201D; --celadon:#5E8C78; --shadow:#48705E; --amber:#C08A3E;
  --lapis:#2E4A7D; --paper:#F4F1EA; --bg:#ffffff; --muted:#6c6862;
  --rule:#e3e0d8; --codebg:#f7f5f0; --sidebar:#fbfaf7;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --ink:#e9e6df; --celadon:#8fc0a9; --shadow:#a9cdb9; --amber:#e0ae66;
    --lapis:#9db9e8; --paper:#26241f; --bg:#1b1a17; --muted:#a9a49b;
    --rule:#37342e; --codebg:#232119; --sidebar:#201e1a;
  }
}
:root[data-theme="dark"]{
  --ink:#e9e6df; --celadon:#8fc0a9; --shadow:#a9cdb9; --amber:#e0ae66;
  --lapis:#9db9e8; --paper:#26241f; --bg:#1b1a17; --muted:#a9a49b;
  --rule:#37342e; --codebg:#232119; --sidebar:#201e1a;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth; scroll-padding-top:1.5rem}
body{margin:0;background:var(--bg);color:var(--ink);
 font:16px/1.65 "Palatino Linotype",Palatino,"Book Antiqua","Noto Serif",Georgia,serif;}
code,pre,.k,.v,.f,.c,.e{font-family:ui-monospace,"SF Mono",Menlo,Consolas,"DejaVu Sans Mono",monospace}

/* layout */
#wrap{display:grid;grid-template-columns:19rem minmax(0,1fr);min-height:100vh}
#side{background:var(--sidebar);border-right:1px solid var(--rule);
 position:sticky;top:0;height:100vh;overflow:auto;padding:1.1rem 0 3rem}
main{padding:2.2rem 2.6rem 8rem;max-width:56rem}
@media (max-width:950px){
  #wrap{grid-template-columns:minmax(0,1fr)}
  #side{position:fixed;inset:0 auto 0 0;width:19rem;z-index:40;
   transform:translateX(-100%);transition:transform .2s ease}
  #side.open{transform:none}
  main{padding:1.2rem 1.1rem 6rem}
  #burger{display:flex !important}
}

/* sidebar */
.brand{display:flex;align-items:center;gap:.6rem;padding:0 1.1rem 1rem;
 border-bottom:1px solid var(--rule);margin-bottom:.8rem}
.brand svg{width:30px;height:30px;flex:none}
.brand a.lang{display:inline-block;margin-top:.25rem;font:400 .72rem/1.3 ui-sans-serif,system-ui,sans-serif;color:var(--muted)}
.brand b{font:600 1.15rem/1 ui-sans-serif,system-ui,sans-serif;letter-spacing:.01em}
.brand span{display:block;font:400 .72rem/1.3 ui-sans-serif,system-ui,sans-serif;
 color:var(--muted);margin-top:.18rem}
#side nav{padding:0 .6rem}
#side .part{font:600 .68rem/1 ui-sans-serif,system-ui,sans-serif;
 text-transform:uppercase;letter-spacing:.09em;color:var(--muted);
 margin:1.1rem .5rem .35rem}
#side a{display:block;padding:.2rem .5rem;border-radius:4px;color:var(--ink);
 text-decoration:none;font:400 .85rem/1.45 ui-sans-serif,system-ui,sans-serif}
#side a:hover{background:var(--paper)}
#side a.ch{font-weight:600}
#side a.s2{padding-left:1.6rem;font-size:.8rem;color:var(--muted);display:none}
#side a.s2.show{display:block}
#side a.cur{background:var(--celadon);color:#fff}
#side a .n{color:var(--muted);margin-right:.35rem;font-variant-numeric:tabular-nums}
#side a.cur .n{color:#fff;opacity:.8}

/* search */
.sbox{padding:0 1.1rem .7rem;position:sticky;top:0;background:var(--sidebar);z-index:5}
#q{width:100%;padding:.48rem .6rem;border:1px solid var(--rule);border-radius:6px;
 background:var(--bg);color:var(--ink);font:400 .86rem/1.3 ui-sans-serif,system-ui,sans-serif}
#q:focus{outline:2px solid var(--celadon);outline-offset:-1px;border-color:transparent}
.hint{font:400 .68rem/1.3 ui-sans-serif,system-ui,sans-serif;color:var(--muted);
 margin:.35rem .1rem 0}
#results{padding:.2rem .6rem}
#results .r{display:block;padding:.45rem .5rem;border-radius:5px;text-decoration:none;
 color:var(--ink);border-bottom:1px solid var(--rule)}
#results .r:hover{background:var(--paper)}
#results .rt{font:600 .82rem/1.3 ui-sans-serif,system-ui,sans-serif}
#results .rc{font:400 .7rem/1.25 ui-sans-serif,system-ui,sans-serif;color:var(--muted)}
#results .rx{font:400 .76rem/1.4 ui-sans-serif,system-ui,sans-serif;color:var(--muted);
 margin-top:.18rem}
#results mark{background:var(--amber);color:#1b1a17;border-radius:2px;padding:0 .1em}
#results .none{padding:.6rem .5rem;color:var(--muted);
 font:400 .8rem ui-sans-serif,system-ui,sans-serif}

/* content */
h1,h2,h3,h4{font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;
 font-weight:650;line-height:1.25;color:var(--shadow);scroll-margin-top:1rem}
h1{font-size:2rem;margin:0 0 .3rem;padding-bottom:.5rem;border-bottom:2px solid var(--rule)}
h1 .n{color:var(--celadon);margin-right:.5rem}
h2{font-size:1.32rem;margin:2.6rem 0 .7rem}
h3{font-size:1.06rem;margin:1.8rem 0 .5rem;color:var(--ink)}
h4{font-size:.97rem;margin:1.3rem 0 .4rem;color:var(--ink);font-style:italic;font-weight:500}
h2 .n,h3 .n,h4 .n{color:var(--muted);margin-right:.45rem;font-weight:500;
 font-variant-numeric:tabular-nums}
.chapter{display:none}
.chapter.on{display:block}
p{margin:0 0 1rem}
p.runin b{color:var(--shadow)}
strong{font-weight:650}
ul,ol{margin:0 0 1rem;padding-left:1.5rem}
li{margin:.3rem 0}
dl{margin:0 0 1rem}
dt{font-weight:650;color:var(--shadow);margin-top:.7rem}
dd{margin:.15rem 0 0 1.3rem}
blockquote{margin:1.1rem 0;padding:.1rem 0 .1rem 1.1rem;
 border-left:3px solid var(--celadon);color:var(--ink)}
a.xref{color:var(--lapis);text-decoration:none;border-bottom:1px solid transparent}
a.xref:hover{border-bottom-color:currentColor}
a.ext{color:var(--shadow);word-break:break-all}

/* code voice */
code{font-size:.87em;background:var(--codebg);padding:.1em .32em;border-radius:4px;
 border:1px solid var(--rule)}
code.k{color:var(--lapis)} code.f{color:var(--shadow)}
code.req{background:none;border:none;color:var(--amber);
 font:600 .74rem ui-sans-serif,system-ui,sans-serif}
pre.code{background:var(--codebg);border:1px solid var(--rule);
 border-left:3px solid var(--celadon);border-radius:0 6px 6px 0;
 padding:.8rem 1rem;overflow:auto;margin:1.1rem 0}
pre.code code{background:none;border:none;padding:0;font-size:.82rem;line-height:1.55}

/* boxes */
.box{margin:1.3rem 0;padding:.85rem 1rem .1rem;border-radius:0 6px 6px 0;
 border-left:3px solid;position:relative}
.box-t{font:650 .78rem/1.3 ui-sans-serif,system-ui,sans-serif;
 margin-bottom:.5rem;letter-spacing:.01em}
.box.warn{border-color:var(--amber);background:color-mix(in srgb,var(--amber) 8%,transparent)}
.box.warn .box-t{color:var(--amber)}
.box.note{border-color:var(--celadon);background:color-mix(in srgb,var(--celadon) 8%,transparent)}
.box.note .box-t{color:var(--shadow)}
.box.tip{border-color:var(--lapis);background:color-mix(in srgb,var(--lapis) 8%,transparent)}
.box.tip .box-t{color:var(--lapis)}
.box p:last-child{margin-bottom:.85rem}

/* tables and figures */
figure{margin:1.6rem 0}
figcaption{font:400 .82rem/1.5 ui-sans-serif,system-ui,sans-serif;color:var(--muted);
 margin-bottom:.5rem}
figcaption b{color:var(--shadow)}
.fig figcaption{margin:.6rem 0 0}
.scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:6px}
table{border-collapse:collapse;width:100%;font-size:.84rem;line-height:1.45}
th{background:var(--paper);text-align:left;font:650 .78rem/1.4 ui-sans-serif,system-ui,sans-serif;
 color:var(--shadow);padding:.5rem .7rem;border-bottom:1px solid var(--rule);
 white-space:nowrap}
td{padding:.45rem .7rem;border-bottom:1px solid var(--rule);vertical-align:top}
tr:last-child td{border-bottom:none}
td code{font-size:.85em}
.figsvg{width:100%;height:auto;display:block;background:#fff;border-radius:6px;
 border:1px solid var(--rule)}

/* maths */
.eq{margin:1.3rem 0;text-align:center;position:relative;overflow-x:auto}
.eqno{position:absolute;right:0;top:50%;transform:translateY(-50%);
 color:var(--muted);font-size:.85rem}
mjx-container[display]{margin:0 !important}
svg[data-mjx]{vertical-align:-0.2ex}
.mjx-c,mjx-container svg{fill:currentColor;stroke:currentColor}

/* footnotes */
sup.fn{color:var(--amber);font-size:.72em;padding-left:.1em}
.fnlist{margin-top:2.5rem;padding-top:.8rem;border-top:1px solid var(--rule);
 font-size:.83rem;color:var(--muted)}
.fnlist li{margin:.35rem 0}

/* chrome */
#burger{display:none;position:fixed;left:.7rem;top:.7rem;z-index:50;
 width:2.4rem;height:2.4rem;border-radius:7px;border:1px solid var(--rule);
 background:var(--bg);color:var(--ink);align-items:center;justify-content:center;
 cursor:pointer;font-size:1.1rem}
#theme{position:fixed;right:1rem;top:.8rem;z-index:30;border:1px solid var(--rule);
 background:var(--bg);color:var(--muted);border-radius:7px;cursor:pointer;
 padding:.3rem .6rem;font:400 .74rem ui-sans-serif,system-ui,sans-serif}
.pager{display:flex;justify-content:space-between;gap:1rem;margin-top:3.5rem;
 padding-top:1.2rem;border-top:1px solid var(--rule)}
.pager a{color:var(--lapis);text-decoration:none;
 font:400 .85rem ui-sans-serif,system-ui,sans-serif;max-width:45%}
.pager a.next{text-align:right;margin-left:auto}
.pager a b{display:block;color:var(--muted);font-size:.7rem;text-transform:uppercase;
 letter-spacing:.08em;font-weight:600}
"""

JS = r"""
(function(){
  var side = document.getElementById('side');
  var q    = document.getElementById('q');
  var res  = document.getElementById('results');
  var nav  = document.getElementById('nav');
  var chapters = Array.prototype.slice.call(document.querySelectorAll('.chapter'));

  // ---- show one chapter at a time, driven by the hash ---------------------
  function idOf(el){ return el.getAttribute('data-ch'); }
  function showChapter(ch, anchor){
    chapters.forEach(function(c){ c.classList.toggle('on', idOf(c) === ch); });
    var curLink = null;
    document.querySelectorAll('#nav a').forEach(function(a){
      var mine = a.getAttribute('data-ch') === ch;
      if(a.classList.contains('s2')){
        a.classList.toggle('show', mine);
      } else {
        a.classList.toggle('cur', mine);
        if(mine) curLink = a;
      }
    });
    if(curLink && nav.style.display !== 'none'){
      var r = curLink.getBoundingClientRect();
      if(r.top < 120 || r.bottom > window.innerHeight - 40){
        curLink.scrollIntoView({block:'center'});
      }
    }
    if(anchor){
      var t = document.getElementById(anchor);
      if(t){ t.scrollIntoView(); return; }
    }
    window.scrollTo(0,0);
  }
  function chapterOfAnchor(id){
    var el = document.getElementById(id);
    while(el && !el.classList.contains('chapter')) el = el.parentElement;
    return el ? idOf(el) : null;
  }
  function route(){
    var h = decodeURIComponent(location.hash.replace(/^#/,''));
    if(!h){ showChapter(idOf(chapters[0])); return; }
    var ch = chapterOfAnchor(h);
    if(ch) showChapter(ch, h); else showChapter(idOf(chapters[0]));
    side.classList.remove('open');
  }
  window.addEventListener('hashchange', route);

  // ---- search -------------------------------------------------------------
  var IDX = window.__IDX__;
  function esc(s){ return s.replace(/[&<>]/g, function(c){
    return {'&':'&amp;','<':'&lt;','>':'&gt;'}[c]; }); }

  IDX.forEach(function(e){ e.lb = e.b.toLowerCase(); });

  function score(entry, terms){
    var s = 0, ok = 0;
    for(var i=0;i<terms.length;i++){
      var t = terms[i];
      var inT = entry.t.toLowerCase().indexOf(t);
      var inB = entry.lb.indexOf(t);
      if(inT < 0 && inB < 0) return 0;
      ok++;
      if(inT >= 0) s += (inT === 0 ? 40 : 25);
      if(inB >= 0) s += 6;
      // a whole-word hit in the body counts for more than a substring
      if(new RegExp('\\b' + t.replace(/[.*+?^${}()|[\]\\]/g,'\\$&')).test(entry.lb)) s += 4;
    }
    return ok === terms.length ? s : 0;
  }

  function snippet(body, terms){
    var lo = body.toLowerCase(), at = -1;
    for(var i=0;i<terms.length && at<0;i++) at = lo.indexOf(terms[i]);
    if(at < 0) at = 0;
    var from = Math.max(0, at - 60), to = Math.min(body.length, at + 140);
    var txt = (from > 0 ? '…' : '') + body.slice(from, to) +
              (to < body.length ? '…' : '');
    txt = esc(txt);
    terms.forEach(function(t){
      txt = txt.replace(new RegExp('(' + t.replace(/[.*+?^${}()|[\]\\]/g,'\\$&') + ')', 'ig'),
                        '<mark>$1</mark>');
    });
    return txt;
  }

  function search(){
    var v = q.value.trim().toLowerCase();
    if(v.length < 2){ res.innerHTML = ''; nav.style.display = ''; return; }
    nav.style.display = 'none';
    var terms = v.split(/\s+/);
    var hits = [];
    for(var i=0;i<IDX.length;i++){
      var s = score(IDX[i], terms);
      if(s) hits.push([s, IDX[i]]);
    }
    hits.sort(function(a,b){ return b[0]-a[0]; });
    if(!hits.length){ res.innerHTML = '<div class="none">Nothing found.</div>'; return; }
    res.innerHTML = hits.slice(0,60).map(function(h){
      var e = h[1];
      return '<a class="r" href="#' + e.i + '">' +
             '<div class="rt">' + (e.n ? '<span style="opacity:.6">' + e.n + '</span> ' : '') +
             esc(e.t) + '</div>' +
             '<div class="rc">' + esc(e.c) + '</div>' +
             '<div class="rx">' + snippet(e.b, terms) + '</div></a>';
    }).join('');
  }
  q.addEventListener('input', search);
  q.addEventListener('keydown', function(e){
    if(e.key === 'Escape'){ q.value=''; search(); q.blur(); }
    if(e.key === 'Enter'){ var a = res.querySelector('a.r'); if(a) a.click(); }
  });
  document.addEventListener('keydown', function(e){
    if((e.key === '/' || (e.key === 'k' && (e.metaKey || e.ctrlKey))) &&
       document.activeElement !== q){
      e.preventDefault(); side.classList.add('open'); q.focus(); q.select();
    }
  });
  res.addEventListener('click', function(e){
    var a = e.target.closest('a.r'); if(a){ side.classList.remove('open'); }
  });

  // ---- chrome -------------------------------------------------------------
  document.getElementById('burger').addEventListener('click', function(){
    side.classList.toggle('open');
  });
  var tb = document.getElementById('theme');
  function applyTheme(t){
    if(t) document.documentElement.setAttribute('data-theme', t);
    else document.documentElement.removeAttribute('data-theme');
    tb.textContent = t === 'dark' ? 'Light' : (t === 'light' ? 'Dark' : 'Dark');
  }
  var saved = null;
  try { saved = localStorage.getItem('yuzen-theme'); } catch(err){}
  applyTheme(saved);
  tb.addEventListener('click', function(){
    var cur = document.documentElement.getAttribute('data-theme');
    var next = cur === 'dark' ? 'light' : 'dark';
    applyTheme(next);
    try { localStorage.setItem('yuzen-theme', next); } catch(err){}
  });

  route();
})();
"""


def clean_index(t):
    """Index text: drop the maths sentinels and tidy the spacing.

    Case is preserved, because the snippet a reader sees is cut from this
    string and key names such as NewtonEpsilonRel are case-significant; matching
    folds case separately, at query time.
    """
    t = re.sub("\x02M?\\d*\x02?", " ", t)
    t = "".join(c for c in t if c >= " ")
    t = re.sub(r"\s+", " ", t)
    # the tag stripper leaves a space before punctuation; close it up
    t = re.sub(r"\s+([,.;:)\]])", r"\1", t)
    t = re.sub(r"([(\[])\s+", r"\1", t)
    return t.strip()


def text_of(h):
    """Plain text of an HTML fragment, for the search index."""
    h = re.sub(r"<svg.*?</svg>", " ", h, flags=re.S)
    h = re.sub(r"<[^>]+>", " ", h)
    h = html.unescape(h)
    return re.sub(r"\s+", " ", h).strip()


def main():
    if "--check-install" in sys.argv[1:]:
        sys.exit(report_install())

    items = read_structure()
    nums, titles = number_everything(items)
    figs = load_figures()
    ctx = {"nums": nums, "titles": titles, "figs": figs,
           "missing": set(), "unhandled": set(), "outline": [], "footnotes": []}

    body_parts, nav_parts, index = [], [], []
    last_part = None
    order = [it for it in items if not it.get("skip")]

    for n, it in enumerate(order):
        src = open(os.path.join(CH, it["file"] + ".tex"), encoding="utf-8").read()
        src = strip_comments(src)
        m = re.search(r"\\chapter\*?\s*\{", src)
        _, j = take_group(src, m.end() - 1)
        rest = src[j:]
        rest = re.sub(r"\\label\{[^}]*\}", "", rest, count=1)
        rest = re.sub(r"\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", rest)
        rest = re.sub(r"\\markboth\{[^}]*\}\{[^}]*\}", "", rest)

        ctx["outline"] = []
        ctx["footnotes"] = []
        ctx["counter"] = [0, 0, 0]
        ctx["auto"] = 0
        ctx["chnum"] = it["num"]
        ctx["chlabel"] = it["label"]
        html_body = "".join(blocks(rest, ctx))
        if ctx["footnotes"]:
            html_body += '<ol class="fnlist">%s</ol>' % "".join(
                "<li>%s</li>" % inline(f, ctx) for f in ctx["footnotes"])

        anchor = it["label"]
        kind = WORDS["appendix"] if it["appendix"] else WORDS["chapter"]
        head = '<h1><span class="n">%s</span> %s</h1>' % (it["num"], html.escape(it["title"])) \
            if it["num"] else '<h1>%s</h1>' % html.escape(it["title"])

        prev_a = order[n-1] if n else None
        next_a = order[n+1] if n+1 < len(order) else None
        pager = '<div class="pager">'
        if prev_a:
            pager += '<a href="#%s"><b>%s</b>%s</a>' % (
                prev_a["label"], WORDS["prev"], html.escape(prev_a["title"]))
        if next_a:
            pager += '<a class="next" href="#%s"><b>%s</b>%s</a>' % (
                next_a["label"], WORDS["next"], html.escape(next_a["title"]))
        pager += "</div>"

        body_parts.append('<section class="chapter" data-ch="%s" id="%s">%s%s%s</section>'
                          % (anchor, anchor, head, html_body, pager))

        # sidebar
        if it["part"] != last_part:
            nav_parts.append('<div class="part">%s</div>' % html.escape(it["part"] or ""))
            last_part = it["part"]
        nav_parts.append('<a class="ch" href="#%s" data-ch="%s">'
                         '<span class="n">%s</span>%s</a>'
                         % (anchor, anchor, it["num"], html.escape(it["title"])))
        for lvl, num, title, lab in ctx["outline"]:
            if lvl != 2 or not lab:
                continue
            nav_parts.append('<a class="s2" href="#%s" data-ch="%s">'
                             '<span class="n">%s</span>%s</a>'
                             % (lab, anchor, num, title))

        # Search index: one entry per chapter and per section, built from the
        # text of the rendered body with the maths stripped, so that a section
        # is indexed in full rather than for as many characters of raw HTML as
        # its first equation happens to leave room for.
        marked = re.sub(r'<span data-sec="([^"]*)"></span>',
                        lambda m: "\x05%s\x05" % m.group(1), html_body)
        plain = clean_index(text_of(marked))
        segs, pos = {}, plain.split("\x05")
        # segs: label -> text from that heading up to the next one
        i = 1
        while i + 1 < len(pos) + 1 and i < len(pos):
            lbl = pos[i].strip()
            body_txt = pos[i+1] if i + 1 < len(pos) else ""
            segs[lbl] = body_txt.strip()
            i += 2
        if LANG == "ja":
            chap_label = (kind % it["num"]) if it["num"] else ""
        else:
            chap_label = "%s %s" % (kind, it["num"]) if it["num"] else kind
        index.append({"i": anchor, "t": it["title"], "n": it["num"],
                      "c": chap_label, "b": clean_index(plain.replace("\x05", " "))})
        for lvl, num, title, lab in ctx["outline"]:
            if not lab:
                continue
            index.append({"i": lab, "t": text_of(title), "n": num,
                          "c": "%s \u00b7 %s" % (chap_label, it["title"]),
                          "b": segs.get(lab, "")})

    if ctx["missing"]:
        print("  undefined references:", ", ".join(sorted(ctx["missing"])[:20]),
              file=sys.stderr)
    ctx["unhandled"].update(UNKNOWN_UNITS)
    if ctx["unhandled"]:
        print("  unhandled LaTeX commands:",
              ", ".join(sorted(ctx["unhandled"])), file=sys.stderr)
        sys.exit("refusing to write an HTML manual with unconverted markup")

    svgs = render_all_math()
    page = "".join(body_parts)
    nav_html = "".join(nav_parts)
    #  A formula that occurs more than once is rendered once, but each copy of
    #  its SVG carries the glyph definitions under the same ids. Give every
    #  copy after the first its own id prefix (MJX-12- -> MJX-12-1-), so that
    #  ids stay unique in the page.
    uses = {}

    def put(m):
        k = int(m.group(1))
        n = uses[k] = uses.get(k, -1) + 1
        if n == 0:
            return svgs[k]
        return re.sub(r"MJX-(\d+)-", lambda g: "MJX-%s-%d-" % (g.group(1), n), svgs[k])
    sub = lambda t: re.sub(r"\x02M(\d+)\x02", put, t)
    page, nav_html = sub(page), sub(nav_html)
    for e in index:
        e["t"] = re.sub(r"\x02M\d+\x02", "", e["t"]).strip()

    mark = open(os.path.join(HERE, "..", "logo", "yuzen-mark.svg"),
                encoding="utf-8").read()
    mark = re.sub(r"<\?xml[^>]*\?>", "", mark).strip()
    mark = re.sub(r'(<svg[^>]*?)\s(width|height)="[^"]*"', r"\1", mark)

    W = WORDS
    js = (JS.replace("Nothing found.", html.escape(W["none"]))
            .replace("t === 'dark' ? 'Light' : (t === 'light' ? 'Dark' : 'Dark')",
                     "t === 'dark' ? %s : %s" % (json.dumps(W["light"]), json.dumps(W["dark"]))))
    out = """<!doctype html>
<html lang="%s"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<meta name="description" content="%s">
<style>%s</style>
</head><body>
<button id="burger" aria-label="Menu">&#9776;</button>
<button id="theme">%s</button>
<div id="wrap">
<aside id="side">
  <div class="brand">%s<div><b>Yuzen</b><span>%s</span><a class="lang" href="%s">%s</a></div></div>
  <div class="sbox">
    <input id="q" type="search" placeholder="%s"
           autocomplete="off" spellcheck="false">
    <div class="hint">%s</div>
  </div>
  <div id="results"></div>
  <nav id="nav">%s</nav>
</aside>
<main>%s</main>
</div>
<script>window.__IDX__=%s;</script>
<script>%s</script>
</body></html>
""" % (W["lang"], W["title"], W["desc"], CSS, W["dark"], mark, W["sub"], W["other_href"], W["other"],
       W["search"], W["hint"], nav_html, page, json.dumps(index, ensure_ascii=False), js)

    dest = os.path.join(HERE, OUT_HTML)
    open(dest, "w", encoding="utf-8").write(out)
    print("wrote %s  (%.1f MB, %d chapters, %d search entries, %d maths)"
          % (dest, len(out.encode()) / 1e6, len(order), len(index), len(svgs)))


if __name__ == "__main__":
    main()
