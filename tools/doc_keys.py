#!/usr/bin/env python3
"""Check that the documentation of the input keys says what the parser does.

    python3 tools/doc_keys.py                      # source and documents only
    python3 tools/doc_keys.py --exe build/bin/yuzen --workdir /tmp/dk
    python3 tools/doc_keys.py --list               # print the keys the parser reads

Every key Yuzen reads is read in yuzen/src/IO/Input.cpp, with one getin() per
key: with a default or without one (required), followed where there is one
by its range check or its list of accepted values. That file is read here,
with the member types of Input.h, and becomes a table of every block path
(required, optional or repeatable) and every key (type, default, required,
accepted values, range). It is then held against

  * the key tables of the reference chapters (keytable, and the yzntables
    with Key/Type/Default columns): every key read is in a table and every
    key in a table is read, with the same type, default, required mark, the
    values it accepts and the range it states;
  * the counts and the required-key and top-level-block statements of
    ch-inputformat.tex;
  * every \\key{...} of the manual and every dotted key path of README.md and
    docs/*.md, which must name a block or a key that exists;
  * every TOML listing of the manual and every file in cases/, whose keys
    must all be read by the parser;
  * the table of keys no shipped case uses (ch-ref-numerics.tex).

With --exe the binary is asked as well: a file holding every block and key
with its default must be read with no unknown entry and no range error, and
every case and every complete TOML listing of the manual must pass
`yuzen --check-input`, which reads and checks a file and runs nothing.

Only the Python standard library is used (tomllib, Python 3.11 or later; the
TOML checks are skipped on an older Python).

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys

try:
    import tomllib
except ImportError:          # Python < 3.11
    tomllib = None

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "yuzen", "src")
MANUAL = os.path.join(ROOT, "docs", "manual", "chapters")


# ---------------------------------------------------------------------------
#  Small helpers
# ---------------------------------------------------------------------------
def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def split_top(s, sep=","):
    """Split s at sep where no bracket or quote is open."""
    out, depth, cur, q = [], 0, "", None
    for c in s:
        if q:
            cur += c
            if c == q:
                q = None
            continue
        if c in "\"'":
            q = c
        elif c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        if c == sep and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += c
    if cur.strip():
        out.append(cur.strip())
    return out


def strip_comments_cpp(text):
    text = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def norm_target(t):
    return re.sub(r"\[[^\]]*\]", "", t.strip())


# ---------------------------------------------------------------------------
#  The code: Input.h (types) and Input.cpp (blocks, keys, defaults, checks)
# ---------------------------------------------------------------------------
class ClassNode:
    def __init__(self, name, parent=None):
        self.name, self.parent = name, parent
        self.members = {}      # member name -> type text
        self.classes = {}      # nested class name -> ClassNode

    def find_class(self, name):
        node = self
        while node is not None:
            if name in node.classes:
                return node.classes[name]
            node = node.parent
        return None


def parse_input_h(path):
    text = strip_comments_cpp(read(path))
    root = None
    stack = []
    pending = None
    for line in text.splitlines():
        s = line.strip()
        m = re.match(r"class\s+(\w+)\s*(?::.*)?$", s)
        if m:
            pending = m.group(1)
            continue
        if s.startswith("{") and pending:
            node = ClassNode(pending, stack[-1] if stack else None)
            if stack:
                stack[-1].classes[pending] = node
            else:
                root = node
            stack.append(node)
            pending = None
            continue
        if s.startswith("};") and stack:
            stack.pop()
            continue
        if not stack or "(" in s:
            continue
        m = re.match(r"((?:std::)?(?:vector|list)\s*<\s*[\w:]+\s*>|long\s+long|[\w:]+)\s+(\w+)\s*(?:=[^;]*)?;", s)
        if m:
            stack[-1].members[m.group(2)] = re.sub(r"\s+", "", m.group(1))
    return root


def member_type(root, dotted):
    node = root
    parts = dotted.split(".")
    for i, p in enumerate(parts):
        t = node.members.get(p)
        if t is None:
            return None
        if i == len(parts) - 1:
            return t
        inner = re.sub(r"^(?:std::)?vector<(.*)>$", r"\1", t)
        nxt = node.find_class(inner)
        if nxt is None:
            return None
        node = nxt
    return None


def doc_type(cpp_type):
    t = cpp_type.replace("std::", "")
    return {"string": "string", "double": "double", "int": "int", "longlong": "int", "bool": "bool",
            "vector<double>": "array", "vector<int>": "array",
            "vector<string>": "array"}.get(t, t)


def literal(text):
    """A C++ literal as a Python value, or None when it is not one."""
    t = text.strip()
    m = re.fullmatch(r'(?:std::string\s*\(\s*)?"([^"]*)"\s*\)?', t)
    if m:
        return m.group(1)
    if t in ("true", "false"):
        return t == "true"
    try:
        return float(t)
    except ValueError:
        return None


def version_series():
    text = read(os.path.join(SRC, "Version.h"))
    ma = re.search(r"#define\s+YUZEN_VERSION_Major\s+(\d+)", text).group(1)
    mi = re.search(r"#define\s+YUZEN_VERSION_Minor\s+(\d+)", text).group(1)
    return "%s.%s" % (ma, mi)


def scheme_keys():
    text = read(os.path.join(SRC, "numerics", "FiniteDifferences.cpp"))
    body = re.search(r"kKeys\[\]\s*=\s*\{(.*?)\};", text, re.S).group(1)
    return re.findall(r'\{\s*"(\w+)"', body)


class Key:
    def __init__(self, block, name, ctype, default, required, line):
        self.block, self.name, self.ctype = block, name, ctype
        self.default, self.required, self.line = default, required, line
        self.absent = None     # value when the whole (optional) block is absent
        self.enum = None
        self.range = None      # (lo, hi, lo_open, hi_open)

    @property
    def path(self):
        return self.block + "." + self.name if self.block else self.name

    @property
    def type(self):
        return doc_type(self.ctype)


def parse_input_cpp(path, header):
    raw = read(path)
    text = strip_comments_cpp(raw)
    lines = text.splitlines()
    consts = dict(re.findall(r"const\s+(?:int|long long|double)\s+(\w+)\s*=\s*([-\d.eE+]+)L*L*\s*;", text))

    blocks = {"": "required"}  # path -> required | optional | repeatable
    keys = {}                  # path -> Key
    by_target = {}             # normalized C++ target -> Key
    var_path = {}
    cand = []
    pushes = {}                # target -> [values] pushed before its getin
    assigns = {}               # target -> value from .assign(n, v)

    chain_re = r'((?:root_block|\w+)(?:->getBlock\(\s*"\w+"(?:\s*,\s*\w+)?\s*\))*)'

    def chain_path(chain):
        head = re.match(r"(\w+)", chain).group(1)
        base = "" if head == "root_block" else var_path.get(head, "")
        names = re.findall(r'getBlock\(\s*"(\w+)"', chain)
        return ".".join([p for p in [base] + names if p])

    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        if re.match(r"^(?:void|int|bool)\s+Input::\w+\s*\(", s):
            var_path, cand, pushes = {"sb": ""}, [], {}

        # --- blocks asked for ---
        for m in re.finditer(chain_re + r'->getBlockSize\(\s*"(\w+)"\s*\)', s):
            parent = chain_path(m.group(1))
            bpath = (parent + "." if parent else "") + m.group(2)
            kind = "optional"
            for la in lines[i:i + 14]:
                if re.search(r"bnum\s*!=\s*1|\w+\s*!=\s*1\)", la):
                    kind = "required"
                    break
                if ".resize(" in la:
                    kind = "repeatable"
                    break
                if re.search(r">\s*1\)", la):
                    kind = "optional"
                    break
            blocks[bpath] = kind

        # --- which block a variable points at ---
        m = re.search(r"(?:Block\s*\*\s*)?\b(\w+)\s*=\s*(?:\([^?;]*\)\s*\?\s*)?" + chain_re + r"\s*(?:;|:)", s)
        if m and "getBlock(" in m.group(2):
            var_path[m.group(1)] = chain_path(m.group(2))
        m = re.search(r"\b(\w+)\s*=\s*root_block\.get\(\)", s)
        if m:
            var_path[m.group(1)] = ""

        # --- default values of vectors, set up before their getin ---
        m = re.match(r"([\w.\[\]]+)\.push_back\(\s*([^)]*)\)\s*;", s)
        if m and m.group(1) != "cand":
            pushes.setdefault(norm_target(m.group(1)), []).append(float(m.group(2)))
        m = re.match(r"([\w.\[\]]+)\.assign\(\s*(\d+)\s*,\s*([^)]*)\)\s*;", s)
        if m:
            assigns[norm_target(m.group(1))] = [float(m.group(3))] * int(m.group(2))

        # --- the keys ---
        m = re.match(r"(\w+)->getin\(\s*\"(\w+)\"\s*,(.*)\)\s*;$", s)
        if m:
            var, name, args = m.group(1), m.group(2), split_top(m.group(3))
            target = norm_target(args[0])
            block = var_path.get(var, "")
            ctype = member_type(header, target)
            if ctype is None:   # a local variable
                d = re.search(r"([\w:]+)\s+" + re.escape(target) + r"\s*;", text)
                ctype = d.group(1) if d else "?"
            required = len(args) == 1
            default = None
            if not required:
                expr = args[1]
                if norm_target(expr) == target:
                    default = pushes.get(target) or assigns.get(target)
                else:
                    default = literal(expr)
                    if default is None and expr in consts:
                        default = float(consts[expr])
                if doc_type(ctype) == "bool" and isinstance(default, float):
                    default = bool(default)
            k = Key(block, name, ctype, default, required, i + 1)
            keys[k.path] = k
            by_target[target] = k
            continue

        # --- an effective default put in after the getin (if(!given) x = v) ---
        m = re.match(r"([\w.]+)\s*=\s*(\"[^\"]*\"|[-\w.+]+)\s*;$", s)
        if m and m.group(1) in by_target:
            prev = next((lines[j].strip() for j in range(i - 1, -1, -1) if lines[j].strip()), "")
            val = literal(m.group(2))
            k = by_target[m.group(1)]
            if re.match(r"if\s*\(\s*!", prev) and val is not None:
                k.default = val
            elif re.match(r"if\s*\(\s*\w+\s*==\s*0\s*\)", prev) and val is not None:
                k.absent = val

        # --- accepted values ---
        if re.search(r"\bcand\.clear\(\)", s):
            cand = []
        m = re.search(r'cand\.push_back\(\s*"([^"]*)"\s*(\+\s*C_Yuzen::VersionSeries\(\))?\s*\)', s)
        if m:
            cand.append(m.group(1) + (version_series() if m.group(2) else ""))
        m = re.search(r"CheckCandidates\(\s*([\w.\[\]]+)\s*,", s)
        if m and norm_target(m.group(1)) in by_target:
            by_target[norm_target(m.group(1))].enum = list(cand)
        m = re.search(r"C_Interpolation_base4::FromKey\(\s*([\w.]+)\s*,", s)
        if m and m.group(1) in by_target:
            by_target[m.group(1)].enum = scheme_keys()
        m = re.search(r'if\s*\(\s*(\w+)\s*!=\s*"([^"]*)"\s*\)', s)
        if m:
            k = next((k for t, k in by_target.items() if t == m.group(1)), None)
            if k is not None:
                k.enum = [m.group(2)]

        # --- ranges ---
        m = re.search(r"CheckRange\(\s*([\w.\[\]]+)\s*,\s*([^,]+),\s*([^)]+)\)", s)
        if m and norm_target(m.group(1)) in by_target:
            lo, hi = (float(consts.get(v.strip(), v.strip())) for v in (m.group(2), m.group(3)))
            by_target[norm_target(m.group(1))].range = (lo, hi, False, False)
        m = re.match(r"if\s*\(\s*([\w.]+)\s*(<=|<)\s*([-\d.eE+]+)\s*\)\s*$", s)
        if m and m.group(1) in by_target:
            nxt = lines[i + 1] if i + 1 < len(lines) else ""
            if "Error(" in nxt:
                lo = float(m.group(3))
                by_target[m.group(1)].range = (lo, math.inf, m.group(2) == "<=", False)
        #  a lower bound inside a negated condition that ends in an Error()
        if "Error(" in s:
            cond = " ".join(l.strip() for l in lines[max(0, i - 4):i])
            cond = cond[cond.rfind("if("):] if "if(" in cond else ""
            if "!(" in cond:
                for g in re.finditer(r"([\w.]+)\s*>=\s*([-\d.eE+]+)", cond):
                    if g.group(1) in by_target and by_target[g.group(1)].range is None:
                        by_target[g.group(1)].range = (float(g.group(2)), math.inf, False, False)
        m = re.match(r"if\s*\(\s*!\s*\(\s*([\w.]+)\s*>\s*([-\d.eE+]+)\s*&&\s*\1\s*<\s*([-\d.eE+]+)\s*\)\s*\)", s)
        if m and m.group(1) in by_target:
            by_target[m.group(1)].range = (float(m.group(2)), float(m.group(3)), True, True)

    # the getin calls that this did not understand would be missed silently
    ncalls = len(re.findall(r"->getin\(", text))
    return blocks, keys, ncalls


def code_model():
    header = parse_input_h(os.path.join(SRC, "IO", "Input.h"))
    blocks, keys, ncalls = parse_input_cpp(os.path.join(SRC, "IO", "Input.cpp"), header)
    return blocks, keys, ncalls


# ---------------------------------------------------------------------------
#  The manual
# ---------------------------------------------------------------------------
def tex_arg(text, pos):
    """The brace-balanced argument starting at text[pos] == '{'."""
    assert text[pos] == "{"
    depth = 0
    for j in range(pos, len(text)):
        if text[j] == "{" and text[j - 1] != "\\":
            depth += 1
        elif text[j] == "}" and text[j - 1] != "\\":
            depth -= 1
            if depth == 0:
                return text[pos + 1:j], j + 1
    raise ValueError("unbalanced braces at %d" % pos)


def detex(s):
    s = re.sub(r"\\(?:bk|allowbreak|,|\\)\s*", "", s)
    s = s.replace("\\_", "_").replace("{}", "").replace("~", " ")
    return s.strip()


def macro_args(text, name):
    """Every argument of \\name{...} in text, brace-balanced, de-TeXed."""
    out = []
    for m in re.finditer(r"\\%s\{" % name, text):
        arg, _ = tex_arg(text, m.end() - 1)
        out.append(detex(arg))
    return out


def resolve_block(name, blocks):
    if name in blocks:
        return {name}
    return {b for b in blocks if b.endswith("." + name)}


def table_blocks(caption, blocks):
    found = set()
    for k in macro_args(caption, "key"):
        found |= resolve_block(k, blocks)
    if not found and re.search(r"root level", caption, re.I):
        found = {""}
    # drop a block that only contains the ones the table is about
    return {b for b in found if not any(o != b and o.startswith(b + ".") for o in found)}


def split_rows(body):
    rows = []
    for r in re.split(r"\\\\\s*(?:\n|$)", body):
        r = r.strip()
        if r and "&" in r:
            rows.append(r)
    return rows


def split_cells(row):
    cells, depth, cur, i = [], 0, "", 0
    while i < len(row):
        c = row[i]
        if c == "\\" and i + 1 < len(row):
            cur += row[i:i + 2]
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
        if c == "&" and depth == 0:
            cells.append(cur.strip())
            cur = ""
        else:
            cur += c
        i += 1
    cells.append(cur.strip())
    return cells


class DocRow:
    def __init__(self, file, line, blocks, key, typ, default, meaning, section):
        self.file, self.line, self.blocks = file, line, blocks
        self.key, self.type, self.default = key, typ, default
        self.meaning, self.section = meaning, section

    def where(self):
        return "%s:%d" % (os.path.relpath(self.file, ROOT), self.line)


def section_text(text, pos):
    starts = [m.start() for m in re.finditer(r"\\(?:section|chapter)\{", text)]
    lo = max([s for s in starts if s <= pos] or [0])
    hi = min([s for s in starts if s > pos] or [len(text)])
    return text[lo:hi]


def doc_tables(blocks):
    rows = []
    for f in sorted(glob.glob(os.path.join(MANUAL, "*.tex"))):
        text = read(f)
        for m in re.finditer(r"\\begin\{(keytable|yzntable)\}", text):
            env = m.group(1)
            caption, pos = tex_arg(text, m.end())
            if env == "yzntable":
                _, pos = tex_arg(text, pos)        # the column spec
            end = text.index("\\end{%s}" % env, pos)
            body = text[pos:end]
            if env == "yzntable":
                head = body.split("\\midrule", 1)
                if len(head) < 2 or not all(w in head[0] for w in ("Key", "Type", "Default")):
                    continue
                body = head[1]
            tb = table_blocks(caption, blocks)
            sec = section_text(text, m.start())
            line0 = text.count("\n", 0, pos) + 1
            off = 0
            for r in split_rows(body):
                cells = split_cells(r)
                if len(cells) < 4 or not cells[0].startswith("\\key{"):
                    continue
                key = detex(tex_arg(cells[0], 4)[0])
                deflt = cells[2]
                line = line0 + body.count("\n", 0, body.find(r, off))
                off = body.find(r, off) + 1
                rows.append(DocRow(f, line, tb, key, detex(cells[1]), deflt,
                                   " & ".join(cells[3:]), sec))
    return rows


def doc_value(cell, typ):
    """The default column as a Python value; 'required' for \\reqd."""
    if "\\reqd" in cell:
        return "required"
    v = macro_args(cell, "val")
    if not v:
        return None
    s = v[0].strip()
    if typ == "string":
        return s.strip('"')
    if typ == "bool":
        return s == "true"
    if typ == "array":
        return [float(x) for x in s.strip("[]").split(",") if x.strip()]
    try:
        return float(s)
    except ValueError:
        return s


def same_value(a, b):
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same_value(x, y) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    if isinstance(a, float) and isinstance(b, float):
        return a == b or abs(a - b) <= 1e-12 * max(abs(a), abs(b))
    return a == b


def fmt(v):
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return "%g" % v
    if isinstance(v, list):
        return "[" + ", ".join(fmt(x) for x in v) + "]"
    if isinstance(v, str) and v != "required":
        return '"%s"' % v
    return str(v)


def doc_ranges(meaning):
    """The ranges a meaning column promises the parser checks.

    'range $[a:b]$', 'Range $(a:b)$' and 'at least N' are promises; 'at least
    N in practice', 'positive', 'non-negative' and the like describe what
    makes sense and are not (ch-ref-domain.tex says so)."""
    out = []
    for m in re.finditer(r"[Rr]ange\s*\$([\[(])\s*([^:$]+?)\s*:\s*([^$]+?)\s*([\])])\$", meaning):
        lo, hi = m.group(2), m.group(3)
        out.append((float(lo), float(hi), m.group(1) == "(", m.group(4) == ")"))
    for m in re.finditer(r"\b[Aa]t least\s+(-?[\d.]+)(?!\s+in practice)", meaning):
        out.append((float(m.group(1)), math.inf, False, False))
    return out


def check_tables(blocks, keys, problems):
    rows = doc_tables(blocks)
    documented = {}
    for r in rows:
        if not r.blocks:
            problems.append("%s: table of key %s: the caption names no block the parser reads"
                            % (r.where(), r.key))
            continue
        for b in r.blocks:
            path = b + "." + r.key if b else r.key
            documented.setdefault(path, []).append(r)
            k = keys.get(path)
            if k is None:
                problems.append("%s: %s is documented but the parser does not read it"
                                % (r.where(), path))
                continue
            if r.type != k.type:
                problems.append("%s: %s: type is %s in the code, %s in the manual"
                                % (r.where(), path, k.type, r.type))
            dv = doc_value(r.default, k.type)
            if k.required != (dv == "required"):
                problems.append("%s: %s: %s in the code, %s in the manual"
                                % (r.where(), path, "required" if k.required else "optional",
                                   "required" if dv == "required" else "optional"))
            elif not k.required and not same_value(dv, k.default):
                problems.append("%s: %s: default is %s in the code, %s in the manual"
                                % (r.where(), path, fmt(k.default), fmt(dv)))
            if k.enum:
                said = set(v.strip('"') for v in macro_args(r.meaning + " " + r.default, "val"))
                if not set(k.enum) <= said:
                    said |= set(v.strip('"') for v in macro_args(r.section, "val"))
                missing = [e for e in k.enum if e not in said]
                if missing:
                    problems.append("%s: %s: accepted values %s are not listed in the manual"
                                    % (r.where(), path, ", ".join(missing)))
            stated = doc_ranges(r.meaning)
            if k.range and not stated:
                problems.append("%s: %s: the code checks the range %s, the manual states none"
                                % (r.where(), path, fmt_range(k.range)))
            for st in stated:
                if not k.range or not same_range(st, k.range):
                    problems.append("%s: %s: the manual states the range %s, the code checks %s"
                                    % (r.where(), path, fmt_range(st),
                                       fmt_range(k.range) if k.range else "none"))
    for path, k in sorted(keys.items()):
        if path not in documented:
            problems.append("Input.cpp:%d: %s is read by the parser but in no key table of the manual"
                            % (k.line, path))
        elif len(documented[path]) > 1:
            problems.append("%s: %s is documented in more than one key table"
                            % (", ".join(r.where() for r in documented[path]), path))
    return rows


def fmt_range(r):
    lo, hi, lo_open, hi_open = r
    return "%s%s:%s%s" % ("(" if lo_open else "[", fmt(lo), fmt(hi), ")" if hi_open else "]")


def same_range(a, b):
    return (same_value(a[0], b[0]) and same_value(a[1], b[1])
            and a[2] == b[2] and a[3] == b[3])


# ---------------------------------------------------------------------------
#  ch-inputformat.tex: counts, required keys, the map of the top level
# ---------------------------------------------------------------------------
NUMBER_WORDS = {w: i for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen "
    "fourteen fifteen sixteen seventeen eighteen nineteen twenty".split())}


def number(word):
    w = word.lower()
    return int(w) if w.isdigit() else NUMBER_WORDS.get(w)


def check_inputformat(blocks, keys, problems):
    f = os.path.join(MANUAL, "ch-inputformat.tex")
    text = read(f)
    rel = os.path.relpath(f, ROOT)
    nblocks = len(blocks)
    with_keys = len({k.block for k in keys.values()})
    nkeys = len(keys)
    m = re.search(r"(\d+)\s+block\s+paths,\s+(\d+)\s+of\s+which\s+carry\s+keys,\s+and\s+(\d+)\s+key\s+entries", text)
    if not m:
        problems.append("%s: the sentence giving the number of block paths and keys is gone" % rel)
    else:
        said = tuple(int(x) for x in m.groups())
        if said != (nblocks, with_keys, nkeys):
            problems.append("%s: counts %d block paths, %d with keys, %d keys; the parser has %d, %d, %d"
                            % ((rel,) + said + (nblocks, with_keys, nkeys)))

    req = sorted(p for p, k in keys.items() if k.required)
    m = re.search(r"(\w+) keys are required outright", text)
    if not m or number(m.group(1)) != len(req):
        problems.append("%s: the number of required keys is %d (%s)" % (rel, len(req), ", ".join(req)))
    para = text[m.start():text.find("\n\n", m.start())] if m else ""
    for p in req:
        if p.split(".")[-1] not in macro_args(para, "key"):
            problems.append("%s: the required key %s is not in the list of required keys" % (rel, p))

    shapes_i = [b for b, k in blocks.items() if k == "repeatable" and b.startswith("Initial.SurfaceShape.")]
    shapes_s = [b for b, k in blocks.items() if k == "repeatable" and b.startswith("Substrate.")]
    m = re.search(r"These are the (\w+) shape blocks,\s*(\w+) under\s+\\key\{Initial\.SurfaceShape\}\s+"
                  r"and (\w+) under \\key\{Substrate\}", text)
    if not m or (number(m.group(1)), number(m.group(2)), number(m.group(3))) != \
            (len(shapes_i) + len(shapes_s), len(shapes_i), len(shapes_s)):
        problems.append("%s: the parser has %d shape blocks, %d under Initial.SurfaceShape and %d under Substrate"
                        % (rel, len(shapes_i) + len(shapes_s), len(shapes_i), len(shapes_s)))
    others = sorted(b for b, k in blocks.items() if k == "repeatable"
                    and b not in shapes_i and b not in shapes_s)
    sec = section_text(text, text.find("label{sec:inputformat:repeat}"))
    for b in others:
        if b not in macro_args(sec, "key"):
            problems.append("%s: the repeatable block %s is not named in the section on repeated blocks"
                            % (rel, b))

    # the map of the top level
    m = re.search(r"\\begin\{yzntable\}\{Every top-level block", text)
    if m:
        body = text[m.start():text.index("\\end{yzntable}", m.start())].split("\\midrule", 1)[1]
        seen = set()
        for r in split_rows(body):
            c = split_cells(r)
            if not c[0].startswith("\\key{"):
                continue
            b = detex(tex_arg(c[0], 4)[0])
            seen.add(b)
            kind = blocks.get(b)
            if kind is None:
                problems.append("%s: the map names a block %s the parser does not read" % (rel, b))
                continue
            said_req = c[1].strip() == "yes"
            said_rep = c[2].strip() == "yes"
            if said_req != (kind == "required") or said_rep != (kind == "repeatable"):
                problems.append("%s: map of the file: %s is %s in the code, required=%s repeats=%s in the manual"
                                % (rel, b, kind, c[1].strip(), c[2].strip()))
        for b in blocks:
            if b and "." not in b and b not in seen:
                problems.append("%s: the top-level block %s is missing from the map" % (rel, b))


# ---------------------------------------------------------------------------
#  Every mention of a key or block in the manual and in the Markdown documents
# ---------------------------------------------------------------------------
def known_names(blocks, keys):
    full = set(blocks) | set(keys)
    leaf = {p.split(".")[-1] for p in full if p}
    return full, leaf


def resolves(name, full, leaf):
    if name in full or name in leaf:
        return True
    return any(p.endswith("." + name) for p in full)


def check_mentions(blocks, keys, problems):
    full, leaf = known_names(blocks, keys)
    for f in sorted(glob.glob(os.path.join(MANUAL, "*.tex"))):
        text = read(f)
        for m in re.finditer(r"\\key\{", text):
            arg, _ = tex_arg(text, m.end() - 1)
            name = detex(arg)
            name = re.split(r"\s*=", name)[0].strip()
            name = re.sub(r"\[\d*\]$", "", name)
            if not name or re.search(r"[<*\\$ ]", name) or name.startswith("["):
                continue
            if not resolves(name, full, leaf):
                problems.append("%s:%d: \\key{%s} names no block or key the parser reads"
                                % (os.path.relpath(f, ROOT), text.count("\n", 0, m.start()) + 1, name))
    tops = sorted({b.split(".")[0] for b in blocks if b})
    pat = re.compile(r"(?<![\w./-])((?:%s)(?:\.[A-Za-z_]\w*)+)" % "|".join(tops))
    mds = [os.path.join(ROOT, "README.md"), os.path.join(ROOT, "CONTRIBUTING.md"),
           os.path.join(ROOT, "cases", "README.md")] + sorted(glob.glob(os.path.join(ROOT, "docs", "*.md")))
    for f in mds:
        if not os.path.exists(f):
            continue
        for n, line in enumerate(read(f).splitlines(), 1):
            for m in pat.finditer(line):
                name = m.group(1)
                if name.endswith((".toml", ".cpp", ".h", ".py", ".tex", ".md", ".dat")):
                    continue
                if not resolves(name, full, leaf):
                    problems.append("%s:%d: %s names no block or key the parser reads"
                                    % (os.path.relpath(f, ROOT), n, name))


# ---------------------------------------------------------------------------
#  The key tables of the Markdown documents (docs/numerics.md) and the map
#  of the top-level blocks in README.md
# ---------------------------------------------------------------------------
def md_value(cell, typ):
    if "required" in cell:
        return "required"
    m = re.match(r"\s*`([^`]*)`", cell)
    if not m:
        return None
    v = m.group(1)
    return doc_value("\\val{%s}" % v, typ)


def check_markdown(blocks, keys, problems):
    for f in sorted(glob.glob(os.path.join(ROOT, "docs", "*.md"))):
        rel = os.path.relpath(f, ROOT)
        lines = read(f).splitlines()
        current = []
        for n, line in enumerate(lines, 1):
            if line.startswith("#"):
                current = [b for b in re.findall(r"`([\w.]+)`", line) if b in blocks]
                continue
            if not current or not line.startswith("|"):
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) < 3 or not re.fullmatch(r"`\w+`", cells[0]):
                continue
            name = cells[0].strip("`")
            for b in current:
                path = b + "." + name
                k = keys.get(path)
                where = "%s:%d" % (rel, n)
                if k is None:
                    problems.append("%s: %s is documented but the parser does not read it" % (where, path))
                    continue
                dv = md_value(cells[1], k.type)
                if k.required != (dv == "required"):
                    problems.append("%s: %s: %s in the code, %s in the document"
                                    % (where, path, "required" if k.required else "optional",
                                       "required" if dv == "required" else "optional"))
                elif not k.required and dv is not None and not same_value(dv, k.default):
                    problems.append("%s: %s: default is %s in the code, %s in the document"
                                    % (where, path, fmt(k.default), fmt(dv)))
                desc = " ".join(cells[2:])
                stated = [(float(m.group(2)), float(m.group(3)), m.group(1) == "(", m.group(4) == ")")
                          for m in re.finditer(r"[Rr]ange `([\[(])([^:`]+):([^`]+)([\])])`", desc)]
                stated += [(float(m.group(1)), math.inf, False, False)
                           for m in re.finditer(r"\b[Aa]t least\s+(-?[\d.]+)(?!\s+in practice)", desc)]
                if k.range and not stated:
                    problems.append("%s: %s: the code checks the range %s, the document states none"
                                    % (where, path, fmt_range(k.range)))
                for st in stated:
                    if not k.range or not same_range(st, k.range):
                        problems.append("%s: %s: the document states the range %s, the code checks %s"
                                        % (where, path, fmt_range(st),
                                           fmt_range(k.range) if k.range else "none"))
                if k.enum and len(k.enum) > 1:
                    said = set(re.findall(r"`\"?([^`\"]+)\"?`", desc + " " + cells[1]))
                    missing = [e for e in k.enum if e not in said]
                    if missing:
                        problems.append("%s: %s: accepted values %s are not listed in the document"
                                        % (where, path, ", ".join(missing)))

    f = os.path.join(ROOT, "README.md")
    text = read(f)
    m = re.search(r"^\| Block \| Contents \|\n\|[-|]+\|\n((?:\|.*\|\n)+)", text, re.M)
    if not m:
        problems.append("README.md: the table of the top-level blocks is gone")
        return
    rows = {r.split("|")[1].strip().strip("`"): r for r in m.group(1).splitlines()}
    for b in blocks:
        if b and "." not in b and b not in rows:
            problems.append("README.md: the top-level block %s is missing from the table" % b)
    for b in rows:
        if b not in blocks:
            problems.append("README.md: the table names a block %s the parser does not read" % b)
    if "Boundary" in rows:
        for t in keys["Boundary.XminBC.Type"].enum:
            if "`%s`" % t not in rows["Boundary"]:
                problems.append("README.md: the boundary type %s is missing from the table" % t)


# ---------------------------------------------------------------------------
#  TOML: the listings of the manual and the shipped cases
# ---------------------------------------------------------------------------
def toml_paths(doc, prefix=""):
    """(path, is_block) for every key and table of a parsed TOML document."""
    for k, v in doc.items():
        p = prefix + "." + k if prefix else k
        if isinstance(v, dict):
            yield p, True
            yield from toml_paths(v, p)
        elif isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
            yield p, True
            for x in v:
                yield from toml_paths(x, p)
        else:
            yield p, False


def check_toml_doc(where, doc, blocks, keys, problems):
    for p, is_block in toml_paths(doc):
        if is_block and p not in blocks:
            problems.append("%s: [%s] is not a block the parser reads" % (where, p))
        elif not is_block and p not in keys:
            problems.append("%s: %s is not a key the parser reads" % (where, p))
        elif is_block and blocks[p] != "repeatable":
            node = doc
            for part in p.split("."):
                node = node[part] if isinstance(node, dict) else None
            if isinstance(node, list):
                problems.append("%s: [[%s]] repeats a block that may appear only once" % (where, p))


def is_complete(doc, blocks):
    """Whether a listing is a whole input file rather than a fragment: it has
    every block that is required and sits under required blocks only."""
    present = {p for p, is_block in toml_paths(doc) if is_block}

    def always_needed(b):
        parts = b.split(".")
        return all(blocks.get(".".join(parts[:i])) == "required" for i in range(1, len(parts) + 1))

    return "SolverName" in doc and all(b in present for b in blocks if b and always_needed(b))


def manual_listings():
    out = []
    for f in sorted(glob.glob(os.path.join(MANUAL, "*.tex"))):
        text = read(f)
        for m in re.finditer(r"\\begin\{toml\}(?:\[[^\]]*\])?\n(.*?)\\end\{toml\}", text, re.S):
            out.append(("%s:%d" % (os.path.relpath(f, ROOT), text.count("\n", 0, m.start()) + 1),
                        m.group(1)))
    return out


def case_files():
    return sorted(glob.glob(os.path.join(ROOT, "cases", "**", "*.toml"), recursive=True))


def check_toml(blocks, keys, problems):
    if tomllib is None:
        print("doc_keys: Python %s has no tomllib; the TOML checks are skipped" % sys.version.split()[0])
        return [], []
    complete = []
    for where, src in manual_listings():
        try:
            doc = tomllib.loads(src)
        except tomllib.TOMLDecodeError as e:
            problems.append("%s: the TOML listing does not parse: %s" % (where, e))
            continue
        check_toml_doc(where, doc, blocks, keys, problems)
        if is_complete(doc, blocks):
            complete.append((where, src))
    used = set()
    for f in case_files():
        rel = os.path.relpath(f, ROOT)
        try:
            with open(f, "rb") as fh:
                doc = tomllib.load(fh)
        except tomllib.TOMLDecodeError as e:
            problems.append("%s: does not parse: %s" % (rel, e))
            continue
        check_toml_doc(rel, doc, blocks, keys, problems)
        used |= {p for p, ib in toml_paths(doc)}
    return complete, used


def check_unused_table(blocks, keys, used, problems):
    """The table of what no shipped case uses must name exactly that."""
    f = os.path.join(MANUAL, "ch-ref-numerics.tex")
    text = read(f)
    m = re.search(r"\\begin\{yzntable\}\{Parser features not exercised", text)
    if not m or not used:
        return
    body = text[m.start():text.index("\\end{yzntable}", m.start())].split("\\midrule", 1)[1]
    named = set()
    for r in split_rows(body):
        c = split_cells(r)
        bl = c[0]
        if bl.strip() == "root":
            bpaths = [""]
        elif "Both" in bl and "LinearSolver" in bl:
            bpaths = sorted(b for b in blocks if b.endswith(".LinearSolver"))
        else:
            bpaths = [detex(tex_arg(bl, bl.index("{"))[0])]
        for b in bpaths:
            listed = macro_args(re.split(r";|although", c[1])[0], "key")
            if not listed and "whole block" in c[1]:
                named.add(b)
            for k in listed:
                p = (b + "." if b else "") + k
                if p in keys or p in blocks:
                    named.add(p)
                else:
                    sub = [q for q in list(keys) + list(blocks) if q.startswith(b + ".") and q.endswith("." + k)]
                    named |= set(sub) if sub else {p}
    #  a full path named in the prose of the section counts as listed too
    tail = section_text(text, m.start())
    tail = tail[tail.index("\\end{yzntable}") if "\\end{yzntable}" in tail else 0:]
    named |= {k for k in macro_args(tail, "key") if k in keys}
    rel = os.path.relpath(f, ROOT)
    for p in sorted(named):
        if p in used:
            problems.append("%s: table of unused keys: %s is used by a case in cases/" % (rel, p))
    unused = set()
    for p in list(keys) + list(blocks):
        if not p or p in used:
            continue
        parent = p.rsplit(".", 1)[0] if "." in p else ""
        if parent and parent not in used:
            continue   # naming a block no case uses covers what is inside it
        unused.add(p)
    #  the keys of an unused block are covered by naming the block
    for p in sorted(unused):
        covered = p in named or any(p.startswith(n + ".") for n in named)
        if not covered:
            problems.append("%s: table of unused keys: %s is used by no case but not listed" % (rel, p))
    m = re.search(r"every block path and key used in\s+the\s+(\d+)\s+input\s+files", text)
    ncases = len(case_files())
    if not m or int(m.group(1)) != ncases:
        problems.append("%s: the number of input files the cross-check covers is not given as %d"
                        % (rel, ncases))


# ---------------------------------------------------------------------------
#  The binary
# ---------------------------------------------------------------------------
def toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, list):
        return "[" + ", ".join(toml_value(x) for x in v) + "]"
    return '"%s"' % v


REQUIRED_VALUES = {"SolverName": None, "JobName": "everykey", "Grid.Dimension": "xh",
                   "Grid.Nx": 8.0, "Grid.Xmin": 0.0, "Grid.Xmax": 1e-3,
                   "Initial.Concentration.Value": 0.5, "Physics.ConcField.Model": "WM",
                   "Numerics.InitialDt": 1e-6, "Output.ProbeList.Probe.X": 0.0,
                   "Output.FieldData.LogSchedule.Start": 1e-6, "Output.FieldData.LogSchedule.End": 1.0,
                   "Output.FieldData.LogSchedule.PerDecade": 10.0}
#  Defaults that, all present at once, make a combination CheckParameters
#  refuses: an EdgeBump (at Xmax by default) needs a pinned edge there.
EVERY_KEY_OVERRIDES = {"Boundary.XmaxBC.Type": "Pinning"}
#  Blocks that are an alternative to keys given at their default: a
#  LogSchedule replaces Output.FieldData.TimeInterval, and may not be given
#  with it (CheckParameters)
EVERY_KEY_SKIP = {"Output.FieldData.LogSchedule"}


def every_key_file(blocks, keys):
    lines = []
    by_block = {}
    for k in keys.values():
        by_block.setdefault(k.block, []).append(k)
    for b in sorted(blocks, key=lambda b: (b != "", b)):
        if b in EVERY_KEY_SKIP:
            continue
        ks = by_block.get(b, [])
        if b:
            br = "[[%s]]" if blocks[b] == "repeatable" else "[%s]"
            lines.append("\n" + br % b)
        for k in sorted(ks, key=lambda k: k.name):
            if k.required:
                v = REQUIRED_VALUES.get(k.path, k.enum[0] if k.enum else 0.0)
                if v is None:
                    v = k.enum[0]
            else:
                v = k.default
            v = EVERY_KEY_OVERRIDES.get(k.path, v)
            if k.type == "int" and isinstance(v, float):
                v = int(v)
            if v is None or v == []:
                continue        # an empty list is not given, but left out
            lines.append("%s = %s" % (k.name, toml_value(v)))
    return "\n".join(lines) + "\n"


def run_check(exe, path, cwd):
    p = subprocess.run([exe, "--check-input", path], cwd=cwd, capture_output=True, text=True,
                       timeout=120)
    return p.returncode, p.stdout + p.stderr


def check_binary(exe, workdir, blocks, keys, complete, problems):
    os.makedirs(workdir, exist_ok=True)
    f = os.path.join(workdir, "everykey.toml")
    with open(f, "w") as fh:
        fh.write(every_key_file(blocks, keys))
    rc, out = run_check(exe, f, workdir)
    if rc != 0 or "Input check passed" not in out:
        problems.append("the binary refuses a file with every block and key at its default:\n"
                        + out[-2000:])
    else:
        #  ... and it must be the binary that reads them: one misspelled key more
        #  has to be refused, or the check above proves nothing
        g = os.path.join(workdir, "misspelled.toml")
        with open(g, "w") as fh:
            fh.write(every_key_file(blocks, keys).replace("\nNx = ", "\nNxx = 1\nNx = ", 1))
        rc, out = run_check(exe, g, workdir)
        if rc == 0 or "Grid.Nxx" not in out:
            problems.append("the binary does not refuse an unknown key Grid.Nxx:\n" + out[-2000:])
    for c in case_files():
        rc, out = run_check(exe, c, workdir)
        if rc != 0:
            problems.append("%s: yuzen --check-input fails:\n%s" % (os.path.relpath(c, ROOT), out[-1500:]))
    for n, (where, src) in enumerate(complete):
        g = os.path.join(workdir, "listing%02d.toml" % n)
        with open(g, "w") as fh:
            fh.write(src)
        rc, out = run_check(exe, g, workdir)
        if rc != 0:
            problems.append("%s: the listing fails yuzen --check-input:\n%s" % (where, out[-1500:]))


# ---------------------------------------------------------------------------
def print_keys(blocks, keys):
    for b in sorted(blocks):
        print("[%s]  %s" % (b or "root", blocks[b]))
        for k in sorted((k for k in keys.values() if k.block == b), key=lambda k: k.line):
            extra = []
            if k.enum:
                extra.append("one of " + "|".join(k.enum))
            if k.range:
                extra.append("range " + fmt_range(k.range))
            if k.absent is not None:
                extra.append("absent block: " + fmt(k.absent))
            print("    %-22s %-7s %-16s %s" % (k.name, k.type,
                                               "required" if k.required else fmt(k.default),
                                               "; ".join(extra)))
    print("%d block paths, %d with keys, %d keys"
          % (len(blocks), len({k.block for k in keys.values()}), len(keys)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", help="the yuzen binary; without it the binary is not asked")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--list", action="store_true", help="print the keys the parser reads and stop")
    a = ap.parse_args()

    blocks, keys, ncalls = code_model()
    if a.list:
        print_keys(blocks, keys)
        return 0
    problems = []
    if ncalls != len(keys):
        problems.append("Input.cpp has %d getin() calls but %d keys were understood; extend doc_keys.py"
                        % (ncalls, len(keys)))
    for k in keys.values():
        if k.type not in ("string", "double", "int", "bool", "array"):
            problems.append("Input.cpp:%d: type of %s not understood (%s)" % (k.line, k.path, k.ctype))
        if not k.required and k.default is None:
            problems.append("Input.cpp:%d: default of %s not understood" % (k.line, k.path))
        if k.range and not k.required and isinstance(k.default, float):
            lo, hi, lo_open, hi_open = k.range
            ok = (k.default > lo if lo_open else k.default >= lo) and \
                 (k.default < hi if hi_open else k.default <= hi)
            if not ok:
                problems.append("Input.cpp:%d: the default of %s is outside its own range" % (k.line, k.path))
        if k.enum and not k.required and k.default not in k.enum:
            problems.append("Input.cpp:%d: the default of %s is not one of its accepted values"
                            % (k.line, k.path))

    check_tables(blocks, keys, problems)
    check_inputformat(blocks, keys, problems)
    check_mentions(blocks, keys, problems)
    check_markdown(blocks, keys, problems)
    complete, used = check_toml(blocks, keys, problems)
    check_unused_table(blocks, keys, used, problems)
    if a.exe:
        wd = a.workdir or os.path.join(os.getcwd(), "doc_keys")
        shutil.rmtree(wd, ignore_errors=True)
        check_binary(os.path.abspath(a.exe), wd, blocks, keys, complete, problems)

    for p in problems:
        print("MISMATCH " + p)
    print("doc_keys: %d block paths, %d keys, %d cases, %d manual listings; %d mismatch%s"
          % (len(blocks), len(keys), len(case_files()), len(manual_listings()),
             len(problems), "" if len(problems) == 1 else "es"))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
