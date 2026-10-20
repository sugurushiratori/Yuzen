# The Yuzen manual

Source for the full manual — getting started, the complete input-file
reference, the theory, the numerical method, and the verification record — in
two editions built from **the same chapters**, so they cannot drift:

| | |
|---|---|
| `yuzen-manual.pdf` | LuaLaTeX, A4, for reading and printing |
| `yuzen-manual.html` | one self-contained file, full-text search, dark mode, no network needed |
| `yuzen-manual-ja.pdf`, `yuzen-manual-ja.html` | the same two in Japanese, from `chapters-ja/` (STYLE.md, "The Japanese edition") |

**Neither built file is kept in git** — the HTML alone is close to six
megabytes, and committing it would grow the history by that much every time a
chapter changes. `.github/workflows/release.yml` builds both when a `v*` tag is
pushed and attaches them to the GitHub release; that workflow also runs
`make check`, so a broken cross-reference stops the release rather than
shipping in it. To get the files without a tag, run `make` here, or run the
workflow by hand from the Actions tab and take the artifact.

## Building

```sh
make                    # both editions, English and Japanese
make en / make ja       # one language
make html               # the HTML edition only, no TeX needed
make check              # fail on undefined references and overfull lines
make watch              # rebuild the PDF on every change
make doctor             # which node and mathjax-full will be used, and whether
                        # that install is intact
```

`make` needs a TeX installation and `make html` needs Node; each says which
tool is missing, and how to install it, if it cannot find one. There is no need
to install either merely to read the manual — the PDF and the HTML edition are
attached to every release.

The HTML edition pre-renders every formula to SVG with MathJax, which runs
under Node:

```sh
sudo apt install nodejs npm         # or: brew install node
cd docs/manual && npm install mathjax-full
```

`mathjax-full` is found whether it is installed here, at the top of the
repository, or globally.

**On WSL, do not `npm install` into a repository on a Windows drive.** A
checkout under `/mnt/c/...` is on DrvFs, where npm regularly writes package
files incompletely — a truncated `.js`, or a `package.json` that is not JSON at
all. npm reports success, and node then fails deep inside MathJax or inside its
own module resolver. Reinstalling in the same place usually fails the same way.
`tohtml.py` checks the package before it renders anything and says so plainly
if it is damaged.

`make doctor` reports which install would be used and whether it is intact.

Install on the Linux filesystem and point node at it:

```sh
(cd ~ && npm install mathjax-full)
NODE_PATH=~/node_modules make html
```

Setting `NODE_PATH` once in `~/.bashrc` saves repeating it.

This works even when a damaged `node_modules` is still sitting in
`docs/manual`: the script uses the first install that passes its check and
loads it by absolute path, so the broken one nearer the script cannot win.
Keeping the whole repository on the Linux filesystem avoids the problem, and
is much faster to build on besides.

**LuaLaTeX is required** for the PDF: the preamble uses `fontspec` and
`unicode-math`, so `pdflatex` and `xelatex` will not build it. The HTML edition
needs Python 3 and Node with `mathjax-full`; regenerating the three shape
figures additionally needs matplotlib, and the checked-in SVGs mean that is
only required when a figure changes.

### What you need

| | |
|---|---|
| Engine | LuaLaTeX (TeX Live 2021 or later) |
| Packages | `koma-script`, `fontspec`, `unicode-math`, `amsmath`, `mathtools`, `siunitx`, `booktabs`, `longtable`, `tabularx`, `enumitem`, `caption`, `listings`, `tcolorbox`, `microtype`, `ragged2e`, `hyperref`, `cleveref`, `geometry` |
| Fonts | TeX Gyre Pagella and Pagella Math, TeX Gyre Heros, Inconsolata |

On Debian and Ubuntu:

```sh
sudo apt install texlive-luatex texlive-latex-recommended texlive-latex-extra \
                 texlive-science texlive-fonts-extra fonts-texgyre fonts-texgyre-math \
                 fonts-inconsolata latexmk
sudo apt install texlive-lang-japanese     # for the Japanese edition (LuaTeX-ja)
```

`fonts-texgyre` carries TeX Gyre Pagella and Heros. `apt` installs it only as
a recommendation of `texlive-fonts-recommended`, so with
`--no-install-recommends` it has to be named, or LuaLaTeX stops with
``The font "TeXGyrePagella" cannot be found``.

`texlive-latex-recommended` is what carries `scrbook.cls`, `fontspec` and
`unicode-math`; without it the build stops at ``File `scrbook.cls' not found``.

Use `latexmk` or `make`, not `pdflatex`. The document requires LuaTeX and says
so if another engine is used.

The two logo figures come from `../logo/` via `cairosvg` (`pip install cairosvg`) and
`make figures` regenerates them; the three shape figures come from
`figures/make_shapes.py` and need matplotlib. The checked-in PDFs mean you need
neither unless you are changing a figure.

## Layout

```
yuzen-manual.tex     the main file: document class, parts, \include list
yuzenman.sty         everything else: fonts, colors, macros, environments
STYLE.md             the writing contract — read this before editing a chapter
chapters/*.tex       one file per chapter, body text only
tohtml.py            the HTML builder: the same chapters, a different renderer
figures/             logo PDFs (from ../logo/*.svg) and the shape figures,
                     generated by figures/make_shapes.py
```

## Editing

`STYLE.md` is the contract: the semantic markup (`\key`, `\val`, `\file`,
`\code`, `\env`), the `keytable` and `yzntable` environments, the `warn` /
`note` / `tip` boxes, the label scheme, and the maths shorthands. Follow it —
uniform markup is what lets the whole reference be restyled or re-checked
mechanically.

Two rules that are easy to get wrong:

* Put the `\label` **inside** the caption argument,
  `\begin{keytable}{Keys of ...\label{tab:...}}`. A `\label` as the first thing
  in the table body pushes the first cell down by one line.
* Long key paths do not break. Insert `\bk` at a sensible point:
  `\key{Output.\bk FieldData.\bk TimeInterval}`.

Run `make check` before committing, and `make html` if you changed a chapter —
`tohtml.py` understands exactly the markup `STYLE.md` defines and **stops with
an error** rather than silently dropping a construct it does not know. If it
complains about an unhandled command, either use the documented markup or teach
the converter about the new one.

### How the HTML edition is built

`tohtml.py` reads `yuzen-manual.tex` for the chapter order, numbers everything
the way LaTeX would, and renders the chapters itself. Two choices are worth
knowing:

* **Mathematics is pre-rendered to SVG** with MathJax at build time, using the
  same macro definitions as `yuzenman.sty`. The page therefore shows every
  equation correctly with JavaScript disabled and with no network access, and
  the SVG inherits `currentColor`, so equations follow the light/dark theme.
* **Figures are inlined as SVG** from `figures/*.svg`, which
  `figures/make_shapes.py` writes alongside the PDFs with text converted to
  outlines. Converting the PDF instead would work but loses the text.

The only JavaScript in the page is the chapter router and the full-text search,
whose index is built at the same time, one entry per chapter and per section.

## Where the content comes from

The manual is the authoritative description. `../theory.md`, `../numerics.md`
and `../validation.md` are shorter companions that follow it; **if you change
the physics or the numerics, change the manual and, where they touch the same
point, the companions too**, or they will drift. The input-file reference in
Part II was read out of the parser source (`yuzen/src/IO/`), not out of the
prose documents or the example cases, and where the two disagreed the manual
follows the source and says so in a footnote.
