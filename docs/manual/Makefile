# ---------------------------------------------------------------------------
#  The Yuzen manual.  LuaLaTeX only -- the preamble uses fontspec and
#  unicode-math, so pdflatex and xelatex will not build it.
#
#    make          build both languages: yuzen-manual{,-ja}.{pdf,html}
#    make en       the English edition only (PDF and HTML)
#    make ja       the Japanese edition only (PDF and HTML)
#    make watch    rebuild the English PDF on every change
#    make html     build the searchable HTML editions only
#    make figures  regenerate the logo PDFs from docs/logo/*.svg
#    make clean    remove build by-products, keep the PDF
#    make distclean  remove the PDF and the HTML too
#
#  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
#  SPDX-License-Identifier: Apache-2.0
# ---------------------------------------------------------------------------
MAIN    := yuzen-manual
MAINJA  := yuzen-manual-ja
LATEXMK := latexmk -lualatex -interaction=nonstopmode -file-line-error
LOGO    := ../logo
FIGS    := figures/yuzen-mark.pdf figures/yuzen-lockup-h.pdf
SVGS    := figures/shapes-1d.svg figures/shapes-createtype.svg \
           figures/film-geometry.svg
SCALING := figures/parallel-scaling.svg
#  The example-result figures live in docs/figures/ because the repository
#  README shows them too. They are committed, and regenerating them needs the
#  solver and a set of finished runs, so they are dependencies of the manual
#  and never a rule: see docs/figures/README.md.
EXFIGS  := $(wildcard ../figures/*.pdf)
EXSVGS  := $(wildcard ../figures/*.svg)

.PHONY: all en ja html watch figures clean distclean check doctor need-latexmk need-python

all: en ja
en: $(MAIN).pdf $(MAIN).html
ja: $(MAINJA).pdf $(MAINJA).html

#  Say which tool is missing and how to get it. Without this the build stops
#  at "latexmk: No such file or directory" and an exit code, which does not
#  tell the reader that a TeX installation is what they are missing.
#  These sit after `all` on purpose: the first target in the file is what a
#  bare `make` builds.
need-latexmk:
	@command -v latexmk >/dev/null 2>&1 || { \
	  echo "make: latexmk was not found."; \
	  echo "  The PDF needs a LuaLaTeX installation. On Debian and Ubuntu:"; \
	  echo "    sudo apt install texlive-luatex texlive-latex-recommended \\"; \
	  echo "                     texlive-latex-extra texlive-science \\"; \
	  echo "                     texlive-fonts-extra fonts-texgyre fonts-texgyre-math \\"; \
	  echo "                     fonts-inconsolata latexmk"; \
	  echo "  and, for the Japanese edition, texlive-lang-japanese."; \
	  echo "  On macOS: MacTeX. On Windows: TeX Live or MiKTeX."; \
	  echo "  See README.md. To build only the HTML edition, run 'make html'."; \
	  exit 1; }
	@command -v lualatex >/dev/null 2>&1 || { \
	  echo "make: latexmk is installed but lualatex is not."; \
	  echo "  Install texlive-luatex (Debian/Ubuntu) or the full TeX distribution."; \
	  exit 1; }

#  Report which node and which mathjax-full install would be used, and
#  whether that install is intact. The first thing to run when make html
#  fails inside node.
doctor:
	@python3 tohtml.py --check-install

need-python:
	@command -v python3 >/dev/null 2>&1 || { \
	  echo "make: python3 was not found; it builds the HTML edition."; exit 1; }

$(MAIN).pdf: $(MAIN).tex yuzenman.sty $(wildcard chapters/*.tex) $(FIGS) $(SCALING) \
             $(EXFIGS) | need-latexmk
	$(LATEXMK) $(MAIN).tex

$(MAINJA).pdf: $(MAINJA).tex yuzenman.sty $(wildcard chapters-ja/*.tex) $(FIGS) $(SCALING) \
               $(EXFIGS) | need-latexmk
	$(LATEXMK) $(MAINJA).tex

# The HTML edition is generated from the same chapters, so the two cannot
# drift. Maths is pre-rendered to SVG with MathJax, so the page needs no
# network and no JavaScript to show an equation; the JavaScript it does carry
# provides the full-text search.
html: $(MAIN).html $(MAINJA).html

$(MAIN).html: tohtml.py $(MAIN).tex $(wildcard chapters/*.tex) $(SVGS) $(SCALING) \
              $(EXSVGS) | need-python
	python3 tohtml.py

$(MAINJA).html: tohtml.py $(MAINJA).tex $(wildcard chapters-ja/*.tex) $(SVGS) $(SCALING) \
                $(EXSVGS) | need-python
	python3 tohtml.py --lang ja

$(SVGS): figures/make_shapes.py
	python3 figures/make_shapes.py

$(SCALING): figures/make_scaling.py ../../tools/scaling/results/scaling-mc076.tsv
	python3 figures/make_scaling.py

watch: $(FIGS) | need-latexmk
	$(LATEXMK) -pvc $(MAIN).tex

figures: $(FIGS)

figures/%.pdf: $(LOGO)/%.svg
	@command -v cairosvg >/dev/null 2>&1 || { \
	  echo "make: cairosvg was not found (pip install cairosvg)."; \
	  echo "  It converts the logo SVGs to PDF. The converted PDFs are"; \
	  echo "  committed, so you only need it when changing a logo."; \
	  exit 1; }
	@mkdir -p figures
	cairosvg $< -o $@

# Fail loudly on anything a reader would notice: undefined references and
# overfull lines wider than 10pt, in either edition; and the Japanese chapters
# must carry the same labels, references, keys, values, file names, numbers,
# maths and code as the English ones (../../tools/manual_ja_check.py).
check: $(MAIN).pdf $(MAINJA).pdf
	@for m in $(MAIN) $(MAINJA); do \
	  test -f $$m.log || { echo "check: $$m.log is missing -- run 'make' first"; exit 1; }; \
	  ! grep -q "LaTeX Warning: Reference" $$m.log || \
	     { grep "LaTeX Warning: Reference" $$m.log; exit 1; }; \
	  grep -o 'Overfull \\hbox ([0-9.]*pt' $$m.log \
	     | sed 's/.*(\([0-9.]*\)pt/\1/' \
	     | awk -v m=$$m '$$1 > 10 { print m ": overfull: " $$1 "pt"; bad = 1 } END { exit bad }' || exit 1; \
	done
	@python3 ../../tools/manual_ja_check.py
	@echo "check: no undefined references, no overfull box wider than 10pt, the two editions agree"

clean:
	$(LATEXMK) -c $(MAIN).tex
	$(LATEXMK) -c $(MAINJA).tex

distclean:
	$(LATEXMK) -C $(MAIN).tex
	$(LATEXMK) -C $(MAINJA).tex
	rm -f $(MAIN).html $(MAINJA).html
