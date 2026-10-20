#!/usr/bin/env python3
"""Build every Yuzen logo asset from one geometry.

    python3 docs/logo/make_logo.py            # the SVGs only
    python3 docs/logo/make_logo.py --raster   # also png/, favicon.ico, the GIF
                                              # and the manual's PDFs

The mark is a Y with a celadon glaze on it, lit from the upper left: a light
gradient from the top of the left arm to the foot of the stem, a highlight
along the left arm, and a darker run of thicker glaze down the right arm and
the right side of the stem. Its outline is built here from straight edges and
circular arcs (the PARAMS below); the outline of the run is fixed (RUN).

Every file this writes is generated: change the geometry or a color here and
run it again, rather than editing an SVG by hand. The rasters need cairosvg
and Pillow (pip install cairosvg pillow), and the GIF rsvg-convert (librsvg).

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MANUAL_FIGS = os.path.join(HERE, "..", "manual", "figures")

# ---------------------------------------------------------------------------
# Geometry, in the units of the mark's 640 x 640 canvas
# ---------------------------------------------------------------------------
SIZE = 640.0
DX, DY = 42.5, 35.0          # places the drawing so the clear space is even

#  The four edges of the two arms, x = a*y + b (before the DX, DY shift), and
#  the side of each on which the letter lies (+1: to the right of it).
#  The arms taper towards the junction; the right one is a little steeper.
EDGES = {
    "Lo": (0.8332, -42.06, +1),    # left arm, outer edge
    "Li": (0.7484, 53.83, -1),     # left arm, edge facing the crotch
    "Ri": (-0.8751, 506.45, +1),   # right arm, edge facing the crotch
    "Ro": (-0.9451, 605.73, -1),   # right arm, outer edge
}
#  radius of the end of the left arm, of the right arm, of the crotch, of the
#  flare of the left and of the right arm into the stem; the stem's two edges
#  and its foot
PARAMS = dict(r_left=35.49, r_right=35.09, r_crotch=60.9, r_flare_l=50.0,
              r_flare_r=50.0, stem_l=203.08, stem_r=330.4, foot=542.73)

#  The run of thicker glaze (already shifted by DX, DY)
RUN = ("M319.8,554.8 322,555 324.1,554.8 326,555 328.1,554.8 329.5,555 332.1,553.6 335.9,551.1 "
       "339.4,548.3 341.1,546.8 344.2,543.4 345.6,541.6 348.2,537.9 349.3,535.9 351.3,531.8 "
       "352.1,529.7 353.5,525.3 354.4,520.9 354.8,516.4 354.9,513.8 354.9,343.8 355.2,340.4 "
       "355.6,336.5 356.2,333.2 356.9,329.9 357.9,326.2 359,323 360.2,319.8 361.8,316.3 "
       "364.2,311.9 366.8,307.8 369.7,303.8 372.9,300.1 553.4,109.1 554.7,107.7 555.6,106.2 "
       "556.5,104.7 557.2,103.1 557.7,101.5 558,99.8 558.2,97.2 558.1,95.5 557.9,93.8 557.4,92.1 "
       "556.8,90.5 556.1,89 555.2,87.5 554.1,86.1 552.9,84.9 551.6,83.8 549.6,82.4 545.9,85.3 "
       "537.2,93.2 534.4,95.8 531.1,99.4 524.1,107.6 508.3,127.7 497.3,140.8 487.8,151.5 "
       "459.6,183 413.7,235.5 402.1,248.4 390.8,260.4 381.4,270.2 362.1,289.6 356.7,295.2 "
       "352.2,300.2 346.1,307.5 340.5,315 336,322 333.8,326 331.7,330 329.1,335.8 326.5,342.6 "
       "324.2,349.5 322.4,356.6 320.8,363.9 319.5,372.1 318.7,379.6 318,388.1 317.7,399.5 "
       "317.9,412.1 318.5,426.7 320.3,455.1 321,467.7 321.3,480.1 321.2,491.5 320.6,502.8 "
       "319.3,513.2 317.6,522.9 313.5,540.9 312.8,545.1 312.6,548.6 312.8,550.5 313.2,552.1 "
       "313.9,553.3 314.6,554 315.8,554.7 316.9,554.9 318.1,555Z")
#  The highlight along the left arm: a stroke
SHINE = "M108.5,90 L256.5,270 Q267.5,284 278.5,287"
SHINE_W = 10.0

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
#  glaze: the gradient from the top left to the foot; run: the darker glaze,
#  from the top of the right arm down to the foot; shine: the highlight
PALETTES = {
    "celadon": dict(glaze=("#B5D7C5", "#1F5644"), run=("#618E81", "#346A69"),
                    shine=("#FFFFFF", 0.20)),
    "reversed": dict(glaze=("#C4E2D2", "#3E7A66"), run=("#6F9C8E", "#437C7A"),
                     shine=("#FFFFFF", 0.22)),
    "amber": dict(glaze=("#EED5A6", "#6E4310"), run=("#B8884A", "#7E4F1C"),
                  shine=("#FFFFFF", 0.22)),
    "lapis": dict(glaze=("#B9CBEA", "#1C3263"), run=("#6E86B5", "#2E4777"),
                  shine=("#FFFFFF", 0.20)),
    "mono": dict(glaze=("#D3D2CF", "#3A3936"), run=("#8E8C88", "#55534F"),
                 shine=("#FFFFFF", 0.22)),
}
SUMI, CELADON, SHADOW = "#22201D", "#5E8C78", "#3C5F4E"
PAPER, KILN = "#F4F1EA", "#2B211B"

# ---------------------------------------------------------------------------
# The outline
# ---------------------------------------------------------------------------
def _normal(a, s):
    h = math.hypot(1.0, a)
    return (s / h, -s * a / h)


def _offset(e, d):
    """The edge e moved by d towards the letter (d < 0: away from it)."""
    a, b, s = e
    nx, ny = _normal(a, s)
    x0, y0 = b + d * nx, d * ny
    return (a, x0 - a * y0, s)


def _meet(e1, e2):
    (a1, b1), (a2, b2) = e1[:2], e2[:2]
    y = (b2 - b1) / (a1 - a2)
    return (a1 * y + b1, y)


def _foot(e, c):
    a, b = e[:2]
    h = math.hypot(a, 1.0)
    ux, uy = a / h, 1.0 / h
    t = (c[0] - b) * ux + c[1] * uy
    return (b + t * ux, t * uy)


def outline_path():
    p = PARAMS
    E = EDGES
    stem_l, stem_r = (0.0, p["stem_l"], +1), (0.0, p["stem_r"], -1)
    rb = (p["stem_r"] - p["stem_l"]) / 2
    yb = p["foot"] - rb
    cL = _meet(_offset(E["Lo"], p["r_left"]), _offset(E["Li"], p["r_left"]))
    cR = _meet(_offset(E["Ri"], p["r_right"]), _offset(E["Ro"], p["r_right"]))
    cC = _meet(_offset(E["Li"], -p["r_crotch"]), _offset(E["Ri"], -p["r_crotch"]))
    cFL = _meet(_offset(E["Lo"], -p["r_flare_l"]), _offset(stem_l, -p["r_flare_l"]))
    cFR = _meet(_offset(E["Ro"], -p["r_flare_r"]), _offset(stem_r, -p["r_flare_r"]))
    pt = lambda q: "%s,%s" % (_f(q[0] + DX), _f(q[1] + DY))
    arc = lambda r, large, sweep, q: "A%s,%s 0 %d %d %s" % (_f(r), _f(r), large, sweep, pt(q))
    d = "M" + pt(_foot(E["Lo"], cL))
    d += " L" + pt(_foot(E["Lo"], cFL))
    d += " " + arc(p["r_flare_l"], 0, 1, _foot(stem_l, cFL))
    d += " L" + pt((p["stem_l"], yb))
    d += " " + arc(rb, 0, 0, (p["stem_r"], yb))
    d += " L" + pt(_foot(stem_r, cFR))
    d += " " + arc(p["r_flare_r"], 0, 1, _foot(E["Ro"], cFR))
    d += " L" + pt(_foot(E["Ro"], cR))
    d += " " + arc(p["r_right"], 1, 0, _foot(E["Ri"], cR))
    d += " L" + pt(_foot(E["Ri"], cC))
    d += " " + arc(p["r_crotch"], 0, 1, _foot(E["Li"], cC))
    d += " L" + pt(_foot(E["Li"], cL))
    d += " " + arc(p["r_left"], 1, 0, _foot(E["Lo"], cL))
    return d + "Z"


def _f(v):
    s = "%.2f" % v
    return s.rstrip("0").rstrip(".")


# ---------------------------------------------------------------------------
# SVG pieces
# ---------------------------------------------------------------------------
#  The extent of the letter on the canvas (for the gradient)
X0, Y0, X1, Y1 = 21.36 + DX, 27.14 + DY, 533.72 + DX, 542.71 + DY
RUN_Y0, RUN_Y1 = 60 + DY, 520 + DY


def mark_group(pal, uid, dash=None, shine_opacity=None):
    """The mark's elements and their defs. uid keeps the ids apart when
    several marks share a page; pal is a palette name, or one fill color.
    dash = (percent of the run hidden, SMIL children of the mask's stroke)
    draws the run as far as it has flowed."""
    O = outline_path()
    if pal not in PALETTES:
        return '<path d="%s" fill="%s"/>' % (O, pal), ""
    P = PALETTES[pal]
    defs = ('<linearGradient id="%s-g" gradientUnits="userSpaceOnUse" x1="%s" y1="%s" x2="%s" y2="%s">'
            '<stop offset="0" stop-color="%s"/><stop offset="1" stop-color="%s"/></linearGradient>'
            '<linearGradient id="%s-r" gradientUnits="userSpaceOnUse" x1="0" y1="%s" x2="0" y2="%s">'
            '<stop offset="0" stop-color="%s"/><stop offset="1" stop-color="%s"/></linearGradient>'
            '<clipPath id="%s-c"><path d="%s"/></clipPath>'
            % (uid, _f(X0), _f(Y0), _f(X1), _f(Y1), P["glaze"][0], P["glaze"][1],
               uid, _f(RUN_Y0), _f(RUN_Y1), P["run"][0], P["run"][1], uid, O))
    run = '<path d="%s" fill="url(#%s-r)"/>' % (RUN, uid)
    if dash is not None:
        defs += ('<mask id="%s-m"><path d="%s" fill="none" stroke="#fff" stroke-width="80" '
                 'stroke-linecap="round" stroke-linejoin="round" '
                 'stroke-dasharray="%s %s" stroke-dashoffset="%s">%s</path></mask>'
                 % (uid, RUN_SPINE, _f(SPINE_LEN), _f(SPINE_LEN), _f(dash[0] * SPINE_LEN / 100), dash[1]))
        run = '<g mask="url(#%s-m)">%s</g>' % (uid, run)
    so = P["shine"][1] if shine_opacity is None else shine_opacity
    shine = ('<path d="%s" fill="none" stroke="%s" stroke-opacity="%s" stroke-width="%s" '
             'stroke-linecap="round" stroke-linejoin="round"/>' % (SHINE, P["shine"][0], _f(so), _f(SHINE_W)))
    body = ('<path d="%s" fill="url(#%s-g)"/><g clip-path="url(#%s-c)">%s%s</g>'
            % (O, uid, uid, run, shine))
    return body, defs


#  The line the run follows, for the animation (shifted)
RUN_SPINE = "M578,70 L549,93 L528,120 L492,160 L455,200 L418,240 L381,280 L349,320 L339,360 L336,400 L337,480 L334,520 L318,575"


def _spine_length():
    import re
    p = [tuple(map(float, q)) for q in re.findall(r"([\d.]+),([\d.]+)", RUN_SPINE)]
    return sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(p, p[1:]))


SPINE_LEN = _spine_length()


def flat_group(base, run):
    O = outline_path()
    return ('<defs><clipPath id="f-c"><path d="%s"/></clipPath></defs><path d="%s" fill="%s"/>'
            '<g clip-path="url(#f-c)"><path d="%s" fill="%s"/></g>' % (O, O, base, RUN, run))


def svg(w, h, title, body, defs="", bg=None, px=None):
    pw, ph = px if px else (w / 2.0, h / 2.0)
    out = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %s %s" width="%s" height="%s">'
           '<title>%s</title>' % (_f(w), _f(h), _f(pw), _f(ph), title))
    if defs:
        out += "<defs>%s</defs>" % defs
    if bg:
        out += '<rect width="%s" height="%s" fill="%s"/>' % (_f(w), _f(h), bg)
    return out + body + "</svg>\n"


# ---------------------------------------------------------------------------
# The logotype: Bitstream Charter Bold, as outlines (font units, 1000/em)
# ---------------------------------------------------------------------------
GLYPHS = [
    (0, "M-11 671V622L53 612L245 264V69L229 55L162 50V0H462V50L394 55L378 69V288L568 611L627 622V671H409V623L469 612L477 596L334 351L197 605L206 617L262 622V671Z"),
    (610.0, "M164 499 23 451V412H88V162C88 47 132 -9 225 -9C278 -9 330 14 393 55V0H582V45L530 49L516 64V499H470L328 451V412H393V93C360 83 332 77 307 77C238 77 212 111 212 189V499Z"),
    (1210.0, "M46 0H465V164H397L385 72L370 56H185L450 434V486H54V335H119L129 418L142 432H304L46 62Z"),
    (1705.0, "M491 124C433 87 379 69 326 69C227 69 165 137 165 252C165 253 165 254 165 255H492C492 261 492 268 492 274C492 418 415 500 284 500C140 500 37 395 37 242C37 91 135 -9 284 -9C354 -9 424 11 491 50ZM361 308H170C173 395 214 445 274 445C329 445 361 402 361 319C361 315 361 312 361 308Z"),
    (2229.0, "M31 0H278V45L234 49L220 64V394C258 405 290 411 315 411C378 411 405 377 405 303V64L391 49L348 45V0H593V45L540 49L527 64V332C527 443 483 498 391 498C347 498 291 476 220 437V496H172L28 449V409H96V64L82 49L31 45Z"),
]
TEXT_LEFT, TEXT_RIGHT, CAP = -11.0, 2229.0 + 593.0, 671.0
#  Cap height of the logotype against the height of the Y, and the gap
#  between them against that height (the proportions of the first lockup)
CAP_RATIO, GAP_RATIO, BASE_RATIO = 0.429, 0.183, 0.675


def logotype(x_left, baseline, scale, color):
    tx = x_left - TEXT_LEFT * scale
    return "".join('<path transform="translate(%s %s) scale(%s %s)" d="%s" fill="%s"/>'
                   % (_f(tx + adv * scale), _f(baseline), _f(scale), _f(-scale), d, color)
                   for adv, d in GLYPHS)


def lockup_h(pal, text, bg=None, uid="lh"):
    yh = Y1 - Y0
    s = CAP_RATIO * yh / CAP
    x_text = X1 + GAP_RATIO * yh
    base = Y0 + BASE_RATIO * yh
    w = x_text + (TEXT_RIGHT - TEXT_LEFT) * s + X0
    body, defs = mark_group(pal, uid)
    return svg(w, SIZE, "Yuzen — lockup", body + logotype(x_text, base, s, text), defs, bg,
               px=(w / SIZE * 140, 140))


def lockup_v(pal, text, bg=None, uid="lv"):
    yh = Y1 - Y0
    s = CAP_RATIO * yh / CAP
    tw = (TEXT_RIGHT - TEXT_LEFT) * s
    w = max(SIZE, tw + 2 * X0)
    shift = (w - SIZE) / 2
    base = Y1 + GAP_RATIO * yh + CAP * s
    h = base + Y0
    body, defs = mark_group(pal, uid)
    body = '<g transform="translate(%s 0)">%s</g>' % (_f(shift), body)
    return svg(w, h, "Yuzen — lockup", body + logotype((w - tw) / 2, base, s, text), defs, bg,
               px=(w / h * 300, 300))


# ---------------------------------------------------------------------------
def write(name, text):
    with open(os.path.join(HERE, name), "w", encoding="utf-8") as f:
        f.write(text)
    return os.path.join(HERE, name)


ANIM_DUR = 4.0


def build_svgs():
    W = SIZE
    out = {}
    b, d = mark_group("celadon", "ym")
    out["yuzen-mark.svg"] = svg(W, W, "Yuzen — mark (celadon)", b, d)
    b, d = mark_group("reversed", "yr")
    out["yuzen-mark-reversed.svg"] = svg(W, W, "Yuzen — reversed (dark ground)", b, d, bg=KILN)
    b, d = mark_group("amber", "ya")
    out["yuzen-mark-amber.svg"] = svg(W, W, "Yuzen — amber", b, d)
    b, d = mark_group("lapis", "yl")
    out["yuzen-mark-lapis.svg"] = svg(W, W, "Yuzen — lapis", b, d)
    b, d = mark_group("mono", "yg")
    out["yuzen-mark-mono.svg"] = svg(W, W, "Yuzen — monochrome", b, d)
    out["yuzen-mark-flat.svg"] = svg(W, W, "Yuzen — flat (two tones, no gradient)",
                                     flat_group(CELADON, SHADOW))
    b, _ = mark_group(CELADON, "ys")
    out["yuzen-mark-solid.svg"] = svg(W, W, "Yuzen — solid (small sizes)", b)
    b, _ = mark_group(SUMI, "y1")
    out["yuzen-mark-1color.svg"] = svg(W, W, "Yuzen — single color (print and engraving)", b)
    L = _f(SPINE_LEN)
    anim = ('<animate attributeName="stroke-dashoffset" values="%s;%s;0;0" '
            'keyTimes="0;0.15;0.55;1" dur="%gs" repeatCount="indefinite"/>' % (L, L, ANIM_DUR))
    b, d = mark_group("celadon", "yn", dash=(0, anim))
    out["yuzen-mark-animated.svg"] = svg(W, W, "Yuzen — animated", b, d)
    out["yuzen-lockup-h.svg"] = lockup_h("celadon", SUMI)
    out["yuzen-lockup-h-reversed.svg"] = lockup_h("reversed", PAPER, bg=KILN, uid="lhr")
    out["yuzen-lockup-v.svg"] = lockup_v("celadon", SUMI)
    for k, v in out.items():
        write(k, v)
    return out


def build_spec():
    """yuzen-logo-spec.svg: the mark, its clear space, the colors and the
    variants on one sheet."""
    W, H = 1600, 1040
    defs, body = [], []
    font = 'font-family="DejaVu Sans, Helvetica, Arial, sans-serif"'
    t = lambda x, y, txt, size=22, fill=SUMI, extra="": body.append(
        '<text x="%s" y="%s" %s font-size="%s" fill="%s" %s>%s</text>'
        % (x, y, "" if "font-family" in extra else font, size, fill, extra, txt))

    def place(pal, uid, x, y, px, bg=None):
        k = px / SIZE
        if bg:
            body.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s"/>' % (x, y, px, px, bg))
        if pal == "flat":
            g = flat_group(CELADON, SHADOW).replace('id="f-c"', 'id="%s-fc"' % uid).replace(
                "url(#f-c)", "url(#%s-fc)" % uid)
            body.append('<g transform="translate(%s %s) scale(%s)">%s</g>' % (x, y, _f(k), g))
            return
        b, d = mark_group(pal, uid)
        defs.append(d)
        body.append('<g transform="translate(%s %s) scale(%s)">%s</g>' % (x, y, _f(k), b))

    body.append('<rect width="%d" height="%d" fill="%s"/>' % (W, H, PAPER))
    t(60, 70, "Yuzen — logo specification", 34, extra='font-weight="bold"')
    t(60, 104, "A glazed Y: celadon lit from the upper left, a highlight on the left arm, "
               "a darker run of thicker glaze down the right arm and stem.", 19, "#55504A")
    # the mark with its clear space
    x0, y0, px = 60, 140, 520
    k = px / SIZE
    body.append('<rect x="%s" y="%s" width="%s" height="%s" fill="#FFFFFF" stroke="#C9C3B8" '
                'stroke-dasharray="6 5"/>' % (x0, y0, px, px))
    bx, by = x0 + X0 * k, y0 + Y0 * k
    body.append('<rect x="%s" y="%s" width="%s" height="%s" fill="none" stroke="#C08A3E" '
                'stroke-width="1.5"/>' % (_f(bx), _f(by), _f((X1 - X0) * k), _f((Y1 - Y0) * k)))
    place("celadon", "s0", x0, y0, px)
    t(x0, y0 + px + 34, "Clear space: the dashed square, built into every SVG viewBox", 18, "#55504A")
    t(x0, y0 + px + 60, "(about one arm width on every side of the letter).", 18, "#55504A")
    t(x0, y0 + px + 96, "Minimum size: 32 px for the full mark;", 18, "#55504A")
    t(x0, y0 + px + 122, "below that, the solid mark.", 18, "#55504A")
    t(x0, y0 + px + 158, "The highlight sits on the left arm and the run on", 18, "#55504A")
    t(x0, y0 + px + 184, "the right: never mirror the mark.", 18, "#55504A")
    # colors
    cx, cy = 660, 160
    t(cx, cy, "Colors", 24, extra='font-weight="bold"')
    swatches = [
        ("Glaze, light", PALETTES["celadon"]["glaze"][0]), ("Glaze, deep", PALETTES["celadon"]["glaze"][1]),
        ("Run, top", PALETTES["celadon"]["run"][0]), ("Run, foot", PALETTES["celadon"]["run"][1]),
        ("Celadon (solid)", CELADON), ("Shadow (flat run)", SHADOW),
        ("Sumi", SUMI), ("Paper", PAPER), ("Kiln", KILN)]
    for i, (name, col) in enumerate(swatches):
        sx, sy = cx + (i % 3) * 300, cy + 24 + (i // 3) * 92
        body.append('<rect x="%d" y="%d" width="64" height="64" rx="6" fill="%s" stroke="#C9C3B8"/>' % (sx, sy, col))
        t(sx + 78, sy + 28, name, 18)
        t(sx + 78, sy + 54, col, 18, "#55504A", 'font-family="DejaVu Sans Mono, monospace"')
    t(cx, cy + 330, "Glaze: linear, from the top left of the letter to the foot of the stem.", 17, "#55504A")
    t(cx, cy + 354, "Run: vertical, top to foot. Highlight: white at 20 % opacity.", 17, "#55504A")
    # variants
    vy = 560
    t(cx, vy, "Variants", 24, extra='font-weight="bold"')
    names = [("celadon", "mark"), ("reversed", "reversed"), ("amber", "amber"), ("lapis", "lapis"),
             ("mono", "mono"), ("flat", "flat"), (CELADON, "solid"), (SUMI, "1color")]
    for i, (pal, label) in enumerate(names):
        vx, vyy = cx + (i % 4) * 225, vy + 20 + (i // 4) * 215
        place(pal, "v%d" % i, vx, vyy, 170, bg=KILN if pal == "reversed" else "#FFFFFF")
        t(vx, vyy + 194, "yuzen-mark%s.svg" % ("" if label == "mark" else "-" + label), 15, "#55504A",
          'font-family="DejaVu Sans Mono, monospace"')
    out = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d">'
           '<title>Yuzen — logo specification</title><defs>%s</defs>%s</svg>\n'
           % (W, H, W, H, "".join(defs), "".join(body)))
    write("yuzen-logo-spec.svg", out)


def anim_frame(t):
    """The run's mask offset at time t (s), as the SMIL animation has it."""
    u = (t % ANIM_DUR) / ANIM_DUR
    if u < 0.15:
        return 100.0
    if u < 0.55:
        x = (u - 0.15) / 0.40
        return 100.0 * (1 - x)
    return 0.0


def build_rasters():
    import io
    import cairosvg
    from PIL import Image
    png = os.path.join(HERE, "png")
    os.makedirs(png, exist_ok=True)

    def render(name, px_w, out=None, data=None):
        if data:
            img = cairosvg.svg2png(bytestring=data, output_width=px_w)
        else:
            img = cairosvg.svg2png(url=os.path.join(HERE, name), output_width=px_w)
        if out:
            open(os.path.join(png, out), "wb").write(img)
        return Image.open(io.BytesIO(img)).convert("RGBA")

    for n in (64, 128, 256, 512, 1024):
        render("yuzen-mark.svg", n, "yuzen-mark-%d.png" % n)
        render("yuzen-mark-reversed.svg", n, "yuzen-mark-reversed-%d.png" % n)
    # icons: the clear space cropped to a thin margin, the solid mark below 32px
    icons = {}
    side = max(X1 - X0, Y1 - Y0) + 2 * 16
    vb = "%s %s %s %s" % (_f((X0 + X1 - side) / 2), _f((Y0 + Y1 - side) / 2), _f(side), _f(side))
    for n in (16, 24, 32, 48):
        b, d = mark_group(CELADON if n < 32 else "celadon", "ic")
        data = svg(SIZE, SIZE, "icon", b, d).replace('viewBox="0 0 640 640"', 'viewBox="%s"' % vb)
        icons[n] = render(None, n, "yuzen-icon-%d.png" % n, data=data.encode())
    for n in (400, 800, 1600):
        render("yuzen-lockup-h.svg", n, "yuzen-lockup-h-%d.png" % n)
    render("yuzen-lockup-v.svg", 800, "yuzen-lockup-v-800.png")
    icons[48].save(os.path.join(HERE, "favicon.ico"), sizes=[(16, 16), (32, 32), (48, 48)],
                   append_images=[icons[16], icons[32]])
    # the GIF: the run flowing down, on paper. cairosvg ignores masks, so the
    # frames are drawn by librsvg (rsvg-convert).
    import shutil
    import subprocess
    if not shutil.which("rsvg-convert"):
        raise SystemExit("make_logo.py: rsvg-convert (librsvg) is needed for the GIF")
    frames = []
    fps = 15
    for i in range(int(ANIM_DUR * fps)):
        b, d = mark_group("celadon", "gf", dash=(anim_frame(i / fps), ""))
        data = svg(SIZE, SIZE, "frame", b, d, bg=PAPER).encode()
        im = Image.open(io.BytesIO(subprocess.run(["rsvg-convert", "-w", "256"], input=data,
                                                  capture_output=True, check=True).stdout))
        frames.append(im.convert("RGB").quantize(colors=128, method=Image.Quantize.MEDIANCUT))
    frames[0].save(os.path.join(HERE, "yuzen-mark-animated.gif"), save_all=True,
                   append_images=frames[1:], duration=int(1000 / fps), loop=0, optimize=True)
    cairosvg.svg2png(url=os.path.join(HERE, "yuzen-logo-spec.svg"),
                     write_to=os.path.join(HERE, "yuzen-logo-spec.png"), output_width=1600)
    # the manual's title page and running heads
    os.makedirs(MANUAL_FIGS, exist_ok=True)
    for name in ("yuzen-mark", "yuzen-lockup-h"):
        cairosvg.svg2pdf(url=os.path.join(HERE, name + ".svg"),
                         write_to=os.path.join(MANUAL_FIGS, name + ".pdf"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raster", action="store_true",
                    help="also write png/, favicon.ico, the GIF and the manual's PDFs")
    a = ap.parse_args()
    build_svgs()
    build_spec()
    if a.raster:
        build_rasters()


if __name__ == "__main__":
    main()
