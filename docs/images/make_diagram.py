#!/usr/bin/env python3
"""MSTFlow workflow as a metro map, in light and dark SVG.

Ancient and modern libraries enter on their own lines, are trimmed on their own
branch, then run side by side through the shared trunk. The two lines never merge:
each keeps its own fragment semantics through classification.

Conventions follow transit maps: ticks for single-line stops, capsules for stops
shared by both lines, labels above, roles below. Background is transparent so the
file sits on any GitHub theme. Edit this file, not the SVGs, and re-run it.
"""
from math import hypot
from pathlib import Path

W, H = 870, 213
LW = 6                                   # line width
A_TOP, M_BOT = 59.5, 160.5               # the two entry lines
A_Y, M_Y = 105.5, 114.5                  # the same lines on the trunk
X0, XT, XB = 122, 196, 246               # line start, trim stop, bend start
BEND = A_Y - A_TOP                       # 45 degree bend: dx equals dy
XM = XB + BEND                           # where the lines reach the trunk
STOPS = [("BWA", "host removal", 384), ("Kraken2", "classification", 482),
         ("Bracken", "abundance", 580)]
XEND, XOPT = 678, 794

THEME = {
    "light": dict(ink="#1f2328", muted="#6a737d", st_fill="#ffffff", st_edge="#1f2328",
                  ancient="#c4552f", modern="#2367a6", opt="#b4bcc4"),
    "dark": dict(ink="#e6edf3", muted="#8d96a0", st_fill="#e6edf3", st_edge="#0d1117",
                 ancient="#ec8a62", modern="#6aaee6", opt="#59626c"),
}
FONT = "Inter,ui-sans-serif,-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif"


def rounded(pts, r=22):
    d = [f"M{pts[0][0]} {pts[0][1]}"]
    for (x0, y0), (x1, y1), (x2, y2) in zip(pts, pts[1:], pts[2:]):
        a, b = hypot(x1 - x0, y1 - y0), hypot(x2 - x1, y2 - y1)
        p = (x1 - (x1 - x0) * r / a, y1 - (y1 - y0) * r / a)
        q = (x1 + (x2 - x1) * r / b, y1 + (y2 - y1) * r / b)
        d.append(f"L{p[0]:.2f} {p[1]:.2f} Q{x1} {y1} {q[0]:.2f} {q[1]:.2f}")
    d.append(f"L{pts[-1][0]} {pts[-1][1]}")
    return " ".join(d)


def label(x, y, s, t, size=13, weight=500, fill=None, anchor="middle"):
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
            f'font-weight="{weight}" letter-spacing="-0.1" fill="{fill or t["ink"]}">{s}</text>')


def capsule(x, t):
    top = A_Y - LW / 2 - 4.5
    h = (M_Y - A_Y) + LW + 9
    return (f'<rect x="{x - 8}" y="{top}" width="16" height="{h}" rx="8" '
            f'fill="{t["st_fill"]}" stroke="{t["st_edge"]}" stroke-width="2.4"/>')


def tick(x, y, colour):
    return (f'<path d="M{x} {y - LW / 2 + 0.5} V{y - LW / 2 - 8}" stroke="{colour}" '
            f'stroke-width="4.5" fill="none"/>')


def build(name):
    t = THEME[name]
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
         f'height="{H}" font-family="{FONT}">']
    mid = (A_Y + M_Y) / 2
    # optional spur first, so the terminus capsule covers its start
    s.append(f'<path d="M{XEND} {mid} H{XOPT}" stroke="{t["opt"]}" stroke-width="{LW - 1}" '
             f'stroke-dasharray="0.1 9" stroke-linecap="round" fill="none"/>')
    lines = {"ancient": [(X0, A_TOP), (XB, A_TOP), (XM, A_Y), (XEND, A_Y)],
             "modern": [(X0, M_BOT), (XB, M_BOT), (XM, M_Y), (XEND, M_Y)]}
    for key, pts in lines.items():
        s.append(f'<path d="{rounded(pts)}" stroke="{t[key]}" stroke-width="{LW}" '
                 f'stroke-linecap="round" stroke-linejoin="round" fill="none"/>')
    # line names at the origins
    s.append(label(X0 - 16, A_TOP + 4.5, "Ancient DNA", t, size=12, weight=600,
                   fill=t["ancient"], anchor="end"))
    s.append(label(X0 - 16, M_BOT + 4.5, "Modern DNA", t, size=12, weight=600,
                   fill=t["modern"], anchor="end"))
    # trimming, one stop on each line
    for y, key, tool, role in ((A_TOP, "ancient", "AdapterRemoval", "trim and merge"),
                               (M_BOT, "modern", "fastp", "trim")):
        s += [tick(XT, y, t[key]), label(XT, y - 17, tool, t),
              label(XT, y + LW / 2 + 14, role, t, size=11, weight=400, fill=t["muted"])]
    # shared trunk
    below = M_Y + LW / 2 + 4.5 + 15
    for tool, role, x in STOPS:
        s += [capsule(x, t), label(x, A_Y - LW / 2 - 4.5 - 10, tool, t),
              label(x, below, role, t, size=11, weight=400, fill=t["muted"])]
    s += [capsule(XEND, t), label(XEND, A_Y - LW / 2 - 4.5 - 10, "Genus table", t, weight=650)]
    # optional destination
    s += [f'<circle cx="{XOPT}" cy="{mid}" r="6.5" fill="{t["st_fill"]}" stroke="{t["opt"]}" '
          f'stroke-width="2.4"/>',
          label(XOPT, A_Y - LW / 2 - 4.5 - 10, "Source attribution", t, fill=t["muted"]),
          label(XOPT, below, "optional, 5 methods", t, size=11, weight=400, fill=t["muted"])]
    s.append("</svg>")
    return "\n".join(s) + "\n"


here = Path(__file__).parent
for name in THEME:
    (here / f"mstflow-workflow-{name}.svg").write_text(build(name))
    print("wrote", here / f"mstflow-workflow-{name}.svg")
