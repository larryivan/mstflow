#!/usr/bin/env python3
"""MSTFlow workflow as a metro map, in light and dark SVG.

Ancient and modern libraries enter on their own lines, are trimmed on their own
branch, then run side by side. After host removal the pair splits, like a line with
two branches through the centre of town: Kraken2 on one, KrakenUniq on the other,
and rejoins at Bracken, which follows either, then the filters and the count table. The source attribution
line leaves the table, fans out into one track per method and gathers again for the
comparison. A thin QC line runs underneath.

Conventions follow transit maps: ticks for single-line stops, capsules for stops
shared by both lines, small rings on the method tracks; labels above, roles below;
line names in their line's colour. Where the pair splits and rejoins, its lines cross
under a station. Background is transparent so the file sits on any GitHub theme.
Edit this file, not the SVGs, and re-run it.
"""
from math import hypot
from pathlib import Path

W, H = 1012, 268
LW = 6                                   # line width
GAP = 9                                  # centre-to-centre spacing of the two DNA lines
A_TOP, M_BOT = 52.5, 168.5               # the two entry lines
MID = 110.5                              # centre of the pair once it runs together
A_Y, M_Y = MID - GAP / 2, MID + GAP / 2
X0, XT, XB = 118, 168, 220               # line start, trim stop, bend start

X_BWA, X_BRACKEN, X_FILTERS, X_TABLE = 352, 532, 648, 730
SPLIT = 36                               # how far each branch moves off the centre line
UPPER = [("Kraken2", "classification", 442)]
LOWER = [("KrakenUniq", "or instead", 442)]

X_FAN, X_GATHER, X_END = 754, 936, 966   # the method tracks spread and gather here
METHODS = ["SourceTracker2", "FEAST", "FastST", "SourceID-NMF", "SourcePredict"]
FAN_STEP = 23

Y_QC = 238
QC_STOPS = [("FastQC", 146), ("FastQC", X_BWA), ("sample summary", 614), ("MultiQC", X_TABLE)]

THEME = {
    "light": dict(ink="#1f2328", muted="#6a737d", st_fill="#ffffff", st_edge="#1f2328",
                  ancient="#c4552f", modern="#2367a6", attr="#7d5bbe", qc="#c3cad2"),
    "dark": dict(ink="#e6edf3", muted="#8d96a0", st_fill="#e6edf3", st_edge="#0d1117",
                 ancient="#ec8a62", modern="#6aaee6", attr="#ad91e6", qc="#3d444d"),
}
FONT = "Inter,ui-sans-serif,-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif"


def rounded(pts, radii=22):
    """A polyline with its corners rounded; `radii` is one radius or one per corner."""
    radii = radii if isinstance(radii, (list, tuple)) else [radii] * (len(pts) - 2)
    d = [f"M{pts[0][0]:.2f} {pts[0][1]:.2f}"]
    for (x0, y0), (x1, y1), (x2, y2), r in zip(pts, pts[1:], pts[2:], radii):
        a, b = hypot(x1 - x0, y1 - y0), hypot(x2 - x1, y2 - y1)
        r = min(r, a / 2, b / 2)
        p = (x1 - (x1 - x0) * r / a, y1 - (y1 - y0) * r / a)
        q = (x1 + (x2 - x1) * r / b, y1 + (y2 - y1) * r / b)
        d.append(f"L{p[0]:.2f} {p[1]:.2f} Q{x1:.2f} {y1:.2f} {q[0]:.2f} {q[1]:.2f}")
    d.append(f"L{pts[-1][0]:.2f} {pts[-1][1]:.2f}")
    return " ".join(d)


def offset(pts, dist):
    """The polyline `dist` above (left of) `pts`, mitred, so paired lines stay
    parallel through 45 degree bends."""
    def normal(p, q):
        dx, dy = q[0] - p[0], q[1] - p[1]
        n = hypot(dx, dy)
        return dy / n, -dx / n
    out = []
    for i, p in enumerate(pts):
        normals = [normal(pts[j], pts[j + 1]) for j in (i - 1, i) if 0 <= j < len(pts) - 1]
        nx, ny = sum(n[0] for n in normals), sum(n[1] for n in normals)
        k = hypot(nx, ny)
        nx, ny = nx / k, ny / k
        cos = nx * normals[0][0] + ny * normals[0][1]
        out.append((p[0] + nx * dist / cos, p[1] + ny * dist / cos))
    return out


def path(d, colour, width=LW):
    return (f'<path d="{d}" stroke="{colour}" stroke-width="{width}" stroke-linecap="round" '
            f'stroke-linejoin="round" fill="none"/>')


def text(x, y, s, t, size=13, weight=500, fill=None, anchor="middle"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-size="{size}" '
            f'font-weight="{weight}" letter-spacing="-0.1" fill="{fill or t["ink"]}">{s}</text>')


def role(x, y, s, t, anchor="middle"):
    return text(x, y, s, t, size=11, weight=400, fill=t["muted"], anchor=anchor)


def capsule(x, centre, t):
    """A stop shared by both DNA lines, centred on `centre`."""
    h = GAP + LW + 9
    return (f'<rect x="{x - 8}" y="{centre - h / 2:.2f}" width="16" height="{h}" rx="8" '
            f'fill="{t["st_fill"]}" stroke="{t["st_edge"]}" stroke-width="2.4"/>')


def tick(x, y, colour, up=True, length=8, width=4.5, line=LW):
    y0 = y - line / 2 + 0.5 if up else y + line / 2 - 0.5
    y1 = y0 - length if up else y0 + length
    return f'<path d="M{x} {y0} V{y1}" stroke="{colour}" stroke-width="{width}" fill="none"/>'


def ring(x, y, colour, t, r=4.2, width=2.2):
    return (f'<circle cx="{x}" cy="{y}" r="{r}" fill="{t["st_fill"]}" stroke="{colour}" '
            f'stroke-width="{width}"/>')


def branch(direction):
    """Centre line of one branch between host removal and Bracken."""
    y = MID + direction * SPLIT
    return [(X_BWA - 12, MID), (X_BWA - 2, MID), (X_BWA - 2 + SPLIT, y),
            (X_BRACKEN + 2 - SPLIT, y), (X_BRACKEN + 2, MID), (X_BRACKEN + 12, MID)]


def pair(centre, t, radii):
    """The ancient and modern lines, GAP apart, along `centre`."""
    return [path(rounded(offset(centre, side * GAP / 2), [r - side * GAP / 2 if r > 8 else r for r in radii]),
                 t[key]) for key, side in (("ancient", 1), ("modern", -1))]


def build(name):
    t = THEME[name]
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
         f'height="{H}" font-family="{FONT}">']

    # QC runs underneath everything
    s.append(path(f"M{X0} {Y_QC} H{QC_STOPS[-1][1]}", t["qc"], width=2.5))
    s.append(text(X0 - 14, Y_QC + 4.5, "QC", t, size=12, weight=600, fill=t["muted"], anchor="end"))
    for label, x in QC_STOPS:
        s += [tick(x, Y_QC, t["qc"], up=False, length=6, width=2.5, line=2.5), role(x, Y_QC + 22, label, t)]

    # source attribution: one track per method, gathered again for the comparison
    x_label = (X_FAN + X_GATHER) / 2
    for i, method in enumerate(METHODS):
        dy = (i - (len(METHODS) - 1) / 2) * FAN_STEP
        y, lift = MID + dy, abs(dy)
        if lift == 0:
            s.append(path(f"M{X_FAN} {MID} H{X_GATHER}", t["attr"], width=2.6))
        else:
            track = [(X_FAN - 10, MID), (X_FAN, MID), (X_FAN + lift, y), (X_GATHER - lift, y), (X_GATHER, MID),
                     (X_GATHER + 10, MID)]
            s.append(path(rounded(track, [8, 14, 14, 8]), t["attr"], width=2.6))
        s.append(text(x_label, y - 6, method, t, size=11, weight=500))
    s.append(path(f"M{X_TABLE} {MID} H{X_FAN}", t["attr"]))
    s.append(path(f"M{X_GATHER} {MID} H{X_END}", t["attr"]))
    bottom = MID + (len(METHODS) - 1) / 2 * FAN_STEP
    s.append(text(x_label, bottom + 24, "Source attribution", t, size=12, weight=600, fill=t["attr"]))
    s.append(role(x_label, bottom + 38, "optional", t))

    # the DNA lines: entry, trim, host removal, then the two branches, then one table
    entry = {"ancient": [(X0, A_TOP), (XB, A_TOP), (XB + A_Y - A_TOP, A_Y), (X_BWA, A_Y)],
             "modern": [(X0, M_BOT), (XB, M_BOT), (XB + M_BOT - M_Y, M_Y), (X_BWA, M_Y)]}
    for key, pts in entry.items():
        s.append(path(rounded(pts), t[key]))
    for direction in (-1, 1):
        s += pair(branch(direction), t, [6, 26, 26, 6])
    for key, y in (("ancient", A_Y), ("modern", M_Y)):
        s.append(path(f"M{X_BRACKEN} {y} H{X_TABLE}", t[key]))

    s.append(text(X0 - 14, A_TOP + 4.5, "Ancient DNA", t, size=12, weight=600, fill=t["ancient"], anchor="end"))
    s.append(text(X0 - 14, M_BOT + 4.5, "Modern DNA", t, size=12, weight=600, fill=t["modern"], anchor="end"))
    for y, key, tool, what, up in ((A_TOP, "ancient", "AdapterRemoval", "trim, merge pairs", True),
                                   (M_BOT, "modern", "fastp", "trim", True)):
        s += [tick(XT, y, t[key]), text(XT, y - 17, tool, t), role(XT, y + LW / 2 + 14, what, t)]

    above, below = -(GAP / 2 + LW / 2 + 4.5 + 10), GAP / 2 + LW / 2 + 4.5 + 15
    # where the branches split and rejoin, names sit on the side with no branch
    # (the branch above leaves at 45 degrees, so the name can sit closer than the role)
    for tool, what, x, anchor, side in (("BWA", "host removal", X_BWA, "end", -1),
                                        ("Bracken", "abundance", X_BRACKEN, "start", 1)):
        s += [capsule(x, MID, t), text(x - side * 8, MID + above, tool, t, anchor=anchor),
              role(x + side * 12, MID + below, what, t, anchor=anchor)]
    s += [capsule(X_FILTERS, MID, t), text(X_FILTERS, MID + above, "Filters", t),
          role(X_FILTERS, MID + below, "k-mers, taxa", t)]
    for stops, direction in ((UPPER, -1), (LOWER, 1)):
        y = MID + direction * SPLIT
        for tool, what, x in stops:
            s += [capsule(x, y, t), text(x, y + above, tool, t) if direction < 0 else text(x, y + below, tool, t),
                  role(x, y + below, what, t) if direction < 0 else role(x, y + below + 14, what, t)]
    s += [capsule(X_TABLE, MID, t), text(X_TABLE, MID + above, "Count table", t, weight=650)]

    # the end of the attribution line
    s += [f'<rect x="{X_END - 8}" y="{MID - 8}" width="16" height="16" rx="8" fill="{t["st_fill"]}" '
          f'stroke="{t["attr"]}" stroke-width="2.4"/>',
          text(X_END, MID + above, "Comparison", t, weight=650),
          role(X_END, MID + below - 4, "cross-validated", t),
          role(X_END, MID + below + 10, "by study", t)]

    s.append("</svg>")
    return "\n".join(s) + "\n"


here = Path(__file__).parent
for name in THEME:
    (here / f"mstflow-workflow-{name}.svg").write_text(build(name))
    print("wrote", here / f"mstflow-workflow-{name}.svg")
