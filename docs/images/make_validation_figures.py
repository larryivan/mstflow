#!/usr/bin/env python3
"""Figures for docs/validation.md, from the tables in docs/validation/.

Each figure is written as a light and a dark SVG for the documentation, and as a
light PDF for print. Environments keep one colour in every figure (oral, gut, soil,
laboratory, in the validated categorical order; the unknown source in neutral grey);
where a figure compares conditions rather than environments, marker shape carries
the condition. Needs matplotlib. Edit this file, not the figures, and re-run it.
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "validation"

THEME = {
    "light": dict(surface="#ffffff", ink="#0b0b0b", secondary="#52514e", muted="#898781", grid="#e1e0d9",
                  axis="#c3c2b7", unknown="#b9b8b0",
                  env={"oral": "#2a78d6", "gut": "#eb6834", "soil": "#1baf7a", "laboratory": "#eda100"}),
    "dark": dict(surface="#0d1117", ink="#f0efec", secondary="#c3c2b7", muted="#898781", grid="#2c2c2a",
                 axis="#383835", unknown="#5d5c57",
                 env={"oral": "#3987e5", "gut": "#d95926", "soil": "#199e70", "laboratory": "#c98500"}),
}
ENVS = ["oral", "gut", "soil", "laboratory", "Unknown"]
METHOD_NAMES = {"expected": "Known composition", "sourcetracker2": "SourceTracker2", "feast": "FEAST",
                "fastst": "FastST", "sourceid_nmf": "SourceID-NMF"}


def read(name):
    with open(DATA / name, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    return None if value in ("NA", "") else float(value)


def style(t):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 7,
        "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0,
        "xtick.major.size": 2.5, "ytick.major.size": 0, "svg.fonttype": "path", "svg.hashsalt": "mstflow",
        "pdf.fonttype": 42,
        "figure.facecolor": t["surface"], "axes.facecolor": t["surface"], "savefig.facecolor": t["surface"],
        "text.color": t["ink"], "axes.labelcolor": t["secondary"], "xtick.color": t["muted"],
        "ytick.color": t["ink"], "axes.edgecolor": t["axis"], "axes.titlecolor": t["ink"],
    })


def clean(ax, t, xgrid=True):
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(t["axis"])
    if xgrid:
        ax.grid(axis="x", color=t["grid"], linewidth=0.5)
        ax.set_axisbelow(True)


def letter(ax, text):
    ax.text(-0.02, 1.08, text, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")


def env_colour(t, env):
    return t["unknown"] if env == "Unknown" else t["env"][env]


def mixtures(t):
    rows = read("mixtures.tsv")
    panels = [("mix_oral_gut", "Plaque + faeces, 1 : 1"), ("mix_oral_gut_ancient", "The same, made ancient-like"),
              ("mix_gut_soil", "Faeces + meadow soil, 1 : 9")]
    order = ["expected", "sourcetracker2", "feast", "fastst", "sourceid_nmf"]
    share = {(r["mixture"], r["row"], r["environment"]): float(r["share"]) for r in rows}
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.3), sharey=True)
    for index, (ax, (mixture, title)) in enumerate(zip(axes, panels)):
        for y, row in enumerate(order):
            left = 0.0
            for env in ENVS:
                width = share.get((mixture, row, env), 0.0)
                if width > 0:
                    ax.barh(y, width, left=left, height=0.62, color=env_colour(t, env), edgecolor=t["surface"],
                            linewidth=1.2)
                left += width
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xticklabels(["0", "", "0.5", "", "1"])
        ax.set_title(title, loc="left", pad=6)
        clean(ax, t, xgrid=False)
        letter(ax, "abc"[index])
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels([METHOD_NAMES[row] for row in order])
    axes[0].invert_yaxis()
    axes[0].get_yticklabels()[0].set_fontweight("bold")
    axes[1].set_xlabel("Share of the sink's counts")
    fig.legend(handles=[Patch(color=env_colour(t, env), label=env if env != "Unknown" else "unknown") for env in ENVS],
               loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.55, 1.02), handlelength=1.0,
               columnspacing=1.4)
    fig.subplots_adjust(left=0.14, right=0.99, bottom=0.2, top=0.8, wspace=0.12)
    return fig


def ancient_like(t):
    rows = read("ancient_like.tsv")
    names = [r["library"] for r in rows]
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.2, 2.5), sharey=True)
    for y, r in enumerate(rows):
        colour = t["env"][r["environment"]]
        cut, between = float(r["bray_curtis_cut_damaged"]), float(r["bray_curtis_same_environment"])
        left.plot([cut, between], [y, y], color=t["axis"], linewidth=0.8, zorder=1)
        left.scatter(cut, y, s=22, color=colour, edgecolor=t["surface"], linewidth=0.8, zorder=3)
        left.scatter(between, y, s=22, facecolor=t["surface"], edgecolor=colour, linewidth=1.1, zorder=3)
        values = [number(r["own_share_modern"]), number(r["own_share_cut"]), number(r["own_share_cut_damaged"])]
        present = [v for v in values if v is not None]
        right.plot([min(present), max(present)], [y, y], color=t["axis"], linewidth=0.8, zorder=1)
        # The three conditions sit a little apart within the row, so equal values stay visible.
        for value, marker, filled, offset in zip(values, ("o", "o", "^"), (False, True, True), (0, -0.22, 0.22)):
            if value is None:
                continue
            right.scatter(value, y + offset, s=24, marker=marker, zorder=3, linewidth=1.1 if not filled else 0.8,
                          facecolor=colour if filled else t["surface"], edgecolor=colour if not filled else t["surface"])
        if None in values:
            right.text(min(present) - 0.015, y, "cut: too few counts", color=t["muted"], va="center", ha="right",
                       fontsize=6)
    left.set_yticks(range(len(names)))
    left.set_yticklabels(names)
    left.invert_yaxis()
    left.set_xlim(0, 1)
    left.set_xlabel("Bray–Curtis dissimilarity of genus profiles")
    left.set_title("Change from cutting, against the spread of one environment", loc="left", pad=6)
    right.set_xlim(0.5, 1)
    right.set_xlabel("SourceTracker2 share of the library's own environment")
    right.set_title("Attribution with the library's study left out", loc="left", pad=6)
    for ax, text in ((left, "a"), (right, "b")):
        clean(ax, t)
        letter(ax, text)
    ink = t["secondary"]
    left.legend(handles=[Line2D([], [], marker="o", linestyle="", color=ink, markersize=4.5,
                                label="cut to ancient lengths, damaged"),
                         Line2D([], [], marker="o", linestyle="", markerfacecolor=t["surface"], markeredgecolor=ink,
                                markersize=4.5, label="another library, same environment")],
                loc="lower right", frameon=False, handletextpad=0.2)
    right.legend(handles=[Line2D([], [], marker="o", linestyle="", markerfacecolor=t["surface"], markeredgecolor=ink,
                                 markersize=4.5, label="modern"),
                          Line2D([], [], marker="o", linestyle="", color=ink, markersize=4.5, label="cut"),
                          Line2D([], [], marker="^", linestyle="", color=ink, markersize=4.5, label="cut, damaged")],
                 loc="lower left", frameon=False, handletextpad=0.2)
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.18, top=0.86, wspace=0.08)
    return fig


def counting(t):
    rows = read("counting.tsv")
    labels = ["Bracken", "≥ 10 reads", "≥ 5 reads", "every read"]
    panels = [("sourcetracker2_correct", "SourceTracker2", "Right on held-out sources", (0, 1.05)),
              ("feast_correct", "FEAST", None, (0, 1.05)), ("fastst_correct", "FastST", None, (0, 1.05)),
              ("sourcetracker2_mixture_error", "SourceTracker2, mixtures", "Largest error on the mixtures",
               (0, 0.13))]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.1))
    x = range(len(rows))
    for index, (ax, (column, title, ylabel, ylim)) in enumerate(zip(axes, panels)):
        values = [float(r[column]) for r in rows]
        ax.plot(x, values, color=t["axis"], linewidth=0.9, zorder=1)
        ax.scatter(x, values, s=22, color=t["ink"], edgecolor=t["surface"], linewidth=0.8, zorder=3)
        ax.scatter([len(rows) - 1], [values[-1]], s=46, facecolor="none", edgecolor=t["ink"], linewidth=0.8,
                   zorder=4)
        if column.endswith("error"):
            ax.axhline(0.10, color=t["muted"], linewidth=0.6, zorder=0)
            ax.text(len(rows) - 1, 0.103, "check fails above", color=t["muted"], fontsize=6, va="bottom", ha="right")
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylim(*ylim)
        ax.set_title(title, loc="left", pad=6)
        if ylabel:
            ax.set_ylabel(ylabel)
        ax.tick_params(axis="y", labelsize=6.5, colors=t["muted"])
        clean(ax, t, xgrid=False)
        ax.grid(axis="y", color=t["grid"], linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        letter(ax, "abcd"[index])
    fig.text(0.5, 0.01, "Counts: Bracken, or the classifier's own, with a minimum of reads per taxon and sample "
             "or none (circled: MSTFlow's choice)", ha="center", color=t["secondary"], fontsize=6.5)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.36, top=0.86, wspace=0.42)
    return fig


def alpha2(t):
    rows = read("alpha2.tsv")
    names = {"calculus": "calculus (ancient)", "sediment": "paleosoil", "coprolite": "coprolite", "bone": "bone",
             "plaque": "plaque (modern)", "faeces": "faeces (modern)", "mix_oral_gut": "plaque + faeces",
             "mix_oral_gut_ancient": "the same, ancient-like", "mix_gut_soil": "faeces + meadow soil"}
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for y, r in enumerate(rows):
        low, high = float(r["unknown_alpha2_0.001"]), float(r["unknown_alpha2_0.1"])
        ax.plot([high, low], [y, y], color=t["axis"], linewidth=0.8, zorder=1)
        ax.scatter(low, y, s=22, facecolor=t["surface"], edgecolor=t["ink"], linewidth=1.0, zorder=3)
        ax.scatter(high, y, s=22, color=t["ink"], edgecolor=t["surface"], linewidth=0.8, zorder=3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([names.get(r["sink"], r["sink"]) for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("SourceTracker2 share given to the unknown source")
    clean(ax, t)
    ink = t["secondary"]
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", markerfacecolor=t["surface"], markeredgecolor=ink,
                              markersize=4.5, label="alpha2 = 0.001 (SourceTracker2 2.0.1)"),
                       Line2D([], [], marker="o", linestyle="", color=ink, markersize=4.5,
                              label="alpha2 = 0.1 (original SourceTracker)")],
              loc="upper center", bbox_to_anchor=(0.4, 1.2), ncol=1, frameon=False, handletextpad=0.2)
    fig.subplots_adjust(left=0.33, right=0.97, bottom=0.17, top=0.82)
    return fig


def main():
    figures = {"mixtures": mixtures, "ancient-like": ancient_like, "counting": counting, "alpha2": alpha2}
    for mode, t in THEME.items():
        style(t)
        for name, draw in figures.items():
            fig = draw(t)
            # No dates in the files, so redrawing unchanged data leaves them byte-identical.
            fig.savefig(HERE / f"validation-{name}-{mode}.svg", metadata={"Date": None})
            if mode == "light":
                fig.savefig(HERE / f"validation-{name}.pdf", metadata={"CreationDate": None})
            plt.close(fig)
            print("wrote", f"validation-{name}-{mode}.svg")


if __name__ == "__main__":
    main()
