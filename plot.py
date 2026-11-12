"""Publication-style tree figures.

A rectangular phylogram drawn directly with matplotlib so that query isolates
can be highlighted, support values placed on branches, and a substitutions-per-
site scale bar added.
"""
from __future__ import annotations

import io
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo


def _load(newick: str):
    return Phylo.read(io.StringIO(newick if newick.endswith(";") else newick + ";"),
                      "newick")


def plot_tree(newick: str, out_path: str | Path, query_labels: set[str] | None = None,
              title: str = "", support_threshold: float = 70.0,
              support_label: str = "UFBoot", scale_label: str = "substitutions/site",
              width: float = 9.0, row_height: float = 0.26, dpi: int = 300):
    """Draw a rooted phylogram; returns the output path.

    Tips whose label is in `query_labels` are coloured and bold - these are the
    isolates under study; everything else is a reference sequence.
    """
    tree = _load(newick)
    query_labels = set(query_labels or ())
    terminals = tree.get_terminals()
    n_tips = len(terminals)

    # ---- layout ---------------------------------------------------------
    ypos, xpos = {}, {}
    for i, t in enumerate(terminals):
        ypos[t] = n_tips - i

    def depth(clade, acc=0.0):
        xpos[clade] = acc
        for ch in clade.clades:
            depth(ch, acc + (ch.branch_length or 0.0))
    depth(tree.root, 0.0)

    def yof(clade):
        if clade in ypos:
            return ypos[clade]
        y = sum(yof(c) for c in clade.clades) / len(clade.clades)
        ypos[clade] = y
        return y
    yof(tree.root)

    fig_h = max(2.6, n_tips * row_height + 1.2)
    fig, ax = plt.subplots(figsize=(width, fig_h))

    for clade in tree.find_clades(order="level"):
        x0, y0 = xpos[clade], ypos[clade]
        for ch in clade.clades:
            x1, y1 = xpos[ch], ypos[ch]
            ax.plot([x0, x0], [y0, y1], color="0.25", lw=1.0, solid_capstyle="round")
            ax.plot([x0, x1], [y1, y1], color="0.25", lw=1.0, solid_capstyle="round")

    xmax = max(xpos.values()) or 1.0
    pad = xmax * 0.012
    for t in terminals:
        is_q = t.name in query_labels
        ax.text(xpos[t] + pad, ypos[t], str(t.name).replace("_", " "),
                va="center", ha="left", fontsize=8,
                color="#0b5394" if is_q else "0.15",
                fontweight="bold" if is_q else "normal")

    for clade in tree.get_nonterminals():
        sup = clade.confidence
        if sup is None and clade.name:
            try:
                sup = float(clade.name)
            except ValueError:
                sup = None
        if sup is None or sup < support_threshold or clade is tree.root:
            continue
        ax.text(xpos[clade] - pad * 0.6, ypos[clade] + 0.22, f"{sup:.0f}",
                ha="right", va="bottom", fontsize=6.5, color="#a33")

    # ---- scale bar ------------------------------------------------------
    import math
    raw = xmax / 5
    step = 10 ** math.floor(math.log10(raw)) if raw > 0 else 0.01
    bar = max(step, round(raw / step) * step)
    y0 = 0.25
    ax.plot([0, bar], [y0, y0], color="0.15", lw=1.4)
    ax.text(bar / 2, y0 - 0.45, f"{bar:g} {scale_label}", ha="center",
            va="top", fontsize=7.5, color="0.15")

    ax.set_xlim(-xmax * 0.02, xmax * 1.55)
    ax.set_ylim(-0.6, n_tips + 1.0)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    foot = f"bold blue = query isolates; numbers on branches = {support_label} >= {support_threshold:g}%"
    fig.text(0.01, 0.005, foot, fontsize=7, color="0.35")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path
