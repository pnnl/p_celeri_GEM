"""
Supplementary Fig. S2 for the iPiccel954 manuscript (supports Fig. 2, gene association and the comparison without fatty
acid reactions).

Usage
    python make_figS2.py iPiccel954.xml FigS2.png
where iPiccel954.xml is the final SBML model, or edit MODEL_PATH and OUT_PATH below.
Requires cobra, matplotlib, numpy and pandas. Nothing is cached and the run takes about a minute.
Values for iCO1515, iCre1355, iCZ843 and iAJ526 are taken from Table 1 of Soleymani et al. 2025 (New Phytologist 248,
1572-1583). Values for iPiccel954 are computed from the SBML file.

Panels
    A  gene association by pathway, with the number of reactions in brackets
    B  reactions with a gene association (%), iPiccel954 against the published models
    C  transport reactions (%)
    D  reversible reactions (%)
    E  reactions per gene
"""
import re
import sys
import warnings
import xml.etree.ElementTree as ET

import cobra
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

MODEL_PATH = "iPiccel954.xml"
OUT_PATH = "FigS2.png"
model_path = sys.argv[1] if len(sys.argv) > 1 else MODEL_PATH
out_path = sys.argv[2] if len(sys.argv) > 2 else OUT_PATH

BIO = "Biomass_Picre_Auto"
M = cobra.io.read_sbml_model(model_path)

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
avail = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in avail), "DejaVu Sans")
plt.rcParams.update({
    "font.family": family, "font.size": 11.0, "axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.major.size": 4, "ytick.major.size": 4, "xtick.labelsize": 10.2, "ytick.labelsize": 10.2,
    "axes.labelsize": 10.7, "legend.frameon": False, "legend.fontsize": 9.9, "text.color": "#000000",
    "axes.labelcolor": "#000000", "xtick.color": "#000000", "ytick.color": "#000000", "mathtext.fontset": "stix",
})
INK, AXIS, GREY = "#000000", "#7a7a7a", "#b9b8b0"
BLUE, ORANGE, GREEN, LGREEN = "#2a78d6", "#eb6834", "#0b8a0b", "#8fd18f"


def short(p):
    return (p.replace("Glycolysis / Gluconeogenesis", "Glycolysis").replace("Glycine, serine and threonine metabolism", "Gly, Ser, Thr")
            .replace(" metabolism", "").replace("Phenylalanine, tyrosine and tryptophan biosynthesis", "Aromatic amino acids")
            .replace("Oxidative phosphorylation", "Oxid. phosphorylation").replace("Cofactor and Prosthetic Group Biosynthesis", "Cofactor biosynthesis")
            .replace("Valine, leucine and isoleucine", "Val, Leu, Ile").replace("Pentose phosphate pathway", "Pentose phosphate")
            .replace("Alanine and aspartate", "Ala and Asp").replace("Glutamate and arginine", "Glu and Arg").strip())


def style(ax, left=True):
    ax.minorticks_off()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK, pad=1.5)
    if not left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=2)


def dot(c):
    return Line2D([], [], marker="o", ls="", ms=9.3, mfc=c, mec=c)


def hbar(ax, y, v, col, lw=8.0, b=0.0):
    ax.plot([b, v], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=3)


def val(ax, x, y, s_, dx=4, dy=0, **kw):
    kw.setdefault("fontsize", 8.8)
    ax.annotate(s_, (x, y), xytext=(dx, dy), textcoords="offset points", annotation_clip=False, **kw)


HEADS = []


def head(ax, letter, title):
    HEADS.append((ax, letter, title))


def finish(fig, path):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for ax, letter, title in HEADS:
        bb = ax.get_tightbbox(rend).transformed(fig.transFigure.inverted())
        p = ax.get_position()
        x0 = max(0.006, bb.x0 - 0.004)
        fig.text(x0, p.y1 + 0.012, letter, fontsize=15.3, fontweight="bold", va="bottom", ha="left")
        fig.text(x0 + 0.017, p.y1 + 0.013, title, fontsize=11.2, va="bottom", ha="left")
    fig.savefig(path, dpi=400, facecolor="white")
    plt.close(fig)
    print("saved", path, flush=True)


# ---------------------------------------------------------------------------
# Model statistics from the final SBML
# ---------------------------------------------------------------------------
grp = {}
for g_ in M.groups:
    for mem in g_.members:
        grp.setdefault(mem.id, []).append(g_.name)
sub_of = {k: v[0] for k, v in grp.items()}
FA = {r for r, s_ in sub_of.items() if s_ == "Fatty acid metabolism"}

is_bio = lambda r: r.id == BIO
is_tr = lambda r: (not r.boundary) and (not is_bio(r)) and len(r.compartments) > 1

_dec = lambda i: re.sub(r"__(\d+)__", lambda mo: chr(int(mo.group(1))), i[2:] if i.startswith("R_") else i)
REV_ATTR = {_dec(e.get("id")): e.get("reversible") == "true" for _, e in ET.iterparse(model_path) if e.tag.endswith("}reaction")}


def model_stats(rxns):
    rxns = list(rxns)
    genes = {g.id for r in rxns for g in r.genes}
    mets = {x.id for r in rxns for x in r.metabolites}
    return dict(rxn=len(rxns), genes=len(genes), met=len(mets), gpr=sum(bool(r.gene_reaction_rule) for r in rxns),
                transport=sum(is_tr(r) for r in rxns),
                reversible=sum(REV_ATTR.get(r.id, r.lower_bound < 0 < r.upper_bound) for r in rxns))


THIS = model_stats(M.reactions)
NOFA = model_stats([r for r in M.reactions if r.id not in FA])
print("iPiccel954", THIS, "\nwithout fatty acid reactions", NOFA, flush=True)

PUBLISHED = {
    "iCO1515": dict(genes=1515, rxn=2516, gpr=1975, transport=306, reversible=196, met=3517),
    "iCre1355": dict(genes=1460, rxn=2394, gpr=1924, transport=449, reversible=723, met=1845),
    "iCZ843": dict(genes=843, rxn=2294, gpr=1838, transport=378, reversible=626, met=1770),
    "iAJ526": dict(genes=531, rxn=1455, gpr=856, transport=317, reversible=561, met=1236),
}
ROWS = [("iPiccel954", THIS, GREEN), ("iPiccel954\nno fatty acid", NOFA, LGREEN)] + [(k, v, BLUE) for k, v in PUBLISHED.items()]

fig = plt.figure(figsize=(11.9, 9.0), dpi=100, layout="constrained")
fig.get_layout_engine().set(w_pad=0.04, h_pad=0.06, wspace=0.12, hspace=0.14, rect=(0.0, 0.0, 1.0, 0.955))
gs = fig.add_gridspec(2, 6)

# ---------------------------------------------------------------------------
# A. Gene association by pathway
# ---------------------------------------------------------------------------
axA = fig.add_subplot(gs[:, 0:2])
head(axA, "A", "Gene association by pathway")
rows = []
for g_ in M.groups:
    mem = list(g_.members)
    if not mem:
        continue
    rows.append((g_.name, len(mem), 100 * sum(bool(x.gene_reaction_rule) for x in mem) / len(mem)))
rows = sorted(rows, key=lambda t: -t[1])[:16]
for i, (nm, n_, p_) in enumerate(rows):
    y = len(rows) - 1 - i
    hbar(axA, y, p_, ORANGE if nm == "Fatty acid metabolism" else BLUE, 7.0)
    val(axA, p_, y, f"{p_:.0f}", va="center", ha="left", fontsize=8.2)
axA.set_yticks(range(len(rows)))
axA.set_yticklabels([f"{short(nm)} ({n_:,})" for nm, n_, _ in rows][::-1], fontsize=8.2)
axA.set_ylim(-0.7, len(rows) - 0.3)
axA.set_xlim(0, 125)
axA.set_xlabel("Reactions with a gene association (%)", labelpad=2)
style(axA, left=False)

# ---------------------------------------------------------------------------
# B to E. iPiccel954 against the published models
# ---------------------------------------------------------------------------
METRICS = [("Reactions with GPR (%)", lambda d_: 100 * d_["gpr"] / d_["rxn"], "{:.0f}"),
           ("Transport (%)", lambda d_: 100 * d_["transport"] / d_["rxn"], "{:.0f}"),
           ("Reversible (%)", lambda d_: 100 * d_["reversible"] / d_["rxn"], "{:.0f}"),
           ("Reactions per gene", lambda d_: d_["rxn"] / d_["genes"], "{:.1f}")]
pos = [(0, 2), (0, 4), (1, 2), (1, 4)]
for k, (title, f_, fmt) in enumerate(METRICS):
    ax = fig.add_subplot(gs[pos[k][0], pos[k][1]:pos[k][1] + 2])
    head(ax, "BCDE"[k], title)
    vals = [f_(d_) for _, d_, _ in ROWS]
    n = len(ROWS)
    for i, (v, (_, _, c)) in enumerate(zip(vals, ROWS)):
        y = n - 1 - i
        hbar(ax, y, v, c, 8.5)
        val(ax, v, y, fmt.format(v), va="center", ha="left", fontsize=8.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels([n_ for n_, _, _ in ROWS][::-1] if pos[k][1] == 2 else [""] * n, fontsize=8.3)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, 1.3 * max(vals))
    style(ax, left=False)

_ent = [([GREEN], "iPiccel954"), ([LGREEN], "iPiccel954 without fatty acid reactions"), ([BLUE], "Published algal models")]
fig.legend([tuple(dot(c) for c in cols) for cols, _ in _ent], [t for _, t in _ent],
           handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)}, loc="outside lower center", ncol=3, fontsize=9.5,
           frameon=False)

finish(fig, out_path)