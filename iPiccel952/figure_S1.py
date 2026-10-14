"""
Supplementary Fig. S1 for the iPiccel952 manuscript (supports Fig. 1, model content and quality).

Usage
    python make_figS1.py iPiccel952.xml FigS1.png
where iPiccel952.xml is the final SBML model, or edit MODEL_PATH and OUT_PATH below.
Requires cobra, matplotlib, numpy and pandas. Nothing is cached and the run takes about a minute.

Panels
    A  largest biomass components, by mass
    B  residual mass and charge imbalance of the final model, by element
    C  reactions per pathway
    D  share of reactions with a gene association, by reaction class
    E  genes per reaction and reactions per gene
    F  most connected metabolites
"""
import sys
import warnings

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

MODEL_PATH = "iPiccel952.xml"
OUT_PATH = "FigS1.png"
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
    "font.family": family, "font.size": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
    "axes.labelsize": 6.3, "legend.frameon": False, "legend.fontsize": 5.8, "text.color": "#000000",
    "axes.labelcolor": "#000000", "xtick.color": "#000000", "ytick.color": "#000000", "mathtext.fontset": "stix",
})
INK, AXIS, GREY = "#000000", "#7a7a7a", "#b9b8b0"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#0b8a0b"


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
    return Line2D([], [], marker="o", ls="", ms=5.5, mfc=c, mec=c)


def legend(ax, entries, **kw):
    hs = [tuple(dot(c) for c in cols) for cols, _ in entries]
    ax.legend(hs, [t for _, t in entries], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
              handletextpad=0.4, labelspacing=0.35, borderpad=0, fontsize=5.6, **kw)


def hbar(ax, y, v, col, lw=5.0, b=0.0):
    ax.plot([b, v], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=3)


def vbar(ax, x, v, col, lw=7.0, b=0.0):
    ax.plot([x, x], [b, v], color=col, lw=lw, solid_capstyle="round", zorder=3)


def val(ax, x, y, s_, dx=4, dy=0, **kw):
    kw.setdefault("fontsize", 5.0)
    ax.annotate(s_, (x, y), xytext=(dx, dy), textcoords="offset points", annotation_clip=False, **kw)


def hpanel(ax, labels, vals, cols, fmt="{:.0f}", xmax=None, lw=5.0, fs=5.2, extra=None):
    n = len(labels)
    for i, (v, c) in enumerate(zip(vals, cols)):
        y = n - 1 - i
        hbar(ax, y, v, c, lw)
        val(ax, v, y, fmt.format(v) + (extra[i] if extra else ""), va="center", ha="left")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels[::-1], fontsize=fs)
    ax.set_ylim(-0.7, n - 0.3)
    if xmax:
        ax.set_xlim(0, xmax)
    style(ax, left=False)


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
        fig.text(x0, p.y1 + 0.012, letter, fontsize=9, fontweight="bold", va="bottom", ha="left")
        fig.text(x0 + 0.021, p.y1 + 0.013, title, fontsize=6.6, va="bottom", ha="left")
    fig.savefig(path, dpi=400, facecolor="white")
    plt.close(fig)
    print("saved", path, flush=True)


# ---------------------------------------------------------------------------
# Reaction classes
# ---------------------------------------------------------------------------
is_bio = lambda r: r.id == BIO
is_ex = lambda r: r.boundary
is_tr = lambda r: (not r.boundary) and (not is_bio(r)) and len(r.compartments) > 1

fig = plt.figure(figsize=(7.0, 5.9), dpi=100, layout="constrained")
fig.get_layout_engine().set(w_pad=0.04, h_pad=0.06, wspace=0.12, hspace=0.14, rect=(0.0, 0.0, 1.0, 0.955))
gs = fig.add_gridspec(2, 3)

# ---------------------------------------------------------------------------
# A. Largest biomass components
# ---------------------------------------------------------------------------
b = M.reactions.get_by_id(BIO)
bm = sorted(((abs(c) * x.formula_weight, x.name) for x, c in b.metabolites.items()
             if c < 0 and x.formula and x.name not in ("ATP", "H2O")), reverse=True)
topb = bm[:12]
axA = fig.add_subplot(gs[0, 0])
head(axA, "A", "Largest biomass components (no ATP, H$_2$O)")
hpanel(axA, [n_[:24] for _, n_ in topb], [v for v, _ in topb], [BLUE] * len(topb), fmt="{:.0f}",
       xmax=1.25 * topb[0][0], lw=4.6, fs=4.9)
axA.set_xlabel("Mass in biomass (mg/gDW)", labelpad=2)

# ---------------------------------------------------------------------------
# B. Residual mass and charge imbalance of the final model
# ---------------------------------------------------------------------------
imb = {"C": 0, "H": 0, "O": 0, "N": 0, "P": 0, "S": 0, "Charge": 0, "Other": 0}
tot_imb = 0
for r in M.reactions:
    if r.boundary or is_bio(r):
        continue
    try:
        bal = r.check_mass_balance()
    except Exception:
        continue
    if not bal:
        continue
    tot_imb += 1
    for k in bal:
        kk = "Charge" if k == "charge" else k if k in imb else "Other"
        imb[kk] += 1
axB = fig.add_subplot(gs[0, 1])
head(axB, "B", "Residual imbalance by element")
ks = list(imb)
for i, k in enumerate(ks):
    vbar(axB, i, imb[k], BLUE, 8.0)
    val(axB, i, imb[k], f"{imb[k]}", dx=0, dy=4, ha="center", va="bottom")
axB.set_xticks(range(len(ks)))
axB.set_xticklabels(ks, fontsize=5.0, rotation=35, ha="right")
axB.set_ylim(0, 1.2 * max(max(imb.values()), 1))
axB.set_xlim(-0.6, len(ks) - 0.4)
axB.set_ylabel(f"Reactions ({tot_imb} imbalanced)", labelpad=2)
style(axB)
print("imbalanced reactions in the final model", tot_imb, imb, flush=True)

# ---------------------------------------------------------------------------
# C. Reactions per pathway
# ---------------------------------------------------------------------------
cnt = pd.Series({g_.name: len(g_.members) for g_ in M.groups}).sort_values(ascending=False)
top = cnt.head(14)
axC = fig.add_subplot(gs[0, 2])
head(axC, "C", "Reactions per pathway")
for i, (k, v) in enumerate(top.items()):
    y = len(top) - 1 - i
    hbar(axC, y, np.log10(v), ORANGE if k == "Fatty acid metabolism" else BLUE, 4.6)
    val(axC, np.log10(v), y, f"{v:,}", va="center", ha="left", fontsize=4.8)
axC.set_yticks(range(len(top)))
axC.set_yticklabels([short(k) for k in top.index][::-1], fontsize=4.9)
axC.set_ylim(-0.7, len(top) - 0.3)
axC.set_xlim(0, 3.9)
axC.set_xticks([0, 1, 2, 3])
axC.set_xticklabels(["1", "10", "100", "1000"])
axC.set_xlabel("Reactions (log scale)", labelpad=2)
style(axC, left=False)

# ---------------------------------------------------------------------------
# D. Gene associations by reaction class
# ---------------------------------------------------------------------------
cls = {"Metabolic": [], "Transport": [], "Exchange and demand": [], "Biomass": []}
for r in M.reactions:
    k = "Biomass" if is_bio(r) else "Exchange and demand" if is_ex(r) else "Transport" if is_tr(r) else "Metabolic"
    cls[k].append(bool(r.gene_reaction_rule))
axD = fig.add_subplot(gs[1, 0])
head(axD, "D", "Reactions with a gene association")
pct = [100 * np.mean(v) if v else 0 for v in cls.values()]
hpanel(axD, list(cls), pct, [BLUE, GREEN, GREY, GREY], fmt="{:.0f}", xmax=118, lw=7.0,
       extra=[f"%  (n = {len(v):,})" for v in cls.values()])
axD.set_xlabel("Share of reactions (%)", labelpad=2)

# ---------------------------------------------------------------------------
# E. Genes per reaction and reactions per gene
# ---------------------------------------------------------------------------
gpr = [len(r.genes) for r in M.reactions if r.genes]
rpg = [len(g.reactions) for g in M.genes]
bins = [1, 2, 3, 4, 5]
f1 = [100 * sum((x == b_) if b_ < 5 else (x >= 5) for x in gpr) / len(gpr) for b_ in bins]
f2 = [100 * sum((x == b_) if b_ < 5 else (x >= 5) for x in rpg) / len(rpg) for b_ in bins]
axE = fig.add_subplot(gs[1, 1])
head(axE, "E", "Gene and reaction multiplicity")
for i in range(5):
    vbar(axE, i - 0.17, f1[i], BLUE, 5.0)
    vbar(axE, i + 0.17, f2[i], ORANGE, 5.0)
    val(axE, i - 0.17, f1[i], f"{f1[i]:.0f}", dx=0, dy=3, ha="center", va="bottom", fontsize=4.6)
    val(axE, i + 0.17, f2[i], f"{f2[i]:.0f}", dx=0, dy=3, ha="center", va="bottom", fontsize=4.6)
axE.set_xticks(range(5))
axE.set_xticklabels(["1", "2", "3", "4", "5 or more"], fontsize=5.2)
axE.set_ylim(0, 1.3 * max(f1 + f2))
axE.set_xlim(-0.6, 4.6)
axE.set_ylabel("Share (%)", labelpad=2)
axE.set_xlabel("Number per reaction or per gene", labelpad=2)
style(axE)
legend(axE, [([BLUE], "Genes per reaction"), ([ORANGE], "Reactions per gene")], loc="upper right")

# ---------------------------------------------------------------------------
# F. Most connected metabolites
# ---------------------------------------------------------------------------
deg = sorted(((len(x.reactions), x.name, x.compartment) for x in M.metabolites if x.name not in ("H+", "H2O")),
             reverse=True)[:12]
axF = fig.add_subplot(gs[1, 2])
head(axF, "F", "Most connected metabolites")
hpanel(axF, [f"{n_[:20]} [{c}]" for _, n_, c in deg], [d_ for d_, _, _ in deg], [BLUE] * 12, fmt="{:.0f}",
       xmax=1.25 * deg[0][0], lw=4.6, fs=4.9)
axF.set_xlabel("Reactions per metabolite", labelpad=2)

finish(fig, out_path)