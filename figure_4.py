"""
Fig. 4 for the iPiccel954 manuscript (nitrogen and carbon coupling, WT-Air and WT-CO2, urea as the nitrogen source).

Usage
    python make_fig4.py ecGEMs.zip protein_abundance_g_per_gCDW.xlsx iPceleri.fasta Fig4_NC_coupling.png

Inputs
    ecGEMs.zip                         models, reaction upper bounds and kcat values of the two conditions
    protein_abundance_g_per_gCDW.xlsx  iBAQ proteomics, columns Wt_CO2 and Wt_Air (g/gCDW), one row per gene
    iPceleri.fasta                     protein sequences and annotations (molecular weights, gene check)

Everything is computed here. The coupling model is rebuilt from the files and solved with cobrapy.
  1. Urea is the nitrogen source, every other nitrogen-containing import is closed. Valine import stays open because
     the model needs it to grow. The two gene-less urea shortcuts (uahCv, ATAH_REV) are blocked and GOGAT is NADH dependent.
  2. The upper bounds of urea carboxylase (URCB) and allophanate hydrolase (ALPH) are rebuilt from the genes that carry
     those annotations in the fasta file, bound = kcat x measured abundance / molecular weight. The kcat values are the
     TurNuP values stored for the two reactions in kcat_mw_*.xlsx.
Requires cobra, matplotlib, numpy, pandas and openpyxl.
"""
import os
import sys
import tempfile
import warnings
import zipfile

import cobra
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.flux_analysis import pfba
from matplotlib import font_manager
from matplotlib.colors import to_rgb
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")
zip_path = sys.argv[1] if len(sys.argv) > 1 else "ecGEMs.zip"
prot_path = sys.argv[2] if len(sys.argv) > 2 else "protein_abundance_g_per_gCDW.xlsx"
fasta_path = sys.argv[3] if len(sys.argv) > 3 else "iPceleri.fasta"
out_path = sys.argv[4] if len(sys.argv) > 4 else "Fig4_NC_coupling.png"

BIO = "Biomass_Picre_Auto"
COND = ["WT-Air", "WT-CO2"]
FOLDER = {"WT-Air": ("WT_Air", "Wt_Air"), "WT-CO2": ("Wt_CO2", "Wt_CO2")}
PCOL = {"WT-Air": "Wt_Air", "WT-CO2": "Wt_CO2"}
EXP_GROWTH = {"WT-Air": 0.028, "WT-CO2": 0.324941}
ANNOT = {"URCB": "GENE4066869", "ALPH": "GENE4066871"}  # fasta annotation urea carboxylase, allophanate hydrolase

# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------
tmp = tempfile.mkdtemp(prefix="ecgems_")
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(tmp)
root = tmp
for r_, dirs, _ in os.walk(tmp):
    if "WT_Air" in dirs:
        root = r_
        break

names, seqs, cur = {}, {}, None
for line in open(fasta_path):
    line = line.strip()
    if line.startswith(">"):
        h = line[1:].split(" ", 1)
        cur = h[0].split("_")[0]
        names[cur] = h[1] if len(h) > 1 else ""
        seqs[cur] = ""
    elif cur:
        seqs[cur] += line
AA = {"A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388, "E": 129.1155, "Q": 128.1307,
      "G": 57.0519, "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766,
      "P": 97.1167, "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326}


def mw(gene):
    return sum(AA.get(c, 110.0) for c in seqs[gene]) + 18.02


pa = pd.read_excel(prot_path).iloc[:, :3].set_index("Protein")
for rid, g in ANNOT.items():
    print(f"{rid} uses {g} ({names[g]}), MW {mw(g):.0f} Da")
    assert "urea carboxylase" in names[g].lower() or "allophanate hydrolase" in names[g].lower()

kcat = {}
for c in COND:
    f, t = FOLDER[c]
    kc = pd.read_excel(f"{root}/{f}/relaxation_algorithm/kcat_mw_{t}.xlsx", sheet_name="reaction_protein_cost")
    kcat[c] = {r: float(kc.loc[kc.Reaction == r, "Kcat"].iloc[0]) for r in ANNOT}  # 1/h


def capacity(c, rid):
    """mmol/gDW/h from kcat (1/h), abundance (g/gCDW) and molecular weight (g/mol)."""
    g = ANNOT[rid]
    return kcat[c][rid] * pa.loc[g, PCOL[c]] / mw(g) * 1000.0


# ---------------------------------------------------------------------------
# Coupling model
# ---------------------------------------------------------------------------
def coupling_model(c, override=None):
    f, t = FOLDER[c]
    m = cobra.io.read_sbml_model(f"{root}/{f}/irrev_iPiccel954.xml")
    ub = pd.read_excel(f"{root}/{f}/relaxation_algorithm/UB_{t}_revised.xlsx", header=None, usecols=[0, 1],
                       names=["r", "u"])
    ids = {r.id for r in m.reactions}
    for r, u in zip(ub.r, ub.u):
        if str(r) in ids and pd.notna(u):
            m.reactions.get_by_id(str(r)).upper_bound = max(0.0, float(u))
    for rid, u in (override or {}).items():
        m.reactions.get_by_id(rid).upper_bound = u
    keep = {"EX_urea(e)_REV", "EX_val-L(e)"}
    for r in m.exchanges:
        if r.id not in keep and any("N" in a.elements for a in r.metabolites):
            r.upper_bound = 0.0
    for rid in ("uahCv", "ATAH_REV"):
        m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    if "GLUS(nadph)" in m.reactions:      # older model versions, make GOGAT NADH dependent
        nadph, nadp, nadh, nad = [m.metabolites.get_by_id(x) for x in ("nadph_h", "nadp_h", "nadh_h", "nad_h")]
        m.reactions.get_by_id("GLUS(nadph)").add_metabolites({nadph: 1, nadp: -1, nadh: -1, nad: 1})
    m.objective = BIO
    return m


FL = {}      # pFBA fluxes with the annotated-gene bounds
GROW = {}    # growth under the four capacity scenarios
CAPUSED, CN, UBV = {}, {}, {}
for c in COND:
    ov = {rid: capacity(c, rid) for rid in ANNOT}
    UBV[c] = ov
    m = coupling_model(c, ov)
    m.optimize()
    p = pfba(m)
    FL[c] = p.fluxes
    CAPUSED[c] = {rid: 100.0 * p.fluxes[rid] / ov[rid] for rid in ANNOT}
    Cn = Nn = 0.0
    for r in m.exchanges:
        v = p.fluxes[r.id]
        for a, co in r.metabolites.items():
            Cn += co * v * a.elements.get("C", 0)
            Nn += co * v * a.elements.get("N", 0)
    CN[c] = Cn / Nn
    GROW[c] = []
    for o in (ov, {**ov, "URCB": 1000.0}, {**ov, "ALPH": 1000.0}, {**ov, "URCB": 1000.0, "ALPH": 1000.0}):
        GROW[c].append(coupling_model(c, o).optimize().objective_value)
    print(c, "growth", [f"{g:.5f}" for g in GROW[c]], "C:N", round(CN[c], 2), "capacity used",
          {k: round(v) for k, v in CAPUSED[c].items()})

ENZ = [  # label, short, gene(s), node
    ("Urea carboxylase", "URC", ["GENE4066869"], 1),
    ("Allophanate hydrolase", "ALPH", ["GENE4066871"], 2),
    ("Pyruvate carboxylase", "PC", ["GENE4066235"], 1),
    ("PEP carboxylase", "PEPC", ["GENE4062117"], 2),
    ("Glutamate synthase (NADH)", "GOGAT", ["GENE4062926"], 3),
    ("Aspartate transaminase", "ASPAT", ["GENE4066933"], 3),
]
abund = {c: [1000.0 * sum(pa.loc[g, PCOL[c]] for g in genes) for _, _, genes, _ in ENZ] for c in COND}
fold = [abund["WT-CO2"][k] / abund["WT-Air"][k] for k in range(len(ENZ))]
g_exp = EXP_GROWTH["WT-CO2"] / EXP_GROWTH["WT-Air"]
g_mod = GROW["WT-CO2"][0] / GROW["WT-Air"][0]
print("fold abundance", [round(x, 2) for x in fold], "growth fold experiment", round(g_exp, 1), "coupling model",
      round(g_mod, 1))

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
avail = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in avail), "DejaVu Sans")
plt.rcParams.update({
    "font.family": family, "font.size": 11, "axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.major.size": 4.2, "ytick.major.size": 4.2, "xtick.labelsize": 10.2, "ytick.labelsize": 10.2,
    "axes.labelsize": 10.7, "legend.frameon": False, "legend.fontsize": 9.9, "text.color": "#000000",
    "axes.labelcolor": "#000000", "axes.edgecolor": "#000000", "xtick.color": "#000000", "ytick.color": "#000000",
})
INK, AXIS, GREY, TRACK = "#000000", "#000000", "#b9b8b0", "#eef0f3"
BLUE, ORANGE = "#2a78d6", "#eb6834"
N1, N2, N3 = "#eda100", "#1baf7a", "#e87ba4"
NODE = {1: N1, 2: N2, 3: N3}
CC = {"WT-Air": BLUE, "WT-CO2": ORANGE}
BAND = "#fbf1d3"


def tint(c, a):
    r_, g_, b_ = to_rgb(c)
    return (1 - a) + a * r_, (1 - a) + a * g_, (1 - a) + a * b_


fig = plt.figure(figsize=(11.9, 9.9), dpi=100)
gs = fig.add_gridspec(3, 12, height_ratios=[1.0, 1.0, 1.05], left=0.075, right=0.985, top=0.965, bottom=0.065,
                      hspace=0.62, wspace=2.6)
HEADS = []


def style(ax):
    ax.minorticks_off()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK, pad=1.5)


def dot(c):
    return Line2D([], [], marker="o", ls="", ms=9.3, mfc=c, mec=c)


def legend(ax, entries, **kw):
    hs = [tuple(dot(c) for c in cols) for cols, _ in entries]
    ax.legend(hs, [t for _, t in entries], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
              handletextpad=0.4, labelspacing=0.35, borderpad=0, fontsize=9.5, **kw)


def vbar(ax, x, v, col, lw=20.2, base=0.0):
    ax.plot([x, x], [base, v], color=col, lw=lw, solid_capstyle="round", zorder=3)


def hbar(ax, y, v, col, lw=10.2, base=0.0):
    ax.plot([base, v], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=3)


# A. internal C:N
axA = fig.add_subplot(gs[0, 0:4])
HEADS.append((axA, "A", "Proteome-derived capacity (mmol/gDW/h)"))
for gi, rid in enumerate(["URCB", "ALPH"]):
    for ci, c in enumerate(COND):
        x = gi + (ci - 0.5) * 0.36
        v = UBV[c][rid]
        vbar(axA, x, v, CC[c], lw=15.3)
        axA.annotate(f"{v:.3f}", (x, v), xytext=(0, 10.2), textcoords="offset points", ha="center", va="bottom", fontsize=9.2)
axA.set_xlim(-0.55, 1.55)
axA.set_ylim(0, 0.62)
axA.set_xticks([0, 1])
axA.set_xticklabels(["Urea carboxylase", "Allophanate hydrolase"], fontsize=9.5)
axA.set_yticks([0, 0.2, 0.4, 0.6])
style(axA)
legend(axA, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left", bbox_to_anchor=(0.0, 1.04))

# B. measured abundance
axB = fig.add_subplot(gs[0, 4:8])
HEADS.append((axB, "B", "Measured abundance (mg/gCDW)"))
n = len(ENZ)
for k in range(n):
    for ci, c in enumerate(COND):
        x = k + (ci - 0.5) * 0.36
        v = abund[c][k]
        vbar(axB, x, v, CC[c], lw=7.1)
        axB.annotate(f"{v:.2f}", (x, v), xytext=(0, 9.3), textcoords="offset points", ha="center", va="bottom", fontsize=7.5, rotation=90)
axB.set_xlim(-0.5, n - 0.5)
axB.set_ylim(0, 1.55)
axB.set_xticks(range(n))
axB.set_xticklabels([e[1] for e in ENZ], fontsize=8.8)
axB.set_yticks([0, 0.5, 1.0])
style(axB)

# C. fold change
axC = fig.add_subplot(gs[0, 8:12])
HEADS.append((axC, "C", "WT-CO$_2$ / WT-Air (fold)"))
rows = [("Growth, experiment", g_exp, INK), ("Growth, model", g_mod, "#6b6b6b")]
rows += [(e[1], fold[k], NODE[e[3]]) for k, e in enumerate(ENZ)]
m_ = len(rows)
for i, (lab, v, col) in enumerate(rows):
    y = m_ - 1 - i
    hbar(axC, y, v, col, lw=9.5)
    axC.annotate(f"{v:.1f}", (v, y), xytext=(11.9, 0), textcoords="offset points", va="center", fontsize=9.2)
axC.set_yticks(range(m_))
axC.set_yticklabels([r[0] for r in rows][::-1], fontsize=9.2)
axC.set_xlim(0, 15)
axC.set_ylim(-0.6, m_ - 0.35)
axC.set_xticks([0, 4, 8, 12])
style(axC)
axC.spines["left"].set_visible(False)
axC.tick_params(axis="y", length=0, pad=2)
legend(axC, [([N1], "Node 1"), ([N2], "Node 2"), ([N3], "Node 3")], loc="lower right", bbox_to_anchor=(1.02, 0.0))
p_ = axC.get_position()
axC.set_position([p_.x0 + 0.035, p_.y0, p_.width - 0.035, p_.height])


# D to F. node fluxes
def node_panel(ax, spec, ymax, yticks, letter_, title):
    HEADS.append((ax, letter_, title))
    ng = len(spec)
    for gi, (lab, key) in enumerate(spec):
        for ci, c in enumerate(COND):
            x = gi + (ci - 0.5) * 0.36
            v = float(FL[c][key])
            vbar(ax, x, v, CC[c], lw=10.2)
            ax.annotate(f"{v:.3f}" if v >= 0.0005 else "<0.001", (x, v), xytext=(0, 9.3), textcoords="offset points", ha="center", va="bottom", fontsize=8.2, rotation=90)
    ax.set_xlim(-0.55, ng - 0.45)
    ax.set_ylim(0, ymax)
    ax.set_xticks(range(ng))
    ax.set_xticklabels([s[0] for s in spec], fontsize=9.2)
    ax.set_yticks(yticks)
    ax.set_ylabel("mmol/gDW/h", labelpad=3, fontsize=9.9)
    style(ax)


axD = fig.add_subplot(gs[1, 0:4])
node_panel(axD, [("URCB\nN", "URCB"), ("PCr\nC", "PCr"), ("ACCOAC\nlipid", "accoaciCv")], 0.50, [0, 0.1, 0.2, 0.3], "D",
           "Node 1, bicarbonate pool")
axE = fig.add_subplot(gs[1, 4:8])
node_panel(axE, [("ALPH\nN release", "ALPH"), ("PPC\nre-fix", "PPC"), ("CO$_2$\nefflux", "CO2t_REV")], 0.46,
           [0, 0.1, 0.2, 0.3], "E", "Node 2, CO$_2$ recapture")
axF = fig.add_subplot(gs[1, 8:12])
node_panel(axF, [("GOGAT\nconsumes", "GLUS(nadh)"), ("ASPATh\nsupplies", "ASPATh_FWD"), ("TYRTAh\nsupplies", "TYRTAh")],
           1.30, [0, 0.3, 0.6, 0.9], "F", "Node 3, 2-oxoglutarate")

# G. capacity used
axG = fig.add_subplot(gs[2, 0:4])
HEADS.append((axG, "G", "Enzyme capacity used (%)"))
rows_g = [("URCB", "Urea\ncarboxylase"), ("ALPH", "Allophanate\nhydrolase")]
for ri, (rid, lab) in enumerate(rows_g):
    y = 1 - ri
    for ci, c in enumerate(COND):
        yy = y + (0.19 if ci == 0 else -0.19)
        v = min(CAPUSED[c][rid], 100.0)
        axG.plot([0, 100], [yy, yy], color=TRACK, lw=12.8, solid_capstyle="round", zorder=1)
        hbar(axG, yy, max(v, 1.5), CC[c], lw=12.8)
        axG.text(104, yy, f"{v:.0f}%", va="center", fontsize=9.9)
axG.set_yticks([1, 0])
axG.set_yticklabels([r[1] for r in rows_g])
axG.set_ylim(-0.7, 1.5)
axG.set_xlim(-3, 128)
axG.set_xticks([0, 50, 100])
style(axG)
axG.spines["left"].set_visible(False)
axG.tick_params(axis="y", length=0, pad=2)

# H. growth when capacity is relaxed
axH = fig.add_subplot(gs[2, 4:12])
HEADS.append((axH, "H", "Coupling-model growth when capacity is relaxed (1/h)"))
groups = ["Measured\ncapacity", "Urea carboxylase\nrelaxed", "Allophanate hydrolase\nrelaxed", "Both\nrelaxed", "Experiment"]
YB = 5e-5
for gi in range(5):
    for ci, c in enumerate(COND):
        x = gi + (ci - 0.5) * 0.34
        v = GROW[c][gi] if gi < 4 else EXP_GROWTH[c]
        col = CC[c] if gi < 4 else tint(CC[c], 0.45)
        axH.plot([x, x], [YB * 1.2, v], color=col, lw=15.3, solid_capstyle="round", zorder=3)
        axH.annotate((f"{v:.2g}" if v < 0.01 else f"{v:.3f}"), (x, v), xytext=(0, 11), textcoords="offset points",
                     ha="center", va="bottom", fontsize=8.5, rotation=90, annotation_clip=False)
axH.set_yscale("log")
axH.set_ylim(YB, 3.0)
axH.set_xlim(-0.5, 4.5)
axH.set_xticks(range(5))
axH.set_xticklabels(groups, fontsize=9.2)
axH.set_yticks([1e-4, 1e-3, 1e-2, 1e-1])
axH.set_yticklabels(["10$^{-4}$", "10$^{-3}$", "10$^{-2}$", "10$^{-1}$"])
style(axH)
legend(axH, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left", bbox_to_anchor=(0.0, 1.03))

# letters and titles
for ax, txt, title in HEADS:
    p = ax.get_position()
    x0 = 0.008 if p.x0 < 0.2 else p.x0 - 0.062
    fig.text(x0, p.y1 + 0.012, txt, fontsize=15.3, fontweight="bold", va="bottom", ha="left")
    fig.text(x0 + 0.022, p.y1 + 0.013, title, fontsize=11.2, va="bottom", ha="left")

fig.savefig(out_path, dpi=450, facecolor="white")
print("saved", out_path)