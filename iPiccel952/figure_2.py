"""
Fig. 2 for the iPiccel952 manuscript (structural comparison with other algal GEMs).

Usage
    python make_fig2.py iPiccel952.xml Fig2.png
or edit XML_PATH and OUT_PATH below and run the file without arguments.

Requires matplotlib and numpy. Values for iCO1515, iCre1355, iCZ843 and iAJ526 are taken from Table 1 of
Soleymani et al. 2025 (New Phytologist 248, 1572-1583). Values for iPiccel952 are computed from the SBML file.
"""
import sys
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

# Set these two paths when you run the script from an editor without command-line arguments.
# On Windows use a raw string, for example r"C:\Users\name\Figures\iPiccel952.xml".
XML_PATH = "iPiccel952.xml"
OUT_PATH = "Fig2.png"
xml_path = sys.argv[1] if len(sys.argv) > 1 else XML_PATH
out_path = sys.argv[2] if len(sys.argv) > 2 else OUT_PATH

# ---------------------------------------------------------------------------
# Published models (Soleymani et al. 2025, Table 1)
# ---------------------------------------------------------------------------
PUBLISHED = {
    "iCO1515": dict(genes=1515, rxn=2516, gpr=1975, transport=306, reversible=196, comps=8, met=3517, unique=2539),
    "iCre1355": dict(genes=1460, rxn=2394, gpr=1924, transport=449, reversible=723, comps=11, met=1845, unique=1132),
    "iCZ843": dict(genes=843, rxn=2294, gpr=1838, transport=378, reversible=626, comps=6, met=1770, unique=457),
    "iAJ526": dict(genes=531, rxn=1455, gpr=856, transport=317, reversible=561, comps=6, met=1236, unique=770),
}

# ---------------------------------------------------------------------------
# iPiccel952 statistics from the SBML file
# ---------------------------------------------------------------------------
tree = ET.parse(xml_path)
root = tree.getroot()
SB = "{http://www.sbml.org/sbml/level3/version1/core}"
FBC = None
for _, (prefix, uri) in ET.iterparse(xml_path, events=["start-ns"]):
    if prefix == "fbc":
        FBC = "{" + uri + "}"
model = root.find(SB + "model")
GRP = "{http://www.sbml.org/sbml/level3/version1/groups/version1}"

species = {}
for s in model.find(SB + "listOfSpecies"):
    sid, comp = s.get("id"), s.get("compartment")
    base = sid[:-(len(comp) + 1)] if sid.endswith("_" + comp) else sid
    species[sid] = (comp, base)

n_rxn = n_gpr = n_transport = n_reversible = 0
has_gpr = {}
for r in model.find(SB + "listOfReactions"):
    n_rxn += 1
    has_gpr[r.get("id")] = r.find(FBC + "geneProductAssociation") is not None
    if r.find(FBC + "geneProductAssociation") is not None:
        n_gpr += 1
    if r.get("reversible") == "true":
        n_reversible += 1
    lr, lp = r.find(SB + "listOfReactants"), r.find(SB + "listOfProducts")
    reac = [x.get("species") for x in (lr if lr is not None else [])]
    prod = [x.get("species") for x in (lp if lp is not None else [])]
    if reac and prod and "biomass" not in r.get("id").lower():
        if len({species[i][0] for i in reac + prod}) > 1:
            n_transport += 1

genes = {g.get(FBC + "id") for g in model.find(FBC + "listOfGeneProducts")}
THIS = dict(genes=len(genes), rxn=n_rxn, gpr=n_gpr, transport=n_transport, reversible=n_reversible,
            comps=len(model.find(SB + "listOfCompartments")), met=len(species),
            unique=len({b for _, b in species.values()}))
print("iPiccel952", THIS)

MODELS = ["iPiccel952"] + list(PUBLISHED)
DATA = {"iPiccel952": THIS, **PUBLISHED}
for m in MODELS:
    d = DATA[m]
    d["gpr_pct"] = 100 * d["gpr"] / d["rxn"]
    d["tr_pct"] = 100 * d["transport"] / d["rxn"]
    d["rev_pct"] = 100 * d["reversible"] / d["rxn"]
    d["rpg"] = d["rxn"] / d["genes"]
    d["uniq_pct"] = 100 * d["unique"] / d["met"]

# Subsystems (SBML groups): reactions per subsystem and share with a gene association
SUBS = []
for g in model.iter(GRP + "group"):
    ids = [x.get(GRP + "idRef") for x in g.iter(GRP + "member") if x.get(GRP + "idRef") in has_gpr]
    SUBS.append((g.get(GRP + "name"), len(ids), sum(has_gpr[i] for i in ids)))
SUBS.sort(key=lambda t: -t[1])

# ---------------------------------------------------------------------------
# Photoautotrophic growth, active reactions and essentiality
# iPiccel952 values come from run_model.py and run_model2.py (cobrapy, model as shipped, GLPK).
# Values for the other models come from Soleymani et al. 2025 (text and Fig. 4). Items marked approximate were
# derived from reported fold changes or read from the bars of Fig. 4c.
# ---------------------------------------------------------------------------
EXTRA = {
    "iPiccel952": dict(growth=0.354, active=1827, ess_g=103, ess_r=1256),
    "iCO1515": dict(growth=1.711, active=715, ess_g=237, ess_r=670),
    "iCre1355": dict(growth=1.711 / 8, active=1350, ess_g=242, ess_r=437),
    "iCZ843": dict(growth=1.711 / 67, active=1280, ess_g=188, ess_r=388),
    "iAJ526": dict(growth=1.711 / 4, active=735, ess_g=134, ess_r=268),
}
# Same model with the 1,700 reactions of the Fatty acid metabolism group removed from the counts (fa_check.py, fa_genes.py).
NOFA = dict(active=510, ess_r=272, ess_g=103, n_rxn=1825, n_genes=931)
for m in MODELS:
    DATA[m].update(EXTRA[m])
    DATA[m]["active_pct"] = 100 * DATA[m]["active"] / DATA[m]["rxn"]
    DATA[m]["essg_pct"] = 100 * DATA[m]["ess_g"] / DATA[m]["genes"]
    DATA[m]["essr_pct"] = 100 * DATA[m]["ess_r"] / DATA[m]["rxn"]
DATA["iPiccel952 excl. FA"] = dict(active=NOFA["active"], ess_g=NOFA["ess_g"], ess_r=NOFA["ess_r"],
                                   active_pct=100 * NOFA["active"] / NOFA["n_rxn"],
                                   essg_pct=100 * NOFA["ess_g"] / NOFA["n_genes"],
                                   essr_pct=100 * NOFA["ess_r"] / NOFA["n_rxn"])

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
avail = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in avail), "DejaVu Sans")
plt.rcParams.update({
    "font.family": family, "font.size": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.3, "legend.frameon": False,
})
INK, INK2, MUTED, AXIS = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
BLUE, ORANGE, YELLOW, MAG, GREEN = "#2a78d6", "#eb6834", "#eda100", "#e87ba4", "#008300"
MODEL_COLOR = {"iPiccel952 excl. FA": "#8fd18f", "iPiccel952": GREEN, "iCO1515": BLUE, "iCre1355": ORANGE, "iCZ843": MAG, "iAJ526": YELLOW}

fig = plt.figure(figsize=(7.0, 5.5), dpi=100)
gs = fig.add_gridspec(3, 8, left=0.105, right=0.99, top=0.945, bottom=0.085, hspace=0.42, wspace=0.55)

OTHER = "#7fb0ea"
SHARED = "#bcd6f6"
TRACK = "#eef0f3"
BAND = "#eaf5ea"
PANELS = [
    ("A", "Genes", "genes", "{:,.0f}"),
    ("B", "Reactions", "rxn", "{:,.0f}"),
    ("C", "Metabolites", None, None),
    ("D", "Compartments", "comps", "{:.0f}"),
    ("E", "Reactions with GPR (%)", "gpr_pct", "{:.0f}"),
    ("F", "Transport (%)", "tr_pct", "{:.0f}"),
    ("G", "Reversible (%)", "rev_pct", "{:.0f}"),
    ("H", "Reactions per gene", "rpg", "{:.1f}"),
    ("I", "Photoautotrophic growth (1/h)", "growth", "{:.2f}"),
]
ypos = np.arange(len(MODELS))[::-1]
LW = 8.0


def rbar(ax, y, x0, x1, color, z=3):
    """Horizontal bar with rounded ends, drawn as a thick round-capped line."""
    ax.plot([x0, x1], [y, y], color=color, lw=LW, solid_capstyle="round", zorder=z)


axes = []
for k, (letter, title, key, fmt) in enumerate(PANELS):
    ax = fig.add_subplot(gs[k // 4, (k % 4) * 2:(k % 4) * 2 + 2]) if k < 8 else fig.add_subplot(gs[2, 0:2])
    ax.patch.set_alpha(0)
    ax.set_zorder(3)
    axes.append(ax)
    if key is None:
        mx = max(DATA[m]["met"] for m in MODELS)
        xmax = mx * 1.42
        vals = {m: DATA[m]["met"] for m in MODELS}
    else:
        mx = max(DATA[m][key] for m in MODELS)
        xmax = mx * 1.3
        vals = {m: DATA[m][key] for m in MODELS}
    pad = xmax * 0.025
    for y, m in zip(ypos, MODELS):
        foc = m == "iPiccel952"
        rbar(ax, y, pad, xmax * 0.97, TRACK, z=1)
        if key is None:
            d = DATA[m]
            rbar(ax, y, pad, d["met"], "#8fd18f" if foc else SHARED, z=2)
            rbar(ax, y, pad, max(d["unique"], pad * 1.01), GREEN if foc else BLUE, z=3)
        else:
            rbar(ax, y, pad, max(vals[m], pad * 1.01), GREEN if foc else OTHER, z=3)
        ax.text(vals[m] + xmax * 0.065, y, fmt.format(vals[m]) if key else f"{vals[m]:,}", va="center",
                fontsize=5.8, color=INK, fontweight="bold" if foc else "normal", zorder=5)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.55, len(MODELS) - 0.45)
    ax.set_yticks(ypos)
    if k % 4 == 0:
        ax.set_yticklabels(MODELS, fontsize=6.2)
        ax.get_yticklabels()[0].set_color(GREEN)
        ax.get_yticklabels()[0].set_fontweight("bold")
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=0, pad=3)
    ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0, 1.1, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom", color=INK)
    ax.text(0.12, 1.115, title, transform=ax.transAxes, fontsize=6.6, ha="left", va="bottom", color=INK2)

# J and K. Grouped vertical bars, one group per measure and one bar per model
def grouped(ax, letter, title, keys, labels, fmt, ymax, ylab_pad):
    ng = len(keys)
    w = 0.15
    order = list(MODELS)
    for gi, key in enumerate(keys):
        for mi, m in enumerate(order):
            x = gi + (mi - 2) * (w + 0.02)
            v = DATA[m][key]
            ax.plot([x, x], [ymax * 0.012, max(v, ymax * 0.02)], color=MODEL_COLOR[m], lw=6.2,
                    solid_capstyle="round", zorder=3)
            ax.text(x, v + ymax * 0.035, fmt.format(v), rotation=90, ha="center", va="bottom", fontsize=5.0,
                    color=INK, fontweight="bold" if m == "iPiccel952" else "normal")
    ax.set_xlim(-0.55, ng - 0.45)
    ax.set_ylim(0, ymax)
    ax.set_xticks(range(ng))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_yticks([])
    ax.tick_params(length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.plot([-0.5, ng - 0.5], [0, 0], color=AXIS, lw=0.6, zorder=1)
    ax.text(0, 1.1, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left", va="bottom", color=INK)
    ax.text(0.07, 1.115, title, transform=ax.transAxes, fontsize=6.6, ha="left", va="bottom", color=INK2)


axJ = fig.add_subplot(gs[2, 2:5])
grouped(axJ, "J", "Active reactions and essentiality (counts)",
        ["active", "ess_g", "ess_r"], ["Active reactions", "Essential genes", "Essential reactions"], "{:,.0f}",
        2400, 0)
axK = fig.add_subplot(gs[2, 5:8])
grouped(axK, "K", "Active reactions and essentiality (% of model)", ["active_pct", "essg_pct", "essr_pct"],
        ["Active reactions", "Essential genes", "Essential reactions"], "{:.0f}", 62, 0)

# Soft band behind the iPiccel952 row of the horizontal panels
from matplotlib.patches import Patch, Rectangle as Rect
for r, last in ((0, 3), (1, 7), (2, 8)):
    first = 4 * r if r < 2 else 8
    a0, a3 = axes[first], axes[last]
    p0, p3 = a0.get_position(), a3.get_position()
    rowh = p0.height / len(MODELS)
    y0 = p0.y1 - rowh
    x1 = 0.997 if r < 2 else p3.x1 + 0.012
    fig.patches.append(Rect((0.003, y0), x1 - 0.003, rowh, transform=fig.transFigure, facecolor=BAND,
                            edgecolor="none", zorder=0))

handles = [Patch(color="#555555", label="Unique metabolites (C, dark)"),
           Patch(color="#c9c9c9", label="Shared metabolites (C, light)")]
handles += [Patch(color=MODEL_COLOR[m], label=m) for m in MODELS]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, -0.006), ncol=7, fontsize=5.6,
           handlelength=1, handleheight=0.8, columnspacing=0.9, handletextpad=0.35)

fig.savefig(out_path, dpi=500, facecolor="white")
print("saved", out_path)