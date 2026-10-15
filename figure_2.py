"""
Fig. 2 for the iPiccel954 manuscript (structural comparison with other algal GEMs).

Usage
    python make_fig2.py iPiccel954.xml Fig2.png
or edit XML_PATH and OUT_PATH below and run the file without arguments.

Requires matplotlib and numpy. Values for iCO1515, iCre1355, iCZ843 and iAJ526 are taken from Table 1 of
Soleymani et al. 2025 (New Phytologist 248, 1572-1583). Values for iPiccel954 are computed from the SBML file.
"""
import sys
import json
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Patch

# Set these two paths when you run the script from an editor without command-line arguments.
# On Windows use a raw string, for example r"C:\Users\name\Figures\iPiccel954.xml".
XML_PATH = "iPiccel954.xml"
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
# iPiccel954 statistics from the SBML file
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
    if "biomass" in r.get("id").lower():
        continue
    if not reac and not prod:
        n_transport += 1          # Htu is kept with an empty stoichiometry and counted as transport, as in Fig. 1
    elif reac and prod and len({species[i][0] for i in reac + prod}) > 1:
        n_transport += 1

genes = {g.get(FBC + "id") for g in model.find(FBC + "listOfGeneProducts")}
THIS = dict(genes=len(genes), rxn=n_rxn, gpr=n_gpr, transport=n_transport, reversible=n_reversible,
            comps=len(model.find(SB + "listOfCompartments")), met=len(species),
            unique=len({b for _, b in species.values()}))
print("iPiccel954", THIS)

MODELS = ["iPiccel954"] + list(PUBLISHED)
DATA = {"iPiccel954": THIS, **PUBLISHED}
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
# iPiccel954 values come from cobrapy with GLPK on the model as shipped (photon uptake bound at -1000; active reactions by FVA at 100% of optimal growth; essential by single deletion below 1% of wild-type growth).
# Values for the other models come from Soleymani et al. 2025 (text and Fig. 4). Items marked approximate were
# derived from reported fold changes or read from the bars of Fig. 4c.
# ---------------------------------------------------------------------------
EXTRA = {
    "iPiccel954": dict(growth=0.337, active=1860, ess_g=116, ess_r=1258),
    "iCO1515": dict(growth=1.711, active=715, ess_g=237, ess_r=670),
    "iCre1355": dict(growth=1.711 / 8, active=1350, ess_g=242, ess_r=437),
    "iCZ843": dict(growth=1.711 / 67, active=1280, ess_g=188, ess_r=388),
    "iAJ526": dict(growth=1.711 / 4, active=735, ess_g=134, ess_r=268),
}
# Same model with the 1,700 reactions of the Fatty acid metabolism group removed from the counts (fa_check.py, fa_genes.py).
NOFA = dict(active=531, ess_r=274, ess_g=116, n_rxn=1825, n_genes=931)
for m in MODELS:
    DATA[m].update(EXTRA[m])
    DATA[m]["active_pct"] = 100 * DATA[m]["active"] / DATA[m]["rxn"]
    DATA[m]["essg_pct"] = 100 * DATA[m]["ess_g"] / DATA[m]["genes"]
    DATA[m]["essr_pct"] = 100 * DATA[m]["ess_r"] / DATA[m]["rxn"]
DATA["iPiccel954 excl. FA"] = dict(active=NOFA["active"], ess_g=NOFA["ess_g"], ess_r=NOFA["ess_r"],
                                   active_pct=100 * NOFA["active"] / NOFA["n_rxn"],
                                   essg_pct=100 * NOFA["ess_g"] / NOFA["n_genes"],
                                   essr_pct=100 * NOFA["ess_r"] / NOFA["n_rxn"])

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
avail = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in avail), "DejaVu Sans")
plt.rcParams.update({
    "font.family": family, "font.size": 9.5, "axes.linewidth": 0.8,
    "text.color": "#000000", "axes.labelcolor": "#000000", "axes.edgecolor": "#000000",
    "xtick.color": "#000000", "ytick.color": "#000000",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5, "legend.frameon": False,
})
INK, INK2, MUTED, AXIS = "#000000", "#000000", "#898781", "#000000"   # all text black
BLUE, ORANGE, YELLOW, MAG, GREEN = "#2a78d6", "#eb6834", "#eda100", "#e87ba4", "#008300"
MODEL_COLOR = {"iPiccel954 excl. FA": "#8fd18f", "iPiccel954": GREEN, "iCO1515": BLUE, "iCre1355": ORANGE, "iCZ843": MAG, "iAJ526": YELLOW}

fig = plt.figure(figsize=(11.5, 13.0), dpi=100)
gs = fig.add_gridspec(5, 8, left=0.125, right=0.99, top=0.96, bottom=0.075, hspace=0.70, wspace=0.50)

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
LW = 12.0


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
        foc = m == "iPiccel954"
        rbar(ax, y, pad, xmax * 0.97, TRACK, z=1)
        if key is None:
            d = DATA[m]
            rbar(ax, y, pad, d["met"], "#8fd18f" if foc else SHARED, z=2)
            rbar(ax, y, pad, max(d["unique"], pad * 1.01), GREEN if foc else BLUE, z=3)
        else:
            rbar(ax, y, pad, max(vals[m], pad * 1.01), GREEN if foc else OTHER, z=3)
        ax.text(vals[m] + xmax * 0.065, y, fmt.format(vals[m]) if key else f"{vals[m]:,}", va="center",
                fontsize=9.5, color=INK, fontweight="bold" if foc else "normal", zorder=5)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.55, len(MODELS) - 0.45)
    ax.set_yticks(ypos)
    if k % 4 == 0:
        ax.set_yticklabels(MODELS, fontsize=10)
        ax.get_yticklabels()[0].set_color(GREEN)
        ax.get_yticklabels()[0].set_fontweight("bold")
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=0, pad=3)
    ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0, 1.1, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", ha="left", va="bottom", color=INK)
    ax.text(0.12, 1.115, title, transform=ax.transAxes, fontsize=10.5, ha="left", va="bottom", color=INK2)

# J and K. Grouped vertical bars, one group per measure and one bar per model
def grouped(ax, letter, title, keys, labels, fmt, ymax, ylab_pad):
    ng = len(keys)
    w = 0.15
    order = list(MODELS)
    for gi, key in enumerate(keys):
        for mi, m in enumerate(order):
            x = gi + (mi - 2) * (w + 0.02)
            v = DATA[m][key]
            ax.plot([x, x], [ymax * 0.012, max(v, ymax * 0.02)], color=MODEL_COLOR[m], lw=9.5,
                    solid_capstyle="round", zorder=3)
            ax.text(x, v + ymax * 0.035, fmt.format(v), rotation=90, ha="center", va="bottom", fontsize=8.5,
                    color=INK, fontweight="bold" if m == "iPiccel954" else "normal")
    ax.set_xlim(-0.55, ng - 0.45)
    ax.set_ylim(0, ymax)
    ax.set_xticks(range(ng))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_yticks([])
    ax.tick_params(length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.plot([-0.5, ng - 0.5], [0, 0], color=AXIS, lw=0.6, zorder=1)
    ax.text(0, 1.1, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", ha="left", va="bottom", color=INK)
    ax.text(0.07, 1.115, title, transform=ax.transAxes, fontsize=10.5, ha="left", va="bottom", color=INK2)


axJ = fig.add_subplot(gs[2, 2:5])
grouped(axJ, "J", "Active reactions and essentiality (counts)",
        ["active", "ess_g", "ess_r"], ["Active\nreactions", "Essential\ngenes", "Essential\nreactions"], "{:,.0f}",
        2400, 0)
axK = fig.add_subplot(gs[2, 5:8])
grouped(axK, "K", "Active reactions and essentiality (% of model)", ["active_pct", "essg_pct", "essr_pct"],
        ["Active\nreactions", "Essential\ngenes", "Essential\nreactions"], "{:.0f}", 62, 0)


# L and M. Validation against published data (values from val_compute.py to val_compute4.py, saved in val_cache.json)
VAL = json.load(open("val_cache.json"))
MEAS = "#9a9a9a"        # measured values
LIGHTG = "#8fd18f"      # model at the default photon bound
FLUX = [("rubisco_c", "RuBisCO\ncarboxylation"), ("rubisco_o", "RuBisCO\noxygenation"), ("pga_red", "PGA\nreduction"),
        ("sbp", "SBPase"), ("starch", "Starch\nsynthesis"), ("pdh_h", "PDH\nplastid"), ("enol", "Enolase"),
        ("pepc", "PEP\ncarboxylase"), ("pdh_m", "PDH\nmitochondrion"), ("cs_m", "Citrate\nsynthase"), ("me_m", "Malic\nenzyme")]
axL = fig.add_subplot(gs[3, 0:8])
FLOOR = 0.1
wL = 0.26
for gi, (k, lab) in enumerate(FLUX):
    meas = VAL["flux_meas"][k]
    lo, hi = VAL["fva"][k]
    mod = VAL["flux_model"][k]
    inside = lo - 1e-6 <= meas <= hi + 1e-6
    for x, v, col in ((gi - wL / 1.6, meas, MEAS), (gi + wL / 1.6, mod, GREEN)):
        axL.plot([x, x], [FLOOR, max(v, FLOOR * 1.3)], color=col, lw=11, solid_capstyle="round", zorder=3)
    xm = gi + wL / 1.6
    axL.plot([xm, xm], [max(lo, FLOOR), hi], color="#000000", lw=1.3, zorder=4)
    for yy in (max(lo, FLOOR), hi):
        axL.plot([xm - 0.06, xm + 0.06], [yy, yy], color="#000000", lw=1.3, zorder=4)
    axL.text(gi - wL / 1.6, max(meas, FLOOR * 1.3) * 1.25, f"{meas:.3g}", ha="center", va="bottom", fontsize=8.3, rotation=90, color=INK)
axL.set_yscale("log"); axL.set_ylim(FLOOR, 1500); axL.set_xlim(-0.6, len(FLUX) - 0.4)
axL.set_xticks(range(len(FLUX))); axL.set_xticklabels([f[1] for f in FLUX], fontsize=9.5, linespacing=0.95)
axL.set_yticks([]); axL.minorticks_off(); axL.tick_params(length=0, pad=2)
for sp in axL.spines.values():
    sp.set_visible(False)
axL.plot([-0.55, len(FLUX) - 0.45], [FLOOR, FLOOR], color=AXIS, lw=0.6, zorder=1)
axL.text(0, 1.1, "L", transform=axL.transAxes, fontsize=14, fontweight="bold", ha="left", va="bottom", color=INK)
axL.text(0.035, 1.115, "Fluxes per 100 net CO$_2$ uptake, model FVA range and measured (log scale)", transform=axL.transAxes, fontsize=10.5, ha="left", va="bottom", color=INK2)
axL.legend([Patch(color=MEAS), Patch(color=GREEN), plt.Line2D([], [], color="#000000", lw=1.3)],
           ["Steichen et al., 2024", "iPiccel954", "FVA range"], loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=9.5, handlelength=1, columnspacing=1.0, handletextpad=0.4)

axM = fig.add_subplot(gs[4, 0:5])
from_cache = VAL
MU_MAX = VAL["mu_max"]
CURVE = np.array(VAL["curve"]); CFV = VAL["cfrac"]
XG = np.array([60 * 0.34 * 3.6, 2000 * 0.22 * 3.6]) * CFV
MUMOD = np.interp(XG, CURVE[:, 0], CURVE[:, 1])
GROW = [("Weissman et al.,\n2018, 60", 0.06, 0.008, MUMOD[0], True),
        ("Weissman et al.,\n2018, 2000", 0.26, 0.015, MUMOD[1], True),
        ("Weissman et al.,\n2018, O$_2$", 0.30, 0.03, MU_MAX, False),
        ("Steichen et al.,\n2024", 0.323, 0.043, MU_MAX, False),
        ("Torres-Tiji et al.,\n2026", 0.289, 0.058, MU_MAX, False)]
wM = 0.24
for gi, (lab, v, e, mod, lit) in enumerate(GROW):
    xs = (gi - wM / 1.6, gi + wM / 1.6)
    axM.plot([xs[0]] * 2, [0.003, v], color=MEAS, lw=9.5, solid_capstyle="round", zorder=3)
    axM.plot([xs[0]] * 2, [v - e, v + e], color="#000000", lw=1.2, zorder=4)
    axM.plot([xs[1]] * 2, [0.003, mod], color=GREEN if lit else LIGHTG, lw=9.5, solid_capstyle="round", zorder=3)
    axM.text(xs[0] - 0.02, v + e + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=8.5, rotation=0, color=INK)
    axM.text(xs[1], mod + 0.012, f"{mod:.2f}" if mod >= 0.1 else f"{mod:.3f}", ha="center", va="bottom", fontsize=8.5, rotation=0, color=INK,
             fontweight="bold" if lit else "normal")
axM.set_xlim(-0.6, len(GROW) - 0.4); axM.set_ylim(0, 0.42)
axM.set_xticks(range(len(GROW))); axM.set_xticklabels([g[0] for g in GROW], fontsize=9.5, linespacing=0.95)
axM.set_yticks([]); axM.tick_params(length=0, pad=2)
for sp in axM.spines.values():
    sp.set_visible(False)
axM.plot([-0.55, len(GROW) - 0.45], [0, 0], color=AXIS, lw=0.6, zorder=1)
axM.text(0, 1.1, "M", transform=axM.transAxes, fontsize=14, fontweight="bold", ha="left", va="bottom", color=INK)
axM.text(0.06, 1.115, "Growth rate (1/h)", transform=axM.transAxes, fontsize=10.5, ha="left", va="bottom", color=INK2)
axM.legend([Patch(color=MEAS), Patch(color=GREEN), Patch(color=LIGHTG)], ["Measured", "Model at study light", "Model at default photon bound"],
           loc="center left", bbox_to_anchor=(1.02, 0.55), fontsize=9.5, handlelength=1, handletextpad=0.4, labelspacing=0.4)

# Soft band behind the iPiccel954 row of the horizontal panels
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
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, -0.004), ncol=4, fontsize=9.5,
           handlelength=1, handleheight=0.8, columnspacing=0.9, handletextpad=0.35)

fig.savefig(out_path, dpi=500, facecolor="white")
print("saved", out_path)