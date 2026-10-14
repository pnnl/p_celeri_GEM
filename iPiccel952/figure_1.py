"""
Fig. 1 for the iPiccel952 manuscript (reconstruction and curation).

Usage
    python make_fig1.py iPiccel952.xml Fig1.png
or edit XML_PATH and OUT_PATH below and run the file without arguments.

Requires matplotlib only. Counts from the final SBML model are computed here. Numbers from
the curation report and the MEMOTE report are set as constants in the block marked REPORT VALUES.
"""
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# Set these two paths when you run the script from an editor without command-line arguments.
# On Windows use a raw string, for example r"C:\Users\name\Figures\iPiccel952.xml".
XML_PATH = "iPiccel952.xml"
OUT_PATH = "Fig1.png"
xml_path = sys.argv[1] if len(sys.argv) > 1 else XML_PATH
out_path = sys.argv[2] if len(sys.argv) > 2 else OUT_PATH

# ---------------------------------------------------------------------------
# REPORT VALUES (edit here if the reports change)
# ---------------------------------------------------------------------------
DRAFT_COUNTS = dict(rxn=4344, met=3235, genes=1084)          # curation report
FINAL_GENES = 954                                              # MEMOTE report
BIOMASS_BEFORE = [9.5, 22.2, 57.4, 8.1, 0.4, 2.5]              # % non-maintenance mass, curation report
BIOMASS_AFTER = [60.6, 14.4, 13.0, 6.1, 0.7, 5.5]
BIOMASS_CLASSES = ["Protein", "Carbohydrate", "Lipid", "RNA", "DNA", "Pigments"]
UNBALANCED_STAGES = [613, 127]                                 # draft, after hydrogen correction
CYCLE_REACTIONS = [1437, 1246, 2]                              # draft, after cycle removal, final (MEMOTE)
MEMOTE = [("Total", 71.5), ("Consistency", 56.9), ("Metabolite\nannotation", 71.2),
          ("Reaction\nannotation", 70.1), ("Gene\nannotation", 33.3), ("SBO terms", 90.9)]

# ---------------------------------------------------------------------------
# Statistics from the SBML file
# ---------------------------------------------------------------------------
tree = ET.parse(xml_path)
root = tree.getroot()
SB = "{http://www.sbml.org/sbml/level3/version1/core}"
FBC = None
for _, (prefix, uri) in ET.iterparse(xml_path, events=["start-ns"]):
    if prefix == "fbc":
        FBC = "{" + uri + "}"
model = root.find(SB + "model")

species = {}
for s in model.find(SB + "listOfSpecies"):
    species[s.get("id")] = dict(comp=s.get("compartment"), formula=s.get(FBC + "chemicalFormula"))

reactions = []
for r in model.find(SB + "listOfReactions"):
    lr, lp = r.find(SB + "listOfReactants"), r.find(SB + "listOfProducts")
    reac = [(x.get("species"), float(x.get("stoichiometry"))) for x in (lr if lr is not None else [])]
    prod = [(x.get("species"), float(x.get("stoichiometry"))) for x in (lp if lp is not None else [])]
    reactions.append(dict(id=r.get("id"), reac=reac, prod=prod))


def parse_formula(f):
    d = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", f or ""):
        d[el] += float(n) if n else 1.0
    return d


def is_boundary(r):
    return len(r["reac"]) == 0 or len(r["prod"]) == 0


categories = Counter()
metabolic_by_comp = Counter()
unbalanced = []
for r in reactions:
    if "biomass" in r["id"].lower():
        categories["biomass"] += 1
        continue
    if is_boundary(r):
        categories["exchange"] += 1
        continue
    comps = {species[s]["comp"] for s, _ in r["reac"] + r["prod"]}
    if len(comps) > 1:
        categories["transport"] += 1
    else:
        categories["metabolic"] += 1
        metabolic_by_comp[next(iter(comps))] += 1
    tot = Counter()
    ok = True
    for sign, side in ((-1, r["reac"]), (1, r["prod"])):
        for s, n in side:
            f = species[s]["formula"]
            ok = ok and bool(f)
            for el, k in parse_formula(f).items():
                tot[el] += sign * n * k
    if not (ok and all(abs(v) < 1e-6 for v in tot.values())):
        unbalanced.append(r)


def imbalance_kind(r):
    ids = [s for s, _ in r["reac"] + r["prod"]]
    forms = [species[s]["formula"] or "" for s in ids]
    if any("Z" in f for f in forms):
        return "Photon placeholder"
    if r["id"].upper().replace("R_", "", 1).startswith("STARCH300") or any("starch300" in s for s in ids):
        return "Starch"
    return "Ceramide, R group, other"


kinds = Counter(imbalance_kind(r) for r in unbalanced)
n_rxn, n_met, n_unbalanced = len(reactions), len(species), len(unbalanced)
print("reactions", n_rxn, "metabolites", n_met, "unbalanced", n_unbalanced, dict(kinds))

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
    "legend.fontsize": 6.5, "legend.frameon": False,
})

INK, INK2, MUTED, AXIS = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW, MAG, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"
LBLUE, LGRAY = "#cde2fb", "#ecebe6"

fig = plt.figure(figsize=(7.2, 6.6), dpi=100)
gsA = fig.add_gridspec(1, 1, left=0.075, right=0.985, top=0.975, bottom=0.765)
gs = fig.add_gridspec(2, 3, hspace=0.58, wspace=0.55, left=0.075, right=0.985, top=0.70, bottom=0.085)


def panel_label(ax, letter, dx=-0.16, dy=1.08):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left", color=INK)


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2)


# ---------------------------------------------------------------------------
# Panel A. Reconstruction pipeline
# ---------------------------------------------------------------------------
axA = fig.add_subplot(gsA[0, 0])
axA.set_xlim(0, 100)
axA.set_ylim(0.5, 42.5)
axA.axis("off")
panel_label(axA, "A", dx=-0.045, dy=0.97)


def box(cx, cy, w, h, title, sub=None, fc=LGRAY, tc=INK, ec=AXIS, bold=False):
    axA.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=1.2",
                                 fc=fc, ec=ec, lw=0.7))
    weight = "bold" if bold else "normal"
    if sub:
        two = "\n" in sub
        axA.text(cx, cy + (h * 0.24 if two else 1.3), title, ha="center", va="center", fontsize=6.6,
                 color=tc, fontweight=weight, linespacing=1.15)
        axA.text(cx, cy - (h * 0.17 if two else 2.6), sub, ha="center", va="center", fontsize=5.8,
                 color=tc, linespacing=1.15)
    else:
        axA.text(cx, cy, title, ha="center", va="center", fontsize=6.6, color=tc,
                 fontweight=weight, linespacing=1.15)


def arrow(p0, p1, color=INK2):
    axA.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=6, lw=0.7, color=color,
                                  shrinkA=0, shrinkB=0))


H = 9.0
box(7, 35, 12.5, H, "iCre1355 GEM\n(literature)")
box(7, 9, 12.5, H, "P. celeri\ngenome")
box(21.5, 22, 9.5, H, "RAPS", fc=LBLUE, ec=BLUE)
box(36, 22, 12, H, "Draft iPcel\nmodel", fc=LBLUE, ec=BLUE)
box(51, 22, 13, H, "Homology to\niPren GEM", fc=LBLUE, ec=BLUE)
box(51, 37, 13, H - 1, "iPren GEM\n(literature)")
box(65, 35, 11, H, "Common to\niPcel, iPren")
box(65, 22, 11, H, "Unique to\niPcel")
box(65, 9, 11, H, "Unique to\niPren")
box(83.5, 22, 12, 12, "Final iPcel GEM", "4344 rxn\n1084 genes", fc=BLUE, tc="white", ec=BLUE, bold=True)
box(95, 22, 9.5, 12, "Curation", "7 stages", fc=LBLUE, ec=BLUE)
box(95, 6.5, 9.5, 11, "iPiccel952", "3525 rxn\n954 genes", fc=GREEN, tc="white", ec=GREEN, bold=True)

# The RAPS box spans x 16.75 to 26.25, so the two input arrows stop just outside its left edge.
arrow((13.6, 31.2), (16.45, 25.0))
arrow((13.6, 12.8), (16.45, 19.0))
arrow((26.5, 22), (29.8, 22))
arrow((42.2, 22), (44.3, 22))
arrow((51, 32.5), (51, 26.6))
arrow((57.6, 25), (59.3, 31.5))
arrow((57.6, 22), (59.3, 22))
arrow((57.6, 19), (59.3, 12.5))
arrow((70.7, 31.5), (77.6, 26.9))
arrow((70.7, 22), (77.4, 22))
arrow((70.7, 12.5), (77.6, 17.1))
arrow((89.6, 22), (90.2, 22))
arrow((95, 16.0), (95, 12.2))
axA.text(73.4, 32.4, "all added", fontsize=5.8, color=INK2, ha="left", style="italic")
axA.text(74.0, 23.4, "all added", fontsize=5.8, color=INK2, ha="center", style="italic")
axA.text(73.8, 10.4, "added only\nas needed", fontsize=5.8, color=INK2, ha="left", va="center",
         style="italic", linespacing=1.1)

# ---------------------------------------------------------------------------
# Panel B. Biomass composition
# ---------------------------------------------------------------------------
axB = fig.add_subplot(gs[0, 0])
panel_label(axB, "B")
colors = [BLUE, ORANGE, AQUA, YELLOW, MAG, GREEN]
ys = [1.55, 0]
for yi, vals in zip(ys, (BIOMASS_BEFORE, BIOMASS_AFTER)):
    left = 0
    for v, c in zip(vals, colors):
        axB.barh(yi, v, left=left, color=c, edgecolor="white", linewidth=0.8, height=0.62)
        if v >= 8:
            axB.text(left + v / 2, yi, f"{v:.1f}", ha="center", va="center", fontsize=6.3, fontweight="bold",
                     color="white" if c in (BLUE, GREEN, AQUA, ORANGE) else INK)
        left += v
axB.set_yticks(ys)
axB.set_yticklabels(["Draft", "Curated"])
axB.set_xlim(0, 100)
axB.set_ylim(-0.45, 2.0)
axB.set_xlabel("Share of biomass mass excluding maintenance (%)")
style(axB)
axB.spines["left"].set_visible(False)
axB.tick_params(axis="y", length=0)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
axB.legend(handles, BIOMASS_CLASSES, ncol=3, loc="center", bbox_to_anchor=(0.5, 0.47),
           columnspacing=1.0, handlelength=0.9, handletextpad=0.4)

# ---------------------------------------------------------------------------
# Panel C. Mass-imbalanced reactions
# ---------------------------------------------------------------------------
axC = fig.add_subplot(gs[0, 1])
panel_label(axC, "C")
vals = UNBALANCED_STAGES + [n_unbalanced]
axC.bar(range(3), vals, color=[ORANGE, ORANGE, BLUE], width=0.58, edgecolor="white", linewidth=0.8)
for i, v in enumerate(vals):
    axC.text(i, v + 15, str(v), ha="center", va="bottom", fontsize=6.6, color=INK, fontweight="bold")
axC.set_xticks(range(3))
axC.set_xticklabels(["Draft", "After H\ncorrection", "Final"])
axC.set_ylabel("Mass-imbalanced reactions")
axC.set_ylim(0, 720)
axC.text(2.38, 330,
         f"Final model\n{kinds.get('Photon placeholder', 0)} photon placeholder\n"
         f"{kinds.get('Ceramide, R group, other', 0)} ceramide, R group, other\n{kinds.get('Starch', 0)} starch",
         fontsize=6.0, color=INK2, ha="right", va="center", linespacing=1.25)
style(axC)

# ---------------------------------------------------------------------------
# Panel D. Balanced cycles
# ---------------------------------------------------------------------------
axD = fig.add_subplot(gs[0, 2])
panel_label(axD, "D")
axD.bar(range(3), CYCLE_REACTIONS, color=[ORANGE, ORANGE, BLUE], width=0.58, edgecolor="white", linewidth=0.8)
for i, v in enumerate(CYCLE_REACTIONS):
    axD.text(i, v + 25, f"{v:,}", ha="center", va="bottom", fontsize=6.6, color=INK, fontweight="bold")
axD.set_xticks(range(3))
axD.set_xticklabels(["Draft", "After cycle\nremoval", "Final"])
axD.set_ylabel("Reactions in balanced cycles")
axD.set_ylim(0, 1650)
axD.text(2.45, 1330, "4 energy-generating\ncycles removed\nUnbounded flux\n36.6% to 0.27%", fontsize=6.0,
         color=INK2, ha="right", va="center", linespacing=1.25)
style(axD)

# ---------------------------------------------------------------------------
# Panel E. Model size
# ---------------------------------------------------------------------------
axE = fig.add_subplot(gs[1, 0])
panel_label(axE, "E")
draft = [DRAFT_COUNTS["rxn"], DRAFT_COUNTS["met"], DRAFT_COUNTS["genes"]]
final = [n_rxn, n_met, FINAL_GENES]
x = range(3)
w = 0.36
axE.bar([i - w / 2 for i in x], draft, w, color=ORANGE, edgecolor="white", linewidth=0.8, label="Draft")
axE.bar([i + w / 2 for i in x], final, w, color=BLUE, edgecolor="white", linewidth=0.8, label="Curated")
for i in x:
    for xpos, val in ((i - w / 2, draft[i]), (i + w / 2, final[i])):
        axE.text(xpos, val + 70, f"{val:,}", ha="center", va="bottom", rotation=90, fontsize=6.0, color=INK2)
axE.set_xticks(list(x))
axE.set_xticklabels(["Reactions", "Metabolites", "Genes"])
axE.set_ylabel("Count")
axE.set_ylim(0, 5600)
axE.legend(loc="upper right", handlelength=0.9, bbox_to_anchor=(1.02, 1.04))
style(axE)

# ---------------------------------------------------------------------------
# Panel F. MEMOTE scores
# ---------------------------------------------------------------------------
axF = fig.add_subplot(gs[1, 1])
panel_label(axF, "F")
names = [m[0] for m in MEMOTE]
scores = [m[1] for m in MEMOTE]
ypos = list(range(len(names)))[::-1]
axF.barh(ypos, scores, color=[INK2] + [BLUE] * (len(names) - 1), height=0.62, edgecolor="white", linewidth=0.8)
for y, s in zip(ypos, scores):
    axF.text(s + 1.5, y, f"{s:.1f}", va="center", fontsize=6.3, color=INK)
axF.set_yticks(ypos)
axF.set_yticklabels(names)
axF.set_xlim(0, 108)
axF.set_xlabel("MEMOTE score (%)")
style(axF)
axF.spines["left"].set_visible(False)
axF.tick_params(axis="y", length=0)

# ---------------------------------------------------------------------------
# Panel G. Reactions by compartment and type
# ---------------------------------------------------------------------------
axG = fig.add_subplot(gs[1, 2])
panel_label(axG, "G")
cname = {"c": "Cytosol", "h": "Chloroplast", "m": "Mitochondrion", "x": "Peroxisome", "u": "Thylakoid",
         "r": "ER", "n": "Nucleus", "e": "Extracellular", "g": "Golgi"}
items = [(cname.get(k, k), v, BLUE) for k, v in metabolic_by_comp.items()]
items.append(("Transport", categories["transport"], AQUA))
items.append(("Exchange, demand", categories["exchange"], MUTED))
items.sort(key=lambda t: t[1])
axG.barh(range(len(items)), [t[1] for t in items], color=[t[2] for t in items], height=0.62,
         edgecolor="white", linewidth=0.8)
for i, t in enumerate(items):
    axG.text(t[1] + 25, i, f"{t[1]:,}", va="center", fontsize=6.3, color=INK)
axG.set_yticks(range(len(items)))
axG.set_yticklabels([t[0] for t in items])
axG.set_xlim(0, 1750)
axG.set_xlabel("Reactions")
style(axG)
axG.spines["left"].set_visible(False)
axG.tick_params(axis="y", length=0)

fig.savefig(out_path, dpi=500, facecolor="white")
print("saved", out_path)