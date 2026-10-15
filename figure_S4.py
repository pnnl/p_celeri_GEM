"""
Supplementary Fig. S4 for the eciPiccel954 manuscript (supports Fig. 4, sensitivity of the nitrogen and carbon coupling).

Usage
    python figure_S4.py ecGEMs.zip protein_abundance_g_per_gCDW.xlsx iPceleri.fasta iPiccel954.xml FigS4.png [coupling_sens.json]
    ecGEMs.zip                          zip file with the WT_Air and Wt_CO2 folders
    protein_abundance_g_per_gCDW.xlsx   measured proteome (g/gCDW), columns Protein, Wt_Air, Wt_CO2
    iPceleri.fasta                      protein sequences (molecular weights of urea carboxylase and allophanate hydrolase)
    iPiccel954.xml                      original SBML model, used only for reaction names
    coupling_sens.json                  cache of the sensitivity runs (default coupling_sens.json in the current folder).
                                        If the file exists it is reused, delete it to recompute.
Requires cobra, matplotlib, numpy, pandas and openpyxl. The first run solves a few hundred LPs with GLPK and takes about a minute.

The coupling model is the one of Fig. 4. Only the nitrogen exchanges of urea and valine are open, the shortcuts uahCv and
ATAH_REV are blocked, glutamate synthase uses NADH, there is no protein pool, and the bounds of urea carboxylase (URCB)
and allophanate hydrolase (ALPH) are kcat times measured abundance over molecular weight.

Panels
    A to C  growth when the capacity of URCB, of ALPH or of both is changed relative to the measured value
    D       model assumptions (NADH or NADPH glutamate synthase, original bound file)
    E       growth fold WT-CO2 over WT-Air for several capacities, against the experiment
    F       reactions whose ten-fold bound relaxation raises growth by more than 1%
"""
import json
import os
import re
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
from matplotlib import font_manager
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")
args = sys.argv[1:]
if len(args) < 5:
    sys.exit(__doc__)
zip_path, prot_path, fasta_path, model_path, out_path = args[:5]
cache_path = args[5] if len(args) > 5 else "coupling_sens.json"

BIO = "Biomass_Picre_Auto"
COND = {"air": ("WT_Air", "Wt_Air", "WT-Air"), "co2": ("Wt_CO2", "Wt_CO2", "WT-CO2")}
EXP_GROWTH = {"air": 0.028, "co2": 0.324941}
ANNOT = {"URCB": "GENE4066869", "ALPH": "GENE4066871"}
MULT = [0.5, 1, 2, 5, 10, None]

tmp = tempfile.mkdtemp(prefix="ecgems_")
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(tmp)
root = tmp
for r_, dirs, _ in os.walk(tmp):
    if "WT_Air" in dirs:
        root = r_
        break

base = lambda r: re.sub(r"_(FWD|REV)$", "", r)
orig = cobra.io.read_sbml_model(model_path)
names = {r.id: r.name for r in orig.reactions}


def rname(r, n=30):
    s_ = names.get(r) or names.get(base(r)) or r
    s_ = s_.replace("_", " ")
    return s_ if len(s_) <= n else s_[: n - 1] + "."


# ---------------------------------------------------------------------------
# Capacities of urea carboxylase and allophanate hydrolase, kcat x abundance / MW
# ---------------------------------------------------------------------------
seqs, cur = {}, None
for line in open(fasta_path):
    line = line.strip()
    if line.startswith(">"):
        cur = line[1:].split(" ", 1)[0].split("_")[0]
        seqs[cur] = ""
    elif cur:
        seqs[cur] += line
AA = {"A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388, "E": 129.1155, "Q": 128.1307,
      "G": 57.0519, "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766,
      "P": 97.1167, "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326}
pa = pd.read_excel(prot_path).iloc[:, :3].set_index("Protein")
mw = lambda g: sum(AA.get(c, 110.0) for c in seqs[g]) + 18.02
CAP = {}
for cond, (folder, tag, _) in COND.items():
    kc = pd.read_excel(f"{root}/{folder}/relaxation_algorithm/kcat_mw_{tag}.xlsx", sheet_name="reaction_protein_cost")
    CAP[cond] = {rid: float(kc.loc[kc.Reaction == rid, "Kcat"].iloc[0]) * pa.loc[g, tag] / mw(g) * 1000.0
                 for rid, g in ANNOT.items()}
    print(cond, "capacity (mmol/gDW/h)", {k: round(float(v), 4) for k, v in CAP[cond].items()}, flush=True)


# ---------------------------------------------------------------------------
# Coupling model and sensitivity runs
# ---------------------------------------------------------------------------
def coupling_model(cond, nadh=True, orig_ub=False):
    folder, tag, _ = COND[cond]
    m = cobra.io.read_sbml_model(f"{root}/{folder}/irrev_iPiccel954.xml")
    m.solver = "glpk"
    ubf = f"UB_{tag}.xlsx" if orig_ub else f"UB_{tag}_revised.xlsx"
    ub = pd.read_excel(f"{root}/{folder}/relaxation_algorithm/{ubf}", header=None, usecols=[0, 1], names=["r", "u"])
    ids = {r.id for r in m.reactions}
    for r, u in zip(ub.r, ub.u):
        u = pd.to_numeric(u, errors="coerce")
        if str(r) in ids and pd.notna(u):
            m.reactions.get_by_id(str(r)).upper_bound = max(0.0, float(u))
    for rid, u in CAP[cond].items():
        m.reactions.get_by_id(rid).upper_bound = u
    keep = {"EX_urea(e)_REV", "EX_val-L(e)"}
    for r in m.exchanges:
        if r.id not in keep and any("N" in a.elements for a in r.metabolites):
            r.upper_bound = 0.0
    for rid in ("uahCv", "ATAH_REV"):
        m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    nadph, nadp, nadh_, nad = [m.metabolites.get_by_id(x) for x in ("nadph_h", "nadp_h", "nadh_h", "nad_h")]
    if nadh and "GLUS(nadph)" in m.reactions:      # older model versions, make GOGAT NADH dependent
        m.reactions.get_by_id("GLUS(nadph)").add_metabolites({nadph: 1, nadp: -1, nadh_: -1, nad: 1})
    if not nadh and "GLUS(nadh)" in m.reactions:   # final model, test an NADPH dependent GOGAT
        m.reactions.get_by_id("GLUS(nadh)").add_metabolites({nadh_: 1, nad: -1, nadph: -1, nadp: 1})
    m.objective = BIO
    return m


def coupling_sensitivity():
    if os.path.exists(cache_path):
        print("reusing", cache_path, flush=True)
        return json.load(open(cache_path))
    res = {}
    for cond in COND:
        m = coupling_model(cond)
        r_ = {}
        for key, rids in (("URCB", ["URCB"]), ("ALPH", ["ALPH"]), ("both", ["URCB", "ALPH"])):
            vals = []
            for mu in MULT:
                with m:
                    for rid in rids:
                        m.reactions.get_by_id(rid).upper_bound = 1000.0 if mu is None else CAP[cond][rid] * mu
                    vals.append(m.optimize().objective_value)
            r_[key] = vals
        r_["nadh"] = m.optimize().objective_value
        r_["nadph"] = coupling_model(cond, nadh=False).optimize().objective_value
        r_["orig_ub"] = coupling_model(cond, orig_ub=True).optimize().objective_value
        s0 = m.optimize()
        g0, fl = s0.objective_value, s0.fluxes
        gains = {}
        for rx in m.reactions:
            if 0 < rx.upper_bound < 1000 and abs(fl[rx.id]) > 1e-9 and rx.id != BIO:
                with m:
                    rx.upper_bound = min(1000.0, rx.upper_bound * 10)
                    g = m.optimize().objective_value
                gains[rx.id] = (g - g0) / g0
        r_["gains"] = gains
        r_["g0"] = g0
        res[cond] = r_
        print(f"coupling sensitivity {cond} done, baseline growth {g0:.5f}", flush=True)
    json.dump(res, open(cache_path, "w"))
    return res


CS = coupling_sensitivity()

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
BLUE, ORANGE = "#2a78d6", "#eb6834"


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


def legend(ax, entries, **kw):
    hs = [tuple(dot(c) for c in cols) for cols, _ in entries]
    ax.legend(hs, [t for _, t in entries], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
              handletextpad=0.4, labelspacing=0.35, borderpad=0, fontsize=9.5, **kw)


def hbar(ax, y, v, col, lw=8.0, b=0.0):
    ax.plot([b, v], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=3)


def vbar(ax, x, v, col, lw=11.0, b=0.0):
    ax.plot([x, x], [b, v], color=col, lw=lw, solid_capstyle="round", zorder=3)


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


fig = plt.figure(figsize=(11.9, 9.6), dpi=100, layout="constrained")
fig.get_layout_engine().set(w_pad=0.04, h_pad=0.06, wspace=0.12, hspace=0.14, rect=(0.0, 0.0, 1.0, 0.955))
gs = fig.add_gridspec(2, 3)
XL = ["0.5x", "1x", "2x", "5x", "10x", "Free"]


def gfmt(v):
    return f"{v:.2f}" if v >= 0.01 else f"{v:.1e}".replace("e-0", "e-")


# ---------------------------------------------------------------------------
# A to C. Growth against the capacity of the two urea enzymes
# ---------------------------------------------------------------------------
for k, (key, ttl) in enumerate((("URCB", "Urea carboxylase capacity"), ("ALPH", "Allophanate hydrolase capacity"),
                                ("both", "Both capacities"))):
    ax = fig.add_subplot(gs[0, k])
    head(ax, "ABC"[k], ttl)
    for i in range(len(MULT)):
        for j, (cond, col) in enumerate((("air", BLUE), ("co2", ORANGE))):
            v = CS[cond][key][i]
            x = i + (-0.19 if j == 0 else 0.19)
            vbar(ax, x, v, col, 7.5, b=1e-6)
            val(ax, x, v, gfmt(v), dx=0, dy=6, ha="center", va="bottom", fontsize=7.1, rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(3e-5, 3.0)
    ax.set_xlim(-0.6, len(MULT) - 0.4)
    ax.set_xticks(range(len(MULT)))
    ax.set_xticklabels(XL, fontsize=8.8)
    ax.set_xlabel("Capacity relative to the measured value", labelpad=2)
    ax.set_ylabel("Coupling-model growth (1/h)", labelpad=2)
    style(ax)
    if k == 0:
        legend(ax, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left")

# ---------------------------------------------------------------------------
# D. Model assumptions
# ---------------------------------------------------------------------------
axD = fig.add_subplot(gs[1, 0])
head(axD, "D", "Model assumptions")
grp_ = [("nadh", "Used"), ("nadph", "NADPH\nglutamate\nsynthase"), ("orig_ub", "Original\nbound file")]
for i, (key, lab) in enumerate(grp_):
    for j, (cond, col) in enumerate((("air", BLUE), ("co2", ORANGE))):
        v = CS[cond][key]
        x = i + (-0.19 if j == 0 else 0.19)
        vbar(axD, x, v, col, 11.0, b=1e-6)
        val(axD, x, v, gfmt(v), dx=0, dy=6, ha="center", va="bottom", fontsize=7.8, rotation=90)
axD.set_yscale("log")
axD.set_ylim(3e-5, 0.05)
axD.set_xticks(range(len(grp_)))
axD.set_xticklabels([l for _, l in grp_], fontsize=8.5)
axD.set_xlim(-0.6, len(grp_) - 0.4)
axD.set_ylabel("Coupling-model growth (1/h)", labelpad=2)
style(axD)
legend(axD, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left")

# ---------------------------------------------------------------------------
# E. Fold change of growth, WT-CO2 over WT-Air
# ---------------------------------------------------------------------------
axE = fig.add_subplot(gs[1, 1])
head(axE, "E", "Growth fold, WT-CO$_2$ / WT-Air")
fold_items = [("Measured\ncapacity", CS["co2"]["both"][1] / CS["air"]["both"][1]),
              ("Urea carb.\n10x", CS["co2"]["URCB"][4] / CS["air"]["URCB"][4]),
              ("Both\n10x", CS["co2"]["both"][4] / CS["air"]["both"][4]),
              ("Both\nfree", CS["co2"]["both"][5] / CS["air"]["both"][5]),
              ("Experiment", EXP_GROWTH["co2"] / EXP_GROWTH["air"])]
for i, (lab, v) in enumerate(fold_items):
    vbar(axE, i, v, GREY if lab == "Experiment" else BLUE, 14.0)
    val(axE, i, v, f"{v:.1f}", dx=0, dy=7, ha="center", va="bottom", fontsize=8.5)
axE.set_xticks(range(len(fold_items)))
axE.set_xticklabels([l for l, _ in fold_items], fontsize=8.2)
axE.set_xlim(-0.6, len(fold_items) - 0.4)
axE.set_ylim(0, 1.25 * max(v for _, v in fold_items))
axE.set_ylabel("Fold", labelpad=2)
style(axE)

# ---------------------------------------------------------------------------
# F. Ten-fold relaxation of one bound
# ---------------------------------------------------------------------------
axF = fig.add_subplot(gs[1, 2])
head(axF, "F", "Ten-fold relaxation of one bound")
allr = {}
for cond in COND:
    for r, g_ in CS[cond]["gains"].items():
        allr[r] = max(allr.get(r, 0), g_)
topr = [r for r in sorted(allr, key=lambda r: -allr[r]) if allr[r] > 0.01]
nact = len(allr)
n = max(len(topr), 1)
for i, r in enumerate(topr):
    y = n - 1 - i
    for off, cond, col in ((0.18, "air", BLUE), (-0.18, "co2", ORANGE)):
        v = 100 * CS[cond]["gains"].get(r, 0)
        hbar(axF, y + off, v, col, 8.0)
        val(axF, v, y + off, f"{v:.0f}%", va="center", ha="left", fontsize=8.5)
axF.set_yticks(range(n))
axF.set_yticklabels([rname(r, 24) for r in topr][::-1], fontsize=8.8)
axF.set_ylim(-0.7, n - 0.3)
axF.set_xlim(0, 100 * max(CS[c]["gains"].get(topr[0], 0) for c in COND) * 1.35)
axF.set_xlabel("Growth gain (%)", labelpad=2)
style(axF, left=False)
axF.annotate(f"The other {nact - len(topr)} active bounded reactions\nchange growth by less than 1%", (0.02, 0.97),
             xycoords="axes fraction", ha="left", va="top", fontsize=8.5)
legend(axF, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="lower right")

finish(fig, out_path)