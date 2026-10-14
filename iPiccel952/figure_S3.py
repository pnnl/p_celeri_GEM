"""
Supplementary Fig. S3 for the eciPiccel952 manuscript (supports Fig. 3, protein pool calibration).

Usage
    python make_figS3.py ecGEMs.zip iPiccel952.xml FigS3.png
    ecGEMs.zip      zip file with the WT_Air and Wt_CO2 folders (relaxation algorithm results and revised files)
    iPiccel952.xml  original SBML model, used only for reaction names
or edit ZIP_PATH, MODEL_PATH and OUT_PATH below.
Requires cobra, matplotlib, numpy, pandas and openpyxl. The run solves two LPs and takes about a minute.

Panels
    A, B  protein pool demand in each relaxation round, WT-Air and WT-CO2
    C     protein cost of the ten most costly boosted reactions in WT-CO2, before and after the boost
    D     range of the TurNuP kcat values of the enzyme-constrained reactions
    E     the ten largest consumers of the protein pool in WT-Air and their share in WT-CO2
    F     reactions with a protein cost, with a kcat-derived bound and carrying flux
"""
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

ZIP_PATH = "ecGEMs.zip"
MODEL_PATH = "iPiccel952.xml"
OUT_PATH = "FigS3.png"
zip_path = sys.argv[1] if len(sys.argv) > 1 else ZIP_PATH
model_path = sys.argv[2] if len(sys.argv) > 2 else MODEL_PATH
out_path = sys.argv[3] if len(sys.argv) > 3 else OUT_PATH

BIO = "Biomass_Picre_Auto"
COND = {"air": ("WT_Air", "Wt_Air", 0.077759, "WT-Air"), "co2": ("Wt_CO2", "Wt_CO2", 0.091287, "WT-CO2")}

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
# Calibrated ecGEMs, baseline solution under the measured pool
# ---------------------------------------------------------------------------
def baseline(cond):
    """Reaction table of the calibrated ecGEM, columns reaction, cost, ub, flux and pool_use."""
    folder, tag, pool, _ = COND[cond]
    d = f"{root}/{folder}"
    m = cobra.io.read_sbml_model(f"{d}/irrev_iPiccel952.xml")
    ub = pd.read_excel(f"{d}/UB_{tag}_revised.xlsx", usecols=[0, 1], names=["r", "u"])
    ids = {r.id: r for r in m.reactions}
    applied = {}
    for r, u in zip(ub.r, ub.u):
        r = str(r)
        if r in ids and pd.notna(u):
            ids[r].upper_bound = max(0.0, float(u))
            applied[r] = float(u)
    cost = pd.read_excel(f"{d}/kcat_mw_{tag}_revised.xlsx", sheet_name="reaction_protein_cost")
    pc = {str(r): float(c) for r, c in zip(cost.Reaction, cost.Protein_cost) if str(r) in ids and pd.notna(c)}
    m.add_cons_vars(m.problem.Constraint(sum((c / 1000.0) * ids[r].flux_expression for r, c in pc.items()),
                                         lb=0, ub=pool, name="pool"))
    m.solver.update()
    m.objective = BIO
    sol = m.optimize()
    print(f"{cond} growth {sol.objective_value:.5f}", flush=True)
    rows = [(r, pc[r], applied.get(r, np.nan), sol.fluxes[r]) for r in pc]
    t = pd.DataFrame(rows, columns=["reaction", "cost", "ub", "flux"])
    t["pool_use"] = t.cost * t.flux / 1000.0
    return t


ca, cc = baseline("air"), baseline("co2")
rx_air = pd.read_excel(f"{root}/WT_Air/relaxation_algorithm/iterative_kcat_fix_results_WtAir.xlsx", sheet_name=None)
rx_co2 = pd.read_excel(f"{root}/Wt_CO2/relaxation_algorithm/iterative_kcat_fix_results_Wt_CO2.xlsx", sheet_name=None)

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


fig = plt.figure(figsize=(7.0, 5.9), dpi=100, layout="constrained")
fig.get_layout_engine().set(w_pad=0.04, h_pad=0.06, wspace=0.12, hspace=0.14, rect=(0.0, 0.0, 1.0, 0.955))
gs = fig.add_gridspec(2, 3)

# ---------------------------------------------------------------------------
# A and B. Protein pool demand in each relaxation round
# ---------------------------------------------------------------------------
def rounds_panel(ax, letter, d_, label):
    rd = d_["rounds"]
    tgt = float(rd.Target_pool.iloc[0])
    n = len(rd)
    for i, v in enumerate(rd.Total_pool):
        vbar(ax, rd.Round.iloc[i], v, GREEN if v <= tgt * 1.02 else ORANGE, 4.6 if n < 25 else 2.4, b=0.01)
    ax.set_yscale("log")
    ax.set_ylim(0.03, 1.6 * rd.Total_pool.max())
    ax.set_xlim(0, n + 1)
    ax.set_xlabel("Relaxation round", labelpad=2)
    ax.set_ylabel("Protein pool demand (g/gDW)", labelpad=2)
    ax.annotate(f"{n} rounds, {int(d_['summary'].N_reactions_boosted_step2.iloc[0])} boosted\ntarget pool {tgt:.3f} g/gDW",
                (0.97, 0.95), xycoords="axes fraction", ha="right", va="top", fontsize=5.4)
    style(ax)
    head(ax, letter, label)


axA = fig.add_subplot(gs[0, 0])
rounds_panel(axA, "A", rx_air, "WT-Air calibration")
axB = fig.add_subplot(gs[0, 1])
rounds_panel(axB, "B", rx_co2, "WT-CO$_2$ calibration")
legend(axB, [([ORANGE], "Above target"), ([GREEN], "Target reached")], loc="upper right", bbox_to_anchor=(1.0, 0.80))

# ---------------------------------------------------------------------------
# C. Most costly boosted reactions, WT-CO2
# ---------------------------------------------------------------------------
bo = rx_co2["boosted_reactions"].sort_values("Old_cost", ascending=False).head(10)
axC = fig.add_subplot(gs[0, 2])
head(axC, "C", "Boosted reactions, WT-CO$_2$")
n = len(bo)
for i, (_, r) in enumerate(bo.iterrows()):
    y = n - 1 - i
    for off, v, c in ((0.18, r.Old_cost, ORANGE), (-0.18, r.New_cost, BLUE)):
        hbar(axC, y + off, np.log10(v) + 3, c, 3.4)
axC.set_yticks(range(n))
axC.set_yticklabels([rname(r) for r in bo.Reaction][::-1], fontsize=4.8)
axC.set_ylim(-0.7, n - 0.3)
axC.set_xlim(0, 5.8)
axC.set_xticks([0, 1, 2, 3, 4])
axC.set_xticklabels(["0.001", "0.01", "0.1", "1", "10"])
axC.set_xlabel("Protein cost (g h/mmol)", labelpad=2)
style(axC, left=False)
legend(axC, [([ORANGE], "Before"), ([BLUE], "After")], loc="lower right")

# ---------------------------------------------------------------------------
# D. Range of kcat values
# ---------------------------------------------------------------------------
kc = pd.read_excel(f"{root}/WT_Air/kcat_mw_Wt_Air_revised.xlsx", sheet_name="reaction_protein_cost")
lk = np.log10(kc.Kcat.astype(float) / 3600.0)
edges = np.arange(-4, 5, 1)
cnts = [int(((lk >= a_) & (lk < a_ + 1)).sum()) for a_ in edges[:-1]]
axD = fig.add_subplot(gs[1, 0])
head(axD, "D", "Range of kcat values")
for i, c in enumerate(cnts):
    vbar(axD, i, c, BLUE, 8.0)
    val(axD, i, c, f"{c}", dx=0, dy=3, ha="center", va="bottom", fontsize=4.8)
axD.set_xticks(range(len(cnts)))
axD.set_xticklabels([f"{int(a_)}" for a_ in edges[:-1]], fontsize=5.4)
axD.set_ylim(0, 1.2 * max(cnts))
axD.set_xlim(-0.6, len(cnts) - 0.4)
axD.set_xlabel("log$_{10}$ kcat (1/s), lower bin edge", labelpad=2)
axD.set_ylabel("Enzyme-constrained reactions", labelpad=2)
style(axD)

# ---------------------------------------------------------------------------
# E. Pool allocation, top consumers
# ---------------------------------------------------------------------------
topu = ca.sort_values("pool_use", ascending=False).head(10)
axE = fig.add_subplot(gs[1, 1])
head(axE, "E", "Pool allocation, top consumers")
n = len(topu)
sa_, sc_ = ca.pool_use.sum(), cc.pool_use.sum()
cmap = cc.set_index("reaction").pool_use
for i, (_, r) in enumerate(topu.iterrows()):
    y = n - 1 - i
    va_, vc_ = 100 * r.pool_use / sa_, 100 * cmap.get(r.reaction, 0) / sc_
    hbar(axE, y + 0.18, va_, BLUE, 3.4)
    hbar(axE, y - 0.18, vc_, ORANGE, 3.4)
    val(axE, va_, y + 0.18, f"{va_:.1f}", va="center", ha="left", fontsize=4.4)
    val(axE, vc_, y - 0.18, f"{vc_:.1f}", va="center", ha="left", fontsize=4.4)
axE.set_yticks(range(n))
axE.set_yticklabels([rname(r) for r in topu.reaction][::-1], fontsize=4.8)
axE.set_ylim(-0.7, n - 0.3)
axE.set_xlabel("Share of the protein pool (%)", labelpad=2)
axE.set_xlim(0, 9.5)
style(axE, left=False)
legend(axE, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="lower right", bbox_to_anchor=(1.0, 0.0))

# ---------------------------------------------------------------------------
# F. Constrained reactions
# ---------------------------------------------------------------------------
axF = fig.add_subplot(gs[1, 2])
head(axF, "F", "Constrained reactions")
cats = ["Cost\nassigned", "Bound\nassigned", "Active"]
for k, (d_, col) in enumerate(((ca, BLUE), (cc, ORANGE))):
    vals = [len(d_), int(d_.ub.notna().sum()), int((d_.flux > 1e-9).sum())]
    for i, v in enumerate(vals):
        x = i + (-0.19 if k == 0 else 0.19)
        vbar(axF, x, v, col, 7.0)
        val(axF, x, v, f"{v:,}", dx=0, dy=4, ha="center", va="bottom", fontsize=4.4, rotation=90)
axF.set_xticks(range(3))
axF.set_xticklabels(cats, fontsize=5.2)
axF.set_ylim(0, 1.25 * len(ca))
axF.set_xlim(-0.6, 2.6)
axF.set_ylabel("Reactions", labelpad=2)
style(axF)
legend(axF, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper right")

finish(fig, out_path)