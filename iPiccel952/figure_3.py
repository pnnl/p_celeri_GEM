"""
Fig. 3 for the iPiccel952 manuscript (enzyme-constrained model, relaxation algorithm and calibration).

Usage
    python make_fig3.py ecGEMs.zip Fig3_ecGEM.png leave_out_results.xlsx
where ecGEMs.zip is the zip file with the WT_Air and Wt_CO2 folders and leave_out_results.xlsx is written by
leave_out_test.py (panel I). Or edit ZIP_PATH, OUT_PATH and LEAVE_PATH below.
Requires matplotlib, numpy, pandas and openpyxl. Everything is read from the files you generated. The experimental growth
rates and measured proteome pools are the targets stored in the relaxation_algorithm.py parameter blocks.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ZIP_PATH = "ecGEMs.zip"
OUT_PATH = "Fig3_ecGEM.png"
LEAVE_PATH = "leave_out_results.xlsx"
zip_path = sys.argv[1] if len(sys.argv) > 1 else ZIP_PATH
out_path = sys.argv[2] if len(sys.argv) > 2 else OUT_PATH
leave_path = sys.argv[3] if len(sys.argv) > 3 else LEAVE_PATH

# Read ecGEMs.zip only, extract to a temporary folder and locate WT_Air inside it
import os
import tempfile
import zipfile
tmp = tempfile.mkdtemp(prefix="ecgems_")
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(tmp)
zip_dir = tmp
for root, dirs, _ in os.walk(tmp):
    if "WT_Air" in dirs:
        zip_dir = root
        break

COND = {
    "WT-Air": dict(dir="WT_Air", tag="Wt_Air", fix="iterative_kcat_fix_results_WtAir.xlsx", pfba="ecGEM_pfba_Wt_Air.xlsx",
                   growth_target=0.028, pool_target=0.077759),
    "WT-CO2": dict(dir="Wt_CO2", tag="Wt_CO2", fix="iterative_kcat_fix_results_Wt_CO2.xlsx", pfba="ecGEM_pfba_Wt_CO2.xlsx",
                   growth_target=0.324941, pool_target=0.091287),
}
BIOMASS = "Biomass_Picre_Auto"

# ---------------------------------------------------------------------------
# Read the files
# ---------------------------------------------------------------------------
D = {}
for name, c in COND.items():
    ra = f"{zip_dir}/{c['dir']}/relaxation_algorithm"
    xl = pd.ExcelFile(f"{ra}/{c['fix']}")
    rounds, boosted, summary = xl.parse("rounds"), xl.parse("boosted_reactions"), xl.parse("summary")
    pf = pd.ExcelFile(f"{zip_dir}/{c['dir']}/{c['pfba']}")
    flux, prot = pf.parse("flux"), pf.parse("protein_requirements")
    kc = pd.read_excel(f"{ra}/kcat_mw_{c['tag']}.xlsx", sheet_name="reaction_protein_cost")
    cum = boosted.groupby("Reaction")["Fold_applied"].prod()
    first_round = boosted.groupby("Reaction")["Round"].min()
    D[name] = dict(
        rounds=rounds, cum=cum, first=first_round, kcat=kc,
        n_stage1=int(summary["N_reactions_relaxed_step1"].iloc[0]), n_stage2=len(cum),
        growth=float(flux.loc[flux.Reaction == BIOMASS, "Flux"].iloc[0]),
        pool=float(prot.iloc[:, 1].sum()),
        pool_start=float(rounds["Total_pool"].iloc[0]), n_rounds=len(rounds),
    )
    print(name, {k: v for k, v in D[name].items() if k not in ("rounds", "cum", "first", "kcat")})

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
avail = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in avail), "DejaVu Sans")
plt.rcParams.update({
    "font.family": family, "font.size": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.labelsize": 6.2, "ytick.labelsize": 6.2, "axes.labelsize": 6.5, "legend.frameon": False,
    "legend.fontsize": 6.0, "text.color": "#000000", "axes.labelcolor": "#000000",
    "xtick.color": "#000000", "ytick.color": "#000000",
})
INK, INK2, MUTED, AXIS = "#000000", "#000000", "#000000", "#7a7a7a"
BLUE, ORANGE, AQUA, YELLOW, MAG, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"
CCOL = {"WT-Air": BLUE, "WT-CO2": ORANGE}
TRACK, GREY = "#eef0f3", "#b9b8b0"

fig = plt.figure(figsize=(7.0, 8.9), dpi=100)
gs = fig.add_gridspec(5, 6, height_ratios=[0.34, 1.22, 0.9, 1.5, 1.2], left=0.075, right=0.985, top=0.983, bottom=0.045,
                      hspace=0.46, wspace=1.15)


HEADS = []


def letter(ax, txt, dx=0, dy=0, title=None, wide=False):
    HEADS.append((ax, txt, title, wide))


def style(ax):
    ax.minorticks_off()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, pad=1.5)


def box(ax, x, y, w, h, lines, fc, ec=None, tc=INK, fs=5.8, bold_first=True, rounding=1.2):
    ec = ec or fc
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rounding}", fc=fc, ec=ec, lw=0.8,
                                zorder=2))
    n = len(lines)
    for i, t in enumerate(lines):
        yy = y + h / 2 + (n - 1) * 0.5 * fs * 0.62 - i * fs * 0.62 * 1.18
        ax.text(x + w / 2, yy, t, ha="center", va="center", fontsize=fs, color=tc, zorder=3,
                fontweight="bold" if (bold_first and i == 0) else "normal")


def arrow(ax, p0, p1, color=INK2, rad=0.0, lw=0.9, style_="-|>"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style_, mutation_scale=6, lw=lw, color=color,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1, shrinkA=0, shrinkB=0))


# ---------------------------------------------------------------------------
# A. ecGEM pipeline
# ---------------------------------------------------------------------------
axA = fig.add_subplot(gs[0, :])
axA.set_xlim(0, 100)
axA.set_ylim(2.5, 20.5)
axA.axis("off")
letter(axA, "A", wide=True)
bw, bh, y0, gap = 14.4, 15, 4.5, 2.7
xs = [i * (bw + gap) for i in range(6)]
specs = [
    (["iPiccel952", "reactions split into", "irreversible pairs"], "#dde9f8", INK),
    (["TurNuP", "kcat predicted from", "sequence and substrates"], "#d3efe4", INK),
    (["Protein cost", "molecular weight / kcat", "and kcat-derived bounds"], "#d3efe4", INK),
    (["Proteome pool", "measured abundance,", "WT-Air and WT-CO2"], "#fbe3c0", INK),
    (["Relaxation algorithm", "bounds and kcat", "corrected (B)"], ORANGE, "white"),
    (["Calibrated ecGEM", "growth and pool match", "measurements"], GREEN, "white"),
]
for x, (lines, fc, tc) in zip(xs, specs):
    box(axA, x, y0, bw, bh, lines, fc, tc=tc)
for i in range(5):
    arrow(axA, (xs[i] + bw + 0.2, y0 + bh / 2), (xs[i + 1] - 0.2, y0 + bh / 2))
axA.text(xs[1] + bw / 2, y0 - 1.2, "1,282 reactions", ha="center", va="top", fontsize=5.6, color=MUTED)
axA.text(xs[4] + bw / 2, y0 - 1.2, "two stages", ha="center", va="top", fontsize=5.6, color=MUTED)

# ---------------------------------------------------------------------------
# B. Relaxation algorithm
# ---------------------------------------------------------------------------
axB = fig.add_subplot(gs[1, :])
axB.set_xlim(0, 100)
axB.set_ylim(-0.5, 67.2)
axB.axis("off")
letter(axB, "B", wide=True)
plt.rcParams["mathtext.fontset"] = "stix"
SERIF = ["Liberation Serif", "STIXGeneral", "DejaVu Serif"]
FS = 6.4


def txt(x, y, s_, bold=False, color=INK, fs=FS, ha="left"):
    axB.text(x, y, s_, fontsize=fs, family=SERIF, fontweight="bold" if bold else "normal", color=color, ha=ha,
             va="center")


def code_line(col_x, y, num, indent, s_):
    txt(col_x, y, f"{num}", color=MUTED, fs=5.6)
    txt(col_x + 3.0 + 2.6 * indent, y, s_)


axB.plot([0.5, 99.5], [66.5, 66.5], color=INK, lw=0.9)
txt(1.0, 64.2, "Algorithm 1", bold=True)
txt(13.0, 64.2, "Two-stage relaxation of reaction upper bounds and kcat values", bold=False)
axB.plot([0.5, 99.5], [62.2, 62.2], color=INK, lw=0.5)
txt(1.0, 59.8, "Input", bold=True)
txt(7.0, 59.8, r"$S\,v=0$,  upper bounds $UB_i$,  protein costs $c_i=MW_i/k_{cat,i}$,  targets $\mu^{*}$ and $P^{*}$,  "
    r"threshold $\tau=0.2$,  relaxed bound $U=1000$,  fold $f=100$,  $R_{max}=100$")
txt(1.0, 55.8, "Output", bold=True)
txt(7.0, 55.8, r"relaxed reaction set $B$,  revised protein costs $c_i'$")
axB.plot([0.5, 99.5], [53.4, 53.4], color=INK, lw=0.5)
axB.plot([50, 50], [1.6, 52.8], color="#d8d7cf", lw=0.5)

XL, XR = 1.0, 51.5
txt(XL, 51.0, "Stage 1", bold=True)
txt(XL + 8.0, 51.0, "bound relaxation, pool constraint removed", color=INK2)
st = 5.3
L1 = [
    (0, r"$v^{free}\leftarrow$ FBA with every $UB_i=U$"),
    (0, r"$F\leftarrow\{\,i:(v^{free}_i-UB_i)/UB_i>\tau\,\}$, sorted by deviation, descending"),
    (0, r"$B\leftarrow\emptyset$,   $\mu\leftarrow$ FBA with the original $UB$"),
    (0, r"$\mathbf{for}$ $i\in F$ $\mathbf{do}$"),
    (1, r"$\mathbf{if}$ $\mu\geq\mu^{*}$ $\mathbf{then}$ $\mathbf{break}$"),
    (1, r"$UB_i\leftarrow U$,   $B\leftarrow B\cup\{i\}$"),
    (1, r"$\mu\leftarrow$ FBA"),
    (0, r"$\mathbf{end\ for}$"),
    (0, r"$\mu_f\leftarrow\min(\mu,\mu^{*})$"),
]
for k, (ind, t) in enumerate(L1):
    code_line(XL, 46.4 - k * st, k + 1, ind, t)
txt(XR, 51.0, "Stage 2", bold=True)
txt(XR + 8.0, 51.0, "kcat relaxation at fixed growth", color=INK2)
L2 = [
    (0, r"add the constraint $v_{biomass}\geq\mu_f$"),
    (0, r"$\mathbf{for}$ $r=1,\ldots,R_{max}$ $\mathbf{do}$"),
    (1, r"$P\leftarrow\min\ \sum_{i\notin B}c_i\,v_i$   subject to $S\,v=0$"),
    (1, r"$\mathbf{if}$ $P\leq P^{*}$ $\mathbf{then}$ $\mathbf{break}$"),
    (1, r"$k\leftarrow\arg\max_i\ c_i\,v_i$   (top consumer)"),
    (1, r"$c_k\leftarrow c_k/f$   ($k_{cat,k}\leftarrow f\,k_{cat,k}$)"),
    (0, r"$\mathbf{end\ for}$"),
    (0, r"$\mathbf{return}$ $B$, $c'$"),
]
for k, (ind, t) in enumerate(L2):
    code_line(XR, 46.4 - k * st, 10 + k, ind, t)
axB.plot([0.5, 99.5], [1.0, 1.0], color=INK, lw=0.9)

# ---------------------------------------------------------------------------
# C and D. Growth rate and proteome pool, experiment against ecGEM
# ---------------------------------------------------------------------------
def rbars(ax, items, lw=7.0):
    for x, v, col in items:
        ax.plot([x, x], [0, v], color=col, lw=lw, solid_capstyle="round", zorder=3)


def paired(ax, key_t, key_m, ymax, yticks, off):
    for i, name in enumerate(COND):
        exp = COND[name][key_t]
        mod = D[name][key_m]
        rbars(ax, [(i - 0.17, exp, GREY), (i + 0.17, mod, CCOL[name])])
        ax.text(i - 0.17, exp + off, f"{exp:.3f}", ha="center", va="bottom", fontsize=5.2, color=INK2, rotation=90)
        ax.text(i + 0.17, mod + off, f"{mod:.3f}", ha="center", va="bottom", fontsize=5.2, color=INK, rotation=90)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(list(COND))
    ax.set_yticks(yticks)
    style(ax)


axC = fig.add_subplot(gs[2, 0:2])
letter(axC, "C", title="Growth rate (1/h)")
paired(axC, "growth_target", "growth", 0.42, [0, 0.1, 0.2, 0.3], 0.02)
def dot(c):
    return Line2D([], [], marker="o", ls="", ms=5.5, mfc=c, mec=c)


def pair_legend(ax, entries, **kw):
    hs = [tuple(dot(c) for c in cols) for cols, _ in entries]
    ax.legend(hs, [t for _, t in entries], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
              handletextpad=0.4, labelspacing=0.35, borderpad=0, fontsize=5.6, **kw)


pair_legend(axC, [([GREY], "Experiment"), ([BLUE, ORANGE], "ecGEM")], loc="upper left")

axD = fig.add_subplot(gs[2, 2:4])
letter(axD, "D", title="Proteome pool (g/gDCW)")
paired(axD, "pool_target", "pool", 0.125, [0, 0.05, 0.1], 0.006)

# ---------------------------------------------------------------------------
# E. Reactions corrected per stage
# ---------------------------------------------------------------------------
axE = fig.add_subplot(gs[2, 4:6])
letter(axE, "E", title="Reactions corrected")
for i, name in enumerate(COND):
    y = 1 - i
    s1, s2 = D[name]["n_stage1"], D[name]["n_stage2"]
    axE.plot([0, 60], [y, y], color=TRACK, lw=9, solid_capstyle="round", zorder=1)
    x = 0
    if s1:
        axE.plot([0, s1], [y, y], color=YELLOW, lw=9, solid_capstyle="round", zorder=3)
        x = s1
    axE.plot([x + 1.0, x + s2], [y, y], color=AQUA, lw=9, solid_capstyle="round", zorder=4)
    axE.text(63, y, f"{s1} + {s2}", va="center", fontsize=5.8, color=INK)
axE.set_yticks([1, 0])
axE.set_yticklabels(list(COND))
axE.set_ylim(-0.6, 1.6)
axE.set_xlim(-2, 72)
axE.set_xticks([])
for sp in axE.spines.values():
    sp.set_visible(False)
axE.tick_params(length=0)
axE.scatter([], [], s=18, color=YELLOW, label="Stage 1, bound")
axE.scatter([], [], s=18, color=AQUA, label="Stage 2, kcat")
axE.legend(loc="lower right", bbox_to_anchor=(1.0, -0.12), fontsize=5.6, handletextpad=0.2, labelspacing=0.3,
           borderpad=0)

# ---------------------------------------------------------------------------
# F. Minimized pool at the start, at the end and the measured pool
# ---------------------------------------------------------------------------
def tint(c, a):
    from matplotlib.colors import to_rgb
    r_, g_, b_ = to_rgb(c)
    return (1 - a) + a * r_, (1 - a) + a * g_, (1 - a) + a * b_


axF = fig.add_subplot(gs[3, 0:2])
letter(axF, "F", title="Proteome pool during relaxation")
YB = 0.03
for gi, name in enumerate(COND):
    vals = [D[name]["pool_start"], D[name]["rounds"]["Total_pool"].iloc[-1], COND[name]["pool_target"]]
    cols = [tint(CCOL[name], 0.4), CCOL[name], GREY]
    for k, (v, col) in enumerate(zip(vals, cols)):
        x = gi + (k - 1) * 0.27
        axF.plot([x, x], [YB * 1.15, v], color=col, lw=7.0, solid_capstyle="round", zorder=3)
        axF.text(x, v * 1.25, f"{v:.3g}", ha="center", va="bottom", fontsize=4.8, color=INK, rotation=90)
axF.set_yscale("log")
axF.set_ylim(YB, 40)
axF.set_xlim(-0.55, 1.55)
axF.set_xticks([0, 1])
axF.set_xticklabels([f"{n}\n{D[n]['n_rounds']} rounds" for n in COND])
axF.set_yticks([0.1, 1, 10]); axF.minorticks_off()
axF.set_yticklabels(["0.1", "1", "10"])
axF.set_ylabel("g/gDCW", labelpad=1)
style(axF)
pair_legend(axF, [([tint(BLUE, 0.4), tint(ORANGE, 0.4)], "Start"), ([BLUE, ORANGE], "Final"), ([GREY], "Measured")],
            loc="upper left", bbox_to_anchor=(0.0, 1.02))

# ---------------------------------------------------------------------------
# G. Median TurNuP kcat of all reactions and of the corrected reactions, before and after
# ---------------------------------------------------------------------------
axG = fig.add_subplot(gs[3, 2:4])
letter(axG, "G", title="Median kcat (1/s)")
kc_all = D["WT-CO2"]["kcat"]
med_all = float(np.median(kc_all["Kcat"]) / 3600.0)
items = [(0, med_all, GREY, "All")]
x = 1.2
labels = [("All 1,282", 0.0)]
groups = [(0, med_all, GREY)]
for gi, name in enumerate(COND):
    kcn = D[name]["kcat"]
    cm = D[name]["cum"]
    sub = kcn[kcn["Reaction"].isin(cm.index)].copy()
    before = float(np.median(sub["Kcat"]) / 3600.0)
    after = float(np.median(sub["Kcat"] * sub["Reaction"].map(cm)) / 3600.0)
    groups.append((1.25 + gi * 1.0, before, tint(CCOL[name], 0.4)))
    groups.append((1.55 + gi * 1.0, after, CCOL[name]))
YG = 3.0
for xg, v, col in groups:
    axG.plot([xg, xg], [YG * 1.12, v], color=col, lw=7.0, solid_capstyle="round", zorder=3)
    axG.text(xg, v * 1.3, f"{v:,.0f}", ha="center", va="bottom", fontsize=5.4, color=INK)
axG.set_yscale("log")
axG.set_ylim(YG, 14000)
axG.set_xlim(-0.4, 3.4)
axG.set_xticks([0, 1.4, 2.4])
axG.set_xticklabels(["All 1,282\nreactions", "WT-Air\ncorrected", "WT-CO2\ncorrected"], fontsize=5.8)
axG.set_yticks([10, 100, 1000, 10000])
axG.set_yticklabels(["10", "100", "1,000", "10,000"]); axG.minorticks_off()
style(axG)
pair_legend(axG, [([tint(BLUE, 0.4), tint(ORANGE, 0.4)], "Before"), ([BLUE, ORANGE], "After")], loc="upper left")

# ---------------------------------------------------------------------------
# H. The 16 reactions corrected in both conditions, round in which each was chosen
# ---------------------------------------------------------------------------
axH = fig.add_subplot(gs[3, 4:6])
letter(axH, "H", title="Reactions corrected in both conditions")
shared = [r for r in D["WT-Air"]["cum"].index if r in D["WT-CO2"]["cum"].index]
shared.sort(key=lambda r: (D["WT-Air"]["first"][r], D["WT-CO2"]["first"][r]))
n = len(shared)
for i, rid in enumerate(shared):
    y = n - 1 - i
    rnd = D["WT-CO2"]["first"][rid]
    col = ORANGE if D["WT-CO2"]["cum"][rid] > 100 else AQUA
    axH.plot([0.4, 20], [y, y], color=TRACK, lw=5.2, solid_capstyle="round", zorder=1)
    axH.plot([0.4, rnd], [y, y], color=col, lw=5.2, solid_capstyle="round", zorder=3)
    axH.text(rnd + 0.9, y, f"{int(rnd)}", va="center", fontsize=5.2, color=INK)
axH.set_yticks(range(n))
axH.set_yticklabels(shared[::-1], fontsize=5.0)
axH.set_xlim(-0.5, 33)
axH.set_ylim(-0.7, n - 0.3)
axH.set_xlabel("Round chosen", labelpad=1)
style(axH)
axH.spines["left"].set_visible(False)
axH.tick_params(axis="y", length=0, pad=1)
axH.scatter([], [], s=16, color=AQUA, label="100-fold")
axH.scatter([], [], s=16, color=ORANGE, label="10,000-fold")
axH.legend(loc="lower right", bbox_to_anchor=(1.04, 0.0), handletextpad=0.1, labelspacing=0.2, borderpad=0,
           fontsize=5.2)

# ---------------------------------------------------------------------------
# I. Leave-out test, growth predicted under the measured pool with kcat calibrated on the other condition
# ---------------------------------------------------------------------------
def shade(c, a):
    from matplotlib.colors import to_rgb
    r_, g_, b_ = to_rgb(c)
    return (1 - a) * r_, (1 - a) * g_, (1 - a) * b_


lo = pd.read_excel(leave_path, sheet_name="leave_out")
axI = fig.add_subplot(gs[4, :])
letter(axI, "I", title="Growth under the measured pool with kcat calibrated on the other condition (1/h)")
YI = 1e-4
for gi, name in enumerate(COND):
    base = gi * 5.0
    sub = lo[lo.target == name].set_index("calibration")["growth"]
    bars = [(sub["measured"], GREY), (sub["none"], tint(CCOL[name], 0.4)), (sub["other"], CCOL[name]),
            (sub["self"], shade(CCOL[name], 0.35))]
    for k, (v, col) in enumerate(bars):
        x = base + k
        top = max(v, YI * 1.3)
        axI.plot([x, x], [YI, top], color=col, lw=14, solid_capstyle="butt", zorder=3)
        label = f"{v:.3g}"
        if k == 2:
            label += f"\n{100 * v / sub['measured']:.0f}% of measured"
        axI.text(x, top * 1.12, label, ha="center", va="bottom", fontsize=5.6, color=INK, linespacing=1.15)
axI.set_yscale("log")
axI.set_ylim(YI, 20)
axI.set_xlim(-0.8, 8.8)
axI.set_xticks([1.5, 6.5])
axI.set_xticklabels(["WT-Air", "WT-CO2"])
axI.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
axI.set_yticklabels(["0.0001", "0.001", "0.01", "0.1", "1"])
axI.minorticks_off()
style(axI)
pair_legend(axI, [([GREY], "Measured"), ([tint(BLUE, 0.4), tint(ORANGE, 0.4)], "Uncalibrated"),
                  ([BLUE, ORANGE], "Calibrated on the other condition"),
                  ([shade(BLUE, 0.35), shade(ORANGE, 0.35)], "Calibrated on itself")],
            loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=4, columnspacing=1.4)

# ---------------------------------------------------------------------------
# Panel letters and titles, placed from the final axes positions
# ---------------------------------------------------------------------------
for ax, txt, title, wide in HEADS:
    p = ax.get_position()
    if wide:
        fig.text(0.008, p.y1 + 0.004, txt, fontsize=9, fontweight="bold", va="bottom", ha="left", color=INK)
        continue
    x0 = p.x0 - 0.062 if txt in ("C", "D", "E", "F", "G") else p.x0 - 0.062
    x0 = 0.008 if p.x0 < 0.2 else p.x0 - 0.068
    fig.text(x0, p.y1 + 0.012, txt, fontsize=9, fontweight="bold", va="bottom", ha="left", color=INK)
    fig.text(x0 + 0.022, p.y1 + 0.013, title, fontsize=6.6, color=INK2, va="bottom", ha="left")

fig.savefig(out_path, dpi=500, facecolor="white")
print("saved", out_path)