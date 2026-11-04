"""
Supplementary Fig. S5 for the eciPiccel954 manuscript (supports Fig. 5, sensitivity of the bottleneck analysis).

Usage
    python figure_S5.py ecGEMs.zip iPiccel954.xml work_dir FigS5.png [CO2_DRAWS] [SEED] [robust_only]
    ecGEMs.zip      zip file with the WT_Air and Wt_CO2 folders
    iPiccel954.xml  original SBML model, used only for reaction names
    work_dir        folder for the cached results. Every file in it is reused when it exists, delete a file to recompute
                    it. The files are control_air.csv, control_co2.csv, control_air_1.05.csv, control_air_1.20.csv,
                    scan_*.csv, pool_scaling.csv, robust_<seed>.csv (WT-Air draws) and robust_co2_<seed>.csv (WT-CO2 draws).
                    control_air.csv, control_co2.csv and robust_1.csv can be copied from the work folder of Fig. 5.
    CO2_DRAWS       number of kcat perturbation draws for WT-CO2 (default 0, which only uses the cached draws)
    SEED            random seed of the draws (default 1). Several seeds can be run in parallel with robust_only as the
                    last argument and are combined when the figure is drawn.
    ONLY            optional environment variable, comma separated, to run only some computations and stop. The names
                    are c110, c105, c120, scan_cost10, scan_cost100, scan_cost0, pool and rco2.
Requires cobra, matplotlib, numpy, pandas, scipy and openpyxl. The first run solves several thousand LPs. With GLPK
the control scans take about 6 minutes each and each kcat perturbation draw takes one to two minutes.

Pool-limited analysis. As in Fig. 5 there is no growth target and no kcat-derived bound, growth is limited by the protein
pool only, and the protein costs are the calibrated ones.

kcat error. Every draw multiplies each kcat by exp(eps), eps ~ N(0, sigma^2) with sigma = ln 2 (two-fold typical error).
This changes the protein costs, and the control coefficients of the reactions that carry flux are recomputed. The pool is
not recalibrated, so the growth of a draw is reported with it.

Panels
    A  number of bottleneck reactions against the control coefficient threshold
    B  largest growth gain when the protein cost of one reaction is lowered at a time
    C  growth against the size of the protein pool
    D  step size of the kcat increase
    E  frequency of the WT-Air bottlenecks across kcat draws
    F  rank agreement of the control coefficients with the calibrated ones, WT-Air
    G  number of bottlenecks in the calibrated model and in every kcat draw, WT-Air and WT-CO2
    H  frequency of the WT-CO2 bottlenecks across kcat draws
"""
import glob
import os
import re
import sys
import tempfile
import time
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
from scipy import stats

warnings.filterwarnings("ignore")
cobra.Configuration().processes = 1
args = sys.argv[1:]
if len(args) < 4:
    sys.exit(__doc__)
zip_path, model_path, work, out_path = args[:4]
NDRAW = int(args[4]) if len(args) > 4 else 0
SEED = int(args[5]) if len(args) > 5 else 1
ROBUST_ONLY = len(args) > 6 and args[6] == "robust_only"
os.makedirs(work, exist_ok=True)

BIO = "Biomass_Picre_Auto"
THR = 0.01
POOL_ONLY = True      # control analysis at the pool-limited optimum, as in Fig. 5
ONLY = [x for x in os.environ.get("ONLY", "").split(",") if x]   # run only these computations, then stop
want = lambda k: (not ONLY) or k in ONLY
SIGMA = np.log(2.0)
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
t_start = time.time()
orig = cobra.io.read_sbml_model(model_path)
names = {r.id: r.name for r in orig.reactions}


def rname(r, n=40, width=28):
    """Full enzyme name, wrapped over two lines when it is long."""
    import textwrap
    s_ = names.get(r) or names.get(base(r)) or r
    s_ = s_.replace("_", " ").replace(":", ": ").replace(":  ", ": ")
    if len(s_) > n:
        s_ = s_[: n - 1].rstrip(" ,;:-(/") + "."
    return "\n".join(textwrap.wrap(s_, width, break_long_words=True))


# ---------------------------------------------------------------------------
# Computations (cached)
# ---------------------------------------------------------------------------
def build(cond):
    folder, tag, pool, _ = COND[cond]
    d = f"{root}/{folder}"
    m = cobra.io.read_sbml_model(f"{d}/irrev_iPiccel954.xml")
    ub = pd.read_excel(f"{d}/UB_{tag}_revised.xlsx", usecols=[0, 1], names=["r", "u"])
    m.solver = "glpk"
    ids = {r.id: r for r in m.reactions}
    applied = {}
    # Pool-limited analysis. The growth target and the kcat-derived bounds are dropped, so growth is limited by the
    # protein pool only. The protein costs are the calibrated (revised) ones.
    if not POOL_ONLY:
        for r, u in zip(ub.r, ub.u):
            r = str(r)
            if r in ids and pd.notna(u):
                ids[r].upper_bound = max(0.0, float(u))
                applied[r] = float(u)
    cost = pd.read_excel(f"{d}/kcat_mw_{tag}_revised.xlsx", sheet_name="reaction_protein_cost")
    pc = {str(r): float(c) for r, c in zip(cost.Reaction, cost.Protein_cost) if str(r) in ids and pd.notna(c)}
    con = m.problem.Constraint(sum((c / 1000.0) * ids[r].flux_expression for r, c in pc.items()), lb=0, ub=pool, name="pool")
    m.add_cons_vars(con)
    m.solver.update()
    m.objective = BIO
    return m, ids, applied, pc, con


def control_coefficients(cond, step=1.1):
    """C_i = (mu_i - mu_0) / (delta mu_0), kcat of reaction i raised by delta = step - 1."""
    tag = "" if abs(step - 1.1) < 1e-9 else f"_{step:.2f}"
    path = f"{work}/control_{cond}{tag}.csv"
    if os.path.exists(path):
        return
    m, ids, applied, pc, con = build(cond)
    s0 = m.optimize()
    g0, fl = s0.objective_value, s0.fluxes
    print(f"control coefficients {cond} step {step}  baseline growth {g0:.5f}", flush=True)
    rows, t0 = [], time.time()
    for i, r in enumerate(pc):
        rx = ids[r]
        v = rx.forward_variable
        oldc, oldu = con.get_linear_coefficients([v])[v], rx.upper_bound
        con.set_linear_coefficients({v: oldc / step})
        if oldu < 1000:
            rx.upper_bound = min(1000.0, oldu * step)
        s = m.optimize()
        g = s.objective_value if s.status == "optimal" else np.nan
        con.set_linear_coefficients({v: oldc})
        rx.upper_bound = oldu
        rows.append((r, pc[r], applied.get(r, np.nan), fl[r], (g - g0) / g0 / (step - 1)))
        if i % 300 == 0:
            print(f"   {i}/{len(pc)}  {time.time() - t0:.0f} s", flush=True)
    d_ = pd.DataFrame(rows, columns=["reaction", "cost", "ub", "flux", "C"])
    d_["pool_use"] = d_.cost * d_.flux / 1000.0
    d_.to_csv(path, index=False)


def robustness(cond, ndraw, seed=1, step=1.1):
    """kcat error analysis for one condition, written to robust_<seed>.csv (WT-Air) or robust_co2_<seed>.csv."""
    path = f"{work}/robust_{seed}.csv" if cond == "air" else f"{work}/robust_{cond}_{seed}.csv"
    done = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame(columns=["seed", "draw", "reaction", "growth0", "C"])
    k0 = int(done.draw.max()) + 1 if len(done) else 0
    if k0 >= ndraw:
        return
    m, ids, applied, pc, con = build(cond)
    base_cost = {r: con.get_linear_coefficients([ids[r].forward_variable])[ids[r].forward_variable] for r in pc}
    rng = np.random.default_rng(seed)
    out = [tuple(x) for x in done.itertuples(index=False)]
    for k in range(ndraw):
        f = np.exp(rng.normal(0, SIGMA, len(pc)))
        if k < k0:
            continue
        for (r, c), fi in zip(base_cost.items(), f):
            con.set_linear_coefficients({ids[r].forward_variable: c / fi})
            if applied.get(r, 1000) < 1000:
                ids[r].upper_bound = min(1000.0, applied[r] * fi)
        s0 = m.optimize()
        g0, fl = s0.objective_value, s0.fluxes
        act_ = [r for r in pc if fl[r] > 1e-9]
        for r in act_:
            v = ids[r].forward_variable
            c, u = con.get_linear_coefficients([v])[v], ids[r].upper_bound
            con.set_linear_coefficients({v: c / step})
            if u < 1000:
                ids[r].upper_bound = min(1000.0, u * step)
            g = m.optimize().objective_value
            con.set_linear_coefficients({v: c})
            ids[r].upper_bound = u
            out.append((seed, k, r, g0, (g - g0) / g0 / (step - 1)))
        print(f"robustness {cond} seed {seed} draw {k + 1}/{ndraw}  growth {g0:.5f}  active reactions {len(act_)}", flush=True)
        pd.DataFrame(out, columns=["seed", "draw", "reaction", "growth0", "C"]).to_csv(path, index=False)
        for r, c in base_cost.items():
            con.set_linear_coefficients({ids[r].forward_variable: c})
            if r in applied:
                ids[r].upper_bound = applied[r]


def single_scan(mode):
    """Lower the protein cost of one reaction at a time (WT-Air) and record the relative growth gain.
    cost10    protein cost divided by 10
    cost100   protein cost divided by 100
    cost0     protein cost removed"""
    path = f"{work}/scan_{mode}.csv"
    if os.path.exists(path):
        return
    f = {"cost10": 10.0, "cost100": 100.0, "cost0": np.inf}[mode]
    m, ids, applied, pc, con = build("air")
    g0 = m.optimize().objective_value
    print(f"single-reaction scan {mode}  baseline {g0:.5f}", flush=True)
    rows = []
    for i, r in enumerate(pc):
        v = ids[r].forward_variable
        oldc = con.get_linear_coefficients([v])[v]
        con.set_linear_coefficients({v: oldc / f})
        s_ = m.optimize()
        g = s_.objective_value if s_.status == "optimal" else np.nan
        con.set_linear_coefficients({v: oldc})
        rows.append((r, oldc, g, (g - g0) / g0))
        if i % 300 == 0:
            print(f"   {i}", flush=True)
    pd.DataFrame(rows, columns=["reaction", "cost", "growth", "rel_gain"]).to_csv(path, index=False)


def pool_scaling():
    path = f"{work}/pool_scaling.csv"
    if os.path.exists(path):
        return
    rows = []
    for cond in COND:
        m, ids, applied, pc, con = build(cond)
        P = COND[cond][2]
        for k in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0):
            con.ub = P * k
            rows.append((cond, k, m.optimize().objective_value))
        con.ub = 1e6
        rows.append((cond, 1e6, m.optimize().objective_value))
    pd.DataFrame(rows, columns=["cond", "k", "growth"]).to_csv(path, index=False)


if ROBUST_ONLY:
    robustness("co2", NDRAW, SEED)
    sys.exit(0)
if want("c110"):
    for c_ in COND:
        control_coefficients(c_)
if want("c105"):
    control_coefficients("air", 1.05)
if want("c120"):
    control_coefficients("air", 1.2)
for mode in ("cost10", "cost100", "cost0"):
    if want("scan_" + mode):
        single_scan(mode)
if want("pool"):
    pool_scaling()
if NDRAW > 0 and want("rco2"):
    robustness("co2", NDRAW, SEED)
if ONLY:
    sys.exit(0)
print(f"computations ready after {time.time() - t_start:.0f} s", flush=True)

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
        fig.text(x0, p.y1 + 0.007, letter, fontsize=15.3, fontweight="bold", va="bottom", ha="left")
        fig.text(x0 + 0.017, p.y1 + 0.008, title, fontsize=11.2, va="bottom", ha="left")
    fig.savefig(path, dpi=400, facecolor="white")
    plt.close(fig)
    print("saved", path, flush=True)


ca, cc = pd.read_csv(f"{work}/control_air.csv"), pd.read_csv(f"{work}/control_co2.csv")
ca12 = pd.read_csv(f"{work}/control_air_1.20.csv")
ca05 = pd.read_csv(f"{work}/control_air_1.05.csv")

fig = plt.figure(figsize=(12.6, 14.6), dpi=100)
gs = fig.add_gridspec(3, 3, left=0.075, right=0.985, top=0.955, bottom=0.05, wspace=0.78, hspace=0.36,
                      width_ratios=[1.0, 1.0, 1.0], height_ratios=[1.0, 1.0, 1.3])

# ---------------------------------------------------------------------------
# A. Threshold
# ---------------------------------------------------------------------------
taus = [0.001, 0.003, 0.01, 0.03, 0.1]
axA = fig.add_subplot(gs[0, 0])
head(axA, "A", "Bottlenecks against threshold")
for i, tau in enumerate(taus):
    for j, (d_, col) in enumerate(((ca, BLUE), (cc, ORANGE))):
        v = int((d_.C >= tau).sum())
        x = i + (-0.19 if j == 0 else 0.19)
        vbar(axA, x, max(v, 0.3), col, 10.0, b=0.3)
        val(axA, x, max(v, 0.3), f"{v}", dx=0, dy=6, ha="center", va="bottom", fontsize=8.2)
axA.set_yscale("log")
axA.set_ylim(0.3, 300)
axA.set_xticks(range(len(taus)))
axA.set_xticklabels([f"{t:g}" for t in taus], fontsize=9.2)
axA.set_xlim(-0.6, len(taus) - 0.4)
axA.set_xlabel("Control coefficient threshold", labelpad=2)
axA.set_ylabel("Bottleneck reactions", labelpad=2)
style(axA)
legend(axA, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper right")

# ---------------------------------------------------------------------------
# B. Single relaxations
# ---------------------------------------------------------------------------
modes = [("cost10", "Protein cost\nlowered 10-fold"), ("cost100", "Protein cost\nlowered 100-fold"),
         ("cost0", "Protein cost\nremoved")]
axB = fig.add_subplot(gs[0, 1])
head(axB, "B", "One reaction relaxed at a time")
gmax, n10 = [], []
for mname, _ in modes:
    d_ = pd.read_csv(f"{work}/scan_{mname}.csv")
    gmax.append(100 * d_.rel_gain.max())
    n10.append(int((d_.rel_gain >= 0.1).sum()))
for i, v in enumerate(gmax):
    y = len(modes) - 1 - i
    hbar(axB, y, max(v, 0.0), BLUE, 11.0)
    val(axB, max(v, 0.0), y, f"{max(v, 0.0):.1f}%\n{n10[i]} reaction with gain of 10% or more" if n10[i] == 1 else f"{max(v, 0.0):.1f}%\n{n10[i]} reactions with gain of 10% or more", va="center", ha="left", fontsize=8.5, linespacing=1.25)
axB.set_yticks(range(len(modes)))
axB.set_yticklabels([m_ for _, m_ in modes][::-1], fontsize=8.5)
axB.set_ylim(-0.7, len(modes) - 0.3)
axB.set_xlim(0, 150); axB.set_xticks([0, 20, 40, 60, 80])
axB.set_xlabel("Largest growth gain (%)", labelpad=2)
style(axB, left=False)

# ---------------------------------------------------------------------------
# C. Pool size
# ---------------------------------------------------------------------------
ps = pd.read_csv(f"{work}/pool_scaling.csv")
KS = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0]
axC = fig.add_subplot(gs[0, 2])
head(axC, "C", "Growth against pool size")
rel_max = 0
for j, (cond, col) in enumerate((("air", BLUE), ("co2", ORANGE))):
    d_ = ps[ps.cond == cond].set_index("k").growth
    ref = d_[1.0]
    for i, k in enumerate(KS):
        v = d_[k] / ref
        rel_max = max(rel_max, v)
        x = i + (-0.19 if j == 0 else 0.19)
        vbar(axC, x, v, col, 8.5)
        val(axC, x, v, f"{v:.2f}", dx=0, dy=6, ha="center", va="bottom", fontsize=7.5, rotation=90)
    print(cond, "growth without the pool relative to calibrated", round(d_[1e6] / ref, 2), flush=True)
axC.set_xticks(range(len(KS)))
axC.set_xticklabels([f"{k:g}" for k in KS], fontsize=9)
axC.set_xlim(-0.6, len(KS) - 0.4)
axC.set_ylim(0, 1.4 * rel_max)
axC.set_xlabel("Pool relative to calibrated", labelpad=2)
axC.set_ylabel("Growth relative to calibrated", labelpad=2)
style(axC)
legend(axC, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left")

# ---------------------------------------------------------------------------
# D. Step size
# ---------------------------------------------------------------------------
sets = {"5%": set(ca05.reaction[ca05.C >= THR]), "10%": set(ca.reaction[ca.C >= THR]), "20%": set(ca12.reaction[ca12.C >= THR])}
ref = sets["10%"]
axD = fig.add_subplot(gs[1, 0])
head(axD, "D", "Step size of the kcat increase")
for i, (k, s_) in enumerate(sets.items()):
    vbar(axD, i, len(s_), GREY, 15.0)
    vbar(axD, i, len(s_ & ref), BLUE, 15.0)
    val(axD, i, len(s_), f"{len(s_)}\n{len(s_ & ref)} shared", dx=0, dy=8, ha="center", va="bottom", fontsize=8.2)
axD.set_xticks(range(3))
axD.set_xticklabels(list(sets), fontsize=9.2)
axD.set_xlim(-0.6, 2.6)
axD.set_ylim(0, 1.3 * max(len(s_) for s_ in sets.values()))
axD.set_xlabel("kcat increase", labelpad=2)
axD.set_ylabel("Bottleneck reactions, WT-Air", labelpad=2)
style(axD)
legend(axD, [([BLUE], "Shared with 10%"), ([GREY], "Other")], loc="upper right")
rho_ = stats.spearmanr(ca.set_index("reaction").C, ca05.set_index("reaction").C.reindex(ca.reaction).values, nan_policy="omit")[0]
print("Spearman of control coefficients, 10% against 5%", round(rho_, 3), flush=True)


# ---------------------------------------------------------------------------
# kcat perturbation draws
# ---------------------------------------------------------------------------
def load_draws(pattern):
    files = sorted(glob.glob(f"{work}/{pattern}"))
    if not files:
        return None
    rb = pd.concat([pd.read_csv(f) for f in files])
    rb["run"] = rb.seed.astype(str) + "_" + rb.draw.astype(str)
    return rb


rb_air = load_draws("robust_[0-9]*.csv")
rb_co2 = load_draws("robust_co2_*.csv")

# E and F. WT-Air
if rb_air is not None:
    nrun = rb_air.run.nunique()
    freq = (rb_air[rb_air.C >= THR].groupby("reaction").run.nunique() / nrun).sort_values(ascending=False)
    bot = list(ca.sort_values("C", ascending=False).reaction.head(12))
    axE = fig.add_subplot(gs[1, 1])
    head(axE, "E", f"WT-Air frequency, {nrun} draws")
    n = len(bot)
    for i, r in enumerate(bot):
        y = n - 1 - i
        v = 100 * freq.get(r, 0.0)
        hbar(axE, y, max(v, 0.6), BLUE, 7.5)
        val(axE, max(v, 0.6), y, f"{v:.0f}", va="center", ha="left", fontsize=8.2)
    axE.set_yticks(range(n))
    axE.set_yticklabels([rname(r) for r in bot][::-1], fontsize=7.8)
    axE.set_ylim(-0.7, n - 0.3)
    axE.set_xlim(0, 118)
    axE.set_xlabel("Draws as a bottleneck (%)", labelpad=2)
    style(axE, left=False)

    axF = fig.add_subplot(gs[1, 2])
    head(axF, "F", "Rank agreement with C")
    c0 = ca.set_index("reaction").C
    rhos = []
    for run, g_ in rb_air.groupby("run", sort=False):
        s_ = g_.set_index("reaction").C
        common = s_.index.intersection(c0.index)
        rhos.append(stats.spearmanr(c0[common], s_[common])[0])
    xs_ = np.arange(len(rhos))
    med_ = float(np.median(rhos))
    axF.axhline(med_, color=AXIS, lw=1.0, ls=(0, (4, 3)), zorder=1)
    axF.scatter(xs_, rhos, s=46, color=BLUE, zorder=3, linewidths=0)
    axF.text(0.03, 0.04, f"median {med_:.2f}", transform=axF.transAxes, ha="left", va="bottom", fontsize=9.5)
    axF.set_xticks(list(range(0, len(rhos), 5)) + ([len(rhos) - 1] if (len(rhos) - 1) % 5 else []))
    axF.set_xticklabels([str(i + 1) for i in axF.get_xticks()], fontsize=9.2)
    axF.set_xlim(-0.8, len(rhos) - 0.2)
    axF.set_ylim(max(0.0, np.floor((min(rhos) - 0.03) * 20) / 20), 1.0)
    axF.set_xlabel("kcat perturbation draw", labelpad=2)
    axF.set_ylabel("Spearman correlation with C", labelpad=2)
    style(axF)
    print("WT-Air rank correlation", [round(v, 2) for v in rhos], flush=True)

# G and H. WT-CO2 and the comparison of the two conditions
if rb_air is not None and rb_co2 is not None:
    def n_bottle(rb):
        return [int((g_.C >= THR).sum()) for _, g_ in rb.groupby("run", sort=False)]

    na, nc = n_bottle(rb_air), n_bottle(rb_co2)
    na0, nc0 = int((ca.C >= THR).sum()), int((cc.C >= THR).sum())
    axG = fig.add_subplot(gs[2, 0:2])
    head(axG, "G", "Number of bottlenecks, calibrated model and kcat draws")
    ndr = max(len(na), len(nc))
    for j, (vals, col, v0) in enumerate(((na, BLUE, na0), (nc, ORANGE, nc0))):
        off = -0.19 if j == 0 else 0.19
        for i, v in enumerate([v0] + vals):
            vv = max(v, 0.3)
            vbar(axG, i + off, vv, col, 9.0 if ndr < 16 else 5.0, b=0.3)
            val(axG, i + off, vv, f"{v}", dx=0, dy=6, ha="center", va="bottom", fontsize=7.8)
    axG.set_yscale("log")
    axG.set_ylim(0.3, 300)
    axG.set_xticks(range(ndr + 1))
    axG.set_xticklabels(["Calibrated"] + [str(i + 1) for i in range(ndr)], fontsize=8.5)
    axG.set_xlim(-0.6, ndr + 0.6)
    axG.set_xlabel("kcat perturbation draw", labelpad=2)
    axG.set_ylabel("Bottleneck reactions", labelpad=2)
    style(axG)
    legend(axG, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper right")
    print("bottlenecks per draw WT-Air", na, "WT-CO2", nc, flush=True)

    top_cal = cc.sort_values("C", ascending=False).reaction.iloc[0]
    nrun2 = rb_co2.run.nunique()
    freq2 = (rb_co2[rb_co2.C >= THR].groupby("reaction").run.nunique() / nrun2)
    tops = [g_.loc[g_.C.idxmax(), "reaction"] for _, g_ in rb_co2.groupby("run", sort=False)]
    n_first = sum(t == top_cal for t in tops)
    order = freq2.sort_values(ascending=False).index.tolist()
    if top_cal not in order[:12]:
        order = order[:11] + [top_cal]
    order = order[:12]
    axH = fig.add_subplot(gs[2, 2])
    head(axH, "H", f"WT-CO$_2$ frequency, {nrun2} draws")
    n = len(order)
    for i, r in enumerate(order):
        y = n - 1 - i
        v = 100 * freq2.get(r, 0.0)
        hbar(axH, y, max(v, 0.6), ORANGE if r == top_cal else BLUE, 7.5)
        val(axH, max(v, 0.6), y, f"{v:.0f}", va="center", ha="left", fontsize=8.2)
    axH.set_yticks(range(n))
    axH.set_yticklabels([rname(r) for r in order][::-1], fontsize=7.8)
    axH.set_ylim(-0.7, n - 0.3)
    axH.set_xlim(0, 118)
    axH.set_xlabel("Draws as a bottleneck (%)", labelpad=2)
    style(axH, left=False)
    print(f"{top_cal} is the top bottleneck in {n_first} of {nrun2} WT-CO2 draws; "
          f"it is a bottleneck in {100 * freq2.get(top_cal, 0):.0f}% of draws", flush=True)

finish(fig, out_path)
print(f"all done in {time.time() - t_start:.0f} s")