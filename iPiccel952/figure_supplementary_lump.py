"""
Supplementary Figures S1 to S5 for the eciPiccel952 manuscript, all from one script.

Usage
    python make_supp_figs.py ecGEMs.zip protein_abundance_g_per_gCDW.xlsx iPiccel952.xml iPceleri.fasta out_dir [DRAWS]

    ecGEMs.zip                          zip with the WT_Air and Wt_CO2 folders
    protein_abundance_g_per_gCDW.xlsx   iBAQ proteomics (g/gCDW), columns Protein, Wt_Air, Wt_CO2
    iPiccel952.xml                      original SBML model (pathway groups, formulas, genes)
    iPceleri.fasta                      protein sequences (molecular weights of urea carboxylase and allophanate hydrolase)
    out_dir                             folder for FigS1.png to FigS5.png (heavy results are cached in out_dir/supp_work)
    DRAWS                               kcat perturbation draws for the robustness panels (default 8)
    SEED                                random seed of the draws (default 1). Draws of several seeds can be run in parallel
                                        (add robust_only as the last argument) and are combined when the figures are drawn.

Figure to figure map
    FigS1 supports Fig. 1  model content and quality (pathways, gene associations, biomass, balance, connectivity)
    FigS2 supports Fig. 2  why the gene-association share is low and how the comparison changes without fatty acid reactions
    FigS3 supports Fig. 3  protein pool calibration (relaxation algorithm rounds, boosted reactions, kcat range, pool allocation)
    FigS4 supports Fig. 4  sensitivity of the nitrogen and carbon coupling result
    FigS5 supports Fig. 5  sensitivity of the bottleneck analysis (threshold, step, pool size, single relaxations, kcat error)

The first run takes about one hour on two CPUs (control scans, kcat perturbation draws, single-reaction scans).
Cached files are reused, so a second run only redraws the figures. Delete a file in supp_work to recompute it.
Requires cobra, matplotlib, numpy, pandas, scipy and openpyxl.
"""
import glob
import json
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
if len(args) < 5:
    sys.exit(__doc__)
zip_path, prot_path, orig_xml, fasta_path, out_dir = args[:5]
NDRAW = int(args[5]) if len(args) > 5 else 8
SEED = int(args[6]) if len(args) > 6 else 1
ROBUST_ONLY = len(args) > 7 and args[7] == "robust_only"
work = os.path.join(out_dir, "supp_work")
os.makedirs(work, exist_ok=True)

BIO = "Biomass_Picre_Auto"
THR = 0.01
SIGMA = np.log(2.0)
COND = {"air": ("WT_Air", "Wt_Air", 0.077759, "WT-Air"), "co2": ("Wt_CO2", "Wt_CO2", 0.091287, "WT-CO2")}
EXP_GROWTH = {"air": 0.028, "co2": 0.324941}
ANNOT = {"URCB": "GENE4066869", "ALPH": "GENE4066871"}

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
print("reading the original model", flush=True)
orig = cobra.io.read_sbml_model(orig_xml)
names = {r.id: r.name for r in orig.reactions}


def rname(r, n=30):
    s_ = names.get(r) or names.get(base(r)) or r
    s_ = s_.replace("_", " ")
    return s_ if len(s_) <= n else s_[: n - 1] + "."


# ---------------------------------------------------------------------------
# Heavy computations (cached)
# ---------------------------------------------------------------------------
def build(cond):
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


def robustness(ndraw, seed=1, step=1.1):
    """kcat error analysis. Every draw multiplies each kcat by exp(eps), eps ~ N(0, sigma^2) with sigma = ln 2
    (two-fold typical error), which changes the protein costs and the kcat-derived bounds together, and recomputes the
    control coefficients of the reactions that carry flux. The pool is not recalibrated, so growth of a draw is reported."""
    path = f"{work}/robust_{seed}.csv"
    done = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame(columns=["seed", "draw", "reaction", "growth0", "C"])
    k0 = int(done.draw.max()) + 1 if len(done) else 0
    if k0 >= ndraw:
        return
    m, ids, applied, pc, con = build("air")
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
        print(f"robustness seed {seed} draw {k + 1}/{ndraw}  growth {g0:.5f}  active reactions {len(act_)}", flush=True)
        pd.DataFrame(out, columns=["seed", "draw", "reaction", "growth0", "C"]).to_csv(path, index=False)
        for r, c in base_cost.items():
            con.set_linear_coefficients({ids[r].forward_variable: c})
            if r in applied:
                ids[r].upper_bound = applied[r]


def single_scan(mode):
    """Relax the bound of one reaction at a time (WT-Air) and record the relative growth gain.
    nopool_ub       original bounds and kcat files, no protein pool
    pool_ub         calibrated files with the pool, bound raised to 1000
    pool_ub_cost    as pool_ub and the protein cost of the reaction divided by 100"""
    path = f"{work}/scan_{mode}.csv"
    if os.path.exists(path):
        return
    d = f"{root}/WT_Air"
    m = cobra.io.read_sbml_model(f"{d}/irrev_iPiccel952.xml")
    pooled = mode != "nopool_ub"
    ubf = "UB_Wt_Air_revised.xlsx" if pooled else "relaxation_algorithm/UB_Wt_Air.xlsx"
    kcf = "kcat_mw_Wt_Air_revised.xlsx" if pooled else "relaxation_algorithm/kcat_mw_Wt_Air.xlsx"
    ub = pd.read_excel(f"{d}/{ubf}", usecols=[0, 1], names=["r", "u"])
    ids = {r.id: r for r in m.reactions}
    applied = {}
    for r, u in zip(ub.r, ub.u):
        r = str(r)
        if r in ids and pd.notna(u):
            ids[r].upper_bound = max(0.0, float(u))
            applied[r] = float(u)
    cost = pd.read_excel(f"{d}/{kcf}", sheet_name="reaction_protein_cost")
    pc = {str(r): float(c) for r, c in zip(cost.Reaction, cost.Protein_cost) if str(r) in ids and pd.notna(c)}
    m.objective = BIO
    if pooled:
        con = m.problem.Constraint(sum((c / 1000.0) * ids[r].flux_expression for r, c in pc.items()), lb=0, ub=0.077759, name="pool")
        m.add_cons_vars(con)
        m.solver.update()
    g0 = m.optimize().objective_value
    print(f"single-reaction scan {mode}  baseline {g0:.5f}", flush=True)
    rows = []
    for i, r in enumerate([r for r, u in applied.items() if u < 1000 and r != BIO]):
        rx = ids[r]
        old = rx.upper_bound
        rx.upper_bound = 1000.0
        oldc = None
        if mode == "pool_ub_cost" and r in pc:
            oldc = con.get_linear_coefficients([rx.forward_variable])[rx.forward_variable]
            con.set_linear_coefficients({rx.forward_variable: oldc / 100.0})
        s = m.optimize()
        g = s.objective_value if s.status == "optimal" else np.nan
        rx.upper_bound = old
        if oldc is not None:
            con.set_linear_coefficients({rx.forward_variable: oldc})
        rows.append((r, old, g, (g - g0) / g0))
        if i % 300 == 0:
            print(f"   {i}", flush=True)
    pd.DataFrame(rows, columns=["reaction", "ub", "growth", "rel_gain"]).to_csv(path, index=False)


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


# nitrogen and carbon coupling model (Fig. 4)
seqs = {}
cur = None
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
for cond, (folder, tag, _, _) in COND.items():
    kc = pd.read_excel(f"{root}/{folder}/relaxation_algorithm/kcat_mw_{tag}.xlsx", sheet_name="reaction_protein_cost")
    CAP[cond] = {rid: float(kc.loc[kc.Reaction == rid, "Kcat"].iloc[0]) * pa.loc[g, tag] / mw(g) * 1000.0
                 for rid, g in ANNOT.items()}


def coupling_model(cond, nadh=True, orig_ub=False):
    folder, tag, _, _ = COND[cond]
    m = cobra.io.read_sbml_model(f"{root}/{folder}/irrev_iPiccel952.xml")
    ubf = f"UB_{tag}.xlsx" if orig_ub else f"UB_{tag}_revised.xlsx"
    ub = pd.read_excel(f"{root}/{folder}/relaxation_algorithm/{ubf}", header=None, usecols=[0, 1], names=["r", "u"])
    ids = {r.id for r in m.reactions}
    for r, u in zip(ub.r, ub.u):
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
    if nadh:
        nadph, nadp, nadh_, nad = [m.metabolites.get_by_id(x) for x in ("nadph_h", "nadp_h", "nadh_h", "nad_h")]
        m.reactions.get_by_id("GLUS(nadph)").add_metabolites({nadph: 1, nadp: -1, nadh_: -1, nad: 1})
    m.objective = BIO
    return m


MULT = [0.5, 1, 2, 5, 10, None]


def coupling_sensitivity():
    path = f"{work}/coupling_sens.json"
    if os.path.exists(path):
        return
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
        # relax every active bounded reaction ten-fold, one at a time
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
        print(f"coupling sensitivity {cond} done, baseline {g0:.5f}", flush=True)
    json.dump(res, open(path, "w"))


# seed the cache with what the main pipeline needs
if ROBUST_ONLY:
    robustness(NDRAW, SEED)
    sys.exit(0)
for c_ in COND:
    control_coefficients(c_)
control_coefficients("air", 1.05)
control_coefficients("air", 1.2)
for mode in ("pool_ub", "nopool_ub", "pool_ub_cost"):
    single_scan(mode)
pool_scaling()
coupling_sensitivity()
if NDRAW > 0:
    robustness(NDRAW, SEED)
print(f"computations ready after {time.time() - t_start:.0f} s", flush=True)

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
BLUE, ORANGE, GREEN, LGREEN, PINK, PURPLE = "#2a78d6", "#eb6834", "#0b8a0b", "#8fd18f", "#e87ba4", "#9b6bd3"


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


def new_fig(w, h, nr, nc, **kw):
    fig = plt.figure(figsize=(w, h), dpi=100, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.06, wspace=0.12, hspace=0.14, rect=(0.0, 0.0, 1.0, 0.955))
    gs = fig.add_gridspec(nr, nc, **{k: v for k, v in kw.items() if k in ("height_ratios", "width_ratios")})
    return fig, gs


HEADS = []


def head(ax, letter, title):
    HEADS.append((ax, letter, title))


def finish(fig, name):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for ax, letter, title in HEADS:
        bb = ax.get_tightbbox(rend).transformed(fig.transFigure.inverted())
        p = ax.get_position()
        x0 = max(0.006, bb.x0 - 0.004)
        fig.text(x0, p.y1 + 0.012, letter, fontsize=9, fontweight="bold", va="bottom", ha="left")
        fig.text(x0 + 0.021, p.y1 + 0.013, title, fontsize=6.6, va="bottom", ha="left")
    HEADS.clear()
    fig.savefig(f"{out_dir}/{name}.png", dpi=400, facecolor="white")
    plt.close(fig)
    print("saved", f"{out_dir}/{name}.png", flush=True)


# ---------------------------------------------------------------------------
# Model statistics from the final SBML
# ---------------------------------------------------------------------------
M = orig
grp = {}
for g_ in M.groups:
    for mem in g_.members:
        grp.setdefault(mem.id, []).append(g_.name)
sub_of = {k: v[0] for k, v in grp.items()}
sub_path = os.path.join(work, "subsys_full.json")
FA = {r for r, s_ in sub_of.items() if s_ == "Fatty acid metabolism"}

is_bio = lambda r: r.id == BIO
is_ex = lambda r: r.boundary
is_tr = lambda r: (not r.boundary) and (not is_bio(r)) and len(r.compartments) > 1


import xml.etree.ElementTree as ET
_dec = lambda i: re.sub(r"__(\d+)__", lambda mo: chr(int(mo.group(1))), i[2:] if i.startswith("R_") else i)
REV_ATTR = {_dec(e.get("id")): e.get("reversible") == "true" for _, e in ET.iterparse(orig_xml) if e.tag.endswith("}reaction")}


def model_stats(rxns):
    rxns = list(rxns)
    genes = {g.id for r in rxns for g in r.genes}
    mets = {x.id for r in rxns for x in r.metabolites}
    n = len(rxns)
    return dict(rxn=n, genes=len(genes), met=len(mets), gpr=sum(bool(r.gene_reaction_rule) for r in rxns),
                transport=sum(is_tr(r) for r in rxns), reversible=sum(REV_ATTR.get(r.id, r.lower_bound < 0 < r.upper_bound) for r in rxns))


THIS = model_stats(M.reactions)
NOFA = model_stats([r for r in M.reactions if r.id not in FA])
print("iPiccel952", THIS, "\nwithout fatty acid reactions", NOFA, flush=True)

# ===========================================================================
# FigS1  model content and quality (supports Fig. 1)
# ===========================================================================
fig, gs = new_fig(7.0, 5.9, 2, 3, hspace=0.62, wspace=1.05)

# A pathways
cnt = pd.Series({g_.name: len(g_.members) for g_ in M.groups}).sort_values(ascending=False)
top = cnt.head(14)
axA = fig.add_subplot(gs[0, 2])
head(axA, "C", "Reactions per pathway")
for i, (k, v) in enumerate(top.items()):
    y = len(top) - 1 - i
    hbar(axA, y, np.log10(v), ORANGE if k == "Fatty acid metabolism" else BLUE, 4.6)
    val(axA, np.log10(v), y, f"{v:,}", va="center", ha="left", fontsize=4.8)
axA.set_yticks(range(len(top)))
axA.set_yticklabels([short(k) for k in top.index][::-1], fontsize=4.9)
axA.set_ylim(-0.7, len(top) - 0.3)
axA.set_xlim(0, 3.9)
axA.set_xticks([0, 1, 2, 3])
axA.set_xticklabels(["1", "10", "100", "1000"])
axA.set_xlabel("Reactions (log scale)", labelpad=2)
style(axA, left=False)

# B gene associations by class
cls = {"Metabolic": [], "Transport": [], "Exchange and demand": [], "Biomass": []}
for r in M.reactions:
    k = "Biomass" if is_bio(r) else "Exchange and demand" if is_ex(r) else "Transport" if is_tr(r) else "Metabolic"
    cls[k].append(bool(r.gene_reaction_rule))
axB = fig.add_subplot(gs[1, 0])
head(axB, "D", "Reactions with a gene association")
labs = list(cls)
pct = [100 * np.mean(v) if v else 0 for v in cls.values()]
hpanel(axB, labs, pct, [BLUE, GREEN, GREY, GREY], fmt="{:.0f}", xmax=118, lw=7.0, extra=[f"%  (n = {len(v):,})" for v in cls.values()])
axB.set_xlabel("Share of reactions (%)", labelpad=2)

# C biomass precursors
b = M.reactions.get_by_id(BIO)
bm = sorted(((abs(c) * x.formula_weight, x.name) for x, c in b.metabolites.items() if c < 0 and x.formula and x.name not in ("ATP", "H2O")), reverse=True)
topb = bm[:12]
axC = fig.add_subplot(gs[0, 0])
head(axC, "A", "Largest biomass components (no ATP, H$_2$O)")
hpanel(axC, [n_[:24] for _, n_ in topb], [v for v, _ in topb], [BLUE] * len(topb), fmt="{:.0f}", xmax=1.25 * topb[0][0], lw=4.6, fs=4.9)
axC.set_xlabel("Mass in biomass (mg/gDW)", labelpad=2)

# D residual mass imbalance in the final model
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
axD = fig.add_subplot(gs[0, 1])
head(axD, "B", "Residual imbalance by element")
ks = list(imb)
for i, k in enumerate(ks):
    vbar(axD, i, imb[k], BLUE, 8.0)
    val(axD, i, imb[k], f"{imb[k]}", dx=0, dy=4, ha="center", va="bottom")
axD.set_xticks(range(len(ks)))
axD.set_xticklabels(ks, fontsize=5.0, rotation=35, ha="right")
axD.set_ylim(0, 1.2 * max(max(imb.values()), 1))
axD.set_xlim(-0.6, len(ks) - 0.4)
axD.set_ylabel(f"Reactions ({tot_imb} imbalanced)", labelpad=2)
style(axD)
print("imbalanced reactions in the final model", tot_imb, imb, flush=True)

# E genes per reaction and reactions per gene
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

# F connectivity hubs
deg = sorted(((len(x.reactions), x.name, x.compartment) for x in M.metabolites if x.name not in ("H+", "H2O")), reverse=True)[:12]
axF = fig.add_subplot(gs[1, 2])
head(axF, "F", "Most connected metabolites")
hpanel(axF, [f"{n_[:20]} [{c}]" for _, n_, c in deg], [d_ for d_, _, _ in deg], [BLUE] * 12, fmt="{:.0f}", xmax=1.25 * deg[0][0], lw=4.6, fs=4.9)
axF.set_xlabel("Reactions per metabolite", labelpad=2)
finish(fig, "FigS1")

# ===========================================================================
# FigS2  gene-association share and the comparison without fatty acid reactions (supports Fig. 2)
# ===========================================================================
PUBLISHED = {
    "iCO1515": dict(genes=1515, rxn=2516, gpr=1975, transport=306, reversible=196, met=3517),
    "iCre1355": dict(genes=1460, rxn=2394, gpr=1924, transport=449, reversible=723, met=1845),
    "iCZ843": dict(genes=843, rxn=2294, gpr=1838, transport=378, reversible=626, met=1770),
    "iAJ526": dict(genes=531, rxn=1455, gpr=856, transport=317, reversible=561, met=1236),
}
ROWS = [("iPiccel952", THIS, GREEN), ("iPiccel952\nno fatty acid", NOFA, LGREEN)] + [(k, v, BLUE) for k, v in PUBLISHED.items()]
fig, gs = new_fig(7.0, 5.6, 2, 6, hspace=0.62, wspace=2.6)

axA = fig.add_subplot(gs[:, 0:2])
head(axA, "A", "Gene association by pathway")
rows = []
for g_ in M.groups:
    mem = [x for x in g_.members]
    if not mem:
        continue
    rows.append((g_.name, len(mem), 100 * sum(bool(x.gene_reaction_rule) for x in mem) / len(mem)))
rows = sorted(rows, key=lambda t: -t[1])[:16]
for i, (nm, n_, p_) in enumerate(rows):
    y = len(rows) - 1 - i
    hbar(axA, y, p_, ORANGE if nm == "Fatty acid metabolism" else BLUE, 4.4)
    val(axA, p_, y, f"{p_:.0f}", va="center", ha="left", fontsize=4.8)
axA.set_yticks(range(len(rows)))
axA.set_yticklabels([f"{short(nm)} ({n_:,})" for nm, n_, _ in rows][::-1], fontsize=4.8)
axA.set_ylim(-0.7, len(rows) - 0.3)
axA.set_xlim(0, 118)
axA.set_xlabel("Reactions with a gene association (%)", labelpad=2)
style(axA, left=False)

METRICS = [("Reactions with GPR (%)", lambda d_: 100 * d_["gpr"] / d_["rxn"], "{:.0f}"),
           ("Transport (%)", lambda d_: 100 * d_["transport"] / d_["rxn"], "{:.0f}"),
           ("Reversible (%)", lambda d_: 100 * d_["reversible"] / d_["rxn"], "{:.0f}"),
           ("Reactions per gene", lambda d_: d_["rxn"] / d_["genes"], "{:.1f}"),
           ("Reactions per metabolite", lambda d_: d_["rxn"] / d_["met"], "{:.2f}"),
           ("Genes per 100 reactions", lambda d_: 100 * d_["genes"] / d_["rxn"], "{:.0f}")]
pos = [(0, 2), (0, 4), (1, 2), (1, 4)]
for k, (title, f_, fmt) in enumerate(METRICS[:4]):
    ax = fig.add_subplot(gs[pos[k][0], pos[k][1]:pos[k][1] + 2])
    head(ax, "BCDE"[k], title)
    vals = [f_(d_) for _, d_, _ in ROWS]
    labs = [n_ if (pos[k][1] == 2) else "" for n_, _, _ in ROWS]
    n = len(ROWS)
    for i, (v, (_, _, c)) in enumerate(zip(vals, ROWS)):
        y = n - 1 - i
        hbar(ax, y, v, c, 5.2)
        val(ax, v, y, fmt.format(v), va="center", ha="left", fontsize=5.0)
    ax.set_yticks(range(n))
    ax.set_yticklabels([n_ for n_, _, _ in ROWS][::-1] if pos[k][1] == 2 else [""] * n, fontsize=4.9)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, 1.3 * max(vals))
    style(ax, left=False)
_ent = [([GREEN], "iPiccel952"), ([LGREEN], "iPiccel952 without fatty acid reactions"), ([BLUE], "Published algal models")]
fig.legend([tuple(dot(c) for c in cols) for cols, _ in _ent], [t for _, t in _ent],
           handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)}, loc="outside lower center", ncol=3, fontsize=5.6, frameon=False)
finish(fig, "FigS2")

# ===========================================================================
# FigS3  protein pool calibration (supports Fig. 3)
# ===========================================================================
rx_air = pd.read_excel(f"{root}/WT_Air/relaxation_algorithm/iterative_kcat_fix_results_WtAir.xlsx", sheet_name=None)
rx_co2 = pd.read_excel(f"{root}/Wt_CO2/relaxation_algorithm/iterative_kcat_fix_results_Wt_CO2.xlsx", sheet_name=None)
fig, gs = new_fig(7.0, 5.9, 2, 3, hspace=0.62, wspace=1.05)


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
    ax.annotate(f"{n} rounds, {int(d_['summary'].N_reactions_boosted_step2.iloc[0])} boosted\ntarget pool {tgt:.3f} g/gDW", (0.97, 0.95),
                xycoords="axes fraction", ha="right", va="top", fontsize=5.4)
    style(ax)
    head(ax, letter, label)


axA = fig.add_subplot(gs[0, 0])
rounds_panel(axA, "A", rx_air, "WT-Air calibration")
axB = fig.add_subplot(gs[0, 1])
rounds_panel(axB, "B", rx_co2, "WT-CO$_2$ calibration")
legend(axB, [([ORANGE], "Above target"), ([GREEN], "Target reached")], loc="upper right", bbox_to_anchor=(1.0, 0.80))

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

ca, cc = pd.read_csv(f"{work}/control_air.csv"), pd.read_csv(f"{work}/control_co2.csv")
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
finish(fig, "FigS3")

# ===========================================================================
# FigS4  sensitivity of the nitrogen and carbon coupling result (supports Fig. 4)
# ===========================================================================
CS = json.load(open(f"{work}/coupling_sens.json"))
fig, gs = new_fig(7.0, 5.7, 2, 3)
XL = ["0.5x", "1x", "2x", "5x", "10x", "Free"]


def gfmt(v):
    return f"{v:.2f}" if v >= 0.01 else f"{v:.1e}".replace("e-0", "e-")


for k, (key, ttl) in enumerate((("URCB", "Urea carboxylase capacity"), ("ALPH", "Allophanate hydrolase capacity"),
                                ("both", "Both capacities"))):
    ax = fig.add_subplot(gs[0, k])
    head(ax, "ABC"[k], ttl)
    for i in range(len(MULT)):
        for j, (cond, col) in enumerate((("air", BLUE), ("co2", ORANGE))):
            v = CS[cond][key][i]
            x = i + (-0.19 if j == 0 else 0.19)
            vbar(ax, x, v, col, 4.6, b=1e-6)
            val(ax, x, v, gfmt(v), dx=0, dy=2, ha="center", va="bottom", fontsize=4.2, rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(3e-5, 3.0)
    ax.set_xlim(-0.6, len(MULT) - 0.4)
    ax.set_xticks(range(len(MULT)))
    ax.set_xticklabels(XL, fontsize=5.2)
    ax.set_xlabel("Capacity relative to the measured value", labelpad=2)
    ax.set_ylabel("Coupling-model growth (1/h)", labelpad=2)
    style(ax)
    if k == 0:
        legend(ax, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left")

# D model assumptions
axD = fig.add_subplot(gs[1, 0])
head(axD, "D", "Model assumptions")
grp_ = [("nadh", "Used"), ("nadph", "NADPH\nglutamate\nsynthase"), ("orig_ub", "Original\nbound file")]
for i, (key, lab) in enumerate(grp_):
    for j, (cond, col) in enumerate((("air", BLUE), ("co2", ORANGE))):
        v = CS[cond][key]
        x = i + (-0.19 if j == 0 else 0.19)
        vbar(axD, x, v, col, 7.0, b=1e-6)
        val(axD, x, v, gfmt(v), dx=0, dy=2, ha="center", va="bottom", fontsize=4.6, rotation=90)
axD.set_yscale("log")
axD.set_ylim(3e-5, 0.05)
axD.set_xticks(range(len(grp_)))
axD.set_xticklabels([l for _, l in grp_], fontsize=5.0)
axD.set_xlim(-0.6, len(grp_) - 0.4)
axD.set_ylabel("Coupling-model growth (1/h)", labelpad=2)
style(axD)
legend(axD, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left")

# E fold change of growth, WT-CO2 over WT-Air
axE = fig.add_subplot(gs[1, 1])
head(axE, "E", "Growth fold, WT-CO$_2$ / WT-Air")
fold_items = [("Measured\ncapacity", CS["co2"]["both"][1] / CS["air"]["both"][1]),
              ("Urea carb.\n10x", CS["co2"]["URCB"][4] / CS["air"]["URCB"][4]),
              ("Both\n10x", CS["co2"]["both"][4] / CS["air"]["both"][4]),
              ("Both\nfree", CS["co2"]["both"][5] / CS["air"]["both"][5]),
              ("Experiment", EXP_GROWTH["co2"] / EXP_GROWTH["air"])]
for i, (lab, v) in enumerate(fold_items):
    vbar(axE, i, v, GREY if lab == "Experiment" else BLUE, 9.0)
    val(axE, i, v, f"{v:.1f}", dx=0, dy=3, ha="center", va="bottom", fontsize=5.0)
axE.set_xticks(range(len(fold_items)))
axE.set_xticklabels([l for l, _ in fold_items], fontsize=4.8)
axE.set_xlim(-0.6, len(fold_items) - 0.4)
axE.set_ylim(0, 1.25 * max(v for _, v in fold_items))
axE.set_ylabel("Fold", labelpad=2)
style(axE)

# F single bound relaxations
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
        hbar(axF, y + off, v, col, 5.0)
        val(axF, v, y + off, f"{v:.0f}%", va="center", ha="left", fontsize=5.0)
axF.set_yticks(range(n))
axF.set_yticklabels([rname(r, 24) for r in topr][::-1], fontsize=5.2)
axF.set_ylim(-0.7, n - 0.3)
axF.set_xlim(0, 100 * max(CS[c]["gains"].get(topr[0], 0) for c in COND) * 1.35)
axF.set_xlabel("Growth gain (%)", labelpad=2)
style(axF, left=False)
axF.annotate(f"The other {nact - len(topr)} active bounded reactions\nchange growth by less than 1%", (0.02, 0.97), xycoords="axes fraction",
             ha="left", va="top", fontsize=5.0)
legend(axF, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="lower right")
finish(fig, "FigS4")

# ===========================================================================
# FigS5  sensitivity of the bottleneck analysis (supports Fig. 5)
# ===========================================================================
fig, gs = new_fig(7.0, 5.7, 2, 3)
ca12 = pd.read_csv(f"{work}/control_air_1.20.csv")
ca05 = pd.read_csv(f"{work}/control_air_1.05.csv")

# A threshold
taus = [0.001, 0.003, 0.01, 0.03, 0.1]
axA = fig.add_subplot(gs[0, 0])
head(axA, "A", "Bottlenecks against threshold")
for i, tau in enumerate(taus):
    for j, (d_, col) in enumerate(((ca, BLUE), (cc, ORANGE))):
        v = int((d_.C >= tau).sum())
        x = i + (-0.19 if j == 0 else 0.19)
        vbar(axA, x, max(v, 0.3), col, 6.0, b=0.3)
        val(axA, x, max(v, 0.3), f"{v}", dx=0, dy=3, ha="center", va="bottom", fontsize=4.8)
axA.set_yscale("log")
axA.set_ylim(0.3, 300)
axA.set_xticks(range(len(taus)))
axA.set_xticklabels([f"{t:g}" for t in taus], fontsize=5.4)
axA.set_xlim(-0.6, len(taus) - 0.4)
axA.set_xlabel("Control coefficient threshold", labelpad=2)
axA.set_ylabel("Bottleneck reactions", labelpad=2)
style(axA)
legend(axA, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper right")

# B single relaxations
modes = [("nopool_ub", "No protein pool,\nbound relaxed"), ("pool_ub", "Protein pool,\nbound relaxed"),
         ("pool_ub_cost", "Protein pool, bound\nand cost relaxed")]
axB = fig.add_subplot(gs[0, 1])
head(axB, "B", "One reaction relaxed at a time")
gmax, n10 = [], []
for mname, _ in modes:
    d_ = pd.read_csv(f"{work}/scan_{mname}.csv")
    gmax.append(100 * d_.rel_gain.max())
    n10.append(int((d_.rel_gain >= 0.1).sum()))
for i, v in enumerate(gmax):
    y = len(modes) - 1 - i
    hbar(axB, y, max(v, 0.0), BLUE, 7.0)
    val(axB, max(v, 0.0), y, f"max {max(v, 0.0):.1f}%, {n10[i]} with 10% or more", va="center", ha="left", fontsize=5.0)
axB.set_yticks(range(len(modes)))
axB.set_yticklabels([m_ for _, m_ in modes][::-1], fontsize=5.0)
axB.set_ylim(-0.7, len(modes) - 0.3)
axB.set_xlim(0, 140)
axB.set_xlabel("Largest growth gain (%)", labelpad=2)
style(axB, left=False)

# C pool size
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
        vbar(axC, x, v, col, 5.2)
        val(axC, x, v, f"{v:.2f}", dx=0, dy=2, ha="center", va="bottom", fontsize=4.4, rotation=90)
    print(cond, "growth without the pool relative to calibrated", round(d_[1e6] / ref, 2), flush=True)
axC.set_xticks(range(len(KS)))
axC.set_xticklabels([f"{k:g}x" for k in KS], fontsize=5.2)
axC.set_xlim(-0.6, len(KS) - 0.4)
axC.set_ylim(0, 1.4 * rel_max)
axC.set_xlabel("Pool relative to calibrated", labelpad=2)
axC.set_ylabel("Growth relative to calibrated", labelpad=2)
style(axC)
legend(axC, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper left")

# D step size
sets = {"5%": set(ca05.reaction[ca05.C >= THR]), "10%": set(ca.reaction[ca.C >= THR]), "20%": set(ca12.reaction[ca12.C >= THR])}
ref = sets["10%"]
axD = fig.add_subplot(gs[1, 0])
head(axD, "D", "Step size of the kcat increase")
for i, (k, s_) in enumerate(sets.items()):
    vbar(axD, i, len(s_), GREY, 9.0)
    vbar(axD, i, len(s_ & ref), BLUE, 9.0)
    val(axD, i, len(s_), f"{len(s_)}\n({len(s_ & ref)} shared)", dx=0, dy=3, ha="center", va="bottom", fontsize=4.8)
axD.set_xticks(range(3))
axD.set_xticklabels(list(sets), fontsize=5.4)
axD.set_xlim(-0.6, 2.6)
axD.set_ylim(0, 1.3 * max(len(s_) for s_ in sets.values()))
axD.set_xlabel("kcat increase", labelpad=2)
axD.set_ylabel("Bottleneck reactions, WT-Air", labelpad=2)
style(axD)
legend(axD, [([BLUE], "Shared with 10%"), ([GREY], "Other")], loc="upper right")
rho_ = stats.spearmanr(ca.set_index("reaction").C, ca05.set_index("reaction").C.reindex(ca.reaction).values, nan_policy="omit")[0]
print("Spearman of control coefficients, 10% against 5%", round(rho_, 3), flush=True)

# E and F robustness
files = sorted(glob.glob(f"{work}/robust_*.csv"))
if files:
    rb = pd.concat([pd.read_csv(f) for f in files])
    rb["run"] = rb.seed.astype(str) + "_" + rb.draw.astype(str)
    nrun = rb.run.nunique()
    freq = (rb[rb.C >= THR].groupby("reaction").run.nunique() / nrun).sort_values(ascending=False)
    bot = list(ca.sort_values("C", ascending=False).reaction.head(12))
    axE = fig.add_subplot(gs[1, 1])
    head(axE, "E", f"Bottleneck frequency, {nrun} kcat draws")
    n = len(bot)
    for i, r in enumerate(bot):
        y = n - 1 - i
        v = 100 * freq.get(r, 0.0)
        hbar(axE, y, max(v, 0.6), BLUE, 4.6)
        val(axE, max(v, 0.6), y, f"{v:.0f}", va="center", ha="left", fontsize=4.8)
    axE.set_yticks(range(n))
    axE.set_yticklabels([rname(r, 24) for r in bot][::-1], fontsize=4.8)
    axE.set_ylim(-0.7, n - 0.3)
    axE.set_xlim(0, 118)
    axE.set_xlabel("Draws in which the reaction is a bottleneck (%)", labelpad=2)
    style(axE, left=False)

    axF = fig.add_subplot(gs[1, 2])
    head(axF, "F", "Rank agreement with the calibrated C")
    c0 = ca.set_index("reaction").C
    rhos = []
    for run, g_ in rb.groupby("run"):
        s_ = g_.set_index("reaction").C
        common = s_.index.intersection(c0.index)
        rhos.append(stats.spearmanr(c0[common], s_[common])[0])
    for i, v in enumerate(rhos):
        vbar(axF, i, v, BLUE, 8.0 if len(rhos) < 12 else 5.0)
        val(axF, i, v, f"{v:.2f}", dx=0, dy=3 if v >= 0 else -3, ha="center", va="bottom" if v >= 0 else "top", fontsize=4.8)
    axF.set_xticks(range(len(rhos)))
    axF.set_xticklabels([str(i + 1) for i in range(len(rhos))], fontsize=5.4)
    axF.set_xlim(-0.6, len(rhos) - 0.4)
    axF.set_ylim(-0.25, 1.2 * max(max(rhos), 0.1))
    axF.axhline(0, color=AXIS, lw=0.6, zorder=2)
    axF.set_xlabel("kcat perturbation draw", labelpad=2)
    axF.set_ylabel("Spearman correlation of C", labelpad=2)
    style(axF)
    print("robustness rank correlation", [round(v, 2) for v in rhos], flush=True)
finish(fig, "FigS5")
print(f"all done in {time.time() - t_start:.0f} s")