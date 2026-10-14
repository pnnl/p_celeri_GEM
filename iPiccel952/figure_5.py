"""
Fig. 5 of the eciPiccel952 manuscript in one script.
Kcat control analysis of the enzyme-constrained model, pathway assignment, robustness to kcat error,
proteomics tests and the figure.

Usage
    python fig5_all.py ecGEMs.zip protein_abundance_g_per_gCDW.xlsx iPiccel952.xml Fig5_bottlenecks.png [DRAWS]

    ecGEMs.zip                          zip with WT_Air and Wt_CO2 folders (irrev model, revised UB and kcat_mw files)
    protein_abundance_g_per_gCDW.xlsx   iBAQ proteomics (g/gCDW), columns Protein, Wt_Air, Wt_CO2
    iPiccel952.xml                      original SBML model with the pathway groups
    DRAWS                               number of kcat perturbation draws for the robustness whiskers (default 6, 0 skips)

Intermediate results are stored in ./fig5_work and reused on a second run, so only the figure is redrawn.
Delete a file in that folder to recompute it.

Method in short
    Control coefficient of reaction i  C_i = (mu_i - mu_0) / (delta mu_0), mu_i being the growth rate after raising kcat_i
    by delta = 10 percent (its kcat-derived bound rises and its protein cost falls by the same factor), with the
    protein pool constraint in place. Bottleneck reactions have C_i >= 0.01. Pathway control is the sum of C_i over the
    pathway, tested by permutation. Robustness repeats the analysis after lognormal kcat errors (sigma = ln 2, two-fold typical error).
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
from matplotlib.colors import to_rgb
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from scipy import stats

warnings.filterwarnings("ignore")
args = sys.argv[1:]
zip_path, prot_path, orig_xml, out_path = args[:4] if len(args) >= 4 else (
    "ecGEMs.zip", "protein_abundance_g_per_gCDW.xlsx", "iPiccel952.xml", "Fig5_bottlenecks.png")
NDRAW = int(args[4]) if len(args) > 4 else 6
work = "fig5_work"
os.makedirs(work, exist_ok=True)
THR = 0.01
NPERM = 5000
STEP = 1.1
SIGMA = np.log(2.0)
BIO = "Biomass_Picre_Auto"
COND = {"air": ("WT_Air", "Wt_Air", 0.077759), "co2": ("Wt_CO2", "Wt_CO2", 0.091287)}

tmp = tempfile.mkdtemp(prefix="ecgems_")
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(tmp)
root = tmp
for r_, dirs, _ in os.walk(tmp):
    if "WT_Air" in dirs:
        root = r_
        break


# ---------------------------------------------------------------------------
# Step 1. Pathway of every reaction (SBML groups, then rule-based assignment for ungrouped reactions)
# ---------------------------------------------------------------------------
RULES = [
    (r'cytochrome|photosystem|plastoquinone|plastocyanin|ferredoxin|thylakoid', 'Photosynthesis'),
    (r'ATP synthase|transhydrogenase|diphosphatase|NADH dehydrogenase|ubiquinol|ubiquinone|complex', 'Oxidative phosphorylation'),
    (r'transport|antiporter|symport|transporter|uniport|diffusion|\bvia\b|exchange', 'TRANSPORT'),
    (r'purine|nucleoside|nucleotide|pyrimidine|uridine|uracil|adenyl|guanyl|\bAMP\b|cytidine|thymidine|ribonucleotide|urate|allant|dihydropyrimidine|dcamp|hypoxanthine|inosine|xanthine|\bGMP\b|\bIMP\b|nucleoside', 'Nucleotide metabolism'),
    (r'enolase|phosphoglycerate|pyruvate kinase|hexokinase|phosphofructo|aldolase|triosephosphate|phosphopyruvate hydratase|glyceraldehyde|glucose-6-phosphate|phosphoglucomutase|2-o-phosphotransferase', 'Glycolysis / Gluconeogenesis'),
    (r'carboxylase|ribulose|phosphoribulokinase|sedoheptulose|pyruvate, phosphate dikinase', 'Carbon fixation'),
    (r'ketol-acid|acetolactate|isopropylmalate|dihydroxy-acid|methylpropanoyl|methylbutanoyl|branched|isovaleryl|valine|leucine|isoleucine', 'Valine, leucine and isoleucine metabolism'),
    (r'pyruvate dehydrogenase|aldehyde dehydrogenase|alcohol dehydrogenase|acetaldehyde|ethanol|butanol|lactate|acetate|acetyl-CoA synthetase', 'Pyruvate metabolism'),
    (r'glycolate|glyoxylate|phosphoglycolate', 'Glyoxylate metabolism'),
    (r'glutamyl-trna|porphyr|chlorophyll|heme|protoporph|ferrochelatase', 'Porphyrin and chlorophyll metabolism'),
    (r'glutathione|glutaredoxin|thioredoxin', 'Glutathione metabolism'),
    (r'anthranilate|chorismate|shikimate|tryptophan|phenylalanine|tyrosine', 'Phenylalanine, tyrosine and tryptophan biosynthesis'),
    (r'glycerol|glycerophospho|phosphatidyl|lipid', 'Glycerolipid metabolism'),
    (r'pantoth|folate|riboflavin|thiamin|biotin|NAD\b|lipoyl|lipoate|cofactor|coenzyme|pyridox', 'Cofactor and Prosthetic Group Biosynthesis'),
    (r'acyl|fatty|enoyl|hydroxyacyl|oxoacyl|acetyl-coa c-acetyl|long-chain|very-long|acp\b|sulfoquinovos', 'Fatty acid metabolism'),
    (r'prephenate|chorismate|anthranil', 'Phenylalanine, tyrosine and tryptophan biosynthesis'),
    (r'lycopene|carotene|zeaxanthin|violaxanthin|phytoene|xanthophyll', 'Carotenoid biosynthesis'),
    (r'tocopherol|homogentisate|geranylgeranyl|hydroxy-3-methylbut', 'Biosynthesis of steroids'),
    (r'isocitrate|oxoglutarate dehydrogenase|succinate|fumarase|malate dehydrogenase|citrate synthase|aconitase', 'TCA cycle'),
    (r"thymidylate|5'-nucleotidase", 'Nucleotide metabolism'),
    (r'glycerate|hydroxyisobutyrate', 'Glycine, serine and threonine metabolism'),
    (r'oxidoreductase|peroxidase|dioxygenase', 'Miscellaneous'),
    (r'catalase|peroxide', 'Miscellaneous'),
    (r'dihydroorotate|aspartate|asparag|glutam|glutamine|ammonia|urea|arginine|ornithine|proline', 'Glutamate and arginine metabolism'),
    (r'serine|glycine|threonine', 'Glycine, serine and threonine metabolism'),
    (r'cysteine|methionine|homocysteine|sulfate|sulfur', 'Cysteine and methionine metabolism'),
    (r'methylerythritol|terpen|isopentenyl|prenyl|carotenoid|steroid|sterol|squalene|farnesyl', 'Biosynthesis of steroids'),
]


def assign_rule(rx):
    nm = rx.name or ""
    for pat, s_ in RULES:
        if re.search(pat, nm, re.I):
            if s_ == "TRANSPORT":
                if re.search("mitochond", nm, re.I) or rx.id.endswith(("tm", "_tm")):
                    return "Transport, mitochondrial"
                if re.search("chloroplast|thylakoid", nm, re.I) or re.search(r"(_th|_thr|th)$", rx.id):
                    return "Transport, chloroplast"
                return "Transport, extracellular"
            return s_
    return None


sub_path = f"{work}/subsys_full.json"
if not os.path.exists(sub_path):
    orig = cobra.io.read_sbml_model(orig_xml)
    grp = {}
    for g_ in orig.groups:
        for mem in g_.members:
            grp.setdefault(mem.id, []).append(g_.name)
    out_ = {}
    n_rule = n_ungrouped = 0
    for rx in orig.reactions:
        if grp.get(rx.id):
            out_[rx.id] = grp[rx.id][0]
            continue
        n_ungrouped += 1
        a_ = assign_rule(rx)
        if a_:
            out_[rx.id] = a_
            n_rule += 1
    print("pathways | ungrouped reactions", n_ungrouped, "assigned by rules", n_rule, "still missing", n_ungrouped - n_rule)
    json.dump(out_, open(sub_path, "w"))
sub = json.load(open(sub_path))
base = lambda r: re.sub(r"_(FWD|REV)$", "", r)


# ---------------------------------------------------------------------------
# Step 2. Enzyme-constrained model and kcat control coefficients
# ---------------------------------------------------------------------------
def build(cond):
    folder, tag, pool = COND[cond]
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


def control_coefficients(cond):
    path = f"{work}/control_{cond}.csv"
    if os.path.exists(path):
        return
    m, ids, applied, pc, con = build(cond)
    s0 = m.optimize()
    base_g, fl = s0.objective_value, s0.fluxes
    print(cond, "baseline growth", round(base_g, 5), flush=True)
    rows, t0 = [], time.time()
    for i, r in enumerate(pc):
        rx = ids[r]
        v = rx.forward_variable
        oldc, oldu = con.get_linear_coefficients([v])[v], rx.upper_bound
        con.set_linear_coefficients({v: oldc / STEP})
        if oldu < 1000:
            rx.upper_bound = min(1000.0, oldu * STEP)
        s = m.optimize()
        g = s.objective_value if s.status == "optimal" else np.nan
        con.set_linear_coefficients({v: oldc})
        rx.upper_bound = oldu
        rows.append((r, pc[r], applied.get(r, np.nan), fl[r], (g - base_g) / base_g / (STEP - 1)))
        if i % 300 == 0:
            print(f"  {cond} control {i}/{len(pc)}  {time.time() - t0:.0f} s", flush=True)
    df_ = pd.DataFrame(rows, columns=["reaction", "cost", "ub", "flux", "C"])
    df_["pool_use"] = df_.cost * df_.flux / 1000.0
    df_.to_csv(path, index=False)


def robustness(ndraw, seed=1, step=STEP):
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


control_coefficients("air")
control_coefficients("co2")
if NDRAW > 0:
    robustness(NDRAW)

# ---------------------------------------------------------------------------
# Step 3. Analysis and figure
# ---------------------------------------------------------------------------
model = cobra.io.read_sbml_model(f"{root}/WT_Air/irrev_iPiccel952.xml")
pa = pd.read_excel(prot_path).iloc[:, :3].set_index("Protein")

def enzyme_abundance(rid, col):
    """GPR evaluation, isoenzymes (or) add and complex subunits (and) take the minimum. Result in g/gCDW."""
    rule = model.reactions.get_by_id(rid).gene_reaction_rule
    if not rule:
        return np.nan

    def split(s, op):
        depth, parts, cur = 0, [], ""
        for t in re.split(r"(\(|\)|\s+and\s+|\s+or\s+)", s):
            if t == "(":
                depth += 1
                cur += t
            elif t == ")":
                depth -= 1
                cur += t
            elif depth == 0 and t.strip() == op:
                parts.append(cur)
                cur = ""
            else:
                cur += t
        parts.append(cur)
        return parts

    def ev(e):
        e = e.strip()
        ors = split(e, "or")
        if len(ors) > 1:
            return sum(ev(x) for x in ors)
        ands = split(e, "and")
        if len(ands) > 1:
            return min(ev(x) for x in ands)
        if e.startswith("(") and e.endswith(")"):
            return ev(e[1:-1])
        g = e.replace("G_", "")
        return float(pa.loc[g, col]) if g in pa.index else 0.0

    return ev(rule)


def load(cond):
    d = pd.read_csv(f"{work}/control_{cond}.csv")
    d["subsystem"] = [sub.get(r) or sub.get(base(r)) or "Unassigned" for r in d.reaction]
    d["ab_air"] = [enzyme_abundance(r, "Wt_Air") for r in d.reaction]
    d["ab_co2"] = [enzyme_abundance(r, "Wt_CO2") for r in d.reaction]
    return d


air, co2 = load("air"), load("co2")
names = {r.id: r.name for r in model.reactions}
a = air[air.subsystem != "Unassigned"].copy()
print("reactions", len(air), "unassigned", int((air.subsystem == "Unassigned").sum()), "sum C", round(air.C.sum(), 3),
      "bottlenecks", int((air.C >= THR).sum()), "| WT-CO2 sum C", round(co2.C.sum(), 3))

# pathway control with permutation test (reactions of pathways with at least 5 members)
rng = np.random.default_rng(1)
Cv, N = a.C.values, len(a)
rows = []
for s, g in a.groupby("subsystem"):
    K = len(g)
    if K < 5:
        continue
    obs = g.C.sum()
    perm = np.array([Cv[rng.choice(N, K, replace=False)].sum() for _ in range(NPERM)])
    rows.append((s, K, obs, obs / Cv.sum(), (1 + (perm >= obs).sum()) / (NPERM + 1)))
E = pd.DataFrame(rows, columns=["sub", "K", "C", "share", "p"])
o = np.argsort(E.p.values)
q = E.p.values[o] * len(E) / (np.arange(len(E)) + 1)
q = np.minimum.accumulate(q[::-1])[::-1]
E["q"] = np.nan
E.loc[E.index[o], "q"] = np.minimum(q, 1)
E = E.sort_values("C", ascending=False).reset_index(drop=True)
print(E.head(10).to_string())

# robustness to kcat error
rob_files = sorted(glob.glob(f"{work}/robust_*.csv"))
rob = None
if rob_files:
    rb = pd.concat([pd.read_csv(f) for f in rob_files])
    rb["subsystem"] = [sub.get(r) or sub.get(base(r)) or "Unassigned" for r in rb.reaction]
    rb["run"] = rb.seed.astype(str) + "_" + rb.draw.astype(str)
    ps = rb[rb.subsystem != "Unassigned"].pivot_table(index="run", columns="subsystem", values="C", aggfunc="sum").fillna(0)
    rob = ps.div(ps.sum(axis=1), axis=0)
    print("perturbation draws", len(rob))

# proteomics tests
act = a[(a.flux > 1e-9) & a.ab_air.notna() & (a.ab_air > 0)].copy()
bm = act.C >= THR
p_bn = stats.mannwhitneyu(act.ab_air[bm], act.ab_air[~bm]).pvalue
rho, p_rho = stats.spearmanr(act.pool_use, act.ab_air)
print("bottleneck median abundance", act.ab_air[bm].median() * 1000, "others", act.ab_air[~bm].median() * 1000, "p", p_bn)
print("predicted demand vs measured, Spearman", rho, p_rho, "n", len(act))
share = a.groupby("subsystem")[["ab_air", "ab_co2"]].sum()
share = share / share.sum()

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
INK, AXIS, GREY, TRACK = "#000000", "#7a7a7a", "#b9b8b0", "#eef0f3"
BLUE, ORANGE = "#2a78d6", "#eb6834"
PALETTE = ["#eda100", "#1baf7a", "#e87ba4", "#2a78d6", "#9b6bd3", "#eb6834", "#14a3b8", "#8a9a1b"]
top_paths = list(E.sort_values("C", ascending=False)["sub"].head(8))
PCOL = {p: PALETTE[i] for i, p in enumerate(top_paths)}
pcol = lambda p: PCOL.get(p, GREY)


def tint(c, al):
    r_, g_, b_ = to_rgb(c)
    return (1 - al) + al * r_, (1 - al) + al * g_, (1 - al) + al * b_


def short(p):
    return (p.replace("Glycolysis / Gluconeogenesis", "Glycolysis").replace("Glycine, serine and threonine metabolism", "Gly, Ser, Thr")
            .replace("metabolism", "").replace("Transport, ", "Transport ").replace("Pentose phosphate pathway", "Pentose phosphate")
            .replace("Phenylalanine, tyrosine and tryptophan biosynthesis", "Aromatic amino acids")
            .replace("Oxidative phosphorylation", "Oxid. phosphorylation").strip())


fig = plt.figure(figsize=(7.0, 6.8), dpi=100)
gs = fig.add_gridspec(3, 12, height_ratios=[1.45, 1.25, 1.0], left=0.075, right=0.985, top=0.985, bottom=0.07,
                      hspace=0.36, wspace=3.4)
HEADS = []


def style(ax):
    ax.minorticks_off()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK, pad=1.5)


def dot(c):
    return Line2D([], [], marker="o", ls="", ms=5.5, mfc=c, mec=c)


def legend(ax, entries, **kw):
    hs = [tuple(dot(c) for c in cols) for cols, _ in entries]
    ax.legend(hs, [t for _, t in entries], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
              handletextpad=0.4, labelspacing=0.35, borderpad=0, fontsize=5.6, **kw)


def hbar(ax, y, v, col, lw=6.0, base_=0.0):
    ax.plot([base_, v], [y, y], color=col, lw=lw, solid_capstyle="round", zorder=3)


def vbar(ax, x, v, col, lw=7.0, base_=0.0):
    ax.plot([x, x], [base_, v], color=col, lw=lw, solid_capstyle="round", zorder=3)


def val(ax, x, y, s_, dx=7, dy=0, **kw):
    ax.annotate(s_, (x, y), xytext=(dx, dy), textcoords="offset points", annotation_clip=False, **kw)


# ---------------------------------------------------------------------------
# A. Algorithm
# ---------------------------------------------------------------------------
axA = fig.add_subplot(gs[0, :])
axA.set_xlim(0, 100)
axA.set_ylim(4.0, 67.2)
axA.axis("off")
HEADS.append((axA, "A", None, True))
SERIF = ["Liberation Serif", "STIXGeneral", "DejaVu Serif"]
FS = 6.4


def txt(x, y, s_, bold=False, color=INK, fs=FS):
    axA.text(x, y, s_, fontsize=fs, family=SERIF, fontweight="bold" if bold else "normal", color=color, ha="left",
             va="center")


def code_line(col_x, y, num, indent, s_):
    txt(col_x, y, f"{num}", fs=5.6)
    txt(col_x + 3.0 + 2.6 * indent, y, s_)


axA.plot([0.5, 99.5], [66.5, 66.5], color=INK, lw=0.9)
txt(1.0, 64.2, "Algorithm 2", bold=True)
txt(13.0, 64.2, "Bottleneck detection by kcat control analysis of an enzyme-constrained model")
axA.plot([0.5, 99.5], [62.2, 62.2], color=INK, lw=0.5)
txt(1.0, 59.8, "Input", bold=True)
txt(7.0, 59.8, r"ecGEM ($S\,v=0$, bounds $UB_i$, costs $c_i$, pool $P^{*}$),  enzyme-constrained reactions $E$,  step $\delta=0.1$,  "
    r"threshold $\tau=0.01$,  pathway map $\Pi$,  proteomics $A$")
txt(1.0, 55.8, "Output", bold=True)
txt(7.0, 55.8, r"bottleneck set $B$,  pathway control $C_p$ with $q_p$,  robustness frequency $f_i$")
axA.plot([0.5, 99.5], [53.4, 53.4], color=INK, lw=0.5)
axA.plot([50, 50], [6.4, 52.8], color="#d8d7cf", lw=0.5)
XL, XR = 1.0, 51.5
txt(XL, 51.0, "Stage 1", bold=True)
txt(XL + 8.0, 51.0, "control coefficients", color=INK)
st = 5.3
L1 = [
    (0, r"$\mu_0\leftarrow\max\ v_{biomass}$  s.t.  $S\,v=0$,  $v_i\leq UB_i$,  $\sum_i c_i v_i\leq P^{*}$"),
    (0, r"$\mathbf{for}$ $i\in E$ $\mathbf{do}$"),
    (1, r"$k_{cat,i}\leftarrow(1+\delta)\,k_{cat,i}$   ($UB_i\leftarrow(1+\delta)UB_i$,  $c_i\leftarrow c_i/(1+\delta)$)"),
    (1, r"$\mu_i\leftarrow$ re-solve the model"),
    (1, r"$C_i\leftarrow(\mu_i-\mu_0)/(\delta\,\mu_0)$"),
    (1, r"restore $k_{cat,i}$"),
    (0, r"$\mathbf{end\ for}$"),
    (0, r"$B\leftarrow\{\,i\in E:\ C_i\geq\tau\,\}$"),
]
for k, (ind, t) in enumerate(L1):
    code_line(XL, 46.4 - k * st, k + 1, ind, t)
txt(XR, 51.0, "Stage 2", bold=True)
txt(XR + 8.0, 51.0, "pathways, robustness and proteomics", color=INK)
L2 = [
    (0, r"$\mathbf{for}$ each pathway $p\in\Pi$ $\mathbf{do}$   $C_p\leftarrow\sum_{i\in p}C_i$"),
    (1, r"$q_p\leftarrow$ FDR of the permutation test of $C_p$  (5000 random sets)"),
    (0, r"$\mathbf{for}$ $r=1,\ldots,R$ $\mathbf{do}$   $k_{cat,i}\leftarrow k_{cat,i}\,e^{\varepsilon_i}$,  $\varepsilon_i\sim\mathcal{N}(0,\sigma^{2})$"),
    (1, r"repeat Stage 1 and record the pathway shares $C_p/\sum C_p$"),
    (0, r"$f_i\leftarrow$ fraction of repeats with $i\in B$"),
    (0, r"$\mathbf{test}$ measured $A_i$ of $i\in B$ against $i\notin B$  (Mann-Whitney)"),
    (0, r"$\mathbf{test}$ predicted $c_i v_i$ against measured $A_i$  (Spearman)"),
    (0, r"$\mathbf{return}$ $B$, $C_p$, $q_p$, $f_i$"),
]
for k, (ind, t) in enumerate(L2):
    code_line(XR, 46.4 - k * st, 9 + k, ind, t)
axA.plot([0.5, 99.5], [5.4, 5.4], color=INK, lw=0.9)

# ---------------------------------------------------------------------------
# B. Top reactions
# ---------------------------------------------------------------------------
NTOP = 14
tb = a.sort_values("C", ascending=False).head(NTOP)
axB = fig.add_subplot(gs[1, 0:4])
HEADS.append((axB, "B", "Control coefficient, WT-Air"))
for i, (_, r) in enumerate(tb.iterrows()):
    y = NTOP - 1 - i
    hbar(axB, y, r.C, pcol(r.subsystem), lw=5.0)
    val(axB, r.C, y, f"{r.C:.3f}", va="center", ha="left", fontsize=5.0)
DISPLAY = {
    "PPC": "PEP carboxylase", "NanoG1457": "Purine nucleoside phosph.", "THRAL_FWD": "Threonine ammonia-lyase",
    "NanoG0606": "ATP synthase", "ATP_3h_tm": "ADP/ATP carrier", "AM6PT": "Mannose 6-P transferase",
    "3PG_pi_thr": "3-PGA/Pi translocator", "SERth_FWD": "Serine transporter", "PGKh_FWD": "Phosphoglycerate kinase",
    "RPEh_REV": "Ribulose-5-P epimerase", "TKT2h_REV": "Transketolase 2", "TKT1h_REV": "Transketolase 1",
    "ATUDm": "UDP phosphotransferase", "ARGSL": "Argininosuccinate lyase", "FBA3hi": "SBP aldolase",
    "GBDP_FWD": "ppGpp hydrolase",
}
labs = [DISPLAY.get(r.reaction, names.get(r.reaction, r.reaction)[:26]) for _, r in tb.iterrows()]
axB.set_yticks(range(NTOP))
axB.set_yticklabels(labs[::-1], fontsize=5.0)
axB.set_xlim(0, 0.105)
axB.set_ylim(-0.7, NTOP - 0.3)
axB.set_xticks([0, 0.03, 0.06])
p_ = axB.get_position()
axB.set_position([p_.x0 + 0.075, p_.y0, p_.width - 0.075, p_.height])
style(axB)
axB.spines["left"].set_visible(False)
axB.tick_params(axis="y", length=0, pad=2)

# ---------------------------------------------------------------------------
# C. Pathway control share with robustness interval
# ---------------------------------------------------------------------------
NP = 9
tp = E.head(NP)
axC = fig.add_subplot(gs[1, 4:8])
HEADS.append((axC, "C", "Pathway share of control (%)"))
for i, (_, r) in enumerate(tp.iterrows()):
    y = NP - 1 - i
    hbar(axC, y, 100 * r.share, pcol(r["sub"]), lw=5.4)
    if rob is not None and r["sub"] in rob.columns:
        lo, hi = np.percentile(100 * rob[r["sub"]], [10, 90])
        axC.plot([lo, hi], [y + 0.0, y + 0.0], color=INK, lw=0.7, zorder=4)
        axC.plot([lo, lo], [y - 0.12, y + 0.12], color=INK, lw=0.7, zorder=4)
        axC.plot([hi, hi], [y - 0.12, y + 0.12], color=INK, lw=0.7, zorder=4)
        xr = max(hi, 100 * r.share)
    else:
        xr = 100 * r.share
    lab = f"{100 * r.share:.0f}"
    val(axC, xr, y, lab, va="center", ha="left", fontsize=5.0, dx=4)
axC.set_yticks(range(NP))
axC.set_yticklabels([short(s) for s in tp["sub"]][::-1], fontsize=5.2)
axC.set_xlim(0, 36)
axC.set_ylim(-0.7, NP - 0.3)
axC.set_xticks([0, 10, 20, 30])
style(axC)
axC.spines["left"].set_visible(False)
axC.tick_params(axis="y", length=0, pad=2)

# ---------------------------------------------------------------------------
# D. WT-Air against WT-CO2
# ---------------------------------------------------------------------------
sc = co2[co2.subsystem != "Unassigned"].groupby("subsystem").C.sum()
sa = a.groupby("subsystem").C.sum()
order = list(E["sub"].head(5)) + [p for p in sc.sort_values(ascending=False).index[:1] if p not in list(E["sub"].head(5))]
axD = fig.add_subplot(gs[1, 8:12])
HEADS.append((axD, "D", "Control per pathway"))
for i, p in enumerate(order):
    y = len(order) - 1 - i
    for ci, (v, col) in enumerate([(sa.get(p, 0.0), BLUE), (sc.get(p, 0.0), ORANGE)]):
        yy = y + (0.19 if ci == 0 else -0.19)
        hbar(axD, yy, max(v, 0.0015), col, lw=5.4)
        val(axD, max(v, 0.0015), yy, f"{max(v, 0.0):.2f}", va="center", ha="left", fontsize=5.0, dx=4)
axD.set_yticks(range(len(order)))
axD.set_yticklabels([short(s) for s in order][::-1], fontsize=5.2)
axD.set_xlim(0, 0.34)
axD.set_ylim(-0.7, len(order) - 0.3)
axD.set_xticks([0, 0.1, 0.2, 0.3])
style(axD)
axD.spines["left"].set_visible(False)
axD.tick_params(axis="y", length=0, pad=2)
legend(axD, [([BLUE], "WT-Air"), ([ORANGE], "WT-CO$_2$")], loc="upper right", bbox_to_anchor=(1.02, 1.03))

# ---------------------------------------------------------------------------
# E. Measured abundance of bottleneck enzymes
# ---------------------------------------------------------------------------
axE = fig.add_subplot(gs[2, 0:4])
HEADS.append((axE, "E", "Measured abundance (mg/gCDW)"))
YE = 0.003
for k, (m_, lab) in enumerate([(act.ab_air[bm], "Bottleneck"), (act.ab_air[~bm], "Other")]):
    v = 1000 * m_.median()
    vbar(axE, k, v, [BLUE, GREY][k], lw=15)
    lo, hi = 1000 * m_.quantile(0.25), 1000 * m_.quantile(0.75)
    axE.plot([k, k], [lo, hi], color=INK, lw=0.8, zorder=4)
    axE.plot([k - 0.06, k + 0.06], [lo, lo], color=INK, lw=0.8, zorder=4)
    axE.plot([k - 0.06, k + 0.06], [hi, hi], color=INK, lw=0.8, zorder=4)
    val(axE, k, hi, f"median {v:.2f}", dx=0, dy=5, ha="center", va="bottom", fontsize=5.4)
axE.text(0.5, 0.96, f"p = {p_bn:.3f}", transform=axE.transAxes, ha="center", va="top", fontsize=5.8)
axE.set_xlim(-0.6, 1.6)
axE.set_ylim(0, 1.5 * 1000 * act.ab_air[bm].quantile(0.75))
axE.set_xticks([0, 1])
axE.set_xticklabels([f"Bottleneck\n(n = {int(bm.sum())})", f"Other\n(n = {int((~bm).sum())})"])
style(axE)

# ---------------------------------------------------------------------------
# F. Predicted demand against measured abundance
# ---------------------------------------------------------------------------
axF = fig.add_subplot(gs[2, 4:8])
HEADS.append((axF, "F", "Predicted demand and measured abundance"))
act["bin"] = pd.qcut(act.pool_use.rank(method="first"), 5, labels=False)
meds = [1000 * act.ab_air[act.bin == b].median() for b in range(5)]
for b in range(5):
    vbar(axF, b, meds[b], tint(BLUE, 0.35 + 0.13 * b), lw=10)
    val(axF, b, meds[b], f"{meds[b]:.2f}", dx=0, dy=6, ha="center", va="bottom", fontsize=5.2)
axF.text(0.03, 0.96, f"Spearman r = {rho:.2f}\np = {p_rho:.0e}", transform=axF.transAxes, ha="left", va="top", fontsize=5.8)
axF.set_xlim(-0.6, 4.6)
axF.set_ylim(0, 1.4 * max(meds))
axF.set_xticks(range(5))
axF.set_xticklabels(["1", "2", "3", "4", "5"])
axF.set_xlabel("Predicted demand quintile (low to high)", labelpad=2, fontsize=5.8)
style(axF)

# ---------------------------------------------------------------------------
# G. Proteome share of the control pathways, WT-CO2 over WT-Air
# ---------------------------------------------------------------------------
axG = fig.add_subplot(gs[2, 8:12])
HEADS.append((axG, "G", "Pathway proteome, WT-CO$_2$ / WT-Air"))
gp = list(E["sub"].head(6))
fl_ = [share.loc[p, "ab_co2"] / share.loc[p, "ab_air"] for p in gp]
for i, (p, f) in enumerate(zip(gp, fl_)):
    y = len(gp) - 1 - i
    hbar(axG, y, f, pcol(p), lw=5.4)
    val(axG, f, y, f"{f:.2f}", va="center", ha="left", fontsize=5.0, dx=4)
axG.set_yticks(range(len(gp)))
axG.set_yticklabels([short(p) for p in gp][::-1], fontsize=5.2)
axG.set_xlim(0, 2.0)
axG.set_ylim(-0.7, len(gp) - 0.3)
axG.set_xticks([0, 0.5, 1.0, 1.5])
p_ = axG.get_position()
axG.set_position([p_.x0 + 0.03, p_.y0, p_.width - 0.03, p_.height])
style(axG)
axG.spines["left"].set_visible(False)
axG.tick_params(axis="y", length=0, pad=2)

# ---------------------------------------------------------------------------
for h in HEADS:
    ax, txt_, title = h[:3]
    wide = len(h) == 4 and h[3]
    p = ax.get_position()
    if wide:
        fig.text(0.008, p.y1 - 0.004, txt_, fontsize=9, fontweight="bold", va="top", ha="left")
        continue
    x0 = 0.008 if p.x0 < 0.2 else p.x0 - 0.062
    fig.text(x0, p.y1 + 0.012, txt_, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text(x0 + 0.022, p.y1 + 0.013, title, fontsize=6.6, va="bottom", ha="left")

fig.savefig(out_path, dpi=500, facecolor="white")
print("saved", out_path)