"""
Leave-out test for the iPiccel952 ecGEM calibration.

Question
    Do the kcat corrections found for one condition predict growth in the other condition when only the measured protein
    pool of that condition is imposed?

Design
    For each target condition (WT-Air, WT-CO2) the growth rate is maximized under the measured protein pool of that
    condition, using the condition's own upper bounds and TurNuP protein costs, with three versions of the calibration.
      none   no calibration (original upper bounds and original protein costs)
      other  calibration found for the other condition (leave-out)
      self   calibration found for the condition itself (reproduces the fitted growth)
    A calibration is the set of upper bounds relaxed in Stage 1 (these reactions also leave the pool accounting) and the
    cumulative kcat fold applied to each reaction in Stage 2 (protein cost divided by the fold). Both are read from the
    files written by relaxation_algorithm.py.

Usage
    python leave_out_test.py ecGEMs.zip leave_out_results.xlsx
Requires cobra, pandas and openpyxl. Gurobi is used when it is available, otherwise the default cobra solver.
"""
import os
import sys
import tempfile
import zipfile

import cobra
import pandas as pd

ZIP_PATH = "ecGEMs.zip"
OUT_PATH = "leave_out_results.xlsx"
zip_path = sys.argv[1] if len(sys.argv) > 1 else ZIP_PATH
out_path = sys.argv[2] if len(sys.argv) > 2 else OUT_PATH

BIOMASS = "Biomass_Picre_Auto"
RELAXED_UB = 1000.0
PROTEIN_COST_DIVISOR = 1000.0

COND = {
    "WT-Air": dict(dir="WT_Air", tag="Wt_Air", fix="iterative_kcat_fix_results_WtAir.xlsx",
                   growth=0.028, pool=0.077759),
    "WT-CO2": dict(dir="Wt_CO2", tag="Wt_CO2", fix="iterative_kcat_fix_results_Wt_CO2.xlsx",
                   growth=0.324941, pool=0.091287),
}

# Locate the condition folders inside the zip file
tmp = tempfile.mkdtemp(prefix="ecgems_")
with zipfile.ZipFile(zip_path) as zf:
    zf.extractall(tmp)
zip_dir = tmp
for root, dirs, _ in os.walk(tmp):
    if "WT_Air" in dirs:
        zip_dir = root
        break


def ra_dir(name):
    return f"{zip_dir}/{COND[name]['dir']}/relaxation_algorithm"


def read_ub(path):
    """Upper bounds from the two column UB files. The original files carry a header row and the revised ones do not,
    so non-numeric values are dropped."""
    df = pd.read_excel(path, header=None, usecols=[0, 1], names=["Reaction", "Upper_Bound"])
    df["Upper_Bound"] = pd.to_numeric(df["Upper_Bound"], errors="coerce")
    df = df.dropna(subset=["Upper_Bound"])
    return {str(r): float(u) for r, u in zip(df["Reaction"], df["Upper_Bound"])}


def read_calibration(name):
    """Relaxed reactions (Stage 1) and cumulative kcat folds (Stage 2) written by the relaxation algorithm."""
    c, ra = COND[name], ra_dir(name)
    boosted = pd.read_excel(f"{ra}/{c['fix']}", sheet_name="boosted_reactions")
    fold = {}
    for rid, f in zip(boosted["Reaction"].astype(str), boosted["Fold_applied"].astype(float)):
        fold[rid] = fold.get(rid, 1.0) * f
    ub_orig = read_ub(f"{ra}/UB_{c['tag']}.xlsx")
    ub_rev = read_ub(f"{ra}/UB_{c['tag']}_revised.xlsx")
    relaxed = {rid for rid, u in ub_rev.items() if u == RELAXED_UB and ub_orig.get(rid) != RELAXED_UB}
    print(f"{name} calibration  relaxed bounds {len(relaxed)}  boosted reactions {len(fold)}")
    return dict(relaxed=relaxed, fold=fold)


MODELS = {}


def load_model(name):
    if name not in MODELS:
        m = cobra.io.read_sbml_model(f"{ra_dir(name)}/irrev_iPiccel952.xml")
        try:
            m.solver = "gurobi"
        except Exception:
            pass
        MODELS[name] = m
    return MODELS[name]


def predict_growth(target, calib):
    """Maximum growth of the target condition under its measured pool, with the given calibration (None = none)."""
    c, ra = COND[target], ra_dir(target)
    model = load_model(target)
    ub = read_ub(f"{ra}/UB_{c['tag']}.xlsx")
    cost_df = pd.read_excel(f"{ra}/kcat_mw_{c['tag']}.xlsx", sheet_name="reaction_protein_cost")
    relaxed = calib["relaxed"] if calib else set()
    fold = calib["fold"] if calib else {}
    with model as m:
        rxns = {r.id: r for r in m.reactions}
        for rid, u in ub.items():
            if rid in rxns:
                rxns[rid].upper_bound = max(0.0, u)
        for rid in relaxed:
            if rid in rxns:
                rxns[rid].upper_bound = RELAXED_UB
        cost = {}
        for rid, pc in zip(cost_df["Reaction"].astype(str), cost_df["Protein_cost"]):
            if rid in rxns and pd.notna(pc) and rid not in relaxed:
                cost[rid] = float(pc) / fold.get(rid, 1.0)
        expr = sum((cost[rid] / PROTEIN_COST_DIVISOR) * rxns[rid].flux_expression for rid in cost)
        m.add_cons_vars(m.problem.Constraint(expr, lb=0, ub=c["pool"], name="protein_capacity"))
        m.objective = BIOMASS
        sol = m.optimize()
        return float(sol.objective_value) if sol.status == "optimal" else 0.0


cal = {n: read_calibration(n) for n in COND}
other = {"WT-Air": "WT-CO2", "WT-CO2": "WT-Air"}

rows = []
for tgt in COND:
    rows.append(dict(target=tgt, calibration="measured", growth=COND[tgt]["growth"]))
    rows.append(dict(target=tgt, calibration="none", growth=predict_growth(tgt, None)))
    rows.append(dict(target=tgt, calibration="other", growth=predict_growth(tgt, cal[other[tgt]])))
    rows.append(dict(target=tgt, calibration="self", growth=predict_growth(tgt, cal[tgt])))

res = pd.DataFrame(rows)
res.to_excel(out_path, sheet_name="leave_out", index=False)
print(res.to_string(index=False))
print("saved", out_path)