"""Build the enzyme-constrained P. celeri models (WT_Air, Wt_CO2) from the patched iPiccel954
(iPiccel954_patched.xml, made by patch_iPiccel954.py), recalibrate them with the unchanged PNNL relaxation
algorithm, and test the chloroplast PGM (GENE4067541) knockout.

Inputs (read only): /home/yasin/Claudio/GEM/p_celeri_GEM/ecGEMs.zip (Niaz's irreversible models, unrelaxed UB and
kcat/MW tables, relaxation_algorithm.py), /home/yasin/Claudio/GEM/p_celeri_GEM/iPceleri.fasta,
/home/yasin/Claudio/GEM/p_celeri_GEM/protein_abundance_g_per_gCDW.xlsx, /home/yasin/Claudio/GEM/iPiccel954_patched.xml.

Steps
  1. Irreversible split of the patched model with Niaz's convention (ChatGEM 48: reversible -> _FWD + _REV,
     reverse-only -> _FWD with flipped stoichiometry), keeping the patched model's bound magnitudes.
  2. Niaz's unrelaxed protein cost (kcat, MW) and upper bound (kcat x abundance / MW x 1000) rows are carried over
     by reaction AND direction (matched on stoichiometry with protons ignored), because several patches change a
     reaction's direction. Rows are never applied to a reaction the patch closed.
  3. Reactions that gained genes in the patch get a cost row: kcat = the gene's DeepEnzyme prediction
     (kcat_predictions.tsv from predict_enzyme_kcats.py, geometric mean over the reaction substrates run; the three
     DGATs in dgat_kcat_predictions.tsv), or the median of Niaz's TurNuP kcats where no prediction exists (genes
     without a resolvable substrate), MW from iPceleri.fasta of the gene used, gene = the most abundant measured gene
     for OR rules and the least abundant subunit for AND rules; their upper bound uses the same formula with the
     condition's measured abundance (left at the model bound when no gene was measured). An OR rule whose genes
     carry predicted kcats gets the sum over its measured isozymes (each at its own kcat), since the single-gene
     convention would discard the expressed paralogs (the three DGATs were the first case).
  4. Niaz's relaxation_algorithm.py is run unchanged on each condition (growth targets as in his runs, pool target = summed measured abundance of the model enzymes).
  5. Tests: wild type, PGM knockout (starch removed from biomass, starch synthase fixed at zero), and a urea-cap
     scan on each calibrated model.
Outputs (GEM/): kcat_substrate_requests.tsv (genes still at the median kcat, with their reaction substrates; with
--substrates-only the script stops after writing it), irrev_iPiccel954_patched_<cond>.xml, UB_<tag>_patched_revised.xlsx, kcat_mw_<tag>_patched_revised.xlsx,
relaxation_<tag>_patched_results.xlsx, ecGEM_patch_tests.tsv
Run:  uv run --project /home/yasin/tools/ChatGEM python -I /home/yasin/Claudio/GEM/patch_ecGEMs.py
"""
import os
import re
import shutil
import tempfile
import zipfile

import sys

import cobra
import numpy as np
import pandas as pd

GEM = "/home/yasin/Claudio/GEM"
ZIP = f"{GEM}/p_celeri_GEM/ecGEMs.zip"
FASTA = f"{GEM}/p_celeri_GEM/iPceleri.fasta"
ABUND = f"{GEM}/p_celeri_GEM/protein_abundance_g_per_gCDW.xlsx"
# Our own measured absolute proteome (Total Protein Approach on the integrated MS1 precursor areas; protein per cell and AFDW per cell from the amino-acid analysis;
# AFDW/CDW 0.85 assumed): g protein per g cell dry weight, 3 vessels per condition. The PNNL file above is an earlier export of the same data on another scaling
# (total protein 0.604 g/gCDW at high CO2 and 0.480 in air against 0.3145 and 0.4174 measured), so its abundances, and the bounds built on them, were too high.
MEASURED_ABUND = "/home/yasin/Claudio/Gen_Trans/final_analysis/GEM_share/protein_abundance_g_per_gCDW_TPA_current.csv"
MEASURED_GROUPS = "/home/yasin/Claudio/Gen_Trans/final_analysis/GEM_share/protein_group_mapping.csv"
_SUFFIX = re.compile(r"_(?:p0x?|p1x?|chloroplast|mitochondria)$")


def measured_abundance():
    """DataFrame(Protein = model gene id, Wt_CO2, Wt_Air) in g/gCDW: the mean of the three vessels, each protein group's abundance given to every member
    gene (the group stands for all of them; NA means not quantified and is never replaced by zero)."""
    new = pd.read_csv(MEASURED_ABUND).set_index("Protein")
    grp = pd.read_csv(MEASURED_GROUPS)
    out = {}
    for cond in ("Wt_CO2", "Wt_Air"):
        v = new[[f"{cond}_R{i}" for i in (1, 2, 3)]].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)
        d = {}
        for rep, members in zip(grp.Protein, grp["Group members"]):
            a = v.get(rep)
            if a is None or pd.isna(a):
                continue
            for mem in str(members).split(";"):
                d.setdefault(_SUFFIX.sub("", mem.strip()), a)
        out[cond] = d
    genes = sorted(set(out["Wt_CO2"]) | set(out["Wt_Air"]))
    return pd.DataFrame({"Protein": genes, "Wt_CO2": [out["Wt_CO2"].get(g) for g in genes], "Wt_Air": [out["Wt_Air"].get(g) for g in genes]})
PATCHED = f"{GEM}/iPiccel954_patched.xml"
BIO = "Biomass_Picre_Auto"
PGM_GENE = "GENE4067541"
CONDS = {"WT_Air": "Wt_Air", "Wt_CO2": "Wt_CO2"}
AA = {"A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388, "E": 129.1155, "Q": 128.1307,
      "G": 57.0519, "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766,
      "P": 97.1167, "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326}
cobra.Configuration().solver = "glpk"
# GLPK's default simplex cycles forever on a degenerate LP met in Niaz's relaxation (one solve ran for hours at full CPU); with presolve the
# same relaxation finishes in 9 s. The optimum is unchanged, so every model built or solved through this module uses presolve.
import cobra.core.model as _cm
_orig_optimize, _orig_slim = _cm.Model.optimize, _cm.Model.slim_optimize


def _presolved(fn):
    def wrapper(self, *a, **k):
        try:
            self.solver.configuration.presolve = True
        except Exception:  # noqa: BLE001
            pass
        return fn(self, *a, **k)
    return wrapper


_cm.Model.optimize, _cm.Model.slim_optimize = _presolved(_orig_optimize), _presolved(_orig_slim)
# Per-gene kcat predictions (h-1) that replace the median placeholder: DeepEnzyme (predict_enzyme_kcats.py) on every
# gene that otherwise carries the median, run on the reaction substrates listed in kcat_substrate_requests.tsv (written
# by this script with --substrates-only); the three DGATs came first (dgat_kcat_predictions.tsv). Per gene: geometric
# mean over the substrates run.
KCAT = {}
for _kp in (f"{GEM}/dgat_kcat_predictions.tsv", f"{GEM}/kcat_predictions.tsv"):
    if os.path.exists(_kp):
        _k = pd.read_csv(_kp, sep="\t")
        _k = _k[_k.gene.str.startswith("GENE") & (_k.kcat_per_h > 0)]
        _k["gene"] = _k.gene.str.split("_").str[0]
        KCAT.update(_k.groupby("gene").kcat_per_h.apply(lambda s: float(np.exp(np.log(s).mean()))).to_dict())
CURRENCY = {"h2o", "h", "atp", "adp", "amp", "pi", "ppi", "nad", "nadh", "nadp", "nadph", "coa", "o2", "co2", "hco3", "nh4",
            "h2o2", "fad", "fadh2", "gtp", "gdp", "ctp", "cdp", "utp", "udp", "ump", "gmp", "cmp", "so4", "so3", "na1", "k",
            "cl", "fe2", "fe3", "mg2", "ca2", "zn2", "mn2", "q8", "q8h2", "pq", "pqh2", "pcrd", "pcox", "fdxrd", "fdxox",
            "trdrd", "trdox", "ficytc", "focytc", "ficytb5", "focytb5", "ACP", "amet", "ahcys", "thf", "mlthf", "glu__L",
            "akg", "e_", "photon"}
SUBSTRATE_REQUESTS = f"{GEM}/kcat_substrate_requests.tsv"


def fasta_mw(path):
    seqs, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0].split("_")[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip()
    return {g: sum(AA.get(a, 0) for a in s.rstrip("*")) + 18.015 for g, s in seqs.items()}


def to_irrev(model):
    m = model.copy()
    add, drop = [], []
    for r in m.reactions:
        lb, ub = r.lower_bound, r.upper_bound
        if lb < 0 < ub:
            f = cobra.Reaction(f"{r.id}_FWD", name=f"{r.name} (fwd)", lower_bound=0, upper_bound=ub, subsystem=r.subsystem)
            b = cobra.Reaction(f"{r.id}_REV", name=f"{r.name} (rev)", lower_bound=0, upper_bound=-lb, subsystem=r.subsystem)
            f.add_metabolites(dict(r.metabolites))
            b.add_metabolites({x: -c for x, c in r.metabolites.items()})
            f.gene_reaction_rule = b.gene_reaction_rule = r.gene_reaction_rule
            add += [f, b]
            drop.append(r)
        elif lb < 0 and ub <= 0:
            f = cobra.Reaction(f"{r.id}_FWD", name=f"{r.name} (conv)", lower_bound=0, upper_bound=-lb, subsystem=r.subsystem)
            f.add_metabolites({x: -c for x, c in r.metabolites.items()})
            f.gene_reaction_rule = r.gene_reaction_rule
            add.append(f)
            drop.append(r)
    m.remove_reactions(drop)
    m.add_reactions(add)
    for r in m.reactions:
        r.lower_bound = max(0.0, r.lower_bound)
    m.objective = BIO
    return m


def base(rid):
    return re.sub(r"_(FWD|REV)$", "", rid)


def skey(r):
    return frozenset((x.id, round(c, 6)) for x, c in r.metabolites.items() if not x.id.startswith("h_"))


def orient_match(new, old_by_base):
    """Find Niaz's irreversible reaction with the same base ID and the same direction."""
    cands = old_by_base.get(base(new.id), [])
    k = skey(new)
    for o in cands:
        if skey(o) == k:
            return o
    heavy = [x.id for x, c in new.metabolites.items() if not x.id.startswith("h_")]
    for o in cands:   # stoichiometry changed by a patch (e.g. FNORh, FACOAL lumps): match direction on a shared metabolite
        om = {x.id: c for x, c in o.metabolites.items()}
        shared = [h for h in heavy if h in om]
        if shared and (om[shared[0]] > 0) == (new.metabolites[new.model.metabolites.get_by_id(shared[0])] > 0):
            return o
    return None


GPR_CORRECTED = {"BTNC", "ACCOAC", "ACCOAhi"}   # plastid acetyl-CoA carboxylase, gene rules corrected in patch_iPiccel954.py step N


def prepare(cond, tag, src, mw, ab, ab_old):
    ra = f"{src}/ecGEMs/{cond}/relaxation_algorithm"
    old = cobra.io.read_sbml_model(f"{ra}/irrev_iPiccel954.xml")
    old_by_base = {}
    for r in old.reactions:
        old_by_base.setdefault(base(r.id), []).append(r)
    rev = cobra.io.read_sbml_model(PATCHED)
    closed = {r.id for r in rev.reactions if r.bounds == (0, 0)}
    m = to_irrev(rev)
    # Niaz's relaxation relaxes bounds one by one until growth reaches the target; with the forced maintenance flux the model is infeasible
    # under the measured-abundance bounds (they sustain 0.008 of the 0.62 mmol ATP gDW-1 h-1 needed) and the algorithm then loops through
    # thousands of useless relaxations. It therefore calibrates without maintenance, and build() adds it back and relaxes what it needs.
    if "NGAM" in m.reactions:
        m.reactions.NGAM.lower_bound = 0.0
    ub_old = pd.read_excel(f"{ra}/UB_{tag}.xlsx")
    ub_old = dict(zip(ub_old.Reaction.astype(str), ub_old.Upper_Bound))
    xl = pd.ExcelFile(f"{ra}/kcat_mw_{tag}.xlsx")
    cost_old = xl.parse("reaction_protein_cost")
    pool = xl.parse("proteome_pool")
    cost_by = {str(r.Reaction): r for _, r in cost_old.iterrows()}
    kmed = float(cost_old.Kcat.median())
    A = dict(zip(ab.Protein, ab[tag]))
    A_old = dict(zip(ab_old.Protein, pd.to_numeric(ab_old[tag], errors="coerce")))
    # PNNL's upper bounds are exactly kcat x PNNL abundance / MW x 1000 (1178 of 1282 rows within 1 %): each carried-over bound is rescaled by the ratio of our measured
    # abundance to the PNNL abundance of the same gene; a gene that was not quantified keeps its bound times the median ratio of the condition.
    ratios = [A[g] / A_old[g] for g in A if g in A_old and pd.notna(A[g]) and pd.notna(A_old[g]) and A[g] > 0 and A_old[g] > 0]
    med_ratio = float(pd.Series(ratios).median())
    n_rescaled = n_median = 0
    ub_rows, cost_rows, n_new, n_pred = [], [], 0, 0
    pool_iso = set()   # measured isozymes whose capacity is summed into an OR-rule bound: they belong to the pool as well
    requests = []
    for r in m.reactions:
        if base(r.id) in closed:
            r.bounds = (0, 0)
            continue
        o = orient_match(r, old_by_base)
        if o is not None and base(r.id) in GPR_CORRECTED:   # step N rewrote the gene rule: the carried-over cost row names the wrong gene, so the reaction is treated as new
            o = None
        if o is not None and o.id in ub_old and pd.notna(ub_old[o.id]):
            ub = ub_old[o.id]
            c0 = cost_by.get(o.id)
            if c0 is not None and "GENE" in r.gene_reaction_rule:
                g0 = c0.Gene
                if g0 in A and pd.notna(A[g0]) and A[g0] > 0 and g0 in A_old and pd.notna(A_old[g0]) and A_old[g0] > 0:
                    ub = ub * A[g0] / A_old[g0]; n_rescaled += 1
                else:
                    ub = ub * med_ratio; n_median += 1
            ub_rows.append((r.id, ub))
        if o is not None and o.id in cost_by and "GENE" in r.gene_reaction_rule:
            c = cost_by[o.id]
            cost_rows.append(dict(Reaction=r.id, Gene=c.Gene, Kcat=c.Kcat, MW=c.MW, Protein_cost=c.Protein_cost))
        elif "GENE" in r.gene_reaction_rule and (o is None or o.id not in cost_by):
            genes = [g.id for g in r.genes if g.id in mw]
            if not genes:
                continue
            meas = {g: A[g] for g in genes if g in A and pd.notna(A[g]) and A[g] > 0}
            if " and " in r.gene_reaction_rule and meas:
                g = min(meas, key=meas.get)
            elif meas:
                g = max(meas, key=meas.get)
            else:
                g = sorted(genes)[0]
            kc = KCAT.get(g, kmed)
            if meas and " and " not in r.gene_reaction_rule and any(x in KCAT for x in meas):
                # predicted per-gene kcats available (DeepEnzyme): the OR rule's capacity is the sum over the measured
                # isozymes, each at its own kcat (median where none was predicted)
                cost_rows.append(dict(Reaction=r.id, Gene=g, Kcat=kc, MW=mw[g], Protein_cost=mw[g] / kc))
                ub_rows.append((r.id, sum(KCAT.get(x, kmed) * a / mw[x] for x, a in meas.items()) * 1000))
                pool_iso.update(meas)
                n_pred += 1
            else:
                cost_rows.append(dict(Reaction=r.id, Gene=g, Kcat=kc, MW=mw[g], Protein_cost=mw[g] / kc))
                if g in meas:
                    ub_rows.append((r.id, kc * meas[g] / mw[g] * 1000))
                if g in KCAT:
                    n_pred += 1
            # substrates DeepEnzyme would need for every gene of this rule still at the median kcat
            for x in genes:
                if x in KCAT:
                    continue
                for s, c in r.metabolites.items():
                    if c < 0 and re.sub(r"_[a-z]$", "", s.id) not in CURRENCY:
                        ch = s.annotation.get("chebi", "")
                        kg = s.annotation.get("kegg.compound", "")
                        requests.append(dict(condition=cond, reaction=r.id, gene=x, metabolite=s.id, name=s.name,
                                             formula=s.formula or "", charge=s.charge if s.charge is not None else "",
                                             kegg=kg[0] if isinstance(kg, list) else kg, chebi=ch[0] if isinstance(ch, list) else ch,
                                             inchikey=s.annotation.get("inchikey", "")))
            n_new += 1
    work = tempfile.mkdtemp()
    for f in os.listdir(ra):
        if f.endswith(".py"):
            shutil.copy(f"{ra}/{f}", work)
    cobra.io.write_sbml_model(m, f"{work}/irrev_iPiccel954.xml")
    pd.DataFrame(ub_rows, columns=["Reaction", "Upper_Bound"]).to_excel(f"{work}/UB_{tag}.xlsx", index=False)
    with pd.ExcelWriter(f"{work}/kcat_mw_{tag}.xlsx", engine="openpyxl") as w:
        pd.DataFrame(cost_rows).to_excel(w, sheet_name="reaction_protein_cost", index=False)
        # the pool is the summed measured abundance of the enzymes in the model (Niaz's definition), now on the measured TPA scale; the PNNL pool belonged to the PNNL abundance scale
        genes_in = {g for r_ in cost_rows for g in str(r_["Gene"]).replace("(", " ").replace(")", " ").replace(" and ", " ").replace(" or ", " ").split()}
        genes_in |= pool_iso
        pool_meas = float(sum(A[g] for g in genes_in if g in A and pd.notna(A[g])))
        print(f"  {cond}: protein pool {float(pool.iloc[0, 1]):.4f} (PNNL scale) -> {pool_meas:.4f} g/gDW (measured, {len(genes_in)} enzyme genes)", flush=True)
        pool = pool.copy()
        pool.iloc[0, 1] = pool_meas
        pool.to_excel(w, sheet_name="proteome_pool", index=False)
    ra_py = f"{work}/relaxation_algorithm.py"
    if os.path.exists(ra_py):   # the algorithm's own target pool is a constant in the script: set it to the measured pool as well
        import re as _re
        txt = open(ra_py).read()
        txt2 = _re.sub(r'("TARGET_POOL_g_per_gDW":\s*)[0-9.eE+-]+', lambda mm: f"{mm.group(1)}{pool_meas:.6f}", txt)
        assert txt2 != txt, "TARGET_POOL_g_per_gDW not found in relaxation_algorithm.py"
        open(ra_py, "w").write(txt2)
    print(f"  {cond}: {n_pred} new reactions carry predicted (DeepEnzyme) kcats; {n_new - n_pred} keep the median; "
          f"{len({q['gene'] for q in requests})} genes still want a kcat", flush=True)
    print(f"{cond}: {n_rescaled} carried-over bounds rescaled to the measured abundance, {n_median} by the median ratio {med_ratio:.3f}", flush=True)
    return work, n_new, len(cost_rows), len(ub_rows), requests


def recalibrate(cond, tag, work):
    cwd = os.getcwd()
    os.chdir(work)
    try:
        exec(compile(open("relaxation_algorithm.py").read(), "relaxation_algorithm.py", "exec"), {"__name__": "__main__"})
    finally:
        os.chdir(cwd)
    res = [f for f in os.listdir(work) if f.startswith("iterative_kcat_fix_results")][0]
    shutil.copy(f"{work}/UB_{tag}_revised.xlsx", f"{GEM}/UB_{tag}_patched_revised.xlsx")
    shutil.copy(f"{work}/kcat_mw_{tag}_revised.xlsx", f"{GEM}/kcat_mw_{tag}_patched_revised.xlsx")
    shutil.copy(f"{work}/{res}", f"{GEM}/relaxation_{tag}_patched_results.xlsx")
    shutil.copy(f"{work}/irrev_iPiccel954.xml", f"{GEM}/irrev_iPiccel954_patched_{cond}.xml")
    s = pd.read_excel(f"{work}/{res}", sheet_name="summary").iloc[0].to_dict()
    shutil.rmtree(work)
    return s


NGAM_LB = 0.62
MAINT_RELAXED = []


def ensure_maintenance(m, max_rounds=60):
    """Put the forced non-growth ATP maintenance back and relax, one binding bound at a time (largest reduced cost first), only what the
    measured-abundance bounds need to sustain it. The relaxed reactions are kept in MAINT_RELAXED and printed."""
    if "NGAM" not in m.reactions:
        return
    r = m.reactions.NGAM
    with m:
        r.lower_bound = 0.0
        m.objective = r
        for k in range(max_rounds):
            sol = m.optimize()
            if sol.status != "optimal":
                break
            if sol.objective_value >= NGAM_LB - 1e-9:
                break
            bind = sorted(((abs(sol.reduced_costs[x.id]), x.id) for x in m.reactions
                           if 0 < x.upper_bound < 999 and abs(sol.fluxes[x.id] - x.upper_bound) < 1e-7 and abs(sol.reduced_costs[x.id]) > 1e-9), reverse=True)
            if not bind:
                break
            top = m.reactions.get_by_id(bind[0][1])
            MAINT_RELAXED.append(top.id)
            top.upper_bound = 1000.0
        keep = {x.id: x.upper_bound for x in m.reactions if x.id in set(MAINT_RELAXED)}
    for rid, ub in keep.items():
        m.reactions.get_by_id(rid).upper_bound = ub
    r.lower_bound = NGAM_LB


def build(model_path, ub_path, kcat_path):
    m = cobra.io.read_sbml_model(model_path)
    rx = {r.id: r for r in m.reactions}
    ub = pd.read_excel(ub_path, usecols=[0, 1], names=["R", "U"])
    for _, row in ub.iterrows():
        r = rx.get(str(row.R))
        if r is not None and pd.notna(row.U) and r.bounds != (0, 0):
            r.upper_bound = max(0.0, float(row.U))
    ensure_maintenance(m)
    cost = pd.read_excel(kcat_path, sheet_name="reaction_protein_cost")
    pool = float(pd.read_excel(kcat_path, sheet_name="proteome_pool").iloc[0, 1])
    expr = sum((float(c.Protein_cost) / 1000.0) * rx[str(c.Reaction)].flux_expression
               for _, c in cost.iterrows() if str(c.Reaction) in rx and pd.notna(c.Protein_cost))
    m.add_cons_vars(m.problem.Constraint(expr, lb=0, ub=pool, name="protein_capacity"))
    m.objective = BIO
    return m


def ko(m):
    m.genes.get_by_id(PGM_GENE).knock_out()
    b = m.reactions.get_by_id(BIO)
    s = m.metabolites.get_by_id("starch300_h")
    b.add_metabolites({s: -b.metabolites[s]})
    m.reactions.get_by_id("STARCH300S").bounds = (0, 0)


if __name__ == "__main__":
    mw = fasta_mw(FASTA)
    ab = measured_abundance()
    ab_old = pd.read_excel(ABUND)
    tests, all_requests = [], []
    substrates_only = "--substrates-only" in sys.argv
    with tempfile.TemporaryDirectory() as tmp:
        zipfile.ZipFile(ZIP).extractall(tmp)
        for cond, tag in CONDS.items():
            work, n_new, n_cost, n_ub, requests = prepare(cond, tag, tmp, mw, ab, ab_old)
            all_requests += requests
            if substrates_only:
                shutil.rmtree(work)
                continue
            s = recalibrate(cond, tag, work)
            print(f"{cond}: {n_cost} cost rows ({n_new} new, median kcat), {n_ub} upper-bound rows | relaxation: step-1 "
                  f"{s['Step1_target_reached']}, relaxed {s['N_reactions_relaxed_step1']}, boosted {s['N_reactions_boosted_step2']}, "
                  f"growth {s['Final_growth_achieved']:.4f} (target {s['Target_growth']}), pool {s['Final_pool_achieved']:.4f}")
            m = build(f"{GEM}/irrev_iPiccel954_patched_{cond}.xml", f"{GEM}/UB_{tag}_patched_revised.xlsx",
                      f"{GEM}/kcat_mw_{tag}_patched_revised.xlsx")
            wt = m.slim_optimize(error_value=0.0)
            ue = m.metabolites.get_by_id("urea_e")
            urea = [r for r in m.reactions if r.id.startswith("EX_urea(e)") and r.metabolites.get(ue, 0) > 0]  # uptake copies
            with m:
                ko(m)
                k = m.slim_optimize(error_value=0.0)
            f = m.optimize().fluxes
            u0 = sum(f[r.id] for r in urea)
            tests.append(dict(condition=cond, urea_cap="replete", WT=wt, KO=k, KO_over_WT=k / wt if wt else None))
            for frac in (0.8, 0.6, 0.4, 0.2):
                with m:
                    for r in urea:
                        r.upper_bound = min(r.upper_bound, frac * abs(u0)) if abs(u0) > 0 else r.upper_bound
                    w2 = m.slim_optimize(error_value=0.0)
                    ko(m)
                    k2 = m.slim_optimize(error_value=0.0)
                tests.append(dict(condition=cond, urea_cap=frac, WT=w2, KO=k2, KO_over_WT=k2 / w2 if w2 else None))
            print(f"  WT {wt:.4f} | PGM KO {k:.4f} | KO/WT {k / wt if wt else float('nan'):.3f} | replete urea uptake {abs(u0):.3f}")
    q = pd.DataFrame(all_requests)
    q.to_csv(SUBSTRATE_REQUESTS, sep="\t", index=False)
    print(f"substrate requests: {len(q)} rows, {q.gene.nunique() if len(q) else 0} genes, "
          f"{q.reaction.nunique() if len(q) else 0} reactions -> {SUBSTRATE_REQUESTS}")
    if substrates_only:
        sys.exit(0)
    t = pd.DataFrame(tests)
    t.to_csv(f"{GEM}/ecGEM_patch_tests.tsv", sep="\t", index=False)
    print(t.to_string(index=False))
