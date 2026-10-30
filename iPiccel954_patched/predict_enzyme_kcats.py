"""DeepEnzyme turnover-number prediction for the P. celeri enzymes that the ecGEMs otherwise cost at the median kcat.

Why: patch_ecGEMs.py gives every reaction that gained a gene in the patch (steps E-G) a cost row with Niaz's median
TurNuP kcat (78,662 h-1 = 21.9 s-1) as a placeholder. This script predicts a per-enzyme kcat from sequence + ESM3
structure with the published DeepEnzyme model (Wang et al. 2024, Brief Bioinform), reusing the vendored inference
wrapper from PGM/ml_variant_design/scripts/deepenzyme_runner.py. The three DGAT genes were the first case
(dgat_kcat_predictions.tsv, 2026-10-17); this is the same procedure for every remaining gene.

Which genes and substrates: kcat_substrate_requests.tsv, written by `patch_ecGEMs.py --substrates-only`, lists every
(reaction, gene) still at the median with the reaction's non-currency substrates. Per gene, up to MAX_SUBSTRATES distinct
substrates are run (ranked by the number of reactions they appear in, one species per lipid class), and the gene's kcat
is the geometric mean over them (patch_ecGEMs.py).

Substrate SMILES, in this order (recorded per metabolite in kcat_substrate_smiles.tsv):
  1. KEGG compound ID -> KEGG REST molfile -> RDKit;
  2. InChIKey -> PubChem;
  3. ChEBI ID -> ChEBI web service;
  4. acyl species without a database ID, built from the iRC1080 id (acyl code CCD + positions, e.g. acoa1819Z = 18:1
     delta-9Z CoA) and the model formula: acyl-CoA, free fatty acid, acyl-ACP (the ACP arm represented by
     N-acetylcysteamine, the standard ACP mimic), 1-acyl-G3P; glycerolipid classes with two or three acyl chains
     (PA, PGP, CDP-DAG, DAG, PE, PG, PC, PI, MGDG, DGDG, SQDG, TAG, ceramides) by one representative species per class
     (16:0/18:1 or dipalmitoyl), a stated simplification: the model sees one species of the class;
  5. metabolite name -> PubChem.
Control: the authors' example (human creatine kinase with creatine) is rerun first.
Resumable: every prediction is appended to kcat_predictions.jsonl at once and skipped on restart.

Inputs: kcat_substrate_requests.tsv, genome_annotation/genome/Pceleri_fullproteome_nonredundant_v3.fasta,
        genome_annotation/structures/all_pdbs_flat/<gene>_<suffix>.pdb
Outputs (GEM/): kcat_substrate_smiles.tsv, kcat_predictions.jsonl, kcat_predictions.tsv
Run:    eval "$(/home/yasin/miniconda3/bin/conda shell.bash hook)" && conda activate deepenzyme && \
        python -I /home/yasin/Claudio/GEM/predict_enzyme_kcats.py [--limit N]
"""
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

sys.path.insert(0, "/home/yasin/Claudio/PGM/ml_variant_design/scripts")
from deepenzyme_runner import Predictor, REPO  # noqa: E402

RDLogger.DisableLog("rdApp.*")
GEM = Path("/home/yasin/Claudio/GEM")
FASTA = GEM / "genome_annotation/genome/Pceleri_fullproteome_nonredundant_v3.fasta"
PDB = GEM / "genome_annotation/structures/all_pdbs_flat"
REQ = GEM / "kcat_substrate_requests.tsv"
SMI = GEM / "kcat_substrate_smiles.tsv"
OUT_JSONL = GEM / "kcat_predictions.jsonl"
OUT_TSV = GEM / "kcat_predictions.tsv"
MAX_SUBSTRATES = 4
CTRL_SEQ = ("MPFGNTHNKFKLNYKPEEEYPDLSKHNNHMAKVLTLELYKKLRDKETPSGFTVDDVIQTGVDNPGHPFIMTVGCVAGDEESYEVFKELFDPIISDRHGGYKPTDKHKTDLNHENLKGGDDLDPNYVLSSRVRTGRSIKGYTLPPHCSRGERRAVEKLSVEALNSLTGEFKGKYYPLKSMTEKEQQQLIDDHFLFDKPVSPLLLASGMARDWPDARGIWHNDNKSFLVWVNEEDHLRVISMEKGGNMKEVFRRFCVGLQKIEEIFKKAGHPFMWNQHLGYVLTCPSNLGTGLRGGVHVKLAHLSTHPKFEEILTRLRLQKRGTGGVDTAAVGSVFDVSNADRLGSSEVEQVQLVVDGVKLMVEMEKKLEKGQSIDDMIPAQK")
CTRL_PDB = REPO / "Data/Input/2.7.3.2_Homo_sapiens_Seq_5476_relaxed_rank_1_model_3.pdb"

# ---------------------------------------------------------------- acyl-chain SMILES from the iRC1080 id
COA = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
SNAC = "SCCNC(C)=O"  # N-acetylcysteamine thioester arm, the ACP mimic
G3P = "OC[C@H](O)COP(=O)(O)O"
CHAIN_NAMES = {"but": 4, "buta": 4, "penta": 5, "hex": 6, "hexa": 6, "hx": 6, "hepta": 7, "oc": 8, "octa": 8, "nona": 9, "non": 9,
               "dca": 10, "dec": 10, "deca": 10, "undeca": 11, "undec": 11, "ddca": 12, "dodeca": 12, "trideca": 13,
               "ttdca": 14, "myrs": 14, "tetradeca": 14, "pentadeca": 15, "hdca": 16, "palm": 16, "hdec": 16, "hexadeca": 16,
               "hpdca": 17, "heptadeca": 17, "heptdec": 17, "ocdca": 18, "ocdc": 18, "st": 18, "octadeca": 18, "octadec": 18,
               "nondca": 19, "nonadeca": 19, "nonadec": 19, "arach": 20, "eicosa": 20, "heneca": 21, "henecosa": 21,
               "heneicosa": 21, "docosa": 22, "docosac": 22, "docos": 22, "tricosa": 23, "tricosac": 23, "ttcosa": 24,
               "lgnc": 24, "tetracosa": 24, "pentacosa": 25, "pentacosac": 25, "hexacosa": 26, "heptacosa": 27,
               "octacosa": 28, "nonacosa": 29, "triaconta": 30}


def acyl(n, positions, hydroxy2=False):
    """Acyl group SMILES written from the carbonyl: C(=O) then n-1 carbons with cis double bonds at the positions."""
    if n is None or n < 2:
        return None
    s = "C(=O)"
    if hydroxy2:
        s += "C(O)"
        start = 3
    else:
        start = 2
    i = start
    while i <= n:
        if i in positions and i + 1 <= n:
            s += "/C=C\\"
            i += 2
        else:
            s += "C"
            i += 1
    return s


def code_chain(code):
    """'1819Z' -> (18, [9]); '160' -> (16, []); '2239Z12Z15Z' -> (22, [9, 12, 15])."""
    mm = re.fullmatch(r"(\d{2})(\d)((?:\d+[ZE])*)", code)
    if not mm:
        return None, []
    n, d, pos = int(mm.group(1)), int(mm.group(2)), [int(p) for p in re.findall(r"(\d+)[ZE]", mm.group(3))]
    if d and not pos:  # positions not written: methylene-interrupted from delta-9
        pos = [9 + 3 * k for k in range(d)]
    return n, pos


def chain_from_formula(formula, c_offset, h_sat_offset):
    """Chain length and double bonds from the neutral formula: C = n + c_offset, H_sat = 2n + h_sat_offset."""
    mc = re.search(r"C(\d+)", formula)
    mh = re.search(r"H(\d+)", formula)
    if not (mc and mh):
        return None, []
    n = int(mc.group(1)) - c_offset
    d = (2 * n + h_sat_offset - int(mh.group(1))) // 2
    return n, [9 + 3 * k for k in range(max(d, 0))]


REP = {  # one representative species per glycerolipid / sphingolipid class (16:0 at sn-1, 18:1(9Z) at sn-2)
    "pa": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "pgp": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@@H](O)COP(=O)(O)O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "pg": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@@H](O)CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "cdp12dgr": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2ccc(N)nc2=O)[C@H](O)[C@@H]1O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "cdp": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2ccc(N)nc2=O)[C@H](O)[C@@H]1O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "12dgr": "CCCCCCCCCCCCCCCC(=O)OC[C@H](CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "pe": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "pchol": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "pail": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "mgdg": "CCCCCCCCCCCCCCCC(=O)OC[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "dgdg": "CCCCCCCCCCCCCCCC(=O)OC[C@H](CO[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "sqdg": "CCCCCCCCCCCCCCCC(=O)OC[C@H](CO[C@@H]1O[C@H](CS(=O)(=O)O)[C@H](O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "tag": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "clpn": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@@H](O)COP(=O)(O)OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "crm": "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)/C=C/CCCCCCCCCCCCC",
    "dhcrm": "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",
    "glucrm": "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)/C=C/CCCCCCCCCCCCC",
    "sphgn": "CCCCCCCCCCCCCCC[C@@H](O)[C@@H](N)CO",
}


def lipid_class(stem):
    """Class prefix of a species-expanded lipid id (12dgr1601819Z -> 12dgr, pa18111Z160 -> pa); '' when none."""
    s = re.sub(r"\(.*?\)", "", stem)
    for k in sorted(REP, key=len, reverse=True):
        if s.startswith(k) and re.match(r"\d", s[len(k):] or "0"):
            return k
    return re.sub(r"\d.*$", "", s)


def _join(a, tail, src):
    return (a + tail, src) if a else (None, None)


def build_smiles(stem, formula):
    """SMILES for an acyl species without a database ID, from the iRC1080 id and the model formula; None if not a class."""
    hydroxy = "(2OH)" in stem
    s = stem.replace("(2OH)", "")
    m1 = re.fullmatch(r"acoa(\d{2}\d(?:\d+[ZE])*)", s)
    if m1:
        n, pos = code_chain(m1.group(1))
        return _join(acyl(n, pos, hydroxy), COA, "built: acyl-CoA from id")
    m2 = re.fullmatch(r"ac(\d{2}\d(?:\d+[ZE])*)", s)
    if m2:
        n, pos = code_chain(m2.group(1))
        return _join(acyl(n, pos, hydroxy), "O", "built: fatty acid from id")
    m3 = re.fullmatch(r"([a-z]+?)(?:e|di|tri)?(?:\((\d+Z(?:\d+Z)*)\))?coa", s)
    if m3 and "N7O17P3S" in formula.replace("O18", "O17"):
        n, pos = chain_from_formula(formula, 21, 34 + (2 if hydroxy else 0))
        return _join(acyl(n, pos, hydroxy), COA, "built: acyl-CoA from formula")
    m4 = re.fullmatch(r"([a-z0-9]+?)(?:ACP)", s)
    if m4 and formula.endswith("SR"):
        mc, mh = re.search(r"C(\d+)", formula), re.search(r"H(\d+)", formula)
        n = int(mc.group(1))
        d = (2 * n - 1 - int(mh.group(1))) // 2
        pos = [9 + 3 * k for k in range(max(d, 0))]
        return _join(acyl(n, pos, hydroxy), SNAC, "built: acyl-ACP from formula, N-acetylcysteamine arm")
    m5 = re.fullmatch(r"1([a-z]+?)(?:e|di|tri)?(?:\((\d+Z(?:\d+Z)*)\))?g3p", s)
    if m5:
        n, pos = chain_from_formula(formula, 3, 5)  # model 1-acyl-G3P (charge -2) C(n+3)H(2n+5)O7P
        if m5.group(2):
            pos = [int(p) for p in re.findall(r"(\d+)Z", m5.group(2))]
        return _join(acyl(n, pos, hydroxy), G3P, "built: 1-acyl-G3P from formula")
    if re.fullmatch(r"C\d+H\d+O[23]", formula):
        for suf in ("", "ac", "a", "c"):
            base = s[:-len(suf)] if suf else s
            if base in CHAIN_NAMES:
                n = CHAIN_NAMES[base]
                mh = int(re.search(r"H(\d+)", formula).group(1))
                d = (2 * n - mh) // 2
                return _join(acyl(n, [9 + 3 * k for k in range(max(d, 0))], hydroxy), "O", "built: fatty acid from name")
    cls = lipid_class(stem)
    if cls in REP and REP[cls]:
        return REP[cls], f"built: representative {cls} species 16:0/18:1"
    return None, None


# ---------------------------------------------------------------- database lookups
def fetch(url, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as h:
                return h.read().decode()
        except urllib.error.HTTPError as e:
            err = e
            if e.code == 404:
                break  # not found is final, no retry
            time.sleep(2 + 3 * k)
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(2 + 3 * k)
    print("  fetch failed:", url, err)
    return None


def smiles_kegg(cid):
    t = fetch(f"https://rest.kegg.jp/get/{cid}/mol")
    if t:
        mol = Chem.MolFromMolBlock(t)
        if mol is not None:
            return Chem.MolToSmiles(mol)
    return None


def smiles_pubchem_inchikey(key):
    for prop in ("SMILES", "IsomericSMILES", "ConnectivitySMILES"):
        t = fetch(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{key}/property/{prop}/TXT")
        if t and not t.startswith("Status"):
            return t.strip().splitlines()[0]
    return None


def smiles_pubchem_name(name):
    q = urllib.parse.quote(name)
    for prop in ("SMILES", "IsomericSMILES"):
        t = fetch(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/{prop}/TXT")
        if t and not t.startswith("Status"):
            return t.strip().splitlines()[0]
    return None


def smiles_chebi(cid):
    t = fetch(f"https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId={cid}")
    if t:
        mm = re.search(r"<smiles>([^<]+)</smiles>", t)
        if mm:
            return mm.group(1)
    return None


def valid(smi):
    return smi is not None and Chem.MolFromSmiles(smi) is not None


def resolve(stem, name, formula, kegg, chebi, inchikey, cache):
    if stem in cache:
        return cache[stem]
    smi, src = None, None
    if kegg:
        smi = smiles_kegg(kegg)
        src = f"KEGG {kegg}"
    if not valid(smi) and inchikey:
        smi, src = smiles_pubchem_inchikey(inchikey), f"PubChem {inchikey}"
    if not valid(smi) and chebi:
        smi, src = smiles_chebi(chebi), f"ChEBI {chebi}"
    if not valid(smi):
        smi, src = build_smiles(stem, formula)
    if not valid(smi) and name and not re.search(r"\[[a-z]\]$", name):
        smi, src = smiles_pubchem_name(name), f"PubChem name '{name}'"
    if not valid(smi):
        return None, "unresolved"  # not cached, so a transient network failure is retried on the next run
    cache[stem] = (smi, src)
    return smi, src


def read_fasta(path):
    seqs, cur = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip()
    return seqs


def main():
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else None
    q = pd.read_csv(REQ, sep="\t", dtype=str).fillna("")
    q["stem"] = q.metabolite.str.replace(r"_[a-z]$", "", regex=True)
    q["cls"] = q.stem.map(lipid_class)
    q["has_id"] = (q.kegg != "") | (q.inchikey != "") | (q.chebi != "")
    seqs = read_fasta(FASTA)
    gene_key = {k.split("_")[0]: k for k in seqs}
    pdbs = {p.name.split("_")[0]: p for p in PDB.glob("GENE*.pdb")}
    cache = {}
    if SMI.exists() and SMI.stat().st_size > 10:
        for _, r in pd.read_csv(SMI, sep="\t", dtype=str).fillna("").iterrows():
            if r.smiles:
                cache[r.stem] = (r.smiles, r.source)
    done = set()
    if OUT_JSONL.exists():
        for line in OUT_JSONL.read_text().splitlines():
            if line.strip():
                d = json.loads(line)
                done.add((d["gene"], d["substrate"]))
    P = Predictor()
    if ("control", "creatine") not in done and CTRL_PDB.exists():
        r = P.predict(CTRL_SEQ, "CN(CC(=O)O)C(=N)N", CTRL_PDB)
        rec = dict(gene="control", substrate="creatine", substrate_name="control: human creatine kinase (authors' example)",
                   smiles_source="example", smiles="CN(CC(=O)O)C(=N)N", length=len(CTRL_SEQ), n_reactions=0,
                   kcat_per_s=r["kcat"], log2_kcat=r["log2_kcat"], oov_seq_ngrams=r["oov_seq_ngrams"], oov_fingerprints=r["oov_fingerprints"])
        with open(OUT_JSONL, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"control creatine kinase: {r['kcat']:.1f} s-1 (published example ~118)")
    genes = sorted(q.gene.unique())
    if "--genes" in sys.argv:
        want = set(sys.argv[sys.argv.index("--genes") + 1].split(","))
        genes = [g for g in genes if g in want]
    if limit:
        genes = genes[:limit]
    n_pred = n_skip = 0
    unresolved = {}
    t0 = time.time()
    for gi, g in enumerate(genes):
        sub = q[q.gene == g]
        ranked = (sub.groupby(["stem", "cls"]).agg(n=("reaction", "nunique"), has_id=("has_id", "max"), name=("name", "first"),
                                                     formula=("formula", "first"), kegg=("kegg", "first"), chebi=("chebi", "first"),
                                                     inchikey=("inchikey", "first")).reset_index()
                  .sort_values(["has_id", "n"], ascending=[False, False]))
        chosen, seen_cls = [], set()
        for _, s in ranked.iterrows():
            if len(chosen) >= MAX_SUBSTRATES:
                break
            if not s.has_id and s.cls in seen_cls:
                continue  # one built species per lipid class
            smi, src = resolve(s.stem, s["name"], s.formula, s.kegg, s.chebi, s.inchikey, cache)
            if smi is None:
                unresolved[s.stem] = s["name"]
                continue
            chosen.append((s, smi, src))
            if not s.has_id:
                seen_cls.add(s.cls)
        key = gene_key.get(g)
        pdb = pdbs.get(g)
        if key is None or pdb is None:
            print(f"{g}: no sequence/structure, skipped")
            continue
        for s, smi, src in chosen:
            if (g, s.stem) in done:
                n_skip += 1
                continue
            try:
                r = P.predict(seqs[key], smi, pdb)
            except Exception as e:  # noqa: BLE001
                print(f"{g} {s.stem}: prediction failed: {e}")
                continue
            rec = dict(gene=key, substrate=s.stem, substrate_name=s["name"], smiles_source=src, smiles=smi, length=len(seqs[key]),
                       n_reactions=int(s.n), kcat_per_s=r["kcat"], log2_kcat=r["log2_kcat"], oov_seq_ngrams=r["oov_seq_ngrams"],
                       oov_fingerprints=r["oov_fingerprints"])
            with open(OUT_JSONL, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            done.add((g, s.stem))
            n_pred += 1
        if gi % 10 == 0:
            print(f"[{gi + 1}/{len(genes)}] {g}: {len(chosen)} substrates | predicted {n_pred}, skipped {n_skip} | {time.time() - t0:.0f} s", flush=True)
        pd.DataFrame([dict(stem=k, smiles=v[0] or "", source=v[1]) for k, v in cache.items()]).to_csv(SMI, sep="\t", index=False)
    rows = [json.loads(l) for l in OUT_JSONL.read_text().splitlines() if l.strip()]
    d = pd.DataFrame(rows)
    d["kcat_per_h"] = d.kcat_per_s * 3600
    d.to_csv(OUT_TSV, sep="\t", index=False)
    g = d[d.gene.str.startswith("GENE")]
    print(f"\n{len(g)} predictions for {g.gene.nunique()} genes; kcat s-1 median {g.kcat_per_s.median():.2f}, "
          f"IQR {g.kcat_per_s.quantile(.25):.2f}-{g.kcat_per_s.quantile(.75):.2f}, range {g.kcat_per_s.min():.3f}-{g.kcat_per_s.max():.1f}")
    print("substrate sources:", d.smiles_source.str.split(" ").str[0].value_counts().to_dict())
    print(f"unresolved substrates ({len(unresolved)}):", "; ".join(f"{k} ({v})" for k, v in list(unresolved.items())[:40]))
    print("wrote", OUT_TSV)


if __name__ == "__main__":
    main()
