"""Patch a COPY of the PNNL P. celeri GEM (iPiccel954) so it is mass balanced, free of energy-generating cycles,
and carries the reactions and genes that the P. celeri genome, the two source models (iPicre1230, RAPS template
iCre1355) and the reference algal GEMs support.

Input  (read only, never modified): /home/yasin/Claudio/GEM/p_celeri_GEM/iPiccel954.xml
Output: /home/yasin/Claudio/GEM/iPiccel954_patched.xml
Evidence tables: iPicre1230_vs_iPiccel954_reaction_map.xlsx, iCre1355_template_vs_iPiccel954.xlsx,
iPicre1230_GPR_transfer.xlsx (GEM/), and the P. celeri annotation (curated layer used where present).

A. Carbon partitioning (kept from the first round)
  A1 triose-phosphate/phosphate translocator GENE4065742: 3PG_pi_thr reversible + DHAP and GAP antiport
     (iPicre1230 DHAP(pi)thr / G3P(pi)thr, all three reference GEMs).
  A2 cytosolic aldolase fbaCv reversible on GENE4062223 / GENE4068431.
  A3 cytosolic FBPase FBPc on GENE4066543 (F16P2; iPicre1230 BFBP).
  A4 lumped fatty-acid elongation FACOAL*_ACP* uses 2 NADPH instead of 4 H (KAR and ENR steps).
  A5 CoA exchange COAth reversible (iRC1080, iCre1355).
  A6 AGAKI (glucose-fructose isomerase on a ribokinase gene) closed.
  Plastid G1P import is LEFT OPEN (it was closed in round 1): iRC1080 and iCre1355 keep it and Fettke 2011
  measured G1P uptake. The PGM knockout instead fixes starch synthesis at zero (pgm_phenotype.py).
B. Energy-generating cycles (Fritzemeier 2017; checked with mass-balanced drains, energy_cycle_check below)
  B1 ppGpp reactions GBDP, GDTP, GTPDPK forward only (the source forced them backwards).
  B2 mitochondrial ATP export ATP_3h_tm carries 1 H+ (adp_c + atp_m + h_c -> adp_m + atp_c + h_m), the usual
     cost of AAC + phosphate carrier; the source imported 3 H+ per ATP while the synthase used 2.
  B3 mitochondrial ATP synthase NanoG0606 3 H+ per ATP.
  B4 CY_focytc_h (cytochrome c6 reduced by a bare proton) closed.
  B5 proton export allowed (EX_h upper bound 1000, as iCre1355).
C. Reactions the genome supports but the model lacked (gene from our annotation; precedent in iPicre1230)
  C1 cytosolic and mitochondrial pyrophosphatase PPAc GENE4065198, PPAm GENE4068945 (iPicre1230 PPA).
  C2 cytosolic NDPK for GDP/UDP/CDP, GENE4065673 (iPicre1230 ATGD).
  C3 mitochondrial complex I (mitochondrial nad genes + nuclear subunits) and complex II to ubiquinone
     (GENE4064131, GENE4066270) (iPicre1230 NADHOR, SUCDH(q8)m).
  C4 ADP-dependent pyruvate kinase PYK_c on the three pyruvate kinase genes (iPicre1230 and iCre1355 PYK; the
     source's pyruvate kinase reactions use CDP or GDP only).
  C5 plastid NADP-glyceraldehyde-3-phosphate dehydrogenase on GapA GENE4065210 (EC 1.2.1.13; iCre1355
     GAPDH_nadp_hi): the source's only plastid GAPDH uses NADH.
  C6 plastid NADP-malate dehydrogenase on GENE4062923 (EC 1.1.1.82; iCre1355 MDH_nadp_hi), the malate valve.
  C7 water and CO2 diffusion cytosol-chloroplast and cytosol-mitochondrion (iPicre1230 H2Oth/H2Otm/CO2th/CO2tm,
     all reference GEMs); without them the stroma had to make its own water and import every carbon with a proton.
  C8 Calvin-cycle hexose branch: plastid aldolase FBAh reversible (iRC1080) and plastid FBPase FBPh on GENE4063806
     (FBPase 1, chloroplastic). Without them plastid FBP was a dead end and starch carbon had to come from the cytosol.
D. Mass balance and electron units
  D1 FNORh balanced (fdxrd_u + nadp_h + 2 h_h -> fdxox_u + nadph_h, the iPicre1230 form); CBFC releases the 2 H
     its plastocyanin carries. Both were made unbalanced in Niaz's curation stage 2 (his Table S5G).
  D2 FDXuh links the 2-electron thylakoid ferredoxin to the 1-electron stromal pool 1:2 (the 1:1 join lost an electron).
  D3 starch formula corrected to C1800H3002O1501 (300 glucose units); STARCH300S no longer releases 300 H+, which had
     made starch synthesis require nitrogen through proton-coupled amino-acid export.
E. Gene rules
  E1 duplicates: where two reactions have the same equation (protons ignored) and only one carries P. celeri genes,
     the gene-less copy is closed when its allowed directions are a subset of the gene copy's, so a knockout acts.
  E2 gene-less reactions receive P. celeri genes when the P. renovo proteome transfer (EC-checked) names the gene,
     the gene's effective EC (curated layer first, retractions honoured) matches the reaction's EC, and its
     predicted localization matches the reaction's compartment.
F. Storage oil and electron sinks (genes found by Pfam domain, 2026-10-17)
  F1 the 143 TAG-synthase reactions (TG*) receive the three DGAT genes GENE4063758, GENE4067271, GENE4071740
     (EC 2.3.1.20, tier 1). Before this no TAG-assembly reaction had a gene, so TAG was a free sink in the ecGEM
     (iRC1080, iCre1355 and iPicre1230 all carry genes on every TAG-assembly reaction).
  F2 flavodiiron O2 photoreduction FLVh on GENE4067072 and GENE4067074 (plastid, PF01613).
  F3 plastid terminal oxidase PTOXu on GENE4067969 (plastid, PF01786).
     F2/F3 give the model costed ways to dispose of excess electrons, so that the free photon demand reactions
     (DM_photon*, an iRC1080 convention kept by every algal GEM compared) can be closed in a simulation.
G. Gene rules recovered for the gene-less reactions (2026-10-18; 1,730 open internal reactions had no gene, 1,728 of
   them exist in iPicre1230 and 1,311 carry a P. renovo rule there, lost in the RAPS build)
  G1 the iPicre1230 rule translated through the P. renovo gene key (Pren_gene_key.tsv, reciprocal best hits): genes
     kept when their effective EC matches the reaction's EC, or, for reactions without EC, when the source rule
     names at most 6 genes (the generic transporter rules with dozens of genes are not translated).
  G2 a sibling reaction of the same family (same id prefix, same iPicre1230 rule, same EC) already carries a P. celeri
     rule: that rule is copied (the source model used one enzyme set for every acyl species of the family).
  G3 reactions with EC and no route above: P. celeri genes whose high-confidence level-4 EC matches (at most 8 genes;
     not for transport reactions; generic ATPase ECs 3.6.3.*/7.* excluded).
  G4 lipid species-expanded reactions without EC (iRC1080 convention): the EC is read from the reaction family
     (1*G3P_ACYLTRANSFERASE 2.3.1.15, AGPAT 2.3.1.51, PGPS 2.7.8.5, PGPP 3.1.3.27, CDPDAGS 2.7.7.41, CDIPT 2.7.8.11,
     MGDG 2.4.1.46, DGDG 2.4.1.241, ETHAPT 2.7.8.1, GLUCRM 2.4.1.80, DHCER 2.3.1.24, CER 1.14.19.17, DAG*_ACHYD
     3.1.1.23, CLPN*HYDROLASE 3.1.1.4) and matched as in G3. Desaturases (DESAT*, all 1.14.19.1 in the source) are
     assigned by the position of the new double bond: delta-9 GENE4063770/GENE4067065 (acyl-ACP 9-desaturases, the
     only delta-9 enzymes found), delta-12 GENE4062097/GENE4066813, delta-7 GENE4066008, omega-3 GENE4067276;
     positions without a P. celeri enzyme (delta-3, -4, -5, -6, -18 and longer than C22) stay gene-less.
  G5 (2026-10-18, second pass over the 876 reactions still gene-less; 342 of them have a source gene the key never
     placed, 236 carry an EC absent from the P. celeri annotation at any confidence, 225 carry no EC, so the key cannot
     be extended for them and the family is the only route) more families read to an EC: FACOAL*_ACP (plastid
     acyl-ACP elongation lumps, condensation 2.3.1.179, localization enforced: KAS I/II GENE4067847/GENE4064434),
     FACOAL*c/h (acyl-CoA ligase 6.2.1.3), FA*(2OH)ACPHi (acyl-ACP thioesterase 3.1.2.14), TGTAG/SQDGATPA (acyl-CoA
     acyltransferase 2.3.1.20, the DGAT genes of F1, same source rule as the TAG synthases), PETH*MT (phosphoethanolamine
     N-methyltransferase 2.1.1.103), mmcoaeCv 5.1.99.1, pcoacCv 6.4.1.3, DEDOLDPS 2.5.1.87, DEDOLDPSAT 1.3.1.94,
     BUP 3.5.1.6. Explicit rules: PCHLDR gets ChlL/ChlN/ChlB (GENE2038481 and GENE2038482 and GENE2038489,
     chloroplast-encoded, EC 1.3.7.7, a complex); PCHOL*_SYN (2.7.8.2) gets GENE4064348, annotated 2.7.8.1 but a
     member of InterPro IPR014472 "Choline/ethanolamine phosphotransferase" (the dual-specificity CEPT family), so
     the one P. celeri aminoalcohol phosphotransferase is given both activities.
     Still without a P. celeri enzyme after G5 (searched: EC field at any confidence, Pfam PF01066 / PF04116, gene
     names): cardiolipin synthase 2.7.8.41 (CLPN*_SYN; only Chopt, PGS1 and PIS carry the CDP-alcohol
     phosphatidyltransferase domain), beta-carotene hydroxylase (BCAROH/BCRPTXANH; no CrtR-b/BCH domain gene, only
     CYP97B3 GENE4066520 of unassigned function), 144 desaturations at positions without an enzyme, and the
     reactions whose EC is absent from the annotation (listed with their family in the report).
  Localization: for lipid subsystems the compartment is not checked (iRC1080 draws lipid reactions in the plastid
  while the enzymes are ER or cytosolic, the F1 precedent); elsewhere the gene's predicted localization must match a
  compartment of the reaction. Every rule is logged with its route in patch_iPiccel954_changes.tsv; the reactions still
  gene-less are listed with the reason in geneless_remaining.tsv.
H. Corpus-driven changes (2026-10-20; sources in GEM/corpus_*.tsv and Pceleri_corpus audit):
  H1 close the reactions of pathways reported ABSENT in Picochlorum: the glyoxylate cycle (isocitrate lyase ICL, malate
     synthase MALSm/MALSx), pyruvate formate lyase (PFLh, PFLACTm) and arginase (argCv). Evidence: Bec et al. 2025
     (Environ Microbiol; 31 trebouxiophyte genomes, ICL/MS/PFL absent in Picochlorum), Xiong et al. 2025 (arginase nearly
     universally missing), and our own search of the P. celeri annotation by Pfam (PF00463, PF01274, PF02901, PF01228,
     PF00491), EC (4.1.3.1, 2.3.3.9, 2.3.1.54, 3.5.3.1) and name, which found no such gene ("not found by that method",
     not proof of loss). The two genes ICL carried (GENE4063068, GENE4068608) are PEP phosphomutase-superfamily proteins
     (PF13714, annotated 2.7.8.23 and 4.1.3.32), not isocitrate lyases; malate synthase had no gene. An earlier attempt
     to close the glyoxylate cycle alone cut the light-limited biomass yield from 1.08 to 0.45 g per mol photons, because
     the model could not move mitochondrial ammonium out (step I); with step I the cycle closes at a 4 % yield cost.
  I. Nitrogen-assimilation module (2026-10-20), found by tracing why amino acids cost 23 photons per carbon against 9 for
     starch once malate synthase was closed. The patched model lacked every ammonium route between compartments
     (cytosol-plastid, cytosol-mitochondrion), plastid glutamine synthetase and ferredoxin-GOGAT, so mitochondrial NH4+
     from glycine decarboxylase and glutamate dehydrogenase could not leave, and glycine, pyruvate and ethanol had to be
     exported. iRC1080, iLB1027 and iPicre1230 carry NH4th (cytosol-plastid) and NH4tm (cytosol-mitochondrion) as
     diffusion reactions; iPicre1230, iRC1080 and iCre1355 carry ferredoxin-GOGAT and iPicre1230 plastid GS.
     I1 NH4th and NH4tm, reversible NH3/NH4+ diffusion, no gene (as in the sources).
     I2 plastid glutamine synthetase GALh on GENE4068506, the GS2 ortholog (reciprocal best hits to Q42689, 69.3 %
        identity, annotated Cytoplasm but plastidic by orthology, see project note of 2026-08-25).
     I3 ferredoxin-dependent glutamate synthase GLUS(ferr) on GENE4068118 (EC 1.4.7.1, plastid, DeepLoc 0.986 with
        transit peptide); the NADH-dependent GLUS(nadh) in the model was gene-less and now carries GENE4062926
        (glutamate synthase 1 [NADH], chloroplastic, EC 1.4.1.14).
     I4 the 2-oxoglutarate/malate translocator AKGMALth of iPicre1230 (the route that imports 2-oxoglutarate for plastid
        GS/GOGAT); the patched model had only the glutamate-malate exchange MALGLUth. It takes the MALGLUth gene rule
        (the iPicre1230 rules of the two reactions are identical, one dicarboxylate-translocator family). Without it
        glutamate cost 17 photons per carbon against 9 for 2-oxoglutarate; with it amino acids cost 9 to 13.
     I1 is written proton-free on purpose: the iRC1080/iPicre1230 form (H+ and NH4+ moved together) let the model run a
     urea -> NH4+ export cycle on a free proton shuttle (yield 2.15 g/mol, negative Ci uptake, 30-fold urea need), which
     the drain test of step energy_cycle_check does not catch; a Ci-uptake and NH4-export check now follows the patch.
  H2 non-growth ATP maintenance NGAM (ATP + H2O -> ADP + Pi, the model's proton convention, lower bound 0.62 mmol gDW-1 h-1). The patched model
     had none (growth-associated ATP only, 39.2 mmol gDW-1). Derived from the one Picochlorum maintenance measurement,
     Barten et al. 2022 Biotechnol Bioeng (BPE23, light-limited chemostats): maintenance rate 0.079 d-1 over a yield on
     light of 1.38 g mol-1 gives 0.0572 mol photons gDW-1 d-1 = 2.385 mmol photons gDW-1 h-1; the model needs 3.85
     photons per ATP at the optimum, so NGAM = 2.385 / 3.85 = 0.62 mmol ATP gDW-1 h-1. The Chlamydomonas value
     (2.85, Kliphuis 2012, Imam 2015) is the upper sensitivity. The energy-cycle test sets NGAM to zero.
  H3 condition-specific biomass reactions on the measured composition of P. celeri (LaPanse et al. 2024 Fig 2, same
     strain and reactor regime; cross-checked against Huesemann 2022 Table 2 and Steichen 2024 Table 2): the native
     biomass (protein 63%, carbohydrate 11%, lipid 13.5%, pigment 5.7%, nucleic 5.2% of core mass) is rescaled class by
     class, starch excluded as the free storage product. Biomass_replete_LaPanse: protein 51.3, carbohydrate 12.2,
     lipid 11.0, pigment 5.3, nucleic 5.7 (% AFDW, renormalised). Biomass_N40_LaPanse: 19.9, 10.8, 13.2, 3.1, 2.2 (pigment = chlorophyll 4.4 / 2.1 x 1.203, see CORE).
     The default objective Biomass_Picre_Auto is unchanged.
  J. Biomass micro-composition on the lab's own and published measurements (2026-10-20). The native biomass carried an
     amino-acid profile with arginine at 0.09 mol % (measured 5.8 CO2 / 4.0 air), asparagine+aspartate 15.6 (measured 9.3),
     proline 5.1 (8.5); chlorophyll a/b 3.0 (4.3 to 4.6 measured) and 5 % carotenoids in the pigments (about 20 % of
     chlorophyll mass measured); and a fatty-acid profile in which C12 to C14 and C19 to C40 chains made 40 % of the acyl
     groups (C18:3 10 %, against 30 % measured).
  J1 amino-acid profile: the 18 hydrolysis groups of the lab's amino-acid analysis (aaa_composition_mol_pct.tsv, high-CO2
     mean of four bottles; Asx split asparagine:aspartate and Glx split glutamine:glutamate in the native ratio) replace the
     native profile; cysteine and tryptophan keep their native share because acid hydrolysis destroys or oxidises them;
     the protein mass is unchanged.
  J2 pigments: chlorophyll a/b 4.4 (Cano 2024 4.3 to 4.6; Karns 2024 HL 5.5), carotenoids 0.25 of chlorophyll mass
     (Krishnan 2023 0.25 to 0.32, Cano 2024 about 0.20), lutein 10.3 % and beta-carotene 5 % of chlorophyll mass (Cano 2024
     Fig 1D), neoxanthin 3.4 % (Krishnan 2023: 0.12 % AFDW); the rest of the carotenoids keep their native proportions;
     the pigment mass is unchanged.
  J3 lipid fatty acids and classes: species with a chain outside C16/C18 are removed (the lab's FAME profile and LaPanse 2024 find
     no C12 to C14 or C19+ acyl groups of note; 219 of 299 species), and the species of the classes already in the biomass,
     including galactolipid species the native biomass did not use, are re-weighted by iterative proportional fitting (minimum
     change from the native weights) to the measured acyl profile: lab GC replete wild type mean of 10 am and 2 pm, C18:3 30.0,
     C16:0 19.2, C18:2 10.7, C18:1 1.3, C18:0 0.6 % of integrated peak area, LaPanse 2024 day 0 C16:0 about 23 and C18:3 about
     32; the ~33 % of the GC profile that is unassigned (probably C16:1 to C16:4; LaPanse lists C16:3 among the six main FAMEs)
     is one group not constrained within itself. Class shares: MGDG 22.5, PC 8.1 and TAG 5.1 % of lipid mass are pulled to
     Picochlorum sp. D3 (Grubisic 2022 Table 3, another strain and medium) and free DAG is capped at 3 % (an assumption), the other
     classes fill the rest pro rata, which inflates PI (33 %) and PG (21 %): no class composition exists for P. celeri (open
     item; DGDG and SQDG are absent from the native biomass and were not added). The C16 polyunsaturate split is open as well
     (the model carries mostly C16:1 because few galactolipid species with C16:3 survive the class shares; Delta-4 for C16:4 has no enzyme).
  J4 the default biomass Biomass_Picre_Auto becomes the replete measured composition (LaPanse 2024 classes, starch term kept);
     the pre-curation reaction is kept closed as Biomass_Niaz_native.
Open, not changed: PGCM (cytosolic PGM) keeps its "NF or GENE4067541" rule. P. celeri has two genes with the
three phosphoglucomutase domains (PF02878, PF02879, PF02880): the eukaryotic-type PGMP1 GENE4067541 (alleles p0/p1
GENE4070978, both predicted plastid by DeepLoc, 0.99/0.85) and a bacterial-type phosphohexomutase GENE4064687
(allele GENE4071452; Swiss-Prot hits split between PGM, phosphomannomutase and phosphoglucosamine mutase at 25-30%;
DeepLoc plastid/mitochondrion; in this model only on phosphomannomutase MPPM). Its substrate is not determinable from
sequence, and the starchless PGMP1 knockout shows it does not replace PGMP1 in the plastid, so no PGM rule is changed.

Run:  uv run --project /home/yasin/tools/ChatGEM python -I /home/yasin/Claudio/GEM/patch_iPiccel954.py
"""
import collections
import re

import cobra
import pandas as pd

GEM = "/home/yasin/Claudio/GEM"
SRC = f"{GEM}/p_celeri_GEM/iPiccel954.xml"
OUT = f"{GEM}/iPiccel954_patched.xml"
ANNOT = f"{GEM}/genome_annotation/annotation/Picochlorum_celeri_annotation.tsv"
TRANSFER = f"{GEM}/iPicre1230_GPR_transfer.xlsx"
REC = f"{GEM}/iPicre1230_reaction_EC.tsv"
LOG = f"{GEM}/patch_iPiccel954_changes.tsv"
BIO = "Biomass_Picre_Auto"
PGM_GENE = "GENE4067541"  # PGMP1, phosphoglucomutase 1, chloroplastic
COMPLEX_I = ("GENE2038411 and GENE2038413 and GENE2038414 and GENE2038417 and GENE2038418 and GENE2038424 and "
             "GENE2038425 and GENE2038438 and GENE4065393 and GENE4063537")
changes = []


def note(step, rid, what):
    changes.append(dict(step=step, reaction=rid, change=what))


def add(m, rid, stoich, gpr="", lb=0.0, ub=1000.0, name="", step=""):
    M = m.metabolites
    r = cobra.Reaction(rid, name=name, lower_bound=lb, upper_bound=ub)
    m.add_reactions([r])
    r.add_metabolites({M.get_by_id(k): v for k, v in stoich.items()})
    r.gene_reaction_rule = gpr
    assert not r.check_mass_balance(), (rid, r.check_mass_balance())
    note(step, rid, f"added: {r.reaction} [{gpr or 'no gene'}]")
    return r


def patch_A(m):
    M = m.metabolites
    t = m.reactions.get_by_id("3PG_pi_thr")
    t.bounds, t.gene_reaction_rule = (-1000, 1000), "GENE4065742"
    note("A1", t.id, "reversible, GENE4065742")
    for x in ["dhap", "g3p"]:
        add(m, f"TPT_{x}", {f"{x}_h": -1, "pi_c": -1, f"{x}_c": 1, "pi_h": 1}, "GENE4065742", -1000,
            name=f"triose phosphate/phosphate translocator ({x})", step="A1")
    f = m.reactions.get_by_id("fbaCv")
    f.bounds, f.gene_reaction_rule = (-1000, 1000), "GENE4062223 or GENE4068431"
    note("A2", f.id, "reversible, cytosolic aldolase genes")
    add(m, "FBPc", {"fdp-B_c": -1, "h2o_c": -1, "f6p-B_c": 1, "pi_c": 1}, "GENE4066543",
        name="fructose-1,6-bisphosphatase, cytosolic", step="A3")
    n = 0
    for r in m.reactions:
        if not (r.id.startswith("FACOAL") and "_ACP" in r.id):
            continue
        acyl = [x for x in r.metabolites if "ACP" in x.id]
        hs = [x for x in r.metabolites if x.id in ("h_c", "h_h", "h_m") and r.metabolites[x] == -4.0]
        if not hs:
            continue
        before = r.check_mass_balance()
        comp = acyl[0].compartment
        r.add_metabolites({hs[0]: 4.0, M.get_by_id(f"nadph_{comp}"): -2.0, M.get_by_id(f"nadp_{comp}"): 2.0})
        assert r.check_mass_balance() == before
        n += 1
    note("A4", "FACOAL*_ACP*", f"{n} elongation lumps use 2 NADPH instead of 4 H")
    m.reactions.COAth.bounds = (-1000, 1000)
    note("A5", "COAth", "reversible")
    m.reactions.AGAKI.bounds = (0, 0)
    note("A6", "AGAKI", "closed")


def patch_B(m):
    M = m.metabolites
    for rid in ["GBDP", "GDTP", "GTPDPK"]:
        m.reactions.get_by_id(rid).bounds = (0, 1000)
        note("B1", rid, "forward only")
    a = m.reactions.ATP_3h_tm
    # source: adp_m + atp_c + 3 h_m --> adp_c + atp_m + 3 h_c (ATP import); keep one proton, allow export only
    assert a.metabolites[M.h_m] == -3 and a.metabolites[M.h_c] == 3 and a.metabolites[M.atp_c] == -1, a.reaction
    a.add_metabolites({M.h_m: 2, M.h_c: -2})
    a.bounds = (-1000, 0)
    note("B2", a.id, f"{a.reaction} {a.bounds}")
    s = m.reactions.NanoG0606
    s.add_metabolites({M.h_c: -1, M.h_m: 1})
    assert not s.check_mass_balance()
    note("B3", s.id, s.reaction)
    m.reactions.CY_focytc_h.bounds = (0, 0)
    note("B4", "CY_focytc_h", "closed")
    m.reactions.get_by_id("EX_h(e)").upper_bound = 1000
    note("B5", "EX_h(e)", "export allowed")


def patch_C(m):
    add(m, "PPAc", {"ppi_c": -1, "h2o_c": -1, "pi_c": 2}, "GENE4065198", name="inorganic pyrophosphatase, cytosolic", step="C1")
    add(m, "PPAm", {"ppi_m": -1, "h2o_m": -1, "pi_m": 2}, "GENE4068945", name="inorganic pyrophosphatase, mitochondrial", step="C1")
    for n in ("g", "u", "c"):
        add(m, f"NDPK_{n}_c", {"atp_c": -1, f"{n}dp_c": -1, "adp_c": 1, f"{n}tp_c": 1}, "GENE4065673", -1000,
            name=f"nucleoside-diphosphate kinase ({n.upper()}DP), cytosolic", step="C2")
    add(m, "NADHOR_m", {"nadh_m": -1, "q8_m": -1, "h_m": -4, "nad_m": 1, "q8h2_m": 1, "h_c": 4}, COMPLEX_I,
        name="NADH:ubiquinone oxidoreductase (complex I)", step="C3")
    add(m, "SUCDH_q8_m", {"succ_m": -1, "q8_m": -1, "fum_m": 1, "q8h2_m": 1}, "GENE4064131 and GENE4066270", -1000,
        name="succinate dehydrogenase (complex II)", step="C3")
    add(m, "PYK_c", {"pep_c": -1, "adp_c": -1, "pyr_c": 1, "atp_c": 1}, "GENE4062482 or GENE4064386 or GENE4065717",
        name="pyruvate kinase (ADP), cytosol", step="C4")
    add(m, "GAPDH_nadp_h", {"13dpg_h": -1, "nadph_h": -1, "g3p_h": 1, "nadp_h": 1, "pi_h": 1}, "GENE4065210", -1000,
        name="glyceraldehyde-3-phosphate dehydrogenase (NADP), chloroplast", step="C5")
    add(m, "MDH_nadp_h", {"oaa_h": -1, "nadph_h": -1, "mal__L_h": 1, "nadp_h": 1}, "GENE4062923", -1000,
        name="malate dehydrogenase (NADP), chloroplast", step="C6")
    # C8 Calvin-cycle hexose branch. In the source (and iPicre1230) plastid FBP is a dead end: FBAh only cleaves it
    # and no plastid FBPase exists, so the stroma could not make hexose phosphate from trioses and starch carbon had
    # to come from the cytosol. FBAh is made reversible (iRC1080 FBAh is reversible; its rule already carries the
    # chloroplastic aldolase GENE4067854) and plastid FBPase is added on GENE4063806 (FBPase 1, chloroplastic).
    fa = m.reactions.FBAh
    fa.bounds = (-1000, 1000)
    note("C8", "FBAh", "reversible (aldolase condensation in the Calvin cycle)")
    add(m, "FBPh", {"fdp-B_h": -1, "h2o_h": -1, "f6p-B_h": 1, "pi_h": 1}, "GENE4063806",
        name="fructose-1,6-bisphosphatase, chloroplast", step="C8")
    for rid, a_, b_ in [("H2Oth", "h2o_c", "h2o_h"), ("H2Otm", "h2o_c", "h2o_m"),
                        ("CO2th", "co2_c", "co2_h"), ("CO2tm", "co2_c", "co2_m")]:
        add(m, rid, {a_: -1, b_: 1}, "", -1000, name=f"diffusion {a_} <-> {b_}", step="C7")


def patch_D(m):
    M = m.metabolites
    f = m.reactions.FNORh
    f.add_metabolites({M.h_h: -3})
    c = m.reactions.CBFC
    c.add_metabolites({M.h_u: 2})
    assert not f.check_mass_balance() and not c.check_mass_balance()
    note("D1", "FNORh", f.reaction)
    note("D1", "CBFC", c.reaction)
    add(m, "FDXuh", {"fdxrd_u": -1, "fdxox_h": -2, "fdxox_u": 1, "fdxrd_h": 2}, "", -1000,
        name="ferredoxin, thylakoid 2-electron pool to stromal 1-electron pool", step="D2")
    # D3 starch formula. starch300_h was C1800H2702O1501; a 300-glucose polymer is (C6H10O5)300 + H2O = C1800H3002O1501.
    # The 300 missing H were balanced by STARCH300S releasing 300 H+ into the stroma, and with no free stromal proton
    # exit those protons could only leave on proton-coupled amino-acid export, which made starch synthesis require
    # nitrogen (no starch at all without urea). Formula corrected, the 300 H removed from STARCH300S, and the starch
    # degradation reactions re-balanced on their protons.
    st = M.starch300_h
    st.formula = "C1800H3002O1501"
    s = m.reactions.STARCH300S
    s.add_metabolites({M.h_h: -s.metabolites[M.h_h]})
    assert not s.check_mass_balance(), s.check_mass_balance()
    note("D3", "starch300_h", "formula C1800H2702O1501 -> C1800H3002O1501")
    note("D3", "STARCH300S", s.reaction)
    for r in st.reactions:
        if r.id.startswith("STARCH300DEGR"):
            dh = r.check_mass_balance().get("H", 0)
            if dh:
                r.add_metabolites({M.h_h: -dh})
            assert not r.check_mass_balance(), (r.id, r.check_mass_balance())
            note("D3", r.id, r.reaction)


def eq_key(r):
    d = {}
    for x, v in r.metabolites.items():
        if x.id.rsplit("_", 1)[0] == "h":
            continue
        d[x.id] = round(v, 6)
    items = sorted(d.items())
    if not items:
        return None, 1
    sign = 1 if items[0][1] > 0 else -1
    return tuple((k, v * sign) for k, v in items), sign


def patch_E(m):
    groups = collections.defaultdict(list)
    for r in m.reactions:
        if r.boundary:
            continue
        k, s = eq_key(r)
        if k:
            groups[k].append((r, s))
    closed = 0
    for k, mem in groups.items():
        if len(mem) < 2:
            continue
        gene = [(r, s) for r, s in mem if "GENE" in r.gene_reaction_rule]
        bare = [(r, s) for r, s in mem if "GENE" not in r.gene_reaction_rule]
        if not gene or not bare:
            continue
        lo = min(min(r.lower_bound * s, r.upper_bound * s) for r, s in gene)
        hi = max(max(r.lower_bound * s, r.upper_bound * s) for r, s in gene)
        for r, s in bare:
            a, b = sorted((r.lower_bound * s, r.upper_bound * s))
            if lo <= a and b <= hi and r.bounds != (0, 0):
                r.bounds = (0, 0)
                closed += 1
                note("E1", r.id, f"closed: duplicates gene-bearing {', '.join(x.id for x, _ in gene)}")
    ann = pd.read_csv(ANNOT, sep="\t", dtype=str).fillna("")
    ann["gene"] = ann["gene_id"].str.split("_").str[0]
    ann = ann.drop_duplicates("gene").set_index("gene")
    eff = {g: ("" if r.EC_curated == "RETRACTED" else (r.EC_curated or r.EC)) for g, r in ann.iterrows()}
    loc = ann["localization"].to_dict()
    rec = dict(pd.read_csv(REC, sep="\t", dtype=str).fillna("")[["reaction_id", "EC"]].values)
    comp_loc = {"c": "Cytoplasm", "h": "Plastid", "u": "Plastid", "m": "Mitochondrion", "x": "Peroxisome",
                "r": "Endoplasmic", "g": "Golgi", "e": "Cell membrane", "n": "Nucleus"}
    tr = pd.read_excel(TRANSFER, sheet_name="GPR_transfer").fillna("")
    tr = tr[tr.verdict.str.startswith("shared, gene-less")]
    added = 0
    for _, row in tr.iterrows():
        r = m.reactions.get_by_id(row.reaction_id)
        if r.bounds == (0, 0) or "GENE" in r.gene_reaction_rule:
            continue
        e = r.annotation.get("ec-code")
        recs = set(e if isinstance(e, list) else [e] if e else []) or {x for x in rec.get(r.id, "").split(";") if x}
        if not recs:
            continue
        comps = {x.compartment for x in r.metabolites}
        keep = sorted({g for g in re.findall(r"GENE\d+", row.P_celeri_genes_via_P_renovo_proteome)
                       if g in m.genes or g in eff
                       if eff.get(g) and (eff[g] in recs or any(eff[g].rstrip(".-") == x.rstrip(".-") for x in recs))
                       and any(comp_loc.get(c) and comp_loc[c] in loc.get(g, "") for c in comps)})
        if keep:
            r.gene_reaction_rule = " or ".join(keep)
            added += 1
            note("E2", r.id, f"gene rule {r.gene_reaction_rule} (EC {';'.join(sorted(recs))})")
    return closed, added


def patch_F(m):
    """Storage-oil assembly and alternative electron sinks, all on P. celeri genes found by Pfam domain."""
    # F1 TAG synthase (DGAT, EC 2.3.1.20) on the three DGAT-family genes: GENE4063758 (DGAT2A, PF03982), GENE4067271
    # (DGAT1, PF03062 MBOAT) and GENE4071740 (DGAT2, PF03982, p1-only allele). EC tier 1, 4-5 independent evidence
    # classes each. The 143 TG* reactions sit in the plastid (iRC1080 convention) while DeepLoc places the three
    # enzymes in the ER (0.80-0.87); the rule is assigned regardless, because the enzyme-capacity bound does not
    # depend on where the model draws the compartment. GENE4068199 (LPLAT, 2.3.1.23), GENE4068792 (MGAT-like,
    # 2.3.1.-) and GENE4066711 (GUP1 MBOAT) carry the domain but not the DGAT EC and are not used.
    dgat = "GENE4063758 or GENE4067271 or GENE4071740"
    n = 0
    for r in m.reactions:
        if re.match(r"^TG\d", r.id) and not r.gene_reaction_rule:
            r.gene_reaction_rule = dgat
            n += 1
    note("F1", f"{n} TG* TAG-synthase reactions", f"gene rule {dgat} (EC 2.3.1.20, Pfam PF03982/PF03062)")
    # F2 flavodiiron protein (FLV) O2 photoreduction: 4 e- from thylakoid ferredoxin reduce O2 to water. Two plastid
    # FLV genes (PF01613, 'Diflavin flavoprotein A' GENE4067072 and GENE4067074); FLVA/FLVB act as a heterodimer, so
    # AND. Measured abundance 3.1e-4 and 2.2e-4 g gDW-1 in high CO2.
    add(m, "FLVh", {"fdxrd_u": -2, "o2_h": -1, "h_h": -4, "fdxox_u": 2, "h2o_h": 2},
        "GENE4067072 and GENE4067074", name="flavodiiron protein, ferredoxin:O2 oxidoreductase (Mehler-like)", step="F2")
    # F3 plastid terminal oxidase (PTOX, PF01786): plastoquinol + O2 -> plastoquinone + water. GENE4067969 (ubiquinol
    # oxidase, DeepLoc plastid); abundance 4.3e-5 (high CO2) vs 3.1e-4 g gDW-1 (air).
    add(m, "PTOXu", {"pqh2_u": -2, "o2_h": -1, "pq_u": 2, "h2o_h": 2}, "GENE4067969",
        name="plastid terminal oxidase, plastoquinol:O2 oxidoreductase", step="F3")


REF = f"{GEM}/p_celeri_GEM/iPicre1230.xml"
KEY = f"{GEM}/Pren_gene_key.tsv"
LIPID_SUBSYSTEMS = {"Fatty acid metabolism", "Glycerolipid metabolism", "Glycerophospholipid metabolism",
                    "Sphingolipid metabolism", "Biosynthesis of steroids", "Terpenoid metabolism", "Carotenoid biosynthesis"}
FAMILY_EC = [(r"^1[A-Z].*G3P_ACYLTRANSFERASE$", "2.3.1.15"), (r"^AGPAT", "2.3.1.51"), (r"^PGPS", "2.7.8.5"),
             (r"^PGPP", "3.1.3.27"), (r"^CDPDAGS", "2.7.7.41"), (r"^CDIPT", "2.7.8.11"), (r"^MGDG", "2.4.1.46"),
             (r"^DGDG", "2.4.1.241"), (r"^ETHAPT", "2.7.8.1"), (r"^GLUCRM", "2.4.1.80"), (r"^DHCER", "2.3.1.24"),
             (r"^CER\d.*_SYN$", "1.14.19.17"), (r"^DAG.*_ACHYD$", "3.1.1.23"), (r"^CLPN.*HYDROLASE$", "3.1.1.4"),
             (r"^CLPN.*_SYN$", "2.7.8.41"), (r"^PCHOL.*_SYN$", "2.7.8.2"),
             # G5: a third element True enforces the localization check even in a lipid subsystem
             (r"^FACOAL\d+(\(2OH\))?_ACP$", "2.3.1.179", True), (r"^FACOAL\d+(\(2OH\))?[ch]$", "6.2.1.3"),
             (r"^FA\d+\(2OH\)ACPHi", "3.1.2.14"), (r"^(TGTAG|SQDGATPA)", "2.3.1.20"), (r"^PETH[DT]?MT", "2.1.1.103"),
             (r"^mmcoaeCv$", "5.1.99.1"), (r"^pcoacCv$", "6.4.1.3"), (r"^DEDOLDPS$", "2.5.1.87"),
             (r"^DEDOLDPSAT$", "1.3.1.94"), (r"^BUP\d*$", "3.5.1.6")]
# G5 explicit rules (complexes and a dual-specificity enzyme; evidence in the docstring)
FAMILY_RULE = [(r"^PCHLDR$", "GENE2038481 and GENE2038482 and GENE2038489",
                "light-independent protochlorophyllide reductase ChlL/ChlN/ChlB, chloroplast-encoded, EC 1.3.7.7"),
               (r"^PCHOL.*_SYN$", "GENE4064348",
                "InterPro IPR014472 choline/ethanolamine phosphotransferase family (dual 2.7.8.1/2.7.8.2), the only member")]
DESATURASE = {9: "GENE4063770 or GENE4067065", 12: "GENE4062097 or GENE4066813", 7: "GENE4066008", "omega3": "GENE4067276"}


def _db_positions(met_id):
    """Double-bond positions encoded in an acyl id such as acoa2053Z6Z9Z12Z15Z or ac1419Z: chain, [positions]."""
    s = re.sub(r"_[a-z]$", "", met_id)
    mm = re.search(r"(\d{2})(\d)((?:\d+[ZE])+)$", s)
    if not mm:
        return None, []
    return int(mm.group(1)), [int(p) for p in re.findall(r"(\d+)[ZE]", mm.group(3))]


def patch_G(m):
    ref = cobra.io.read_sbml_model(REF)
    key = pd.read_csv(KEY, sep="\t", dtype=str).fillna("")
    p2c = {r.Pren_gene: re.findall(r"GENE\d+", r.P_celeri_RBH) for _, r in key.iterrows() if "GENE" in r.P_celeri_RBH}
    ann = pd.read_csv(ANNOT, sep="\t", dtype=str).fillna("")
    ann["gene"] = ann["gene_id"].str.split("_").str[0]
    ann = ann.drop_duplicates("gene").set_index("gene")
    eff = {g: ("" if r.EC_curated == "RETRACTED" else (r.EC_curated or r.EC)) for g, r in ann.iterrows()}
    conf, loc = ann["EC_confidence"].to_dict(), ann["localization"].to_dict()
    rec = dict(pd.read_csv(REC, sep="\t", dtype=str).fillna("")[["reaction_id", "EC"]].values)
    comp_loc = {"c": "Cytoplasm", "h": "Plastid", "u": "Plastid", "m": "Mitochondrion", "x": "Peroxisome",
                "r": "Endoplasmic", "g": "Golgi", "e": "Cell membrane", "n": "Nucleus"}
    ec2g = collections.defaultdict(set)
    for g, e in eff.items():
        if conf.get(g) == "High":
            for x in re.split(r"[;,\s]+", e):
                if x and x.count(".") == 3 and "-" not in x:
                    ec2g[x].add(g)

    def recs_of(r):
        e = r.annotation.get("ec-code")
        return set(e if isinstance(e, list) else [e] if e else []) or {x for x in rec.get(r.id, "").split(";") if x}

    def ec_ok(g, recs):
        e = eff.get(g, "")
        return bool(e) and any(e == x or e.rstrip(".-") == x.rstrip(".-") for x in recs)

    def loc_ok(g, r, force=False):
        if r.subsystem in LIPID_SUBSYSTEMS and not force:
            return True
        comps = {x.compartment for x in r.metabolites}
        return any(comp_loc.get(c) and comp_loc[c] in loc.get(g, "") for c in comps)

    def is_transport(r):
        comps = {x.compartment for x in r.metabolites}
        stems = {re.sub(r"_[a-z]$", "", x.id) for x in r.metabolites}
        return len(comps) > 1 and len(stems) < len(r.metabolites)

    fam = lambda rid: re.sub(r"[\d(\[].*$", "", rid)
    rule2cel = collections.defaultdict(set)
    for r in m.reactions:
        if "GENE" in r.gene_reaction_rule and "NF" not in r.gene_reaction_rule and r.id in ref.reactions:
            rr = ref.reactions.get_by_id(r.id).gene_reaction_rule
            if rr and len(fam(r.id)) >= 3:
                rule2cel[(fam(r.id), rr, frozenset(recs_of(r)))].add(r.gene_reaction_rule)
    counts = collections.Counter()
    for r in list(m.reactions):
        if r.boundary or r.gene_reaction_rule or r.bounds == (0, 0):
            continue
        recs = recs_of(r)
        if re.match(r"^12DGR.*G3P_ACYLTRANSFERASE$", r.id):
            # the source EC on these 1-acyl-G3P + acyl-CoA -> phosphatidate steps is 3.1.3.4 (the phosphatase), which
            # sent the PAP2 gene GENE4063154 onto the acyltransferase; the reaction is an LPAAT (2.3.1.51)
            recs = {"2.3.1.51"}
        route, force = "G3", False
        # G5: explicit rules for complexes and the dual-specificity enzyme
        for pat, rule, why in FAMILY_RULE:
            if re.match(pat, r.id):
                r.gene_reaction_rule = rule
                counts["G5"] += 1
                note("G5", r.id, f"gene rule {rule} ({why})")
                break
        if r.gene_reaction_rule:
            continue
        if not recs:  # lipid species-expanded reaction without EC: the family gives the EC, so every route is EC-checked
            for pat, ec, *flags in FAMILY_EC:
                if re.match(pat, r.id):
                    recs, route, force = {ec}, "G4", bool(flags and flags[0])
                    break
        ref_rule = ref.reactions.get_by_id(r.id).gene_reaction_rule if r.id in ref.reactions else ""
        pren = re.findall(r"Pren_g\d+", ref_rule)
        # G1: translate the source rule through the key
        if pren and (recs or len(set(pren)) <= 6):
            cel = sorted({g for p in pren for g in p2c.get(p, []) if (not recs or ec_ok(g, recs)) and loc_ok(g, r, force)})
            if cel:
                r.gene_reaction_rule = " or ".join(cel)
                counts["G1"] += 1
                note("G1", r.id, f"gene rule {r.gene_reaction_rule} (iPicre1230 rule via P. renovo key, "
                                 f"{'EC ' + ';'.join(sorted(recs)) + (' read from the reaction family' if route == 'G4' else '') if recs else 'orthology only'})")
                continue
        # G2: a sibling of the same family with the same source rule and EC already carries a P. celeri rule
        sib = rule2cel.get((fam(r.id), ref_rule, frozenset(recs))) if ref_rule and len(fam(r.id)) >= 3 else None
        if sib and len(sib) == 1:
            r.gene_reaction_rule = next(iter(sib))
            counts["G2"] += 1
            note("G2", r.id, f"gene rule {r.gene_reaction_rule} (copied from sibling {fam(r.id)}* with the same iPicre1230 rule)")
            continue
        # G3 / G4a: EC match in the annotation (G4 when the EC was read from the lipid reaction family)
        if recs and not is_transport(r) and not all(x.startswith(("3.6.3.", "7.")) for x in recs):
            cand = sorted({g for ec in recs for g in ec2g.get(ec, ()) if loc_ok(g, r, force)})
            if 0 < len(cand) <= 8:
                r.gene_reaction_rule = " or ".join(cand)
                counts[route] += 1
                note(route, r.id, f"gene rule {r.gene_reaction_rule} (EC {';'.join(sorted(recs))}"
                                  f"{' read from the reaction family' if route == 'G4' else ''}, high-confidence annotation)")
                continue
        # G4b: desaturases by the position of the new double bond
        if r.id.startswith("DESAT"):
            subs = [x for x, c in r.metabolites.items() if c < 0 and _db_positions(x.id)[0]]
            prods = [x for x, c in r.metabolites.items() if c > 0 and _db_positions(x.id)[0]]
            chain, new = None, []
            if prods:
                chain, pp = _db_positions(prods[0].id)
                ps = _db_positions(subs[0].id)[1] if subs else []
                new = [p for p in pp if p not in ps]
            if chain and len(new) == 1 and chain <= 22:
                p = new[0]
                rule = DESATURASE.get(p) or (DESATURASE["omega3"] if chain - p == 3 and chain >= 16 else None)
                if rule:
                    r.gene_reaction_rule = rule
                    counts["G4"] += 1
                    note("G4", r.id, f"gene rule {rule} (desaturase, new double bond delta-{p} on C{chain})")
                    continue
            counts["desaturase without P. celeri enzyme"] += 1
            continue
        counts["still gene-less"] += 1
    print("G: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    # what keeps the remaining reactions gene-less (any-confidence EC index, for the report)
    ec_any = collections.defaultdict(set)
    for g, e in eff.items():
        for x in re.split(r"[;,\s]+", e):
            if x and x.count(".") == 3 and "-" not in x:
                ec_any[x].add(g)
    rows = []
    for r in m.reactions:
        if r.boundary or r.gene_reaction_rule or r.bounds == (0, 0):
            continue
        recs = recs_of(r)
        for pat, ec, *flags in FAMILY_EC:
            if not recs and re.match(pat, r.id):
                recs = {ec}
                break
        if is_transport(r):
            cat = "transport reaction (no rule in the source model)"
        elif r.id.startswith("DESAT"):
            cat = "desaturation at a position without a P. celeri enzyme"
        elif not recs:
            cat = "no EC on the reaction or its family"
        elif not any(ec_any.get(x) for x in recs):
            cat = "EC absent from the P. celeri annotation at any confidence"
        elif not any(ec2g.get(x) for x in recs):
            cat = "EC only at low confidence in the annotation"
        else:
            cat = "high-confidence EC genes excluded (localization, generic ATPase EC, or more than 8 candidates)"
        rows.append(dict(reaction=r.id, subsystem=r.subsystem, family=fam(r.id), ECs=";".join(sorted(recs)), category=cat))
    pd.DataFrame(rows).to_csv(f"{GEM}/geneless_remaining.tsv", sep="\t", index=False)
    print("remaining gene-less by category: " + "; ".join(f"{k} {v}" for k, v in collections.Counter(x["category"] for x in rows).most_common()))
    return counts


H_CLOSE = {"ICL": "isocitrate lyase", "MALSm": "malate synthase", "MALSx": "malate synthase",
           "PFLACTm": "pyruvate formate lyase", "PFLh": "pyruvate formate lyase", "argCv": "arginase"}
AA20 = ("ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile", "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val")
# pigment class = chlorophyll (LaPanse 2024 Fig 2: 4.4 replete, 2.1 at 40 ppm) x (1 + 0.203 carotenoids per chlorophyll, Cano 2024 Fig 1D): replete 5.3, 40 ppm 2.5
# (the carotenoid:chlorophyll ratio is assumed unchanged at 40 ppm, not measured; Cano's HPLC chlorophyll of 6.2 % AFDW and total pigment of 7.5 % are on a different basis
# from LaPanse's spectrophotometric chlorophyll and are not mixed with it)
CORE = {"Biomass_replete_LaPanse": dict(protein=51.3, carbohydrate=12.2, lipid=11.0, pigment=5.3, nucleic=5.7),
        "Biomass_N40_LaPanse": dict(protein=19.9, carbohydrate=10.8, lipid=13.2, pigment=2.5, nucleic=2.2)}
NGAM_LB = 0.62


def biomass_class(x):
    i = x.id
    if re.match(r"^(%s)(__L)?_[a-z]$" % "|".join(AA20), i):
        return "protein"
    if re.match(r"^d?[acgtu]tp_", i):
        return "nucleic"
    if re.match(r"^(12dgr|tag|mag|mgdg|dgdg|sqdg|pg|pe|pchol|pc|dgts|pail|crm)\d", i):
        return "lipid"
    if i in ("man_c", "gal_c", "rib-D_c"):
        return "carbohydrate"
    if re.match(r"^(chla|chlb|acaro|anxan|caro|loroxan|lut|neoxan|vioxan|zaxan)_", i):
        return "pigment"
    return None


def patch_H(m):
    for rid, what in H_CLOSE.items():
        r = m.reactions.get_by_id(rid)
        note("H1", rid, f"closed (was {r.bounds}, rule '{r.gene_reaction_rule}'): {what} reported absent in Picochlorum (Bec 2025, "
                        "Xiong 2025) and not found in the P. celeri annotation by Pfam, EC or name")
        r.bounds = (0, 0)
    add(m, "NGAM", {"atp_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1}, lb=NGAM_LB, ub=1000.0,
        name="non-growth ATP maintenance (Picochlorum, Barten 2022 maintenance rate converted to ATP)", step="H2")
    b = m.reactions.get_by_id(BIO)
    st = m.metabolites.get_by_id("starch300_h")
    gam = {m.metabolites.get_by_id(k) for k in ("atp_c", "h2o_c")}
    cls = {x: biomass_class(x) for x, c in b.metabolites.items() if c < 0 and x not in gam and x is not st}
    g = {}
    for x, k in cls.items():
        if k:
            g[k] = g.get(k, 0) + (-b.metabolites[x] * x.formula_weight / 1000)
    for rid, frac in CORE.items():
        tot = sum(frac.values())
        target = {k: v / tot for k, v in frac.items()}
        core_g = sum(g.values())
        stoich = {x: c for x, c in b.metabolites.items() if x is not st}
        for x, k in cls.items():
            if k in target:
                stoich[x] = b.metabolites[x] * target[k] * core_g / g[k]
        r = cobra.Reaction(rid, name=f"biomass on the measured composition: {frac}", lower_bound=0.0, upper_bound=1000.0)
        m.add_reactions([r])
        r.add_metabolites(stoich)
        note("H3", rid, f"added: biomass with class fractions {frac} (% AFDW before renormalisation), starch excluded, "
                        "growth-associated ATP and non-class components as in the native biomass")


def patch_I(m):
    """Ammonium routes between compartments, plastid glutamine synthetase and ferredoxin-GOGAT (see docstring, section I)."""
    # proton-free: the iRC1080/iPicre1230 form (H+ + NH4+ moved together) lets the model run a urea -> NH4+ export cycle on a
    # free proton shuttle (biomass yield 2.15 g/mol, negative net Ci uptake), so NH3/NH4+ diffusion is written without H+
    add(m, "NH4th", {"nh4_c": -1, "nh4_h": 1}, lb=-1000.0, name="ammonium (NH3) diffusion, cytosol-plastid", step="I1")
    add(m, "NH4tm", {"nh4_c": -1, "nh4_m": 1}, lb=-1000.0, name="ammonium (NH3) diffusion, cytosol-mitochondrion", step="I1")
    add(m, "AKGMALth", {"akg_c": -1, "mal__L_h": -1, "akg_h": 1, "mal__L_c": 1}, m.reactions.get_by_id("MALGLUth").gene_reaction_rule, lb=-1000.0,
        name="2-oxoglutarate/malate translocator, plastid envelope (iPicre1230 AKGMALth, same transporter family rule as MALGLUth)", step="I4")
    add(m, "GALh", {"atp_h": -1, "glu__L_h": -1, "nh4_h": -1, "adp_h": 1, "gln__L_h": 1, "pi_h": 1}, "GENE4068506",
        name="glutamine synthetase, chloroplastic (GS2)", step="I2")
    add(m, "GLUS(ferr)", {"akg_h": -1, "fdxrd_h": -2, "gln__L_h": -1, "h_h": -2, "fdxox_h": 2, "glu__L_h": 2}, "GENE4068118",
        name="glutamate synthase (ferredoxin), chloroplastic", step="I3")
    r = m.reactions.get_by_id("GLUS(nadh)")
    note("I3", r.id, f"gene rule '{r.gene_reaction_rule}' -> GENE4062926 (glutamate synthase 1 [NADH], chloroplastic, EC 1.4.1.14)")
    r.gene_reaction_rule = "GENE4062926"


AAA_FILE = "/home/yasin/Claudio/Gen_Trans/final_analysis/Bioreactor_data/Protein/aaa_composition_mol_pct.tsv"
AA_GROUPS = {"Asx": ["asp__L_c", "asn__L_c"], "Thr": ["thr__L_c"], "Ser": ["ser__L_c"], "Glx": ["glu__L_c", "gln__L_c"], "Pro": ["pro__L_c"],
             "Gly": ["gly_c"], "Ala": ["ala__L_c"], "Cys": ["cys__L_c"], "Val": ["val__L_c"], "Met": ["met__L_c"], "Ile": ["ile__L_c"],
             "Leu": ["leu__L_c"], "Tyr": ["tyr__L_c"], "Phe": ["phe__L_c"], "His": ["his__L_c"], "Lys": ["lys__L_c"], "Trp": ["trp__L_c"],
             "Arg": ["arg__L_c"]}
PIG = dict(chl_ab=4.4, car_over_chl=0.203, lutein=0.103, bcar=0.048, neox=0.020)   # Cano 2024 Fig 1D/E (TG2 WT, ALGiSIM, HPLC): total carotenoids 20.3 % of total chlorophyll mass, lutein 10.3, beta-carotene 4.8, neoxanthin 2.0 % of chlorophyll; chlorophyll itself is taken from LaPanse 2024 (4.4 % replete, 2.1 % at 40 ppm), the same basis as the rest of the composition
# acyl-chain targets, % of acyl groups (docstring J3): measured C16:0, C18:x; C16:2 3.2 (LaPanse 2024), C16:1 >= C16:2 (LaPanse ordering, set equal);
# the rest of the lab GC profile that is unassigned (about 27 %, probably C16:3 and C16:4) goes to C16:3, the only C16 polyunsaturate with a
# complete P. celeri desaturase set (delta-4, needed for C16:4, has no enzyme; open item)
# PE and PI from the lab's own EMSL lipidomics of the WT in high CO2 (Gen_Trans/final_analysis/lipidome/data/61507-Posewitz_Lipidomics_Results.xlsx, negative mode, mean of vessels B1/B3/B4):
# PG 67.7, PE 23.1, PI 8.1 % of the signal, so PE/PG 0.341 and PI/PG 0.120, anchored on the measured PG of 21.9 % of total lipid (Grubisic 2022): PE 7.5, PI 2.6 % of lipid mass.
# DGDG 17.7 = DGDG/MGDG 0.79 (positive mode, same [M+NH4]+ family) x MGDG 22.5; SQDG 5.6 = SQDG/MGDG 0.25 (less reliable: SQDG ionises best in negative mode); DGTS 1.2 = DGTS/PC 0.146 x PC 8.1
# (DGTS and PC share the quaternary-ammonium positive-mode response).
# mag capped at 1 % (a minor pool, assumption). 11 of 36 DGTS species remain unproducible, so their
# mass is distributed pro rata over the modelled classes by the final normalisation of the class targets.
# Ratios inside one ionisation family only (phospholipids deprotonated in negative mode); positive-mode class shares are not comparable across classes (PC 71 % of the positive signal).
CLASS_TARGET = {"mgdg": 22.5, "pg": 21.9, "pchol": 8.1, "tag": 5.1, "12dgr": 3.0, "pe": 7.5, "pail": 2.6, "mag": 1.0, "dgdg": 17.7, "sqdg": 5.6, "dgts": 1.2}   # % of lipid mass: MGDG, PG (21.93), PC, TAG from Picochlorum sp. D3 (Grubisic 2022 Section 2.4: MGDG 22.48, PG 21.93, PC 8.08 % of total lipid by TLC-FID; TAG from Table 3; other strain and medium); free DAG capped at 3 % (assumption: a minor pool, the pro-rata fill had made it 24 %)
FA_TARGET = {"C16:0": 22.0, "C18:3": 31.0, "C18:2": 11.0, "C18:1": 2.0, "C18:0": 1.0, "C16:1": 18.1, "C16:2": 3.2, "C16:3": 11.7}   # the unassigned third of the GC profile (33 %): C16:3 11.7 % of FAME from Picochlorum atomus (von Alvensleben 2013 Table 4, replete), C16:2 3.2 % (LaPanse 2024), C16:1 the remainder
_POS = re.compile(r"\d+[ZE]")


def lipid_chains(sid):
    """(class, [chain types]) of a lipid species id such as pchol1601826Z9Z_c; (None, None) if not a glycerolipid.
    Each acyl code is the chain length (2 digits, or 1 digit for C4 to C9), 1 digit of double bonds d, then exactly d positions
    (digits + Z/E), the only way to split codes such as 1839Z12Z15Z2249Z12Z15Z18Z into C18:3 and C22:4 or 401001133Z6Z9Z into
    C4:0, C10:0 and C11:3. Backtracks over the 1- and 2-digit chain lengths and returns the first complete parse."""
    s_ = re.sub(r"_[a-z]$", "", sid)
    mm = re.match(r"^(12dgr|tag|mgdg|dgdg|sqdg|pg|pe|pchol|pail|mag|dgts)(\d.*)$", s_)
    if not mm:
        return None, None
    cls, rest = mm.groups()

    def go(i):
        if i == len(rest):
            return []
        for width in (2, 1):
            if i + width + 1 > len(rest):
                continue
            if not (rest[i:i + width].isdigit() and rest[i + width].isdigit()):
                continue
            ch, d, j = int(rest[i:i + width]), int(rest[i + width]), i + width + 1
            if not 3 <= ch <= 40 or d > 8:
                continue
            ok, pos = True, []
            for _ in range(d):
                t = _POS.match(rest, j)
                if not t or int(t.group(0)[:-1]) > ch:
                    ok = False
                    break
                pos.append(t.group(0))
                j = t.end()
            if not ok:
                continue
            tail = go(j)
            if tail is not None:
                return [f"C{ch}:{d}"] + tail
        return None

    out = go(0)
    return (cls, out) if out is not None else (cls, None)


def patch_L(m):
    """L. Galactolipid and sulfolipid synthesis for the lipid classes the lab's EMSL lipidomics detects in the P. celeri WT (DGDG 24 species, SQDG 10, DGTS 22) but the
    model could not make (0 of 17 DGDG and 0 of 7 SQDG species producible). Genes found by EC, KO and name in the annotation (Pceleri_annotation.tsv): DGD1 digalactosyldiacylglycerol
    synthase (EC 2.4.1.241, K09480) GENE4062052 and GENE4065431; UDP-sulfoquinovose synthase SQD1 GENE4065148 (UDPSQS already in the model); SQDG synthase SQD2 GENE4068253
    (K06119; the two putative copies GENE4063927 and GENE4066131 are unpaired partial models and are not used); DGTS synthase BTA1 GENE4064750 (K13621; its reactions exist).
    For every producible MGDG species a DGDG species of the same acyl code is created (formula = MGDG + UDP-galactose - UDP) with MGDG + UDP-galactose -> DGDG + UDP;
    for every producible C16/C18 diacylglycerol without one an SQDG synthase reaction DAG + UDP-sulfoquinovose -> SQDG + UDP is added the same way."""
    M = m.metabolites
    ph = m.reactions.get_by_id("EX_photonVis(e)")

    def producible(x):
        with m:
            ph.lower_bound = -200.0
            d = m.add_boundary(x, type="demand")
            m.objective = d
            return m.slim_optimize(error_value=0.0) > 1e-6

    def formula_sum(parts):
        el = {}
        for x, sign in parts:
            for k, v in x.elements.items():
                el[k] = el.get(k, 0) + sign * v
        return "".join(f"{k}{int(v)}" for k, v in sorted(el.items()) if v), sum(sign * x.charge for x, sign in parts)

    # L0 prokaryotic-type plastid diacylglycerol. The plastid PA with 18:x at sn-1 and 16:0 at sn-2 (ATS1 then ATS2 order; AGPAT1819Z160h, AGPAT18111Z160h) fed CDP-DAG only and
    # had no phosphatase, so the 18/16 diacylglycerols (the main MGDG species of a prokaryotic-pathway chloroplast) could not be made and 18 of 36 MGDG species
    # and the DGDG made from them were unproducible. The plastid PAP (same gene and stoichiometry as the existing 12dgr*_PHOSPHATIDATEPHOSPHATASE1 reactions) is added.
    tpap = m.reactions.get_by_id("12dgr1601819Z_PHOSPHATIDATEPHOSPHATASE1")
    n_pap = 0
    for x in [y for y in M if re.match(r"^pa(1819Z|18111Z)160_h$", y.id)]:
        code = x.id[2:-2]
        d = f"12dgr{code}_h"
        if d in M and not any(r.id.startswith("PAPH") for r in M.get_by_id(d).reactions) and f"PAPH{code}" not in m.reactions:
            add(m, f"PAPH{code}", {x.id: -1, "h2o_h": -1, d: 1, "pi_h": 1}, tpap.gene_reaction_rule, 0.0,
                name=f"plastid phosphatidate phosphatase ({code})", step="L0")
            n_pap += 1
    note("L", "L0", f"{n_pap} plastid phosphatidate phosphatase reactions added for 18:x/16:0 PA")
    udpgal, udp, udpsq = M.get_by_id("udpgal_h"), M.get_by_id("udp_h"), M.get_by_id("udpsq_h")
    n_dg = n_sq = 0
    for x in [y for y in M if re.match(r"^mgdg.+_h$", y.id)]:
        code = x.id[4:-2]
        if not producible(x):
            continue
        did = f"dgdg{code}_h"
        if did not in M:
            f, ch = formula_sum([(x, 1), (udpgal, 1), (udp, -1)])
            d = cobra.Metabolite(did, name=f"Digalactosyldiacylglycerol ({code})", formula=f, charge=ch, compartment="h")
            m.add_metabolites([d])
        if f"DGDGS{code}" not in m.reactions:
            add(m, f"DGDGS{code}", {x.id: -1, "udpgal_h": -1, did: 1, "udp_h": 1}, "GENE4062052 or GENE4065431", 0.0,
                name=f"digalactosyldiacylglycerol synthase DGD1 ({code})", step="L1")
            n_dg += 1
    for x in [y for y in M if re.match(r"^12dgr.+_h$", y.id)]:
        cls_, ch_ = lipid_chains(x.id)
        if not ch_ or not all(t in ("C16:0", "C16:1", "C16:2", "C16:3", "C18:0", "C18:1", "C18:2", "C18:3") for t in ch_):
            continue
        code = x.id[5:-2]
        if f"SQDGS{code}" in m.reactions or not producible(x):
            continue
        sid = f"sqdg{code}_h"
        if sid not in M:
            f, ch = formula_sum([(x, 1), (udpsq, 1), (udp, -1)])
            m.add_metabolites([cobra.Metabolite(sid, name=f"Sulfoquinovosyldiacylglycerol ({code})", formula=f, charge=ch, compartment="h")])
        add(m, f"SQDGS{code}", {x.id: -1, "udpsq_h": -1, sid: 1, "udp_h": 1}, "GENE4068253", 0.0,
            name=f"sulfoquinovosyldiacylglycerol synthase SQD2 ({code})", step="L2")
        n_sq += 1
    # L3 DGTS: BTA1 (GENE4064750) makes DGTS from DAG + SAM (diacylglyceryl-homoserine, DGHS) and three further SAM (trimethylation). The native model carried
    # these two steps for six 16/18 diacylglycerols only; they are added for every producible cytosolic C16/C18 diacylglycerol, formulas from the native DGHS and DGTS.
    n_ts = 0
    tdghs, tdgts = M.get_by_id("dghs1601819Z_c"), M.get_by_id("dgts1601819Z_c")
    tdag = M.get_by_id("12dgr1601819Z_c")
    for x in [y for y in M if re.match(r"^12dgr.+_c$", y.id)]:
        cls_, ch_ = lipid_chains(x.id)
        if not ch_ or not all(t in ("C16:0", "C16:1", "C16:2", "C16:3", "C18:0", "C18:1", "C18:2", "C18:3") for t in ch_):
            continue
        code = x.id[5:-2]
        if f"ACPT{code}" in m.reactions or not producible(x):
            continue
        dh, dt = f"dghs{code}_c", f"dgts{code}_c"
        for mid, tmpl, nm in ((dh, tdghs, "Diacylglyceryl-homoserine"), (dt, tdgts, "Diacylglyceryl-trimethylhomoserine")):
            if mid not in M:
                el = dict(tmpl.elements)
                for k, v in x.elements.items():
                    el[k] = el.get(k, 0) + v
                for k, v in tdag.elements.items():
                    el[k] = el.get(k, 0) - v
                m.add_metabolites([cobra.Metabolite(mid, name=f"{nm} ({code})", formula="".join(f"{k}{int(v)}" for k, v in sorted(el.items()) if v),
                                                    charge=tmpl.charge, compartment="c")])
        add(m, f"ACPT{code}", {x.id: -1, "amet_c": -1, "5mta_c": 1, dh: 1}, "GENE4064750", 0.0, name=f"BTA1 diacylglyceryl homoserine synthase ({code})", step="L3")
        add(m, f"TM{code}", {dh: -1, "amet_c": -3, "ahcys_c": 3, dt: 1}, "GENE4064750", 0.0, name=f"BTA1 DGHS trimethylation ({code})", step="L3")
        n_ts += 1
    note("L", "summary", f"{n_ts} DGTS (BTA1) species added; DGTS species producible afterwards: {sum(producible(y) for y in M if re.match(r'^dgts.+_c$', y.id))}")
    note("L", "summary", f"{n_dg} DGDG synthase (DGD1) and {n_sq} SQDG synthase (SQD2) reactions added; DGDG and SQDG species producible afterwards: "
                          f"{sum(producible(y) for y in M if re.match(r'^dgdg.+_h$', y.id))} and {sum(producible(y) for y in M if re.match(r'^sqdg.+_h$', y.id))}")


def _flux_producible(m, x):
    """True when the network can actually make metabolite x (a demand reaction carries flux with light available); a producing reaction alone is not enough."""
    with m:
        m.reactions.get_by_id("EX_photonVis(e)").lower_bound = -200.0
        d = m.add_boundary(x, type="demand")
        m.objective = d
        return m.slim_optimize(error_value=0.0) > 1e-6


def patch_J(m):
    import numpy as np
    b = m.reactions.get_by_id(BIO)
    native = {x: c for x, c in b.metabolites.items()}
    ghost = cobra.Reaction("Biomass_Niaz_native", name="biomass before the 2026-10-20 curation (kept closed)", lower_bound=0.0, upper_bound=0.0)
    m.add_reactions([ghost])
    ghost.add_metabolites(native)
    note("J4", ghost.id, "added closed: the native biomass stoichiometry before step J")
    mw = lambda x: x.formula_weight / 1000.0
    # ---- J1 amino acids
    meas = pd.read_csv(AAA_FILE, sep="\t")
    meas = meas[meas.condition == "CO2"][list(AA_GROUPS)].mean()
    mol = {k: sum(-native.get(m.metabolites.get_by_id(i), 0) for i in ids) for k, ids in AA_GROUPS.items()}
    keep = ("Cys", "Trp")
    free_share = 1 - sum(mol[k] for k in keep) / sum(mol.values())
    new = {k: (mol[k] / sum(mol.values())) if k in keep else meas[k] / sum(v for kk, v in meas.items() if kk not in keep) * free_share for k in AA_GROUPS}
    grams = sum(-native[m.metabolites.get_by_id(i)] * mw(m.metabolites.get_by_id(i)) for ids in AA_GROUPS.values() for i in ids if m.metabolites.get_by_id(i) in native)
    per_mol = {}
    coef = {}
    for k, ids in AA_GROUPS.items():
        inner = [-native.get(m.metabolites.get_by_id(i), 0) for i in ids]
        w = [x / sum(inner) if sum(inner) else 1 / len(ids) for x in inner]
        for i, wi in zip(ids, w):
            coef[i] = new[k] * wi
    g_per_unit = sum(coef[i] * mw(m.metabolites.get_by_id(i)) for i in coef)
    scale = grams / g_per_unit
    for i, c in coef.items():
        x = m.metabolites.get_by_id(i)
        b.add_metabolites({x: -(c * scale) - b.metabolites.get(x, 0)}, combine=True)
    arg = -b.metabolites[m.metabolites.arg__L_c]
    note("J1", BIO, f"amino-acid profile set to the measured mol % (Arg {100 * coef['arg__L_c']:.2f} mol %, was {100 * mol['Arg'] / sum(mol.values()):.2f}); protein mass kept; "
                    f"Cys and Trp kept at the native share")
    # ---- J2 pigments
    pig_ids = ["chla_u", "chlb_u", "acaro_h", "anxan_u", "caro_u", "loroxan_u", "lut_u", "neoxan_u", "vioxan_u", "zaxan_u"]
    pm = {i: m.metabolites.get_by_id(i) for i in pig_ids}
    pg_native = {i: -native.get(pm[i], 0) * mw(pm[i]) for i in pig_ids}
    P = sum(pg_native.values())
    chl = P / (1 + PIG["car_over_chl"])
    target = {"chla_u": chl * PIG["chl_ab"] / (PIG["chl_ab"] + 1), "chlb_u": chl / (PIG["chl_ab"] + 1), "lut_u": PIG["lutein"] * chl,
              "caro_u": PIG["bcar"] * chl, "neoxan_u": PIG["neox"] * chl}
    rest = PIG["car_over_chl"] * chl - sum(target[k] for k in ("lut_u", "caro_u", "neoxan_u"))
    others = ["acaro_h", "anxan_u", "loroxan_u", "vioxan_u", "zaxan_u"]
    so = sum(pg_native[i] for i in others)
    for i in others:
        target[i] = rest * pg_native[i] / so
    for i, g_ in target.items():
        b.add_metabolites({pm[i]: -(g_ / mw(pm[i])) - b.metabolites.get(pm[i], 0)}, combine=True)
    note("J2", BIO, f"pigments: Chl a/b {PIG['chl_ab']}, carotenoids {PIG['car_over_chl']:.2f} of chlorophyll mass, lutein {PIG['lutein']}, beta-carotene {PIG['bcar']}, "
                    f"neoxanthin {PIG['neox']} (fractions of chlorophyll mass); pigment mass kept")
    # ---- J3 lipid species, iterative proportional fitting
    sp = {}
    for x, c in b.metabolites.items():
        if c < 0:
            cls, ch = lipid_chains(x.id)
            if cls and ch:
                sp[x] = (cls, ch, -c)
    keepers = {x: v for x, v in sp.items() if all(t.split(":")[0] in ("C16", "C18") for t in v[1])}
    present = {(v[0], x.compartment) for x, v in sp.items()}
    # classes the lab lipidomics detects that the native biomass lacks: their producible species may join the fit (step L made DGDG and SQDG producible)
    for cls_, comp_ in (("dgdg", "h"), ("sqdg", "h"), ("dgts", "c")):
        present.add((cls_, comp_))
    pool_extra = {}
    for x in m.metabolites:
        cls_, ch_ = lipid_chains(x.id)
        if not (ch_ and (cls_, x.compartment) in present and x not in keepers and all(t in FA_TARGET or t in ("C16:1", "C16:2", "C16:3") for t in ch_)):
            continue
        # only species the network can make: some open reaction other than the biomass reactions produces it
        if any(r.metabolites[x] > 0 and r.bounds != (0, 0) and not r.id.startswith("Biomass") for r in x.reactions) and _flux_producible(m, x):
            pool_extra[x] = (cls_, ch_, 1e-6)
    total_g = sum(v[2] * mw(x) for x, v in sp.items())
    dropped_g = total_g - sum(v[2] * mw(x) for x, v in keepers.items())
    cls_g = collections.defaultdict(float)
    for x, v in sp.items():
        cls_g[v[0]] += v[2] * mw(x)
    live = {v[0] for v in keepers.values()}
    lost_classes = sorted(set(cls_g) - live)
    live_g = sum(cls_g[c] for c in live)
    cls_target = {c: cls_g[c] * total_g / live_g for c in live}  # classes without a C16/C18 species hand their mass to the others pro rata
    types = sorted({t for v in keepers.values() for t in v[1]})
    grp = lambda t: t
    tgt = {t: v for t, v in FA_TARGET.items() if t in {grp(u) for u in types}}
    tsum = sum(tgt.values())
    tgt = {t: 100.0 * v / tsum for t, v in tgt.items()}
    keepers.update(pool_extra)
    w = {x: v[2] for x, v in keepers.items()}
    for x, v in keepers.items():
        if not all(grp(t) in tgt for t in v[1]):
            w[x] = 0.0  # a chain with no target cannot be carried
    def class_restore():
        """Soft pull of the three classes with Picochlorum data (MGDG, PC, TAG) towards their measured mass shares; the other classes keep
        their relative native shares and fill the remainder."""
        cg = collections.defaultdict(float)
        for x, v in keepers.items():
            cg[v[0]] += w[x] * mw(x)
        tot_g = sum(cg.values())
        known = {c: CLASS_TARGET[c] / 100.0 * tot_g for c in CLASS_TARGET if cg.get(c, 0) > 0}
        rest_native = sum(cls_target[c] for c in cg if c not in known)
        rest_now = sum(cg[c] for c in cg if c not in known)
        if rest_native <= 0 and known:
            # every present class has a target: the targets are renormalised to the total (the species the model cannot make are spread pro rata)
            kt = sum(known.values())
            known = {c: v * tot_g / kt for c, v in known.items()}
        for x, v in keepers.items():
            c = v[0]
            if cg[c] <= 0:
                continue
            tgt_c = known[c] if c in known else (tot_g - sum(known.values())) * cls_target[c] / rest_native
            w[x] *= (tgt_c / cg[c]) ** 0.5

    def chain_step():
        cnt = collections.defaultdict(float); tot = 0.0
        for x, v in keepers.items():
            for t in v[1]:
                cnt[grp(t)] += w[x]; tot += w[x]
        for x, v in keepers.items():
            if w[x] > 0:
                f = np.mean([(tgt[grp(t)] / 100.0) / max(cnt[grp(t)] / tot, 1e-12) for t in v[1]])
                w[x] *= f ** 0.5
        return cnt, tot

    for it in range(800):
        class_restore()
        chain_step()  # multiplicative updates from the native weights: the least-change (minimum KL) profile that meets the chain targets
    # the last operation is the chain fit, so the acyl profile is what the measurements ask for; class shares may drift a little
    # the multiplicative fit leaves unused species at 1e-19 to 1e-21 mmol per g; coefficients that small make the LPs numerically unsolvable
    # (GLPK cycles on them), so species below 1e-6 mmol gDW-1 (the smallest native coefficient is 6.5e-7) are dropped before rescaling
    pruned = [x for x in keepers if 0 < w[x] < 1e-6 * (total_g / sum(w[y] * mw(y) for y in keepers))]
    for x in pruned:
        w[x] = 0.0
    cgf = collections.defaultdict(float)
    for x, v in keepers.items():
        cgf[v[0]] += w[x] * mw(x)
    gsum = sum(cgf.values())
    for x in keepers:
        w[x] *= total_g / gsum
    cgf = collections.defaultdict(float)
    for x, v in keepers.items():
        cgf[v[0]] += w[x] * mw(x)
    drift = max(abs(cgf.get(c, 0.0) / cls_target[c] - 1) for c in cls_target)
    # write the fitted weights back for the native species AND for the producible species the fit added from the model pool (they were fitted but never written before)
    for x in list(sp) + [y for y in keepers if y not in sp]:
        b.add_metabolites({x: -(w.get(x, 0.0)) - b.metabolites.get(x, 0)}, combine=True)
    cnt = collections.defaultdict(float); tot = 0.0
    for x, v in keepers.items():
        for t in v[1]:
            cnt[t] += w[x]; tot += w[x]
    prof = ", ".join(f"{t} {100 * v / tot:.1f}" for t, v in sorted(cnt.items(), key=lambda z: -z[1]) if v > 0)
    removed = sum(1 for x in sp if w.get(x, 0.0) == 0.0)
    note("J3", BIO, f"lipid species restricted to C16/C18 chains with a measured target ({removed} of {len(sp)} species removed, {100 * dropped_g / total_g:.0f} % of the lipid mass was in other chains) "
                    f"and re-weighted to the measured acyl profile (targets renormalised over the chains the model can carry); result (% of acyl groups): {prof}; class mass shares "
                    f"moved by at most {100 * drift:.0f} %, classes without a usable species ({lost_classes}) hand their mass to the others pro rata")
    print("J3 acyl profile after fit:", prof, "| classes without usable species:", lost_classes)
    return ghost


def patch_J_default(m):
    """J4: Biomass_Picre_Auto becomes the replete measured composition; the starch term of the native reaction is kept."""
    b = m.reactions.get_by_id(BIO)
    rep = m.reactions.get_by_id("Biomass_replete_LaPanse")
    st = m.metabolites.get_by_id("starch300_h")
    starch = b.metabolites.get(st, 0)
    old = dict(b.metabolites)
    b.add_metabolites({x: -c for x, c in old.items()})
    b.add_metabolites(dict(rep.metabolites))
    if starch:
        b.add_metabolites({st: starch})
    note("J4", BIO, f"default biomass is now the replete LaPanse-composition biomass with the curated micro-composition (starch term {starch:.4f} kept)")


def patch_N(m):
    """Plastid acetyl-CoA carboxylase gene rules (step N). The draft rules for the three plastid reactions OR together every
    biotin-dependent carboxylase gene in the genome (pyruvate carboxylase GENE4066235, propionyl-CoA carboxylase GENE4067618,
    urea carboxylase GENE4066869, cytosolic acetyl-CoA carboxylase 1 GENE4066936) with the plastid subunits, so the enzyme
    capacity of the plastid malonyl-CoA step was taken from unrelated enzymes. The annotation has the heteromeric plastid enzyme:
    biotin carboxyl carrier protein GENE4068291, biotin carboxylase GENE4065894 (EC 6.3.4.14), carboxyl transferase alpha
    GENE4065316 and the two chloroplast-encoded carboxyl transferase beta genes GENE2038492 and GENE2038442."""
    rules = {"BTNC": "GENE4065894 and GENE4068291",
             "ACCOAC": "GENE4065316 and (GENE2038492 or GENE2038442)",
             "ACCOAhi": "GENE4068291 and GENE4065894 and GENE4065316 and (GENE2038492 or GENE2038442)"}
    for rid, rule in rules.items():
        if rid in m.reactions:
            old = m.reactions.get_by_id(rid).gene_reaction_rule
            m.reactions.get_by_id(rid).gene_reaction_rule = rule
            note("N1", rid, f"gene rule restricted to the plastid acetyl-CoA carboxylase subunits (was {old[:80]}...)")


def patch_M(m):
    """Plastid fatty acid synthesis and acetyl-CoA supply (step M).
    M1: the draft makes every de novo fatty acid in the cytosol (60 reactions with cytosolic acyl-carrier-protein species
        carrying the chloroplast 3-oxoacyl-ACP synthase II gene GENE4064434) while the plastid copy stays idle; in green algae
        the type II synthase is plastidial, so the cytosolic chain is closed.
    M2: the model has no pyruvate transport into the plastid, so plastid acetyl-CoA cost 110 photons per unit against 23 in
        the cytosol; a proton-neutral exchange is added with the plastid BASS-family transporter genes (GENE4063079,
        GENE4068003), the family member that carries pyruvate in plant plastids (candidate assignment).
    M3: the gene-less plastid pyruvate:ferredoxin oxidoreductase PYRShi is closed, because the P. celeri annotation has no
        EC 1.2.7.1, and plastid pyruvate dehydrogenase is measured (Steichen 2024)."""
    M = m.metabolites
    cyt = [r for r in m.reactions if any(x.id.endswith("ACP_c") for x in r.metabolites) and (r.lower_bound != 0 or r.upper_bound != 0)]
    for r in cyt:
        r.bounds = (0.0, 0.0)
        note("M1", r.id, "closed: cytosolic acyl-carrier-protein fatty acid synthesis duplicates the plastid pathway")
    if "PYRth" not in m.reactions:
        add(m, "PYRth", {"pyr_c": -1, "pyr_h": 1}, gpr="GENE4063079 or GENE4068003", lb=-1000.0, ub=1000.0,
            name="pyruvate transport, cytosol to chloroplast", step="M2")
    if "PYRShi" in m.reactions and (m.reactions.PYRShi.lower_bound != 0 or m.reactions.PYRShi.upper_bound != 0):
        m.reactions.PYRShi.bounds = (0.0, 0.0)
        note("M3", "PYRShi", "closed: gene-less plastid pyruvate:ferredoxin oxidoreductase, no EC 1.2.7.1 in the annotation")


def patch(m):
    patch_A(m)
    patch_B(m)
    patch_C(m)
    patch_D(m)
    out = patch_E(m)
    patch_F(m)
    patch_G(m)
    patch_I(m)
    patch_M(m)
    patch_N(m)
    patch_L(m)
    patch_J(m)
    patch_H(m)
    patch_J_default(m)
    return out


def energy_cycle_check(m):
    """With every uptake closed, no mass-balanced drain of an energy or redox carrier may carry flux."""
    M = m.metabolites
    tests = {}
    for c in "chmxru":
        has = lambda *ids: all(f"{i}_{c}" in M for i in ids)
        for n in ("atp", "gtp", "ctp", "utp"):
            d = n[0] + "dp"
            if has(n, d, "pi", "h2o"):
                tests[f"{n.upper()}_{c}"] = {f"{n}_{c}": -1, f"h2o_{c}": -1, f"{d}_{c}": 1, f"pi_{c}": 1}
        for red, ox in (("nadh", "nad"), ("nadph", "nadp"), ("fadh2", "fad"), ("q8h2", "q8"), ("pqh2", "pq")):
            if has(red, ox, "h"):
                tests[f"{red}_{c}"] = {f"{red}_{c}": -1, f"{ox}_{c}": 1, f"h_{c}": 2}
        if has("fdxrd", "fdxox"):
            tests[f"fdx_{c}"] = {f"fdxrd_{c}": -1, f"fdxox_{c}": 1}
        if has("accoa", "coa", "ac", "h2o"):
            tests[f"accoa_{c}"] = {f"accoa_{c}": -1, f"h2o_{c}": -1, f"coa_{c}": 1, f"ac_{c}": 1}
    out = {}
    with m:
        if "NGAM" in m.reactions:
            m.reactions.NGAM.lower_bound = 0  # maintenance is a forced ATP drain; the test asks about free drains
        for b in m.boundary:
            b.lower_bound = 0
        for name, st in tests.items():
            with m:
                t = cobra.Reaction("EGC_test")
                m.add_reactions([t])
                t.add_metabolites({M.get_by_id(k): v for k, v in st.items()})
                assert not t.check_mass_balance(), (name, t.check_mass_balance())
                m.objective = t
                out[name] = round(m.slim_optimize(error_value=0.0), 6)
    return out


def drop_starch(m):
    b = m.reactions.get_by_id(BIO)
    s = m.metabolites.get_by_id("starch300_h")
    b.add_metabolites({s: -b.metabolites[s]})


if __name__ == "__main__":
    orig = cobra.io.read_sbml_model(SRC)
    orig.solver = "glpk"
    m = orig.copy()
    closed, added = patch(m)
    print(f"E1 closed {closed} gene-less duplicates; E2 gave genes to {added} gene-less reactions")
    for label, model in [("original", orig), ("patched", m)]:
        g = model.slim_optimize(error_value=0.0)
        with model:
            model.genes.get_by_id(PGM_GENE).knock_out()
            drop_starch(model)
            gk = model.slim_optimize(error_value=0.0)
        egc = energy_cycle_check(model)
        bad = {k: v for k, v in egc.items() if v > 1e-6}
        print(f"{label:9s} growth {g:.4f} | PGM KO + no starch {gk:.4f} | energy drains tested {len(egc)}, "
              f"carrying flux {len(bad)} {bad}")
    egc = energy_cycle_check(m)
    assert all(v <= 1e-6 for v in egc.values()), egc
    # light-limited sanity check (catches cycles the drain test cannot: free energy from a transport loop)
    with m:
        m.reactions.NGAM.lower_bound = 0
        for r in m.reactions:
            if r.id.startswith("EX_photonVis") and r.lower_bound < 0:
                r.lower_bound = -30.0
        sol = m.optimize()
        b = m.reactions.get_by_id(BIO)
        gper = sum(-c * x.formula_weight / 1000 for x, c in b.metabolites.items() if c < 0 and x.id not in ("atp_c", "h2o_c"))
        yld = gper * sol.objective_value / 0.03
        ci = -sol.fluxes.get("EX_co2(e)", 0) - sol.fluxes.get("EX_hco3(e)", 0)
        nh4 = sol.fluxes.get("EX_nh4(e)", 0)
        with m:
            m.reactions.get_by_id("EX_nh4(e)").upper_bound = 0  # no ammonium export: the yield must not depend on it
            yld_nx = gper * m.slim_optimize() / 0.03
        print(f"light-limited check at 30 mmol photons: yield {yld:.2f} g/mol ({yld_nx:.2f} with NH4 export blocked), net Ci uptake {ci:.2f}")
        assert ci > 0 and abs(yld - yld_nx) < 0.03 * yld and yld < 2.6, (yld, yld_nx, ci, nh4)
    pd.DataFrame(changes).to_csv(LOG, sep="\t", index=False)
    cobra.io.write_sbml_model(m, OUT)
    print("wrote", OUT, "and", LOG)
