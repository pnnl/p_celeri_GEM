"""Model calculations for Fig. 3 (validation). Writes val_cache.json."""
import json, warnings
import cobra
from cobra.flux_analysis import pfba
warnings.filterwarnings("ignore")
MODEL_PATH = "iPiccel952.xml"
GROWTH_TG2 = 0.323
# measured TG2 INST-MFA fluxes (Steichen et al. 2024, Fig. 3), per 100 units of net CO2 uptake, and model reactions
FLUXES = {
    "rubisco_c": ("RuBisCO carboxylation", 120, {"RBPCh": 1}),
    "rubisco_o": ("RuBisCO oxygenation", 0.4, {"RBCh": 1, "NanoG0460": 1}),
    "pga_red": ("3PGA reduction", 215, {"NanoG0437": 1}),
    "sbp": ("SBPase", 42, {"SBP": 1}),
    "starch": ("Starch synthesis", 1.6, {"GLPThi": 1}),
    "pdh_h": ("Chloroplast pyruvate dehydrogenase", 13, {"NanoG0420": 1}),
    "enol": ("Cytosolic enolase", 30, {"NanoG0405": 1}),
    "pepc": ("PEP carboxylase", 10, {"PPC": 1, "NanoG1707": 1}),
    "pdh_m": ("Mitochondrial pyruvate dehydrogenase", 4, {"NanoG0421": 1}),
    "cs_m": ("Citrate synthase", 3.5, {"CSm": 1}),
    "me_m": ("Malic enzyme", 4, {"MDHCm": 1}),
}
m = cobra.io.read_sbml_model(MODEL_PATH)
m.solver = "glpk"
ph = m.reactions.get_by_id("EX_photonVis(e)")
coef = sum(c for met, c in m.reactions.get_by_id("PRISM_design_growth").metabolites.items()
           if met.id.startswith("photon") and c > 0)
bio = m.reactions.get_by_id("Biomass_Picre_Auto")
GAM = {"atp_c", "h2o_c", "adp_c", "pi_c", "h_c"}
mass = cmass = 0.0
for met, c in bio.metabolites.items():
    if c < 0 and met.id not in GAM and met.formula:
        mass += -c * met.formula_weight
        cmass += -c * met.elements.get("C", 0) * 12.011
out = {"coef": coef, "cfrac": cmass / mass, "mu_max": m.slim_optimize()}
print("C fraction", out["cfrac"], "max growth", out["mu_max"], "photons per photonVis", coef)

# flux state at the measured TG2 growth rate, with the least light, then the least total flux
with m:
    bio.lower_bound = GROWTH_TG2
    m.objective = m.problem.Objective(ph.flux_expression, direction="max")
    s1 = m.optimize()
    ph.lower_bound = s1.fluxes["EX_photonVis(e)"] - 1e-6
    m.objective = bio
    f = pfba(m).fluxes
U = -(f["EX_co2(e)"] + f["EX_hco3(e)"])
out["U"] = U
out["flux_model"] = {k: 100 * sum(c * f[i] for i, c in d.items()) / U for k, (n, v, d) in FLUXES.items()}
out["flux_meas"] = {k: v for k, (n, v, d) in FLUXES.items()}
print("U", U)
for k in FLUXES:
    print(k, FLUXES[k][1], round(out["flux_model"][k], 1))

# growth against absorbed photon flux
curve = []
for P in [10, 25, 50, 100, 200, 300, 450, 600, 800, 1000, 1200, 1400, 1550, 1700, 2000]:
    with m:
        ph.lower_bound = -P / coef
        mu = m.slim_optimize()
        sol = m.optimize()
        g = sol.fluxes
        used = 2 * g["PSIred"] + 4 * g["PSIIred"] + 2 * g["PSIblue"] + 4 * g["PSIIblue"]
        cup = -(g["EX_co2(e)"] + g["EX_hco3(e)"])
    curve.append([P, mu, used, cup])
    print(P, round(mu, 4), round(used, 1), round(cup, 2))
out["curve"] = curve
json.dump(out, open("val_cache.json", "w"))
