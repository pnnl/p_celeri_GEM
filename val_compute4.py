"""Extra validation numbers. Adds 'extra' to val_cache.json.
  net photosynthetic quotient and absorbed photons per net O2 evolved at limiting light (model, min photon then pFBA)
  growth of in silico lutein-free and chlorophyll b-free biomass at the default photon bound
"""
import json, warnings, cobra
from cobra.flux_analysis import pfba
warnings.filterwarnings("ignore")
m = cobra.io.read_sbml_model("iPiccel952.xml"); m.solver = "glpk"
m.reactions.get_by_id("EX_nh4(e)").upper_bound = 0
R = json.load(open("val_cache.json")); coef = R["coef"]
bio = m.reactions.get_by_id("Biomass_Picre_Auto"); ph = m.reactions.get_by_id("EX_photonVis(e)")
with m:
    bio.lower_bound = 0.1
    m.objective = m.problem.Objective(ph.flux_expression, direction="max"); s = m.optimize()
    ph.lower_bound = s.fluxes["EX_photonVis(e)"] * 1.0001
    f = pfba(m).fluxes
absorbed = -f["EX_photonVis(e)"] * coef
net_o2 = f["DM_o2D(u)"] + f["EX_o2(e)"]
cu = -(f["EX_co2(e)"] + f["EX_hco3(e)"])
ex = {"photons_per_net_O2": absorbed / net_o2, "PQ_net": net_o2 / cu, "PQ_gross": f["DM_o2D(u)"] / cu}
# in silico pigment-free biomass
wt = m.slim_optimize()
for key, mid in (("lutein_free", "lut_u"), ("chlb_free", "chlb_u")):
    with m:
        met = m.metabolites.get_by_id(mid); c = bio.metabolites[met]
        bio.add_metabolites({met: -c})
        ex[key] = {"coef_removed": c, "mu": m.slim_optimize(), "mu_wt": wt}
R["extra"] = ex; json.dump(R, open("val_cache.json", "w")); print(ex)
