"""Net inorganic carbon uptake range (FVA) at the measured TG2 growth rate. Adds 'cuptake' to val_cache.json."""
import json, warnings, cobra
warnings.filterwarnings("ignore")
m = cobra.io.read_sbml_model("iPiccel954.xml"); m.solver = "glpk"
m.reactions.get_by_id("EX_nh4(e)").upper_bound = 0
MU = 0.323
m.reactions.get_by_id("Biomass_Picre_Auto").lower_bound = MU
m.reactions.get_by_id("Biomass_Picre_Auto").upper_bound = MU
expr = -(m.reactions.get_by_id("EX_co2(e)").flexpr if False else m.reactions.get_by_id("EX_co2(e)").flux_expression + m.reactions.get_by_id("EX_hco3(e)").flux_expression)
out = {}
for d in ("min", "max"):
    m.objective = m.problem.Objective(expr, direction=d)
    out[d] = m.slim_optimize()
R = json.load(open("val_cache.json"))
R["cuptake"] = {"min": out["min"], "max": out["max"], "pfba": R["U"], "meas_TG2": 9.59, "meas_TG1": 5.96}
json.dump(R, open("val_cache.json", "w"))
print(R["cuptake"])