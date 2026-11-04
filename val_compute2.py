"""FVA of normalized fluxes and minimum quantum requirement for Fig. 3. Adds to val_cache.json."""
import json, warnings, time
import cobra
warnings.filterwarnings("ignore")
exec(open("val_compute.py").read().split("m = cobra.io.read_sbml_model")[0])   # FLUXES, constants
m = cobra.io.read_sbml_model("iPiccel954.xml")
m.solver = "glpk"
R = json.load(open("val_cache.json"))
m.reactions.get_by_id("EX_nh4(e)").upper_bound = 0   # urea is the nitrogen source, no ammonium excretion
coef = R["coef"]
bio = m.reactions.get_by_id("Biomass_Picre_Auto")
ph = m.reactions.get_by_id("EX_photonVis(e)")
co2, hco3 = m.reactions.get_by_id("EX_co2(e)"), m.reactions.get_by_id("EX_hco3(e)")
Uexpr = -(co2.flux_expression + hco3.flux_expression)


DEN = [Uexpr]


def feasible(expr, rho, side):
    """side 'max': is there a solution with expr >= rho*U ; side 'min': with expr <= rho*U."""
    m.objective = m.problem.Objective(expr - rho * DEN[0], direction="max" if side == "max" else "min")
    s = m.optimize()
    if s.status != "optimal":
        return False
    return s.objective_value > 1e-7 if side == "max" else s.objective_value < -1e-7


def ratio_range(expr, hi=30.0):
    # largest ratio
    if feasible(expr, hi, "max"):
        rmax = hi
    else:
        lo, up = 0.0, hi
        for _ in range(22):
            mid = (lo + up) / 2
            lo, up = (mid, up) if feasible(expr, mid, "max") else (lo, mid)
        rmax = lo
    # smallest ratio (flux may be reversible, so allow negative)
    lo, up = -hi, hi
    if not feasible(expr, -hi, "min"):
        for _ in range(22):
            mid = (lo + up) / 2
            lo, up = (lo, mid) if feasible(expr, mid, "min") else (mid, up)
        rmin = up
    else:
        rmin = -hi
    return 100 * rmin, 100 * rmax


rng = {}
t0 = time.time()
for opt_tag, g_lo, g_hi in [("at_0.323", 0.323, 0.323)]:
    with m:
        bio.lower_bound, bio.upper_bound = g_lo, g_hi
        for k, (name, v, d) in FLUXES.items():
            expr = sum(c * m.reactions.get_by_id(i).flux_expression for i, c in d.items())
            rng[k] = ratio_range(expr)
            print(k, v, [round(x, 1) for x in rng[k]], round(time.time() - t0))
R["fva"] = rng


# the same ranges per unit of RuBisCO carboxylation (removes the dependence on how net CO2 uptake is defined)
rng2 = {}
with m:
    bio.lower_bound, bio.upper_bound = 0.323, 0.323
    DEN[0] = m.reactions.get_by_id("RBPCh").flux_expression
    for k, (name, v, d) in FLUXES.items():
        if k == "rubisco_c":
            continue
        expr = sum(c * m.reactions.get_by_id(i).flux_expression for i, c in d.items())
        rng2[k] = ratio_range(expr)
        print("per rubisco", k, round(100 * v / 120, 1), [round(x, 1) for x in rng2[k]])
R["fva_rub"] = rng2
DEN[0] = Uexpr
# growth with unlimited light and the growth-light curve up to saturation
curve = []
for P in [10, 25, 50, 100, 200, 400, 700, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 15000]:
    with m:
        ph.lower_bound = -P / coef
        mu = m.slim_optimize()
        sol = m.optimize()
        g = sol.fluxes
        used = 2 * g["PSIred"] + 4 * g["PSIIred"] + 2 * g["PSIblue"] + 4 * g["PSIIblue"]
    curve.append([P, mu, used])
R["curve2"] = curve
with m:
    ph.lower_bound = -1e5
    R["mu_unlimited"] = m.slim_optimize()
print("unlimited light", R["mu_unlimited"])

# minimum photons per carbon, for starch and for biomass (no maintenance of other products)
with m:
    sd = m.add_boundary(m.metabolites.get_by_id("starch300_h"), type="demand")
    bio.lower_bound = 0
    P = 100.0
    ph.lower_bound = -P / coef
    m.objective = sd
    fs = m.slim_optimize()
    Cstarch = fs * 300 * 6
    R["qr_starch_abs"] = P / Cstarch
    s = m.optimize()
    g = s.fluxes
    used = 2 * g["PSIred"] + 4 * g["PSIIred"] + 2 * g["PSIblue"] + 4 * g["PSIIblue"]
    R["qr_starch_used"] = used / Cstarch
    print("starch flux", fs, "QR absorbed", R["qr_starch_abs"], "QR used", R["qr_starch_used"])
json.dump(R, open("val_cache.json", "w"))