# iPiccel954 patched (copy of iPiccel954 curated against the lab's measurements)

This branch adds one folder and changes nothing else in the repository. `iPiccel954.xml` and all of your scripts, figures and tables on `main` are untouched.

Read first: `iPiccel954_changes_for_Niaz.docx` (one page: what changed, why, and what it means for the draft). Full per-reaction log and every comparison with the draft: `iPiccel954_patched_report.docx`.

| File | What it is |
|---|---|
| `iPiccel954_patched.xml` | patched stoichiometric model (3,641 reactions, 3,040 metabolites, 1,005 genes) |
| `irrev_iPiccel954_patched_Wt_CO2.xml`, `irrev_iPiccel954_patched_WT_Air.xml` | enzyme-constrained models rebuilt with your relaxation algorithm unchanged, on the lab's measured proteome |
| `patch_iPiccel954.py` | every change to the model, step by step (A to N), each logged per reaction |
| `patch_ecGEMs.py` | the ecGEM rebuild |
| `predict_enzyme_kcats.py`, `kcat_predictions.tsv` | DeepEnzyme kcat predictions for genes that carried the median kcat (model estimates, not assays) |

The scripts are included for reference and read inputs from the lab workspace (`/home/yasin/Claudio/GEM`), so they will not run unchanged elsewhere.
