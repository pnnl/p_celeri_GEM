# p_celeri_GEM

Genome-scale metabolic model and enzyme-constrained models of *Picochlorum celeri*, with the data and code for the manuscript "Enzyme-Constrained Metabolic Modeling Reveals CO₂-Dependent Growth Control in *Picochlorum celeri*".

The repository contains iPiccel954, a genome-scale model with 3525 reactions, 2973 metabolites and 954 genes, and eciPiccel954, its enzyme-constrained version. The enzyme-constrained models are built for six conditions, wild type in air, wild type at elevated CO₂, and four carbonic anhydrase lines at elevated CO₂ (CAII, CA3, CA4 and CA7).

## Files

| File | Description |
| --- | --- |
| iPiccel954.xml | The final genome-scale model in SBML format |
| iPicre1230.xml | The *Picochlorum renovo* template model used for reconstruction |
| iPceleri.fasta | Protein sequences of *P. celeri* |
| MEMOTE.html | MEMOTE quality report of the model |
| ecGEMs.zip | Enzyme-constrained models, with one folder per condition and the kcat relaxation algorithm for each |
| protein_abundance_g_per_gCDW.xlsx | Measured protein abundance for each condition |
| leave_out.py and leave_out_results.xlsx | Leave-out test of the kcat calibration and its results |
| val_compute.py to val_compute4.py and val_cache.json | Model validation calculations and their cached results |
| figure_1.py to figure_5.py | Scripts for the main figures |
| figure_S1.py to figure_S5.py | Scripts for the supplementary figures |
| LICENSE.txt and DISCLAIMER.txt | License and disclaimer |

## Requirements

Python 3 with the following packages.

```
pip install cobra pandas openpyxl optlang matplotlib gurobipy
```

The scripts use Gurobi by default, which needs a valid license. Any other solver supported by COBRApy can be used instead.

## How to use

To load the model and run flux balance analysis, use the following.

```python
from cobra import io

model = io.read_sbml_model("iPiccel954.xml")
print(model.optimize().objective_value)
```

To work with the enzyme-constrained models, unzip ecGEMs.zip. Each condition folder contains ecGEM.py, which builds and solves the enzyme-constrained model, and a relaxation_algorithm folder with the kcat relaxation script. Run each script from inside the folder it sits in.

```
python ecGEM.py
```

Figure scripts read the model, the protein abundance file and the unzipped ecGEMs folder, so keep these together in the same directory when running them.

```
python figure_1.py
```

## Data sources

Enzyme turnover numbers were predicted with TurNuP. Protein abundances come from Torres-Tiji et al., 2026.

## Contact

Niaz Bahar Chowdhury, Environmental Molecular Sciences Laboratory, Pacific Northwest National Laboratory
