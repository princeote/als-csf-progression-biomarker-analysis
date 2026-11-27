# Candidate CSF protein panels associated with ALS progression rate

Companion code, results, and figures for *Candidate Cerebrospinal Fluid Protein Panels Associated with ALS Progression Rate: A Proteomic and Machine Learning Analysis* by Prince Otegbulu, California Academy of Mathematics and Science.

This exploratory study screened 793 cerebrospinal fluid (CSF) protein features in 10 fast and 22 slow ALS progressors. It ranked 22 screen-positive proteins with a random forest and evaluated three panels that included age and sex as covariates. The continuous clinical progression rate (CPR) analysis used all 43 ALS participants, including the 11 intermediate progressors.

| Panel | Mean ROC AUC ± SD | Fast sensitivity | Slow specificity |
| --- | ---: | ---: | ---: |
| Top-5 | 0.876 ± 0.029 | 70.0% | 90.9% |
| Top-7 | 0.835 ± 0.045 | 60.0% | 81.8% |
| Top-10 | 0.849 ± 0.038 | 70.0% | 72.7% |

AUC summarizes ten repeats of stratified five-fold cross-validation. Sensitivity and specificity use each participant's vote across those repeats, with classification cutoffs selected inside each outer training split. The protein screen and panel ranking used the full Fast–Slow cohort before cross-validation, so the AUCs are likely optimistic. The progression groups describe prior decline at the sampled visit. The panel still requires independent validation for future decline. No protein reached a 5% false discovery rate in the 793-feature screen.

## Repository contents

| Location | Contents |
| --- | --- |
| [paper/](paper/) | LaTeX manuscript, bibliography, and seven figures. Figures 4 and 6 use aggregate statistics in place of the individual measurements shown in the supplied manuscript. |
| [reproduce/](reproduce/) | Cohort preparation, screening, random forest ranking and cross-validation, and generation of the two aggregate figures. |
| [results/](results/) | Aggregate screen and CPR statistics, Gini rankings, panel membership, performance summaries, per-repeat AUCs, and analysis settings. |

Start with the [panel performance table](results/PANEL_PERFORMANCE_10x5.csv), the [Top-10 protein evidence table](results/FINAL_TOP10_EVIDENCE.csv), or the [full protein screen](results/PROTEIN_SCREENING_AUTHENTIC_793.csv). In the result tables, positive Cohen's d means higher estimated abundance in the Fast group. The CPR correlations are unadjusted Spearman correlations. Their BH q-values correct across the 22 selected proteins.

## Build the paper and figures

The analysis and figure scripts were checked with Python 3.14 and the versions in [requirements.txt](requirements.txt). From the repository root, activate a Python 3.14 virtual environment, then run:

```sh
python -m pip install -r requirements.txt
python reproduce/make_public_figures.py
```

The figure command rebuilds Figures 4 and 6 from the included aggregate CSVs. The other five PNGs are supplied manuscript assets. To compile the paper with Tectonic, run:

```sh
cd paper
tectonic main.tex
```

## Rerun the random forest analysis

The source study's raw mass-spectrometry files are available in [PRIDE PXD024219](https://www.ebi.ac.uk/pride/archive/projects/PXD024219). This analysis uses the matching processed intensity, clinical, and protein annotation tables supplied by the [source study team](https://doi.org/10.3389/fnins.2021.642324). Those inputs and participant-level outputs are restricted under the data use agreement. Contact the source study team for access to the processed exports.

With authorized access, put the matching `clinData2.csv`, `exprData.csv`, and `featData.csv` files in a private directory outside this checkout. The loader verifies that they match the expected exports. Run the following from the repository root using two separate directories outside the checkout:

```sh
python reproduce/rf_age_sex_adjusted_10x5.py --data-dir "/path/to/approved-exports" --output-dir "/path/to/private-rerun"
```

The rerun writes participant-level predictions and other restricted outputs to the private output directory. Keep the full rerun directory outside version control. The public tables and figures can be inspected using the files in this repository.
