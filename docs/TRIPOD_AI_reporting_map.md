# TRIPOD+AI-Style Reporting Map for This Manuscript

This is a manuscript-preparation map aligned to TRIPOD+AI reporting expectations. It is not a replacement for the official TRIPOD+AI checklist; use the official checklist during final submission.

## Title and abstract

- State that this is an exploratory internally cross-validated prediction-model study.
- Identify the target outcome: fast versus slow ALS progression.
- Identify the data source: baseline CSF proteomic and clinical covariate data from PRIDE project PXD024219.
- State the modeling approach: regularized logistic regression, constrained random forest, fixed literature-guided protein panels, and one nested data-driven top-10 feature-selection comparator.
- Avoid claiming clinical validation.

## Introduction

- Explain why ALS progression stratification matters.
- Explain why disease duration and ALSFRS-R are excluded from predictors because CPR/rate of progression is derived from ALSFRS-R and time from onset.
- Justify the use of a predefined mechanism panel in a small cohort where unconstrained discovery is unstable.

## Methods: data source and participants

- Report the baseline cohort size: ALS n=43 and healthy controls n=25.
- Report the binary modeling cohort: n=32, with 10 fast and 22 slow progressors.
- Report that 11 intermediate progressors were excluded from binary modeling.
- Define baseline-visit selection and any sample-level exclusions.

## Methods: outcome

- Define CPR as `(48 - ALSFRS-R) / months from symptom onset`.
- Define fast progressors as CPR >= 1.0 and slow progressors as CPR <= 0.5.
- Explicitly state that CPR depends on ALSFRS-R and time from onset.
- State that intermediate progressors were excluded from binary classification.

## Methods: predictors

- Clinical no-duration predictors: age at visit, sex, bulbar onset, and limb onset.
- Primary predefined mechanism panel: CHIT1, RBP4, AMBP, SERPINA4, CST3, TNR, and NCAN.
- Broader literature-supported comparator panel: the 14 assay-available literature-supported proteins.
- State that disease duration and ALSFRS-R were not used as primary predictors.

## Methods: missing data and preprocessing

- Report log2(x + 1) transformation and sample median centering.
- Report median imputation inside sklearn pipelines.
- Report standardization for logistic-regression pipelines.
- Do not claim batch-effect removal unless it occurred upstream and is separately documented.

## Methods: model development and validation

- Report repeated stratified 5-fold cross-validation with 100 repeats.
- Report balanced and unweighted versions of L2-regularized logistic regression and constrained random forest.
- State that balanced models are primary for the imbalanced cohort and unweighted models are sensitivity analyses.
- State that thresholds were fixed at 0.50 unless Youden's J threshold sensitivity is explicitly enabled.

## Methods: feature selection

- For fixed panels, describe literature/mechanism criteria and assay availability.
- For data-driven models, describe training-fold Mann-Whitney p-value, absolute Cohen's d, combined rank, and exclusion of healthy-control age/sex confounders.
- Emphasize that no full-dataset protein ranking was used for held-out model evaluation.

## Methods: performance metrics

- Primary metric: held-out ROC AUC.
- Secondary metrics: balanced accuracy, sensitivity for fast progressors, specificity for slow progressors, precision, and F1.
- Report repeat-level 95% intervals as internal resampling intervals across repeated cross-validation, not as external-validation confidence intervals.

## Results

- Report sample sizes and class balance.
- Present model comparisons using the cleaned model set only.
- Compare clinical-only, final mechanism 7 protein-only, final mechanism 7 plus clinical, broader 14-protein comparator, and data-driven top-10 models.
- Interpret any data-driven result as exploratory because feature selection is fold-specific and small-sample unstable.

## Discussion

- Interpret the findings as exploratory evidence that literature-guided CSF proteins may contain progression-relevant signal.
- Do not claim clinical readiness.
- Discuss whether adding duration-free clinical covariates improves or weakens performance in this cohort.
- Discuss the NCAN assay limitation: the assay may capture total/intact NCAN abundance rather than a disease-linked cleavage fragment.

## Limitations

- Small binary cohort: n=32 with only 10 fast progressors.
- No external validation cohort.
- Repeated cross-validation reuses the same patients across repeats.
- Multiple model variants were compared.
- Protein measurements depend on assay coverage and feature matching.
- Findings may be sample-size-driven or overfit.

## Transparency and reproducibility

- Provide the analysis code.
- Provide derived output tables and figure data after rerunning the cleaned pipeline.
- State whether raw patient-level data can be shared; if not, explain restrictions.
- Include software dependencies and run commands.
