# Analysis Plan

## Objective

The objective of this exploratory analysis was to evaluate whether CSF proteomic biomarkers improve fast-versus-slow ALS progression classification compared with duration-free clinical covariates and nested data-driven protein selection.

## Cohort

The binary modeling cohort included baseline ALS patients classified as fast or slow progressors. Intermediate progressors were excluded from binary classification.

## Outcome

The outcome was fast-versus-slow ALS progression status, defined using CPR/rate of progression.

## Primary Analysis

The primary biomarker model was the full 14-protein literature-guided CSF panel.

## Secondary Analysis

The secondary biomarker model was the final mechanism 5-10 panel.

## Comparator Models

Comparator models included:

1. A duration-free clinical model using age, sex, and site of onset.
2. A literature-guided protein plus clinical no-duration model.
3. Nested data-driven protein-only models using top-5, top-7, and top-10 proteins selected inside each training fold.

## Leakage-Safety Rules

Disease duration was excluded from all primary models because the outcome depends on time from onset.

ALSFRS-R was excluded from all predictors because the progression outcome depends on functional decline.

Data-driven protein ranking and feature selection were performed inside the outer training fold only.

Signed-score direction estimation, imputation, scaling, and model fitting were performed using training-fold data only.

## Model Evaluation

Models were evaluated using repeated stratified 5-fold cross-validation with 100 repeats.

Performance was summarized using held-out predictions from the outer test folds.

Primary metrics included ROC AUC, balanced accuracy, sensitivity for fast progressors, specificity for slow progressors, precision, and F1 score.

## Interpretation

Because the cohort was small and internally cross-validated, results should be interpreted as exploratory and hypothesis-generating. External validation is required before clinical use.
