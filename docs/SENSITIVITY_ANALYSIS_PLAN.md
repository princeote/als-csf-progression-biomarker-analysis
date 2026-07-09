# Sensitivity Analysis Plan

## Purpose

This file documents the planned sensitivity analyses for the ALS CSF fast-vs-slow progression classifier. The primary manuscript model prioritizes the predefined final mechanism 7 panel, while sensitivity analyses test whether results depend on class weighting, threshold choice, or exploratory data-driven protein selection.

## Class Weighting Sensitivity

The binary ALS modeling cohort is imbalanced, with 10 fast progressors and 22 slow progressors. Primary machine-learning models therefore use `class_weight="balanced"` to reduce bias toward the majority slow-progressor class.

The revised pipeline evaluates four model variants for logistic regression and random forest:

| Algorithm label | Model type | Class weighting role |
|---|---|---|
| `logistic_regression` | L2-regularized logistic regression | Balanced primary model |
| `logistic_regression_unweighted` | L2-regularized logistic regression | Unweighted sensitivity model |
| `random_forest` | Constrained random forest | Balanced primary model |
| `random_forest_unweighted` | Constrained random forest | Unweighted sensitivity model |

Balanced and unweighted models are evaluated using the same cross-validation splits, features, preprocessing, and metrics. This allows direct assessment of whether performance estimates depend on the class-weighting strategy.

## No SMOTE or Synthetic Oversampling

SMOTE and other synthetic oversampling methods are not used. Because the fast-progressor class contains only 10 patients, synthetic sample generation could exaggerate apparent separability and increase overfitting risk. Class weighting is preferred because it changes the training loss without creating artificial patients.

## Data-Driven Protein Panel Sensitivity

The revised pipeline limits data-driven feature-selection models to a single exploratory top-10 protein panel. Feature ranking is performed inside each outer training fold only, using the combined rank of training-fold Mann-Whitney p-value and absolute Cohen's d, after excluding proteins marked as age- or sex-associated in healthy controls.

The top-10 data-driven panel is treated as an exploratory comparator, not as the preferred manuscript panel. The preferred manuscript panel remains the predefined final mechanism 7 panel.

## Youden's J Threshold Sensitivity

The primary analysis uses a fixed classification threshold of 0.50. As a threshold sensitivity analysis, the pipeline can tune thresholds inside the training data using Youden's J statistic:

\[
J = sensitivity + specificity - 1
\]

Run the threshold sensitivity analysis with:

```bash
python code/ALS_nested_model_comparison_pipeline.py --data-dir data --output-dir outputs_youden --n-repeats 100 --tune-thresholds
```

AUC is not changed by threshold tuning. Youden's J affects threshold-dependent metrics such as sensitivity, specificity, balanced accuracy, precision, F1 score, and confusion-matrix counts.

## Manuscript Interpretation

Balanced models should be treated as primary for the imbalanced cohort. Unweighted models, Youden's J threshold tuning, and the data-driven top-10 panel should be described as sensitivity or exploratory comparator analyses. None of these analyses should be interpreted as external clinical validation.
