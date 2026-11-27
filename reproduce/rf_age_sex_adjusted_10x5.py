#!/usr/bin/env python3
"""RF-only ALS progression analysis with age/sex covariate adjustment.

This is the revised analysis requested for the manuscript.  The protein
screen remains the two-part Fast-versus-Slow screen (MWU
p < 0.05 and |Cohen's d| >= 0.8), but the healthy-control significance
filter is not used to remove proteins.  All screen-positive proteins are
available for the age/sex-inclusive Gini ranking, with membership determined by protein importances.  Age and sex are then appended
as adjustment covariates to each Top-5, Top-7, and Top-10 Random Forest fit;
they are not counted as members of a protein panel.  Disease duration is not
used.

The repeated-CV performance calculation uses training-only inner-CV Youden
thresholds.  Protein ranking is a development-stage full-cohort lock, so the
results are internal resampling estimates conditional on that lock rather
than an unbiased estimate of a complete discovery pipeline.  Participant-level
majority-vote ties are resolved by the mean probability-minus-threshold margin.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

from als_data_and_screen import load_data, protein_screen, bh_adjust


SEED = 42
FOLDS = 5
REPEATS = 10
PANEL_SIZES = (5, 7, 10)
ADJUSTMENT_COLUMNS = ("age_visit_years", "sex_male")
RF_PARAMS: dict[str, Any] = {
    "n_estimators": 100,
    "max_depth": 5,
    "min_samples_split": 5,
    "min_samples_leaf": 2,
    "class_weight": {0: 1.0, 1: 3.0},
    "random_state": SEED,
    "n_jobs": 1,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def model_matrix(frame: pd.DataFrame, proteins: list[str]) -> pd.DataFrame:
    """Return protein panel plus age/sex adjustment features only."""
    return frame[proteins + list(ADJUSTMENT_COLUMNS)].apply(pd.to_numeric, errors="coerce")


def new_model(seed: int) -> RandomForestClassifier:
    params = dict(RF_PARAMS)
    params["random_state"] = seed
    return RandomForestClassifier(**params)


def fit_predict(model_x: pd.DataFrame, y: np.ndarray, train: np.ndarray, test: np.ndarray, seed: int) -> np.ndarray:
    imputer = SimpleImputer(strategy="median")
    x_train = imputer.fit_transform(model_x.iloc[train])
    x_test = imputer.transform(model_x.iloc[test])
    model = new_model(seed).fit(x_train, y[train])
    return model.predict_proba(x_test)[:, 1]


def choose_youden(y: np.ndarray, probability: np.ndarray) -> float:
    fpr, tpr, thresholds = roc_curve(y, probability)
    j = tpr - fpr
    return float(thresholds[int(np.nanargmax(j))])


def metrics(y: np.ndarray, probability: np.ndarray, prediction: np.ndarray) -> dict[str, Any]:
    tn, fp, fn, tp = confusion_matrix(y, prediction, labels=[0, 1]).ravel()
    return {
        "auc": float(roc_auc_score(y, probability)),
        "accuracy": float(accuracy_score(y, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(y, prediction)),
        "fast_sensitivity": float(tp / (tp + fn)) if tp + fn else float("nan"),
        "slow_specificity": float(tn / (tn + fp)) if tn + fp else float("nan"),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def inner_threshold(model_x: pd.DataFrame, y: np.ndarray, train: np.ndarray, split_index: int) -> tuple[float, dict[str, Any]]:
    """Choose Youden's cutoff from inner-CV predictions within outer training."""
    inner = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=SEED + split_index)
    x_train = model_x.iloc[train].reset_index(drop=True)
    y_train = y[train]
    probability = np.zeros(len(y_train), dtype=float)
    for inner_index, (inner_train, inner_test) in enumerate(inner.split(x_train, y_train), start=1):
        probability[inner_test] = fit_predict(x_train, y_train, inner_train, inner_test, SEED + split_index * 100 + inner_index)
    threshold = choose_youden(y_train, probability)
    internal = metrics(y_train, probability, (probability >= threshold).astype(int))
    internal["youden_j"] = internal["fast_sensitivity"] + internal["slow_specificity"] - 1.0
    return threshold, internal


def protein_ranking(binary: pd.DataFrame, candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    proteins = candidates["Protein"].tolist()
    adjusted_features = proteins + list(ADJUSTMENT_COLUMNS)
    imputer = SimpleImputer(strategy="median")
    x = imputer.fit_transform(model_matrix(binary, proteins))
    model = new_model(SEED).fit(x, binary["y"].astype(int).to_numpy())
    ranking = candidates.copy()
    importance = pd.Series(model.feature_importances_, index=adjusted_features)
    ranking["RF_Gini_Importance"] = ranking["Protein"].map(importance)
    protein_total = float(ranking["RF_Gini_Importance"].sum())
    ranking["RF_Gini_Share"] = ranking["RF_Gini_Importance"] / protein_total if protein_total else 0.0
    covariate_importance = pd.DataFrame({
        "Feature": list(ADJUSTMENT_COLUMNS),
        "RF_Gini_Importance": [float(importance[column]) for column in ADJUSTMENT_COLUMNS],
        "RF_Gini_Share_All_Features": [float(importance[column]) for column in ADJUSTMENT_COLUMNS],
        "Role": "adjustment_covariate",
    })
    ranking = ranking.sort_values(
        ["RF_Gini_Importance", "Abs_Cohens_d", "Protein"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    ranking.insert(0, "RF_Gini_Rank", np.arange(1, len(ranking) + 1))
    ranking["Top5"] = ranking["RF_Gini_Rank"] <= 5
    ranking["Top7"] = ranking["RF_Gini_Rank"] <= 7
    ranking["Top10"] = ranking["RF_Gini_Rank"] <= 10
    return ranking, covariate_importance


def repeated_panel(binary: pd.DataFrame, proteins: list[str], panel_size: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = binary["y"].astype(int).to_numpy()
    x = model_matrix(binary, proteins)
    cv = RepeatedStratifiedKFold(n_splits=FOLDS, n_repeats=REPEATS, random_state=SEED)
    rows: list[dict[str, Any]] = []
    threshold_rows: list[dict[str, Any]] = []
    for split_index, (train, test) in enumerate(cv.split(np.zeros(len(y)), y), start=1):
        repeat = (split_index - 1) // FOLDS + 1
        fold = (split_index - 1) % FOLDS + 1
        threshold, internal = inner_threshold(x, y, train, split_index)
        probability = fit_predict(x, y, train, test, SEED + split_index)
        threshold_rows.append({
            "panel_size": panel_size,
            "repeat": repeat,
            "fold": fold,
            "split": split_index,
            "threshold": threshold,
            "inner_auc": internal["auc"],
            "inner_fast_sensitivity": internal["fast_sensitivity"],
            "inner_slow_specificity": internal["slow_specificity"],
            "inner_balanced_accuracy": internal["balanced_accuracy"],
            "inner_youden_j": internal["youden_j"],
        })
        for local_index, row_index in enumerate(test):
            row = binary.iloc[row_index]
            prediction = int(probability[local_index] >= threshold)
            rows.append({
                "panel_size": panel_size,
                "repeat": repeat,
                "fold": fold,
                "split": split_index,
                "SubjID": row["SubjID"],
                "Group": row["Group"],
                "y_true": int(y[row_index]),
                "probability_fast": float(probability[local_index]),
                "y_probability_fast": float(probability[local_index]),
                "training_margin": float(probability[local_index] - threshold),
                "training_only_youden_threshold": float(threshold),
                "threshold": float(threshold),
                "training_only_prediction": prediction,
                "y_pred": prediction,
                "selected_features": ";".join(proteins),
            })
    predictions = pd.DataFrame(rows)
    thresholds = pd.DataFrame(threshold_rows)
    repeat_rows = []
    for repeat, sub in predictions.groupby("repeat", sort=True):
        repeat_rows.append({
            "panel_size": panel_size,
            "repeat": int(repeat),
            "auc": float(roc_auc_score(sub["y_true"], sub["probability_fast"])),
            "n_predictions": int(len(sub)),
            "n_fast": int(sub["y_true"].sum()),
            "n_slow": int((sub["y_true"] == 0).sum()),
        })
    repeat_auc = pd.DataFrame(repeat_rows)
    patient = predictions.groupby("SubjID", as_index=False).agg(
        panel_size=("panel_size", "first"),
        Group=("Group", "first"),
        y_true=("y_true", "first"),
        probability_fast=("probability_fast", "mean"),
        probability_fast_sd=("probability_fast", lambda s: float(s.std(ddof=1))),
        probability_fast_min=("probability_fast", "min"),
        probability_fast_max=("probability_fast", "max"),
        mean_training_margin=("training_margin", "mean"),
        training_only_positive_votes=("training_only_prediction", "sum"),
        n_repeats=("repeat", "nunique"),
        mean_training_threshold=("training_only_youden_threshold", "mean"),
    )
    patient["training_only_fast_repeat_percentage"] = patient["training_only_positive_votes"] / patient["n_repeats"]
    strict_majority = patient["training_only_positive_votes"] > patient["n_repeats"] / 2
    tie = patient["training_only_positive_votes"] == patient["n_repeats"] / 2
    patient["training_only_prediction"] = (strict_majority | (tie & (patient["mean_training_margin"] > 0))).astype(int)
    patient["y_pred"] = patient["training_only_prediction"]
    return predictions, thresholds, repeat_auc, patient


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path, help="Directory with the approved cohort exports")
    parser.add_argument("--output-dir", required=True, type=Path, help="Private output directory outside this checkout")
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    out = args.output_dir.resolve()
    checkout = Path(__file__).resolve().parents[1]
    if any(Path(os.path.normcase(path)).is_relative_to(Path(os.path.normcase(checkout))) for path in (data_dir, out)):
        parser.error("Keep the cohort inputs and restricted outputs outside this checkout.")
    if out == data_dir:
        parser.error("Use a separate output directory from the cohort inputs.")
    out.mkdir(parents=True, exist_ok=True)

    data = load_data(data_dir)
    binary: pd.DataFrame = data["binary"]
    als: pd.DataFrame = data["als"]
    proteins: list[str] = data["protein_cols"]
    screen = protein_screen(binary, proteins)
    screen.to_csv(out / "PROTEIN_SCREENING_AUTHENTIC_793.csv", index=False)

    candidates = screen.loc[screen["Pass_nominal_p_lt_0.05_and_abs_d_ge_0.8"]].copy().reset_index(drop=True)
    if len(candidates) != 22:
        raise RuntimeError(f"Expected 22 dual-screen-positive candidates, found {len(candidates)}")
    candidates["HC_demographic_filter_applied"] = False
    candidates["Adjustment_covariates"] = "age_visit_years;sex_male"
    candidates.to_csv(out / "FILTERED_CANDIDATES.csv", index=False)

    ranking, covariate_importance = protein_ranking(binary, candidates)
    ranking.to_csv(out / "RF_GINI_RANKING.csv", index=False)
    covariate_importance.to_csv(out / "RF_GINI_ADJUSTMENT_COVARIATES.csv", index=False)
    panel_rows = []
    for size in PANEL_SIZES:
        for rank, protein in enumerate(ranking.head(size)["Protein"], start=1):
            panel_rows.append({"panel_size": size, "panel_rank": rank, "Protein": protein})
    pd.DataFrame(panel_rows).to_csv(out / "PANEL_DEFINITIONS.csv", index=False)

    all_predictions = []
    all_thresholds = []
    all_repeat = []
    all_patients = []
    all_ties = []
    performance = []
    for size in PANEL_SIZES:
        panel_proteins = ranking.head(size)["Protein"].tolist()
        pred, thresholds, repeat_auc, patient = repeated_panel(binary, panel_proteins, size)
        all_predictions.append(pred)
        all_thresholds.append(thresholds)
        all_repeat.append(repeat_auc)
        all_patients.append(patient)
        ties = patient.loc[patient["training_only_positive_votes"] == patient["n_repeats"] / 2].copy()
        ties["tie_resolved_by"] = np.where(ties["mean_training_margin"] > 0, "mean_margin_fast", "mean_margin_slow")
        all_ties.append(ties[["panel_size", "SubjID", "Group", "y_true", "training_only_positive_votes", "n_repeats", "mean_training_margin", "tie_resolved_by"]])
        patient_metrics = metrics(
            patient["y_true"].to_numpy(int),
            patient["probability_fast"].to_numpy(float),
            patient["training_only_prediction"].to_numpy(int),
        )
        performance.append({
            "panel_size": size,
            "unique_patient_auc": patient_metrics["auc"],
            "participant_averaged_auc_secondary": patient_metrics["auc"],
            "mean_repeat_auc": float(repeat_auc["auc"].mean()),
            "sd_repeat_auc": float(repeat_auc["auc"].std(ddof=1)),
            "min_repeat_auc": float(repeat_auc["auc"].min()),
            "max_repeat_auc": float(repeat_auc["auc"].max()),
            "pooled_repeated_row_auc": float(roc_auc_score(pred["y_true"], pred["probability_fast"])),
            "training_only_fast_sensitivity": patient_metrics["fast_sensitivity"],
            "training_only_slow_specificity": patient_metrics["slow_specificity"],
            "training_only_balanced_accuracy": patient_metrics["balanced_accuracy"],
            "training_only_accuracy": patient_metrics["accuracy"],
            "training_only_tn": patient_metrics["tn"],
            "training_only_fp": patient_metrics["fp"],
            "training_only_fn": patient_metrics["fn"],
            "training_only_tp": patient_metrics["tp"],
            "training_threshold_median": float(thresholds["threshold"].median()),
            "training_threshold_iqr": float(thresholds["threshold"].quantile(.75) - thresholds["threshold"].quantile(.25)),
        })

    predictions = pd.concat(all_predictions, ignore_index=True)
    thresholds = pd.concat(all_thresholds, ignore_index=True)
    repeats = pd.concat(all_repeat, ignore_index=True)
    patients = pd.concat(all_patients, ignore_index=True)
    ties = pd.concat(all_ties, ignore_index=True) if all_ties else pd.DataFrame()
    perf = pd.DataFrame(performance).sort_values("panel_size")
    predictions.to_csv(out / "OOF_PREDICTIONS_10x5.csv", index=False)
    predictions.to_csv(out / "OOF_PREDICTIONS_TRAINING_THRESHOLD.csv", index=False)
    thresholds.to_csv(out / "TRAINING_ONLY_THRESHOLDS.csv", index=False)
    thresholds.to_csv(out / "TRAINING_ONLY_THRESHOLD_RESULTS_10x5.csv", index=False)
    repeats.to_csv(out / "REPEAT_AUCS_10x5.csv", index=False)
    patients.to_csv(out / "PATIENT_LEVEL_PREDICTIONS_10x5.csv", index=False)
    ties.to_csv(out / "TIE_RESOLUTION_SUMMARY.csv", index=False)
    perf.to_csv(out / "PANEL_PERFORMANCE_10x5.csv", index=False)
    perf.to_csv(out / "TRAINING_THRESHOLD_SUMMARY.csv", index=False)

    # CPR is descriptive and includes all 43 ALS participants, including the
    # 11 Intermediate cases.  This does not create a CPR prediction model.
    cpr_rows = []
    for protein in ranking["Protein"]:
        sub = pd.DataFrame({"cpr": pd.to_numeric(als["CPR"], errors="coerce"), "protein": pd.to_numeric(als[protein], errors="coerce")}).dropna()
        from scipy.stats import spearmanr
        rho, p = spearmanr(sub["cpr"], sub["protein"])
        cpr_rows.append({"Protein": protein, "CPR_n": len(sub), "CPR_Spearman_rho": rho, "CPR_p_value": p})
    cpr = pd.DataFrame(cpr_rows)
    cpr["CPR_BH_q_value"] = bh_adjust(cpr["CPR_p_value"])
    cpr["CPR_abs_rho_ge_0.40_and_p_lt_0.05"] = (cpr["CPR_Spearman_rho"].abs() >= .40) & (cpr["CPR_p_value"] < .05)
    cpr = cpr.merge(screen[["Protein", "Pass_nominal_p_lt_0.05_and_abs_d_ge_0.8"]], on="Protein", how="left")
    cpr.to_csv(out / "CPR_CORRELATIONS_WITH_FDR.csv", index=False)

    metadata = {
        "analysis": "Random Forest Top-5/Top-7/Top-10 with age-sex adjustment under 10x5 CV",
        "screen": "all 793 proteins; two-sided MWU p < 0.05 AND |Cohen d| >= 0.8; no healthy-control protein exclusion",
        "candidate_n": int(len(candidates)),
        "covariate_adjustment": {"included": list(ADJUSTMENT_COLUMNS), "duration_included": False, "covariates_not_counted_as_panel_proteins": True},
        "cohort": {"als_first_visit": int(len(als)), "binary_n": int(len(binary)), "fast": int(binary["y"].sum()), "slow": int((binary["y"] == 0).sum()), "intermediate_for_cpr": int(als["Group"].eq("Intermediate").sum())},
        "rf_parameters": RF_PARAMS,
        "ranking": "full-cohort age/sex-adjusted Gini ranking; only protein importances determine Top-5/Top-7/Top-10 membership, while age and sex remain adjustment covariates outside panel counts",
        "cv": {"repeats": REPEATS, "folds": FOLDS, "seed": SEED, "threshold": "inner-CV Youden within each outer training set", "participant_vote_tie_rule": "mean probability-minus-threshold margin; Fast if mean margin > 0, otherwise Slow"},
        "source_files": {name: {"path": str(data_dir / name), "sha256": sha256(data_dir / name)} for name in ["clinData2.csv", "exprData.csv", "featData.csv"]},
        "caveat": "Full-cohort screening and protein ranking are development-stage locks; repeated CV estimates internal performance conditional on those locks.",
    }
    (out / "ANALYSIS_METADATA.json").write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")
    (out / "AGE_SEX_ADJUSTMENT_SUMMARY.csv").write_text(
        "feature,role\nage_visit_years,adjustment_covariate\nsex_male,adjustment_covariate\nduration_months,excluded_from_model\n",
        encoding="utf-8",
    )
    print(json.dumps({"out": str(out), "candidate_n": len(candidates), "ranking": ranking[["RF_Gini_Rank", "Protein"]].to_dict("records"), "performance": performance}, indent=2, default=str))


if __name__ == "__main__":
    main()
