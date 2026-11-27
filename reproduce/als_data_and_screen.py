"""Data loading, sample-wise normalization, and unadjusted dual protein screen.
No healthy-control exclusion or alternative classifier workflow is run here.
Functions are preserved from the verified source used for the adjusted RF results.
"""
from __future__ import annotations
import math,re
import hashlib
from io import BytesIO
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
FAST=1
SLOW=0

# Author-confirmed release: measurement columns correspond to the clinical rows
# and measurement rows correspond to the feature rows, in this exact order.
SOURCE_RDATA_SHA256 = "809dd42a75622a97dc5122420eefda95fcf7a14de9c47d62a27aeb3040f8bc7d"
AUDITED_EXPORTS = {
    "clinData2.csv": "9670c88e04875a5f865243bd0a8d2c083dad290f1dab49565ad76eb27428ae06",
    "exprData.csv": "3a01dd71d844354a1a114b10c1a2a6f016fba22d60547d0fbd4628389b5a58e3",
    "featData.csv": "8a7146cc63a8f35c8ab37d2fcc461e2e719779605dcc8fcc941b0631bb8c3222",
}

def audited_table(data_dir: Path, name: str, **csv_options) -> pd.DataFrame:
    data = (Path(data_dir) / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != AUDITED_EXPORTS[name]:
        raise ValueError(f"Cannot establish sample/feature alignment: {name} differs from the "
                         "source-verified export. Restore the audited release or verify the new "
                         "release against an author-supplied sample and feature key.")
    return pd.read_csv(BytesIO(data), **csv_options)
def clean_name(value: object) -> str:
    text = "Unknown" if pd.isna(value) else str(value).split("OS=")[0].strip()
    text = text.replace(" (Fragment)", "")
    return re.sub(r"\s+", " ", text)

def unique_names(names: Iterable[str]) -> list[str]:
    seen: dict[str, int] = {}
    result: list[str] = []
    for name in names:
        count = seen.get(name, 0)
        result.append(name if count == 0 else f"{name}__dup{count}")
        seen[name] = count + 1
    return result

def bh_adjust(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    n = len(p)
    order = np.argsort(np.where(np.isfinite(p), p, 1.0))
    ranked = np.where(np.isfinite(p[order]), p[order], 1.0)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(q, 1.0)
    return out

def cohen_d(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna().astype(float)
    y = pd.to_numeric(y, errors="coerce").dropna().astype(float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    denom = len(x) + len(y) - 2
    pooled = math.sqrt(max(((len(x) - 1) * sx**2 + (len(y) - 1) * sy**2) / denom, 0.0))
    return float((x.mean() - y.mean()) / pooled) if pooled > 0 else 0.0

def mw_p(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna().astype(float)
    y = pd.to_numeric(y, errors="coerce").dropna().astype(float)
    if len(x) < 3 or len(y) < 3:
        return 1.0
    return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)

def first_visits(clin: pd.DataFrame, group: str) -> pd.DataFrame:
    sub = clin.loc[clin["Patient.Group"].eq(group)].copy()
    sub["original_row_index"] = sub.index
    sub = sub.sort_values(["SubjID", "EncounterNumber"])
    return sub.drop_duplicates("SubjID", keep="first").reset_index(drop=True)

def add_covariates(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["age_visit_years"] = pd.to_numeric(out["ageVisit"], errors="coerce") / 365.25
    out["duration_months"] = pd.to_numeric(out["Duration.in.months"], errors="coerce")
    out["sex_male"] = out["Gender"].astype(str).str.upper().str.strip().eq("M").astype(int)
    return out

def normalize_expression(expr_raw: pd.DataFrame) -> pd.DataFrame:
    # This is the original project transformation: log2(value + 1), followed
    # by sample-wise median centering across proteins.
    log2 = np.log2(expr_raw.astype(float) + 1.0)
    medians = np.nanmedian(log2, axis=0)
    return log2.sub(medians, axis=1)

def load_data(data_dir: Path) -> dict[str, object]:
    clin = audited_table(data_dir, "clinData2.csv")
    expr_raw = audited_table(data_dir, "exprData.csv", header=None)
    feat = audited_table(data_dir, "featData.csv")
    keys = clin[["SubjID", "EncounterNumber"]]
    if keys.isna().any().any() or keys.duplicated().any():
        raise ValueError("Each clinical row requires a unique participant and encounter.")
    if expr_raw.shape[1] != len(clin) or expr_raw.shape[0] != len(feat):
        raise ValueError(f"Data alignment failed: expr={expr_raw.shape}, clin={clin.shape}, feat={feat.shape}")
    expr = normalize_expression(expr_raw)
    clean = [clean_name(x) for x in feat["DESCRIPTION"]]
    proteins = unique_names([f"P{i}_{clean[i]}" for i in range(len(clean))])
    expr_t = expr.T.copy()
    expr_t.columns = proteins
    meta = feat.copy()
    meta.insert(0, "Protein", proteins)
    meta.insert(1, "Clean_Name", clean)

    def attach(base: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for _, row in base.iterrows():
            idx = int(row["original_row_index"])
            rows.append(pd.concat([row.drop(labels=["original_row_index"]), expr_t.iloc[idx]]))
        return pd.DataFrame(rows).reset_index(drop=True)

    als = add_covariates(attach(first_visits(clin, "ALS")))
    hc = add_covariates(attach(first_visits(clin, "HC")))
    als["CPR"] = pd.to_numeric(als["Rate.of.progression"], errors="coerce")
    als["Group"] = np.where(als["CPR"] >= 1.0, "Fast", np.where(als["CPR"] <= 0.5, "Slow", "Intermediate"))
    binary = als.loc[als["Group"].isin(["Fast", "Slow"])].copy().reset_index(drop=True)
    binary["y"] = binary["Group"].eq("Fast").astype(int)
    return {
        "clin": clin,
        "expr_raw": expr_raw,
        "feat": feat,
        "expr_norm": expr,
        "protein_meta": meta,
        "protein_cols": proteins,
        "als": als,
        "hc": hc,
        "binary": binary,
    }

def protein_screen(binary: pd.DataFrame, proteins: list[str]) -> pd.DataFrame:
    fast = binary.loc[binary["y"].eq(FAST)]
    slow = binary.loc[binary["y"].eq(SLOW)]
    rows = []
    for protein in proteins:
        xf, xs = fast[protein], slow[protein]
        p = mw_p(xf, xs)
        d = cohen_d(xf, xs)
        rows.append({
            "Protein": protein,
            "P_value_MWU": p,
            "Cohens_d_fast_minus_slow": d,
            "Abs_Cohens_d": abs(d),
            "Fast_n": int(pd.to_numeric(xf, errors="coerce").notna().sum()),
            "Slow_n": int(pd.to_numeric(xs, errors="coerce").notna().sum()),
            "Mean_Fast": float(pd.to_numeric(xf, errors="coerce").mean()),
            "Mean_Slow": float(pd.to_numeric(xs, errors="coerce").mean()),
            "Median_Fast": float(pd.to_numeric(xf, errors="coerce").median()),
            "Median_Slow": float(pd.to_numeric(xs, errors="coerce").median()),
            "Direction": "higher_in_fast" if d > 0 else "lower_in_fast" if d < 0 else "none",
        })
    out = pd.DataFrame(rows)
    out["BH_q_value"] = bh_adjust(out["P_value_MWU"])
    out["Pass_nominal_p_lt_0.05_and_abs_d_ge_0.8"] = (out["P_value_MWU"] < 0.05) & (out["Abs_Cohens_d"] >= 0.8)
    return out.sort_values(["Pass_nominal_p_lt_0.05_and_abs_d_ge_0.8", "Abs_Cohens_d"], ascending=[False, False]).reset_index(drop=True)
