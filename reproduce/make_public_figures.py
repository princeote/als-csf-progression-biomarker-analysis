#!/usr/bin/env python3
"""Generate the aggregate CPR and protein-effect figures for the public paper."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


TEAL = "#008587"
BLUE = "#2463A6"
INK = "#172334"
GRID = "#E2E8F0"
GRAY = "#52647A"


def plot_signed_bars(ax: plt.Axes, labels: pd.Series, values: pd.Series) -> None:
    colors = [TEAL if value >= 0 else BLUE for value in values]
    ax.barh(labels, values, height=0.62, color=colors, zorder=3)
    ax.invert_yaxis()
    ax.axvline(0, color=INK, linewidth=0.9, zorder=4)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", pad=7, length=0)
    for label, value in zip(labels, values):
        offset = 5 if value >= 0 else -5
        alignment = "left" if value >= 0 else "right"
        ax.annotate(f"{value:+.3f}", (value, label), xytext=(offset, 0),
                    textcoords="offset points", ha=alignment, va="center", color=INK)


def plot_cpr(results: Path, out: Path) -> None:
    cpr = pd.read_csv(results / "CPR_CORRELATIONS_WITH_FDR.csv")
    highlighted = cpr.loc[cpr["CPR_abs_rho_ge_0.40_and_p_lt_0.05"]]
    highlighted = highlighted.sort_values("CPR_Spearman_rho", ascending=False)
    labels = highlighted["Short_Name"]
    fig, (ax, stats) = plt.subplots(
        ncols=2, figsize=(8.1, 3.1), gridspec_kw={"width_ratios": [3.7, 1.65]},
        sharey=True, layout="constrained",
    )
    plot_signed_bars(ax, labels, highlighted["CPR_Spearman_rho"])
    ax.set_xlim(-0.72, 0.72)
    ax.set_title("CPR associations among screen-positive proteins", loc="left")
    ax.set_xlabel("Spearman ρ (43 ALS participants)")
    stats.set_title("Nominal p       BH q", fontweight="normal")
    stats.set_axis_off()
    stats.set_xlim(0, 1)
    for protein in highlighted.itertuples(index=False):
        stats.text(0.30, protein.Short_Name, f"{protein.CPR_p_value:.3g}",
                   ha="center", va="center", color=INK)
        stats.text(0.73, protein.Short_Name, f"{protein.CPR_BH_q_value:.3g}",
                   ha="center", va="center", color=INK)
    fig.savefig(out / "fig4_cpr_summary.png", dpi=300, facecolor="white")
    plt.close(fig)


def plot_effects(results: Path, out: Path) -> None:
    evidence = pd.read_csv(results / "FINAL_TOP10_EVIDENCE.csv").sort_values("Gini_Rank")
    fig, ax = plt.subplots(figsize=(8.1, 4.1), layout="constrained")
    plot_signed_bars(ax, evidence["Short_Name"], evidence["Cohens_d_fast_minus_slow"])
    for cutoff in (-0.8, 0.8):
        ax.axvline(cutoff, color=GRAY, linewidth=0.9, linestyle="--", zorder=2)
    ax.set_xlim(-1.35, 1.42)
    ax.set_title("Fast versus Slow effect sizes for the Top-10 proteins", loc="left")
    ax.set_xlabel("Cohen’s d (Fast minus Slow)")
    fig.savefig(out / "fig6_protein_effects.png", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=root / "results")
    parser.add_argument("--output-dir", type=Path, default=root / "paper" / "figures")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "axes.edgecolor": GRAY,
        "xtick.color": GRAY,
        "ytick.color": INK,
    })
    plot_cpr(args.results_dir, args.output_dir)
    plot_effects(args.results_dir, args.output_dir)


if __name__ == "__main__":
    main()
