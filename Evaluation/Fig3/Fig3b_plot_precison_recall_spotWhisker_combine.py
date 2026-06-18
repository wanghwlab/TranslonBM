#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# ============================================================

LEVEL_CONFIGS = [
    {
        "level": "TIS",
        "base_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_Blocks/",
        "category": "merged_ATG",
        "type": "orf_pred_default_trim",
    },
    {
        "level": "TI-seq",
        "base_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_TISeq/",
        "category": "merged_ATG",
        "type": "orf_pred_default_trim",
    },
    {
        "level": "MaxQuant (protein level)",
        "base_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_Maxquant/",
        "category": "fdr0.1_score_protein",
        "type": "orf_pred_default_trim",
    },
    {
        "level": "MaxQuant (peptide level)",
        "base_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_Maxquant/",
        "category": "fdr0.1_score_peptide",
        "type": "orf_pred_default_trim",
    },
]

OUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_scatter_median_whisker"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PDF = os.path.join(OUT_DIR, "Precision_Recall_combined_all_levels_with_F1.pdf")
OUT_PNG = os.path.join(OUT_DIR, "Precision_Recall_combined_all_levels_with_F1.png")
OUT_CSV = os.path.join(OUT_DIR, "Precision_Recall_combined_all_levels_summary.csv")

ALIGNER_MARKERS = {
    "hisat2": "s",   # square
    "STAR": "o",     # circle
    "tophat2": "^",  # triangle
}

ALIGNER_LABELS = {
    "hisat2": "HISAT2",
    "STAR": "STAR",
    "tophat2": "TopHat2",
}

ALIGNER_ORDER = ["hisat2", "STAR", "tophat2"]

LEVEL_ORDER = [
    "TIS",
    "TI-seq",
    "MaxQuant (protein level)",
    "MaxQuant (peptide level)",
]

LEVEL_SIZE = {
    "TIS": 95,
    "TI-seq": 75,
    "MaxQuant (protein level)": 58,
    "MaxQuant (peptide level)": 42,
}

TOOL_NAME_MAP = {
    "riboseqc": "Ribo-seQC",
    "Ribo-seQC": "Ribo-seQC",
    "gedi": "PRICE",
    "PRICE": "PRICE",
    "ribocode": "RiboCode",
    "RiboCode": "RiboCode",
    "ribotish": "Ribo-TISH",
    "Ribo-TISH": "Ribo-TISH",
    "ribotricer": "Ribotricer",
    "Ribotricer": "Ribotricer",
    "ribotaper": "RiboTaper",
    "RiboTaper": "RiboTaper",
    "ribowave": "RiboWave",
    "RiboWave": "RiboWave",
    "riborf": "RibORF",
    "RibORF": "RibORF",
    "rpbp": "RP-BP",
    "RP-BP": "RP-BP",
    "ORFquant": "ORFquant",
    "orfquant": "ORFquant",
    "orfrater": "ORF-RATER",
    "ORF-RATER": "ORF-RATER",
    "ribohmm": "riboHMM",
    "riboHMM": "riboHMM",
}

ORF_DETECTOR_COLOR = {
    "Ribo-seQC": "#73aaa4",
    "PRICE": "#86BA84",
    "RiboCode": "#53A364",
    "Ribo-TISH": "#C97582",
    "Ribotricer": "#4379B5",
    "RiboTaper": "#7CBCA3",
    "RiboWave": "#7B5B97",
    "RibORF": "#B99CB1",
    "RP-BP": "#699F8C",
    "ORFquant": "#577175",
    "ORF-RATER": "#796653",
    "riboHMM": "#E08B42",
}

TOOL_ORDER = [
    "Ribo-seQC",
    "PRICE",
    "RiboCode",
    "Ribo-TISH",
    "Ribotricer",
    "RiboTaper",
    "RiboWave",
    "RibORF",
    "RP-BP",
    "ORFquant",
    "ORF-RATER",
    "riboHMM",
]

EXCLUDE_SAMPLES = ["SRX1447296"]


# ============================================================

def build_input_path(base_dir, category, type_name):
    category_suffix = category.replace("final_ORFs_", "")
    input_csv_name = f"{type_name}_score_{category_suffix}.csv"
    input_csv_path = os.path.join(base_dir, category_suffix, input_csv_name)
    return input_csv_path


def standardize_tool_name(tool):
    tool = str(tool)
    return TOOL_NAME_MAP.get(tool, tool)


def read_one_level(config):
    level = config["level"]
    base_dir = config["base_dir"]
    category = config["category"]
    type_name = config["type"]

    input_csv = build_input_path(base_dir, category, type_name)

    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    print(f"[READ] {level}: {input_csv}")

    df = pd.read_csv(input_csv)

    required_cols = {"sample", "soft", "tools", "precision", "recall"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{input_csv} 缺少必要列: {missing}")

    df = df.copy()
    df["level"] = level
    df["tools"] = df["tools"].apply(standardize_tool_name)

    if EXCLUDE_SAMPLES:
        df = df[~df["sample"].isin(EXCLUDE_SAMPLES)]

    df = df[df["soft"].isin(ALIGNER_ORDER)]
    df = df.dropna(subset=["precision", "recall"])

    return df


def summarize_precision_recall(df):
    summary = (
        df.groupby(["level", "soft", "tools"], as_index=False)
        .agg(
            recall_median=("recall", "median"),
            precision_median=("precision", "median"),
            recall_min=("recall", "min"),
            recall_max=("recall", "max"),
            precision_min=("precision", "min"),
            precision_max=("precision", "max"),
            n=("recall", "size"),
        )
    )
    return summary


# ============================================================

def add_f1_curves(ax, f1_scores=(0.1, 0.3, 0.5, 0.7, 0.9)):
    """
    在 precision-recall 图中添加 F1-score 等值线。

    F1 = 2 * precision * recall / (precision + recall)

    给定 F1 和 recall，反推 precision：
    precision = F1 * recall / (2 * recall - F1)
    """
    recall = np.linspace(0.001, 1.0, 1000)

    for f1 in f1_scores:
        valid = recall > f1 / 2
        r = recall[valid]
        p = (f1 * r) / (2 * r - f1)

        valid_p = (p >= 0) & (p <= 1)
        r = r[valid_p]
        p = p[valid_p]

        if len(r) == 0:
            continue

        ax.plot(
            r,
            p,
            linestyle="--",
            linewidth=0.9,
            color="gray",
            alpha=0.45,
            zorder=1,
        )

        ax.text(
            r[-1] - 0.015,
            p[-1] + 0.012,
            f"F1={f1}",
            fontsize=8,
            color="gray",
            ha="right",
            va="bottom",
            alpha=0.85,
        )


# ============================================================

def plot_combined(summary_df):
    fig, ax = plt.subplots(figsize=(7, 5.5))

    existing_tools = set(summary_df["tools"].unique())
    tools_order = [x for x in TOOL_ORDER if x in existing_tools]
    remaining_tools = sorted(existing_tools - set(tools_order))
    tools_order.extend(remaining_tools)

    summary_df["tools"] = pd.Categorical(
        summary_df["tools"],
        categories=tools_order,
        ordered=True
    )
    summary_df["level"] = pd.Categorical(
        summary_df["level"],
        categories=LEVEL_ORDER,
        ordered=True
    )
    summary_df = summary_df.sort_values(["level", "soft", "tools"])

    add_f1_curves(ax)

    for _, row in summary_df.iterrows():
        tool = str(row["tools"])
        level = str(row["level"])
        aligner = str(row["soft"])

        color = ORF_DETECTOR_COLOR.get(tool, "#333333")
        marker = ALIGNER_MARKERS.get(aligner, "o")
        size = LEVEL_SIZE.get(level, 70)

        ax.scatter(
            row["recall_median"],
            row["precision_median"],
            s=size,
            c=color,
            marker=marker,
            edgecolors="none",
            linewidths=0,
            alpha=0.9,
            zorder=3
        )

    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Recall", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.35, zorder=0)

    ax.set_title(
        "Precision–Recall median performance of ORF detection tools\n"
        "across TIS, TI-seq, MaxQuant protein-level and peptide-level data",
        fontsize=13,
        fontweight="bold"
    )

    # ========================================================
    tool_handles = []
    for tool in tools_order:
        handle = Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            markerfacecolor=ORF_DETECTOR_COLOR.get(tool, "#333333"),
            markeredgecolor="none",
            markersize=7,
            label=tool
        )
        tool_handles.append(handle)

    legend1 = ax.legend(
        handles=tool_handles,
        title="ORF detection tool (color)",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.00),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
        ncol=1
    )
    ax.add_artist(legend1)

    # ========================================================
    aligner_handles = []
    for aligner in ALIGNER_ORDER:
        handle = Line2D(
            [0], [0],
            marker=ALIGNER_MARKERS[aligner],
            linestyle="",
            markerfacecolor="#666666",
            markeredgecolor="none",
            markersize=8,
            label=ALIGNER_LABELS[aligner]
        )
        aligner_handles.append(handle)

    legend2 = ax.legend(
        handles=aligner_handles,
        title="Aligner (shape)",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.53),
        frameon=False,
        fontsize=9,
        title_fontsize=10
    )
    ax.add_artist(legend2)

    # ========================================================
    level_handles = []
    for level in LEVEL_ORDER:
        ms = (LEVEL_SIZE[level] ** 0.5) * 0.9
        handle = Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            markerfacecolor="#888888",
            markeredgecolor="none",
            markersize=ms,
            label=level
        )
        level_handles.append(handle)

    legend3 = ax.legend(
        handles=level_handles,
        title="Data level (point size)",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.28),
        frameon=False,
        fontsize=9,
        title_fontsize=10
    )
    ax.add_artist(legend3)

    plt.tight_layout(rect=[0, 0, 0.78, 1])

    print(f"[SAVE] {OUT_PDF}")
    plt.savefig(OUT_PDF, dpi=400, bbox_inches="tight")

    print(f"[SAVE] {OUT_PNG}")
    plt.savefig(OUT_PNG, dpi=400, bbox_inches="tight")

    plt.close(fig)


# ============================================================

def main():
    print("=== Start plotting combined precision-recall scatter ===")

    all_df = []
    for config in LEVEL_CONFIGS:
        df = read_one_level(config)
        all_df.append(df)

    all_df = pd.concat(all_df, ignore_index=True)

    summary_df = summarize_precision_recall(all_df)

    print(f"[SAVE] {OUT_CSV}")
    summary_df.to_csv(OUT_CSV, index=False)

    plot_combined(summary_df)

    print("=== Finished ===")


if __name__ == "__main__":
    main()
