#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

OUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_scatter_median_whisker/"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PDF = os.path.join(
    OUT_DIR,
    "Precision_Recall_median_whisker_3x4_with_F1.pdf"
)

OUT_CSV = os.path.join(
    OUT_DIR,
    "Precision_Recall_median_whisker_3x4_summary.csv"
)

ALIGNER_ORDER = ["hisat2", "STAR", "tophat2"]

ALIGNER_LABELS = {
    "hisat2": "HISAT2",
    "STAR": "STAR",
    "tophat2": "TopHat2",
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

ANNOTATE_TOOL_NAME = False


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
    df["category"] = category
    df["type"] = type_name

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
            recall_min=("recall", "min"),
            recall_max=("recall", "max"),
            precision_median=("precision", "median"),
            precision_min=("precision", "min"),
            precision_max=("precision", "max"),
            n=("recall", "size"),
        )
    )
    return summary


# ============================================================

def add_f1_curves(ax, f1_scores=(0.3, 0.5, 0.7, 0.9)):
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
            linewidth=0.8,
            color="gray",
            alpha=0.45,
            zorder=0,
        )

        ax.text(
            r[-1] - 0.02,
            p[-1] + 0.015,
            f"F1={f1}",
            fontsize=7,
            color="gray",
            ha="right",
            va="bottom",
            alpha=0.8,
        )


# ============================================================

def plot_precision_recall(summary_df):
    levels = [x["level"] for x in LEVEL_CONFIGS]
    aligners = ALIGNER_ORDER

    existing_tools = set(summary_df["tools"].unique())
    tools_order = [x for x in TOOL_ORDER if x in existing_tools]
    remaining_tools = sorted(existing_tools - set(tools_order))
    tools_order.extend(remaining_tools)

    n_rows = len(aligners)
    n_cols = len(levels)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.2 * n_cols * 0.75, 4.0 * n_rows * 0.75),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    for i, aligner in enumerate(aligners):
        for j, level in enumerate(levels):
            ax = axes[i][j]

            sub = summary_df[
                (summary_df["level"] == level) &
                (summary_df["soft"] == aligner)
            ].copy()

            add_f1_curves(ax)

            if sub.empty:
                ax.set_title(
                    f"{level} | {ALIGNER_LABELS.get(aligner, aligner)}\nNo data",
                    fontsize=12,
                    fontweight="bold",
                )
                ax.set_xlim(-0.05, 1.05)
                ax.set_ylim(-0.05, 1.05)
                continue

            sub["tools"] = pd.Categorical(
                sub["tools"],
                categories=tools_order,
                ordered=True
            )
            sub = sub.sort_values("tools")

            for _, row in sub.iterrows():
                tool = str(row["tools"])
                color = ORF_DETECTOR_COLOR.get(tool, "#333333")

                x = row["recall_median"]
                y = row["precision_median"]

                xmin = row["recall_min"]
                xmax = row["recall_max"]
                ymin = row["precision_min"]
                ymax = row["precision_max"]

                ax.hlines(
                    y=y,
                    xmin=xmin,
                    xmax=xmax,
                    color=color,
                    linewidth=1.3,
                    alpha=0.9,
                    zorder=2,
                )

                ax.vlines(
                    x=x,
                    ymin=ymin,
                    ymax=ymax,
                    color=color,
                    linewidth=1.3,
                    alpha=0.9,
                    zorder=2,
                )

                # 散点去除黑色边框
                ax.scatter(
                    x,
                    y,
                    s=44,
                    color=color,
                    edgecolors="none",
                    linewidths=0,
                    alpha=0.95,
                    zorder=3,
                )

                if ANNOTATE_TOOL_NAME:
                    ax.text(
                        x + 0.008,
                        y + 0.008,
                        tool,
                        fontsize=6,
                        color=color,
                    )

            ax.set_title(
                f"{level} | {ALIGNER_LABELS.get(aligner, aligner)}",
                fontsize=12,
                fontweight="bold",
            )

            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)

            ax.grid(
                True,
                linestyle="--",
                linewidth=0.4,
                alpha=0.35,
                zorder=0,
            )

            if i == n_rows - 1:
                ax.set_xlabel("Recall", fontsize=11)

            if j == 0:
                ax.set_ylabel("Precision", fontsize=11)

            ax.tick_params(axis="both", labelsize=9)

    legend_handles = []
    legend_labels = []

    for tool in tools_order:
        color = ORF_DETECTOR_COLOR.get(tool, "#333333")

        handle = plt.Line2D(
            [0],
            [0],
            marker="o",
            color=color,
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=6,
            linewidth=0,
        )

        legend_handles.append(handle)
        legend_labels.append(tool)

    fig.legend(
        legend_handles,
        legend_labels,
        title="ORF detection tool",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        fontsize=8,
        title_fontsize=10,
        frameon=False,
    )

    fig.suptitle(
        "Precision–Recall performance of ORF detection tools",
        fontsize=15,
        fontweight="bold",
        y=1.0,
    )

    plt.tight_layout(rect=[0, 0, 0.84, 0.96])

    print(f"[SAVE] {OUT_PDF}")
    plt.savefig(
        OUT_PDF,
        format="pdf",
        dpi=400,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================

def main():
    print("=== Starting precision-recall median-whisker plot ===")

    all_df = []

    for config in LEVEL_CONFIGS:
        df = read_one_level(config)
        all_df.append(df)

    all_df = pd.concat(all_df, ignore_index=True)

    summary_df = summarize_precision_recall(all_df)

    print(f"[SAVE] {OUT_CSV}")
    summary_df.to_csv(OUT_CSV, index=False)

    plot_precision_recall(summary_df)

    print("=== Finished ===")


if __name__ == "__main__":
    main()
