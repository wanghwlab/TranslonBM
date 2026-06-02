#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- 配置区 ---

BASE_DIRS_TO_PROCESS = [
   # "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_TISeq/",
    "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_group/PR_recall_Blocks/",
   # "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_Maxquant/"
]

CATEGORIES_TO_PROCESS = [
    "merged_ATG",
    "fdr0.1_score_peptide",
    "fdr0.1_score_protein",
]

TYPES_TO_PROCESS = [
    "orf_pred_default_untrim",
    "orf_pred_default_trim",
]

# 修改点 1：颜色配置 (STAR: #a0c5b4, HISAT2: #e3ab76, TopHat2: #ca7e7b)
SOFT_ORDER = ['hisat2', 'STAR', 'tophat2']
SOFT_COLORS = {
    "hisat2": "#e3ab76",
    "STAR": "#a0c5b4",
    "tophat2": "#ca7e7b"
}

# --- 核心繪圖函數 ---

def generate_plots_for_csv(input_csv, output_dir, category_name, category_label, trim_label):
    print("-" * 80)
    print(f"Processing CSV file: {input_csv}")
    
    if not os.path.exists(input_csv):
        print(f"[WARNING] Input file not found. Skipping.")
        return

    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(input_csv)

    # 過濾數據
    df = df[df['sample'] != 'SRX1447296']
    df = df[df['soft'].isin(SOFT_ORDER)]
    
    sample_groups = {
        'Replicate': df[df['sample'].str.startswith('SRX7666')],
        'Simulation': df[df['sample'].str.contains('simulation')],
        'Real': df[(df['sample'].str.startswith('SRX')) & (~df['sample'].str.startswith('SRX7666'))]
    }
    
    sns.set_theme(style="ticks")
    plt.rcParams.update({
        'text.color': 'black', 'axes.labelcolor': 'black',
        'xtick.color': 'black', 'ytick.color': 'black',
        'axes.edgecolor': 'black'
    })

    if "non_annotated" in category_name or "annotated" in category_name:
        metrics_to_plot = ['F1-score']
        num_subplots = 1
        fig_height = 4
    else:
        metrics_to_plot = ['precision', 'recall', 'F1-score']
        num_subplots = 3
        fig_height = 9

    for group_name, df_group in sample_groups.items():
        if df_group.empty:
            continue
        print(f"  Processing sample group: {group_name}")
        
        sorted_tools = (
            df_group.groupby('tools')['F1-score']
            .mean()
            .sort_values(ascending=False)
            .index.tolist()
        )

        fig_width = max(6, len(sorted_tools) * 0.4)
        fig, axes = plt.subplots(num_subplots, 1, figsize=(fig_width, fig_height), sharex=True, squeeze=False)
        axes = axes.flatten()

        fig.suptitle(f'Performance Comparison: Aligner Grouping ({group_name})', fontsize=16, y=1.02)

        for i, metric in enumerate(metrics_to_plot):
            ax = axes[i]
            
            sns.boxplot(
                x='tools', 
                y=metric, 
                hue='soft', 
                data=df_group, 
                order=sorted_tools,  
                hue_order=SOFT_ORDER,
                palette=SOFT_COLORS,
                ax=ax, 
                width=0.6,
                showfliers=False, 
                linewidth=1.2,
                boxprops={'edgecolor': 'black'},
                whiskerprops={'color': 'black'},
                capprops={'color': 'black'},
                medianprops={'color': 'black'}
            )
            
            ax.set_ylabel(metric.capitalize(), fontsize=14)
            ax.set_xlabel('')
            ax.tick_params(axis='y', labelsize=12)
            ax.set_ylim(-0.05, 1.05)
            sns.despine(ax=ax)
            
            if ax.get_legend() is not None:
                ax.get_legend().remove()

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, title='Aligner', bbox_to_anchor=(1.01, 0.95), loc='upper left', fontsize=12, title_fontsize=14)
        
        axes[-1].set_xlabel('ORF Prediction Tool', fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.tight_layout(rect=[0, 0.03, 0.9, 0.98])

        output_filename_base = f"boxplot_combined_{category_label}_{trim_label}_{group_name}.pdf"
        output_filename = os.path.join(output_dir, output_filename_base)
        
        print(f"      Saving combined plot to {output_filename}")
        plt.savefig(output_filename, format='pdf', bbox_inches='tight', dpi=400)
        plt.close(fig)

# --- 主執行函數 ---
def main():
    print("--- Starting Combined Plotting Script ---")
    
    for base_input_dir in BASE_DIRS_TO_PROCESS:
        category_label = base_input_dir.rstrip('/').split('_')[-1]
        base_output_dir = base_input_dir.rstrip('/') + '_plots_combined2/'
        print(f"\nProcessing Base Directory: {base_input_dir} (Category: {category_label})")
        
        for category_name in CATEGORIES_TO_PROCESS:
            for type_name in TYPES_TO_PROCESS:
                trim_label = type_name.split('_')[-1]
                category_suffix = category_name.replace("final_ORFs_", "")
                input_csv_name = f"{type_name}_score_{category_suffix}.csv"
                input_csv_path = os.path.join(base_input_dir, category_suffix, input_csv_name)
                output_dir_for_csv = os.path.join(base_output_dir, category_suffix, type_name)
                
                generate_plots_for_csv(input_csv_path, output_dir_for_csv, category_name, category_label, trim_label)

    print("\n--- Script Finished ---")

if __name__ == "__main__":
    main()
