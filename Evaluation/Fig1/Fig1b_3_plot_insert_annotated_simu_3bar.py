#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

def linear_transform(x, pos):
    return x * 100

name_mapping = {
    "simulation_6M_T1": "simulation_6M_T1",
    "simulation_6M_T3": "simulation_6M_T3",
    "simulation_60M_T1": "simulation_60M_T1",
    "simulation_60M_T3": "simulation_60M_T3",

}

color_map_name = 'Set3'
cmap = get_cmap(color_map_name)
colors = [to_hex(cmap(i)) for i in [4, 5, 6, 2, 3]]

# Define the base directory
base_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/tools_overlap/"
#annotation_types = ['merged_ATG/orf_pred_default_untrim', 'merged_NTG/orf_pred_default_untrim']
annotation_types = ['merged_canonical_ATG/orf_pred_default_untrim','merged_non_canonical_ATG/orf_pred_default_untrim','merged_canonical_ATG/orf_pred_default_trim','merged_non_canonical_ATG/orf_pred_default_trim']

softwares_order = ['hisat2', 'STAR', 'tophat2']
soft_colors = {"hisat2": "#E64B35", "STAR": "#4DBBD5", "tophat2": "#00A087"} # 仅作参考，脚本主要使用堆叠颜色

def get_sort_key_from_id(sim_id):
    keys = list(name_mapping.keys())
    for i, key in enumerate(keys):
        if sim_id == key:
            return i
    return float('inf')

def process_dataframe(df, annotation_type):
    if 'NTG' in annotation_type and "predicted by 11 soft ORF num" in df.columns:
        df["predicted by 11 soft ORF num"] = df["predicted by 11 soft ORF num"] - 1

    numeric_columns = [
        "predicted by 1 soft ORF num", "predicted by 2 soft ORF num", "predicted by 3 soft ORF num",
        "predicted by 4 soft ORF num", "predicted by 5 soft ORF num", "predicted by 6 soft ORF num",
        "predicted by 7 soft ORF num", "predicted by 8 soft ORF num", "predicted by 9 soft ORF num",
        "predicted by 10 soft ORF num", "predicted by 11 soft ORF num",
    ]
    for col in numeric_columns:
        if col not in df.columns:
            df[col] = 0
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors="coerce")

    df["2-4 tools"] = df["predicted by 2 soft ORF num"] + df["predicted by 3 soft ORF num"] + df["predicted by 4 soft ORF num"]
    df["5-7 tools"] = df["predicted by 7 soft ORF num"] + df["predicted by 5 soft ORF num"] + df["predicted by 6 soft ORF num"]
    df["8-10 tools"] = df["predicted by 8 soft ORF num"] + df["predicted by 9 soft ORF num"] + df["predicted by 10 soft ORF num"]
    df["11 tools"] = df["predicted by 11 soft ORF num"]
    
    df.rename(columns={"predicted by 1 soft ORF num": "1 tool"}, inplace=True)
    plot_categories = ["1 tool", "2-4 tools", "5-7 tools", "8-10 tools", "11 tools"]
    df["Total"] = df[plot_categories].sum(axis=1)
    
    df['tool_label'] = df.iloc[:, 0].str.split('_').str[-2]
    
    return df, plot_categories

for annotation in annotation_types:
    all_files = glob.glob(f"{base_dir}/{annotation}/simu*_overlap_statistics.csv")
    sim_ids = sorted(list(set([os.path.basename(f).split('_hisat2')[0].split('_STAR')[0].split('_tophat2')[0] for f in all_files])), key=get_sort_key_from_id)
    
    if not sim_ids:
        continue

    global_max_total = 0
    global_max_inset_total = 0
    
    fig, axes = plt.subplots(1, len(sim_ids), figsize=(6 * len(sim_ids), 6), sharey=False)
    if len(sim_ids) == 1: axes = [axes]

    for idx, (sim_id, ax) in enumerate(zip(sim_ids, axes)):
        dataset_name = name_mapping.get(sim_id, sim_id)
        combined_data = {}
        all_tools = []
        
        for soft in softwares_order:
            file_path = glob.glob(f"{base_dir}/{annotation}/{sim_id}_{soft}_overlap_statistics.csv")
            if file_path:
                df, plot_categories = process_dataframe(pd.read_csv(file_path[0]), annotation)
                if soft == 'hisat2':
                    df = df.sort_values(by="Total", ascending=False)
                    all_tools = df['tool_label'].tolist()
                combined_data[soft] = df
        
        if not combined_data: continue

        bar_width = 0.25
        x_base = np.arange(len(all_tools))
        #offsets = [-bar_width, 0, bar_width]
        offsets = [-0.3, 0, 0.3]

        for s_idx, soft in enumerate(softwares_order):
            if soft not in combined_data: continue
            df = combined_data[soft]
            df = df.set_index('tool_label').reindex(all_tools).fillna(0)
            
            bottom = np.zeros(len(all_tools))
            x_pos = x_base + offsets[s_idx]
            
            for c_idx, category in enumerate(plot_categories):
                values = df[category].values
                ax.bar(x_pos, values, width=bar_width, bottom=bottom, color=colors[c_idx], align='center')
                bottom += values
            
            global_max_total = max(global_max_total, df["Total"].max())
            
        ax.set_title(f"{dataset_name}", fontsize=14, fontweight="bold")
        ax.set_xticks(x_base)
        ax.set_xticklabels(all_tools, rotation=45, ha="right")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if idx > 0: ax.set_ylabel("")
        else: ax.set_ylabel("Total number", fontsize=12)

        inset_tools = all_tools[-4:]
        ax_ins = ax.inset_axes([0.5, 0.55, 0.48, 0.42])
        x_ins_base = np.arange(len(inset_tools))
        
        for s_idx, soft in enumerate(softwares_order):
            if soft not in combined_data: continue
            df_ins = combined_data[soft].set_index('tool_label').reindex(inset_tools).fillna(0)
            bottom_ins = np.zeros(len(inset_tools))
            x_ins_pos = x_ins_base + offsets[s_idx]
            
            for c_idx, category in enumerate(plot_categories):
                values = df_ins[category].values
                ax_ins.bar(x_ins_pos, values, width=bar_width, bottom=bottom_ins, color=colors[c_idx])
                bottom_ins += values
            global_max_inset_total = max(global_max_inset_total, df_ins["Total"].max())
            
        ax_ins.set_xticks(x_ins_base)
        ax_ins.set_xticklabels(inset_tools, rotation=45, ha="right", fontsize=8)
        ax_ins.spines['top'].set_visible(False)
        ax_ins.spines['right'].set_visible(False)

    y_lim = global_max_total * 1.1
    y_lim_ins = global_max_inset_total * 1.1
    for ax in axes:
        ax.set_ylim(0, y_lim)
        for child in ax.get_children():
            if isinstance(child, plt.Axes):
                child.set_ylim(0, y_lim_ins)

    fig.legend(
        labels=["unique tool", "2-4 tools", "5-7 tools", "8-10 tools", "11 tools"],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=5,
        frameon=False,
        fontsize=12
    )

    output_base_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap"
    os.makedirs(output_base_dir, exist_ok=True)
    anno_name = annotation.replace("/", "_")
    output_path = os.path.join(output_base_dir, f"grouped_aligners_{anno_name}_real_trim.pdf")
    
    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"Chart saved to {output_path}")
