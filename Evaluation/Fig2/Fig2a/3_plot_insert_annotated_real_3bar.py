#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# 1. 更新后的名称映射（包含仿真数据和 6 套 SRX 实测数据）
name_mapping = {
    #"simulation_6M_T1": "simulation_6M_T1",
    #"simulation_6M_T3": "simulation_6M_T3",
    #"simulation_60M_T1": "simulation_60M_T1",
    #"simulation_60M_T3": "simulation_60M_T3",
    "SRX876063": "Ji et al. (2015)",
    "SRX740748": "Gao et al.(2015)",
    "SRX1254413": "Calviello et al.(2016)",
    "SRX5256543": "Martinez et al.(2020)",
    "SRX5887328": "Chen et al.(2020)",
    "SRX11812007": "Chothani et al.(2022)",
}

# 颜色方案
color_map_name = 'Set3'
cmap = plt.get_cmap(color_map_name)
#colors = [plt.colors.to_hex(cmap(i)) for i in [4, 5, 6, 2, 3]]
colors = [mpl.colors.to_hex(cmap(i)) for i in [4, 5, 6, 2, 3]]

base_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/tools_overlap/"
#annotation_types = ['merged_ATG/orf_pred_default_trim', 'merged_NTG/orf_pred_default_trim']
annotation_types = ['merged_canonical_ATG/orf_pred_default_untrim', 'merged_non_canonical_ATG/orf_pred_default_untrim','merged_canonical_ATG/orf_pred_default_trim', 'merged_non_canonical_ATG/orf_pred_default_trim']
softwares_order = ['hisat2', 'STAR', 'tophat2']

def process_dataframe(df, annotation_type):
    if 'NTG' in annotation_type and "predicted by 11 soft ORF num" in df.columns:
        df["predicted by 11 soft ORF num"] = df["predicted by 11 soft ORF num"] - 1

    numeric_columns = [f"predicted by {i} soft ORF num" for i in range(1, 12)]
    for col in numeric_columns:
        if col not in df.columns:
            df[col] = 0
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors="coerce")

    df["2-4 tools"] = df["predicted by 2 soft ORF num"] + df["predicted by 3 soft ORF num"] + df["predicted by 4 soft ORF num"]
    df["5-7 tools"] = df["predicted by 5 soft ORF num"] + df["predicted by 6 soft ORF num"] + df["predicted by 7 soft ORF num"]
    df["8-10 tools"] = df["predicted by 8 soft ORF num"] + df["predicted by 9 soft ORF num"] + df["predicted by 10 soft ORF num"]
    df["11 tools"] = df["predicted by 11 soft ORF num"]
    df.rename(columns={"predicted by 1 soft ORF num": "1 tool"}, inplace=True)
    
    plot_categories = ["1 tool", "2-4 tools", "5-7 tools", "8-10 tools", "11 tools"]
    df["Total"] = df[plot_categories].sum(axis=1)
    df['tool_label'] = df.iloc[:, 0].str.split('_').str[-2]
    
    return df, plot_categories

for annotation in annotation_types:
    all_files = glob.glob(f"{base_dir}/{annotation}/*_overlap_statistics.csv")
    present_keys = set()
    for f in all_files:
        fname = os.path.basename(f)
        for key in name_mapping.keys():
            if fname.startswith(key + "_"):
                present_keys.add(key)
    
    dataset_ids = sorted(list(present_keys), key=lambda x: list(name_mapping.keys()).index(x))
    
    if not dataset_ids:
        continue

    global_max_total = 0
    global_max_inset_total = 0

    fig, axes = plt.subplots(1, len(dataset_ids), figsize=(5 * len(dataset_ids), 6))
    if len(dataset_ids) == 1: axes = [axes]


    bar_width = 0.22     
    gap = 0.06            
    offsets = [-(bar_width + gap), 0, (bar_width + gap)] 

    for idx, (d_id, ax) in enumerate(zip(dataset_ids, axes)):
        dataset_display_name = name_mapping[d_id]
        combined_data = {}
        all_tools = []
        

        for soft in softwares_order:
            f_pattern = f"{base_dir}/{annotation}/{d_id}*_{soft}_overlap_statistics.csv"
            f_list = glob.glob(f_pattern)
            if f_list:
                df, plot_categories = process_dataframe(pd.read_csv(f_list[0]), annotation)
                if soft == 'hisat2': 
                    df = df.sort_values(by="Total", ascending=False)
                    all_tools = df['tool_label'].tolist()
                combined_data[soft] = df
        
        if not combined_data: continue

        x_base = np.arange(len(all_tools))
        
        for s_idx, soft in enumerate(softwares_order):
            if soft not in combined_data: continue
            df_plot = combined_data[soft].set_index('tool_label').reindex(all_tools).fillna(0)
            
            bottom = np.zeros(len(all_tools))
            x_pos = x_base + offsets[s_idx]
            
            for c_idx, category in enumerate(plot_categories):
                values = df_plot[category].values
                ax.bar(x_pos, values, width=bar_width, bottom=bottom, color=colors[c_idx], align='center')
                bottom += values
            global_max_total = max(global_max_total, df_plot["Total"].max())


        ax.set_title(dataset_display_name, fontsize=12, fontweight="bold")
        ax.set_xticks(x_base)
        ax.set_xticklabels(all_tools, rotation=45, ha="right", fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if idx == 0: ax.set_ylabel("Total number", fontsize=12)


        if len(all_tools) >= 4:
            inset_tools = all_tools[-4:]
            ax_ins = ax.inset_axes([0.5, 0.52, 0.48, 0.45])
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
            ax_ins.set_xticklabels(inset_tools, rotation=45, ha="right", fontsize=7)
            ax_ins.tick_params(axis='y', labelsize=7)
            ax_ins.spines['top'].set_visible(False)
            ax_ins.spines['right'].set_visible(False)


    y_lim_main = global_max_total * 1.1
    y_lim_ins = global_max_inset_total * 1.1
    for ax in axes:
        ax.set_ylim(0, y_lim_main)
        for child in ax.get_children():
            if isinstance(child, plt.Axes):
                child.set_ylim(0, y_lim_ins)


    fig.legend(
        labels=["unique tool", "2-4 tools", "5-7 tools", "8-10 tools", "all (11 tools)"],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=5,
        frameon=False,
        fontsize=11
    )


    output_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap"
    os.makedirs(output_dir, exist_ok=True)
    anno_tag = annotation.replace("/", "_")
    output_path = os.path.join(output_dir, f"grouped_aligners_all_datasets_{anno_tag}.pdf")
    
    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"成功保存至: {output_path}")
