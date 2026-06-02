import pandas as pd
import os
import glob
import itertools
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
import matplotlib

matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# --- 配置区 ---
INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_canonical_ATG/orf_pred_default_untrim/"
OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/Jaccard_Heatmaps/merged_canonical_ATG/orf_pred_default_untrim/"

# 数据集顺序配置
MANUAL_ORDER_REAL = [
    "SRX876063", "SRX740748", "SRX1254413",
    "SRX5256543", "SRX5887328", "SRX11812007",
]
MANUAL_ORDER_REPLICATE = [
    "SRX7666669-73", "SRX7666674-78", "SRX7666679-83",
    "SRX7666684-88", "SRX7666689-93", "SRX7666694-98",
]
MANUAL_ORDER_SIMULATION = [
    "simulation_6M_T1", "simulation_6M_T3",
    "simulation_60M_T1", "simulation_60M_T3",
]

COLORMAP = "YlGnBu"
# 工具顺序
MANUAL_TOOL_ORDER = ['ribohmm', 'rpbp', 'orfrater', 'gedi', 'ORFquant', 'ribotricer', 'ribotish', 'ribotaper', 'ribocode', 'riborf', 'ribowave']

NAME_MAPPING = {
    "SRX876063": "Ji et al. (2015)",
    "SRX740748": "Gao et al.(2015)",
    "SRX1254413": "Calviello et al.(2016)",
    "SRX1447296": "Raj et al.(2016)",
    "SRX5256543": "Martinez et al.(2020)",
    "SRX5887328": "Chen et al.(2020)",
    "SRX11812007": "Chothani et al.(2022)",
}

# --- 核心函数 ---

def parse_filename(filepath):
    """从文件名解析出 sample, aligner, tool。"""
    base = os.path.basename(filepath).replace('_gcoor.tsv.gz', '')
    parts = base.split('_')
    if len(parts) >= 3:
        tool = parts[-1]; aligner = parts[-2]; sample = '_'.join(parts[:-2])
        return sample, aligner, tool
    return None, None, None

def calculate_jaccard_matrix(orf_sets):
    """计算Jaccard矩阵。"""
    tools = [tool for tool in MANUAL_TOOL_ORDER if tool in orf_sets]
    if len(tools) < 2: return None
    matrix = pd.DataFrame(1.0, index=tools, columns=tools)
    for t1, t2 in itertools.combinations(tools, 2):
        set1 = orf_sets.get(t1, set())
        set2 = orf_sets.get(t2, set())
        
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        
        jaccard = intersection / union if union > 0 else 0
        matrix.loc[t1, t2] = matrix.loc[t2, t1] = jaccard
    return matrix

def get_display_name(sample_name, name_mapping):
    for key, study_name in name_mapping.items():
        if sample_name.startswith(key):
            return study_name
    return sample_name

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("--- 正在预加载所有ORF数据... ---")
    all_orf_data = {}
    all_files = glob.glob(os.path.join(INPUT_DIR, "*_gcoor.tsv.gz"))
    
    COLS_TO_USE = [
        'coordinate_id', 'chrom', 'coordinate_0base', 'strand', 
        'ORF_sequence_correct', 'start_codon', 'ORF_sequence_aa'
    ]
    FULL_COLUMN_NAMES = [
        'gene_id', 'coordinate_id', 'chrom',
        'coordinate_0base', 'strand', 'ORF_sequence_correct',
        'start_codon', 'ORF_sequence_aa'
    ]

    for f in all_files:
        sample, aligner, tool = parse_filename(f)
        if not all([sample, aligner, tool]): continue
        
        try:
            df = pd.read_csv(
                f, 
                sep='\t', 
                compression='gzip', 
                low_memory=False,
                
                header=None,
                names=FULL_COLUMN_NAMES,
                
                usecols=lambda c: c in COLS_TO_USE 
            )
            
            df = df[df['start_codon'] == 'ATG']
            
            if df.empty:
                all_orf_data[(sample, aligner, tool)] = set()
                continue
            
            orf_set = {tuple(row) for row in df.itertuples(index=False, name=None)}
            
            all_orf_data[(sample, aligner, tool)] = orf_set
            
        except ValueError as ve:
            try:
                df = pd.read_csv(f, sep='\t', compression='gzip', low_memory=False)
                df = df[df['start_codon'] == 'ATG']
                df = df.drop(columns=['gene_id', 'transcript_id'], errors='ignore')
                orf_set = {tuple(row) for row in df.itertuples(index=False, name=None)}
                all_orf_data[(sample, aligner, tool)] = orf_set
            except Exception as e2:
                print(f"  [错误] 读取文件 {f} 失败: {e2}")
        except Exception as e:
            print(f"  [错误] 读取文件 {f} 失败: {e}")

    all_aligners = sorted(list(set(k[1] for k in all_orf_data.keys())))
    all_samples_found = sorted(list(set(k[0] for k in all_orf_data.keys())))

    for aligner in all_aligners:
        dataset_config = {
            'Simulation': MANUAL_ORDER_SIMULATION,
            'Real': MANUAL_ORDER_REAL,
            'Replicate': MANUAL_ORDER_REPLICATE
        }
        
        for group_name, manual_order_list in dataset_config.items():
            samples_to_plot = []
            for srx_prefix in manual_order_list:
                for full_sample_name in all_samples_found:
                    if full_sample_name.startswith(srx_prefix):
                        samples_to_plot.append(full_sample_name)
                        break
            
            if not samples_to_plot: continue
            
            print(f"\n--- 正在处理: Aligner='{aligner}', Dataset='{group_name}' ---")
            
            num_samples = len(samples_to_plot)
            ncols = 4
            nrows = math.ceil(num_samples / ncols)
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 7, nrows * 6), squeeze=False)
            axes = axes.flatten()
            fig.suptitle(f'ORF Tools Overlap (Jaccard) | Aligner: {aligner} | Dataset: {group_name}', fontsize=20, y=1.02)

            for i, sample in enumerate(samples_to_plot):
                ax = axes[i]
                display_name = get_display_name(sample, NAME_MAPPING)
                ax.set_title(display_name, fontsize=14, pad=10)
                
                sample_orf_sets = {k[2]: v for k, v in all_orf_data.items() if k[0] == sample and k[1] == aligner}
                jaccard_matrix = calculate_jaccard_matrix(sample_orf_sets)

                if jaccard_matrix is None or jaccard_matrix.empty:
                    ax.text(0.5, 0.5, 'Not enough tools', ha='center', va='center', color='gray')
                    ax.set_axis_off()
                    continue

                mask = np.triu(np.ones_like(jaccard_matrix, dtype=bool))
                sns.set_theme(style="white")
                sns.heatmap(jaccard_matrix, mask=mask, annot=True, fmt=".2f", cmap=COLORMAP, linewidths=.5, vmin=0, vmax=1, ax=ax, cbar_kws={"shrink": .8})
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                plt.setp(ax.get_yticklabels(), rotation=0)

            for i in range(num_samples, len(axes)):
                axes[i].set_axis_off()

            plt.tight_layout(rect=[0, 0, 1, 0.96])
            output_filename = os.path.join(OUTPUT_DIR, f"jaccard_{aligner}_{group_name}.pdf")
            print(f"--- 保存整合图: {output_filename} ---")
            plt.savefig(output_filename, format='pdf', bbox_inches='tight')
            plt.close(fig)

    print("\n所有任务已完成！")

if __name__ == '__main__':
    main()
