import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import os
import math


INPUT_CSV = "/home/tangyuewen/ORF_benchmark/final_ORFs_2025.7/plots/Other_calculate/aligner_read_overlap_with_publication.csv"

OUTPUT_PDF = "/home/tangyuewen/ORF_benchmark/final_ORFs_2025.7/plots/Other_calculate/aligner_read_overlap_venn.pdf"

SAMPLE_ORDER_AND_TITLES = {
    "SRX876063": "Ji et al. (2015)",
    "SRX740748": "Gao et al. (2015)",
    "SRX1254413": "Calviello et al. (2016)",
    "SRX1447296": "Raj et al. (2016)",
    "SRX5256543": "Martinez et al. (2020)",
    "SRX5887328": "Chen et al. (2020)",
    "SRX11812007": "Chothani et al. (2022)",
}

def format_label(count, total):
    """一个辅助函数，用于生成带千分位和百分比的标签。"""
    if total == 0:
        percentage = 0.0
    else:
        percentage = (count / total) * 100
    return f"{int(count):,}\n({percentage:.2f}%)"

def main():
    """主执行函数"""
    print(f"--- Starting Venn Diagram Plotting Script ---")
    
    # --- 阶段一: 读取并准备数据 ---
    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"[ERROR] Input CSV file not found at: {INPUT_CSV}. Exiting.")
        return

    # --- 阶段二: 创建绘图画布 ---
    sample_keys_to_plot = list(SAMPLE_ORDER_AND_TITLES.keys())
    num_samples = len(sample_keys_to_plot)
    ncols = 7
    nrows = math.ceil(num_samples / ncols)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 5.5), squeeze=False)
    axes = axes.flatten()

    print(f"Found {num_samples} samples to plot in a {nrows}x{ncols} grid.")

    # --- 阶段三: 按指定顺序循环绘制每个样本的韦恩图 ---
    for i, sample_prefix in enumerate(sample_keys_to_plot):
        ax = axes[i]
        title = SAMPLE_ORDER_AND_TITLES[sample_prefix]
        ax.set_title(title, fontsize=12, pad=15)
        
        sample_data = df[df['sample'].str.startswith(sample_prefix)]
        if sample_data.empty:
            print(f"  - Warning: No data found for sample prefix '{sample_prefix}'. Skipping subplot.")
            ax.text(0.5, 0.5, "No Data", ha='center', va='center', color='gray')
            ax.set_axis_off()
            continue

        sample_row = sample_data.iloc[0]
        
        subsets = (
            sample_row.get('STAR_only_count', 0),
            sample_row.get('hisat2_only_count', 0),
            sample_row.get('STAR_and_hisat2_only_count', 0),
            sample_row.get('tophat2_only_count', 0),
            sample_row.get('STAR_and_tophat2_only_count', 0),
            sample_row.get('hisat2_and_tophat2_only_count', 0),
            sample_row.get('all_three_aligners_count', 0),
        )
        
        total_reads = sample_row.get('total_unique_reads', sum(subsets))

        v = venn3(subsets=subsets, 
                  set_labels=('STAR', 'HISAT2', 'TopHat2'),
                  ax=ax)
        
        subset_map = {
            '100': 0, '010': 1, '110': 2, '001': 3,
            '101': 4, '011': 5, '111': 6
        }
        
        for region_id, data_index in subset_map.items():
            label = v.get_label_by_id(region_id)
            if label is not None:
                count = subsets[data_index]
                label.set_text(format_label(count, total_reads))

        colors = ['#008080', '#FF8C00', '#E9967A'] 
        if v.get_patch_by_id('100'): v.get_patch_by_id('100').set_color(colors[0])
        if v.get_patch_by_id('010'): v.get_patch_by_id('010').set_color(colors[1])
        if v.get_patch_by_id('001'): v.get_patch_by_id('001').set_color(colors[2])
        
        for p_id in ['110', '101', '011', '111']:
             if v.get_patch_by_id(p_id): v.get_patch_by_id(p_id).set_color('dimgrey')
        
        for p in v.patches:
            if p: p.set_alpha(0.6)

        venn3_circles(subsets=subsets, linestyle="-", linewidth=1, color="black", ax=ax)

    # --- 阶段四: 整理并保存PDF ---
    for i in range(num_samples, len(axes)):
        axes[i].set_axis_off()
        
    plt.tight_layout(pad=2.0)
    
    os.makedirs(os.path.dirname(OUTPUT_PDF), exist_ok=True)
    
    print(f"\n--- Saving plot to: {OUTPUT_PDF} ---")
    plt.savefig(OUTPUT_PDF, format='pdf', bbox_inches='tight')
    plt.close(fig)
    
    print("--- Script Finished Successfully ---")

if __name__ == "__main__":
    main()
