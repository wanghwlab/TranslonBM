import pandas as pd

# 1. 定义路径
input_file = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/merged_RLD_stats_final.csv'
output_file_1 = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/proportions_by_sample_tool.csv'
output_file_2 = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/proportions_by_sample_only.csv'

def process_data(file_path):
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        return

    df_filtered = df[(df['read_length'] >= 25) & (df['read_length'] <= 35)].copy()

    # --- 任务 1: 按 sample 和 ORFtools 分组 ---
    df_group1 = df_filtered.copy()
    group1_sums = df_group1.groupby(['sample', 'ORFtools'])['read_count'].transform('sum')
    df_group1['Proportion (%)'] = (df_group1['read_count'] / group1_sums * 100).round(2).astype(str) + '%'
    
    result1 = df_group1.rename(columns={'read_length': 'Reads Length (nt)'})
    result1 = result1.sort_values(['sample', 'ORFtools', 'Reads Length (nt)'])
    result1 = result1[['sample', 'ORFtools', 'Reads Length (nt)', 'Proportion (%)']]
    result1.to_csv(output_file_1, index=False)
    print(f"文件 1 已生成 (按 Sample + Tool 分组): {output_file_1}")

    # --- 任务 2: 仅按 sample 分组 ---
    df_group2 = df_filtered.groupby(['sample', 'read_length'])['read_count'].sum().reset_index()
    group2_sums = df_group2.groupby(['sample'])['read_count'].transform('sum')
    df_group2['Proportion (%)'] = (df_group2['read_count'] / group2_sums * 100).round(2).astype(str) + '%'
    
    result2 = df_group2.rename(columns={'read_length': 'Reads Length (nt)'})
    result2 = result2.sort_values(['sample', 'Reads Length (nt)'])
    result2 = result2[['sample', 'Reads Length (nt)', 'Proportion (%)']]
    result2.to_csv(output_file_2, index=False)
    print(f"文件 2 已生成 (仅按 Sample 分组): {output_file_2}")

if __name__ == "__main__":
    process_data(input_file)
