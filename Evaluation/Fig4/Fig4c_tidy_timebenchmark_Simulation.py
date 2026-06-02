import os
import pandas as pd
import re 

def extract_parameters(path_string, keywords):
    """
    一个辅助函数，用于从路径字符串中查找并连接存在的关键词。
    """
    found_keywords = [k for k in keywords if k in path_string]
    if found_keywords:
        return '_'.join(found_keywords)
    return '-'

def process_benchmark_files(root_directory, output_csv_path):
    """
    遍历指定路径下的所有 .txt 文件，提取信息，并合并成一个 CSV 文件。
    """
    all_data_rows = []

    print(f"开始在 '{root_directory}' 文件夹中搜索 .txt 文件...")

    for root, dirs, files in os.walk(root_directory):
        for filename in files:
            if filename.endswith(".txt"):
                
                if 'time_benchmarks' not in root.split(os.sep):
                    continue 

                full_path = os.path.join(root, filename)
                
                base_name = os.path.splitext(filename)[0]
                try:
                    parts = base_name.rsplit('_', 1)
                    sample = parts[0]
                    aligner = parts[1]
                except IndexError:
                    print(f"警告：文件名格式不符，已跳过：{filename}")
                    continue

                try:
                    data_df = pd.read_csv(full_path, sep='\t')
                    if data_df.empty:
                        print(f"警告：文件为空，已跳过：{full_path}")
                        continue
                    
                    data_df['sample'] = sample
                    data_df['aligner'] = aligner
                    data_df['doc_path'] = full_path

                    all_data_rows.append(data_df)
                except Exception as e:
                    print(f"错误：处理文件 {full_path} 时出错：{e}")

    if not all_data_rows:
        print("未找到任何在 'time_benchmarks' 子文件夹下的 .txt 文件。")
        return

    print("正在合并所有数据...")
    final_df = pd.concat(all_data_rows, ignore_index=True)


    print("正在从路径中提取 software, chr, 和 parameter 信息...")

    regex_pattern = r'(\w+)_(chr[NM])'
    extracted_info = final_df['doc_path'].str.extract(regex_pattern)
    
    final_df['software'] = extracted_info[0].fillna('-')
    final_df['chr'] = extracted_info[1].fillna('-')

    param_keywords = [
        'ribotaper', 'orfrater', 'rpbp', 'orfquant', 'ribowave', 'gedi', 
        'ribotricer', 'ribohmm', 'ribotish', 'ribocode', 'riborf'
    ]
    final_df['parameter'] = final_df['doc_path'].apply(extract_parameters, keywords=param_keywords)

    desired_order = [
        'sample', 'aligner', 'software', 'chr', 'parameter',
        's', 'h:m:s', 'max_rss', 'max_vms', 'max_uss', 
        'max_pss', 'io_in', 'io_out', 'mean_load', 'cpu_time', 'doc_path'
    ]
    
    for col in desired_order:
        if col not in final_df.columns:
            final_df[col] = None
            
    final_df = final_df[desired_order]

    final_df.to_csv(output_csv_path, index=False)
    print(f"\n处理完成！结果已成功保存到：{output_csv_path}")
    print(f"总共处理了 {len(final_df)} 个文件记录。")


if __name__ == '__main__':

    SEARCH_PATH = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu/"
    OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/time_benchmark/ORF_simu_time_benchmark_summary.csv"

    process_benchmark_files(SEARCH_PATH, OUTPUT_FILE)
