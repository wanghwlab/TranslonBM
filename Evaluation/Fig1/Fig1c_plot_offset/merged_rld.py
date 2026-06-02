import pandas as pd
import glob
import os

def extract_rld_from_stats():
    # 1. 设置文件路径
    input_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping/merge_chrN/*.stats"
    
    print(f"正在搜索文件: {input_pattern} ...")
    files = glob.glob(input_pattern)
    
    if not files:
        print("错误：未找到任何 .stats 文件，请检查路径。")
        return

    print(f"共找到 {len(files)} 个文件，开始提取 'RL' 行...")
    
    all_data = []

    for file_path in files:
        try:
            # --- A. 解析文件名 ---
            basename = os.path.basename(file_path)
            filename_no_ext = basename.replace('.stats', '')
            
            if '_' in filename_no_ext:
                sample_name, tool_name = filename_no_ext.rsplit('_', 1)
            else:
                sample_name = filename_no_ext
                tool_name = "Unknown"

            # --- B. 提取内容 (核心修改部分) ---
            file_data = []
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('RL'):
                        parts = line.strip().split()
                        
                        if len(parts) >= 3:
                            try:
                                r_len = int(parts[1])
                                r_count = int(parts[2])
                                
                                file_data.append({
                                    'sample': sample_name,
                                    'ORFtools': tool_name, 
                                    'read_length': r_len,
                                    'read_count': r_count
                                })
                            except ValueError:
                                continue 

            if file_data:
                all_data.extend(file_data)
            else:
                print(f"警告: 文件 {basename} 中没有找到 'RL' 开头的行。")
            
        except Exception as e:
            print(f"处理文件 {basename} 时出错: {e}")

    # --- C. 保存结果 ---
    if all_data:
        final_df = pd.DataFrame(all_data)
        
        cols = ['sample', 'ORFtools', 'read_length', 'read_count']
        final_df = final_df[cols]
        
        final_df = final_df.sort_values(by=['sample', 'ORFtools', 'read_length'])
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/untrim/merged_RLD_stats_final.csv"
        final_df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print("提取完成！")
        print(f"总数据行数: {len(final_df)}")
        print(f"已保存为: {os.path.abspath(output_file)}")
        print("\n数据预览:")
        print(final_df.head())
    else:
        print("未提取到任何数据。")

if __name__ == "__main__":
    extract_rld_from_stats()
