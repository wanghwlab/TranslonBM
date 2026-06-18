import pandas as pd
import glob
import os
import re

def parse_ribotricer_offsets():
    # 1. 定义 Ribotricer 文件的搜索路径
    input_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribotricer/ribotricer_chrN/orf_pred_default/*/*_psite_offsets.txt"
    
    # 2. 定义基准 Offset (Base Offset)
    DEFAULT_BASE_OFFSET = 12 

    files = glob.glob(input_pattern)
    print(f"找到 {len(files)} 个 Ribotricer offset 文件，开始转换...")

    all_data = []

    for file_path in files:
        try:
            # --- 提取元数据 (Sample 和 ORFtools) ---
            filename = os.path.basename(file_path)
            
            # 去掉后缀提取 Sample 名
            sample_name = filename.replace("_psite_offsets.txt", "")
            orf_tool = "ribotricer"

            # --- 解析文件内容 ---
            base_length = None
            offsets_dict = {}

            with open(file_path, 'r') as f:
                lines = f.readlines()

            for line in lines:
                line = line.strip()
                
                # 1. 提取基准长度 (relative lag to base: 29)
                base_match = re.search(r"relative lag to base:\s*(\d+)", line)
                if base_match:
                    base_length = int(base_match.group(1))
                    offsets_dict[base_length] = DEFAULT_BASE_OFFSET
                    continue

                # 2. 提取各长度的 lag (lag of 28: 0)
                lag_match = re.search(r"lag of\s*(\d+):\s*(-?\d+)", line)
                if lag_match:
                    read_len = int(lag_match.group(1))
                    lag = int(lag_match.group(2))
                    
                    # 核心计算公式
                    abs_offset = DEFAULT_BASE_OFFSET + lag
                    offsets_dict[read_len] = abs_offset

            # --- 转换为 DataFrame ---
            if offsets_dict:
                for r_len, off_val in offsets_dict.items():
                    all_data.append({
                        'sample': sample_name,
                        'ORFtools': orf_tool,
                        'read_length': r_len,
                        'offset': off_val
                    })
            else:
                print(f"警告: 文件 {filename} 中未解析到有效数据")

        except Exception as e:
            print(f"处理文件 {file_path} 出错: {e}")

    # --- 保存结果 ---
    if all_data:
        df = pd.DataFrame(all_data)
        
        # 排序：按样本和长度排序
        df = df.sort_values(by=['sample', 'read_length'])
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/ribotricer_converted_offsets.csv"
        df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print(f"转换完成！")
        print(f"保存文件: {os.path.abspath(output_file)}")
        print("\n数据预览 (前5行):")
        print(df.head())
        
        print("\n[验证] 检查 simulation_6M_T3 样本的 26nt offset (应为9):")
        check = df[(df['sample'].str.contains("simulation_6M_T3")) & (df['read_length'] == 26)]
        print(check)
    else:
        print("没有提取到数据。")

if __name__ == "__main__":
    parse_ribotricer_offsets()
