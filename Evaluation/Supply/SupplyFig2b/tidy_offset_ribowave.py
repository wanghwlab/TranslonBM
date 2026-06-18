import glob
import os
import csv

# ================= 配置区域 =================

# 1. 输入文件的匹配模式 
INPUT_PATTERN = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribowave/ribowave_chrN/orf_pred_default/*/P-site/*.psite1nt.txt"

# 2. 输出 CSV 文件路径
OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/ribowave_converted_offsets.csv"

# ================= 主程序 =================
def process_ribowave_offsets():
    files = glob.glob(INPUT_PATTERN)
    
    if not files:
        print("未找到任何文件，请检查路径配置。")
        return

    print(f"找到 {len(files)} 个文件，开始处理...")

    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, delimiter=',') 
        writer.writerow(["sample", "ORFtools", "read_length", "offset"])

        for file_path in files:
            try:
                # --- 1. 获取样本名称 ---
                filename = os.path.basename(file_path)
                sample_name = filename.replace(".psite1nt.txt", "")
                
                # --- 2. 读取文件内容 ---
                with open(file_path, 'r') as infile:
                    for line in infile:
                        line = line.strip()
                        if not line: continue
                        
                        parts = line.split()
                        
                        if len(parts) >= 2:
                            read_len_str = parts[0]
                            psite_pos_str = parts[1]
                            
                            # --- 3. Offset 矫正 ---
                            try:
                                position = int(psite_pos_str)
                                offset = position - 1
                                
                                writer.writerow([sample_name, "ribowave", read_len_str, offset])
                                
                            except ValueError:
                                print(f"警告: 文件 {filename} 中存在非数字行: {line}")

            except Exception as e:
                print(f"处理文件 {file_path} 时出错: {e}")

    print(f"处理完成！结果已保存至: {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    process_ribowave_offsets()
