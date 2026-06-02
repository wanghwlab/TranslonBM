import glob
import os
import csv

# ================= 配置区域 =================

# 1. 输入文件的通配符路径
INPUT_PATTERN = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribocode/ribocode_chrN/P_site_determination/SR*/SRX*_pre_config.txt"

# 2. 输出文件名
OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/ribocode_extracted_offsets.csv"

# ================= 主程序 =================

def extract_ribocode_offsets():
    # 查找所有匹配的文件
    files = glob.glob(INPUT_PATTERN)
    
    if not files:
        print("未找到任何文件，请检查路径配置。")
        return

    print(f"找到 {len(files)} 个文件，开始提取...")

    # 打开 CSV 文件准备写入
    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile) # 默认逗号分隔
        
        # 写入表头
        writer.writerow(["sample", "ORFtools", "read_length", "offset"])

        for file_path in files:
            try:
                # 读取文件所有行
                with open(file_path, 'r') as infile:
                    lines = infile.readlines()
                
                # 寻找非注释行（不以 # 开头且非空的行）
                config_line = None
                for line in lines:
                    stripped_line = line.strip()
                    if stripped_line and not stripped_line.startswith('#'):
                        config_line = stripped_line
                        break
                
                # 如果找到了配置行，进行解析
                if config_line:
                    # 使用空白字符分割 (Tab 或 空格)
                    parts = config_line.split()
                    
                    if len(parts) >= 5:
                        # parts[0] -> SampleName (例如 SRX876063_SRX876069_tophat2)
                        # parts[3] -> ReadLengths (例如 25,28,29)
                        # parts[4] -> Offsets (例如 9,12,12)
                        
                        sample_name = parts[0]
                        lengths_str = parts[3]
                        offsets_str = parts[4]
                        
                        # 将逗号分隔的字符串转为列表
                        lengths = lengths_str.split(',')
                        offsets = offsets_str.split(',')
                        
                        # 确保长度和offset数量一致，然后成对写入
                        if len(lengths) == len(offsets):
                            for l, o in zip(lengths, offsets):
                                # 写入一行: sample, ribocode, length, offset
                                writer.writerow([sample_name, "ribocode", l, o])
                        else:
                            print(f"警告: {sample_name} 的长度与Offset数量不匹配")
                    else:
                        print(f"警告: 文件格式异常 (列数不足): {file_path}")
                else:
                    print(f"警告: 未在文件中找到有效配置行: {file_path}")

            except Exception as e:
                print(f"处理文件出错 {file_path}: {e}")

    print(f"处理完成！结果已保存至: {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    extract_ribocode_offsets()
