#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pandas as pd
import os
import glob

# Step 1: 设置输入和输出目录
input_dir = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_canonical_ATG/orf_pred_default_trim'
output_dir = '/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/tools_overlap/merged_canonical_ATG/orf_pred_default_trim'

# Step 2: 定义需要处理的文件模式列表
file_patterns = [
{"prefix": "SRX11812007_SRX11812008_SRX11812009_hisat2", "pattern": "SRX11812007_SRX11812008_SRX11812009_hisat2*gcoor.tsv.gz"},
{"prefix": "SRX11812007_SRX11812008_SRX11812009_STAR", "pattern": "SRX11812007_SRX11812008_SRX11812009_STAR*gcoor.tsv.gz"},
{"prefix": "SRX11812007_SRX11812008_SRX11812009_tophat2", "pattern": "SRX11812007_SRX11812008_SRX11812009_tophat2*gcoor.tsv.gz"},
{"prefix": "SRX1254413_hisat2", "pattern": "SRX1254413_hisat2*gcoor.tsv.gz"},
{"prefix": "SRX1254413_STAR", "pattern": "SRX1254413_STAR*gcoor.tsv.gz"},
{"prefix": "SRX1254413_tophat2", "pattern": "SRX1254413_tophat2*gcoor.tsv.gz"},
{"prefix": "SRX5256543_SRX5256555_hisat2", "pattern": "SRX5256543_SRX5256555_hisat2*gcoor.tsv.gz"},
{"prefix": "SRX5256543_SRX5256555_STAR", "pattern": "SRX5256543_SRX5256555_STAR*gcoor.tsv.gz"},
{"prefix": "SRX5256543_SRX5256555_tophat2", "pattern": "SRX5256543_SRX5256555_tophat2*gcoor.tsv.gz"},
{"prefix": "SRX5887328_SRX5887329_SRX5887330_hisat2", "pattern": "SRX5887328_SRX5887329_SRX5887330_hisat2*gcoor.tsv.gz"},
{"prefix": "SRX5887328_SRX5887329_SRX5887330_STAR", "pattern": "SRX5887328_SRX5887329_SRX5887330_STAR*gcoor.tsv.gz"},
{"prefix": "SRX5887328_SRX5887329_SRX5887330_tophat2", "pattern": "SRX5887328_SRX5887329_SRX5887330_tophat2*gcoor.tsv.gz"},
{"prefix": "SRX740748_hisat2", "pattern": "SRX740748_hisat2*gcoor.tsv.gz"},
{"prefix": "SRX740748_STAR", "pattern": "SRX740748_STAR*gcoor.tsv.gz"},
{"prefix": "SRX740748_tophat2", "pattern": "SRX740748_tophat2*gcoor.tsv.gz"},
{"prefix": "SRX876063_SRX876069_hisat2", "pattern": "SRX876063_SRX876069_hisat2*gcoor.tsv.gz"},
{"prefix": "SRX876063_SRX876069_STAR", "pattern": "SRX876063_SRX876069_STAR*gcoor.tsv.gz"},
{"prefix": "SRX876063_SRX876069_tophat2", "pattern": "SRX876063_SRX876069_tophat2*gcoor.tsv.gz"},
{"prefix": "simulation_60M_T1_hisat2", "pattern": "simulation_60M_T1_hisat2*gcoor.tsv.gz"},
{"prefix": "simulation_60M_T1_STAR", "pattern": "simulation_60M_T1_STAR*gcoor.tsv.gz"},
{"prefix": "simulation_60M_T1_tophat2", "pattern": "simulation_60M_T1_tophat2*gcoor.tsv.gz"},
{"prefix": "simulation_60M_T3_hisat2", "pattern": "simulation_60M_T3_hisat2*gcoor.tsv.gz"},
{"prefix": "simulation_60M_T3_STAR", "pattern": "simulation_60M_T3_STAR*gcoor.tsv.gz"},
{"prefix": "simulation_60M_T3_tophat2", "pattern": "simulation_60M_T3_tophat2*gcoor.tsv.gz"},
{"prefix": "simulation_6M_T1_hisat2", "pattern": "simulation_6M_T1_hisat2*gcoor.tsv.gz"},
{"prefix": "simulation_6M_T1_STAR", "pattern": "simulation_6M_T1_STAR*gcoor.tsv.gz"},
{"prefix": "simulation_6M_T1_tophat2", "pattern": "simulation_6M_T1_tophat2*gcoor.tsv.gz"},
{"prefix": "simulation_6M_T3_hisat2", "pattern": "simulation_6M_T3_hisat2*gcoor.tsv.gz"},
{"prefix": "simulation_6M_T3_STAR", "pattern": "simulation_6M_T3_STAR*gcoor.tsv.gz"},
{"prefix": "simulation_6M_T3_tophat2", "pattern": "simulation_6M_T3_tophat2*gcoor.tsv.gz"}

]

column_names = ['gene_id', 'coordinate_id', 'chrom', 'coordinate_0base', 'strand', 
                'ORF_sequence_correct', 'start_codon', 'ORF_sequence_aa']

# Step 3: 遍历文件模式进行处理
for file_config in file_patterns:
    prefix = file_config["prefix"]
    file_pattern = file_config["pattern"]
    files = glob.glob(os.path.join(input_dir, file_pattern))

    if not files:
        print(f"No files found for pattern: {file_pattern}")
        continue

    print(f"Processing files for prefix: {prefix}")

    # 保存每个文件的唯一前四列数据（去重列：coordinate_id, chrom, coordinate_0base, strand）
    unique_dataframes = []

    for file in files:
        # Step 4: 读取文件，读取需要的列：coordinate_id, chrom, coordinate_0base, strand, start_codon
        df = pd.read_csv(file, sep='\t', compression='gzip', header=None, names=column_names,
                         usecols=['coordinate_id', 'chrom', 'coordinate_0base', 'strand', 'start_codon'])
        df_unique = df[['coordinate_id', 'chrom', 'coordinate_0base', 'strand']].drop_duplicates()
        unique_dataframes.append(df_unique)

    # Step 5: 合并所有文件四列并直接统计重复次数
    combined_df = pd.concat(unique_dataframes)
    overlap_count = combined_df.value_counts().reset_index(name='count')
    overlap_count.columns = ['coordinate_id', 'chrom', 'coordinate_0base', 'strand', 'count']

    overlap_count['tools'] = overlap_count['count'].apply(lambda x: f"{x}_tools")

    os.makedirs(output_dir, exist_ok=True)
    overlap_file = os.path.join(output_dir, f'{prefix}_overlap_count.txt')
    overlap_count.to_csv(overlap_file, sep='\t', index=False)

    # Step 6: 对每个文件统计重复情况
    result = {}
    tools_range = range(1, overlap_count['count'].max() + 1)

    for file, df_unique in zip(files, unique_dataframes):
        merged_df = pd.merge(df_unique, overlap_count, 
                             on=['coordinate_id', 'chrom', 'coordinate_0base', 'strand'], how='left')
        tools_count = merged_df['count'].value_counts().to_dict()

        result[file] = {f"predicted by {tools} soft ORF num": tools_count.get(tools, 0) for tools in tools_range}

    final_df = pd.DataFrame(result).T.fillna(0).astype(int)

    statistics_file = os.path.join(output_dir, f'{prefix}_overlap_statistics.csv')
    final_df.to_csv(statistics_file)

    print(f"Files processed for prefix: {prefix}. Outputs:\n- {overlap_file}\n- {statistics_file}")

print("All patterns processed successfully.")
