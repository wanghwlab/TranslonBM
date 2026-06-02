import pandas as pd
import os
import itertools
import re
from multiprocessing import Pool
import sys

# --- 配置区 (请根据您的环境修改) ---

WORKDIR = "/home/tangyuewen/ORF_benchmark/Maxquant_2026.1"
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_group/PR_recall_Maxquant/"

REF_BLOCK_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.block_Trans_Mod.txt"

SAMPLES = [ "SRX1254413", "SRX11812007_SRX11812008_SRX11812009", "SRX5887328_SRX5887329_SRX5887330"]
ALIGNERS = ["tophat2", "STAR", "hisat2"]
PRED_FLAGS = ["orf_pred_default"] 
TRIM_FLAGS = ["trim","untrim"]
FDR = [""]

NUM_PROCESSES = 16

ANNOTATED_CDS_SET = set()

# --- 辅助功能函数 ---
def normalize_coords_string(coords_str: str) -> str:
    try:
        blocks = coords_str.strip().split(',')
        block_tuples = []
        for b in blocks:
            start, end = map(int, b.split('-'))
            block_tuples.append((start, end))
        block_tuples.sort(key=lambda x: x[0])
        return ",".join([f"{s}-{e}" for s, e in block_tuples])
    except Exception:
        return coords_str

def load_annotated_cds(ref_path: str):
    global ANNOTATED_CDS_SET
    print(f"Loading annotated CDS blocks from: {ref_path} ...")
    try:
        with open(ref_path, 'r') as f:
            for line in f:
                if not line.strip(): continue
                fields = line.strip().split()
                if len(fields) >= 5:
                    chrom, raw_blocks, strand = fields[1].strip(), fields[2].strip(), fields[3].strip()
                    norm_blocks = normalize_coords_string(raw_blocks)
                    ANNOTATED_CDS_SET.add((chrom, strand, norm_blocks))
        print(f"Successfully loaded {len(ANNOTATED_CDS_SET)} annotated CDS records.")
    except Exception as e:
        print(f"[ERROR] Failed to load annotated CDS file: {e}")
        sys.exit(1)

def is_all_seq(protein_id_str):
    if not isinstance(protein_id_str, str): return False
    main_ids = protein_id_str.split('|')[0]
    return all(id_val.startswith('seq_') for id_val in main_ids.split(';'))

def has_any_unique_peptides(count_str):
    if not isinstance(count_str, str): return False
    try:
        return any(int(c) > 0 for c in count_str.split(';'))
    except (ValueError, TypeError):
        return False

def extract_seq_ids(id_string):
    if not isinstance(id_string, str): return []
    return re.findall(r'seq_\d+', id_string)

# --- 核心工作函数 ---
def process_level(params):
    sample, trim_flag, pred_flag, aligner, fdr, level = params
    print(f"开始处理 [{level.capitalize()} Level]: Sample={sample}, Aligner={aligner}")

    try:
        output_subdir = os.path.join(BASE_OUTPUT_DIR, trim_flag, pred_flag)
        stats_path = os.path.join(output_subdir, f"{sample}_MS_valid_{level}_{aligner}_score.tsv")
        
        ref_path = os.path.join(WORKDIR, "maxquant_ref", f"{pred_flag}_{trim_flag}", f"{sample}_ref.tsv")
        rerun_subdir = f"{pred_flag}_{trim_flag}"
        
        if level == 'peptide':
            ms_data_path = os.path.join(WORKDIR, "maxquant_rerun", rerun_subdir, fdr, aligner, sample, "combined", "txt", "peptides.txt")
            ms_df = pd.read_csv(ms_data_path, sep='\t', low_memory=False)
            valid_items = ms_df[(ms_df['Unique (Proteins)'] == 'yes') & (ms_df['Potential contaminant'].isna()) & (ms_df['Reverse'].isna())].copy()
            id_column = 'Proteins'
        else: # protein level
            ms_data_path = os.path.join(WORKDIR, "maxquant_rerun", rerun_subdir, fdr, aligner, sample, "combined", "txt", "proteinGroups.txt")
            ms_df = pd.read_csv(ms_data_path, sep='\t', low_memory=False)
            has_mask = ms_df['Peptide counts (unique)'].apply(has_any_unique_peptides)
            valid_items = ms_df[(ms_df['Potential contaminant'].isna()) & (ms_df['Reverse'].isna()) & (has_mask)].copy()
            valid_items = valid_items[valid_items['Protein IDs'].apply(is_all_seq)]
            id_column = 'Protein IDs'

        # 1. 加载完整的 reference 文件，并清理表头空白字符
        full_ref_df = pd.read_csv(ref_path, sep='\t')
        full_ref_df.rename(columns=lambda x: x.strip(), inplace=True)
        
        ref_df = full_ref_df[(full_ref_df['samples'] == sample) & (full_ref_df['aligner'] == aligner)].copy()
        
        if ref_df.empty:
            print(f"  -> [警告] 在 ref 文件中未找到 Sample={sample}, Aligner={aligner} 的数据。跳过。")
            return

        # ★★★ 核心修复 1：提取当前 aligner 的所有合法 ID 集合（宇宙集合） ★★★
        aligner_seq_ids = set(ref_df['sequence_id'].astype(str).str.strip())

        # 2. 划分 Novel ORF 集合
        novel_seq_ids = set()
        if 'chrom' in ref_df.columns and 'coordinate_0base' in ref_df.columns:
            for _, row in ref_df.iterrows():
                try:
                    chrom = str(row['chrom']).strip()
                    strand = str(row['strand']).strip()
                    raw_coords = str(row['coordinate_0base']).strip()
                    norm_blocks = normalize_coords_string(raw_coords)
                    
                    if (chrom, strand, norm_blocks) not in ANNOTATED_CDS_SET:
                        novel_seq_ids.add(str(row['sequence_id']).strip())
                except Exception:
                    continue
        else:
            print(f"  -> [错误] ref 文件中缺少坐标列。当前列名: {list(ref_df.columns)}")
            return
            
        print(f"  -> [统计] {sample} ({aligner}) 在 _ref.tsv 中包含的 Novel ORF 数量为: {len(novel_seq_ids)}")

        results = []
        
        # 3. 从质谱结果中提取所有被验证的 seq_id
        raw_validated_seq_ids = set(valid_items[id_column].apply(extract_seq_ids).explode().dropna())
        
        # ★★★ 核心修复 2：将质谱检测到的全局集合与当前 aligner 集合取交集 ★★★
        validated_seq_ids = raw_validated_seq_ids.intersection(aligner_seq_ids)
        ms_valid_orfs_num = len(validated_seq_ids)
        
        # 计算 MS_novel 集合 (MS集合 ∩ Novel集合)
        ms_novel_seq_ids = validated_seq_ids.intersection(novel_seq_ids)
        ms_novel_orfs_num = len(ms_novel_seq_ids)

        # 4. 按 orf_detector 分组进行计算
        for detector, group in ref_df.groupby('orf_detector'):
            predicted_seq_ids = set(group['sequence_id'].astype(str).str.strip().unique())
            pred_orfs_num = len(predicted_seq_ids)
            
            # 总体指标
            valid_orfs_num = len(predicted_seq_ids.intersection(validated_seq_ids))
            validation_rate = valid_orfs_num / pred_orfs_num if pred_orfs_num > 0 else 0
            recall = valid_orfs_num / ms_valid_orfs_num if ms_valid_orfs_num > 0 else 0
            f1 = (2 * validation_rate * recall) / (validation_rate + recall) if (validation_rate + recall) > 0 else 0
            
            # Novel 指标
            # ★★★ 必须先定义 pred_novel_seq_ids (预测出的序列与novel全集的交集) ★★★
            pred_novel_seq_ids = predicted_seq_ids.intersection(novel_seq_ids)
            # 1. 预测的 novel 总数 (计算 Validation Rate 的分母)
            pred_novel_orfs_num = len(pred_novel_seq_ids)

            # 2. 被质谱验证的 novel 数量 (分子)
            valid_novel_orfs_num = len(pred_novel_seq_ids.intersection(ms_novel_seq_ids))

            # 3. 计算 Novel Recall
            novel_recall = valid_novel_orfs_num / ms_novel_orfs_num if ms_novel_orfs_num > 0 else 0

            # 4. 计算 Novel Validation Rate (新增)
            novel_validation_rate = valid_novel_orfs_num / pred_novel_orfs_num if pred_novel_orfs_num > 0 else 0

            results.append({
                'aligner': aligner,
                'orf_detector': detector,
                'MS_valid_ORFs_num': ms_valid_orfs_num,
                'pred_ORFs_num': pred_orfs_num,
                'valid_ORFs_num': valid_orfs_num,
                'Validation_Rate': validation_rate,
                'Recall': recall,
                'F1_score': f1,
                'MS_novel_ORFs_num': ms_novel_orfs_num,
                'pred_novel_ORFs_num': pred_novel_orfs_num,         # 输出预测的 novel 数量
                'valid_novel_ORFs_num': valid_novel_orfs_num,
                'Novel_Validation_Rate': novel_validation_rate,     # 输出 Novel Validation Rate
                'Novel_Recall': novel_recall
            })

        output_df = pd.DataFrame(results)
        os.makedirs(os.path.dirname(stats_path), exist_ok=True)
        output_df.to_csv(stats_path, sep='\t', index=False)
        print(f"完成 [{level.capitalize()} Level]: Sample={sample}, Aligner={aligner}")

    except FileNotFoundError as e:
        print(f"[错误] 文件未找到 for {sample}/{aligner}: {e}. 跳过。")
    except Exception as e:
        print(f"[错误] 处理失败 for {sample}/{aligner}: {e}. 跳过。")

def process_peptide_wrapper(params):
    return process_level(params + ('peptide',))

def process_protein_wrapper(params):
    return process_level(params + ('protein',))

# --- 主程序入口 ---
if __name__ == '__main__':
    load_annotated_cds(REF_BLOCK_FILE)

    task_params = list(itertools.product(SAMPLES, TRIM_FLAGS, PRED_FLAGS, ALIGNERS, FDR))
    print(f"总共生成 {len(task_params) * 2} 个分析任务 (peptide/protein), 使用 {NUM_PROCESSES} 个核心并行处理。")

    print("\n--- 开始 Peptide Level 分析 ---")
    with Pool(processes=NUM_PROCESSES) as pool:
        pool.map(process_peptide_wrapper, task_params)
    
    print("\n--- 开始 Protein Level 分析 ---")
    with Pool(processes=NUM_PROCESSES) as pool:
        pool.map(process_protein_wrapper, task_params)

    print("\n所有分析任务已结束。")
