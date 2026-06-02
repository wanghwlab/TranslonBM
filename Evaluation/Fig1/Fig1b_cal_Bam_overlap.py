import os
import glob
import pysam 
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool

INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/Mapping/merge_chrN/"

OUTPUT_CSV = "/home/tangyuewen/ORF_benchmark/final_ORFs_2025.7/plots/Other_calculate/aligner_read_overlap_with_publication.csv"

NUM_PROCESSES = 16

NAME_MAPPING = {
    "SRX876063": "Ji et al. (2015)",
    "SRX740748": "Gao et al.(2015)",
    "SRX1254413": "Calviello et al.(2016)",
    "SRX1447296": "Raj et al.(2016)",
    "SRX5256543": "Martinez et al.(2020)",
    "SRX5887328": "Chen et al.(2020)",
    "SRX11812007": "Chothani et al.(2022)",
}


def parse_bam_filename(filepath: str) -> tuple:
    """
    一个稳健的函数，用于从BAM文件名中解析出样本名和比对软件名。
    """
    known_aligners = ['tophat2', 'hisat2', 'STAR']
    base = os.path.basename(filepath).replace('.bam', '')
    
    for aligner in known_aligners:
        if base.endswith('_' + aligner):
            sample = base[:-len('_' + aligner)]
            return sample, aligner
            
    return None, None

def extract_read_ids_from_bam(bam_path: str) -> tuple:
    """
    (多进程工作函数)
    从单个BAM文件中读取所有read的ID (query_name)，并返回一个集合。
    """
    print(f"  - Reading: {os.path.basename(bam_path)}...")
    read_ids = set()
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            for read in bamfile:
                read_ids.add(read.query_name)
        print(f"  - Finished: {os.path.basename(bam_path)}, found {len(read_ids)} unique reads.")
        return (bam_path, read_ids)
    except Exception as e:
        print(f"[ERROR] Failed to process {os.path.basename(bam_path)}. Reason: {e}")
        return (bam_path, set())

def main():
    """主执行函数"""
    print("--- Starting Read Overlap Analysis Script ---")

    # --- 阶段一: 查找并按样本分组所有BAM文件 ---
    print("\n--- Phase 1: Finding and grouping BAM files by sample ---")
    
    search_path = os.path.join(INPUT_DIR, "SRX74*.bam")
    all_bam_files = glob.glob(search_path)
    
    grouped_files = defaultdict(dict)
    for f in all_bam_files:
        sample, aligner = parse_bam_filename(f)
        if sample and aligner:
            grouped_files[sample][aligner] = f
            
    if not grouped_files:
        print(f"[ERROR] No valid BAM files matching 'sample_aligner.bam' format found in {INPUT_DIR}.")
        return

    print(f"Found {len(grouped_files)} unique samples to analyze.")

    # --- 阶段二: 并行提取所有文件的Read ID ---
    print(f"\n--- Phase 2: Extracting read IDs in parallel using {NUM_PROCESSES} cores... ---")
    
    files_to_process = [path for sample_data in grouped_files.values() for path in sample_data.values()]
    
    read_sets = {}
    with Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(extract_read_ids_from_bam, files_to_process)
        for path, ids in results:
            if path:
                read_sets[path] = ids
    
    print("\nRead ID extraction complete.")

    # --- 阶段三: 计算重叠并生成结果 ---
    print("\n--- Phase 3: Calculating overlaps for each sample ---")
    
    analysis_results = []
    
    for sample, aligner_paths in grouped_files.items():
        print(f"  - Calculating for sample: {sample}")
        
        star_reads = read_sets.get(aligner_paths.get('STAR'), set())
        hisat_reads = read_sets.get(aligner_paths.get('hisat2'), set())
        tophat_reads = read_sets.get(aligner_paths.get('tophat2'), set())
        
        publication_name = "Unknown"
        for key, value in NAME_MAPPING.items():
            if sample.startswith(key):
                publication_name = value
                break

        total_unique_reads = len(star_reads | hisat_reads | tophat_reads)
        
        all_three_count = len(star_reads & hisat_reads & tophat_reads)
        
        star_hisat_only_count = len((star_reads & hisat_reads) - tophat_reads)
        star_tophat_only_count = len((star_reads & tophat_reads) - hisat_reads)
        hisat_tophat_only_count = len((hisat_reads & tophat_reads) - star_reads)
        
        star_only_count = len(star_reads - (hisat_reads | tophat_reads))
        hisat_only_count = len(hisat_reads - (star_reads | tophat_reads))
        tophat_only_count = len(tophat_reads - (star_reads | hisat_reads))
        
        analysis_results.append({
            'sample': sample,
            'publication': publication_name,
            'STAR_only_count': star_only_count,
            'hisat2_only_count': hisat_only_count,
            'tophat2_only_count': tophat_only_count,
            'STAR_and_hisat2_only_count': star_hisat_only_count,
            'STAR_and_tophat2_only_count': star_tophat_only_count,
            'hisat2_and_tophat2_only_count': hisat_tophat_only_count,
            'all_three_aligners_count': all_three_count,
            'total_unique_reads': total_unique_reads,

            'STAR_only_percent': (star_only_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
            'hisat2_only_percent': (hisat_only_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
            'tophat2_only_percent': (tophat_only_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
            'STAR_and_hisat2_only_percent': (star_hisat_only_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
            'STAR_and_tophat2_only_percent': (star_tophat_only_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
            'hisat2_and_tophat2_only_percent': (hisat_tophat_only_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
            'all_three_aligners_percent': (all_three_count / total_unique_reads * 100) if total_unique_reads > 0 else 0,
        })

    # --- 阶段四: 保存为CSV文件 ---
    if not analysis_results:
        print("[ERROR] No results were generated.")
        return
        
    results_df = pd.DataFrame(analysis_results)
    
    header = [
        'sample', 'publication',
        'STAR_only_count', 'hisat2_only_count', 'tophat2_only_count',
        'STAR_and_hisat2_only_count', 'STAR_and_tophat2_only_count',
        'hisat2_and_tophat2_only_count', 'all_three_aligners_count',
        'total_unique_reads',
        'STAR_only_percent', 'hisat2_only_percent', 'tophat2_only_percent',
        'STAR_and_hisat2_only_percent', 'STAR_and_tophat2_only_percent',
        'hisat2_and_tophat2_only_percent', 'all_three_aligners_percent'
    ]
    results_df = results_df[header]
    
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    results_df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\n--- Analysis complete. Results saved to: {OUTPUT_CSV} ---")

if __name__ == '__main__':
    main()
