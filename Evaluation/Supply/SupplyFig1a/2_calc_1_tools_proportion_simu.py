import os
import pandas as pd

# 1. Define datasets, aligners, and ORF prediction tools
datasets = [
    #"SRX876063_SRX876069", "SRX740748", "SRX1254413", 
    #"SRX5256543_SRX5256555", "SRX5887328_SRX5887329_SRX5887330", 
    #"SRX11812007_SRX11812008_SRX11812009"
    "simulation_6M_T1","simulation_6M_T3","simulation_60M_T1","simulation_60M_T3"
]
aligners = ["tophat2", "hisat2", "STAR"]
tools = [
    "riborf", "ribocode", "ribotricer", "ribohmm", "ribotish", 
    "ribowave", "orfrater", "orfquant", "gedi", "ribotaper", "rpbp"
]

# 2. Define file paths
overlap_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/tools_overlap/merged_ATG/orf_pred_default_untrim/"
pred_dir = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/orf_pred_default_untrim/"
output_csv = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/1tools_percent/1_tools_proportion_summary_simu_untrim.csv"

results_list = []

print("Calculating the software-source proportions for 1_tools...")

# 3. Iterate over all dataset and aligner combinations
for ds in datasets:
    for al in aligners:
        overlap_file = os.path.join(overlap_dir, f"{ds}_{al}_overlap_count.txt")
        
        if not os.path.exists(overlap_file):
            print(f"[Skipped] File not found: {overlap_file}")
            continue
            
        print(f"Processing: {ds} - {al}")
        
        df_overlap = pd.read_csv(overlap_file, sep=r'\s+')
        
        if 'tools' not in df_overlap.columns or 'coordinate_id' not in df_overlap.columns:
            print(f"[Warning] File {overlap_file} is missing required columns; check the format.")
            continue
            
        set_1_tools = set(df_overlap[df_overlap['tools'] == '1_tools']['coordinate_id'])
        total_1_tools = len(set_1_tools)
        
        row_data = {
            'Dataset': ds,
            'Aligner': al,
            'Total_1_tools_ORFs': total_1_tools
        }
        
        if total_1_tools == 0:
            print(f"  -> {ds} {al} contains no ORFs assigned to 1_tools.")
            for tool in tools:
                row_data[f'{tool}_count'] = 0
                row_data[f'{tool}_percent(%)'] = 0.0
            results_list.append(row_data)
            continue
            
        # 4. Check the 11 ORF prediction tools for overlaps
        for tool in tools:
            pred_file = os.path.join(pred_dir, f"{ds}_{al}_{tool}_gcoor.tsv.gz")
            
            if os.path.exists(pred_file):
                try:
                    df_pred = pd.read_csv(pred_file, sep='\t', usecols=['coordinate_id'], compression='gzip')
                    tool_coords = set(df_pred['coordinate_id'])
                    overlap_count = len(tool_coords.intersection(set_1_tools))
                    
                except Exception as e:
                    print(f"  [Error] Reading file {pred_file} failed: {e}")
                    overlap_count = 0
            else:
                overlap_count = 0
                
            percentage = (overlap_count / total_1_tools) * 100
            
            row_data[f'{tool}_count'] = overlap_count
            row_data[f'{tool}_percent(%)'] = round(percentage, 2)
            
        results_list.append(row_data)

# 5. Save results as CSV
df_results = pd.DataFrame(results_list)
df_results.to_csv(output_csv, index=False)

print(f"\nStatistics completed. Results saved in the current directory as: {output_csv}")
