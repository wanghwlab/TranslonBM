import pandas as pd

def offsets_extraction(ribotish_offsets,standard_offsets):
    ribotish_offsets = pd.read_csv(ribotish_offsets,sep = '\t',header = None)
    ribotish_offsets[1] = ribotish_offsets[1] - 1
    ribotish_offsets.to_csv(standard_offsets,sep='\t',index=False,header = 0)

offsets_extraction(snakemake.input[0],snakemake.output[0])