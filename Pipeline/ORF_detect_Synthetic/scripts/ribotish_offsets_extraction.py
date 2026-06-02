import pandas as pd

def offsets_extraction(ribotish_offsets,standard_offsets):
    with open(ribotish_offsets, 'r') as f:
        offsets = f.read()
        loc = locals()        
        exec(offsets)
        offdict = loc['offdict']
        del offdict['m0']
        offsetsdf = pd.DataFrame(list(offdict.items()))
        offsetsdf.to_csv(standard_offsets, index = None, header = None, sep = '\t' )

offsets_extraction(snakemake.input[0],snakemake.output[0])