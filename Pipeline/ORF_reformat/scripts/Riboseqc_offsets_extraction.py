import pandas as pd

def extract_rl(riboseqc_rl,standard_offsets):
    df = pd.read_csv(riboseqc_rl, sep="\t")
    df = df[df['max_coverage']]
    df = df[['read_length', 'cutoff']]
    df = df.sort_values(by=['read_length'])
    df.to_csv(standard_offsets, sep="\t", index=False, header=False)

extract_rl(snakemake.input[1],snakemake.output[0])