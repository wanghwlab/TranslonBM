import pandas as pd

def extract_rl(riboseqc_rl,rl_for_ORFquant):
    df = pd.read_csv(riboseqc_rl, sep="\t") 
    df = df[["read_length", "cutoff", "comp"]]
    df = df.rename(columns={"comp": "compartment"})
    df.to_csv(rl_for_ORFquant, sep="\t", index=False)

extract_rl(snakemake.input[0],snakemake.output[0])
