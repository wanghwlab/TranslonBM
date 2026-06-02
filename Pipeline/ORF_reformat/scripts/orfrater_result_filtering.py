import pandas as pd
import shutil

def filter_result(prediction_table_file,orfrater_08):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep=',')
    #biotype = pd.read_csv(biotype, sep='\t', header=0)
    prediction_table.rename(columns={'tid':'transcript_id', 'tcoord':'ORF_tstart', 'tstop':'ORF_tstop'}, inplace = True)
    #prediction_table = pd.merge(left=prediction_table, right=biotype, on='transcript_id')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    orfrater_result = prediction_table.loc[prediction_table['orfrating'] > 0.8]
    orfrater_result.to_csv(orfrater_08, sep='\t', index=False)

filter_result(snakemake.input[0], snakemake.output[0])