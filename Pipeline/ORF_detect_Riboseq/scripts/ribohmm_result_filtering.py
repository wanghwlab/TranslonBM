## #!/usr/bin/env python
import re 
import pandas as pd
def result(gpd_file, prediction_table_file,ribohmm_filtered):
    gpd = pd.read_csv(gpd_file, sep='\t', header=None, names=[
            'transcript_id', 'chrom', 'strand', 'transcript_gstart',
            'transcript_gstop', 'annotated_gstart',
            'annotated_gstop', 'exon_number',
            'exon_gstart_blocks', 'exon_gstop_blocks', 'score',
            'gene_id', 'cdsStartStat', 'cdsEndStat', 'exonFrames'
        ])
    transcript_id_strand = gpd[['transcript_id', 'strand']]
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep=' ')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    ribohmm_result = prediction_table.loc[prediction_table['posterior'] > 8000]
    ribohmm_result.sort_values(by=['posterior'],ascending=False, inplace=True)
    ribohmm_result.rename(columns={'chromosome':'chrom', 'cdstart':'ORF_gstart', 'cdstop':'ORF_gstop'}, inplace = True)
    ribohmm_result = pd.merge(left=ribohmm_result, right=transcript_id_strand, how='inner', on=[
                           'transcript_id'])
    ribohmm_result = ribohmm_result.loc[ribohmm_result.strand_x == ribohmm_result.strand_y,]
    ribohmm_result.rename(columns={'strand_x':'strand'}, inplace = True)
    ribohmm_result.to_csv(ribohmm_filtered, sep='\t', index=False)

result(snakemake.input[1], snakemake.input[0], snakemake.output[0])