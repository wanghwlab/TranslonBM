.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
#! /usr/bin/env R
library(ORFquant)

## curate ORFquant_result

ORFquant_result <- load(snakemake@input[[1]])
gr <- ORFquant_results$ORFs_tx

df <- data.frame(
    seqnames = seqnames(gr),
    ranges = ranges(gr),
    genomic_region = elementMetadata(gr)[, c("region")],
    pvalue = elementMetadata(gr)[, c("pval")],
    ORF_category_Txe = elementMetadata(gr)[, c("ORF_category_Tx")],
    ORF_category_Tx_compatible = elementMetadata(gr)[, c("ORF_category_Tx_compatible")],
    ORF_category_Gene = elementMetadata(gr)[, c("ORF_category_Gen")]
)

df <- df[, c("seqnames", "ranges.start", "ranges.end", "genomic_region.seqnames", "genomic_region.strand", "pvalue", "genomic_region.start", "genomic_region.end")]
colnames(df) <- c("transcript_id", "ORF_tstart", "ORF_tstop", "chrom", "strand", "pvalue", "ORF_gstart", "ORF_gstop")

write.table(df, snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)