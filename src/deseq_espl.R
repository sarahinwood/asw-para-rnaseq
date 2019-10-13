library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")
library("VennDiagram")

dds_abdo <- readRDS("output/asw_timecourse/deseq2/dds_abdo.rds")
##filter results on alpha
dds_abdo_res <- results(dds_abdo, alpha = 0.1)
##order based off padj
ordered_dds_abdo_res <- dds_abdo_res[order(dds_abdo_res$padj),]
df_all_res <- data.table(data.frame(ordered_dds_abdo_res), keep.rownames = TRUE)

##espl res from long tc
trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt")
espl <-  dplyr::filter(trinotate_report, grepl('enhancer of split', eggnog))
espl_longtc_deseq_res <- merge(df_all_res, espl, by.x="rn", by.y="#gene_id", all.x=FALSE, all.y=TRUE)
fwrite(espl_longtc_deseq_res, "output/espl_deseq/long_tc_espl.csv")

##espl res from asw-mh concat long tc
asw_mh_concat_annots <- fread("data/asw_mh_transcriptome/asw_mh_trinotate_annotation_report.txt")
espl_asw_mh_concat <-  dplyr::filter(asw_mh_concat_annots, grepl('enhancer of split', eggnog))
asw_mh_concat_all_res <- fread("output/asw_mh_timecourse/deseq2/asw/timecourse_all_genes.csv")
espl_asw_mh_concat_deseq_res <- merge(asw_mh_concat_all_res, espl_asw_mh_concat, by.x="rn", by.y="#gene_id", all.x=FALSE, all.y=TRUE)
fwrite(espl_asw_mh_concat_deseq_res, "output/espl_deseq/asw_only_analysis_asw_mh_concat_espl.csv")

##espl from short tc
short_tc_all_res <- fread("short_tc_output/asw_timecourse/deseq2/timecourse_all_genes.csv")
espl_short_tc_deseq_res <- merge(short_tc_all_res, espl, by.x="rn", by.y="#gene_id", all.x=FALSE, all.y=TRUE)
fwrite(espl_short_tc_deseq_res, "output/espl_deseq/short_tc_espl.csv")
