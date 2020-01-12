library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")
library("VennDiagram")

gene2tx <- fread("data/asw_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv")
setkey(sample_data, Sample_name)

##Create dds object and link to sample data
dds_para_status <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##only abdo samples
dds_para_status_group <- copy(dds_para_status)
dds_para_status_group$group <- factor(paste(dds_para_status$Tissue,dds_para_status$Parasitism_status,sep="_"))
##add group to design
design(dds_para_status_group) <- ~group
##run deseq2 and generate results
dds_para_status_group <- DESeq(dds_para_status_group)
saveRDS(dds_para_status_group, file = "output/parasitism_status_pairwise/deseq2/para_status_dds.rds")

resultsNames(dds_para_status_group)

##para vs unpara abdomens
res_group <- results(dds_para_status_group, contrast = c("group", "Abdomen_Para", "Abdomen_Unpara"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
ordered_sig_res_group_table$rn <- tstrsplit(ordered_sig_res_group_table$rn, "ASW_", keep=c(2))
fwrite(ordered_res_group_table, "output/parasitism_status_pairwise/deseq2/res_group.csv")
fwrite(ordered_sig_res_group_table, "output/parasitism_status_pairwise/deseq2/exposed_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

trinotate_report <- fread("data/asw_most_sig_transcript_blastx_hit_for_each_gene.csv")
sig_degs_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_degs_annots, "output/parasitism_status_pairwise/deseq2/sig_degs_annots.csv")

plotCounts(dds_para_status_group, "ASW_TRINITY_DN38122_c0_g1", intgroup = c("group"), main="Bro Frag - TRINITY_DN38122_c0_g1")

counts_matrix <- data.table(data.frame(counts(dds_para_status_group)), keep.rownames = TRUE)
kila_counts <- melt(data.table(counts_matrix[counts_matrix$rn == "ASW_TRINITY_DN1684_c0_g1",]))
kila_counts$h8_fixed <- tstrsplit(kila_counts$variable, "X", keep=c(2))
kila_counts_sample_data <- merge(sample_data, kila_counts, by.x="Sample_name", by.y="variable")
kila_counts_sample_data_h8 <- merge(sample_data, kila_counts, by.x="Sample_name", by.y="h8_fixed")
full_sample_data_kila_counts <-full_join(kila_counts_sample_data, kila_counts_sample_data_h8)

##highest kil in parasitised - 	L2_30m_A
##highest ex - L1_Ex_H


