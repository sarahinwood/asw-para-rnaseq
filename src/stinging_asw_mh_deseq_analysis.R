library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")
library("VennDiagram")

gene2tx <- fread("data/asw_mh_transcriptome/asw_mh_Trinity.fasta.gene_trans_map", header = FALSE)
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
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)

##Select only abdomen samples
dds_abdo <- dds[,dds$Tissue == "Abdomen"]
##convert to factors
time_order <- c("Control", "m30", "m120", "m240", 'm480')
dds_abdo$Treatment <- factor(dds_abdo$Treatment, levels=time_order)
dds_abdo$Wasp_Location <- factor(dds_abdo$Wasp_Location)
dds_abdo$Flow_cell <- factor(dds_abdo$Flow_cell)
##add factors of ineterst to design
design(dds_abdo) <- ~Flow_cell+Wasp_Location+Treatment
##run deseq, must specify reduced model for LRT test
dds_abdo <- DESeq(dds_abdo, test = "LRT", reduced = ~Flow_cell+Wasp_Location)
##filter results on alpha
dds_abdo_res <- results(dds_abdo, alpha = 0.1)
##order based off padj
ordered_dds_abdo_res <- dds_abdo_res[order(dds_abdo_res$padj),]
##save dds as a file for import in clustering timecourse genes script
saveRDS(dds, file = "output/asw_timecourse/deseq2/dds.rds")
##make list of sig genes
sig_genes <- subset(dds_abdo_res, padj < 0.1)
##make list of sig gene names
sig_gene_names <- row.names(sig_genes)
##save list of sig gene names
fwrite(data.table(sig_gene_names), "output/asw_mh_timecourse/deseq2/timecourse_sig_gene_names.csv")

##write list of results for all genes for FGSEA analysis
timecourse_all <- data.table(data.frame(dds_abdo_res), keep.rownames=TRUE)
fwrite(timecourse_all, "output/asw_mh_timecourse/deseq2/timecourse_all_genes.csv")

##Order results based of padj
ordered_sig_degs <- sig_genes[order(sig_genes$padj),]
##make datatable and write to output
ordered_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
fwrite(ordered_degs_table, "output/asw_mh_timecourse/deseq2/timecourse_analysis_sig_degs.csv")

##read in annotated transcriptome
trinotate_report <- fread("data/asw_mh_transcriptome/asw_mh_trinotate_annotation_report.txt")
setnames(ordered_degs_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
sig_w_annots <- merge(ordered_degs_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/asw_mh_timecourse/deseq2/sig_genes_with_annots.csv")