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
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)

##remove samples that aren't what they should be
coldata <- data.frame(colData(dds))
dds <- dds[,-c(1,4,6,7,8,9,11,18)]
coldata_samples_removed <- data.frame(colData(dds))

##read in dds to save rerunning
counts_matrix <- data.frame(counts(dds))
counts_colSums <- data.frame(colSums(counts_matrix, na.rm=TRUE))
fwrite(counts_colSums, "output/asw_timecourse/deseq2/counts_colsums.csv", row.names = TRUE)

##start analysing
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
saveRDS(dds, file = "output/asw_mh_timecourse/deseq2/asw_dds.rds")
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

long_tc_sig <- fread("output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")
long_tc_ids <- long_tc_sig$rn
##compare asw+mh at once to long tc
asw_sig_degs <- dplyr::filter(ordered_degs_table, grepl('ASW_TRINITY', `#gene_id`))
asw_sig_ids <- data.table(asw_sig_degs$`#gene_id`)
asw_sig_ids$V1 <- tstrsplit(asw_sig_ids$V1, "ASW_", keep=c(2))
plot_asw_sig_ids <- asw_sig_ids$V1

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Long Time-course DEGs"=long_tc_ids, "ASW DEGs (ASW-MH concat)"=plot_asw_sig_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

##compare asw-only from asw_mh concat to long tc
asw_sig_degs <- data.table(ordered_degs_table$rn)
asw_sig_degs$V1 <- tstrsplit(asw_sig_degs$V1, "ASW_", keep=c(2))
asw_sig_ids <- asw_sig_degs$V1

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Long Time-course DEGs"=long_tc_ids, "ASW DEGs (ASW-MH concat,
                            ASW-only deseq)"=asw_sig_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

##read in annotated transcriptome
trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_annotation_report.txt")
setnames(ordered_degs_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
sig_w_annots <- merge(ordered_degs_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/asw_mh_timecourse/deseq2/sig_genes_with_annots.csv")


