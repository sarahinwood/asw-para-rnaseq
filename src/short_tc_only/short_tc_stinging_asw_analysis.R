library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")
library("VennDiagram")

##ASW gene to trans map (edited to have ID's match concat transcriptome - ASW_TRINITY)
gene2tx <- fread("data/asw_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files - quantified against concat ASW-MH transcriptome
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
##remove long tc samples
colData(dds)
##do 6 times to remove all 6 8h samples at top
dds <- dds[,-1]
##Save dds
saveRDS(dds, file = "output/asw_mh_timecourse/short_tc_deseq2/dds.rds")

##read in dds to save rerunning
dds <- readRDS("output/asw_mh_timecourse/short_tc_deseq2/dds.rds")
counts_matrix <- data.frame(counts(dds))
counts_colSums <- data.frame(colSums(counts_matrix, na.rm=TRUE))
fwrite(counts_colSums, "output/asw_mh_timecourse/short_tc_deseq2/counts_colsums.csv", row.names = TRUE)

##Select only abdomen samples
dds_abdo <- dds[,dds$Tissue == "Abdomen"]
##convert to factors
time_order <- c("Control", "m30", "m120", "m240")
dds_abdo$Treatment <- factor(dds_abdo$Treatment, levels=time_order)
dds_abdo$Wasp_Location <- factor(dds_abdo$Wasp_Location)
##add factors of ineterst to design
design(dds_abdo) <- ~Wasp_Location+Treatment
##run deseq, must specify reduced model for LRT test
dds_abdo <- DESeq(dds_abdo, test = "LRT", reduced = ~Wasp_Location)
##save dds as a file for import in clustering timecourse genes script
saveRDS(dds_abdo, file = "output/asw_mh_timecourse/short_tc_deseq2/dds_abdo.rds")

dds_abdo <- readRDS("output/asw_mh_timecourse/short_tc_deseq2/dds_abdo.rds")
##filter results on alpha
dds_abdo_res <- results(dds_abdo, alpha = 0.1)
##order based off padj
ordered_dds_abdo_res <- dds_abdo_res[order(dds_abdo_res$padj),]
df_all_res <- data.table(data.frame(ordered_dds_abdo_res), keep.rownames = TRUE)

##make list of sig genes
sig_genes <- subset(dds_abdo_res, padj < 0.1)
##make list of sig gene names
sig_gene_names <- row.names(sig_genes)
##save list of sig gene names
fwrite(data.table(sig_gene_names), "output/asw_mh_timecourse/short_tc_deseq2/timecourse_sig_gene_names.csv")

##write list of results for all genes for FGSEA analysis
timecourse_all <- data.table(data.frame(dds_abdo_res), keep.rownames=TRUE)
timecourse_all$fixed_id <- tstrsplit(timecourse_all$rn, "ASW_", keep=c(2))
fwrite(timecourse_all, "output/asw_mh_timecourse/short_tc_deseq2/timecourse_all_genes.csv")

##Order results based of padj
ordered_sig_degs <- sig_genes[order(sig_genes$padj),]
##make datatable and write to output
ordered_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
ordered_degs_table$fixed_ids <- tstrsplit(ordered_degs_table$rn, "ASW_", keep=c(2))
fwrite(ordered_degs_table, "output/asw_mh_timecourse/short_tc_deseq2/timecourse_analysis_sig_degs.csv")

##plot expression pattern for gene
plot_gene <- plotCounts(dds_abdo, "ASW_TRINITY_DN798_c0_g1", 
                        intgroup = c("Treatment"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = Treatment, y = count)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() + xlab("Parasitism Timepoint") + ylab("Normalized Count")
##plot counts for genes of interest, sub in name
plotCounts(dds_abdo, "ASW_TRINITY_DN1343_c0_g1", intgroup = c("Treatment"))

##make table to view counts for each gene
counts_table <- (data.table(counts(dds_abdo), keep.rownames = TRUE))

##read in annotated transcriptome
trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_annotation_report.txt")
setnames(ordered_degs_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
sig_w_annots <- merge(ordered_degs_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/asw_mh_timecourse/short_tc_deseq2/sig_genes_with_annots.csv")

##sort out DEGs that are new and DEGs that have been lost and take a look at annots

##read back in dedeup degs with annots
dedup_sig_w_annots <- fread("output/asw_mh_timecourse/short_tc_deseq2/dedup_sig_genes_with_annots.csv")
##sum of DEGs with no blastX annotation in transcriptome
sum(dedup_sig_w_annots$sprot_Top_BLASTX_hit==".")
##list of DEGs with no blastX annotation
no_blastx_annot_degs <- dedup_sig_w_annots[dedup_sig_w_annots$sprot_Top_BLASTX_hit == ".",]
##list of DEGs with no blastX OR blastP
no_blast_annot_degs <- no_blastx_annot_degs[no_blastx_annot_degs$sprot_Top_BLASTP_hit == ".",]
##make list of degs with no blast annot.
list_degs_no_annot <- data.table(no_blast_annot_degs$`#gene_id`)
##write list of degs with no annot.
fwrite(list_degs_no_annot, "output/asw_mh_timecourse/short_tc_deseq2/no_annot/degs_with_no_annot.txt")

##look at cellular immune response genes
melanization_annots <- dplyr::filter(trinotate_report, grepl('melanization', gene_ontology_blast))
encapsulation_annots <- dplyr::filter(trinotate_report, grepl('encapsulation of foreign target', gene_ontology_blast))
defense_annots <- dplyr::filter(trinotate_report, grepl('defense response to insect', gene_ontology_blast))
##filter out deseq results for cellular immunity genes
cellular_immunity_annots <- dplyr::filter(trinotate_report,grepl('melanization|encapsulation of foreign target|defense response to insect|melanin', gene_ontology_blast))
cellular_immunity_genes <- data.table(unique(cellular_immunity_annots$`#gene_id`))
cell_imm_rnaseq_res <- merge(cellular_immunity_genes, timecourse_all, by.x="V1", by.y="rn")
fwrite(cell_imm_rnaseq_res, "output/asw_mh_timecourse/short_tc_deseq2/cellular_immune_deseq_res.csv")





