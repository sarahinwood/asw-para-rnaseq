library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")
library("VennDiagram")

##ASW gene to trans map (edited to have ID's match concat transcriptome - ASW_TRINITY)
gene2tx <- fread("data/asw_transcriptome/Trinity.fasta.gene_trans_map", header = FALSE)
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
saveRDS(dds, file = "output/asw_timecourse/deseq2/dds.rds")

##read in dds to save rerunning
dds <- readRDS("output/asw_timecourse/deseq2/dds.rds")

##remove 8h NC2&3 - find sample number in table
colData(dds)
##remove 8hNC2 - 4th in table
dds <- dds[,-4]
##remove 8hNC3 - sample 5 once 8hNC2 removed
dds <- dds[,-5]
##check correct samples were removed
colData(dds)

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
##save dds as a file for import in clustering timecourse genes script
saveRDS(dds_abdo, file = "output/asw_timecourse/deseq2/dds_abdo.rds")

dds_abdo <- readRDS("output/asw_timecourse/deseq2/dds_abdo.rds")
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
fwrite(data.table(sig_gene_names), "output/asw_timecourse/deseq2/timecourse_sig_gene_names.csv")

##write list of results for all genes for FGSEA analysis
timecourse_all <- data.table(data.frame(dds_abdo_res), keep.rownames=TRUE)
timecourse_all$fixed_id <- tstrsplit(timecourse_all$rn, "ASW_", keep=c(2))
fwrite(timecourse_all, "output/asw_timecourse/deseq2/timecourse_all_genes.csv")

##Order results based of padj
ordered_sig_degs <- sig_genes[order(sig_genes$padj),]
##make datatable and write to output
ordered_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
ordered_degs_table$fixed_ids <- tstrsplit(ordered_degs_table$rn, "ASW_", keep=c(2))
fwrite(ordered_degs_table, "output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")

##plot expression pattern for gene
plot_gene <- plotCounts(dds_abdo, "ASW_TRINITY_DN920_c0_g1", 
                        intgroup = c("Treatment"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = Treatment, y = count)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() + xlab("Parasitism Timepoint") + ylab("Normalized Count")
##plot counts for genes of interest, sub in name
plotCounts(dds_abdo, "ASW_TRINITY_DN9264_c0_g1", intgroup = c("Treatment", "Flow_cell"))

##make table to view counts for each gene
counts_table <- (data.table(counts(dds_abdo), keep.rownames = TRUE))

  ##read in annotated transcriptome
trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt")
setnames(ordered_degs_table, old=c("rn"), new=c("#gene_id"))
  ##merge list of sig genes with annotations
sig_w_annots <- merge(ordered_degs_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
  ##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/asw_timecourse/deseq2/sig_genes_with_annots.csv")

##sort out DEGs that are new and DEGs that have been lost and take a look at annots

  ##read back in dedeup degs with annots
dedup_sig_w_annots <- fread("output/asw_timecourse/deseq2/dedup_sig_genes_with_annots.csv")
  ##sum of DEGs with no blastX annotation in transcriptome
sum(dedup_sig_w_annots$sprot_Top_BLASTX_hit==".")
  ##list of DEGs with no blastX annotation
no_blastx_annot_degs <- dedup_sig_w_annots[dedup_sig_w_annots$sprot_Top_BLASTX_hit == ".",]
  ##list of DEGs with no blastX OR blastP
no_blast_annot_degs <- no_blastx_annot_degs[no_blastx_annot_degs$sprot_Top_BLASTP_hit == ".",]
  ##make list of degs with no blast annot.
list_degs_no_annot <- data.table(no_blast_annot_degs$`#gene_id`)
  ##write list of degs with no annot.
fwrite(list_degs_no_annot, "output/asw_timecourse/deseq2/no_annot/degs_with_no_annot.txt")

  ##read in blast results?
timecourse_no_annot_blast <- fread("output/asw_timecourse/no_annot/blastx_titles.outfmt6")
  ##rename columns
setnames(timecourse_no_annot_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("#gene_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
  ##save formatted blast results (then dedup annotations in excel - keep annot with highest e-value)
fwrite(timecourse_no_annot_blast, "output/asw_timecourse/no_annot/blastx_timecourse_annotations.csv")

  ##read in dedup blast results
dedup_timecourse_blast <- fread("output/asw_timecourse/no_annot/dedup_blastx_timecourse_annotations.csv")
dedup_timecourse_blast[,unique(`#gene_id`)]

  ##merge blastx annotations for unann transcripts with transcriptome annotations
all_annots_degs <- merge(dedup_sig_w_annots, dedup_timecourse_blast, by.x="transcript_id", by.y="#gene_id", all = TRUE)
fwrite(all_annots_degs, "output/asw_timecourse/deseq2/degs_trinotate_blastx_annots.csv")

##filter out genes in blastx annotation column that contain "uncharacterized" or "hypothetical"
unchar_or_hypo_annots <- dplyr::filter(all_annots_degs, grepl('uncharacterized|hypothetical', annotation))
##filter out genes with no manual annotation OR trinotate blastx annotation
no_manual_annot <- all_annots_degs %>% filter(is.na(annotation))
no_annot <- no_manual_annot[no_manual_annot$sprot_Top_BLASTX_hit == ".",]
##merge list of genes with no annot OR hypothetical/uncharacterised and save for interproscan
unchar_hypo_ids <- data.table(unchar_or_hypo_annots$transcript_id)
noannot_ids <- data.table(no_annot$transcript_id)
ids_for_interproscan <- merge(unchar_hypo_ids, noannot_ids, all = TRUE)
fwrite(ids_for_interproscan, "output/asw_timecourse/interproscan/interproscan_ids.txt")

##compare to shorter-timecourse (no 8h)
short_tc_degs <- fread("output/asw_timecourse_all_NC/deseq2/timecourse_analysis_sig_degs.csv")
short_tc_degs$fixed_ids <- tstrsplit(short_tc_degs$rn, "ASW_", keep=c(2))
short_tc_ids <- short_tc_degs$fixed_ids
long_sig_id_list <- ordered_degs_table$fixed_ids
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Long TC, all NCs"=short_tc_ids, "Long TC -8hNC2&3"=long_sig_id_list), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

short_tc_only <- data.table(setdiff(short_tc_degs$rn, ordered_degs_table$rn))
short_tc_only_degs <- subset(short_tc_degs, (rn %in% short_tc_only$V1))
annot_short_tc_only <- merge(short_tc_only_degs, dedup_sig_w_annots, by.x ="rn", by.y="#gene_id")
fwrite(annot_short_tc_only, "output/asw_timecourse/deseq2/short_tc_only_degs.csv")

long_tc_only <- data.table(setdiff(ordered_degs_table$rn, short_tc_degs$rn))
long_tc_only_degs <- subset(ordered_degs_table, (rn %in% long_tc_only$V1))
annot_long_tc_only <- merge(long_tc_only_degs, dedup_sig_w_annots, by.x ="rn", by.y="#gene_id")
fwrite(annot_long_tc_only, "output/asw_timecourse/deseq2/long_tc_only_degs.csv")
