library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")

gene2tx <- fread("data/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

  ##Find all salmon quant files
quant_files <- list.files(path="output/salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
  ##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
  ##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
  ##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
  ##Import table describing samples
sample_data <- fread("data/full_sample_key.csv")
setkey(sample_data, Sample_name)

  ##Create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)

  ##Select only abdomen samples
dds_abdo <- dds[,dds$Tissue == "Abdomen"]
  ##convert to factors
dds_abdo$Treatment <- factor(dds_abdo$Treatment)
dds_abdo$Wasp_Location <- factor(dds_abdo$Wasp_Location)
  ##add factors of ineterst to design
design(dds_abdo) <- ~Wasp_Location+Treatment
  ##run deseq, must specify reduced model for LRT test
dds_abdo <- DESeq(dds_abdo, test = "LRT", reduced = ~Wasp_Location)
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
fwrite(data.table(sig_gene_names), "output/asw_timecourse/deseq2/timecourse_sig_gene_names.csv")

##write list of results for all genes for FGSEA analysis
timecourse_all <- data.table(data.frame(dds_abdo_res), keep.rownames=TRUE)
fwrite(timecourse_all, "output/asw_timecourse/deseq2/timecourse_all_genes.csv")

  ##Order results based of padj
ordered_sig_degs <- sig_genes[order(sig_genes$padj),]
  ##make datatable and write to output
ordered_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
fwrite(ordered_degs_table, "output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")

  ##plot counts for genes of interest, sub in name
plotCounts(dds_abdo, "TRINITY_DN13642_c0_g1", intgroup = c("Treatment", "Wasp_Location"))

  ##read in annotated transcriptome
trinotate_report <- fread("data/trinotate_annotation_report.txt")
setnames(ordered_degs_table, old=c("rn"), new=c("#gene_id"))
  ##merge list of sig genes with annotations
sig_w_annots <- merge(ordered_degs_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
  ##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/asw_timecourse/deseq2/sig_genes_with_annots.csv")

  ##read back in dedeup degs with annots
dedup_sig_w_annots <- fread("output/asw_timecourse/deseq2/dedup_sig_genes_with_annots.csv")
  ##sum of DEGs with no blastX annotation in transcriptome
sum(dedup_sig_w_annots$sprot_Top_BLASTX_hit==".")
  ##list of DEGs with no blastX annotation
no_blastx_annot_degs <- dedup_sig_w_annots[dedup_sig_w_annots$sprot_Top_BLASTX_hit == ".",]
  ##list of DEGs with no blastX OR blastP
no_blast_annot_degs <- no_blastx_annot_degs[no_blastx_annot_degs$sprot_Top_BLASTP_hit == ".",]
  ##make list of degs with no blast annot.
list_degs_no_annot <- data.table(no_blast_annot_degs$transcript_id)
  ##write list of degs with no annot.
fwrite(list_degs_no_annot, "output/asw_timecourse/deseq2/degs_with_no_annot.txt")

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

