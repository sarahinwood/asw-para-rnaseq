library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")

gene2tx <- fread("data/asw_Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

  ##Find all salmon quant files
quant_files <- list.files(path="output/salmon/asw", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
  ##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
  ##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
  ##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
  ##Import table describing samples
sample_data <- fread("data/full_sample_key.csv")
setkey(sample_data, Sample_name)

  ##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
  ##create dds object for group analysis
dds_group <- copy(dds)
  ##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
  ##add group to design
design(dds_group) <- ~group
  ##run deseq2 and generate results
dds_group <- DESeq(dds_group)
resultsNames(dds_group)
  ##Make table of results for exposed vs control heads
res_group <- results(dds_group, contrast = c("group", "Head_Exposed", "Head_Control"), lfcThreshold = 1, alpha = 0.1)
  ##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
  ##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/exposed/deseq2/exposed_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

##Sub in any gene of interest to plot counts  
plotCounts(dds_group, "TRINITY_DN2391_c0_g3", intgroup = c("group"))

##PCA plot - first must log transform
##variance stabilising transformations method, impt to set blind=true
vst <- vst(dds, blind=TRUE)
##plot PCA to investigate sample clustering
plotPCA(vst, intgroup=c("Treatment", "Tissue", "Wasp_Location"))

  ##read in annotated transcriptome
trinotate_report <- fread("data/trinotate_annotation_report.txt")
setnames(ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
  ##merge list of sig genes with annotations
sig_w_annots <- merge (ordered_sig_res_group_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
  ##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/exposed/deseq2/sig_genes_with_annots.csv")

  ##read back in dedeup sig genes w/annots
dedup_sig_w_annots <- fread("output/exposed/deseq2/dedup_sig_genes_with_annots.csv")
  ##sum of DEGs with no blastX annotation in transcriptome
sum(dedup_sig_w_annots$sprot_Top_BLASTX_hit==".")
  ##list of DEGs with no blastX annotation
no_blastx_annot_degs <- dedup_sig_w_annots[dedup_sig_w_annots$sprot_Top_BLASTX_hit == ".",]
  ##list of DEGs with no blastX OR blastP
no_blast_annot_degs <- no_blastx_annot_degs[no_blastx_annot_degs$sprot_Top_BLASTP_hit == ".",]
  ##make list of degs with no blast annot.
list_degs_no_annot <- data.table(no_blast_annot_degs$transcript_id)
  ##write list of degs with no annot.
fwrite(list_degs_no_annot, "output/exposed/deseq2/degs_with_no_annot.txt")

  ##dedup annotations in excel and then read in and format blast results for unann degs
blastx_unann_degs <- fread("output/exposed/no_annot/blastx_titles.outfmt6")
setnames(blastx_unann_degs, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("#gene_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
fwrite(blastx_unann_degs, "output/exposed/no_annot/blastx_exposed_results.csv")
##unique gene ids = 121
blastx_unann_degs[,unique(`#gene_id`)]
#read in dedup annotations and check unique id = 121
dedup_blast_exposed <- fread("output/exposed/no_annot/dedup_blastx_exposed_results.csv")
dedup_blast_exposed[,unique(`#gene_id`)]
sig_blastx_trinotate_annots <- merge(dedup_sig_w_annots, dedup_blast_exposed, by.x="annotation_transcript_id", by.y="#gene_id", all = TRUE)
fwrite(sig_blastx_trinotate_annots, "output/exposed/deseq2/sig_genes_trinotate_blastx_annots.csv")
