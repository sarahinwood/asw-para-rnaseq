##lincoln-mh-exposed asw vs all controls

library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")
library("VennDiagram")

gene2tx <- fread("data/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/salmon/asw/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
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
dds_linc_group <- copy(dds)
##create groupings of tissue+treatment
dds_linc_group$group <- factor(paste(dds$Tissue,dds$Treatment,dds$Wasp_Location,sep="_"))
##add group to design
design(dds_linc_group) <- ~group
##run deseq2 and generate results
dds_linc_group <- DESeq(dds_linc_group)

resultsNames(dds_linc_group)
##Make table of results for exposed vs control heads
res_linc_group <- results(dds_linc_group, contrast = c("group", "Head_Exposed_Lincoln", "Head_Control_Lincoln"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_linc_group <- res_linc_group[order(res_linc_group$padj),]
##Make data table and write to output
ordered_res_linc_group_table <- data.table(data.frame(ordered_res_linc_group), keep.rownames = TRUE)
fwrite(ordered_res_linc_group_table, "output/linc_exposed/deseq2/all_degs_linc_exposed_vs_LINC_controls.csv")
ordered_sig_res_linc_group_table <- subset(ordered_res_linc_group_table, padj < 0.05)
fwrite(ordered_sig_res_linc_group_table, "output/linc_exposed/deseq2/linc_exposed_vs_LINC_controls_sig_degs.csv", col.names = TRUE, row.names = FALSE)

plotCounts(dds_linc_group, "TRINITY_DN38122_c0_g1", intgroup = c("group"))

##volcano plot
EnhancedVolcano(ordered_res_linc_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

##read in annotated transcriptome
trinotate_report <- fread("data/trinotate_annotation_report.txt")
setnames(ordered_sig_res_linc_group_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
sig_w_annots <- merge (ordered_sig_res_linc_group_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(sig_w_annots, "output/linc_exposed/deseq2/linc_ex_vs_LINC_control_sig_genes_with_annots.csv")

linc_v_all_controls <- fread("output/linc_exposed/deseq2/linc_exposed_vs_ALL_controls_sig_degs.csv")
linc_v_linc_controls <- fread("output/linc_exposed/deseq2/linc_exposed_vs_LINC_controls_sig_degs.csv")
all_v_all <- fread("output/exposed/deseq2/exposed_analysis_sig_degs.csv")

linc_v_all_names <- linc_v_all_controls$rn
linc_v_linc_names <- linc_v_linc_controls$rn
all_v_all_names <- all_v_all$rn

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Linc-MH-Ex ASW v Linc Controls"=linc_v_linc_names, "Linc-MH-Ex ASW v All Controls"=linc_v_all_names, "Old Analysis With Ruakura-MH-Ex"=all_v_all_names), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)



