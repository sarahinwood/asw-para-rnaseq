library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")
library("VennDiagram")


gene2tx <- fread("data/asw_transcriptome/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files from salmon filtering because res. from filtering with STAR and lost bro that way
quant_files <- list.files(path="output/asw_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv")
setkey(sample_data, Sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds, file = "output/asw_timecourse_pairwise/deseq2/dds.rds")

##create dds object for group analysis
dds_group <- copy(dds)
##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
saveRDS(dds_group, file = "output/asw_timecourse_pairwise/deseq2/dds_group.rds")

resultsNames(dds_group)

##Make table of results for control vs timepoint
res_group <- results(dds_group, contrast = c("group", "Abdomen_m240", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/asw_timecourse_pairwise/deseq2/240m_res_group.csv")
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/asw_timecourse_pairwise/deseq2/240m_timecourse_pairwise_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)
##Sub in any gene of interest to plot counts  
plotCounts(dds_group, "TRINITY_DN35519_c0_g1", intgroup = c("group"), main="")


m30_deg_list <- fread("output/asw_timecourse_pairwise/deseq2/30m_timecourse_pairwise_analysis_sig_degs.csv")
m30_degs <- m30_deg_list$rn

m120_deg_list <- fread("output/asw_timecourse_pairwise/deseq2/120m_timecourse_pairwise_analysis_sig_degs.csv")
m120_degs <- m120_deg_list$rn

m240_deg_list <- fread("output/asw_timecourse_pairwise/deseq2/240m_timecourse_pairwise_analysis_sig_degs.csv")
m240_degs <- m240_deg_list$rn

lrt_degs <- fread("output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")
lrt_ids <- lrt_degs$rn

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("m30"=m30_degs, "m120"=m120_degs, "m240"=m240_degs, "lrt"=lrt_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

annots <- fread("data/asw_transcriptome/trinotate_annotation_report.txt")

m30_deg_annots <- merge(m30_deg_list, annots, by.x="rn", by.y="#gene_id", all.x=TRUE)
m120_deg_annots <- merge(m120_deg_list, annots, by.x="rn", by.y="#gene_id", all.x=TRUE)
m240_deg_annots <- merge(m240_deg_list, annots, by.x="rn", by.y="#gene_id", all.x=TRUE)

