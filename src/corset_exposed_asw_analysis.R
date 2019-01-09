library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(plyr)
library(VennDiagram)

##read in corset output
count_data<-read.delim("output/corset/counts.txt",row.names=1, check.names = FALSE)

##Import table describing samples
sample_data <- fread("data/full_sample_key.csv")
setkey(sample_data, Sample_name)

##make dds object
dds <- DESeqDataSetFromMatrix(count_data, colData = sample_data[colnames(count_data)], design = ~1)
##save dds object
saveRDS(dds, file = "output/exposed_corset/deseq2/dds.rds")
##create dds object for group analysis
dds_group <- copy(dds)
##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
saveRDS(dds_group, file = "output/exposed_corset/deseq2/dds_group.rds")

resultsNames(dds_group)
##Make table of results for exposed vs control heads
res_group <- results(dds_group, contrast = c("group", "Head_Exposed", "Head_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/exposed_corset/deseq2/res_group.csv")
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/exposed_corset/deseq2/exposed_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

##Sub in any gene of interest to plot counts  
plotCounts(dds_group, "Cluster-2682.0", intgroup = c("group"), main="...")

##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

##read in annotations
trinotate_report <- fread("data/trinotate_annotation_report.txt")
##read in corset clusters
cluster_data<-read.delim("output/corset/clusters.txt", header = FALSE)
##Generate counts of transcripts in each cluster
cluster_counts <- count(cluster_data, vars="V2")

##generate table of transcript, annot + cluster allocation
trinotate_clusters<-merge(cluster_data, trinotate_report, by.x="V1", by.y="transcript_id", all.x=TRUE)

##merge annot+clusters with list of DEGs
degs_annots <- merge(trinotate_clusters, ordered_sig_res_group_table, by.x="V2", by.y="rn")
fwrite(degs_annots, "output/exposed_corset/deseq2/sig_clusters_annots.csv")

##look at overlap with original analysis w/out clustering
deg_ids <- data.frame(tstrsplit(degs_annots$V1, "_i", keep=1))
setnames(deg_ids, old=c("c..TRINITY_DN25575_c0_g1....TRINITY_DN39667_c0_g1....TRINITY_DN1916_c0_g1..."), new=c("gene_id"))
deg_iso_ids <- deg_ids$gene_id

##read in longest iso/gene results
no_corset_anal <- fread("output/exposed/deseq2/exposed_analysis_sig_degs.csv")
old_deg_ids <- no_corset_anal$rn
##Venn diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Corset DEGs"=deg_iso_ids, "Longest Isoform/Gene DEGs"=old_deg_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)