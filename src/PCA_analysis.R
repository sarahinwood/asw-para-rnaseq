library("DESeq2")
library("data.table")
library("ggplot2")

##read in dds saved in previous script
dds <- readRDS("output/exposed/deseq2/dds.rds")
dds_group <- readRDS("output/exposed/deseq2/dds_group.rds")

##PCA plot - first must log transform using vst (must set blind = true)
vst <- varianceStabilizingTransformation(dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(vst, intgroup=c("Treatment", "Tissue", "Wasp_Location"))

##pull out results from dds_group
res_group <- results(dds_group, contrast = c("group", "Head_Exposed", "Head_Control"), lfcThreshold = 1, alpha = 0.1)
##only keep genes where padj is not NA
kept_genes <- rownames(subset(res_group, !is.na(padj)))
##create matrix of vst values for only genes where padj didn't = NA
vst_asssay<- assay(vst)[kept_genes,]
##perform PCA on vst data matrix
pc <- prcomp(t(vst_asssay), center = TRUE, scale = TRUE)
##generate data table of results
pc_wide <- data.table(pc$x, keep.rownames = TRUE)
pc_pd <- melt(pc_wide)
fwrite(pc_pd, "output/vst_pca_plot_data.csv")
##plot pcs
ggplot(pc_pd, aes(x=rn, y=value, colour=rn))+
  facet_wrap(~variable)+geom_point()+theme(axis.text.x=element_text(angle = 90))