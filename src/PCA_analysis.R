library("DESeq2")
library("data.table")
library("ggplot2")

sample_data <- fread("data/sample_key.csv")
##read in dds saved in previous script
dds <- readRDS("output/asw_timecourse/deseq2/dds.rds")
dds_abdo <- readRDS("output/asw_timecourse/deseq2/dds_abdo.rds")

##PCA plot - first must log transform using vst (must set blind = true)
vst <- varianceStabilizingTransformation(dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(vst, intgroup=c("Treatment", "Tissue", "Flow_cell"))

##pull out results from dds_group
res_group <- results(dds_abdo, alpha = 0.1)
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

pc_pd_sample_data <- merge(pc_pd, sample_data, by.x="rn", by.y="Sample_name")

##plot pcs
ggplot(pc_pd_sample_data, aes(x=rn, y=value, colour=Treatment, shape=Tissue))+
  facet_wrap(~variable)+geom_point()+theme(axis.text.x=element_text(angle = 90))
