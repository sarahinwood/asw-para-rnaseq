library("DESeq2")
library("data.table")

  ##function to calculate geometric mean
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

  ##read in dds saved in previous script
dds <- readRDS("output/asw_timecourse/deseq2/dds.rds")
  ##filter for only abdo samples
dds_abdo <- dds[,dds$Tissue == "Abdomen"]
  ##read in list of sig gene names
sig_gene_names <- fread("output/asw_timecourse/deseq2/timecourse_sig_gene_names.csv")

  ##vst log transform data - absolute counts now changed so cannot compare between genes (only between samples for 1 gene)
vst <- varianceStabilizingTransformation(dds_abdo, blind = FALSE)
  ##make matrix of transformed data
vst_matrix <- data.table(as.matrix(assay(vst)), keep.rownames = TRUE)
  ##melt to make long table rather than wide
long_vst_data <- melt(vst_matrix, id.vars = "rn", variable.name = "Sample_name", value.name = "vst")
  ##make table of colData from dds object
long_coldata <- data.table(as.data.frame(colData(dds_abdo)))

  ##merge sample data with vst data
merged_exp_coldata <- merge(long_vst_data, long_coldata, all.x = TRUE, all.y = FALSE)
  ##generate table of mean vst values
mean_vst <- merged_exp_coldata[,.(vst_mean = gm_mean(vst)),by=.(rn, Treatment)]
  ##make long table a wide table instead
mean_vst_wide <- dcast(mean_vst, rn~Treatment)
  ##make matrix with gene names as row names
expression_matrix <- as.matrix(data.frame(mean_vst_wide, row.names = "rn"))

  ##linking treatment labels to treatments
pheno_data <- data.frame(row.names = colnames(expression_matrix), treatment = colnames(expression_matrix))



