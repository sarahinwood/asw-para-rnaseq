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

  ##Create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
  ##save dds as a file for import in clustering timecourse genes script
saveRDS(dds, file = "output/asw_timecourse/deseq2/dds.rds")

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
dds_abdo_res <- results(dds_abdo, alpha = 0.05)
 ##order based off padj
ordered_dds_abdo_res <- dds_abdo_res[order(dds_abdo_res$padj),]

  ##make list of sig genes
sig_genes <- subset(dds_abdo_res, padj < 0.05)
  ##make list of sig gene names
sig_gene_names <- row.names(sig_genes)
  ##save list of sig gene names
fwrite(data.table(sig_gene_names), "output/asw_timecourse/deseq2/timecourse_sig_gene_names.csv")

##plot counts for genes of interest, sub in name
plotCounts(dds_abdo, "TRINITY_DN13642_c0_g1", intgroup = c("Treatment", "Wasp_Location"))
