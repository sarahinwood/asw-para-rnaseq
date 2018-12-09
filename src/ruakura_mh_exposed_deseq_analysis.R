##rukaura-mh-exposed asw vs control

library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")

gene2tx <- fread("data/Trinity.fasta.gene_trans_map", header = FALSE)
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
dds_ru_group <- copy(dds)
##create groupings of tissue+treatment
dds_ru_group$group <- factor(paste(dds$Tissue,dds$Treatment,dds$Wasp_Location,sep="_"))
##add group to design
design(dds_ru_group) <- ~group
##run deseq2 and generate results
dds_ru_group <- DESeq(dds_ru_group)

resultsNames(dds_ru_group)
##Make table of results for exposed vs control heads
res_ru_group <- results(dds_ru_group, contrast = c("group", "Head_Exposed_Ruakura", "Head_Control_Ruakura"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_ru_group <- res_ru_group[order(res_ru_group$padj),]
##Make data table and write to output
ordered_res_ru_group_table <- data.table(data.frame(ordered_res_ru_group), keep.rownames = TRUE)
fwrite(ordered_res_ru_group_table, "output/ru_exposed/deseq2/all_degs_ru_ex_v_ru_control.csv")
ordered_sig_res_ru_group_table <- subset(ordered_res_ru_group_table, padj < 0.05)
fwrite(ordered_sig_res_ru_group_table, "output/exposed/deseq2/ru_ex_sig_degs.csv", col.names = TRUE, row.names = FALSE)

plotCounts(dds_ru_group, "TRINITY_DN11440_c1_g1", intgroup = c("group"))

##volcano plot
EnhancedVolcano(ordered_res_ru_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)




