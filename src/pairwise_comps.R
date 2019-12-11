library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("Biostrings")
library("dplyr")
library("VennDiagram")

##ASW gene to trans map (edited to have ID's match concat transcriptome - ASW_TRINITY)
gene2tx <- fread("data/asw_transcriptome/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files - quantified against concat ASW-MH transcriptome
quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv")
setkey(sample_data, Sample_name)

##Create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
#saveRDS(dds, file = "output/asw_timecourse/deseq2/dds.rds")

##read in dds to save rerunning
#dds <- readRDS("output/asw_timecourse/deseq2/dds.rds")

##remove 8h NC2&3 - find sample number in table
colData(dds)
##remove 8hNC2 - 4th in table
dds <- dds[,-4]
##remove 8hNC3 - sample 5 once 8hNC2 removed
dds <- dds[,-5]
##check correct samples were removed
colData(dds)

##Select only abdomen samples
dds_abdo <- dds[,dds$Tissue == "Abdomen"]

##create dds object for group analysis
dds_group <- copy(dds)
##create groupings of tissue+treatment
dds_group$group <- factor(paste(dds$Tissue,dds$Treatment,sep="_"))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
#saveRDS(dds_group, file = "output/exposed/deseq2/dds_group.rds")

resultsNames(dds_group)

res_group <- results(dds_group, contrast = c("group", "Abdomen_m480", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)

fwrite(ordered_res_group_table, "output/exposed/deseq2/res_group.csv")

trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt")

degs_w_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)

plotCounts(dds_group, "ASW_TRINITY_DN9264_c0_g1", intgroup = c("group"), main="")

