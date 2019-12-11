library("data.table")
library("VennDiagram")
library("DESeq2")

sample_names_pot_fixed <- fread("output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")
sample_names_pot_fixed$rn <- tstrsplit(sample_names_pot_fixed$rn, "ASW_", keep=c(2))
sample_names_pot_fixed_ids <- sample_names_pot_fixed$rn

sample_names_pot_swapped <- fread("samplelabelspotswapped_output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")
sample_names_pot_swapped$rn <- tstrsplit(sample_names_pot_swapped$rn, "ASW_", keep=c(2))
sample_names_pot_swapped_ids <- sample_names_pot_swapped$rn

short_tc <- fread("short_tc_output/asw_timecourse/deseq2/timecourse_analysis_sig_degs.csv")
short_tc_ids <- short_tc$rn

vd <- venn.diagram(x = list("sample_names_pot_fixed"=sample_names_pot_fixed_ids, "sample_names_pot_swapped"=sample_names_pot_swapped_ids, "short_tc"=short_tc_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

##comp swapped 8h to 8hNC
dds <- readRDS("output/asw_timecourse/deseq2/dds.rds")
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

##Make table of results for exposed vs control heads
res_group <- results(dds_group, contrast = c("group", "Abdomen_m480", "Abdomen_Control"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
#fwrite(ordered_res_group_table, "output/exposed/deseq2/res_group.csv")
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.1)
###none of the sig genes here are also sig in short tc
###if I try cluster the sample name swapped data it has completely different cluster patterns