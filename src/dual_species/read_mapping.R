library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

##asw reads mapped
dds_concat_group_asw <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
counts_table_asw <- (data.table(counts(dds_concat_group_asw)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("Sample_name", "readpairs_mapped_ASW"))

##mh reads mapped
dds_concat_group_mh <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
counts_table_mh <- (data.table(counts(dds_concat_group_mh)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("Sample_name", "readpairs_mapped_MH"))

##both reads mapped
read_mapping <- merge(counts_colSums_asw, counts_colSums_mh)
bbduk_reads_out <- fread("output/bbduk_trim/bbduk_reads_out.csv")
full_read_mapping <- merge(read_mapping, bbduk_reads_out, by="Sample_name")
full_read_mapping$bbduk_halved <- (full_read_mapping$bbduk_reads_out)/2

##mapping %s
full_read_mapping$total_mapped_reads <- (full_read_mapping$readpairs_mapped_MH + full_read_mapping$readpairs_mapped_ASW)
full_read_mapping$`overall_mapping_%` <- (full_read_mapping$total_mapped_reads)/(full_read_mapping$bbduk_halved)*100
full_read_mapping$`%_ofmapped_ASW` <- (full_read_mapping$readpairs_mapped_ASW/full_read_mapping$total_mapped_reads)*100
full_read_mapping$`%_ofmapped_MH` <- (full_read_mapping$readpairs_mapped_MH/full_read_mapping$total_mapped_reads)*100
full_read_mapping$`%_ofall_ASW` <- (full_read_mapping$readpairs_mapped_ASW/full_read_mapping$bbduk_halved)*100
full_read_mapping$`%_ofall_MH` <- (full_read_mapping$readpairs_mapped_MH/full_read_mapping$bbduk_halved)*100
fwrite(full_read_mapping, "output/deseq2/read_mapping.csv")

##plot or calc means for mh mapping between para and unpara
ggplot(full_read_mapping, aes(x=Parasitism_PCR, y=`%_ofall_MH`))+
  geom_boxplot(aes(fill=Parasitism_PCR), alpha=0.9, show.legend = FALSE)+
  theme_bw()+
  scale_fill_manual(values=c("#440154FF", "#31688EFF", "#FDE725FF", "#35B779FF"))+
  xlab("Parasitism status")+
  scale_colour_viridis_d()+
  ylab(expression(paste("% reads mapping to ", italic("M. hyperodae"), " transcriptome")))

##head or abdo only
head<-subset(full_read_mapping, full_read_mapping$Tissue == "Head")
abdo<-subset(full_read_mapping, full_read_mapping$Tissue == "Abdomen")
##plot or calc means for mh mapping between para and unpara
ggplot(abdo, aes(x=Parasitism_PCR, y=`%_ofall_MH`))+
  geom_boxplot(aes(fill=Parasitism_PCR), alpha=0.9, show.legend = FALSE)+
  theme_bw()+
  scale_fill_manual(values=c("#440154FF", "#FDE725FF"))+
  xlab("Parasitism status")+
  scale_colour_viridis_d()+
  ylab(expression(paste("% reads mapping to ", italic("M. hyperodae"), " transcriptome")))
