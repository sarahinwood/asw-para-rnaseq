##run exposed deseq analysis until results data.table

library("tximport")
library("data.table")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("tidyr")

##read in dedup viral annots
dd_viral <- fread("data/viral/dedup_viral_annots_trinotate.csv")
##split table to include columns of annotation from annotation taxa etc
viral_annots <- dd_viral[,tstrsplit(sprot_Top_BLASTX_hit, ";", fixed=TRUE, keep=c(1,2,3,4)), by = `#gene_id`]

dsDNA_viral_genes <- viral_annots[viral_annots$V3 == " dsDNA viruses, no RNA stage",]
dsDNA_list <- dsDNA_viral_genes$`#gene_id`


##tidy results table from exposed deseq analysis
res <- results(dds_group, tidy=TRUE, contrast=c("group", "Head_Exposed", "Head_Control")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
goi <- res[res$row == dsDNA_list,]
##make table of counts for all dsDNA viral genes of interest
tcounts <- t(log2((counts(dds_group[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds_group), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
##add annotation to tcounts table
annot_full <- viral_annots[,tstrsplit(V1, "Full=", keep=c(2)), by=`#gene_id`]
annot_tcounts <- merge(tcounts, annot_full, by.x="gene", by.y="#gene_id")
##plot counts for dsDNA virus genes
ggplot(annot_tcounts, aes(x=Treatment, y=expression, fill=Tissue)) + 
  geom_boxplot() + 
  facet_wrap(~V1) + 
  labs(x="Treatment", 
       y="Expression (log normalized counts)", 
       fill="(Tissue)")

##merge deseq results with viral annots
viral_deseq_res <- merge(ordered_res_group_table, viral_annots, by.x="rn", by.y="#gene_id", all.y = TRUE)
fwrite(viral_deseq_res, "output/viral/viral_deseq_res.csv")

##plot counts for a gene
plotCounts(dds_group, "TRINITY_DN1259_c0_g1", intgroup = c("group"), main="")