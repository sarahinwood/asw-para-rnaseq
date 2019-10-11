##run exposed deseq analysis until results deseq run on dds_group
##found code at https://rstudio-pubs-static.s3.amazonaws.com/187480_55b8817c975341c0bbc81934d773e768.html

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
##make list of genes from DNA viruses to subset results for plot below
dsDNA_viral_genes <- viral_annots[viral_annots$V3 == " dsDNA viruses, no RNA stage",]
dsDNA_list <- as.character(dsDNA_viral_genes$`#gene_id`)

##make a tidy results table from exposed deseq analysis, arranged by p value
res <- results(dds_group, tidy=TRUE, contrast=c("group", "Head_Exposed", "Head_Control")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
##take counts from dds, and normalise without outlier replacement
##add half a count because next step is to log2() and log2(0)=negative infinity
##now have log-transformed normalized count matrix, and need to transpose so sample names are rownames
##merge with colData(dds) where sample names are rownames
##which gives a wide dataframe with a row for each sample and column for all genes
##gather to have one row per sample per gene
tcounts <- t(log2((counts(dds_group[dsDNA_list, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds_group), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(dsDNA_list)+1):ncol(.))

##add annotation to tcounts table
annot_full <- viral_annots[,tstrsplit(V1, "Full=", keep=c(2)), by=`#gene_id`]
annot_brief <- viral_annots[,tstrsplit(V1, "^", fixed=TRUE, keep=c(1)), by="#gene_id"]
annot_brief$V2 <- paste(annot_brief$"V1", annot_brief$`#gene_id`, sep=" - ")


annot_tcounts <- merge(tcounts, annot_brief, by.x="gene", by.y="#gene_id")


##plot counts for dsDNA virus genes
ggplot(annot_tcounts, aes(x=Treatment, y=expression, colour=Tissue, shape=Wasp_Location)) + 
  geom_point() + 
  facet_wrap(~V2, scales = "free") + 
  labs(x="Treatment", 
       y="Expression (log normalized counts)",
       fill="(Tissue)")

##merge deseq results with viral annots
viral_deseq_res <- merge(ordered_res_group_table, viral_annots, by.x="rn", by.y="#gene_id", all.y = TRUE)
fwrite(viral_deseq_res, "output/viral/viral_deseq_res.csv")

##plot counts for a gene
plotCounts(dds_group, "TRINITY_DN1259_c0_g1", intgroup = c("group"), main="")