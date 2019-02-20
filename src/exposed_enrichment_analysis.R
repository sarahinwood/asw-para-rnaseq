library('data.table')
library('fgsea')
library('ggplot2')

trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_pfam), unique(`#gene_id`)]
res_group <- fread("output/exposed/deseq2/res_group.csv")

go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
term_to_gene <- go_term_table[,list(list(gene_id)), by=accessions]
pathways <- term_to_gene[,V1]
names(pathways) <- term_to_gene[,accessions]

##use stat column from deseq results to rank genes (can change if wanted)
setorder(res_group, stat)
ranks <- res_group[!is.na(stat), stat]
names(ranks) <- res_group[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sum(sorted_fgsea_res$padj<0.05)
fwrite(sorted_fgsea_res, "output/exposed/fgsea/fgsea_exposed_GOtermpfam_deseqstat_res.csv")

##plot enrichment of GO term
plotEnrichment(pathways[["GO:0098586"]], ranks) + labs(title="cellular response to virus")

##subset into only sig terms and merge w/annotations
sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/exposed/fgsea/sig_annot_fgsea_pfam.csv")
##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="biological_process"]
cc_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="cellular_component"]
mf_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="molecular_function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(mf_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Molecular Function GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()
