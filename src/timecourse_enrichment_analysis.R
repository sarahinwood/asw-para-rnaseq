library('data.table')
library('fgsea')
library('ggplot2')

trinotate_report <- fread("data/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_pfam), unique(`#gene_id`)]
res_group <- fread("output/exposed/deseq2/res_group.csv")
res_timecourse <- fread("output/asw_timecourse/deseq2/timecourse_all_genes.csv")

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
setorder(res_timecourse, stat)
ranks <- res_timecourse[!is.na(stat), stat]
names(ranks) <- res_timecourse[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
##43 enriched GO terms
sum(sorted_fgsea_res$padj<0.05)
fwrite(sorted_fgsea_res, "output/asw_timecourse/fgsea/fgsea_timecourse_GOtermpfam_deseqstat_res.csv")

##read in file with functions added to GO terms when padj<0.1 (41 below 0.1)
annot_fgsea_res <- fread("output/asw_timecourse/fgsea/annot_fgsea_timecourse_GOtermpfam_deseqstat_res.csv")
##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_fgsea_res[annot_fgsea_res$pathway_kind=="biological process"]
cc_res <- annot_fgsea_res[annot_fgsea_res$pathway_kind=="cellular component"]
mf_res <- annot_fgsea_res[annot_fgsea_res$pathway_kind=="molecular function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(cc_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Cellular Component GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()

##find genes in GO:signal transduction and look at annots
signal_transduction_genes <- go_term_table[go_term_table$accessions == "GO:0007165"]
sig_trans_annots <- merge(x = signal_transduction_genes, y = trinotate_report, by.x = "gene_id", by.y="#gene_id", all.x = TRUE, all.y = FALSE)
fwrite(sig_trans_annots, "output/asw_timecourse/fgsea/signal_transduction_genes.csv")

####Core members that contribute to ES score (present in list before running sum reaches max.dev. from 0)
sig_trans_res <- fgsea_res[fgsea_res$pathway == "GO:0007165",]
sig_trans_leading_edge <- data.frame(sig_trans_res$leadingEdge)
setnames(sig_trans_leading_edge, old=c("c..TRINITY_DN338_c0_g1....TRINITY_DN8216_c0_g1....TRINITY_DN6471_c0_g1..."), new=c("gene_id"))
sig_trans_leading_annots <- merge(sig_trans_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(sig_trans_leading_annots, "output/asw_timecourse/fgsea/sig_trans_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0007165"]], ranks) + labs(title="signal transduction")

##find ALL genes in GO:regulation of transcription, DNA-templated
transcription_reg_genes <- go_term_table[go_term_table$accessions == "GO:0006355"]
transcr_reg_annots <- merge(x = transcription_reg_genes, y = trinotate_report, by.x = "gene_id", by.y="#gene_id", all.x = TRUE, all.y = FALSE)
fwrite(transcr_reg_annots, "output/asw_timecourse/fgsea/transcription_reg_genes.csv")

####find CORE members that contribute to ES score (present in list before running sum reaches max.dev. from 0)
trans_reg_res <- fgsea_res[fgsea_res$pathway == "GO:0006355",]
trans_reg_leading_edge <- data.frame(trans_reg_res$leadingEdge)
setnames(trans_reg_leading_edge, old=c("c..TRINITY_DN3621_c0_g2....TRINITY_DN4127_c4_g1....TRINITY_DN4199_c0_g1..."), new=c("gene_id"))
trans_leading_annots <- merge(trans_reg_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(trans_leading_annots, "output/asw_timecourse/fgsea/trans_reg_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0006355"]], ranks) + labs(title="regulation of transcription, DNA-templated")

####change to get leading edge genes of any GO pathway
carb_res <- fgsea_res[fgsea_res$pathway == "GO:0005975",]
carb_leading_edge <- data.frame(carb_res$leadingEdge)
setnames(carb_leading_edge, old=c("c..TRINITY_DN902_c0_g1....TRINITY_DN4004_c0_g1....TRINITY_DN11718_c0_g1..."), new=c("gene_id"))
carb_annots <- merge(carb_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(carb_annots, "output/asw_timecourse/fgsea/carb_met_leading_edge_annots.csv")
