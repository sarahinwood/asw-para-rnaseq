library('data.table')
library('fgsea')
library('ggplot2')
library('VennDiagram')

trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_pfam), unique(`#gene_id`)]
res_timecourse <- fread("output/asw_timecourse/deseq2/timecourse_all_genes.csv")

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
setorder(res_timecourse, stat)
ranks <- res_timecourse[!is.na(stat), stat]
names(ranks) <- res_timecourse[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
##43 enriched GO terms
sum(sorted_fgsea_res$padj<0.05)
fwrite(sorted_fgsea_res, "output/asw_timecourse/fgsea/fgsea_timecourse_GOtermpfam_deseqstat_res.csv")

sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/asw_timecourse/fgsea/sig_annot_fgsea_pfam.csv")

##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_sig_fgsea[annot_sig_fgsea$pathway_kind=="biological_process"]
cc_res <- annot_sig_fgsea[annot_sig_fgsea$pathway_kind=="cellular_component"]
mf_res <- annot_sig_fgsea[annot_sig_fgsea$pathway_kind=="molecular_function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(bp_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Biological Process GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()+
  theme(axis.text.y = element_text(size=15))

####Core members that contribute to ES score (present in list before running sum reaches max.dev. from 0)
sig_trans_res <- fgsea_res[fgsea_res$pathway == "GO:0007165",]
sig_trans_leading_edge <- data.frame(sig_trans_res$leadingEdge)
setnames(sig_trans_leading_edge, old=c("c..TRINITY_DN436_c0_g1....TRINITY_DN4549_c1_g1....TRINITY_DN8214_c0_g1..."), new=c("gene_id"))
sig_trans_leading_annots <- merge(sig_trans_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(sig_trans_leading_annots, "output/asw_timecourse/fgsea/sig_trans_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0007165"]], ranks) + labs(title="signal transduction")
##compare to prev.results
old_sig_trans <- fread("nf_output/asw_timecourse/fgsea/dedup_sig_trans_leading_edge_annots.csv")
old_st_id <- old_sig_trans$gene_id
new_st_id <- sig_trans_leading_edge$gene_id
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("NF Signal Transduction"=old_st_id, "F Signal Transduction"=new_st_id), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

####find CORE members that contribute to ES score (present in list before running sum reaches max.dev. from 0)
trans_reg_res <- fgsea_res[fgsea_res$pathway == "GO:0006355",]
trans_reg_leading_edge <- data.frame(trans_reg_res$leadingEdge)
setnames(trans_reg_leading_edge, old=c("c..TRINITY_DN3621_c0_g2....TRINITY_DN4127_c4_g1....TRINITY_DN5168_c0_g1..."), new=c("gene_id"))
trans_leading_annots <- merge(trans_reg_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(trans_leading_annots, "output/asw_timecourse/fgsea/trans_reg_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0006355"]], ranks) + labs(title="regulation of transcription, DNA-templated")
##compare to prev.results
old_trans_reg <- fread("nf_output/asw_timecourse/fgsea/dedup_trans_reg_leading_edge_annots.csv")
old_tr_id <- old_trans_reg$gene_id
new_tr_id <- trans_reg_leading_edge$gene_id
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("NF Transcription Regulation"=old_tr_id, "F Transcription Regulation"=new_tr_id), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)


steroid_res <- fgsea_res[fgsea_res$pathway == "GO:0006694",]
steroid_leading_edge <- data.frame(steroid_res$leadingEdge)
setnames(steroid_leading_edge, old=c("c..TRINITY_DN1464_c0_g2....TRINITY_DN30690_c0_g1....TRINITY_DN4292_c0_g1..."), new=c("gene_id"))
steroid_leading_annots <- merge(steroid_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(steroid_leading_annots, "output/asw_timecourse/fgsea/steroid_leading_edge_annots.csv")

transmembrane_res <- fgsea_res[fgsea_res$pathway == "GO:0055085",]
transmembrane_leading_edge <- data.frame(transmembrane_res$leadingEdge)
setnames(transmembrane_leading_edge, old=c("c..TRINITY_DN5359_c0_g1....TRINITY_DN21633_c0_g1....TRINITY_DN4949_c0_g1..."), new=c("gene_id"))
transmembrane_leading_annots <- merge(transmembrane_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(transmembrane_leading_annots, "output/asw_timecourse/fgsea/transmembrane_leading_edge_annots.csv")
