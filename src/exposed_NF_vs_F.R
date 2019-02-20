library("tximport")
library("data.table")
library("VennDiagram")

nf_sig_degs <- fread("output/exposed/nf_deseq2/dedup_nf_only_deg_annots.csv")
quasi_sig_degs <- fread("output/exposed/deseq2/exposed_analysis_sig_degs.csv")
f_sig_degs <- fread("output/exposed/f_deseq2/exposed_analysis_sig_degs.csv")
annot_f_sig_degs <- fread("output/exposed/deseq2/degs_trinotate_blastx_annots.csv")
asw_trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt", na.strings = ".")
asw_gene_ids <- asw_trinotate_report[!is.na(gene_ontology_pfam), unique(`#gene_id`)]

nf_sig_ids <- nf_sig_degs$rn
f_sig_ids <- f_sig_degs$rn

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Non-Filtered DEGs"=nf_sig_ids, "Filtered DEGs"=f_sig_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)

nf_only_ids <- data.table(setdiff(nf_sig_ids, f_sig_ids))
nf_only_degs <- merge(nf_only_ids, nf_sig_degs, by.x="V1", by.y="rn", all.x=TRUE)
nf_only_annots <- merge(nf_only_degs, asw_trinotate_report, by.x="V1", by.y="#gene_id", all.x=TRUE)
fwrite(nf_only_annots, "output/exposed/nf_deseq2/nf_only_deg_annots.csv")

f_only_ids <- data.table(setdiff(f_sig_ids, nf_sig_ids))
f_only_annots <- merge(f_only_ids, annot_f_sig_degs, by.x="V1", by.y="#gene_id", all.x=TRUE)
fwrite(f_only_annots, "output/exposed/deseq2/f_only_degs.csv")

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}
##Extract ASW GO terms
go_term_list <- lapply(asw_gene_ids, EXTRACT_GO_TERMS, trinotate_report=asw_trinotate_report)
go_term_table <- rbindlist(go_term_list)
##asw nf GO terms for VD
asw_nf_ids<-data.table(asw_nf_deg_dedup$`#gene_id`)
asw_nf_go_terms<-merge(asw_nf_ids, go_term_table, by.x="V1", by.y="gene_id", all.x=TRUE)
asw_nf_go_vd <- na.omit(unique(asw_nf_go_terms$accessions))
##asw f GO terms for VD
asw_f_ids<-data.table(asw_f_deg_dedup$`#gene_id`)
asw_f_go_terms<-merge(asw_f_ids, go_term_table, by.x="V1", by.y="gene_id", all.x=TRUE)
asw_f_go_vd <- na.omit(unique(asw_f_go_terms$accessions))

##extract MH GO terms
mh_trinotate <- fread("data/mh_trinotate_annotation_report.txt", na.strings = ".")
mh_gene_ids <- mh_trinotate[!is.na(gene_ontology_pfam), unique(`#gene_id`)]
mh_go_term_list <- lapply(mh_gene_ids, EXTRACT_GO_TERMS, trinotate_report=mh_trinotate)
mh_go_term_table <- rbindlist(mh_go_term_list)
##mh GO terms for VD
mh_ids <- data.table(mh_sig_degs$`#gene_id`)
mh_go_terms<-merge(mh_ids, mh_go_term_table, by.x="V1", by.y="gene_id", all.x=TRUE)
mh_go_vd <- na.omit(unique(mh_go_terms$accessions))

Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Unique MH DEG GO Terms"=mh_go_vd, "Unique ASW NF DEG GO Terms"=asw_nf_go_vd, "Unique ASW F DEG GO Terms"=asw_f_go_vd), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE, fill=Set1)
grid.newpage()
grid.draw(vd)
##get all that overlap between mh and asw nf
mh_asw_nf_overlap <- intersect(mh_go_vd, asw_nf_go_vd)
##remove those that are in asw f
mh_aswnf_unique_overlap <- data.table(setdiff(mh_asw_nf_overlap, asw_f_go_vd))
fwrite(mh_aswnf_unique_overlap, "output/exposed/nf_deseq2/nf_mh_overlap_GO_terms.csv")
