library("tximport")
library("data.table")
library("VennDiagram")

nf_sig_degs <- fread("output/exposed/nf_deseq2/exposed_analysis_sig_degs.csv")
f_sig_degs <- fread("output/exposed/deseq2/exposed_analysis_sig_degs.csv")

nf_sig_ids <- nf_sig_degs$rn
f_sig_ids <- f_sig_degs$rn
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Non-Filtered DEGs"=nf_sig_ids, "Filtered DEGs"=f_sig_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)
