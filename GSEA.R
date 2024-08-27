library(clusterProfiler)
library(dplyr)
library(enrichplot)
DESeq2 <- DESeq2[DESeq2$padj<0.05,]
DESeq2 <- arrange(DESeq2,desc(log2FoldChange))
##genelist
gene_list <- DESeq2$log2FoldChange
names(gene_list) <- rownames(DESeq2)
gene_list[1:6]
h.all_gmt <- read.gmt("hall.gmt")
h.all_res <- GSEA(gene_list,TERM2GENE = h.all_gmt)
head(h.all_res@result[3:10])
hallmark_all <- h.all_res@result
library(GseaVis)
library(jjAnno)
# change line color
genesets <- h.all_res@result$ID[c(1)]
library(ggsci)
gseaNb(object = h.all_res,
       geneSetID = genesets,
       subPlot = 2,
       rmHt = T,
       termWidth = 35,
       #legend.position = c(0.8,0.8),
       curveCol = pal_npg()(7))
gseaNb(object = h.all_res,
       geneSetID = genesets,
       newGsea = T,
       rmHt = T,
       addPval = T,
       pvalX = 0.9,
       pvalY = 0.6)