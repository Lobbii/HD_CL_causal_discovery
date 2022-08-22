library(dplyr)
library(stringr)
library(fgsea)
library(limma)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(cowplot)

#loadings <- read.csv('loadings.csv', row.names=1)

dim(loadings)

sig <- c(27, 11, 2, 59, 10, 1)
sig <- sort(sig)

set.seed(20220214)

plot.vec <- list()
count <- 1

gsea.res <- list()

for (idx in sig) {
  print(idx)
  
  temp.load <- loadings[,idx]
  #temp.load <- rank(temp.load, ties.method='random')
  names(temp.load) <- row.names(loadings)
  temp.load <- sort(temp.load, decreasing=T)
  #temp.load <- order(temp.load, names(temp.load), decreasing=T)
  
  #print(hist(temp.load))
  
  gsea.res[[count]] <- gseGO(temp.load, ont='BP', OrgDb=org.Hs.eg.db, keyType='SYMBOL', minGSSize=20, maxGSSize=300, pvalueCutoff=0.05)
  print(gsea.res[[count]])
  
  #p <- dotplot(gsea.res, showCategory=10) + ggtitle(paste('LF', idx, 'Top 10 Enriched Gene Sets', sep=' ')) + theme_bw() # + theme(text = element_text(size = 20))
  p <- dotplot(gsea.res[[count]], showCategory=10) + ggtitle(paste('LF', idx, sep=' ')) + theme_bw() #+ theme(text = element_text(size = 20))
  #pdf(paste('LF',idx,'_gsea_dotplot.pdf', sep=''), width=5, height=8)
  #print(p)
  #dev.off()
  plot.vec[[count]] <- p
  count <- count + 1
}

count <- 1
for (idx in sig) {
  p <- dotplot(gsea.res[[count]], showCategory=20, label_format = 20, font.size=14) + ggtitle(paste('LF', idx, sep=' '))# + theme_bw()
  plot.vec[[count]] <- p
  count <- count + 1
}

#plot_grid(plotlist=plot.vec[], ncol=3, nrow=2)

pdf('metabric_gsea_dotplot_ref.pdf', width=25, height=25)
plot_grid(plotlist=plot.vec[], ncol=3, nrow=2)
dev.off()

