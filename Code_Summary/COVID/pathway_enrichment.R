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

# loadings <- read.csv('loadings.csv', row.names=1)
setwd('/Users/mij51/Desktop/EBMF_10run/Covid/')


f.b<-readRDS('d_f.b.rds')
loadings<-f.b$ldf$f

dim(loadings)

sig <- c(3, 23, 2, 5, 20, 36)

set.seed(20220214)

plot.vec <- list()
count <- 1

gsea_list<-c()
for (idx in sig) {
    
    temp.load <- loadings[,idx]
    temp.load <- rank(temp.load, ties.method='random')
    names(temp.load) <- row.names(loadings)
    temp.load <- sort(temp.load, decreasing=T)
    ## temp.load <- order(temp.load, names(temp.load), decreasing=T)

    ## print(hist(temp.load))
 
    gsea.res <- gseGO(temp.load, ont='BP', OrgDb=org.Hs.eg.db, keyType='SYMBOL', minGSSize=20, maxGSSize=300, pvalueCutoff=0.05)
    print(gsea.res)
    gsea_list<-c(gsea_list,gsea.res)
    
    p <- dotplot(gsea.res, showCategory=20) + ggtitle(paste('LF', idx, 'Top 20 Enriched Gene Sets', sep=' ')) + theme_bw()
    pdf(paste('LF',idx,'_gsea_dotplot.pdf', sep=''), width=5, height=8)
    print(p)
    dev.off()
    plot.vec[[count]] <- p
    count <- count + 1
}






pdf('icu_gsea_dotplot.pdf', width=10, height=7)
plot_grid(plotlist=plot.vec[1:2], ncol=2)
dev.off()

pdf('disease_state_gsea_dotplot.pdf', width=20, height=7)
plot_grid(plotlist=plot.vec[3:7], ncol=4, nrow=1)
dev.off()
