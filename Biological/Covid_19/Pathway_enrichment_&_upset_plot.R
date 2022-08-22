library(edgeR)
library(ggpubr)
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
library(UpSetR)

# load data 
count_mat<-read.delim('./GSE157103_genes.ec.tsv',row.names = 1)
Cli_DF<-readRDS('Cli_DF.rds')

d0 <- DGEList(count_mat)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 


Disease_State <- Cli_DF$disease_state
Disease_State<-plyr::mapvalues(Disease_State,
                               from =c("COVID-19","non-COVID-19"),to =c("Yes","No")  )
ICU<-Cli_DF$ICU
norm_d <- voom(d,plot = T)


# Perfrom DGE analysis
mm <- model.matrix(~0 + Disease_State)
fit <- lmFit(norm_d,mm)
head(coef(fit))

contr <- makeContrasts(Disease_StateYes - Disease_StateNo, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "Covid_YesvsNo_DGE.csv", row.names = F, sep = "\t", quote = F)

COVID_sig_gene_tbl<-top.table[top.table$adj.P.Val<0.05,]


mm <- model.matrix(~0 + ICU)
fit <- lmFit(norm_d,mm)
head(coef(fit))

contr <- makeContrasts(ICUYes - ICUNo, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "ICU_YesvsNo_DGE.csv", row.names = F, sep = "\t", quote = F)
ICU_sig_gene_tbl<-top.table[top.table$adj.P.Val<0.05,]




# Perfrom GSEA analysis based on DGE result
covid_T<-COVID_sig_gene_tbl$t
names(covid_T)<-COVID_sig_gene_tbl$Gene
covid_T <- rank(covid_T, ties.method='random')
covid_T<-sort(covid_T,decreasing = T)

gsea.res_covid <- gseGO(covid_T, ont='BP', OrgDb=org.Hs.eg.db, keyType='SYMBOL', minGSSize=20, maxGSSize=300, pvalueCutoff=0.05)
P_covid<-dotplot(gsea.res_covid, showCategory=20) + ggtitle(paste('COVID19',  'Top 10 Enriched Gene Sets', sep=' ')) + theme_bw()
pdf(paste('DGE_COVID_gsea_dotplot.pdf', sep=''), width=8, height=8)
print(P_covid)
dev.off()

ICU_T<-ICU_sig_gene_tbl$t
names(ICU_T)<-ICU_sig_gene_tbl$Gene
ICU_T <- rank(ICU_T, ties.method='random')
ICU_T<-sort(ICU_T,decreasing = T)

gsea.res_ICU <- gseGO(ICU_T, ont='BP', OrgDb=org.Hs.eg.db, keyType='SYMBOL', minGSSize=20, maxGSSize=300, pvalueCutoff=0.05)
P_ICU<-dotplot(gsea.res_ICU, showCategory=10) + ggtitle(paste('ICU',  'Top 10 Enriched Gene Sets', sep=' ')) + theme_bw()
pdf(paste('DGE_ICU_gsea_dotplot.pdf', sep=''), width=8, height=8)
print(P_ICU)
dev.off()


ggarrange(P_covid,P_ICU,ncol = 2)

saveRDS(gsea.res_ICU,'gsea.res_ICU.rds')
saveRDS(gsea.res_covid,'gsea.res_covid.rds')



# Perfrom GSEA analysis based on EBMF latent factors
f.b<-readRDS('d_f.b.rds')
loadings<-f.b$ldf$f

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



# generate Upser Plot
list_COVID<-list(LF2 =gsea_list[[3]]@result$Description,LF5 =gsea_list[[4]]@result$Description,
                 LF20 =gsea_list[[5]]@result$Description,LF36 =gsea_list[[6]]@result$Description,
                 DGE = gsea.res_covid@result$Description)
P_covid<-upset(fromList(list_COVID), order.by = "freq",sets=names(list_COVID),keep.order = TRUE)
pdf(paste('COVID_gsea_UpSet.pdf', sep=''), width=5, height=5)
print(P_covid)
dev.off()


list_ICU<-list(LF3 =gsea_list[[1]]@result$Description, LF23 =gsea_list[[2]]@result$Description,
               DGE = gsea.res_ICU@result$Description)

P_ICU<-upset(fromList(list_ICU),order.by = "freq",sets=names(list_ICU),keep.order = TRUE)
pdf(paste('ICU_gsea_UpSet.pdf', sep=''), width=5, height=5)
print(P_ICU)
dev.off()
