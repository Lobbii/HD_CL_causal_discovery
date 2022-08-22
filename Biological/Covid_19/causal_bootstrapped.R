library(rCausalMGM)
library(dplyr)
library(reshape2)
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(ggpubr)

setwd('Covid/')

set.seed(42)

f.b<-readRDS('d_f.b.rds')
L<-f.b$ldf$l
table(substr(rownames(L),1,1))

group<-substr(rownames(L),1,1)

colnames(L)<-paste0('LF_',1:40)
Cli_DF<-readRDS('Cli_DF_v2.rds')
index<-rownames(L) %in% rownames(Cli_DF)
L<-L[index,]

all(rownames(Cli_DF)==rownames(L))

test_matrix<-cbind(L,Cli_DF,stringsAsFactors=F)
test_matrix<-as.data.frame(test_matrix)

saveRDS(test_matrix, file = 'covid_data.RDS')

# FCI section
g.boot <- bootstrap(test_matrix, method = 'fci-max', alpha=0.05, ensembleMethod = 'majority', verbose = T)
g.boot
g.boot$edges

saveGraph(g.boot, 'boostrap_fci_stable.txt')
write.csv(g.boot$stabilities, 'METABRIC_stabilities.csv', row.names=T)