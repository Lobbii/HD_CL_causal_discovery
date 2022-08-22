setwd('Covid')

count_mat<-read.delim('./GSE157103_genes.ec.tsv',row.names = 1)


library(edgeR)
d0 <- DGEList(count_mat)
d0 <- calcNormFactors(d0)
d0


cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

snames <- colnames(count_mat) # Sample names
snames

cultivar <- substr(snames, 1,  1) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar

group <- interaction(cultivar, time)
group

plotMDS(d, col = as.numeric(factor(cultivar)))

norm_d <- voom(d, plot = T)
norm_d<-t(norm_d$E)
saveRDS(norm_d,'norm_d.rds')


cutoff <- 3
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
norm_d3 <- voom(d, plot = T)
norm_d3<-t(norm_d3$E)
saveRDS(norm_d3,'norm_d3.rds')




# EBMF
setwd('/net/talisker/home/benos/mij51/EBMF/Covid/')
Y<-readRDS('norm_d3.rds')

# srun -p benos --ntasks-per-node=10   --pty /bin/bash -i

library(flashr)
f = flash(Y,Kmax = 1000)
saveRDS(f,'d3_f.rds')
f.b = flash(Y, f_init = f, backfit=TRUE, greedy=FALSE,Kmax = 500)
saveRDS(f.b,'d3_f.b.rds')







