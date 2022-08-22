setwd('')
RNA_data_filtered<-readRDS('./LOVE/RNA_data_filtered.rds') 
metabricdata<-readRDS('./LOVE/metabricdata.rds')
Patient_info<-readRDS('./LOVE/Patient_info.rds')

# metabric RNA normalization
RNA_data_filtered_zscore<-BBmisc::normalize(RNA_data_filtered,margin = 2)


# Run EBMF on metabric dataset
library(flashr)
f = flash(as.matrix(RNA_data_filtered_zscore), Kmax = 2000)
ldf = flash_get_ldf(f)
saveRDS(f,'f.rds')

f.b = flash(as.matrix(RNA_data_filtered_zscore), f_init = f, backfit=TRUE, greedy=FALSE, Kmax = 2000)
saveRDS(f.b,'f.b.rds')


f<-readRDS('f.b.rds')
ldf = flash_get_ldf(f)


library(MASS)
library(ggfortify)

plot(flash_get_fitted_values(f), as.matrix(RNA_data_filtered_zscore),
     main="compare fitted values with true LF'")

brca_L<-flashr::flash_sampler(as.matrix(brca_data_filtered_zscore), f)

# use z-score normalized count matrix 
L_METABRIC_zscore<-ginv(t(ldf$f)  %*% ldf$f) %*% t(ldf$f)%*%t(as.matrix(RNA_data_filtered_zscore))
L_TCGA_zscore<-ginv(t(ldf$f)  %*% ldf$f) %*% t(ldf$f)%*%t(as.matrix(brca_data_filtered_zscore))

# generate PCA plot
pca_res <- prcomp(t(Z_METABRIC_zscore), scale. = TRUE)
autoplot(pca_res, data = Patient_info, colour = 'CLAUDIN_SUBTYPE')+  theme_classic()
pca_res_TCGA <- prcomp(t(Z_TCGA_zscore), scale. = TRUE)
autoplot(pca_res_TCGA, data = brca_clinical, colour = 'subtype')+  theme_classic()

z_merge_zscore<- cbind(Z_METABRIC_zscore,Z_TCGA_zscore)
pca_res_merge_zscore <- prcomp(t(z_merge_zscore), scale. = TRUE)
autoplot(pca_res_merge_zscore,  data = subtuype_info, colour = 'subtype',shape='batch')+ 
  theme_classic()

