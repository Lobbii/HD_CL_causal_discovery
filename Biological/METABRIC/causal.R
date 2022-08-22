library(rCausalMGM)
library(dplyr)

#load('EBMF.Rdata')

f <- readRDS('f.gb_1000.rds')

disc_data <- readRDS('discrete_data.rds')

pat_info <- readRDS('Patient_info.rds')
pat_info$LYMPH_NODE_STATUS <- '4toX'
pat_info$LYMPH_NODE_STATUS[pat_info$LYMPH_NODES_EXAMINED_POSITIVE < 4] <- '1to3'
pat_info$LYMPH_NODE_STATUS[pat_info$LYMPH_NODES_EXAMINED_POSITIVE == 0] <- 'NodeNegative'
pat_info$LYMPH_NODE_STATUS[is.na(pat_info$LYMPH_NODES_EXAMINED_POSITIVE)] <- NA
rownames(pat_info) <- pat_info$PATIENT_ID

samp_info <- readRDS('Sample_info.rds')
rownames(samp_info) <- samp_info$PATIENT_ID

recu_info <- readRDS('recurrent_info.rds')
rownames(recu_info) <- recu_info$METABRIC.ID

clin_data <- pat_info %>% select(CELLULARITY, LYMPH_NODE_STATUS, INFERRED_MENOPAUSAL_STATE,
                                 INTCLUST, AGE_AT_DIAGNOSIS, CLAUDIN_SUBTYPE, )

clin_data <- cbind(clin_data, samp_info %>% select(CANCER_TYPE_DETAILED, ER_STATUS, HER2_STATUS,
                                                   GRADE, PR_STATUS, TUMOR_SIZE))

clin_data <- cbind(clin_data, recu_info %>% select(Distant_relapse, Loco_regional_relapse))

## clin_data$BREAST_SURGERY[clin_data$BREAST_SURGERY==''] <- NA
clin_data$CELLULARITY[clin_data$CELLULARITY==''] <- NA
clin_data$CANCER_TYPE_DETAILED[clin_data$CANCER_TYPE_DETAILED=='Breast'] <- NA
clin_data$AGE_AT_DIAGNOSIS <- as.numeric(clin_data$AGE_AT_DIAGNOSIS)
clin_data$TUMOR_SIZE <- log2(clin_data$TUMOR_SIZE)

CC <- rowSums(is.na(clin_data))==0

latents <- f$ldf$l
colnames(latents) <- paste('Z', seq(1, ncol(latents)), sep='')

loadings <- f$ldf$f
colnames(loadings) <- paste('Z', seq(1, ncol(latents)), sep='')

data <- cbind(latents, clin_data)[CC,]
dim(data)
head(data)

write.csv(data, 'joint_data.csv', row.names=T)

set.seed(42)
ig <- steps(data, maxDiscrete=11, g=0.1,
            lambda=10^(seq(log10(0.5)-1, log10(0.5), length.out=50)), verbose=T)
ig
ig$edges

saveGraph(ig, 'EBMF_bf_mgm.txt')

ig <- loadGraph('EBMF_bf_mgm.txt')
ig

## ggplot(data, aes(x=Z48, y=as.numeric(AGE_AT_DIAGNOSIS))) +
##     geom_point() +
##     geom_smooth(method='lm') + 
##     scale_color_continuous()

set.seed(42)
g <- stars(data, method='pcm', initialGraph=ig, g=0.1,
           params=c(1e-5, 2e-5, 3.5e-5, 5e-5, 7.5e-5, 1e-4, 5e-4, 1e-3, 5e-3),
           maxDiscrete=11, adjacency=F, verbose=T)
## g <- pcMax(data, initialGraph=ig, maxDiscrete=11, alpha=3e-5, verbose=T)
g
g$edges

saveGraph(g, 'EBMF_bf_mgmpcmax.txt')

## ggplot(data, aes(x=Z1, y=Z2, color=factor(CLAUDIN_SUBTYPE))) +
##     geom_point()

## emb <- umap(data[,setdiff(g$markov.blankets$CLAUDIN_SUBTYPE, colnames(clin_data))], n_neighbors=30)

## ggplot(cbind(umap1=emb$layout[,1], umap2=emb$layout[,2], data),
##        aes(x=umap1, y=umap2, color=factor(CLAUDIN_SUBTYPE))) +
##     geom_point()
