## Import libraries
library(rCausalMGM)
library(dplyr)
library(reshape2)
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(ggpubr)

set.seed(42)

## Set directory
setwd('Covid/')

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
alphas <- seq(0.01, 0.1, by=0.005)
for (a in alphas) {
  print(a)
  g <- fciMax(df = test_matrix, maxDiscrete = 5, fdr = T, alpha = a, verbose = F)
  saveGraph(g, paste0('Causal/Graphs/fci_' , a, '.txt'))
}

gr <- rCausalMGM::stars(df = test_matrix, method = 'fcim', maxDiscrete = 5, verbose = F)
gr
gr$edges
saveGraph(gr, 'Causal/Graphs/stars_fci.txt')

#g2 <- fciMax(df = test_matrix, maxDiscrete = 11, fdr = T, verbose = T)
#g2
#g2$edges
#saveGraph(g2, 'Causal/covid_fcicMax.txt')


# Prediction

model_data <- test_matrix
model_data$ICU <- ifelse(test_matrix$ICU=='Yes', 1, 0)
model_data$mechanical_ventilation <- ifelse(test_matrix$mechanical_ventilation=='yes', 1, 0)


#model_ <- lm(ICU ~ LF_38+LF_10+LF_34+LF_14+LF_37+LF_23+LF_3, data = test_matrix)
#model_ <- glm(ICU ~ LF_38+LF_10+LF_34+ventilator_free_days+LF_14+LF_37+LF_23+mechanical_ventilation+LF_3, data = model_data, family = "binomial")
model_ <- glm(ICU ~ LF_38+LF_10+LF_34+LF_14+LF_37+LF_23+LF_3, data = model_data, family = "binomial")

summary(model_)
prob=predict(model_,type=c("response"))
model_data$prob=prob
library(pROC)
g <- plot.roc(ICU ~ prob, data = model_data,
              print.auc=TRUE,
              ci=TRUE,
              thresholds='best',
              print.thres='best',
              main = "ICU Prediction AUC with Markov Blanket")
plot(g)    

# visualize latent factor
gr <- fciMax(df = test_matrix, maxDiscrete = 5, fdr = T, alpha = 0.05, verbose = F)

colnames(L)<-paste0('LF_',1:40)
D_df<-as.data.frame(L[,c(2, 5, 20, 36)])
D_df$Group = substr(rownames(L),1,1)
D_df$Group<-plyr::mapvalues(D_df$Group,
                            from =c("N","C"),to =c("Non-Covid","Covid")  )
D_df$Group<-as.factor(D_df$Group)
df_long <- melt(D_df, id.var = "Group")

my_comparisons <- list(c("Covid","Non-Covid"))
colnames(df_long)[2]<-'LF'
p<-ggplot(df_long,aes(Group,value,fill=Group)) + 
  geom_violin()  + theme_classic() +   
  ylab("Loading") + xlab("EBMF Latent Factors") +
  facet_wrap(~LF,ncol=5,scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')#,label = "p.signif")

cowplot::save_plot('Causal/Graphs/covid_violin_disease.png',p,base_height = 5,base_width = 8)

# ICU

colnames(L)<-paste0('LF_',1:40)
D_df<-as.data.frame(test_matrix[,c(3, 23)])
D_df$ICU = test_matrix$ICU
D_df$ICU<-as.factor(D_df$ICU)
df_long <- melt(D_df, id.var = "ICU")

my_comparisons <- list(c("Yes","No"))
colnames(df_long)[2]<-'LF'
p<-ggplot(df_long,aes(ICU,value,fill=ICU)) + 
  geom_violin()  + theme_classic() +   
  ylab("Loading") + xlab("EBMF Latent Factors") +
  facet_wrap(~LF,ncol=5,scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')#,label = "p.signif")

cowplot::save_plot('Causal/Graphs/covid_violin_icu.png',p,base_height = 5,base_width = 6)

# bio
lfs <- f.b$ldf$f

#lf1 <- as.data.frame(lfs[,1, drop=F])
write.csv(lfs, 'Causal/Biological/loadings.csv', row.names = T)

