library(rCausalMGM)
library(dplyr)
library(reshape2)
library(ggplot2)
library(MetBrewer)
library(cowplot)
library(ggpubr)

L <- f$ldf$l
colnames(L)<-paste0('LF_',1:328)

# ER: 2
d_df<-as.data.frame(L[,2])
colnames(d_df)[1] <- 'LF_2'
d_df$ER_STATUS <- clin_data$ER_STATUS
df_long <- melt(d_df, id.var = "ER_STATUS")

my_comparisons <- list(c("Positive","Negative"))
colnames(df_long)[2]<-'LF'
p1<-ggplot(df_long,aes(ER_STATUS,value,fill=ER_STATUS)) + 
  geom_violin()  + theme_classic() +   
  ylab("Factor Values") + xlab("EBMF Latent Factors") +
  facet_wrap(~LF,ncol=5,scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')#,label = "p.signif")

cowplot::save_plot('Violin/metabric_er.png',p1,base_height = 5,base_width = 3)

# PR: 11, 27
d_df<-as.data.frame(L[,c(11, 27)])
d_df$PR_STATUS <- clin_data$PR_STATUS
df_long <- melt(d_df, id.var = "PR_STATUS")

my_comparisons <- list(c("Positive","Negative"))
colnames(df_long)[2]<-'LF'
p2<-ggplot(df_long,aes(PR_STATUS,value,fill=PR_STATUS)) + 
  geom_violin()  + theme_classic() +   
  ylab("Factor Values") + xlab("EBMF Latent Factors") +
  facet_wrap(~LF,ncol=5,scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')#,label = "p.signif")

cowplot::save_plot('Violin/metabric_pr.png',p2,base_height = 5,base_width = 6)

# Distant Relapse: 1, 10
d_df<-as.data.frame(L[,c(1, 10)])
d_df$Distant_relapse <- clin_data$Distant_relapse
df_long <- melt(d_df, id.var = "Distant_relapse")

my_comparisons <- list(c("Early","Mild","Late","Survivor"))
colnames(df_long)[2]<-'LF'
p3<-ggplot(df_long,aes(Distant_relapse,value,fill=Distant_relapse)) + 
  geom_violin()  + theme_classic() +   
  ylab("Factor Values") + xlab("EBMF Latent Factors") +
  facet_wrap(~LF,ncol=5,scales = "free_y") +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "red",
               aes(color = "Mean")) +
  scale_colour_manual(values = c("red"), # Colors
                      name = "") # Remove the legend title
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')#,label = "p.signif")


#cowplot::save_plot('Violin/metabric_er.png',p,base_height = 5,base_width = 6)

# Tumor Size: 59
d_df<-as.data.frame(L[,c(2)])
d_df$TUMOR_SIZE <- clin_data$TUMOR_SIZE
#df_long <- melt(d_df, id.var = "TUMOR_SIZE")

#my_comparisons <- list(c("Positive","Negative"))
colnames(d_df)[1]<-'LF_2'
p4 <- ggplot(d_df, aes(x = LF_2, y = TUMOR_SIZE)) +
  geom_point() + theme_classic()
#p4<-ggplot(df_long,aes(TUMOR_SIZE,value,fill=TUMOR_SIZE)) + 
 # geom_point() + theme_classic() +   
  #ylab("Loading") + xlab("EBMF Latent Factors") +
#  facet_wrap(~LF,ncol=5,scales = "free_y") #+
  #stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test')#,label = "p.signif")

#cowplot::save_plot('Violin/metabric_er.png',p,base_height = 5,base_width = 6)