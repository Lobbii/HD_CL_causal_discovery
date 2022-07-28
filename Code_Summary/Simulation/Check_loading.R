library(flashr)
library(rCausalMGM)
library(dplyr)
library(grid)
library(gridExtra)


info<-matrix(nrow = 5,ncol = 14)
colnames(info)<-c("Source_edge","EBMF_edge","EBMF_TP","EBMF_FP","EBMF_FN",
                  "FCI_edges","Confounded_edges","FCI_TP","FCI_FP","FCI_FN",
                  "SHD_EBMF","SHD_FCI",
                  "NON <-> edge","NON <-> TP")

filenames<-c("simulation_1","simulation_2","simulation_3","simulation_4","simulation_5")

source('R/Mapping_loading.R')
source('R/generate_adj_matrix.R')
source('R/read_causality_graph.R')

for (i in 1:length(filenames)){
  
  file_name<-filenames[i]
  simulation<-readRDS(paste0(file_name,".rds")) 
  dir.create(paste0("./",file_name))
  Y<-simulation$simulated_data
  
  # perform EBMF 
  f = flash(Y,Kmax = 500)
  saveRDS(f,paste0("./",file_name,"/f.rds"))
  f.b = flash(Y, f_init = f, backfit=TRUE, greedy=FALSE)
  saveRDS(f.b,paste0("./",file_name,"/f.b.rds"))
  
  f.b<-readRDS(paste0("./",file_name,"/f.b.rds"))
  L<-f.b$ldf$f
  colnames(L)<-paste0("EBMF_",seq(ncol(L)))
  LL<-t(simulation$loading)
  colnames(LL)<-paste0("source_",seq(ncol(LL)))
  Mer_mat<-cbind(L,LL)
  cor_mat<-cor(Mer_mat)
  dir.create(paste0("./",file_name,"/f_b"))
  Loading_mat<-Mapping_Loading(L,LL,paste0("./",file_name,"/f_b"))
  
  
  L_norm<-Loading_mat$mapped_matrix/max(abs(Loading_mat$mapped_matrix))
  LL_norm<-LL/max(abs(LL))
  colnames(L_norm)<-paste0("EBMF_",seq(ncol(L_norm)))
  
  p1=pheatmap::pheatmap(L_norm ,cluster_rows = F,cluster_cols = F)
  p2=pheatmap::pheatmap(LL_norm,cluster_rows = F,cluster_cols = F)
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  g <- grid.arrange(grobs=plot_list, ncol=2)
  ggplot2::ggsave(file=paste0("./",file_name,"/LF_heatmap.png"), g) #saves g
  
  # infer edges
  L_norm<-abs(L_norm)
  L_norm[L_norm<=0.2]=0
  L_norm[L_norm>0.2]=1
  LL_norm<-abs(LL_norm)
  LL_norm[LL_norm<=0.2]=0
  LL_norm[LL_norm>0.2]=1
  p1=pheatmap::pheatmap(L_norm ,cluster_rows = F,cluster_cols = F)
  p2=pheatmap::pheatmap(LL_norm,cluster_rows = F,cluster_cols = F)
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  g <- grid.arrange(grobs=plot_list, ncol=2)
  ggplot2::ggsave(file=paste0("./",file_name,"/LF_heatmap_edge.png"), g) #saves g
  
  # adj matrix and SHD 
  L_norm_adj<-gen_adj_mat(L_norm)
  LL_norm_adj<-gen_adj_mat(LL_norm)
  info[i,11]<-sum(abs(L_norm_adj-LL_norm_adj)/2)
  
  # perfrom FCi
  dir.create(paste0("./",file_name,"/fci"))
  test_matrix<-as.data.frame(Y)
  g<-fciStable(test_matrix,alpha = 0.05,fdr = T,verbose = T)
  saveGraph(g, paste0("./",file_name,"/fci/mgm_fci.txt"))
  
  
  net_FCI<-read_causality_graph(paste0("./",file_name,"/fci/mgm_fci.txt"))
  info[i,6]<-nrow(net_FCI$edges) 
  Confounded_Edges<-net_FCI$edges[net_FCI$edges[,3]=="<->",] # "o->" "o-o" "<->" "-->"
  info[i,7]<-nrow(Confounded_Edges)
  
  FCI_adj<-gene_adj_mat_from_FCI(Confounded_Edges,test_matrix)
  pheatmap::pheatmap(FCI_adj,cluster_rows = F,cluster_cols = F)
  info[i,12]<-sum(abs(FCI_adj-LL_norm_adj)/2)
  
  
  
  info[i,1]=sum(LL_norm_adj)/2 # Source_edge
  info[i,2]=sum(L_norm_adj)/2 # EBMF_edge
  info[i,3]=sum(L_norm_adj[LL_norm_adj==1])/2 # "EBMF_TP"
  
  mat_diff<-L_norm_adj - LL_norm_adj
  info[i,4]=sum(mat_diff[mat_diff>0])/2 # "EBMF_FP"
  info[i,5]=abs(sum(mat_diff[mat_diff<0]))/2 # "EBMF_FN"
  
  
  info[i,8]= sum(FCI_adj[LL_norm_adj==1])/2 # FCI_TP
  
  mat_diff<-FCI_adj - LL_norm_adj
  info[i,9]= sum(mat_diff[mat_diff>0])/2 # FCI_FP
  info[i,10]= abs(sum(mat_diff[mat_diff<0]))/2 #FCI_FN
  
  Selected_Edges<-net_FCI$edges[net_FCI$edges[,3]!="<->",]
  info[i,13]<-nrow(Selected_Edges)
  FCI_adj<-gene_adj_mat_from_FCI(Selected_Edges,test_matrix)
  info[i,14]= sum(FCI_adj[LL_norm_adj==1])/2 # FCI_TP
  
}

write.csv(info,"FCI.csv")









