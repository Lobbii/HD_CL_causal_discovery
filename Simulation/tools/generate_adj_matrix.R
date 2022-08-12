gen_adj_mat<-function(Loading_matrix){
  adj_matrix<-matrix(0,nrow=nrow(Loading_matrix),ncol=nrow(Loading_matrix))
  for(i in 1:ncol(Loading_matrix)){
    x_list<-c()
    for(j in 1:nrow(Loading_matrix)){
      if(Loading_matrix[j,i]!=0){
        x_list<-c(x_list,j)
      }
    }
    # print(length(x_list))
    x_pair<-expand.grid(x_list, x_list)
    for (z in 1:nrow(x_pair)){adj_matrix[x_pair[z,1],x_pair[z,2]]<-1}
  }
  diag(adj_matrix) <- 0
  pheatmap::pheatmap(adj_matrix,cluster_rows = F,cluster_cols = F)
  return(adj_matrix)
}
  

gene_adj_mat_from_FCI<-function(Confounded_Edges,test_matrix){
  Confounded_Edges<-as.data.frame(Confounded_Edges)
  Confounded_Edges[,1]<-match(Confounded_Edges[,1],paste0("V",seq(200)) )
  Confounded_Edges[,2]<-match(Confounded_Edges[,2],paste0("V",seq(200)) )
  
  adj_matrix<-matrix(0,nrow=ncol(test_matrix),ncol=ncol(test_matrix))
  
  for(i in 1:nrow(Confounded_Edges)){
    adj_matrix[Confounded_Edges[i,1],Confounded_Edges[i,2]]<-1
    adj_matrix[Confounded_Edges[i,2],Confounded_Edges[i,1]]<-1
  }
  return(adj_matrix)
}
