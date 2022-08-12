Mapping_LF<-function(L,LL,Path){
  Mer_mat<-cbind(LL,L)
  cor_mat<-cor(Mer_mat)
  p<-pheatmap::pheatmap(cor_mat,cluster_rows = F,cluster_cols = F)
  cowplot::save_plot(paste0(Path,"/cor.png"), p,
                     base_aspect_ratio = 1, base_height =15)
  
  
  require(clue)
  A<-abs(cor_mat[1:ncol(LL),(ncol(LL)+1):ncol(Mer_mat)])
  B <- diag(1, nrow(A),ncol(A))
  A<-t(A);B<-t(B)
  
  n <- nrow(A) 
  D <- matrix(NA, nrow(A),nrow(A)) 
  for (i in 1:n) { 
    for (j in 1:n) { 
      D[j, i] <- (sum((B[j, ] - A[i, ])^2)) 
    }
  } 
  vec <- c(solve_LSAP(D)) 
  X<-list(A=A[vec,], pvec=vec) 
  mat<-t(X$A) 
  
  p<-pheatmap::pheatmap(mat,cluster_rows = F,cluster_cols = F)
  cowplot::save_plot(paste0(Path,"/cor_diagonal.png"), p,
                     base_aspect_ratio = 1.5, base_height =15)
  
  L_var<-colnames(mat)[1:ncol(LL)]
  L_var_index<-match(L_var,colnames(L))
  colnames(L)[L_var_index]<-rownames(mat)
  
  return(list('mapped_matrix'=L[,L_var_index],"Latent_Matrix" = L))
}