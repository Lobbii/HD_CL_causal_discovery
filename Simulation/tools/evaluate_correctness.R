evaluate_correctness <- function(net_estimated,net_true) {

  # net_estimated<-read_causality_graph(estimated_graph_path)
  # net_true<-read_causality_graph(true_graph_path)
  adj_precision<-adjacency_precision(net_true,net_estimated) # /n.y.adjs -> prdicted -> y 
  adj_recall<-adjacency_recall(net_true,net_estimated) # n.x.adjs -> true ->x
  adj_MCC<-adjacency_MCC(net_true,net_estimated) # n.x.adjs -> true ->x
  arh_recovery<-arrowhead_graph_recovery(net_true,net_estimated) # arrow_intersect(x, y)/n.y.arrows
 
 # out_list =list('adjacency_precision'=adj_precision,
 #                'adjacency_recall'= adj_recall,
 #                # 'adjacency_MCC' = adj_MCC,
 #                'arrowhead_precision'=arh_recovery$arrowhead_precision,
 #                'arrowhead_recall'=arh_recovery$arrowhead_recall)
 #                # 'arrowhead_MCC'= arh_recovery$arrowhead_MCC
 
 eva<-matrix(NA,nrow = 1,ncol = 4)
 eva[1,]<-c(adj_precision,adj_MCC,arh_recovery$arrowhead_precision,arh_recovery$arrowhead_recall)
 colnames(eva)<-c('AP',"AR","AHP","AHR")
  return(eva)
}



eva_adj<-function(graph_esti,graph_true){
  eva<-matrix(NA,nrow = 1,ncol = 12)
  colnames(eva)<-apply(expand.grid(c('AP',"AR","AHP","AHR"),c("_CC","_CD","_DD")), 1, paste, collapse="")
  if (graph_esti$`Edge Number`$Continuous > 0) {
    adj_precision<-adjacency_precision(graph_true$Continuous,graph_esti$Continuous) # /n.y.adjs -> prdicted -> y 
    adj_recall<-adjacency_recall(graph_true$Continuous,graph_esti$Continuous) # n.x.adjs -> true ->x
    if(adj_precision > 0){
      arh_recovery<-arrowhead_graph_recovery(graph_true$Continuous,graph_esti$Continuous)
      eva[1,][1:4]<-c(adj_precision,adj_recall,arh_recovery$arrowhead_precision,arh_recovery$arrowhead_recall)
    }else{
      eva[1,][1:4]<-c(adj_precision,adj_recall,0,0)
    }
    # print("Continous - Continous egde")
    # print(paste("precision: ",adj_precision ))
    # print(paste("recall: ",adj_recall ))
  } else( 
    eva[1,][1:4]<-c(0,0,0,0)
    # print("No Continous - Continous egde")
  )
  
  if(graph_esti$`Edge Number`$Mixed>0){
    adj_precision<-adjacency_precision(graph_true$Mixed,graph_esti$Mixed) # /n.y.adjs -> prdicted -> y 
    adj_recall<-adjacency_recall(graph_true$Mixed,graph_esti$Mixed) # n.x.adjs -> true ->x
    if(adj_precision > 0){
      arh_recovery<-arrowhead_graph_recovery(graph_true$Mixed,graph_esti$Mixed)
      eva[1,][5:8]<-c(adj_precision,adj_recall,arh_recovery$arrowhead_precision,arh_recovery$arrowhead_recall)
    }else{
      eva[1,][5:8]<-c(adj_precision,adj_recall,0,0)
    }
    # print("Continous - Discrete egde")
    # print(paste("precision: ",adj_precision ))
    # print(paste("recall: ",adj_recall ))
  }else(
    eva[1,][5:8]<-c(0,0,0,0)
    # print("No Continous - Discrete egde")
    )
  
  if(graph_esti$`Edge Number`$Discrete>0){
    adj_precision<-adjacency_precision(graph_true$Discrete,graph_esti$Discrete) # /n.y.adjs -> prdicted -> y 
    adj_recall<-adjacency_recall(graph_true$Discrete,graph_esti$Discrete) # n.x.adjs -> true ->x
    if(adj_precision > 0){
      arh_recovery<-arrowhead_graph_recovery(graph_true$Discrete,graph_esti$Discrete)
      eva[1,][9:12]<-c(adj_precision,adj_recall,arh_recovery$arrowhead_precision,arh_recovery$arrowhead_recall)
    }else{
      eva[1,][9:12]<-c(adj_precision,adj_recall,0,0)
    }
    # print("Discrete - Discrete egde")
    # print(paste("precision: ",adj_precision ))
    # print(paste("recall: ",adj_recall ))
  }else(
    eva[1,][9:12]<-c(0,0,0,0)
    # print("No Discrete - Discrete egde")
    )
  return(eva)
}