edge_intersect <- function(x, y)
{
  edge.same <- c()
  # for each node, calculate the intersection of the node's
  # adjacencies in x and y
  for (node in names(x$adjacencies)) {
    # get the size for intersection of the adjacencies of 'node'
    # in est graph and true graph
    edge_same<-intersect(x$adjacencies[[node]], y$adjacencies[[node]])
    if (length(edge_same) != 0) {
      edge.same <- c(edge.same,paste(node,edge_same))
    }
  }
  return(edge.same)
}

require(data.table)

arrowhead_graph_recovery<-function(x,y){
  
  edge.same<-edge_intersect(x,y) 
  x.edge.same<-as.data.frame(x$edges,stringsAsFactors=F)[paste(x$edges[,1],x$edges[,2],sep = ' ') %in% edge.same,]
  setDT(x.edge.same)[, c('endpoint_1','edge','endpoint_2') := tstrsplit(V3, "")]
  x.edge.same<-rbind(x.edge.same[,c(1,2,6)],x.edge.same[,c(2,1,4)],use.names=FALSE)
  x.edge.same<-as.data.frame(x.edge.same,stringsAsFactors=F)
  
  y.edge.same<-as.data.frame(y$edges,stringsAsFactors=F)[paste(y$edges[,1],y$edges[,2],sep = ' ') %in% edge.same,]
  setDT(y.edge.same)[, c('endpoint_1','edge','endpoint_2') := tstrsplit(V3, "")]
  y.edge.same<-rbind(y.edge.same[,c(1,2,6)],y.edge.same[,c(2,1,4)],use.names=FALSE)
  y.edge.same<-as.data.frame(y.edge.same,stringsAsFactors=F)
  
  
  TP <- 0
  FN <- 0
  FP <- 0
  TN <- 0 
  for (i in 1:nrow(x.edge.same) ){
    index<- match(do.call(paste, x.edge.same[i,c(1,2)]) , do.call(paste, y.edge.same[,c(1,2)]))
    # True edge: A *-> B
    if (x.edge.same[index,3] == '>'){
      if (y.edge.same[i,3] == '>' | y.edge.same[i,3] == '<'){TP = TP + 1 }
      if (y.edge.same[i,3] == '-'){FN = FN + 1 }
    }  
    # Trueedge:A*--B
    if (x.edge.same[index,3] == '-'){
      if (y.edge.same[i,3] == '>' | y.edge.same[i,3] == '<'){FP = FP + 1 }
      if (y.edge.same[i,3] == '-'){TN = TN + 1 }
    }  
    # print(c('true',x.edge.same[index,3],'flase',y.edge.same[i,3],'TP',TP,'FN',FN,'FP',FP,TN,TN))
  }
  
  arrowhead_recall<- TP/(TP+FN)
  arrowhead_precision<- TP/(TP+FP)
  N = length(edge.same)
  S = (TP+FN)/N
  P = (TP+FP)/N
  arrowhead_MCC<- (TP/N-S*P)/sqrt(P*S*(1-S)*(1-P))
  outlist<-list('arrowhead_recall'=arrowhead_recall,'arrowhead_precision'=arrowhead_precision,'arrowhead_MCC'=arrowhead_MCC)
  return(outlist)
}



