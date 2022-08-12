#' Load reference sample causal network
#'
#' this functions loads reference sample causal network (FGES) from TXT file
#' @param file Path to TXT file containing the reference sample network
#' @return A list containing edge, nodes and adjacency information of causal graph learned by FGES
#' @export

read_causality_graph<-function(file){
  if (!file.exists(file)) {
    stop("Cannot find file!")
  }
  tmp <- readLines(file)
  if (tmp[1] != "Graph Nodes:") {    stop("file does not contain a compatible graph.")}

  if (grepl(";", tmp[2])){
    split <- ";"
  }  else {split <- ","}

  nodes <- unlist(strsplit(tmp[2], split = split))
  nodes <- sort(nodes)
  if (tmp[4] != "Graph Edges:") {stop("file does not contain a compatible graph.")}

  tmp <- tmp[-(1:4)]
  edges <- c()
  line <- tmp[1]
  while (!is.na(line) && line != "") {
    line <- sub("[0-9]+\\. ", "", line)
    edge <- unlist(strsplit(line, split = " "))
    if (is.na(sum(match(c(edge[1], edge[3]), nodes)))) {
      print(line)
      stop("file contains a malformed causality graph.")
    }
    edges <- c(edges, c(edge[1], edge[3], edge[2]))
    tmp <- tmp[-1]
    line <- tmp[1]
  }
  edges<-matrix(edges, ncol=3, byrow=TRUE)

  adjacencies <- lapply(nodes, function(node) {
    neighborhood <- c()
    # loop over the edges
    for (i in 1:nrow(edges)) {
      # if a node is in an edge, find its partner (neighbor)
      if( node %in% edges[i, 1:2]) {
        neighbor <- edges[i, c(node != edges[i, 1:2], F)]
        if (!(neighbor %in% neighborhood))
          neighborhood <- c(neighborhood, neighbor)
      }
    }
    return(neighborhood)
  })
  names(adjacencies) <- nodes

  return(list('edges'=edges,'nodes'=nodes,'adjacencies'=adjacencies))
}


Get_adjacencies <- function(nodes,edges) {
  adjacencies <- lapply(nodes, function(node) {
    neighborhood <- c()
    # loop over the edges
    for (i in 1:nrow(edges)) {
      # if a node is in an edge, find its partner (neighbor)
      if( node %in% edges[i, 1:2]) {
        neighbor <- edges[i, c(node != edges[i, 1:2], F)]
        if (!(neighbor %in% neighborhood))
          neighborhood <- c(neighborhood, neighbor)
      }
    }
    return(neighborhood)
  })
  names(adjacencies) <- nodes
  return(adjacencies)
}

Split_graph<- function(graph,con_var){
  
  edges<-graph$edges
  nodes<-graph$nodes
  edges_v1<-as.integer(edges[,1] %in% con_var)
  edges_v2<-as.integer(edges[,2] %in% con_var)
  Type<-edges_v1+edges_v2

  edges_con<-edges[Type==2,,drop=FALSE]
  edges_mix<-edges[Type==1,,drop=FALSE]
  edges_dis<-edges[Type==0,,drop=FALSE]
  
  if(nrow(edges_con)>0){
    adjacencies_con<-Get_adjacencies(nodes,edges_con)
  } else{
    adjacencies_con<-NA
  }
  if(nrow(edges_mix)>0){
    adjacencies_mix<-Get_adjacencies(nodes,edges_mix)
  } else{
    adjacencies_mix<-NA
  }
  if(nrow(edges_dis)>0){
    adjacencies_dis<-Get_adjacencies(nodes,edges_dis)
  } else{
    adjacencies_dis<-NA
  }
  
  net_con<-list('edges'=edges_con,'nodes'=nodes,'adjacencies'=adjacencies_con)
  net_mix<-list('edges'=edges_mix,'nodes'=nodes,'adjacencies'=adjacencies_mix)
  net_dis<-list('edges'=edges_dis,'nodes'=nodes,'adjacencies'=adjacencies_dis)
  
  return(list("Continuous"=net_con,"Mixed"=net_mix,"Discrete"=net_dis,
              "Edge Number" =list("Continuous" = nrow(edges_con),"Mixed" = nrow(edges_mix),"Discrete" = nrow(edges_dis)) ))
 
}

