adjacency_precision <- function(x, y)
{
  n.y.adjs <- sum(lengths(y$adjacencies))
  if (n.y.adjs == 0) {
    warning("y has no adjacencies. Returning NA")
    return(NA)
  }
  # calcluate the intersection of adjacencies over y.adjs (estimated adjacency) and return the ratio
  return(adjacency_intersect(x, y) / n.y.adjs)
}


adjacency_recall <- function(x, y)
{
  n.x.adjs <- sum(lengths(x$adjacencies))
  if (n.x.adjs == 0) {
    warning("x has no adjacencies. Returning NA")
    return(NA)
  }
  # calcluate the intersection of adjacencies over x.adjs (true adjacency) and return the ratio
  return(adjacency_intersect(x, y) / n.x.adjs)
}

# internal function that is used to cacluate the intersection of the adjacencies
# for each node in x and y
adjacency_intersect <- function(x, y)
{
  n.same <- 0
  # for each node, calculate the intersection of the node's
  # adjacencies in x and y
  for (node in names(x$adjacencies)) {
    # get the size for intersection of the adjacencies of 'node'
    # in est graph and true graph
    n.same <- n.same + length(
      intersect(x$adjacencies[[node]], y$adjacencies[[node]]))
  }
  return(n.same)
}


adjacency_MCC <- function(x, y)
{
  n.y.adjs <- sum(lengths(y$adjacencies))
  if (n.y.adjs == 0) {
    warning("y has no adjacencies. Returning NA")
    return(NA)
  }
  # calcluate the intersection of adjacencies over y.adjs (estimated adjacency) and return the ratio
  TP<-adjacency_intersect(x, y)/2
  N <- length(x$nodes)*(length(x$nodes)-1)/2
  S <- nrow(x$edges)/N  # true graph TP+FN /N
  P <- nrow(y$edges)/N # estimated graph TP+FP /N
  MCC <- (TP/N-S*P)/sqrt(S*P*(1-S)*(1-P)) 
  return(MCC)
}
