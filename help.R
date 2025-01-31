##########################################
# Help function to plot Figure 1
##########################################


sbmA <- function(index, numNode, rho, Bmat, block.sizes) {
  graph <- sample_sbm(numNode,
                      rho * Bmat,
                      block.sizes,
                      directed = FALSE,
                      loops = TRUE)
  Amat <- as_adjacency_matrix(graph)
  return(Amat)
}

Afunc <- function(seed,numGraph, numNode, rho, Bmat, block.sizes) {
  set.seed(seed)
  Alist <- lapply(1:numGraph, sbmA,
                  numNode = numNode, rho = rho, 
                  Bmat = Bmat, block.sizes = block.sizes)
  return(Alist)
}


AtruncFunc <- function(k,svdA) {
  if (k == 1) {
    n <- nrow(svdA$u)
    d <- matrix(svdA$d[1], 1, 1)
    u <- matrix(svdA$u[,1], n, 1)
    v <- matrix(svdA$v[,1], n, 1)
  } else {
    d <- diag(svdA$d[1:k])
    u <- svdA$u[,1:k]
    v <- svdA$v[,1:k]
  }
  return(u %*% d %*% t(v))
}
