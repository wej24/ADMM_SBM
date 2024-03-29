###################################
# Simulation Function for HoSBM 
# with deterministic block size 
###################################

sbmA <- function(index, numNode, rho, Bmat, block.sizes) {
  graph <- sample_sbm(numNode,
                      rho * Bmat,
                      block.sizes,
                      directed = FALSE,
                      loops = TRUE)
  Amat <- as_adjacency_matrix(graph)
  return(Amat)
}

## Afunc to generate L=numGraph graphs in HoSBM
Afunc <- function(seed,numGraph, numNode, rho, Bmat, block.sizes) {
  set.seed(seed)
  Alist <- lapply(1:numGraph, sbmA,
                  numNode = numNode, rho = rho, 
                  Bmat = Bmat, block.sizes = block.sizes)
  return(Alist)
}


# dataGene <- function(seed, K, block.sizes, sparseRho, numGraph, numNode, Bmat) {
#   memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
#   Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
#   Z <- matrix(0,numNode,K)
#   Z[cbind(1:numNode,memb_index)] <- 1
#   return(list(Alist=Alist, Z=Z))
# }


# align the clustering membership based on Hamming distance
# The following code from Bias-adjusted paper Github
# requires both vec1 and vec2 to have the same number of unique clusters
align_two_membership_vectors <- function(vec1, vec2, override = F){
  if(override & length(unique(vec1)) > length(unique(vec2))){
    len <- length(unique(vec1))
    vec2[1:len] <- 1:len
  }
  stopifnot(length(unique(vec1)) == length(unique(vec2)))
  
  K <- length(unique(vec1))
  tab <- table(vec1, vec2)
  
  permn_list <- .permn(K)
  similarity_vec <- sapply(permn_list, function(x){
    tab2 <- tab
    tab2 <- tab2[,x]
    sum(diag(tab2))/sum(tab2)
  })
  
  idx <- which.max(similarity_vec)
  
  vec2_new <- vec2
  for(i in 1:K){
    vec2_new[which(vec2 == permn_list[[idx]][i])] <- i
  }
  
  return(vec2_new)
}


## https://github.com/cran/combinat/blob/master/R/permn.R
.permn <- function(x, fun = NULL, ...) {
  if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) x <- seq(x)
  n <- length(x)
  nofun <- is.null(fun)
  out <- vector("list", gamma(n + 1))
  p <- ip <- seqn <- 1:n
  d <- rep(-1, n)
  d[1] <- 0
  m <- n + 1
  p <- c(m, p, m)
  i <- 1
  use <-  - c(1, n + 2)
  
  while(m != 1) {
    out[[i]] <- if(nofun) x[p[use]] else fun(x[p[use]], ...)
    i <- i + 1
    m <- n
    chk <- (p[ip + d + 1] > seqn)
    m <- max(seqn[!chk])
    if(m < n)
      d[(m + 1):n] <-  - d[(m + 1):n]
    index1 <- ip[m] + 1
    index2 <- p[index1] <- p[index1 + d[m]]
    p[index1 + d[m]] <- m
    tmp <- ip[index2]
    ip[index2] <- ip[m]
    ip[m] <- tmp
  }
  
  out
}
