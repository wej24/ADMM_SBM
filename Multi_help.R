##############################################################
# Help function for MultiSBM simulation
##############################################################


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

## A square in bias-ajusted clustering 
A_square <- function(Alist) {
  L <- length(Alist)
  sum_mat <- numeric(0)
  for(i in 1:L){
    tmp <- crossprod(Alist[[i]])
    if(i == 1) {
      sum_mat <- tmp
    } else {
      sum_mat <- sum_mat + tmp
    }
  }
  diag(sum_mat) <- 0
  sum_mat <- as.matrix(sum_mat)
  return(sum_mat)
}


## estimate B on each Layer
simB <- function(A, Zhat) {
  Bhat <- (t(Zhat)/colSums(Zhat)) %*% A %*% t(t(Zhat)/colSums(Zhat))
  return(Bhat[upper.tri(Bhat,diag = TRUE)])
}

## estimate B on each Graph community
## indB is the indices of all graphs in $\tilde{\ell}$-th group
estB <- function(seed, Alist, Z, indB, nfold, lambdas, Winit, sparseRho, B,
                 warm = 0, parallel = TRUE, nCores = 5, 
                 convergence = 1e-10, maxiter = 10000) {
  Abar <- Reduce("+", Alist[indB]) / length(Alist[indB])
  cvRes <- cvSRL2(seed, nfold, Alist[indB], Z, 
                 lambdas, rho = 2.5, Winit = Winit,
                 convergence = convergence, maxiter = maxiter,
                 warm = warm, parallel = parallel, nCores = nCores)
  if (warm == 0) {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas[cvRes$indLambda],rho = 2.5,
                            Winit = Winit,
                            convergence = convergence, maxiter=maxiter,
                            warm = 0)
    estB <- resSRL$estW[[1]] / sparseRho
  } else {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas,rho = 2.5,
                            Winit = Winit,
                            convergence = convergence, maxiter=maxiter,
                            warm = 1)
    estB <- resSRL$estW[[cvRes$indLambda]] / sparseRho
  }
  ourRank <- qr(estB)$rank
  ourMSE <- sqrt(sum((estB - B)^2))
  avgB <- (t(Z)/colSums(Z)) %*% Abar %*% t(t(Z)/colSums(Z))
  avgB <- avgB/sparseRho
  avgRank <- qr(avgB)$rank
  avgMSE <- sqrt(sum((avgB-B)^2))
  
  res <- list(estB = estB, ourRank = ourRank, ourMSE = ourMSE,
              avgB = avgB, avgRank = avgRank, avgMSE = avgMSE)
}



multiSim <- function(seed=seed, K=K, B1=B1, B2=B2, B3=B3, B4=B4, sparseRho=sparseRho,
                     numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
                     lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=0) {
  set.seed(seed)
  memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
  Alist1 <- Afunc(seed, numGraph, numNode, sparseRho, B1, block.sizes)
  Alist2 <- Afunc(seed, numGraph, numNode, sparseRho, B2, block.sizes)
  Alist3 <- Afunc(seed, numGraph, numNode, sparseRho, B3, block.sizes)
  Alist4 <- Afunc(seed, numGraph, numNode, sparseRho, B4, block.sizes)
  Alist <- c(Alist1, Alist2,Alist3,Alist4)
  A_2 <- A_square(Alist)
  svdA_2 <- svdr(A_2, K+1)
  
  ## do clustering to estimate Z
  U2 <- svdA_2$u[,1:K] %*% diag(svdA_2$d[1:K])
  mG <- Mclust(U2, G=K)
  ARI <- adj.rand.index(memb_index, mG$classification)
  memb2 <- align_two_membership_vectors(memb_index, mG$classification)
  Z <- matrix(0,numNode,K)
  Z[cbind(1:numNode,memb2)] <- 1
  
  ## estimate B for each layer and then cluster layer
  uppB <- lapply(Alist,simB, Zhat=Z)
  uppB <- do.call(rbind,uppB)
  mB <- Mclust(uppB, G=4)
  memb_B <- c(rep(1,50),rep(2,50), rep(3,50),rep(4,50))
  ARIB <- adj.rand.index(memb_B, mB$classification)
  memb_B2 <- align_two_membership_vectors(memb_B, mB$classification)
  
  ## indices of graphs in each graph community 
  indB1 <- which(memb_B2 == 1)
  indB2 <- which(memb_B2 == 2)
  indB3 <- which(memb_B2 == 3)
  indB4 <- which(memb_B2 == 4)
  
  ## estimate B for each graph community 
  Winit <- matrix(0,K,K)
  res1 <- estB(seed, Alist, Z, indB1, nfold, lambdas, Winit, sparseRho, B1,
               warm = warm, parallel = parallel, nCores = ncores, 
               convergence = 1e-10, maxiter = 10000)
  res2 <- estB(seed, Alist, Z, indB2, nfold, lambdas, Winit, sparseRho, B2,
               warm = warm, parallel = parallel, nCores = ncores, 
               convergence = 1e-10, maxiter = 10000)
  res3 <- estB(seed, Alist, Z, indB3, nfold, lambdas, Winit, sparseRho, B3,
               warm = warm, parallel = parallel, nCores = ncores, 
               convergence = 1e-10, maxiter = 10000)
  res4 <- estB(seed, Alist, Z, indB4, nfold, lambdas, Winit, sparseRho, B4,
               warm = warm, parallel = parallel, nCores = ncores, 
               convergence = 1e-10, maxiter = 10000)
  return(list(ARI=ARI,ARIB=ARIB,Z=Z,
              indB1=indB1,indB2=indB2,indB3=indB3,indB4=indB4,
              Alist=Alist,
              ourRank=c(res1$ourRank,res2$ourRank,res3$ourRank,res4$ourRank),
              ourMSE=c(res1$ourMSE,res2$ourMSE,res3$ourMSE,res4$ourMSE),
              avgRank=c(res1$avgRank,res2$avgRank,res3$avgRank,res4$avgRank),
              avgMSE=c(res1$avgMSE,res2$avgMSE,res3$avgMSE,res4$avgMSE)))
}



estBcv <- function(fold, nfold, Alist, Zhat, K,sparseRho, B) {
  MSEA <- vector("list", nfold)
  for (i in 1:nfold) {
    indTrain <- which(fold != i)
    indValid <- which(fold == i)
    Avalid <- Reduce("+", Alist[indValid]) / length(indValid)
    Atrain <- Reduce("+", Alist[indTrain]) / length(indTrain)
    estBhat <- (t(Zhat)/colSums(Zhat)) %*% Atrain %*% t(t(Zhat)/colSums(Zhat))
    estBhat <- estBhat / sparseRho
    svdBhat <- svd(estBhat)
    diag(Avalid) <- 0
    for (j in 1:K) {
      if (j == 1) {
        estB <- svdBhat$d[1] * svdBhat$u[,1] %*% t(svdBhat$v[,1])
      } else {
        estB <- svdBhat$u[,1:j] %*% diag(svdBhat$d[1:j]) %*% t(svdBhat$v[,1:j])
      }
      estA <- Zhat %*% (estB*sparseRho) %*% t(Zhat)
      diag(estA) <- 0
      MSEA[[i]] <- c(MSEA[[i]], sum((estA-Avalid)^2)) 
    }
  }
  cv.error <- do.call(rbind, MSEA)
  cvMean <- colMeans(cv.error)
  d <- which.min(cvMean)
  return(d)
}

estBLow <- function(seed,Alist,Zhat,indB,K,B,sparseRho,nfold) {
  set.seed(seed)
  A <- Reduce("+", Alist[indB]) / length(indB)
  Bhat <- (t(Zhat)/colSums(Zhat)) %*% A %*% t(t(Zhat)/colSums(Zhat))
  Bhat <- Bhat/sparseRho
  MSEavg <- sqrt(sum((Bhat-B)^2))
  svdB <- svd(Bhat)
  fold <- foldFunc2(seed, nfold, length(indB))
  d <- estBcv(fold, nfold, Alist[indB], Zhat, K,sparseRho, B)
  if (d == 1) {
    Blow <- svdB$d[1] * svdB$u[,1] %*% t(svdB$v[,1])
  } else {
    Blow <- svdB$u[,1:d] %*% diag(svdB$d[1:d]) %*% t(svdB$v[,1:d])
  }
  MSElow <- sqrt(sum((Blow-B)^2))
  return(list(MSEavg=MSEavg,MSElow=MSElow,
              rd=d,Blow=Blow,Bhat=Bhat))
}

estBheLow <- function(seed, Alist, Zhat, 
                      indB1, indB2, indB3, indB4,
                      B1, B2, B3, B4,
                      sparseRho, nfold) {
  set.seed(seed)
  res1 <- estBLow(seed,Alist,Zhat,indB1,K,B1,sparseRho,nfold)
  res2 <- estBLow(seed,Alist,Zhat,indB2,K,B2,sparseRho,nfold)
  res3 <- estBLow(seed,Alist,Zhat,indB3,K,B3,sparseRho,nfold)
  res4 <- estBLow(seed,Alist,Zhat,indB4,K,B4,sparseRho,nfold)
  return(list(avgMSE=c(res1$MSEavg,res2$MSEavg,res3$MSEavg,res4$MSEavg),
              avgLRMSE=c(res1$MSElow,res2$MSElow,res3$MSElow,res4$MSElow),
              avgLRRank=c(res1$rd,res2$rd,res3$rd,res4$rd)))
}

estBTrun_each <- function(seed,Alist,Zhat,indB,K, d, B,sparseRho) {
  set.seed(seed)
  A <- Reduce("+", Alist[indB]) / length(indB)
  svdA <- svdr(A, K+1)
  Ak <- svdA$u[,1:K] %*% diag(svdA$d[1:K]) %*% t(svdA$v[,1:K])
  if (K != d) {
    if (d == 1) {
      Ad <- svdA$d[1] * svdA$u[,1] %*% t(svdA$v[,1])
    } else {
      Ad <- svdA$u[,1:d] %*% diag(svdA$d[1:d]) %*% t(svdA$v[,1:d])
    }
  } else {
    Ad <- Ak
  }
  Bk <- (t(Zhat)/colSums(Zhat)) %*% Ak %*% t(t(Zhat)/colSums(Zhat))
  Bk <- Bk/sparseRho
  MSEk <- sqrt(sum((Bk-B)^2))
  Bd <- (t(Zhat)/colSums(Zhat)) %*% Ad %*% t(t(Zhat)/colSums(Zhat))
  Bd <- Bd/sparseRho
  MSEd <- sqrt(sum((Bd-B)^2))
  Bhat <- (t(Zhat)/colSums(Zhat)) %*% A %*% t(t(Zhat)/colSums(Zhat))
  Bhat <- Bhat/sparseRho
  MSEavg <- sqrt(sum((Bhat-B)^2))
  return(list(specMSEd=MSEd, specMSEk=MSEk, avgMSE=MSEavg))
}

estBTrun <- function(seed, K, Alist,Zhat,
                     indB1, indB2, indB3, indB4,
                     d1, d2, d3, d4,
                     B1, B2, B3, B4,
                     sparseRho) {
  res1 <- estBTrun_each(seed,Alist,Zhat,indB1, K, d1, B1, sparseRho)
  res2 <- estBTrun_each(seed,Alist,Zhat,indB2, K, d2, B2, sparseRho)
  res3 <- estBTrun_each(seed,Alist,Zhat,indB3, K, d3, B3, sparseRho)
  res4 <- estBTrun_each(seed,Alist,Zhat,indB4, K, d4, B4, sparseRho)
  return(list(avgMSE=c(res1$avgMSE,res2$avgMSE,res3$avgMSE,res4$avgMSE),
              specMSEd=c(res1$specMSEd,res2$specMSEd,res3$specMSEd,res4$specMSEd),
              specMSEk=c(res1$specMSEk,res2$specMSEk,res3$specMSEk,res4$specMSEk)))
}