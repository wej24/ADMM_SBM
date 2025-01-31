#######################################
# Help function for MonoSBM simulation 
#######################################

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


# dataGene <- function(seed, K, block.sizes, sparseRho, numGraph, numNode, Bmat) {
#   memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
#   Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
#   Z <- matrix(0,numNode,K)
#   Z[cbind(1:numNode,memb_index)] <- 1
#   return(list(Alist=Alist, Z=Z))
# }


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




monoSim <- function(seed, K, Bmat, sparseRho, numNode, numGraph, block.sizes, TrueZ=TRUE, 
                    lambdas, nfold=5, parallel=TRUE, ncores=5, warm=0) {
  set.seed(seed)
  memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
  Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
  Abar <- Reduce("+", Alist) / length(Alist)
  if (TrueZ) {
    Z <- matrix(0,numNode,K)
    Z[cbind(1:numNode,memb_index)] <- 1
  } else {
    ## do clustering to estimate Z
    svdA <- svdr(Abar, K+1)
    U <- svdA$u[,1:K] %*% diag(svdA$d[1:K])
    mG <- Mclust(U, G=K)
    ARI <- adj.rand.index(memb_index, mG$classification)
    memb2 <- align_two_membership_vectors(memb_index, mG$classification)
    Z <- matrix(0,numNode,K)
    Z[cbind(1:numNode,memb2)] <- 1
  }
  Winit <- matrix(0, K, K)
  cvRes <- cvSRL2(seed, nfold, Alist, Z, 
                  lambdas, rho = 2.5, Winit = Winit,
                  convergence = 1e-10, maxiter = 10000,
                  warm = warm, parallel = parallel, nCores = ncores)
  if (warm == 0) {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas[cvRes$indLambda],rho = 2.5,
                            Winit = Winit,
                            convergence = 1e-10, maxiter=10000,
                            warm = 0)
    estB <- resSRL$estW[[1]] / sparseRho
  } else {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas,rho = 2.5,
                            Winit = Winit,
                            convergence = 1e-10, maxiter=10000,
                            warm = 1)
    estB <- resSRL$estW[[cvRes$indLambda]] / sparseRho
  }
  ourRank <- qr(estB)$rank
  ourMSE <- sqrt(sum((estB - Bmat)^2))
  avgB <- (t(Z)/colSums(Z)) %*% Abar %*% t(t(Z)/colSums(Z))
  avgB <- avgB/sparseRho
  avgRank <- qr(avgB)$rank
  avgMSE <- sqrt(sum((avgB-Bmat)^2))
  
  res <- list(estB = estB, ourRank = ourRank, ourMSE = ourMSE,
              avgB = avgB, avgRank = avgRank, avgMSE = avgMSE)
  if (!TrueZ) {
    res$ARI <- ARI
    res$Z <- Z
    res$Alist <- Alist
    res$Abar <- Abar
  }
  return(res)
}



specEmbed <- function(d, K, Abar, Zhat, B, sparseRho=1) {
  svdA <- svdr(Abar, d + 1)
  if (d == 1) {
    Ad <- svdA$d[1] * svdA$u[,1] %*% t(svdA$v[,1])
  } else {
    Ad <- svdA$u[,1:d] %*% diag(svdA$d[1:d]) %*% t(svdA$v[,1:d])
  }
  if (K == 1) {
    Ak <- svdA$d[1] * svdA$u[,1] %*% t(svdA$v[,1])
  } else {
    Ak <- svdA$u[,1:K] %*% diag(svdA$d[1:K]) %*% t(svdA$v[,1:K])
  }
  Bk <- (t(Zhat)/colSums(Zhat)) %*% Ak %*% t(t(Zhat)/colSums(Zhat))
  Bk <- Bk / sparseRho
  Bd <- (t(Zhat)/colSums(Zhat)) %*% Ad %*% t(t(Zhat)/colSums(Zhat))
  Bd <- Bd / sparseRho
  specMSEd <- sqrt(sum((Bd-B)^2))
  specMSEk <- sqrt(sum((Bk-B)^2))
  return(list(specMSEd=specMSEd, specMSEk=specMSEk))
}

lowRank <- function(d, Abar, Zhat, sparseRho=1) {
  Bavg <- (t(Zhat)/colSums(Zhat)) %*% Abar %*% t(t(Zhat)/colSums(Zhat))
  Bavg <- Bavg/sparseRho
  
  svdB <- svdr(Bavg, min(d+1,ncol(Zhat)))
  if (d == 1) {
    Blow <- svdB$d[1] * svdB$u[,1] %*% t(svdB$v[,1])
  } else {
    Blow <- svdB$u[,1:d] %*% diag(svdB$d[1:d]) %*% t(svdB$v[,1:d])
  }
  return(Blow)
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


estBLow <- function(seed,Alist,Zhat,K,B,sparseRho,nfold) {
  set.seed(seed)
  A <- Reduce("+", Alist) / length(Alist)
  Bhat <- (t(Zhat)/colSums(Zhat)) %*% A %*% t(t(Zhat)/colSums(Zhat))
  Bhat <- Bhat/sparseRho
  MSEavg <- sqrt(sum((Bhat-B)^2))
  svdB <- svd(Bhat)
  fold <- foldFunc2(seed, nfold, length(Alist))
  d <- estBcv(fold, nfold, Alist, Zhat, K,sparseRho, B)
  if (d == 1) {
    Blow <- svdB$d[1] * svdB$u[,1] %*% t(svdB$v[,1])
  } else {
    Blow <- svdB$u[,1:d] %*% diag(svdB$d[1:d]) %*% t(svdB$v[,1:d])
  }
  MSElow <- sqrt(sum((Blow-B)^2))
  return(list(avgMSE=MSEavg,avgLRMSE=MSElow,avgLRRank=d))
}




geneB <- function(seed, K) {
  set.seed(seed)
  u <- runif(K,0.2,0.9)
  u1 <- runif(K,0.2,0.9)
  B <- u %*% t(u) + u1 %*% t(u1)
  if (max(B) > 1) {
    B <- B / 2
  }
  return(B)
}

reEstZSim <- function(seed, K, Bmat, sparseRho, numNode, numGraph, block.sizes,
                      lambdas, nfold=5, parallel=TRUE, ncores=5, warm=1) {
  set.seed(seed)
  memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
  Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
  A <- Reduce("+", Alist) / length(Alist)
  svdA <- svdr(A, K+1)
  U <- svdA$u[,1:K] %*% diag(svdA$d[1:K])
  mK <- kmeans(U, centers=K,nstart=100)
  ARI <- adj.rand.index(memb_index, mK$cluster)
  memb2 <- align_two_membership_vectors(memb_index, mK$cluster)
  mic <- sum(memb_index != memb2) / length(memb_index)
  Z <- matrix(0,numNode,K)
  Z[cbind(1:numNode,memb2)] <- 1
  Bhat <- (t(Z)/colSums(Z)) %*% A %*% t(t(Z)/colSums(Z))
  Bhat <- Bhat / sparseRho
  avgErr <- sqrt(sum((Bhat - Bmat)^2))
  fold <- foldFunc2(seed, nfold, length(Alist))
  Winit <- matrix(0,K,K)
  res <- gridLambdaSRL(A,A,Z,lambdas,rho = 2.5,
                       Winit = Winit,
                       convergence = 1e-10, maxiter=10000,
                       warm = warm)
  cvRes <- cvSRL2(seed, nfold, Alist, Z, 
                  lambdas, rho = 2.5, Winit = Winit,
                  convergence = 1e-10, maxiter = 10000,
                  warm = warm, parallel = parallel, nCores = ncores)
  estCVB <- res$estW[[cvRes$indLambda]] / sparseRho
  cvRank <- qr(estCVB)$rank
  cvErr <- sqrt(sum((estCVB-Bmat)^2))
  avgRank <- qr(Bhat)$rank
  if (cvRank == 1) {
    U <- svdA$d[1] * svdA$u[,1]
  } else {
    U <- svdA$u[,1:cvRank] %*% diag(svdA$d[1:cvRank])
  }
  mK <- kmeans(U, centers=K,nstart=100)
  ARI_refit <- adj.rand.index(memb_index, mK$cluster)
  memb2 <- align_two_membership_vectors(memb_index, mK$cluster)
  mic_refit <- sum(memb_index != memb2) / length(memb_index)
  Z <- matrix(0,numNode,K)
  Z[cbind(1:numNode,memb2)] <- 1
  Bhat_refit <- (t(Z)/colSums(Z)) %*% A %*% t(t(Z)/colSums(Z))
  Bhat_refit <- Bhat_refit / sparseRho
  avgErr_refit <- sqrt(sum((Bhat_refit - Bmat)^2))
  
  svdB <- svd(Bhat_refit)
  if (cvRank == 1) {
    B2 <- svdB$d[1] * svdB$u[,1] %*% t(svdB$v[,1])
  } else {
    B2 <- svdB$u[,1:cvRank] %*% diag(svdB$d[1:cvRank]) %*% t(svdB$v[,1:cvRank])
  }
  lrErr_refit <- sqrt(sum((B2-Bmat)^2))
  resSummary <- list(mic=mic,mic_refit=mic_refit,
                     err=cvErr, err_refit=lrErr_refit,
                     ourRank=cvRank)
  return(resSummary)
}