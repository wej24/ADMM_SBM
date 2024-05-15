#####################################################################
# Simulation in section 5.3
####################################################################

library(igraph)
library(Matrix)
library(irlba)
library(fossil)
library(mclust)

source("../../Algorithm/cvADMM.R")
source("../simulationFunc.R")


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


seed <- 1
numGraph <- 50
numNode <- 1000
sparseRho <- 0.15
warm <- 1
K <- 10
parallel <- 0
Bmat <- geneB(2,K)

set.seed(seed)
block.sizes <- c(rep(0.15,2),rep(0.1,3), rep(0.08,5)) * numNode
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

nfold <- 5
lambdas <- 10^seq(-4,4,by = 0.005)

fold <- foldFunc2(seed, nfold, length(Alist))
Winit <- matrix(0,K,K)
res <- gridLambdaSRL(A,A,Z,lambdas,rho = 2.5,Winit,convergence = 1e-10,maxiter=10000,warm)

if (parallel == 1) {
  library(doParallel)
  cl <- makeCluster(future::availableCores(),type="FORK")
  registerDoParallel(cl)
  cvRes <- foreach(i = 1:nfold) %dopar% {
    cvFunc2(fold, i, Alist, Z, lambdas, rho=2.5,Winit, warm,convergence = 1e-10, maxiter = 10000)
  }
  stopCluster(cl) 
} else {
  cvRes <- list()
  for (i in 1:nfold) {
    cvRes[[i]] <- cvFunc2(fold, i, Alist, Z, lambdas, rho=2.5,Winit, warm,convergence = 1e-10, maxiter = 10000)
  }
}
cv.error <- do.call(rbind, cvRes)
cvMean <- colMeans(cv.error)
indLambda <- which.min(cvMean)

estCVB <- res$estW[[indLambda]] / sparseRho
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
                   cvErr=cvErr,avgErr=avgErr,
                   avgErr_refit=avgErr_refit,lrErr_refit=lrErr_refit,
                   cvRank=cvRank)

