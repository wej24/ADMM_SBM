## Additional simulation in the supplement
## Single A with true Z

library(igraph)
library(doParallel)
library(Matrix)

source("singleA_func.R")

seed <- 1 ## for one replicate
numNode <- 1000
numGraph <- 1
K <- 2

## d = 1
u1 <- c(0.7, 0.7^2)
Bmat <- u1 %*% t(u1)

## d = 2
# Bmat <- matrix(0.3, K, K)
# diag(Bmat) <- c(0.5, 0.4)

sparseRho <- log(numNode) / numNode
block.sizes <- c(1/4, 3/4) * numNode
memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
Z <- matrix(0,numNode,K)
Z[cbind(1:numNode,memb_index)] <- 1

nfold <- 5
fold <- foldFunc(seed,nfold,K,block.sizes)
Winit <- matrix(0,K,K)
lambdas <- 10^seq(-3, 0, length.out=150)


cl <- makeCluster(future::availableCores(),type="FORK")
registerDoParallel(cl)
cvRes <- foreach(i = 1:nfold) %dopar% {
  cvFunc(fold, i, Alist[[1]], Z,lambdas, Winit)
}
stopCluster(cl)

cv.error <- do.call(rbind, cvRes)
cvMean <- colMeans(cv.error)

indLambda <- which.min(cvMean)

resSRL <- gridLambdaSRL(Alist[[1]],Alist[[1]],Z,lambdas, Winit, maxiter=10000)

estB <- resSRL$estW[[indLambda]] / sparseRho
ourRank <- qr(estB)$rank
ourMSE <- sqrt(sum((estB-Bmat)^2))

Bhat <- (t(Z)/colSums(Z)) %*% Alist[[1]] %*% t(t(Z)/colSums(Z))
Bhat <- Bhat/sparseRho
avgRank <- qr(Bhat)$rank
avgMSE <- sqrt(sum((Bhat-Bmat)^2))

resSummary <- list(ourB=estB, avgB=Bhat,
            ourRank=ourRank, avgRank=avgRank,
            ourMSE=ourMSE, avgMSE=avgMSE)


