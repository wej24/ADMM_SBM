#####################################################################
# Additional simulations in the supplement
# Setting 1: K=2, d=2
# Since we run each replicate on the cluster,
# the following code runs one replicate 
# with sparseRho=log(n) / n, numNode=1000, numGraph=100
# To run other cases, 
# change sparseRho, numNode, numGraph, and seed correspondingly
####################################################################

library(igraph)
library(doParallel)
library(Matrix)
source("../../Algorithm/cvADMM.R")
source("../simulationFunc.R")


## Data generate
seed <- 1 ## seed for one single replication
K <- 2
numGraph <- 100 ## number of graphs L
numNode <- 1000 ## number of nodes n
block.sizes <- c(1/4, 3/4) * numNode
sparseRho <- log(numNode) / numNode ## sparse factor in SBM 
Bmat <- matrix(0.3, K, K)
diag(Bmat) <- c(0.5, 0.4)
memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
Z <- matrix(0,numNode,K)
Z[cbind(1:numNode,memb_index)] <- 1
Abar <- Reduce("+", Alist) / length(Alist)
## run CV
nfold <- 5
lambdas <- 10^seq(-3, 0, length.out=150)
Winit <- matrix(0, K, K)
warm <- 1
cvRes <- cvSRL(seed, nfold, Alist, Z, 
               lambdas, rho = 2.5, Winit = Winit,
               convergence = 1e-10, maxiter = 10000,
               warm = warm, parallel = TRUE, nCores = 10)

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





