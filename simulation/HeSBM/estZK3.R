########################################################################
# Simulation in section 5.4
# Setting: K=3
# Since we run each replicate on the cluster,
# the following code runs one replicate 
# with sparseRho=log(n) / n, numNode=1000, numGraph=50 for each group
# To run other cases, 
# change sparseRho and seed correspondingly
#######################################################################

library(igraph)
library(Matrix)
library(irlba)
library(fossil)
library(mclust)
library(gtools)
source("../../Algorithm/cvADMM.R")
source("../simulationFunc.R")
source("../simulationFuncHeSBM.R")

## Data generate
seed <- 101 ## seed for one single replication
K <- 3
numGraph <- 50 ## number of graph for each graph community
numNode <- 1000 ## number of nodes n 
block.sizes <- c(1/4, 1/4, 1/2) * numNode
memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
sparseRho <- log(numNode) / numNode ## sparse factor in SBM

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))

B1 <- eigen_mat %*% diag(c(1.2, 0.6, -0.7)) %*% t(eigen_mat)
u1 <- 0.8^(1:K)
B2 <- u1 %*% t(u1)
B3 <- eigen_mat %*% diag(c(1.7, 0, -0.6)) %*% t(eigen_mat)
B4 <- eigen_mat %*% diag(c(1.2, 0.6, 0.7)) %*% t(eigen_mat)

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
lambdas <- 10^seq(-4,-1.5,length.out=150)
Winit <- matrix(0,K,K)
nfold <- 5

res1 <- estB(seed, Alist, Z, indB1, nfold, lambdas, Winit, sparseRho, B1,
             warm = 0, parallel = TRUE, nCores = 10, 
             convergence = 1e-10, maxiter = 10000, rho=2.5)
res2 <- estB(seed, Alist, Z, indB2, nfold, lambdas, Winit, sparseRho, B2,
             warm = 0, parallel = TRUE, nCores = 10, 
             convergence = 1e-10, maxiter = 10000, rho=2.5)
res3 <- estB(seed, Alist, Z, indB3, nfold, lambdas, Winit, sparseRho, B3,
             warm = 0, parallel = TRUE, nCores = 10, 
             convergence = 1e-10, maxiter = 10000, rho=2.5)
res4 <- estB(seed, Alist, Z, indB4, nfold, lambdas, Winit, sparseRho, B4,
             warm = 0, parallel = TRUE, nCores = 10, 
             convergence = 1e-10, maxiter = 10000, rho=2.5)




