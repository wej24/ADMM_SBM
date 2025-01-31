###########################################################################
# Simulation for section 5.4 in the main text
# MonoSBM with true memberships, K=3, n=1000, rho=logn/n
# change sparseRho for other cases in the paper
# Note: the following code runs one replicate by setting the seed
# since we run each replicate independently on cluster. 
# If people need to change seed, please ensure ARI >= 0.8
###########################################################################

# required packages
library(igraph)
library(Matrix)
library(irlba)
library(fossil)
library(mclust)
library(gtools)

# load help function
source('cvADMM_multiple.R')
source('Multi_help.R')

#-------------------------------------------------------------------------------
# Generate data
seed <- 101
K <- 3
numGraph <- 50 ## number of graphs for each layer
numNode <- 1000 ## number of nodes n 
block.sizes <- c(1/4, 1/4, 1/2) * numNode
sparseRho <- log(numNode) / numNode ## sparse factor in SBM

# True B matrix for each group
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
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
lambdas <- 10^seq(-4,-1.5,length.out=150) # tuning parameters

# run simulation for one replicate using parallel computing and our method
res <- multiSim(seed=seed, K=K, B1=B1, B2=B2, B3=B3, B4=B4, sparseRho=sparseRho,
                numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
                lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=0)

# check the results
res$ourRank # estimated rank by our method
res$ourMSE # estimated estimation error by our method
res$avgRank # estimated rank by simple averaging method
res$avgMSE # estimated estimation error by simple averaging method
res$ARI # ARI for estimated Z

# run simulation for one replicate using low rank approximation method
res_LR <- estBheLow(seed=seed, Alist=res$Alist, Zhat=res$Z, 
                    indB1=res$indB1, indB2=res$indB2, indB3=res$indB3, indB4=res$indB4,
                    B1=B1, B2=B2, B3=B3, B4=B4,
                    sparseRho=sparseRho, nfold=5)

# check the results
res_LR$avgMSE
res_LR$avgLRMSE
res_LR$avgLRRank

# run simulation for one replicate using spectral embedding method
res_spec <- estBTrun(seed=seed, K=K, Alist=res$Alist, Zhat=res$Z,
                     indB1=res$indB1, indB2=res$indB2, indB3=res$indB3, indB4=res$indB4,
                     d1=3, d2=1, d3=2, d4=3,
                     B1=B1, B2=B2, B3=B3, B4=B4,
                     sparseRho=sparseRho)

# check the results
res_spec$avgMSE
res_spec$specMSEd
res_spec$specMSEk
#-------------------------------------------------------------------------------