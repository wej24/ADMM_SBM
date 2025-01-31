################################################################################
# Simulation for section S3.1 in the supplement
# MonoSBM with true memberships and single A, K=2, d=2, rho=logn/n, n=1000
# change numNode for other cases in the paper
# Note: the following code runs one replicate by setting the seed
# since we run each replicate independently on cluster.
################################################################################

# required packages
library(igraph)
library(doParallel)
library(Matrix)

# load help function
source('Mono_single_help.R')

#-------------------------------------------------------------------------------
# Generate data
seed <- 2
numNode <- 1000
numGraph <- 1
K <- 2
Bmat <- matrix(0.3, K, K)
diag(Bmat) <- c(0.5, 0.4)
sparseRho <- log(numNode) / numNode
block.sizes <- c(1/4, 3/4) * numNode
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
lambdas <- 10^seq(-3, 0, length.out=150)

# run simulation for one replicate using parallel computing and our method
res <- MonoSBM_single(seed=seed, K=K, Bmat=Bmat, sparseRho=sparseRho,
                      numNode=numNode, numGraph=numGraph, block.sizes=block.sizes, 
                      lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5)

# check the results
res$ourRank
res$avgRank
res$ourMSE
res$avgMSE
#-------------------------------------------------------------------------------