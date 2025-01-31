###########################################################################
# Simulation for section S3.2 in the supplement
# MonoSBM with true memberships, K=2, d=1, rho=logn/n, L=100, n=1000
# change sparseRho, numNode, numGraph for other cases in the paper
# Note: the following code runs one replicate by setting the seed
# since we run each replicate independently on cluster. 
###########################################################################

# required packages
library(igraph)
library(doParallel)
library(Matrix)

# load help function
source('cvADMM_multiple.R')
source('Mono_help.R')

#-------------------------------------------------------------------------------
# Generate data
seed <- 1
K <- 2
numGraph <- 100 ## number of graphs L
numNode <- 1000 ## number of nodes n 
block.sizes <- c(1/4, 3/4) * numNode
sparseRho <- log(numNode) / numNode ## sparse factor in SBM
u1 <- c(0.7, 0.7^2)
Bmat <- u1 %*% t(u1) # True B matrix
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Run simulation
lambdas <- 10^seq(-3, 0, length.out=150) ## tuning parameters
# run simulation for one replicate using parallel computing
res <- monoSim(seed=seed, K=K, Bmat=Bmat, sparseRho=sparseRho,
               numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
               TrueZ=TRUE, lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=1)


# check the results
res$ourRank # estimated rank by our method
res$ourMSE # estimated estimation error by our method
res$avgRank # estimated rank by simple averaging method
res$avgMSE # estimated estimation error by simple averaging method
#-------------------------------------------------------------------------------