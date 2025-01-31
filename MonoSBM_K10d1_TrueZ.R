###########################################################################
# Simulation for section 5.1 in the main text
# MonoSBM with true memberships, K=10, d=1, rho=0.1, L=100, n=1000
# change numNode for other cases in the paper
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
K <- 10
numGraph <- 100 ## number of graphs L
numNode <- 1000 ## number of nodes n
block.sizes <- c(rep(0.15,2),rep(0.1,3), rep(0.08,5)) * numNode
sparseRho <- 0.1 ## sparse factor in SBM
u <- 0.9^(1:K)
Bmat <- u %*% t(u) # True B matrix
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
lambdas <- 10^seq(-4,1,length.out=150) ## tuning parameters

# run simulation for one replicate using parallel computing
res <- monoSim(seed=seed, K=K, Bmat=Bmat, sparseRho=sparseRho,
               numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
               TrueZ=TRUE, lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=0)


# check the results
res$ourRank # estimated rank by our method
res$ourMSE # estimated estimation error by our method
res$avgRank # estimated rank by simple averaging method
res$avgMSE # estimated estimation error by simple averaging method
#-------------------------------------------------------------------------------