###########################################################################
# Simulation for section 5.2 in the main text
# MonoSBM with estimated memberships, K=10, d=1, rho=0.1, L=100, n=10000
# Note: the following code runs one replicate by setting the seed
# since we run each replicate independently on cluster.
# If people need to change seed, please ensure ARI >= 0.8
########################################################################### 

# required packages
library(igraph)
library(doParallel)
library(Matrix)
library(irlba)
library(fossil)
library(mclust)
library(gtools)

# load help function
source('cvADMM_multiple.R')
source('Mono_help.R')

#-------------------------------------------------------------------------------
# Generate data
seed <- 1 
K <- 10
numGraph <- 100 # number of graphs L
numNode <- 10000 # number of nodes n
block.sizes <- c(rep(0.15,2),rep(0.1,3), rep(0.08,5)) * numNode
sparseRho <- 0.1 # sparse factor in SBM
u <- 0.9^(1:K)
Bmat <- u %*% t(u) # True B matrix
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
lambdas <- 10^seq(-4,1,length.out=150) # tuning parameters

# run simulation for one replicate using parallel computing and our method
res <- monoSim(seed=seed, K=K, Bmat=Bmat, sparseRho=sparseRho,
               numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
               TrueZ=FALSE, lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=0)

# run simulation for one replicate using low rank approximation method with true d
Blow <- lowRank(d=1, Abar=res$Abar, Zhat=res$Z, sparseRho=sparseRho)
res$avgLRMSE <- sqrt(sum((Blow-Bmat)^2))

# check the results
res$ourRank # estimated rank by our method
res$ourMSE # estimated estimation error by our method
res$avgRank # estimated rank by simple averaging method
res$avgMSE # estimated estimation error by simple averaging method
res$avgLRMSE # estimated estimation error by low rank approximation method with d=1
#-------------------------------------------------------------------------------