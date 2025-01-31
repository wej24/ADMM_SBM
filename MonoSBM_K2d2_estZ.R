###################################################################################
# Simulation for section S3.3 in the supplement
# MonoSBM with estimated memberships, K=2, d=2, rho=sqrt(log(n))/n, L=100, n=1000
# change sparseRho, numGraph for other cases in the paper
# Note: the following code runs one replicate by setting the seed
# since we run each replicate independently on cluster.
# If people need to change seed, please ensure ARI >= 0.8
##################################################################################

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
seed <- 2
K <- 2
numGraph <- 100 ## number of graphs L
numNode <- 1000 ## number of nodes n 
block.sizes <- c(1/4, 3/4) * numNode
sparseRho <- sqrt(log(numNode)) / numNode ## sparse factor in SBM
Bmat <- matrix(0.2, K, K)
diag(Bmat) <- c(0.5, 0.5) # True B matrix
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
lambdas <- 10^seq(-4, -1.5, length.out=150) # tuning parameters

# run simulation for one replicate using parallel computing and our method
res <- monoSim(seed=seed, K=K, Bmat=Bmat, sparseRho=sparseRho,
               numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
               TrueZ=FALSE, lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=0)

# check the results
res$ourRank # estimated rank by our method
res$ourMSE # estimated estimation error by our method
res$avgRank # estimated rank by simple averaging method
res$avgMSE # estimated estimation error by simple averaging method
res$ARI # ARI for estimated Z

# run simulation for one replicate using low rank approximation method
res_LR <- estBLow(seed=seed,Alist=res$Alist,Zhat=res$Z,K=K,B=Bmat,sparseRho=sparseRho,nfold=5)

# check the results
res_LR$avgMSE # estimated estimation error by simple averaging method
res_LR$avgLRMSE # estimated estimation error by low rank approximation method
res_LR$avgLRRank # estimated rank by low rank approximation method

# run simulation for one replicate using spectral embedding method
res_spec <- specEmbed(d=2, K=K, Abar=res$Abar, Zhat=res$Z, B=Bmat, sparseRho=sparseRho)

# check the results
res_spec$specMSEd
res_spec$specMSEk
#-------------------------------------------------------------------------------