###########################################################################
# Simulation for section 5.3 in the main text
# MonoSBM with re-estimated Z, K=10, d=2, rho=0.15, L=50, n=1000
# Note: the following code runs one replicate by setting the seed
# since we run each replicate independently on cluster.
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
numGraph <- 50 # number of graphs L
numNode <- 1000 # number of nodes n
block.sizes <- c(rep(0.15,2),rep(0.1,3), rep(0.08,5)) * numNode
sparseRho <- 0.15 # sparse factor in SBM
Bmat <- geneB(2,K) # True B matrix
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
lambdas <- 10^seq(-4,4,by = 0.005) # tuning parameters

# run simulation for one replicate using parallel computing
res <- reEstZSim(seed=seed, K=K, Bmat=Bmat, sparseRho=sparseRho,
                 numNode=numNode, numGraph=numGraph, block.sizes=block.sizes,
                 lambdas=lambdas, nfold=5, parallel=TRUE, ncores=5, warm=1)


res$mic # misclustering error rate without re-estimation
res$mic_refit # misclustering error rate with re-estimation
res$err # estimation error without re-estimation
res$err_refit # estimation error with re-estimation
res$ourRank # estimated rank by our method
#-------------------------------------------------------------------------------