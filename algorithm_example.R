#########################################
# Simple example for how to use algorithm
#########################################

# required packages for simluating data
library(igraph)

# required packages for the main algorithm
library(doParallel)
library(Matrix)

# load help function for simulating data
source('Mono_help.R')

#-------------------------------------------------------------------------------
# Generate the data
seed <- 1
K <- 10
numGraph <- 100 ## number of graphs L
numNode <- 1000 ## number of nodes n
block.sizes <- c(rep(0.15,2),rep(0.1,3), rep(0.08,5)) * numNode
sparseRho <- 0.1 ## sparse factor in SBM
u <- 0.9^(1:K)
Bmat <- u %*% t(u) # True B matrix
memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
Z <- matrix(0,numNode,K)
Z[cbind(1:numNode,memb_index)] <- 1
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Run the algorithm for the sequence of lambdas using function cvSRL_all

# load help function 
source('cvADMM_multiple.R')

###############################################################################
# Function: cvSRL_all
# Required Input:
#   seed: random seed to reproduce the results
#   nfold: number of folds for CV
#   Alist: a list of adjacency matrices
#   Z: membership matrix
#   sparseRho: sparse factor for SBM, default is 1
#   lambdas: a sequence of tuning parameters
#   rho: a positive constant variable rho1 in Algorithm 1, default is 2.5
#   Winit: initial value for B matrix, default is NULL
#   convergence: tolerance of the algorithm, default is 1e-10
#   maxiter: maximum number of iterations, default is 10000
#   warm: whether to use warm start, 0 is no and 1 is yes
#   parallel: whether to use parallel computing, default is TRUE
#   nCores: number of cores to use for parallel computing
# Output:
#   estB: estimated B matrix
#   lambda: the selected lambda to estimate B matrix
###############################################################################

lambdas <- 10^seq(-4,1,length.out=150) ## tuning parameters
cvRes <- cvSRL_all(seed=seed, nfold=5, Alist=Alist, Z=Z, sparseRho=sparseRho,
                   lambdas=lambdas, rho=2.5, Winit=NULL,
                   convergence = 1e-10, maxiter=10000,
                   warm=0, parallel=TRUE, nCores=5)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Run the algorithm for a single lambda using function cvSRL_all

# load help function 
source('ADMM.R')

###############################################################################
# Function: srl
# Required Input:
#   S: estimated matrix for E(A)
#   X: membership matrix
#   lambda: a single tuning parameter
#   sparseRho: sparse factor in SBM, default is 1
# Output:
#   B: estimated B matrix
###############################################################################

Abar <- Reduce("+", Alist) / length(Alist)
res <- srl(S=Abar, X=Z, lambda=cvRes$lambda, sparseRho=sparseRho)
#-------------------------------------------------------------------------------
