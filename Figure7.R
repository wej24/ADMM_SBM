##########################################################
# Figure 7 in the supplement
##########################################################

# required packages
library(igraph)
library(doParallel)
library(Matrix)
library(irlba)
library(fossil)
library(mclust)
library(gtools)
library(corrplot)

# load help function
source('Mono_help.R')

#-------------------------------------------------------------------------------
# Generate data
K <- 10
numGraph <- 100 # number of graphs L
numNode <- 10000 # number of nodes n
block.sizes <- c(rep(0.15,2),rep(0.1,3), rep(0.08,5)) * numNode
sparseRho <- 0.1 # sparse factor in SBM
u <- 0.9^(1:K)
Bmat <- u %*% t(u) # True B matrix
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
set.seed(115)
memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
Alist <- Afunc(seed, numGraph, numNode, sparseRho, Bmat, block.sizes)
Abar <- Reduce("+", Alist) / length(Alist)
## do clustering to estimate Z
svdA <- svdr(Abar, K+1)
U <- svdA$u[,1:K] %*% diag(svdA$d[1:K])
mG <- Mclust(U, G=K)
ARI <- adj.rand.index(memb_index, mG$classification)
memb2 <- align_two_membership_vectors(memb_index, mG$classification)
Z <- matrix(0,numNode,K)
Z[cbind(1:numNode,memb2)] <- 1

Bavg <- (t(Z)/colSums(Z)) %*% Abar %*% t(t(Z)/colSums(Z))
Bavg <- Bavg/sparseRho

svdB <- svd(Bavg)
BavgLR <- svdB$d[1] * svdB$u[,1] %*% t(svdB$v[,1])

diffavg <- abs(Bavg-Bmat)
diffavgLR <- abs(BavgLR - Bmat)

# scree plot of Bavg
plot(svdB$d, ylab = "Singular Values", ylim = c(0,4.5), pch=16,cex.axis = 2.5, cex.lab = 2.5, cex = 4)

# heatmaps of diffavg
cols <- grey(1-c(0,seq(0.7,1,length=5)))
cols <- grey(c(1,0.9,0.3,0))
corrplot(diffavg,
         is.corr = FALSE, method="color", 
         col = cols,col.lim=c(0,0.7), tl.pos='n',cl.cex = 2.5)
for (i in 1:9) {
  lines(c(0.5 + i, 0.5 + i), c(0.5,10 + 0.5), lwd=5, lty=2, col="coral4")
  lines(c(0.5,10 + 0.5), c(0.5 + i,0.5 + i), lwd=5, lty=2, col="coral4")
}

# heatmaps of diffavgLR
cols <- grey(1-c(0,seq(0.7,1,length=5)))
cols <- grey(c(1,0.9,0.3,0))
corrplot(diffavgLR,
         is.corr = FALSE, method="color", 
         col = cols,col.lim=c(0,0.7), tl.pos='n',cl.cex = 2.5)
for (i in 1:9) {
  lines(c(0.5 + i, 0.5 + i), c(0.5,10 + 0.5), lwd=5, lty=2, col="coral4")
  lines(c(0.5,10 + 0.5), c(0.5 + i,0.5 + i), lwd=5, lty=2, col="coral4")
}
#-------------------------------------------------------------------------------