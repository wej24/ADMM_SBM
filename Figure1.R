##########################################################
# Figure 1 in the main text
# Note: the following code is for number of nodes = 1000
##########################################################

# required packages
library(igraph)

# load help function
source('help.R')

#-------------------------------------------------------------------------------
# Generate data
K <- 2
numGraph <- 1
numNode <- 1000
Bmat <- matrix(sqrt(3)/8,K,K)
diag(Bmat) <- c(3/4, 1/2)
rhoVec <- c(1/numNode, sqrt(log(numNode))/numNode,
            log(numNode)/numNode, 1/sqrt(numNode),1)
block.sizes <- c(1/4, 3/4) * numNode
memb_index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
Z <- matrix(0,numNode,K)
Z[cbind(1:numNode,memb_index)] <- 1
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute relative error
loss <- NULL
for (rho in rhoVec) {
  cur_loss <- matrix(0, 100, numNode)
  for (seed in 1:100) {
    A <- Afunc(seed, numGraph, numNode, rho, Bmat, block.sizes)[[1]]
    svdA <- svd(A)
    for (i in 1:numNode) {
      Atrun <- AtruncFunc(i, svdA)
      Bhat <- (t(Z)/colSums(Z)) %*% Atrun %*% t(t(Z)/colSums(Z))
      Bhat <- Bhat/rho
      cur_loss[seed,i] <- sqrt(sum((Bhat - Bmat)^2) / sum(Bmat^2))
    }
  }
  loss <- rbind(loss, colMeans(cur_loss))
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# plot Figure 1
graphics.off()
png(filename = "reErrNode1000.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
plot(loss[1,], type="l", lty=1, lwd=5,ylim=c(min(loss),max(loss)), 
     col="red",xlab="Number of leading eigenvectors", ylab="Mean relative error",
     cex.lab=2.5,cex.axis=2.5)
lines(loss[2,], type="l", lty=2, lwd=5,col="blue")
lines(loss[3,], type="l", lty=3, lwd=5,col="pink")
lines(loss[4,], type="l", lty=4, lwd=5,col="purple")
lines(loss[5,], type="l", lty=5, lwd=5)
points(which.min(loss[1,]), min(loss[1,]),pch=18,cex=4,lwd=5,col="red")
points(which.min(loss[2,]), min(loss[2,]),pch=18,cex=4,lwd=5,col="blue")
points(which.min(loss[3,]), min(loss[3,]),pch=18,cex=4,lwd=5,col="pink")
points(which.min(loss[4,]), min(loss[4,]),pch=18,cex=4,col="purple")
points(which.min(loss[5,]),min(loss[5,]),pch=18,cex=4,lwd=5)
legend("topright", lty=c(1,2,3,4,5),lwd=5,cex=2,
       c(expression(paste(rho,"=1/n")),
         expression(paste(rho,"=",sqrt("logn"),"/n")),
         expression(paste(rho,"=logn/n")),
         expression(paste(rho,"=1/",sqrt("n"))),
         expression(paste(rho,"=1"))),
       col=c("red","blue","pink","purple","black"),ncol=2)
dev.off()
#-------------------------------------------------------------------------------