##################################
# Analysis for real data
# Section 6 in the main text
##################################

# required packages
library(mclust)
library(corrplot)
library(Matrix)

# load data
load('real_data.RData')
K <- 8

# load help functions
source('Multi_help.R')
source('cvADMM_multiple.R')

#-------------------------------------------------------------------------------
# plot heatmap for adjacency matrix Figure 3 in the main text
# change ind for different layer
ind <- 1
n <- nrow(Z)
cols <- grey(1-c(0,seq(0.2,1,length=25)))
corrplot(adj_list[[ind]],
         is.corr = FALSE, method="color", 
         col = cols,cl.pos= 'n',tl.pos='n')
for (i in n[-8]) {
  lines(c(0.5 + i, 0.5 + i), c(0.5,n[8] + 0.5), lwd=5, lty=2, col="coral4")
  lines(c(0.5,n[8] + 0.5), c(n[8] + 0.5 - i,n[8] + 0.5 - i), lwd=5, lty=2, col="coral4")
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# do clustering for layer group
uppB <- lapply(adj_list,simB, Zhat=Z)
uppB <- do.call(rbind,uppB)
mB <- Mclust(uppB, G=3)
summary(mB)
mB$classification
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# do simple averaging to estimate B for the second group
ind <- c(5,7:10)
Alist <- adj_list[ind]
A <- Reduce("+", Alist) / length(Alist)
B2avg <- (t(Z)/colSums(Z)) %*% A %*% t(t(Z)/colSums(Z))
svdB <- svd(B2avg)

# SVD plot
graphics.off()
png(filename = "svd2.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
plot(svdB$d, ylab = "Singular Values", ylim = c(0,0.8), pch=16,cex.axis = 2.5, cex.lab = 2.5, cex = 4)
dev.off()

# estimate B by our method for the second group
lambdas <- 10^seq(-4,3,by=0.005)
Winit <- matrix(0,K,K)
res2 <- gridLambdaSRL(A,A,Z,lambdas,rho=2.5,Winit,
                     convergence = 1e-10, maxiter=10000, warm=1)
fold_list <- combn(length(ind),2, simplify = FALSE)
nfold <- length(fold_list)
cvRes <- list()
for (i in 1:nfold) {
  cvRes[[i]] <- rep_CVFunc(fold_list, i, Alist, Z, lambdas, rho=2.5, Winit=matrix(0,K,K), warm=1)
}
cvRes <- do.call(rbind, cvRes)
cvMean <- colMeans(cvRes)
indLambda <- which.min(cvMean)
B2our <- res2$estW[[indLambda]]
r2 <- qr(B2our)$rank

# check the rank for each combination of train and validation samples
unlist(lapply(1:length(fold_list), function(x) qr(res2$estW[[which.min(cvRes[x,])]])$rank))
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# do simple averaging to estimate B for the third group
ind <- 1:4
Alist <- adj_list[ind]
A <- Reduce("+", Alist) / length(Alist)
B3avg <- (t(Z)/colSums(Z)) %*% A %*% t(t(Z)/colSums(Z))
svdB <- svd(B3avg)

## SVD plot
graphics.off()
png(filename = "svd3.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
plot(svdB$d, ylab = "Singular Values", ylim = c(0,0.8), pch=16,cex.axis = 2.5, cex.lab = 2.5, cex = 4)
dev.off()


# estimate B by our method for the third group
Winit <- matrix(0,K,K)
lambdas <- 10^seq(-4,3,by=0.005)
res3 <- gridLambdaSRL(A,A,Z,lambdas,rho=2.5,Winit,
                      convergence = 1e-10, maxiter=10000, warm=1)
fold_list <- combn(length(ind),2, simplify = FALSE)
nfold <- length(fold_list)
cvRes <- list()
for (i in 1:nfold) {
  cvRes[[i]] <- rep_CVFunc(fold_list, i, Alist, Z, lambdas, rho=2.5, Winit=matrix(0,K,K), warm=1)
}
cvRes <- do.call(rbind, cvRes)
cvMean <- colMeans(cvRes)
indLambda <- which.min(cvMean)
B3our <- res3$estW[[indLambda]]
r3 <- qr(B3our)$rank

# check the rank for each combination of train and validation samples
unlist(lapply(1:length(fold_list), function(x) qr(res3$estW[[which.min(cvRes[x,])]])$rank))
#-------------------------------------------------------------------------------