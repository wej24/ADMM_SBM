
load("main_analysis_revision.RData")
names(adj_list) <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")


clustering_res2 <- clustering_res
order_vec <- c(3, 7, 1, 2, 6, 8, 4, 5)
for(i in 1:length(order_vec)){
  clustering_res2[which(clustering_res == i)] <- order_vec[i]
}


## estB on each Layer
simB <- function(A, Zhat) {
  Bhat <- (t(Zhat)/colSums(Zhat)) %*% A %*% t(t(Zhat)/colSums(Zhat))
  return(Bhat[upper.tri(Bhat,diag = TRUE)])
}

orderGroup <- order(clustering_res2)

Z <- matrix(0, length(clustering_res2), K)
Z[cbind(1:length(clustering_res2), clustering_res2[orderGroup])] <- 1

adj_list2 <- lapply(1:10,function(x) adj_list[[x]][orderGroup, orderGroup])
uppB <- lapply(adj_list2,simB, Zhat=Z)
uppB <- do.call(rbind,uppB)
plot(svd(uppB)$d)


library(mclust)
mB <- Mclust(uppB, G=3)
summary(mB)
mB$classification

## second group

ind <- c(5,7:10)
Alist <- adj_list2[ind]
A <- Reduce("+", Alist) / length(Alist)
B2avg <- (t(Z)/colSums(Z)) %*% A %*% t(t(Z)/colSums(Z))
svdB <- svd(B2avg)

## SVD plot
graphics.off()
png(filename = "svd1.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
plot(svdB$d, ylab = "Singular Values", ylim = c(0,0.8), pch=16,cex.axis = 2.5, cex.lab = 2.5, cex = 4)
dev.off()


cols <- grey(1-c(0,seq(0.2,1,length=20)))
library(corrplot)
graphics.off()
png(filename = "t2Rank8.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
corrplot(B2avg,is.corr=FALSE,method="color", col=cols, tl.pos='n',cl.cex = 2.5)
lines(c(1.5, 1.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(2.5, 2.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(3.5, 3.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(4.5, 4.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(5.5, 5.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(6.5, 6.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(7.5, 7.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(1.5, 1.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(2.5, 2.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(3.5, 3.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(4.5, 4.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(5.5, 5.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(6.5, 6.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(7.5, 7.5), lwd=5, lty=2, col = "coral4")
dev.off()


## third group
ind <- 1:4
Alist <- adj_list2[ind]
A <- Reduce("+", Alist) / length(Alist)
B3avg <- (t(Z)/colSums(Z)) %*% A %*% t(t(Z)/colSums(Z))

svdB <- svd(B3avg)

## SVD plot
graphics.off()
png(filename = "svd2.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
plot(svdB$d, ylab = "Singular Values", ylim = c(0,0.8), pch=16,cex.axis = 2.5, cex.lab = 2.5, cex = 4)
dev.off()

cols <- grey(1-c(0,seq(0.2,1,length=25)))

graphics.off()
png(filename = "t3Rank8.png", width = 10, height = 10, units = "in", res = 100)
par(pty="s")
corrplot(B3avg,is.corr=FALSE,method="color", col=cols, tl.pos='n',cl.cex = 2.5)
lines(c(1.5, 1.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(2.5, 2.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(3.5, 3.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(4.5, 4.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(5.5, 5.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(6.5, 6.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(7.5, 7.5), c(0.5, 8.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(1.5, 1.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(2.5, 2.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(3.5, 3.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(4.5, 4.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(5.5, 5.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(6.5, 6.5), lwd=5, lty=2, col = "coral4")
lines(c(0.5, 8.5), c(7.5, 7.5), lwd=5, lty=2, col = "coral4")
dev.off()




## estimate B by our method
source("../../Algorithm/cvADMM.R")
## Repeated CV function for real data
cvFunc <- function(fold_list, curFold, Alist, Z, lambdas, rho, Winit, warm,
                   convergence = 1e-10, maxiter = 10000) {
  indTrain <- fold_list[[curFold]]
  Avalid <- Reduce("+", Alist[indTrain]) / length(indTrain)
  Atrain <- Reduce("+", Alist[-indTrain]) / length(Alist[-indTrain])
  res <- gridLambdaSRL(Atrain,Avalid,Z,lambdas, rho, Winit,
                       convergence = convergence, maxiter=maxiter, warm)
  return(res$MSEA)
}

library(Matrix)
time <- 2 ## time = 1 for the second group and 2 for the third group
warm <- 1
if (time == 2) {
  ind <- 1:4
} else {
  ind <- c(5,7:10)
}
Alist <- adj_list2[ind]
A <- Reduce("+", Alist) / length(Alist)
Winit <- matrix(0,K,K)
lambdas <- 10^seq(-4,3,by=0.005)
res <- gridLambdaSRL(A,A,Z,lambdas,rho=2.5,Winit,
                     convergence = 1e-10, maxiter=10000, warm=warm)
fold_list <- combn(length(ind),2, simplify = FALSE)
nfold <- length(fold_list)
cvRes <- list()
for (i in 1:nfold) {
  cvRes[[i]] <- cvFunc(fold_list, i, Alist, Z, lambdas, rho=2.5, Winit=matrix(0,K,K), warm=warm)
}
resSummary <- list(Bour=Bour, r=r, resSRL=res, cvRes=cvRes)
cvMean <- colMeans(cvRes)
which.min(cvMean)
unlist(lapply(1:length(fold_list), function(x) qr(resSummary$resSRL$estW[[ind[x]]])$rank))
apply(cvRes,1, which.min)
## Bour with one combination of train and validation gives rank 2
B3 <- resSummary$resSRL$estW[[688]]

## replace time = 2 with time = 1 will give the estimation for the second group
B2 <- resSummary$resSRL$estW[[794]]
