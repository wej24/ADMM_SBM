##########################################
# Simulation Function for HeSBM
##########################################

## A square
A_square <- function(Alist) {
  L <- length(Alist)
  sum_mat <- numeric(0)
  for(i in 1:L){
    tmp <- crossprod(Alist[[i]])
    if(i == 1) {
      sum_mat <- tmp
    } else {
      sum_mat <- sum_mat + tmp
    }
  }
  diag(sum_mat) <- 0
  sum_mat <- as.matrix(sum_mat)
  return(sum_mat)
}


## estB on each Layer
simB <- function(A, Zhat) {
  Bhat <- (t(Zhat)/colSums(Zhat)) %*% A %*% t(t(Zhat)/colSums(Zhat))
  return(Bhat[upper.tri(Bhat,diag = TRUE)])
}


## estB on each Graph community

estB <- function(seed, Alist, Z, indB, nfold, lambdas, Winit, sparseRho, B,
                 warm = 0, parallel = TRUE, nCores = 10, 
                 convergence = 1e-10, maxiter = 10000, rho) {
  Abar <- Reduce("+", Alist[indB]) / length(Alist[indB])
  cvRes <- cvSRL(seed, nfold, Alist[indB], Z, 
                 lambdas, rho = rho, Winit = Winit,
                 convergence = convergence, maxiter = maxiter,
                 warm = warm, parallel = parallel, nCores = nCores)
  if (warm == 0) {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas[cvRes$indLambda],rho = 2.5,
                            Winit = Winit,
                            convergence = convergence, maxiter=maxiter,
                            warm = 0)
    estB <- resSRL$estW[[1]] / sparseRho
  } else {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas,rho = 2.5,
                            Winit = Winit,
                            convergence = convergence, maxiter=maxiter,
                            warm = 1)
    estB <- resSRL$estW[[cvRes$indLambda]] / sparseRho
  }
  ourRank <- qr(estB)$rank
  ourMSE <- sqrt(sum((estB - B)^2))
  avgB <- (t(Z)/colSums(Z)) %*% Abar %*% t(t(Z)/colSums(Z))
  avgB <- avgB/sparseRho
  avgRank <- qr(avgB)$rank
  avgMSE <- sqrt(sum((avgB-B)^2))
  
  res <- list(estB = estB, ourRank = ourRank, ourMSE = ourMSE,
              avgB = avgB, avgRank = avgRank, avgMSE = avgMSE)
}