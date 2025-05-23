######################################################################
# M-fold cross-validation for multiple A
#####################################################################

# for real data to use repeated cross validation
rep_CVFunc <- function(fold_list, curFold, Alist, Z, lambdas, rho, Winit, warm,
                       convergence = 1e-10, maxiter = 10000) {
  indValid <- fold_list[[curFold]]
  Avalid <- Reduce("+", Alist[indValid]) / length(indValid)
  Atrain <- Reduce("+", Alist[-indValid]) / length(Alist[-indValid])
  res <- gridLambdaSRL(Atrain,Avalid,Z,lambdas, rho, Winit,
                       convergence = convergence, maxiter=maxiter, warm)
  return(res$MSEA)
}


cvSRL_all <- function(seed, nfold, Alist, Z, sparseRho=1,
                      lambdas, rho=2.5, Winit = NULL,
                      convergence = 1e-10, maxiter = 10000,
                      warm=0, parallel=TRUE, nCores = future::availableCores()) {
  set.seed(seed)
  Abar <- Reduce("+", Alist) / length(Alist)
  K <- ncol(Z)
  if (is.null(Winit)) {
    Winit <- matrix(0, K, K) 
  }
  cvRes <- cvSRL2(seed, nfold, Alist, Z, 
                  lambdas, rho = rho, Winit = Winit,
                  convergence = convergence, maxiter = maxiter,
                  warm = warm, parallel = parallel, nCores = nCores)
  if (warm == 0) {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas[cvRes$indLambda],rho = rho,
                            Winit = Winit,
                            convergence = 1e-10, maxiter=10000,
                            warm = 0)
    estB <- resSRL$estW[[1]] / sparseRho
  } else {
    resSRL <- gridLambdaSRL(Abar,Abar,Z,lambdas,rho = rho,
                            Winit = Winit,
                            convergence = 1e-10, maxiter=10000,
                            warm = 1)
    estB <- resSRL$estW[[cvRes$indLambda]] / sparseRho
  }
  return(list(estB=estB, lambda=lambdas[cvRes$indLambda]))
}


cvSRL2 <- function(seed, nfold, Alist, Z, 
                  lambdas, rho, Winit = NULL,
                  convergence = 1e-10, maxiter = 10000,
                  warm, parallel, nCores = future::availableCores()) {
  numGraph <- length(Alist)
  fold <- foldFunc2(seed, nfold, numGraph)
  if (is.null(Winit)==FALSE) {
    Winit = Winit
  } else {
    Winit <- matrix(0, ncol(Z), ncol(Z))
  }
  if (parallel == TRUE) {
    cl <- makeCluster(nCores,type="FORK")
    registerDoParallel(cl)
    cvRes <- foreach(i = 1:nfold) %dopar% {
      cvFunc2(fold, i, Alist, Z, 
              lambdas,rho, Winit, warm,
              convergence = convergence, maxiter = maxiter)
    }
    stopCluster(cl)
  } else {
    cvRes <- list()
    for (i in 1:nfold) {
      cvRes[[i]] <- cvFunc2(fold, i, Alist, Z, 
                            lambdas,rho, Winit, warm,
                            convergence = convergence, maxiter = maxiter)
    }
  }
  cv.error <- do.call(rbind, cvRes)
  cvMean <- colMeans(cv.error)
  indLambda <- which.min(cvMean)
  return(list(cv.error = cv.error, indLambda = indLambda))
}


## foldFunc2 splits number of graphs into M folds
foldFunc2 <- function(seed, nfold, numGraph) {
  set.seed(seed)
  remainder <- numGraph %% nfold
  if (remainder == 0) {
    fold <- sample(rep(1:nfold, numGraph / nfold))
  } else {
    n1 <- numGraph - remainder
    f1 <- sample(rep(1:nfold, n1 / nfold))
    f2 <- sample(1:remainder)
    fold <- c(f1, f2)
  }
  return(fold)
}

## cvFunc2 estimates matrix B in one fold 
## for a sequence of lambda values
cvFunc2 <- function(fold, curFold, Alist, Z, 
                    lambdas,rho, Winit, warm,
                    convergence = 1e-10, maxiter = 10000) {
  indTrain <- which(fold != curFold)
  indValid <- which(fold == curFold)
  Atrain <- Reduce("+", Alist[indTrain]) / length(indTrain)
  Avalid <- Reduce("+", Alist[indValid]) / length(indValid)
  res <- gridLambdaSRL(Atrain,Avalid,Z,lambdas,rho,Winit,
                       convergence = convergence, maxiter=maxiter, warm)
  return(res$MSEA)
}

gridLambdaSRL <- function(S,S2,X,lambdas,rho,Winit,
                          convergence = 1e-10, maxiter=10000,warm) {
  estW <- list()
  diag(S2) <- 0
  MSEA <- NULL
  p <- ncol(X)
  if(is.null(p)) {
    p <- 1
    X <- as.matrix(X)
  }
  T1 <- nrow(S)
  s <- as.vector(S)
  C <- kronecker(X, X)
  Cs <- t(C) %*% s
  D <- kronecker(diag(p),diag(p))
  U <- 2*t(C)%*%C+rho*t(D)%*%D
  solveU <- solve(U)
  
  for(i in 1:length(lambdas)) {
    if(i == 1) {
      init.W <- init.Z <- Winit
      init.B2 <- matrix(0, p, p)
    }
    model <- srl2(s,C,D,solveU,Cs,p,lambdas[i],
                 rho,convergence=convergence,maxiter=maxiter,
                 init.W=init.W,init.Z=init.Z,init.B2=init.B2)
    if (warm == 1) {
      init.W <- model$W
      init.Z <- model$Z
      init.B2 <- model$B2
    } else {
      init.W <- Winit
      init.Z <- Winit
      init.B2 <- matrix(0, p, p)
    }
    estA <- X %*% model$Z %*% t(X)
    diag(estA) <- 0
    MSEA <- c(MSEA, sum((S2 - estA)^2))
    estW[[i]] <- model$Z
  }
  return(list(MSEA = MSEA, estW = estW, ind = which.min(MSEA)))
}



## ADMM algorithm for CV
srl2 <- function(s,C,D,solveU,Cs,dimW,lambda,
                rho,convergence=1e-10,maxiter=10000,
                init.W=NULL,init.Z=NULL,init.B2=NULL){
  
  p <- dimW
  # Variables Initialization
  if(is.null(init.W)==FALSE){
    oldW <- oldZ <- diag(p)
    W <- init.W
    Z <- init.Z
    B2 <- init.B2
  } else {
    W <- Z <-  matrix(0,p,p)
    oldW <- oldZ <- diag(p)
    B2 <- matrix(0,p,p)
  }
  
  criteria <- 1e10
  i <- 1
  
  # While loop for the iterations
  while(criteria > convergence && i <= maxiter){
    
    W <- updateW(s,Cs,D,solveU,Z,B2,rho,p)
    
    Z <- updateZ(W,B2,lambda,rho)
    
    B2 <- B2 + Z - W
    
    criteria <- mean((W-oldW)^2)
    
    oldW <- W
    oldZ <- Z
    i <- i+1
    
  }	
  
  if(i>maxiter){
    warning("The algorithm has not converged 
            by the specified maximum number of iteration")
  }
  
  return(list(W=W,Z=Z,B2=B2,iteration=i))	
}

######################################################
# Update W
######################################################
updateW <- function(s,Cs,D,solveU,Z,B2,rho,p){
  theta <- as.vector(Z)
  b <- as.vector(B2)
  w <- solveU%*%(2*Cs+rho*t(D)%*%(theta+b))
  W <- matrix(w,p,p)
  return(W)
}


######################################################
# Update Z
######################################################
updateZ <- function(W,B2,lambda,rho){
  C1 <- W-B2
  a <- svd(C1)
  if(length(a$d) == 1) {
    D1 <- matrix(a$d,1,1)
  } else {
    D1 <- diag(a$d)
  }
  U1 <- a$u
  V1 <- a$v  
  L <- U1%*%(pmax(D1-lambda/rho,0))%*%t(V1)			
  return(L)
}

