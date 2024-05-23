## SBM model
## with deterministic block size
sbmA <- function(index, numNode, rho, Bmat, block.sizes) {
  graph <- sample_sbm(numNode,
                      rho * Bmat,
                      block.sizes,
                      directed = FALSE,
                      loops = TRUE)
  Amat <- as_adjacency_matrix(graph)
  return(Amat)
}

Afunc <- function(seed,numGraph, numNode, rho, Bmat, block.sizes) {
  set.seed(seed)
  Alist <- lapply(1:numGraph, sbmA,
                  numNode = numNode, rho = rho, 
                  Bmat = Bmat, block.sizes = block.sizes)
  return(Alist)
}




foldFunc <- function(seed, nfold, K, block.sizes) {
  set.seed(seed)
  fold <- lapply(1:K, function(x) {
    sample(rep(1:nfold, block.sizes[x] / nfold))
  })
  fold <- do.call(c, fold)
  return(fold)
}

cvFunc <- function(fold, curFold, A, Z,lambdas, Winit) {
  ind <- which(fold != curFold)
  Atrain <- A[ind, ind]
  Xtrain <- Z[ind,]
  MSEA <- gridLambdaSRLCV(ind,Atrain,Xtrain,A,Z,lambdas,Winit,maxiter=10000)
  return(MSEA)
}


## grid Search for CV
gridLambdaSRLCV <- function(ind,Strain,Xtrain,S,X,lambdas,Winit,maxiter=10000) {
  MSEA <-  NULL
  Stest <- S
  Stest[ind, ind] <- 0
  diag(Stest) <- 0
  p <- ncol(Xtrain)
  if(is.null(p)) {
    p <- 1
    Xtrain <- as.matrix(Xtrain)
  }
  T1 <- nrow(Strain)
  rho <- 2.5
  strain <- as.vector(Strain)
  C <- kronecker(Xtrain, Xtrain)
  Cs <- t(C) %*% strain
  D <- kronecker(diag(p),diag(p))
  U <- 2*t(C)%*%C+rho*t(D)%*%D
  solveU <- solve(U)
  
  for(i in 1:length(lambdas)) {
    if(i == 1) {
      init.W <- init.Z <- Winit
      init.B2 <- matrix(0, p, p)
    }
    model <- srl(strain,C,D,solveU,Cs,p,lambdas[i],
                 rho,convergence=1e-10,maxiter=maxiter,
                 init.W=init.W,init.Z=init.Z,init.B2=init.B2)
    init.W <- model$W
    init.Z <- model$Z
    init.B2 <- model$B2
    estA <- X %*% model$Z %*% t(X)
    estA[ind, ind] <- 0
    diag(estA) <- 0
    MSEA <- c(MSEA, sum((Stest - estA)^2))
  }
  return(MSEA)
}



gridLambdaSRL <- function(S,S2,X,lambdas, Winit, maxiter=10000) {
  estW <- list()
  diag(S2) <- 0
  MSEA <- NULL
  p <- ncol(X)
  if(is.null(p)) {
    p <- 1
    X <- as.matrix(X)
  }
  T1 <- nrow(S)
  rho <- 2.5
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
    model <- srl(s,C,D,solveU,Cs,p,lambdas[i],
                 rho,convergence=1e-10,maxiter=maxiter,
                 init.W=init.W,init.Z=init.Z,init.B2=init.B2)
    init.W <- model$W
    init.Z <- model$Z
    init.B2 <- model$B2
    estA <- X %*% model$Z %*% t(X)
    diag(estA) <- 0
    MSEA <- c(MSEA, sum((S2 - estA)^2))
    estW[[i]] <- model$Z
  }
  return(list(MSEA = MSEA, estW = estW))
}


## ADMM algorithm for CV
srl <- function(s,C,D,solveU,Cs,dimW,lambda,
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

