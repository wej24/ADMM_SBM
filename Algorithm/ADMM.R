############################################################
# ADMM algorithm named as srl (Algorithm 1 in paper) 
# Required Input:
#   S:  n * n matrix for Abar
#   X:  n * K matrix for membership matrix
#   lambda: tuning parameter
#   rho: a positive constant variable rho1 in Algorithm 1
# Output:
#   W and Z are the estimations of probability matrix B
#   Since the nuclear norm is on Z, we use Z as the final 
#   estimation of probability matrix B
###########################################################

srl <- function(S,X,lambda,
                rho,convergence=1e-10,maxiter=10000,
                init.W=NULL,init.Z=NULL,init.B2=NULL){
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

