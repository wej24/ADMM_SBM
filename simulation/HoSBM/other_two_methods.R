## Spectral embedding method and low rank approximation method
## Inputs:
## dhat : desired rank
## Abar: averaged adjacency matrices
## Zhat: estimated membership matrix
## sparseRho: sparse factor in SBM, default is 1

library(irlba)

specEmbed <- function(dhat, Abar, Zhat, sparseRho=1) {
  svdA <- svdr(Abar, dhat + 1)
  if (dhat == 1) {
    Ad <- svdA$d[1] * svdA$u[,1] %*% t(svdA$v[,1])
  } else {
    Ad <- svdA$u[,1:dhat] %*% diag(svdA$d[1:dhat]) %*% t(svdA$v[,1:dhat])
  }
  
  embed_B <- (t(Zhat)/colSums(Zhat)) %*% Ad %*% t(t(Zhat)/colSums(Zhat))
  embed_B <- embed_B / sparseRho
  return(embed_B)
}

lowRank <- function(dhat, Abar, Zhat, sparseRho=1) {
  Bavg <- (t(Zhat)/colSums(Zhat)) %*% Abar %*% t(t(Zhat)/colSums(Zhat))
  Bavg <- Bavg/sparseRho
  
  svdB <- svdr(Bavg, min(dhat+1,ncol(Zhat)))
  if (d == 1) {
    Blow <- svdB$d[1] * svdB$u[,1] %*% t(svdB$v[,1])
  } else {
    Blow <- svdB$u[,1:dhat] %*% diag(svdB$d[1:dhat]) %*% t(svdB$v[,1:dhat])
  }
  return(Blow)
}