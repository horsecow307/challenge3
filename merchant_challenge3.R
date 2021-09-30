sampX <- matrix(nrow = 600,rnorm(90000))
# library(parallel)
# library(foreach)
# library(matlib)
# 
# truth <- matlib::QR(sampX)$Q
# a <- MGSbase(sampX)
# b <- MGS(sampX)
# norm(a-b,type="2") # a and b are off but close
# norm(truth-b,type="2") # b and truth are off
# norm(truth-a,type="2") # a and truth are off
system.time(MGS(sampX))
system.time(MGSbase(sampX))

MGS <- function(x){
  # Computes the orthonormal matrix U with the same column
  # space as X, using the modified Graham Schmidt algorithm
  U <- x
  coln <- ncol(x)
  rown <- nrow(x)
  R = matrix(nrow = coln, ncol = coln, 0)
  Q = matrix(nrow = rown, ncol = coln, 0)

  #iterative steps

  for (i in 1:coln) {
    R[i,i] <- norm(as.matrix(U[,1]), type="2")
    Q[,i] <- U[,i]/R[i,i]
    for (j in (i + 1):coln) {
      if (j > coln) {break}
      else {
        R[i,j] <- t(Q[,i])%*%U[,j]
        U[,j] <- U[,j] - R[i,j]%*%Q[,i]
      }
    }
  }
  return(Q)
}
 
MGSbase <- function(X){
  # Computes the orthonormal matrix U with the same column
  # space as X, using the modified Graham Schmidt algorithm

  U<-X
  coln <- ncol(X)
  rown <- nrow(X)
  R = matrix(nrow = coln, ncol = coln, 0)
  Q = matrix(nrow = rown, ncol = coln, 0)

  #iterative steps
  for (i in 1:coln) {
    U[,i] <- X[,i]
  }
  for (i in 1:coln) {
    R[i,i] <- sqrt(t(U[,i])%*%U[,i])
    Q[,i] <- U[,i]/R[i,i]
    for (j in (i + 1):coln) {
      if (j > coln) {break}
      else {
        R[i,j] <- t(Q[,i])%*%U[,j]
        U[,j] <- U[,j] - R[i,j]%*%Q[,i]
      }
    }
  }
  return(Q)
}

