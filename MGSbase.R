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

