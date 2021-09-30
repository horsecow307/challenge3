# sampX <- matrix(nrow = 600,rnorm(360000))
# library(parallel)
# library(foreach)

# MGS_1 <- function(x){
#   # Computes the orthonormal matrix U with the same column
#   # space as X, using the modified Graham Schmidt algorithm
#   U <- x
#   coln <- ncol(x)
#   rown <- nrow(x)
#   R = matrix(nrow = coln, ncol = coln, 0)
#   Q = matrix(nrow = rown, ncol = coln, 0)
#   
#   #iterative steps
#   
#   foreach (i=1:coln, .combine='c') %dopar% {
#     R[i,i] <- norm(as.matrix(U[,1]), type="2")
#     Q[,i] <- U[,i]/R[i,i]
#     foreach (j=(i + 1):coln, .combine='c') %dopar% {
#       if (j > coln) {break} 
#       else {
#         R[i,j] <- t(Q[,i])%*%U[,j]
#         U[,j] <- U[,j] - R[i,j]%*%Q[,i]
#       }
#     }
#   }
#   return(Q)
# }

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
 
# MGSbase <- function(X){
#   # Computes the orthonormal matrix U with the same column
#   # space as X, using the modified Graham Schmidt algorithm
# 
#   U<-X
#   coln <- ncol(X)
#   rown <- nrow(X)
#   R = matrix(nrow = coln, ncol = coln, 0)
#   Q = matrix(nrow = rown, ncol = coln, 0)
# 
#   #iterative steps
#   for (i in 1:coln) {
#     U[,i] <- X[,i]
#   }
#   for (i in 1:coln) {
#     R[i,i] <- sqrt(t(U[,i])%*%U[,i])
#     Q[,i] <- U[,i]/R[i,i]
#     for (j in (i + 1):coln) {
#       if (j > coln) {break}
#       else {
#         R[i,j] <- t(Q[,i])%*%U[,j]
#         U[,j] <- U[,j] - R[i,j]%*%Q[,i]
#       }
#     }
#   }
#   return(Q%*%R)
# }
# 
# system.time(MGSbase(sampX)) #1.92 sec for 500x500, 11.50 sec for 1000x1000
# R <- parApply(cl, U, 2, FUN = function(x) norm(as.matrix(x), type="2"))
# Q <- mapply(function(x,y) x/y, U,diag(Q))
# 
# apply(cl,U,2, FUN = stepOne)