library(parallel)
library(foreach)

# x <- matrix(nrow=1000,rnorm(50000))
## parallel function ##

MGS<- function(x){
  nnode <- parallel::detectCores()
  cl <- parallel::makeCluster((nnode-1),type="PSOCK")
  doParallel::registerDoParallel(cl = cl)
  # Computes the orthonormal matrix U with the same column
  # space as X, using the modified Graham Schmidt algorithm
  U <- x
  coln <- ncol(x)
  rown <- nrow(x)
  R = matrix(nrow = coln, ncol = coln, 0)
  Q = matrix(nrow = rown, ncol = coln, 0)

  #iterative steps
  
  foreach (i=1:coln, .combine='c') %dopar% {
    R[i,i] <- norm(as.matrix(U[,i]), type="2")
    Q[,i] <- U[,i]/R[i,i]
    for (j in (i+1):coln) {
      if (j > coln) {break} else {
        R[i,j] <- t(Q[,i])%*%U[,j]
        U[,j] <- U[,j] - R[i,j]%*%Q[,i]
      }
    }
  }
  
  return(Q)
}

# MGS(x)
# 
# # Diagnostics
# 
# source("MGSbase.R")
# 
# # test speed difference as N ->
# 
# matrix(nrow=1000, ncol = 6)
# truth = norm(matlib::QR(x)$Q) # reliable matlib QR decomposer
# truth - norm(MGS(x))
# truth - norm(MGSbase(x))
# 
# microbenchmark(MGSbase(x), unit = "us")       
# microbenchmark(MGS(x), unit = "us")
