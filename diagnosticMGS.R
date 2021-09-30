## diagnostics

diagMGS <- function (data, ...) {
  l = list(...)
  numargs <- length(l)
  a <- unlist(lapply(l, FUN = function(x) mean(microbenchmark(x(data))$time)))
  q <- unlist(lapply(l, FUN = function(x) norm(x(data))))
  df <- as.data.frame(rbind("meantime"=a,"norm"=q))
  return(df)
  }

