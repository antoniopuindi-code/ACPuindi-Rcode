# ACPuindi-Rcode
insert.diag <- function(matQ, varianceSeas, start=c(1,1), dir=c(1,1)) {
  
  sq <- seq_along(varianceSeas)-1
  indices <- sapply(1:2, function(i) start[i] + dir[i]*sq)
  stopifnot(all(indices>0))
  stopifnot(all(indices[,1]<=nrow(matQ)))
  stopifnot(all(indices[,2]<=ncol(matQ)))  
  
  matQ[indices] <- varianceSeas
  matQ
}

