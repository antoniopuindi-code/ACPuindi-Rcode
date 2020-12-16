# ACPuindi-Rcode

##
matPhi <- function(smallphi=NULL, SeasPeriod=NULL, tauS) {
  # r code for transition matrix
  Phi <- matrix(c(1,smallphi), nrow=1, ncol=2)
  if(!is.null(SeasPeriod)) {
    zero.tau <- matrix(0, nrow = 1, ncol = tauS)
    Phi <- cbind(Phi, zero.tau)
  }
  # 
  beta.row <- matrix(c(0, smallphi), nrow=1, ncol=2)
  if(!is.null(SeasPeriod)) {
    beta.row <- cbind(beta.row, zero.tau)
  }
  Phi <- rbind(Phi, beta.row)
  #
  if (!is.null(SeasPeriod)) {
    seasonal.row <- t(zero.tau)
    if (!is.null(beta)) {
      seasonal.row <- cbind(seasonal.row, seasonal.row)
    }
    #
    for (i in SeasPeriod) {
      if (i == SeasPeriod[1]) {
        a.row.one <- matrix(0, nrow = 1, ncol = i)
        a.row.one[i] <- 1
        a.row.two <- cbind(diag((i - 1)), matrix(0, nrow = (i - 1), ncol = 1))
        A <- rbind(a.row.one, a.row.two)
      } else {
        old.A.rows <- dim(A)[1]
        old.A.columns <- dim(A)[2]
        a.row.one <- matrix(0, nrow = 1, ncol = i)
        a.row.one[i] <- 1
        a.row.two <- cbind(diag((i - 1)), matrix(0, nrow = (i - 1), ncol = 1))
        Ai <- rbind(a.row.one, a.row.two)
        A <- rbind(A, matrix(0, nrow = dim(Ai)[1], ncol = old.A.columns))
        A <- cbind(A, matrix(0, nrow = dim(A)[1], ncol = dim(Ai)[2]))
        A[((old.A.rows + 1):(old.A.rows + dim(Ai)[1])), ((old.A.columns + 1):(old.A.columns + dim(Ai)[2]))] <- Ai
      }
    }
    seasonal.row <- cbind(seasonal.row, A)
    Phi <- rbind(Phi, seasonal.row)
  }
  return(Phi)
}
##

##
make_A = function(smallphiS=NULL, SPeriodS=NULL){
  # r code for the observation matrix
  adjustPhi = 0
  numCols = 1 
  numSeasonal = 0
  lengthSeasonal = 0
  smallphi = 0
  SeasPeriod = 0
  #
  if(!is.null(smallphiS)) {
    smallphi = as.numeric(smallphiS)
    adjustPhi = 1
    numCols = numCols + 1
  }
  if(!is.null(SPeriodS)) {
    SeasPeriod = as.integer(SPeriodS)
    numSeasonal = length(SPeriodS)
    for(s in 1:numSeasonal) {
      lengthSeasonal = lengthSeasonal + SeasPeriod[s]
    }
    numCols = numCols + lengthSeasonal
  } else {
    lengthSeasonal = 0
  }
  #
  mA = matrix(0, nrow = 1, ncol = numCols)
  #
  if(!is.null(SPeriodS)) {
    position = adjustPhi
    for(s in 1:numSeasonal) {
      position = position + SeasPeriod[s]
      mA[1,position] = 1
    }
  }
  #
  mA[1,1] = 1
  #
  if(adjustPhi == 1) {
    mA[1,2] = smallphi
  }
  return(mA)
}
##
##
insert.diag <- function(matQ, varseas, start=c(1,1), dir=c(1,1)) {
  # r code for diagonal matrix
  sq <- seq_along(varseas)-1
  indices <- sapply(1:2, function(i) start[i] + dir[i]*sq)
  stopifnot(all(indices>0))
  stopifnot(all(indices[,1]<=nrow(matQ)))
  stopifnot(all(indices[,2]<=ncol(matQ)))  
  #
  matQ[indices] <- varseas
  matQ
}

