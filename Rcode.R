##
##
scov <- function (y, input, use.damped.trend = NULL, seasonal.periods = NULL) {
  
  if (any(class(y) %in% c("data.frame", "list", "matrix", "mts"))) 
    stop("y should be a univariate time series")
  seriesname <- deparse(substitute(y))
  origy <- y
  if (is.null(seasonal.periods)) {
    if (any(class(y) == "msts")) 
      seasonal.periods <- attr(y, "msts")
    else if (class(y) == "ts") 
      seasonal.periods <- frequency(y)
    else {
      y <- as.ts(y)
      seasonal.periods <- 1
    }
  }else {
    if (!any(class(y) == "ts")) 
      y <- msts(y, seasonal.periods)
  }
  seasonal.periods <- unique(pmax(seasonal.periods, 1))
  if (all(seasonal.periods == 1)) 
    seasonal.periods <- NULL
  ny <- length(y)
  y <- na.contiguous(y)
  if (ny != length(y)) 
    warning("Missing values encountered")
  
  if (!is.null(seasonal.periods)) {
    seasonal.mask <- (seasonal.periods == 1)
    seasonal.periods <- seasonal.periods[!seasonal.mask]
  }
  
  if (!is.null(use.damped.trend)) {
    use.damped.trend <- TRUE
  }else{
    use.damped.trend <- FALSE
  }
  
  y <- as.numeric(y)
  best.aic <- NULL
  
  for (damping in use.damped.trend) {
    current.model <- filtSCov(y, input, damping = damping, seasonal.periods = seasonal.periods)
    if (!is.null(best.aic)) {
      if (current.model$AIC < best.aic) {
        best.aic <- current.model$AIC
        best.model <- current.model
      }
    }else {
      best.model <- current.model
      best.aic <- best.model$AIC
    }
  }
  best.model$call <- match.call()
  attributes(best.model$fitted.values) <- attributes(best.model$errors) <- attributes(origy)
  best.model$y <- origy
  best.model$series <- seriesname
  best.model$method <- "SCov"
  return(best.model)
}
##
##
filtSCov <- function(y, input, damping, seasonal.periods, force.seasonality=FALSE) {
  
  first.model <- fitSCov(y, input, use.damping = damping, seasonal.periods = seasonal.periods)
  
  if(!is.null(seasonal.periods) && !force.seasonality) {
    non.seasonal.model <- fitSCov(y, input, use.damping = damping, seasonal.periods = NULL)
    
    if (first.model$AIC > non.seasonal.model$AIC) {
      seasonal.periods <- NULL
      first.model <- non.seasonal.model
    }
  }
  starting.params <- first.model$parameters
  
  second.model <- fitSCov(y, input, use.damping = damping, seasonal.periods = seasonal.periods)
  if (second.model$AIC < first.model$AIC) {
    return(second.model)
  } else {
    return(first.model)
  }
}
##
##
fitSCov <- function(y, input, use.damping, seasonal.periods=NULL, starting.params=NULL) {
  # r code for the model estimation process
  if (!is.null(seasonal.periods)) {
    seasonal.periods <- as.integer(sort(seasonal.periods))
  }
  if (is.null(starting.params)) {
    # Calculate starting values:
    p1 = p2 = p3 <- 0.1
    sigv <- 0.1
    sigw1 <- 0.002
    sigw2 <- 0.3
    if (use.damping) {
      small.phi <- .999
    } else {
      small.phi <- 1
    }
    if (!is.null(seasonal.periods)) {
      gamma.v <- rep(.001, length(seasonal.periods))
    } else {
      gamma.v <- NULL
    }
  } else {
    paramz <- unParamSCov(starting.params$vect, starting.params$control)
    p1 <- paramz$p1
    p2 <- paramz$p2
    p3 <- paramz$p3
    sigv <- paramz$sigv
    sigw1 <- paramz$sigw1
    sigw2 <- paramz$sigw2
    small.phi <- paramz$small.phi
    gamma.v <- paramz$gamma.v
  }
  # Optimise the starting values:
  # Make the parameter vector
  param.vector <- paramSCov(p1=p1, p2=p2, p3=p3, sigv=sigv, sigw1=sigw1, sigw2=sigw2, 
                            small.phi = small.phi, gamma.v = gamma.v)
  par.scale <- ParSCov(param.vector$control)
  
  if (!is.null(seasonal.periods)) {
    tau <- sum(seasonal.periods)
  } else {
    tau <- 0
  }
  # Optimise the likelihood function
  if (length(param.vector$vect) > 1) {
    optim.like <<- optim(par = param.vector$vect, fn = LikeSCov, y.data=y, method = "Nelder-Mead", 
                         use.small.phi = use.damping, seasonal.periods = seasonal.periods, 
                         tau = tau, control = list(maxit=(100*length(param.vector$vect)^2), 
                                                   parscale = par.scale))
  } else {
    optim.like <<- optim(par = param.vector$vect, fn = LikeSCov, y.data=y, method = "BFGS", 
                         use.small.phi = use.damping, seasonal.periods = seasonal.periods, 
                         tau = tau, control=list(parscale=par.scale))
  }
  # Get the parameters out of the param.vector
  paramz <- unParamSCov(optim.like$par, param.vector$control)
  p1 <- paramz$p1
  p2 <- paramz$p2
  p3 <- paramz$p3
  sigv <<- paramz$sigv
  sigw1 <<- paramz$sigw1
  sigw2 <<- paramz$sigw2
  small.phi <- paramz$small.phi
  gamma.v <<- paramz$gamma.v
  # One-step-ahead forecast
  fitted.values <- kf$fprev[, , 1:450]
  e <- kf$innov
  # Get the like and the AIC
  likelihood <- optim.like$value
  aic <- likelihood + 2*(length(param.vector$vect) + nrow(kf$xp))
  
  SE = sqrt(diag(solve(optim.like$hessian)))
  parErrors <- round(cbind(estimate=optim.like$par, SE), 6)
  colnames(parErrors) <- c("parameters")
  rownames(parErrors) <- c("cReg1", "cReg2", "cReg3", "sigv", "sigw1", "sigw2", 
                           "phi", "sigSeas")
  parErros <- as.table(parErrors)
  # Listing the objects
  fits <- ts(c(fitted.values))
  fits[1:4] <- y[1:4]
  e <- ts(c(e))
  tsp(fits) <- tsp(e) <- tsp(y)
  
  model.output <- list(y=y,input=input, fitted.values=fits, errors=e,
                       parameters=list(vect=optim.like$par,control=param.vector$control),
                       p1=p1, p2=p2, p3=p3, sigv=sigv, sigw1=sigw1, sigw2=sigw2, 
                       damping.parameter=small.phi, sigSeas=gamma.v, parameters=parErros, 
                       likelihood=likelihood, val.optim=optim.like$convergence, AIC=aic,
                       seasonal.periods=seasonal.periods)
  class(model.output) <- c("SCov")
  return(model.output)
}
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
##
##
matQ = function(pdim, sigw1, sigw2, epsilonS=NULL, SeasPeriodS=NULL){
  # r code for state covariance matrix
  epsilonSLength = 0
  adjustSigw2 = 0
  numCols = 1
  if(!is.null(sigw2)) {
    numCols = numCols + 1
    adjustSigw2 = 1
  }
  #
  if((!is.null(epsilonS))&&(!is.null(SeasPeriodS))) {
    epsilonSVector = as.numeric(epsilonS)
    SeasPeriod = as.integer(SeasPeriodS)
    for(i in 1:length(SeasPeriodS)) {
      epsilonSLength = epsilonSLength + SeasPeriod[i]
    }
    numCols = numCols + epsilonSLength
  } else {
    epsilonSLength = 0
  }
  
  gTranspose = matrix(0, nrow = 1, ncol = numCols)
  gTranspose[1,1] = as.numeric(sigw1)
  if(!is.null(sigw2)) {
    gTranspose[1,2] = as.numeric(sigw2)
  }
  # 
  if((!is.null(epsilonS))&&(!is.null(SeasPeriodS))) {
    position = adjustSigw2 + 1
    gTranspose[1, position] = epsilonSVector
    if(length(epsilonS) > 1) {
      for(s in 1:(length(SeasPeriodS)-1)) {
        position = position + SeasPeriod[s]
        gTranspose[1,position] = epsilonSVector[(s+1)]
      }
    }
  }
  matQ <- diag(as.numeric(gTranspose), pdim)
  return(matQ)
}
##
##
unParamSCov<-function(param.vector, control) {
  
  p1<-param.vector[1]
  p2<-param.vector[2]
  p3<-param.vector[3]
  sigv<-param.vector[4]
  sigw1<-param.vector[5]
  sigw2<-param.vector[6]
  if(control$use.damping) {
    small.phi<-param.vector[7]
    epsilonS.start<-8
  } else {
    small.phi<-1
    epsilonS.start<-7
  }
  
  if(control$length.epsilonS > 0) {
    epsilonS.vector<-param.vector[epsilonS.start:(epsilonS.start+control$length.epsilonS-1)]
    final.epsilonS.pos<-epsilonS.start+control$length.epsilonS-1
  } else {
    epsilonS.vector=NULL
    final.epsilonS.pos<-epsilonS.start-1
  }
  
  return(list(p1=p1, p2=p2, p3=p3, sigv=sigv, sigw1=sigw1, sigw2=sigw2, 
              small.phi=small.phi, epsilonS.v=epsilonS.vector))
}
##
##
paramSCov <- function(p1, p2, p3, sigv, sigw1, sigw2, small.phi=1, epsilonS.v=NULL) {
  
  param.vector<-cbind(p1, p2, p3, sigv, sigw1, sigw2)
  if (is.null(small.phi)) {
    use.damping <- FALSE
  } else if (small.phi != 1) {
    param.vector <- cbind(param.vector, small.phi)
    use.damping <- TRUE
  } else {
    use.damping <- FALSE
  }
  if (!is.null(epsilonS.v)) {
    epsilonS.v <- matrix(epsilonS.v, nrow = 1, ncol = length(epsilonS.v))
    param.vector <- cbind(param.vector, epsilonS.v)
    length.epsilonS <- length(epsilonS.v)
  } else {
    length.epsilonS <- 0
  }
  
  control <- list(use.damping = use.damping, length.epsilonS = length.epsilonS)
  return(list(vect = as.numeric(param.vector), control = control))
}
##
##
ParSCov<-function(control) {
  
  parscale<-c(rep(1e-2,3),rep(1e-2,3))
  if(control$use.damping) {
    parscale<-c(parscale, 1e-2)
  } 
  if(control$length.epsilonS > 0) {
    parscale<-c(parscale, rep(1e-2, control$length.epsilonS))
  }
  return(parscale)
}
##
##
KFSCov<-function(num, y, alpha, A, mu0, Sigma0, Phi, Ups, Gam, cQ, cR, input){
  
  Q0 = t(cQ) %*% cQ
  R0 = t(cR) %*% cR
  Phi = as.matrix(Phi)
  pdim = nrow(Phi)
  y = as.matrix(y)
  qdim = ncol(y)
  rdim = ncol(as.matrix(input))
  
  if (max(abs(Ups)) == 0) 
    Ups = matrix(0, pdim, rdim)
  if (max(abs(Gam)) == 0) 
    Gam = matrix(0, qdim, rdim)
  Ups = as.matrix(Ups)
  Gam = as.matrix(Gam)
  ut = matrix(input, num, rdim)
  xp = array(NA, dim = c(pdim, 1, num))
  Pp = array(NA, dim = c(pdim, pdim, num))
  f = array(NA, dim = c(qdim, 1, num))
  xf = array(NA, dim = c(pdim, 1, num))
  Pf = array(NA, dim = c(pdim, pdim, num))
  innov = array(NA, dim = c(qdim, 1, num))
  res = array(NA, dim = c(qdim, 1, num))
  sig = array(NA, dim = c(qdim, qdim, num))
  Qt = array(NA, dim = c(pdim, pdim, num))
  x00 = as.matrix(mu0, nrow = pdim, ncol = 1)
  P00 = as.matrix(Sigma0, nrow = pdim, ncol = pdim)
  
  xp[, , 1] = Phi %*% x00 + Ups %*% ut[1, ]
  Pp[, , 1] = Phi %*% P00 %*% t(Phi) + Q0
  
  B = matrix(A[, , 1], nrow = qdim, ncol = pdim)
  f[, , 1] = B %*% xp[, , 1] + Gam %*% ut[1, ]
  
  sigtemp = B %*% Pp[, , 1] %*% t(B) + R0
  sig[, , 1] = (t(sigtemp) + sigtemp)/2
  siginv = solve(sig[, , 1])
  K = Pp[, , 1] %*% t(B) %*% siginv
  innov[, , 1] = y[1, ] - B %*% xp[, , 1] - Gam %*% ut[1, ]
  Q = diag(diag(alpha * Q0 + (1-alpha) * (K %*% innov[, , 1] %*% t(innov[, , 1]) %*% t(K))))
  
  xf[, , 1] = xp[, , 1] + K %*% innov[, , 1]
  Pf[, , 1] = Pp[, , 1] - K %*% B %*% Pp[, , 1]
  
  res[, , 1] = y[1, ] - B %*% xf[, , 1] - Gam %*% ut[1, ]
  R = alpha * R0 + (1-alpha) * (res[, , 1] %*% t(res[, , 1]) + B %*% Pp[, , 1] %*% t(B))
  Qt[, , 1] = K %*% (B %*% Pp[, , 1] %*% t(B) + R) %*% t(K)
  sigmat = as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)
  like = log(det(sigmat)) + t(innov[, , 1]) %*% siginv %*% innov[, , 1]
  
  for (i in 2:num) {
    if (num < 2) 
      break
    xp[, , i] = Phi %*% xf[, , i - 1] + Ups %*% ut[i, ]
    Pp[, , i] = Phi %*% Pf[, , i - 1] %*% t(Phi) + Q
    
    B = matrix(A[, , i], nrow = qdim, ncol = pdim)
    f[, , i] = B %*% xp[, , i] + Gam %*% ut[i, ]
    
    siginv = B %*% Pp[, , i] %*% t(B) + R
    sig[, , i] = (t(siginv) + siginv)/2
    siginv = solve(sig[, , i])
    K = Pp[, , i] %*% t(B) %*% siginv
    
    innov[, , i] = y[i, ] - B %*% xp[, , i] - Gam %*% ut[i, ]
    Q = diag(diag(alpha * Qt[, , i-1] + (1-alpha) * (K %*% innov[, , i] %*% t(innov[, , i]) %*% t(K))))
    
    xf[, , i] = xp[, , i] + K %*% innov[, , i]
    Pf[, , i] = Pp[, , i] - K %*% B %*% Pp[, , i]
    
    res[, , i] = y[i, ] - B %*% xf[, , i] - Gam %*% ut[i, ]
    R = alpha * R0 + (1-alpha) * (res[, , i] %*% t(res[, , i]) + B %*% Pp[, , i] %*% t(B))
    Qt[, , i] = K %*% (B %*% Pp[, , i] %*% t(B) + R) %*% t(K)
    sigmat = matrix(sig[, , i], nrow = qdim, ncol = qdim)
    like = like + log(det(sigmat)) + t(innov[, , i]) %*% siginv %*% innov[, , i]
  }
  like = 0.5 * like
  list(xp = xp, Pp = Pp, fprev = f, xf = xf, Pf = Pf, like = like, 
       innov = innov, res = res, sig = sig, Qt = Q, Rt = R, Kn = K,
       Phi=Phi, Gam=Gam)
}
##
##
checkphi<-function(Phi, small.phi=NULL, sigv=NULL, sigw1=NULL, sigw2=NULL,
                   epsilonS.vector=NULL, tau=0) {
  
  if(!is.null(sigv)){
    if(sigv < 0)
      return(FALSE)
  }
  if(!is.null(sigw1)){
    if(sigw1 < 0)
      return(FALSE)
  }
  if(!is.null(sigw2)){
    if(sigw2 < 0)
      return(FALSE)
  }
  if(!is.null(epsilonS.vector)){
    if(epsilonS.vector < 0)
      return(FALSE)
  }
  if(!is.null(small.phi)){
    if(((small.phi < .8) | (small.phi > 1)))
      return(FALSE)
  }
  Phi.eigen.values<-eigen(Phi, symmetric=FALSE, only.values=TRUE)$values
  return(all(abs(Phi.eigen.values) < 1 + 1e-2))
}
##
##
LikeSCov <- function(param.vector, y.data, use.small.phi, seasonal.periods, tau) {
  # r code for the criterion function to be minimized.
  p1 <- param.vector[1]
  p2 <- param.vector[2]
  p3 <- param.vector[3]
  sigv <- param.vector[4]
  sigw1 <- param.vector[5]
  sigw2 <- param.vector[6]
  if (use.small.phi) {
    small.phi <- param.vector[7]
    epsilonS.start <- 8
  } else {
    small.phi <- 1
    epsilonS.start <- 7
  }
  
  if (!is.null(seasonal.periods)) {
    epsilonS.vector <- param.vector[epsilonS.start:(epsilonS.start + length(seasonal.periods) - 1)]
    final.epsilonS.pos <- epsilonS.start + length(epsilonS.vector) - 1
  } else {
    epsilonS.vector <- NULL
    final.epsilonS.pos <- epsilonS.start - 1
  }
  Phi <- matPhi(smallphi=small.phi, SeasPeriod=seasonal.periods, 
                tauS=tau)
  
  A <- array(make_A(smallphiS=small.phi, SPeriodS=seasonal.periods), 
             dim = c(1, nrow(Phi), num))
  #
  pdim <- nrow(Phi) 
  qdim <- ncol(as.matrix(y.data)) 
  rdim <- ncol(as.matrix(input))
  #
  Gam <- matrix(0, nrow = qdim, ncol = rdim) 
  Gam[1,1]<-p1; Gam[1,2]<<-p2; Gam[1,3]<<-p3
  #
  cR <- matrix(0, pdim, qdim)
  cR[1,1] <- sigv
  cQ <- matQ(pdim, sigw1=sigw1, sigw2=sigw2, epsilonS=epsilonS.vector, 
             SeasPeriodS=seasonal.periods)
  mu0 <- matrix(20, pdim, ncol = 1)
  Sigma0 <- diag(80, pdim) 
  alpha <- .3
  kf <- KFSCov(num, y.data, alpha, A, mu0, Sigma0, Phi, Ups=0, Gam, cQ, cR, input)
  
  if (checkphi(Phi, small.phi, sigv, sigw1, sigw2, epsilonS.vector, tau)) {
    return(kf$like)
  } else {
    return(Inf)
  }
}
##
##