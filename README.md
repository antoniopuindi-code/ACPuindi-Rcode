# ACPuindi-Rcode

################################################################################
require(astsa)
require(forecast)
require(nlme)
library(data.table)
# CODIGO DAS FUNÇÕES DO MODELO "SCov"
#===============================================================================
matrizPhi <- function(smallPhi=NULL, seasonalPeriods=NULL, tau_s) {
  # Matriz Phi
  #1. Linha de nivel
  Phi <- matrix(c(1,smallPhi), nrow=1, ncol=2)
  if(!is.null(seasonalPeriods)) {
    #tau <- sum(seasonalPeriods)
    zero.tau <- matrix(0, nrow = 1, ncol = tau_s)
    Phi <- cbind(Phi, zero.tau)
  }
  # 2. Linha da tendencia
  beta.row <- matrix(c(0, smallPhi), nrow=1, ncol=2)
  if(!is.null(seasonalPeriods)) {
    beta.row <- cbind(beta.row, zero.tau)
  }
  Phi <- rbind(Phi, beta.row)
  
  #3. Linha da Sasonalidade
  if (!is.null(seasonalPeriods)) {
    seasonal.row <- t(zero.tau)
    if (!is.null(beta)) {
      seasonal.row <- cbind(seasonal.row, seasonal.row)
    }
    
    # Cria a matriz A
    for (i in seasonalPeriods) {
      if (i == seasonalPeriods[1]) {
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
#phi = .4
#sp = c(1,2,3)
#Phi = matrizPhi(smallPhi=phi, seasonalPeriods=sp)
#===============================================================================
cria_A = function(smallPhi_s=NULL, sPeriods_s=NULL){
  
  adjustPhi = 0
  numCols = 1 
  numSeasonal = 0
  lengthSeasonal = 0
  smallPhi = 0
  seasonalPeriods = 0
  
  if(!is.null(smallPhi_s)) {
    smallPhi = as.numeric(smallPhi_s)
    adjustPhi = 1
    numCols = numCols + 1
  }
  if(!is.null(sPeriods_s)) {
    seasonalPeriods = as.integer(sPeriods_s)
    numSeasonal = length(sPeriods_s)
    for(s in 1:numSeasonal) {
      lengthSeasonal = lengthSeasonal + seasonalPeriods[s]
    }
    numCols = numCols + lengthSeasonal
  } else {
    lengthSeasonal = 0
  }
  
  mA = matrix(0, nrow = 1, ncol = numCols)
  
  if(!is.null(sPeriods_s)) {
    position = adjustPhi
    for(s in 1:numSeasonal) {
      position = position + seasonalPeriods[s]
      mA[1,position] = 1
    }
  }
  
  mA[1,1] = 1
  
  if(adjustPhi == 1) {
    mA[1,2] = smallPhi
  }
  
  return(mA)
}
#phi = .4
#sP = c(1,2,3)
#criaA(smallPhi_s=phi, sPeriods_s=sP)

#===============================================================================
insert.diag <- function(matQ, varianceSeas, start=c(1,1), dir=c(1,1)) {
  
  sq <- seq_along(varianceSeas)-1
  indices <- sapply(1:2, function(i) start[i] + dir[i]*sq)
  stopifnot(all(indices>0))
  stopifnot(all(indices[,1]<=nrow(matQ)))
  stopifnot(all(indices[,2]<=ncol(matQ)))  
  
  matQ[indices] <- varianceSeas
  matQ
}
#-------------------------------------------------------------------------------
matQ = function(pdim, sigw1, sigw2, gammaVector_s=NULL, seasonalPeriods_s=NULL){
  
  gammaLength = 0
  adjustSigw2 = 0
  numCols = 1
  if(!is.null(sigw2)) {
    numCols = numCols + 1
    adjustSigw2 = 1
  }
  # Find the length of the gamma/seasonal bit
  if((!is.null(gammaVector_s))&&(!is.null(seasonalPeriods_s))) {
    gammaVector = as.numeric(gammaVector_s)
    seasonalPeriods = as.integer(seasonalPeriods_s)
    for(i in 1:length(seasonalPeriods_s)) {
      gammaLength = gammaLength + seasonalPeriods[i]
    }
    numCols = numCols + gammaLength
  } else {
    gammaLength = 0
  }
  
  gTranspose = matrix(0, nrow = 1, ncol = numCols)
  gTranspose[1,1] = as.numeric(sigw1)
  if(!is.null(sigw2)) {
    gTranspose[1,2] = as.numeric(sigw2)
  }
  
  # Copy the gamma/seasonal bits
  if((!is.null(gammaVector_s))&&(!is.null(seasonalPeriods_s))) {
    position = adjustSigw2 + 1
    gTranspose[1, position] = gammaVector
    if(length(gammaVector_s) > 1) {
      for(s in 1:(length(seasonalPeriods_s)-1)) {
        position = position + seasonalPeriods[s]
        gTranspose[1, position] = gammaVector[(s+1)]
      }
    }
  }
  #seasonalPeriods = 0
  #gammaVector = 0
  matQ <- diag(as.numeric(gTranspose), pdim)
  return(matQ)
}


#===============================================================================
unParameteriseSMWC<-function(param.vector, control) {
  
  p1<-param.vector[1]
  p2<-param.vector[2]
  p3<-param.vector[3]
  sigv<-param.vector[4]
  sigw1<-param.vector[5]
  sigw2<-param.vector[6]
  if(control$use.damping) {
    small.phi<-param.vector[7]
    gamma.start<-8
  } else {
    small.phi<-1
    gamma.start<-7
  }
  
  if(control$length.gamma > 0) {
    gamma.vector<-param.vector[gamma.start:(gamma.start+control$length.gamma-1)]
    final.gamma.pos<-gamma.start+control$length.gamma-1
  } else {
    gamma.vector=NULL
    final.gamma.pos<-gamma.start-1
  }
  
  return(list(p1=p1, p2=p2, p3=p3, sigv=sigv, sigw1=sigw1, sigw2=sigw2, 
              small.phi=small.phi, gamma.v=gamma.vector))
}

#===============================================================================
parameteriseSMWC <- function(p1, p2, p3, sigv, sigw1, sigw2, small.phi=1, 
                             gamma.v=NULL) {
  
  param.vector<-cbind(p1, p2, p3, sigv, sigw1, sigw2)
  if (is.null(small.phi)) {
    use.damping <- FALSE
  } else if (small.phi != 1) {
    param.vector <- cbind(param.vector, small.phi)
    use.damping <- TRUE
  } else {
    use.damping <- FALSE
  }
  if (!is.null(gamma.v)) {
    gamma.v <- matrix(gamma.v, nrow = 1, ncol = length(gamma.v))
    param.vector <- cbind(param.vector, gamma.v)
    length.gamma <- length(gamma.v)
  } else {
    length.gamma <- 0
  }
  
  control <- list(use.damping = use.damping, length.gamma = length.gamma)
  return(list(vect = as.numeric(param.vector), control = control))
}

#===============================================================================
ParscaleSMWC<-function(control) {
  
  parscale<-c(rep(1e-2,3),rep(1e-2,3))
  
  if(control$use.damping) {
    parscale<-c(parscale, 1e-2)
  } 
  if(control$length.gamma > 0) {
    parscale<-c(parscale, rep(1e-2, control$length.gamma))
  }
  
  return(parscale)
}

#===============================================================================
Kf_setric<-function(num, y, alpha, A, mu0, Sigma0, Phi, Ups, Gam, cQ, cR, input){
  #browser()
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
  
  #xp = array(NA, dim = c(pdim, 1, (num+h)))
  #Pp = array(NA, dim = c(pdim, pdim, (num+h)))
  #f = array(NA, dim = c(qdim, 1, (num+h)))
  #rmse = rep(0,h)
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
    # Conectar o ajuste com a previsão multi-passos
    #for (t in i:num+h) {
    #  xp[, , t] = Phi %*% xf[, , i - 1] + Ups %*% ut[i, ]
    #  Pp[, , t] = Phi %*% Pf[, , i - 1] %*% t(Phi) + Q 
    #  f[, , t] = B %*% xp[, , t] + Gam %*% ut[i, ] 
    #  sigf = B %*% Pp[, , t] %*% t(B) + R
    #  rmse[t] = sqrt(sigf)
    #}
    # Fim da previsão
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
#===============================================================================
checkphi<-function(Phi, small.phi=NULL, sigv=NULL, sigw1=NULL, sigw2=NULL,
                   gamma.vector=NULL, tau=0) {
  
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
  if(!is.null(gamma.vector)){
    if(gamma.vector < 0)
      return(FALSE)
  }
  if(!is.null(small.phi)){
    if(((small.phi < .8) | (small.phi > 1)))
      return(FALSE)
  }
  Phi.eigen.values<-eigen(Phi, symmetric=FALSE, only.values=TRUE)$values
  return(all(abs(Phi.eigen.values) < 1 + 1e-2))
}

#===============================================================================
LinnSMWC <- function(param.vector, y.data, use.small.phi, seasonal.periods, tau=0) {
  
  p1 <- param.vector[1]
  p2 <- param.vector[2]
  p3 <- param.vector[3]
  sigv <- param.vector[4]
  sigw1 <- param.vector[5]
  sigw2 <- param.vector[6]
  if (use.small.phi) {
    small.phi <- param.vector[7]
    gamma.start <- 8
  } else {
    small.phi <- 1
    gamma.start <- 7
  }
  
  if (!is.null(seasonal.periods)) {
    gamma.vector <- param.vector[gamma.start:(gamma.start + length(seasonal.periods) - 1)]
    final.gamma.pos <- gamma.start + length(gamma.vector) - 1
  } else {
    gamma.vector <- NULL
    final.gamma.pos <- gamma.start - 1
  }
  
  Phi <<- matrizPhi(smallPhi=small.phi, seasonalPeriods=seasonal.periods, 
                    tau_s=tau)
  
  A <<- array(cria_A(smallPhi_s=small.phi, sPeriods_s=seasonal.periods), 
              dim = c(1, nrow(Phi), num))
  
  ## Dimensões matriciais
  pdim <- nrow(Phi) 
  qdim <- ncol(as.matrix(y.data)) 
  rdim <- ncol(as.matrix(input))
  
  Gam <<- matrix(0, nrow = qdim, ncol = rdim) 
  Gam[1,1]<<-p1; Gam[1,2]<<-p2; Gam[1,3]<<-p3
  #Ups <<- matrix(0, nrow = pdim, ncol = rdim)
  #Ups[1,1]<<-p1; Ups[2,1]<<-p2
  
  ## Matrizes de variancias
  cR <<- matrix(0, pdim, qdim)
  cR[1,1] <<- sigv
  cQ <<- matQ(pdim, sigw1=sigw1, sigw2=sigw2, gammaVector_s=gamma.vector, 
              seasonalPeriods_s=seasonal.periods)
  #mu0 <<- matrix(0, pdim, ncol = 1)
  #Sigma0 <<- diag(99, pdim)
  mu0 <<- matrix(20, pdim, ncol = 1)
  Sigma0 <<- diag(80, pdim) 
  alpha <<- .99
  
  kf <<- Kf_setric(num, y.data, alpha, A, mu0, Sigma0, Phi, Ups=0, Gam, cQ, cR, 
                   input)
  
  if (checkphi(Phi, small.phi, sigv, sigw1, sigw2, gamma.vector, tau)) {
    return(kf$like)
  } else {
    return(Inf)
  }
}
#===============================================================================

filtroSMWC <- function(y, input, damping, seasonal.periods, force.seasonality=FALSE){
  # Codigo do filtro dos modelos
  if (!damping) {
    return(list(AIC = Inf))
  }
  
  first.model <- fitSMWC(y, input, use.damping = damping, seasonal.periods = seasonal.periods)
  if (!is.null(seasonal.periods) && !force.seasonality) {
    non.seasonal.model <- fitSMWC(y, input, use.damping = damping, seasonal.periods = NULL)
    if (first.model$AIC > non.seasonal.model$AIC) {
      seasonal.periods <- NULL
      first.model <- non.seasonal.model
    }
  }
  starting.params <- first.model$parameters
  return(first.model)
}

============================================================================================
###############################################################################
#===============================================================================
fitSMWC <- function(y, input, use.damping, seasonal.periods=NULL, 
                    starting.params=NULL) {
  
  if (!is.null(seasonal.periods)) {
    seasonal.periods <- as.integer(sort(seasonal.periods))
  }
  
  if (is.null(starting.params)) {
    # Calculate starting values:
    if (sum(seasonal.periods) > 16) {
      p1 = p2 = p3 <- 0.1
      sigv <- 0.1
      sigw1 <- 0.002
      sigw2 <- 0.3
    } else {
      #p1 = p2 = p3 <- 0.01
      #sigv <- 0.001
      #sigw1 <- 0.05
      #sigw2 <- 0.05
      #p1 = p2 = p3 <- 0.1
      #sigv <- 0.1
      #sigw1 <- 0.002
      #sigw2 <- 0.3
      p1 = p2 = p3 <- 1e-6
      sigv <- 1e-10
      sigw1 <- 1e-10
      sigw2 <- 1e-10
    }
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
    paramz <- unParameteriseSMWC(starting.params$vect, starting.params$control)
    p1 <- paramz$p1
    p2 <- paramz$p2
    p3 <- paramz$p3
    sigv <- paramz$sigv
    sigw1 <- paramz$sigw1
    sigw2 <- paramz$sigw2
    small.phi <- paramz$small.phi
    gamma.v <- paramz$gamma.v
  }
  
  ## Optimise the starting values:
  # Make the parameter vector  parameterise
  param.vector <- parameteriseSMWC(p1=p1, p2=p2, p3=p3, sigv=sigv, sigw1=sigw1,
                                   sigw2=sigw2, small.phi = small.phi, 
                                   gamma.v = gamma.v)
  par.scale <- ParscaleSMWC(param.vector$control)
  
  if (!is.null(seasonal.periods)) {
    tau <- sum(seasonal.periods)
  } else {
    tau <- 0
  }
  
  # Optimise the likelihood function
  if (length(param.vector$vect) > 1) {
    optim.like <<- optim(par = param.vector$vect, fn = LinnSMWC, y.data=y,
                         method = "Nelder-Mead", use.small.phi = use.damping, 
                         seasonal.periods = seasonal.periods, tau = tau, 
                         control = list(maxit=(100*length(param.vector$vect)^2), 
                                        parscale = par.scale))
  } else {
    optim.like <<- optim(par = param.vector$vect, fn = LinnSMWC, y.data=y, 
                         method = "BFGS", use.small.phi = use.damping, 
                         seasonal.periods = seasonal.periods, tau = tau, 
                         control=list(parscale=par.scale))
  }
  
  # Get the parameters out of the param.vector
  paramz <- unParameteriseSMWC(optim.like$par, param.vector$control)
  p1 <- paramz$p1
  p2 <- paramz$p2
  p3 <- paramz$p3
  sigv <- paramz$sigv
  sigw1 <- paramz$sigw1
  sigw2 <- paramz$sigw2
  small.phi <- paramz$small.phi
  gamma.v <- paramz$gamma.v
  
  # Previsão um passo a frente da observação e erro de previsão
  fitted.values <- kf$fprev[, , 1:450]
  e <- kf$innov
  #variance <- sum((e * e))/length(y)
  # Obter a verossimilhança e o AIC
  likelihood <- optim.like$value
  aic <- likelihood + 2*(length(param.vector$vect) + nrow(kf$xp))
  
  ## Estimativas e os respetivos erros
  #SE = sqrt(diag(solve(optim.like$hessian)))
  #parErros <- round(cbind(estimate=optim.like$par, SE), 6)
  parErros <- round(cbind(estimate=optim.like$par), 6)
  #colnames(parErros) <- c("parametros")
  #rownames(parErros) <- c("cReg1", "cReg2", "cReg3", "sigv", "sigw1", "sigw2", 
  #                        "phi", "sigSeas")
  parErros <- as.table(parErros)
  
  # Listar dos objectos
  fits <- ts(c(fitted.values))
  fits[1:4] <- y[1:4]
  e <- ts(c(e))
  tsp(fits) <- tsp(e) <- tsp(y)
  
  model.saida <- list(y=y,input=input, fitted.values=fits, errors=e,
                      parameters=list(vect=optim.like$par,control=param.vector$control),
                      p1=p1, p2=p2, p3=p3, sigv=sigv, sigw1=sigw1, sigw2=sigw2, 
                      damping.parameter=small.phi, sigSeas=gamma.v, parameters=parErros, 
                      likelihood=likelihood, val.optim=optim.like$convergence, AIC=aic,
                      seasonal.periods=seasonal.periods)
  class(model.saida) <- c("SMWC")
  return(model.saida)
}

===============================================================================================

# Modelo Estrutural Básico com Covariáveis
SCov <- function (y, input, use.damped.trend = NULL, seasonal.periods = NULL) {
    
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
        current.model <- filtroSMWC(y, input, damping = damping,
                                    seasonal.periods = seasonal.periods)
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
    #if (best.model$optim.return.code != 0) {
    #  warning("optim() did not converge.")
    #}
    attributes(best.model$fitted.values) <- attributes(best.model$errors) <- attributes(origy)
    best.model$y <- origy
    best.model$series <- seriesname
    best.model$method <- "EMWC"
    return(best.model)
}


