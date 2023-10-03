llk <- function(Q, T) {
  valQ <- 0
  signQ <- 0
  valQ <- sum(log(eigen(Q)$val))
  
  trQT <- sum(diag(Q %*% T))
  logL <- 0.5 * (valQ - trQT)
  
  return(logL)
}

Wvert <- function(Wi) {
  nWi <- length(Wi)
  V <- NULL
  
  for (j in 1:nWi) {
    Vcat <- Wi[[j]]
    V <- c(V, Vcat)
  }
  
  return(unique(V))
}

WvertC <- function(Wi, C) {
  nWi <- length(Wi)
  V <- NULL
  
  for (j in 1:nWi) {
    if (j != C) {
      Vcat <- Wi[[j]]
      V <- c(V, Vcat)
    }
  }
  
  return(unique(V))
}

setdifference <- function(Ea, Ec) {
  na <- length(Ea)
  
  for (j in na:1) {
    matchind <- any(Ec == Ea[j])
    if (matchind) {
      Ea <- Ea[-j]
    }
  }
  
  return(Ea)
}

setdifferencev <- function(Ea, Ec) {
  na <- length(Ea)
  
  for (j in na:1) {
    matchind <- any(Ec == Ea[j])
    if (matchind) {
      Ea <- Ea[-j]
    }
  }
  
  return(Ea)
}

sequencen <- function(n) {
  seqn <- 1:n
  
  return(seqn)
}

clique_part_Cpp <- function(W) {
  nW <- length(W)
  Wi <- NULL
  U <- vector("list", nW)
  nWi <- 0
  Ui <- NULL
  
  for (i in 1:nW) {
    Wi <- W[[i]]
    nWi <- length(Wi)
    
    for (j in 1:nWi) {
      Wij <- Wi[[j]]
      Ui <- c(Ui, Wij)
    }
    
    U[[i]] <- unique(Ui)
    Ui <- NULL
  }
  
  return(U)
}

clique_comppart <- function(U, n) {
  nU <- length(U)
  V <- sequencen(n)
  Ubar <- vector("list", nU)
  Ubari <- NULL
  
  for (i in 1:nU) {
    Ui <- U[[i]]
    Ubari <- setdifferencev(V, Ui)
    Ubar[[i]] <- Ubari
  }
  
  return(Ubar)
}

IPS2 <- function(S, n, T, W, maxit, eps) {

U <- clique_part_Cpp(W)
Ubar <- clique_comppart(U, n)

itcount <- 0
prevlogL <- 0
logL <- 0
converged <- FALSE
logLout <- rep(0, maxit)
V <- NULL
nW <- length(W)
Wi <- NULL
nWi <- 0
nEc <- 0
QM <- NULL
QMbar <- NULL
QMMbar <- NULL
Psi <- NULL
LambdaTemp <- NULL
Sc <- NULL
Scinv <- NULL
Qca <- NULL
Qaa <- NULL
Qtemp <- NULL
QT <- NULL
trQT <- NULL

Lambda <- matrix(0, n, n)
Q <- diag(n)

itcount <- 1
while (itcount < maxit) {
  for (i in 1:nW) {
    Wi <- W[[i]]
    V <- Wvert(Wi)
    nWi <- length(Wi)
    
    Uvert <- U[[i]]
    Ubarvert <- Ubar[[i]]
    
    QM <- Q[Uvert, Uvert]
    QMbar <- Q[Ubarvert, Ubarvert]
    QMMbar <- Q[Uvert, Ubarvert]
    if(length(V)==1){
      Psi <- t(QMMbar) %*% solve(QMbar) %*% (QMMbar)
    }else {
      Psi <- QMMbar %*% solve(QMbar) %*% t(QMMbar)
    }
    LambdaTemp <- QM - Psi
    Lambda[V, V] <- LambdaTemp
    
    for (j in 1:nWi) {
      Ec <- Wi[[j]]
      Ea <- WvertC(Wi, j)
      Ea <- setdifference(Ea, Ec)
      
      nEc <- length(Ec)
      Idn <- diag(nEc)
      Onesn <- matrix(1, nEc, nEc)
      
      Sc <- S[Ec, Ec]
      Scinv <- solve(Sc)
      
      Qca <- Lambda[Ec, Ea]
      Qaa <- Lambda[Ea, Ea]
      if(length(Ea)!=0){
      Lambda[Ec, Ec] <- Scinv + Lambda[Ec, Ea] %*% solve(Lambda[Ea, Ea]) %*% Lambda[Ea, Ec]
      } else{
        Lambda[Ec, Ec] <- Scinv
      }
    }
    
    Qtemp <- Lambda[V, V] + Psi
    Q[V, V] <- (Qtemp + t(Qtemp))/2
  }
  
  logL <- llk(Q, T)
  logLout[itcount] <- logL
  itcount <- itcount + 1
  
  if (itcount > 1) {
    if (abs(logL - prevlogL) < eps) {
      converged <- TRUE
      break
    } else {
      if (itcount == maxit) {
        converged <- FALSE
        break
      }
    }
  }
  
  prevlogL <- logL
}

logLout <- logLout[1:itcount-1]

result <- list(
  Q = Q,
  count = itcount,
  converged = converged,
  logL = logLout
)

return(result)
}

clique_part <- function(A){
  rownames(A) <- colnames(A) <- 1:ncol(A)
  A_NELO <- as(A,"matrix")
  Clist  <- getCliques(A_NELO)
  Clist  <- lapply(Clist, function(x) sort(as.numeric(x)))
  
  V      <- sum(lower.tri(A)*A)
  NN     <- length(Clist)
  Nc     <- NN/(1:NN) 
  NcTemp <- Nc[which(Nc %% 1 == 0)]
  if(length(NcTemp)>2){
    M      <- NcTemp[-c(1, length(NcTemp))]
  } else{
    M <- NcTemp
  }
  
  Mcomp  <- sapply(M, clique_select,Clist=Clist, V=V, NN=NN)
  
  nc     <- M[which(Mcomp==min(Mcomp))]
  nC     <- length(Clist)/nc
  W      <- lapply(1:nc, function(i) Clist[(i*nC-nC+1):(i*nC)])
  
  #Need to adjust index since C++ index begins with zero
  #W1  <- lapply(1:nc, function(i) lapply(W[[i]], function(x) x-1))  
  W
}


clique_select <- function(m, Clist, V, NN){
  nC  <- length(Clist)/m
  Wt  <- lapply(1:m, function(i) Clist[(i*nC-nC+1):(i*nC)])
  Umstar <- max(sapply(1:m, function(j) sum(sapply(Wt[[j]], length))))
  log( m*V^3 + NN * Umstar^3 )
}

