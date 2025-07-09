########## FUNCTIONS ##############

#' @title Correlated SBMs model parameter generator
#' @description Generates correlated stochastic block models probability matrices
#' @param A adjacency matrix corresponding to previous temporal snapshot
#' @param B block model probabilities
#' @param C community labels
#' @param rho correlation parameter
#' @return Probability matrix P
#' @export
generatePsbm <- function(A, B, C, rho){
  
  if(nrow(A)!=ncol(A)){
    stop("A must be a square matrix")
  }
  
  if(nrow(B)!=ncol(B)){
    stop("B must be a square matrix")
  }
  
  if(length(unique(C))!=ncol(B)){
    stop("Length of C and dim of B must match.")
  }
  
  if(rho<0 || rho>1){
    stop("Rho must be between 0 and 1.")
  }
  
  n = dim(A)[1]
  K = length(unique(C))
  C <- as.factor(C)
  
  # Membership matrix
  Z <- matrix(model.matrix(~C-1), nrow=n, ncol=K)
  # Probabilities before incorporating correlation
  P <- Z%*%B%*%t(Z)
  
  (P + rho*(1-P))^A * (P*(1-rho))^(1-A)
}


#' @title Correlated SBMs network generator
#' @description Generates correlated stochastic block models adjacency matrices
#' @param n number of nodes
#' @param B block model probabilities
#' @param C community labels
#' @param rho correlation parameter
#' @param TT number of time steps
#' @return Array of TT adjacency matrices
#' @export
generate_SBM <- function(n, B, C, rho, TT){
  
  if(nrow(B)!=ncol(B)){
    stop("B must be a square matrix")
  }
  
  if(length(unique(C))!=ncol(B)){
    stop("Length of C and dim of B must match.")
  }
  
  if(rho<0 || rho>1){
    stop("Rho must be between 0 and 1.")
  }
  
  A <- array(0, c(n, n, TT))
  
  K = length(unique(C))
  C <- as.factor(C)
  
  # Membership matrix
  Z <- matrix(model.matrix(~C-1), nrow=n, ncol=K)
  
  P <- Z%*%B%*%t(Z)
  At <- matrix(rbinom(n*n, 1, P), n, n)
  At[lower.tri(At)] <- 0
  At = At + t(At)
  diag(At) <- 0
  A[,,1] <- At
  
  for(tt in 2:TT){
    P <- generatePsbm(A[,,tt-1], B, C, rho)
    At <- matrix(rbinom(n*n, 1, P),n,n)
    At[lower.tri(At)] <- 0
    At = At + t(At)
    diag(At) <- 0
    A[,,tt] <- At
  }
  
  return(A)
}

#' @title Erdos-Renyi network generator
#' @description Generates adjacency matrix from Erdos-Renyi model
#' @param n number of nodes
#' @param p edge probability
#' @return Adjacency matrix A
#' @export
generateER <- function(n, p){
  A <- matrix(0, n, n)
  A <- matrix(rbinom(n*n, 1, p), n, n)
  A[lower.tri(A)] <- 0
  A = A + t(A)
  
  return(A)
}

#' @title Community detection p-value
#' @description Computes p-value for testing against ER null using Bickel and Sarkar (2016) method
#' @param A adjacency matrix
#' @return p-value for each temporal snapshot
#' @export
spectral.pval <- function(A){
  
  if(nrow(A)!=ncol(A)){
    stop("A must be a square matrix")
  }
  
  n=dim(A)[1]
  TT = dim(A)[3]
  
  pval = numeric(TT)
  
  for(t in 1:TT){
    p.hat <- sum(A[,,t])/(n*(n-1))
    
    P.hat <- p.hat - p.hat*diag(1,n)
    A.prime <- (A[,,t]-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    
    princ.eigen <- RSpectra::eigs_sym(A.prime,1,which="LA")[[1]]
    
    obs.stat <- n^(2/3)*(princ.eigen-2)
    pval[t] <- RMTstat::ptw(obs.stat, beta=1, lower.tail = FALSE)
  }
  
  return(pval)
}

#' @title Community detection p-value (with bootstrap)
#' @description Computes p-value for testing against ER null using Bickel and Sarkar (2016) bootstrap method
#' @param A adjacency matrix
#' @return p-value for each temporal snapshot
#' @export
spectral.adj.pval <- function(A){
  
  if(nrow(A)!=ncol(A)){
    stop("A must be a square matrix")
  }
  
  n=dim(A)[1]
  TT = dim(A)[3]
  
  apply_fun <- function(t){
    
    p.hat <- sum(A[,,t])/(n*(n-1))
    
    P.hat <- p.hat - p.hat*diag(1,n)
    A.prime <- (A[,,t]-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    
    princ.eigen <- RSpectra::eigs_sym(A.prime,1,which="LA")[[1]]
    
    obs.stat <- n^(2/3)*(princ.eigen-2)
    
    mu.tw <- -1.2065335745820 #from wikipedia
    sigma.tw <- sqrt(1.607781034581) #from wikipedia
    
    emp.stats <- numeric(50)
    
    for(i in 1:50){
      A.i <- generateER(n, p.hat)
      A.i.prime <- (A.i-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
      princ.eigen.i <- RSpectra::eigs_sym(A.i.prime,1,which="LA")[[1]]
      emp.stats[i] <- n^(2/3)*(princ.eigen.i-2)
    }
    
    mu.theta <- mean(emp.stats)
    sigma.theta <- sqrt(var(emp.stats))
    
    theta.prime <- mu.tw + (obs.stat-mu.theta)/sigma.theta * sigma.tw
    pval = RMTstat::ptw(theta.prime, beta=1, lower.tail = FALSE)
    
    return(pval)
  }
  
  pval = unlist(parallel::mclapply(1:TT, apply_fun, mc.cores = parallel::detectCores()-1))
  
  return(pval)
}


## Convert p-value to e-value
#' @title p-value to e-value calibrator
#' @description Converts a p-value to an e-value using a calibrator
#' @param p p-values
#' @param method calibration method: int (=integrated); VS (=Vovk-Sellk); kappa (=kappa)
#' @param kappa value to be used if method="kappa" selected
#' @return e-value 
#' @export
p_to_e <- function(p, method=c("int, VS, kappa"), kappa=NULL){
  
  if(method=="int"){
    if(p==1){
      eval = 0.5
    }else{
      eval = (1-p+p*log(p)) / (p*(-log(p))^2)
    }
    
  }else if(method=="VS"){
    if(p > exp(-1)){
      eval = 1
    }else{
      eval = -exp(-1)/(p*log(p))
    }
    
  }else if(method=="kappa"){
    
    if(kappa <= 0 || kappa >= 1){
      stop("kappa must be between 0 and 1 (exclusive).")
    }
    eval = kappa * p^(kappa-1)
    
  }else{
    stop("Calibrator type must be 'int', 'VS', or 'kappa'.")
  }
  
  return(eval)
}

# Convert multiple p-values at once
p_to_e <- Vectorize(p_to_e)