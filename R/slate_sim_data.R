#' Simulate Small Example Data for slate
#'
#' For n-samples and grouped features vector, simulate data for Y (nxp matrix of features),
#' a covariate for survival, and survival time and event status. A matrix Zt for the design
#' of the features, containing indicators for group membership of the features.
#' 
#' @param nsample the number of samples in the experiment
#'
#' @param feature_group index of groupings for genomics/multi-omics features. The full length of the fector is the total number of features to be simulated.
#' 
#' @return a list with items:
#' \describe{
#'   \item{Y}{Numeric matrix of observed features, with subjects in rows and
#'     features in columns.
#'   }
#'   \item{Zt}{ Design matrix for the multivariate feature matrix, the transpose of
#'     the Z matrix
#'   }
#'   \item{surv_covar}{Numeric vector or matrix of survival covariates}
#'   \item{surv_time}{Numeric vector of follow-up or event times}
#'   \item{surv_event}{Numeric vector indicating event status}
#' }
#' 
#' @export
slate_sim_data <- function(nsample=25, feature_group=rep(1:10,rep(10,10))) {
  
  moments_x_given_Y <- function(Y, W, psi, Zt){
    ## conditional moments of X given Y
    
    nlat <- nrow(W)
    imat <- diag(nlat)
    psi_inv = diag(1.0 / psi)
    WZt <- W %*% Zt
    
    inv_lowdim <- ginv(imat + WZt %*% psi_inv %*% t(WZt))
    vinv <-  psi_inv - psi_inv %*% t(WZt) %*% inv_lowdim %*% WZt %*% psi_inv
    
    Ex <-   Y %*% t(WZt %*% vinv)
    var_x_y <-  imat - WZt %*% vinv %*% t(WZt)
    dim(var_x_y) <- NULL
    
    Ex[is.na(Ex)] <- 0   
    
    return(list(mu = Ex, var=var_x_y ))
  }
  
  
  n <- nsample  ##  25
  p <- length(feature_group) ## 100
  Y <- matrix(rnorm(n*p), nrow=nsample)
  z_index <- feature_group #rep(1:10, rep(10,10))
  nz <- length(unique(z_index))  
  if(length(z_index) < 5 | nz <2 ) {
      stop("feature_group has either too few features or inadequate length\n")
  }
  Z <- matrix(0, nrow=p, ncol=nz)
  Z <- 1*(col(Z) == z_index)
  Zt <- t(Z)
  surv_covar <- sample(c(0,1), n,replace=TRUE)
  lambda <- 0.1        # baseline hazard
  
  beta_surv <- 1
  nlatent <- 1
  W <- matrix(0, nrow=nlatent, ncol=nrow(Zt))
  psi <- rep(1, ncol(Y))
  W[1,1] <-  1
  psi[1] <- .9
  mom <- moments_x_given_Y(Y, W, psi, Zt)
  x_sim <- mom$mu + sqrt(mom$var) * rnorm(n)
  
  beta_latent <- .5
  
  # Event times from exponential baseline hazard
  u <- runif(n)
  eta <- surv_covar * beta_surv  + x_sim*beta_latent
  T <- -log(u) / (lambda * exp(eta))
  C <- rexp(n, rate = 0.05)
  
  surv_time <- pmin(T, C)
  surv_event <- as.integer(T <= C)
  
  
  if(is.null(colnames(Y))){
    colnames(Y) <- paste0("Y", 1:ncol(Y))
  }
  
  if(is.null(rownames(Zt))){
    rownames(Zt) <- paste0("Zt", 1:nrow(Zt))
  }
  
  return(list(Y=Y, Zt=Zt, surv_covar=surv_covar, surv_time=surv_time, surv_event=surv_event))
}
