
#' Fit SLATE models over a tuning-parameter grid
#'
#' Fits \code{\link{slate}} models over a two-dimensional grid of elastic-net
#' mixing parameters and penalty strengths.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Design matrix for the multivariate feature matrix. This is transpose of
#'  Z matrix. The Z matrix is easier to conceptualize as a design matrix with
#'  rows representing the different features of matrix Y and columns representing
#'  design elements.
#' @param surv_covar Numeric vector or matrix of survival covariates. If a vector
#'   is supplied, it is converted to a one-column matrix. If surv_covar = NULL,
#'   it is assumed that there is no survival covariate and creates a vector of 0's,
#'   which results in cox regression coefficient forced to be 0.
#' @param surv_time Numeric vector of follow-up or event times.
#' @param surv_event Numeric  vector indicating event status, where
#'   `1` denotes an observed event and `0` denotes censoring.
#' @param alpha_grid Numeric vector of elastic-net mixing parameters.
#' @param lambda_frac Numeric vector of fractions of the maximum penalty value.
#' @param tol Numeric convergence tolerance passed to \code{\link{slate}}.
#'   Default is `1e-5`.
#' @param eps Numeric threshold for treating values as zero. Default is `1e-6`.
#' @param maxiter Integer maximum number of iterations passed to
#'   \code{\link{slate}}. Default is `20`.
#' @param gauss.nnodes Integer number of Gauss-Hermite quadrature nodes passed to
#'   \code{\link{slate}}. Default is `9`.
#' @param verbose Logical; if `TRUE`, prints progress from individual
#'   \code{\link{slate}} fits. Default is `FALSE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{grid}{A matrix of fitted \code{"slate"} objects, with rows
#'     corresponding to `alpha_grid` and columns corresponding to `lambda_frac`.}
#'   \item{alpha}{The input vector `alpha_grid`.}
#'   \item{lambda_frac}{The input vector `lambda_frac`.}
#' }
#'
#' @details
#' The maximum penalty value is computed using \code{compute_max_XYZ}. For each
#' value of `alpha_grid`, the penalty value is calculated as
#' \deqn{
#'   \lambda = \texttt{lambda\_frac} \times
#'   \frac{\texttt{max\_XYZ}}{\alpha}.
#' }
#'
#' Initial values are computed once using \code{\link{slate_init}} and reused for
#' each grid fit. Each fitted model is obtained by calling \code{\link{slate}}.
#'
#' @note
#' Values of `alpha_grid` should be strictly positive, since the maximum penalty
#' is computed as `max_XYZ / alpha`.
#'
#' @seealso
#' \code{\link{slate}}, \code{\link{slate_init}}
#'
#' @examples
#' set.seed(3302)
#' sim_slate <- slate_sim_data(nsample=50, feature_group=rep(1:5, rep(10,5)))
#' alpha_grid <- c((2:1)/2)
#' nlambda <- 2
#' ## in example, limit iterations and alpha/lambda for time constraints
#' ## for actual usage, set to default of 50
#' ## and 4 or more alpha/lambda values each
#' set.iter <- 3 
#' lambda_frac <- exp(seq(log(1), log(.1), length.out = nlambda))
#' fit_grid <- slate_grid(Y=sim_slate$Y, Zt=sim_slate$Zt, surv_covar=sim_slate$surv_covar,
#'                       surv_time=sim_slate$surv_time,
#'                       surv_event=sim_slate$surv_event,
#'                       alpha_grid, lambda_frac, maxiter=set.iter)
#' str(fit_grid)
#' @export
slate_grid <- function(Y, Zt, surv_covar, surv_time, surv_event,
                       alpha_grid, lambda_frac, tol=1e-5, eps=1e-6, maxiter=20,
                       gauss.nnodes=9, verbose=FALSE){

  nlatent <- 1
  
  ## current methods use 1/n  for lnlike
  ## and update for num and denom of W
  ## may need to damp updates to coeff updates, use CCD for complete updates to W
  
  if(is.null(surv_covar)){
    surv_covar <- rep(0, length(surv_time))
  }
  if(is.vector(surv_covar)){
    surv_covar <- matrix(surv_covar, ncol=1)
  }

  nlambda <- length(lambda_frac)
  nalpha <- length(alpha_grid)

  ## note 2-dimensional grid of fits. Each fit is a list
  ## and all fits arranged as matrix (of lists)
  grid <- matrix(list(), nrow=nalpha, ncol=nlambda)


  max_XYZ <- compute_max_XYZ(Y, Zt, nlatent=1)
  if(verbose) {
     cat("Computing initial values\n")
  }
    
  init <- slate_init(Y, Zt, surv_time, surv_event, surv_covar, iter_max=500, tol=1e-6, verbose=verbose)


  for(i in 1:nalpha){
    alpha <- alpha_grid[i]
    max_lambda <- max_XYZ/alpha

    for(j in 1:nlambda){

      lambda <- lambda_frac[j]*max_lambda
      if(verbose) {
          cat("================ alpha = ", alpha, ", lambda.frac = ", lambda_frac[j],
              ", lambda = ", lambda, " ==================\n")
      }
      if(FALSE){
        if(i==1 & j==1){
          init <- NULL
        } else if (i > 1 & j == 1){
          init <- grid[[i-1,j]]
        } else {
          init <- grid[[i,j-1]]
        }
      }

      grid[[i,j]] <- slate(Y, Zt, surv_covar, surv_time, surv_event,
                           alpha=alpha, lambda=lambda, tol=tol, eps=eps,maxiter=maxiter,
                           gauss.nnodes=gauss.nnodes, init=init, verbose=verbose)

    }

  }


  return(list(grid=grid, alpha = alpha_grid, lambda_frac=lambda_frac))
}



#' Fit penalized latent factor models over a tuning grid
#'
#' Fits a collection of penalized latent factor models over a two-dimensional
#' grid of elastic-net mixing parameters and penalty strengths.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric matrix of feature annotations or design covariates,
#'   corresponding to \eqn{Z^\top} in the latent factor model.
#' @param alpha_grid Numeric vector of elastic-net mixing parameters.
#' @param lambda_frac Numeric vector of fractions of the maximum penalty value.
#' @param tol Numeric convergence tolerance passed to
#'   \code{\link{mlm_factor_ccd}}. Default is `1e-5`.
#' @param eps Numeric threshold for treating values as effectively zero.
#'   Included for interface consistency. Default is `1e-6`.
#' @param maxiter Integer maximum number of iterations passed to
#'   \code{\link{mlm_factor_ccd}}. Default is `20`.
#' @param verbose Logical; if `TRUE`, prints progress information from model
#'   fitting routines. Default is `FALSE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{grid}{A matrix of fitted latent factor model objects returned by
#'     \code{\link{mlm_factor_ccd}}. Rows correspond to `alpha_grid` and columns
#'     correspond to `lambda_frac`.}
#'   \item{alpha}{The input vector `alpha_grid`.}
#'   \item{lambda_frac}{The input vector `lambda_frac`.}
#' }
#'
#' @details
#' The function first computes the maximum penalty scaling value using
#' \code{\link{compute_max_XYZ}} and obtains initial parameter estimates using
#' \code{\link{mlm_factor_em}}.
#'
#' For each value of `alpha_grid`, the maximum penalty parameter is computed as:
#' \deqn{
#'   \lambda_{\max} =
#'   \frac{\texttt{max\_XYZ}}{\alpha}.
#' }
#'
#' The actual penalty used for each fit is:
#' \deqn{
#'   \lambda =
#'   \texttt{lambda\_frac} \times \lambda_{\max}.
#' }
#'
#' Each model is then fit using
#' \code{\link{mlm_factor_ccd}} with the corresponding `alpha` and `lambda`
#' values.
#'
#' @note
#' Values in `alpha_grid` should be strictly positive, since the maximum penalty
#' is computed as `max_XYZ / alpha`.
#'
#' @seealso
#' \code{\link{mlm_factor_ccd}}, \code{\link{compute_max_XYZ}},
#' \code{\link{mlm_factor_em}}
#'
#' @export
fit_grid_mlm <- function(Y, Zt,
                     alpha_grid, lambda_frac, tol=1e-5, eps=1e-6, maxiter=20, verbose=FALSE){



  nlatent  <- 1
  
  nlambda <- length(lambda_frac)
  nalpha <- length(alpha_grid)

  ## note 2-dimensional grid of fits. Each fit is a list
  ## and all fits arranged as matrix (of lists)
  grid <- matrix(list(), nrow=nalpha, ncol=nlambda)

  max_XYZ <- compute_max_XYZ(Y, Zt, nlatent=1)
  if(verbose) {   
      cat("Computing initial values\n")
  }
  init <- mlm_factor_em(Y, Zt, nlatent)


  for(i in 1:nalpha){
    alpha <- alpha_grid[i]
    max_lambda <- max_XYZ/alpha

    for(j in 1:nlambda){

      lambda <- lambda_frac[j]*max_lambda


      cat("================ alpha = ", alpha, ", lambda.frac = ", lambda_frac[j], ", lambda = ", lambda, " ==================\n")

      grid[[i,j]] <-  mlm_factor_ccd(Y, Zt, nlatent, init=init, alpha, lambda,
                                     maxiter = 500, tol = 1e-5, verbose = FALSE)
    }

  }


  return(list(grid=grid, alpha = alpha_grid, lambda_frac=lambda_frac))
}


#' Select the Best Model from a SLATE Fit Grid
#'
#' Selects the model with the smallest Bayesian Information Criterion (BIC)
#' from a grid of fitted SLATE models. If multiple models have the same minimum
#' BIC, the model corresponding to the largest lambda penalty (i.e., the
#' largest column index) is returned.
#'
#' @param fit_grid A fitted SLATE grid object containing a component
#'   \code{grid}, where each element is a fitted model with a numeric
#'   \code{bic} value.
#'
#' @return A fitted SLATE model corresponding to the minimum BIC in the grid.
#'
#' @details
#' The function extracts the BIC from each model in the grid, identifies the
#' minimum BIC while ignoring missing values, and returns the corresponding
#' fitted model. When multiple models share the same minimum BIC, the tie is
#' broken by selecting the model with the largest lambda penalty, assumed to
#' correspond to the largest column index of the grid.
#'
#' @examples
#' \dontrun{
#' fit_grid <- slate_fit_grid(...)
#' best_fit <- slate_best_fit(fit_grid)
#' }
#'
#' @export
slate_best_fit <- function(fit_grid){
  nr <- nrow(fit_grid$grid)
  nc <- ncol(fit_grid$grid)
  bic <- matrix(0,nr,nc)  
  
  for(i in 1:nr){
    for(j in 1:nc){
      bic[i,j] <- fit_grid$grid[[i,j]]$bic
    }
  }
  
  
  min_val <- min(bic, na.rm = TRUE)
  index <- which(!is.na(bic) & bic == min_val, arr.ind = TRUE)
  ## if ties,choose index that has largest lambda penalty
  index_row  <- which(index[,2] ==  max(index[,2]), arr.ind = TRUE)[1]
  index <- index[index_row,]
  r <- index[1]
  c <- index[2]
  
  fit_best <- fit_grid$grid[[r,c]]
  return(fit_best)
}


