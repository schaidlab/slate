#' Sparse Latent factor model for Annotated features and Time-to-Event data (SLATE)
#'
#' Fits a survival latent factor model for multivariate features (e.g., omics) that can
#' be annotated and time-to-event
#' outcomes using penalized iterative updates. The current implementation supports
#' only one latent variable.
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
#' @param surv_event Numeric vector indicating event status, where
#'   `1`  denotes an observed event and `0` denotes censoring.
#' @param alpha Elastic-net mixing parameter for penalizing `W`. `alpha = 1`
#'   gives an L1 penalty, `alpha = 0` gives an L2 penalty.
#' @param lambda  nonnegative penalty tuning parameter.
#' @param tol Numeric convergence tolerance. Default is `1e-5`.
#' @param eps Numeric threshold used to determine whether coefficients or loadings
#'   are effectively zero. Default is `1e-6`.
#' @param maxiter Integer maximum number of outer iterations. Default is `50`.
#' @param gauss.nnodes Integer number of Gauss-Hermite quadrature nodes. Default
#'   is `9`.
#' @param init Optional list of initial values. When supplied, it should contain
#'   `W`, `psi`, `coef_surv`, `coef_latent`, `risk_sets`, `haz_subj`,
#'   `cumhaz_subj`, and `cox_event_time`.
#' @param verbose Logical; if `TRUE`, prints progress information during fitting.
#'   Default is `FALSE`.
#'
#' @return An object of class `"slate"`, returned as a list with components:
#' \describe{
#'   \item{W}{Estimated factor loading matrix.}
#'   \item{psi}{Estimated feature-specific residual variances.}
#'   \item{coef_surv}{Estimated coefficients for observed survival covariates.}
#'   \item{coef_latent}{Estimated coefficients for latent factors.}
#'   \item{var_coef}{Estimated covariance matrix for survival coefficients.}
#'   \item{df}{Estimated model degrees of freedom.}
#'   \item{bic}{Bayesian information criterion.}
#'   \item{zstat}{Wald z-statistics for survival coefficients.}
#'   \item{pval}{Two-sided Wald p-values for survival coefficients.}
#'   \item{risk_sets}{Risk-set information at event times.}
#'   \item{haz_subj}{Estimated subject-specific baseline hazard values.}
#'   \item{cumhaz_subj}{Estimated subject-specific cumulative baseline hazards.}
#'   \item{cox_event_time}{Unique event times from the Cox model fit.}
#' }
#'
#' @details
#' The function initializes latent factor parameters either from `init` or by
#' fitting a factor model followed by a Cox proportional hazards model. It then
#' iteratively updates the factor loadings, residual variances, Cox regression
#' coefficients, baseline hazard, and expected latent quantities using and EM
#' algorithm with Gauss-Hermite quadrature to integrate functions that depend
#' on latent variables.
#'
#' Convergence is assessed using the relative change in the penalized objective
#' function:
#' \deqn{
#'   -\ell / n + \alpha \lambda \sum |W| +
#'   (1 - \alpha)\lambda \sum W^2 / 2.
#' }
#'
#' @note
#' This implementation assumes 1 latent variable because numerical
#' integration is currently implemented only for a single latent variable.
#'
#' @seealso
#' \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}}
#'
#' @export
slate <- function(Y, Zt, surv_covar=NULL, surv_time, surv_event,
                            alpha, lambda=NULL, tol=1e-5, eps=1e-6, maxiter=50,
                            gauss.nnodes=9, init=NULL, verbose=FALSE){

  Call <- match.call()
  
  
  
  if(is.null(colnames(Y))){
    colnames(Y) <- paste0("Y", 1:ncol(Y))
  }
  
  if(is.null(rownames(Zt))){
    rownames(Zt) <- paste0("Zt", 1:nrow(Zt))
  }
  
  ## for now only allow 1 latent because numerical integration is setup for that
  ## would need multivariate integration for > 1 latent variable
  
  nlatent <- 1
 
  ## current methods use 1/n  for lnlike
  ## and update for num and denom of W
  ## may need to damp updates to coeff updates, use CCD for complete updates to W

  if(is.null(surv_covar)){
    surv_covar <- rep(0, length(surv_time))
    names(surv_covar) <- "no-covar"
  }

  if(is.vector(surv_covar)){
    surv_covar <- matrix(surv_covar, ncol=1)
  }

  if(is.null(colnames(surv_covar))){
    colnames(surv_covar) <- paste0("cov", 1:ncol(surv_covar))
  }
  
  nsubj <- nrow(Y)
  nfeat <- ncol(Y)
  ncovar <- ncol(surv_covar)
  ncoef <- ncovar + nlatent

  gh.dat <- gauss_hermite_quad(dim=nlatent, nnodes=gauss.nnodes)

  wt <- gh.dat$wt
  node <- gh.dat$node

  if(!is.matrix(node)){
    node <- as.matrix(node)
  }


  ## init values
  if(verbose) {
    cat("estimating initial values\n")
  }


  if(is.null(init)){

    ## try using pca for init W and psi to avoid em for full factor model
    ##fit_pca <- mlm_mle_pca(Y,Zt, nlatent)
    ##fit_factor <- list(W=fit_pca$W, psi=rep(fit_pca$sigma, nfeat))

    ## full factor em takes longer to get init value but reduces number of iterations below
    fit_factor <- mlm_factor_em(Y,Zt, nlatent=nlatent, iter_max = 500, tol=1e-6, verbose=FALSE)
    fit_Ex <- mlm_Ex(Y, Zt, fit_factor$W, fit_factor$psi)
    fit_cox <- coxph(Surv(surv_time, surv_event) ~ surv_covar + fit_Ex$Ex)

    cox_detail <- coxph.detail(fit_cox)
    cox_detail$cum_hazard <- cumsum(cox_detail$hazard)
    cox_event_time <- cox_detail$time

    risk_sets <-vector("list", length=length(cox_event_time))
    index_subj <- 1:length(surv_time)
    for(i in 1:length(cox_event_time)){
      risk_sets[[i]] <- list(subj= index_subj[surv_time >= cox_event_time[i]],
                             event= index_subj[surv_time == cox_event_time[i] & surv_event == 1])
    }


    haz_subj <- expand_hazard(cox_detail$hazard, surv_time, cox_event_time)
    cumhaz_subj <- expand_hazard(cox_detail$cum_hazard, surv_time, cox_event_time)

    coef <- fit_cox$coefficients
    coef[is.na(coef)] <- 0

    coef_surv <- coef[1:ncovar]
    coef_latent <- coef[(ncovar+1):(ncovar+nlatent)]
    W <- fit_factor$W
    psi <- fit_factor$psi
  } else {

    W <- init$W
    psi <- init$psi
    coef_surv <- init$coef_surv
    coef_latent <- init$coef_latent
    coef <- c(coef_surv, coef_latent)
    risk_sets <- init$risk_sets
    haz_subj <- init$haz_subj
    cumhaz_subj <- init$cumhaz_subj
    cox_event_time <- init$cox_event_time

    if(any(is.na(W))){
      cat("init W is NA\n")
      W[is.na(W)] <- 0
    }
    if(any(is.na(coef_latent))){
      cat("init coef is NA\n")
      coef_latent[is.na(coef_latent)] <- 0
    }
  }


  ## expected values of functions of latent variables
  ex.dat <- expected_func_x(Y, W, psi, Zt, surv_covar, coef_surv,coef_latent,
                            surv_event, cumhaz_subj, gh.dat)
  haz <- update_hazard(surv_covar, coef_surv, ex.dat$Exb, risk_sets)
  cumhaz <- cumsum(haz)
  haz_subj <- expand_hazard(haz, surv_time, cox_event_time)
  cumhaz_subj <- expand_hazard(cumhaz, surv_time,cox_event_time)


  lnlike <- mlm_factor_lnlike(Y, Zt, W, psi) +
    sum(surv_event * (log(haz_subj) + surv_covar %*% coef_surv + ex.dat$Ex %*%coef_latent) -
          ex.dat$Exb*exp((surv_covar %*% coef_surv))*cumhaz_subj)

  ##  update penalty function
  penfunc_curr <- -lnlike/(nsubj)  + alpha*lambda*sum(abs(W))  + (1-alpha)*lambda*sum(W^2)/2

  if(is.na(penfunc_curr)) penfunc_curr <- 0
  penfunc_old <- penfunc_curr


  ## begin iterative updates
  converge <- FALSE
  
  for(iter in 1:maxiter){

    ## get expected values
    Ex <- ex.dat$Ex
    Exx_sum <- apply(ex.dat$Exx,c(1, 2), sum)

    ## update W and psi
    fit_W  <- update_W(Y, Zt, W, psi, Ex, Exx_sum, alpha, lambda)
    W <- fit_W$W
    psi <- fit_W$psi

    W[is.na(W)] <- 0
    psi[is.na(psi)] <- 1

    ##  update cox regression coefficients

    coef <- update_cox_coef(coef, surv_covar, surv_event, cumhaz_subj, risk_sets,
                    Y, W, psi, Zt, ex.dat, gh.dat, tol)

    coef[is.na(coef)] <- 0

    coef_surv <- coef[1:ncovar]
    coef_latent <- coef[(ncovar+1):(ncovar + nlatent)]

    ##  update baseline hazard
    haz <- update_hazard(surv_covar, coef_surv, ex.dat$Exb, risk_sets)
    cumhaz <- cumsum(haz)
    haz_subj <- expand_hazard(haz, surv_time, cox_event_time)
    cumhaz_subj <- expand_hazard(cumhaz, surv_time, cox_event_time)

    ## update lnlike
    lnlike <- mlm_factor_lnlike(Y, Zt, W, psi) +
      sum(surv_event * (log(haz_subj) + surv_covar %*% coef_surv + ex.dat$Ex %*%coef_latent) -
      ex.dat$Exb*exp((surv_covar %*% coef_surv))*cumhaz_subj)

    ##  update penalty function
    penfunc_curr <- -lnlike/(nsubj)  + alpha*lambda*sum(abs(W))  + (1-alpha)*lambda*sum(W^2)/2
    if(is.na(penfunc_curr)) penfunc_curr <- 0

    ##  check convergence
    pen_rel_change <- abs(penfunc_curr - penfunc_old) / (abs(penfunc_old) + tol)


    if (verbose) {
      cat("iter =", iter, ", penfunc_curr =", penfunc_curr, ", W range =", range(W),
          ", pen_rel_change = ", pen_rel_change,  "\n")
    }

    if (pen_rel_change < tol) {
      converge <- TRUE
      break
    }

    ## update E[x]
    ex.dat <- expected_func_x(Y, W, psi, Zt, surv_covar, coef_surv,coef_latent,
                              surv_event, cumhaz_subj, gh.dat)
    penfunc_old <- penfunc_curr
  }

  coef <- c(coef_surv, coef_latent)
  coef[is.na(coef)] <- 0

  score_dat <- expected_cox_score(surv_covar, surv_event, coef, cumhaz_subj, risk_sets,
                                  Y, W, psi, Zt, ex.dat, gh.dat)
  var_score <- score_dat$var_score
  info <- score_dat$info

  info_louis <- score_dat$info_louis

  var_coef <- matrix(NA, nrow=ncoef, ncol=ncoef)

  W[is.na(W)] <- 0
  if(all(abs(W) < eps)){
    coef_latent <- 0 *coef_latent
    fit_cox <- coxph(Surv(surv_time, surv_event) ~ surv_covar)
    coef_surv <- fit_cox$coef
    var_coef[1:ncovar, 1:ncovar] <- fit_cox$var
    coef <- c(coef_surv, coef_latent)
  } else{
    var_coef <- ginv(info_louis)
  }

  df_model <- mlm_factor_df(W, psi) + sum(abs(coef) > eps)
  

  bic <- -2*lnlike + df_model*log(nsubj)
  ## if bic is NA, due to lnlike NA, make it Inf to
  if(is.na(bic)) bic <- Inf
  
  zstat <- coef/sqrt(diag(var_coef))
  pval <- 2*(1-pnorm(abs(zstat)))

  
  rownames(W) <- paste0("latent", 1:nlatent)
  colnames(W) <- rownames(Zt)
  covar_names <- c(colnames(surv_covar), paste0("latent", 1:nlatent))
  names(psi) <- colnames(Y)
  coef <- c(coef_surv, coef_latent)
  names(coef) <- covar_names
  rownames(var_coef) <- colnames(var_coef) <- covar_names
  names(zstat) <- covar_names
  names(pval) <- covar_names
  coef <- c(coef_surv, coef_latent)
  names(coef) <- covar_names
  
  fit <- list(W=W, psi=psi, coef_surv=coef_surv,coef_latent=coef_latent,
             var_coef = var_coef,  df=df_model, bic=bic, zstat=zstat, pval=pval,
             coef=coef,
             risk_sets = risk_sets,
             haz_subj = haz_subj,
             cumhaz_subj = cumhaz_subj, cox_event_time=cox_event_time,
             iter=iter, converge=converge)

  oldClass(fit) <- "slate"
  return(fit)
}


#' Initialize parameters for a SLATE model
#'
#' Computes initial values for \code{\link{slate}} by fitting a latent factor
#' model and then a Cox proportional hazards model using the estimated latent
#' factor expectations.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Design matrix for the multivariate feature matrix. This is transpose of
#'  Z matrix. The Z matrix is easier to conceptualize as a design matrix with
#'  rows representing the different features of matrix Y and columns representing
#'  design elements.
#' @param surv_time Numeric vector of follow-up or event times.
#' @param surv_event Numeric vector indicating event status, where
#'   `1`  denotes an observed event and `0`  denotes censoring.
#' @param surv_covar Numeric vector or matrix of survival covariates. If a vector
#'   is supplied, it is converted to a one-column matrix. If surv_covar = NULL,
#'   it is assumed that there is no survival covariate and creates a vector of 0's,
#'   which results in cox regression coefficient forced to be 0.
#' @param iter_max Integer maximum number of EM iterations used by
#'   \code{mlm_factor_em}. Default is `500`.
#' @param tol Numeric convergence tolerance passed to \code{mlm_factor_em}.
#'   Default is `1e-6`.
#' @param verbose Logical; if `TRUE`, prints progress from the factor-model EM
#'   algorithm. Default is `FALSE`.
#'
#' @return A list of initial values suitable for the `init` argument of
#' \code{\link{slate}}, with components:
#' \describe{
#'   \item{W}{Initial factor loading matrix.}
#'   \item{psi}{Initial feature-specific residual variances.}
#'   \item{coef_surv}{Initial Cox coefficients for observed survival covariates.}
#'   \item{coef_latent}{Initial Cox coefficients for latent factor expectations.}
#'   \item{risk_sets}{Risk-set information at Cox event times.}
#'   \item{haz_subj}{Subject-specific baseline hazard values expanded to each subject.}
#'   \item{cumhaz_subj}{Subject-specific cumulative baseline hazard values expanded to each subject.}
#'   \item{cox_event_time}{Event times from the initial Cox model fit.}
#' }
#'
#' @details
#' The function first fits a multivariate latent factor model using
#' \code{mlm_factor_em}. It then computes expected latent factors with
#' \code{mlm_Ex} and fits a Cox proportional hazards model:
#' \deqn{
#'   h(t) = h_0(t)\exp(\gamma^\top Z + \beta^\top E[X \mid Y]).
#' }
#' Risk sets, baseline hazards, and cumulative hazards are extracted from the Cox
#' model and returned for reuse by \code{\link{slate}}.
#'
#' @seealso
#' \code{\link{slate}}, \code{\link[survival]{coxph}},
#' \code{\link[survival]{Surv}}
#'
#' @export
slate_init <- function(Y, Zt, surv_time, surv_event, surv_covar=NULL,
                       iter_max=500, tol=1e-6, verbose=FALSE){

  
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


  ## full factor em takes longer to get init value but reduces number of iterations below
  fit_factor <- mlm_factor_em(Y,Zt, nlatent=nlatent, iter_max = iter_max, tol=tol, verbose=verbose)
  fit_Ex <- mlm_Ex(Y, Zt, fit_factor$W, fit_factor$psi)
  fit_cox <- coxph(Surv(surv_time, surv_event) ~ surv_covar + fit_Ex$Ex)

  cox_detail <- coxph.detail(fit_cox)
  cox_detail$cum_hazard <- cumsum(cox_detail$hazard)
  cox_event_time <- cox_detail$time

  risk_sets <-vector("list", length=length(cox_event_time))
  index_subj <- 1:length(surv_time)
  for(i in 1:length(cox_event_time)){
    risk_sets[[i]] <- list(subj= index_subj[surv_time >= cox_event_time[i]],
                           event= index_subj[surv_time == cox_event_time[i] & surv_event == 1])
  }


  haz_subj <- expand_hazard(cox_detail$hazard, surv_time, cox_event_time)
  cumhaz_subj <- expand_hazard(cox_detail$cum_hazard, surv_time, cox_event_time)
  ncovar <- ncol(surv_covar)
  coef <- fit_cox$coefficients
  coef[is.na(coef)] <- 0

  coef_surv <- coef[1:ncovar]
  coef_latent <- coef[(ncovar+1):(ncovar+nlatent)]
  W <- fit_factor$W
  psi <- fit_factor$psi

  return(list( W=W, psi=psi, coef_surv=coef_surv, coef_latent=coef_latent,
               risk_sets=risk_sets, haz_subj=haz_subj, cumhaz_subj = cumhaz_subj,
               cox_event_time=cox_event_time))
}


#' Summarize a Fitted `slate` Model
#'
#' Prints a summary of a fitted `slate` model, including the estimated
#' coefficients, their standard errors, Wald test statistics, p-values,
#' estimated latent loadings (`W`), convergence information, and the Bayesian
#' Information Criterion (BIC).
#'
#' If all latent loadings are zero (i.e., `W = 0` for every loading), the
#' latent coefficient estimates are not identifiable. In this case, the
#' corresponding coefficient entries are printed as `NA`.
#'
#' @param object A fitted `slate` model object returned by [slate()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return
#' This function is called for its side effect of printing a model summary.
#' It returns its input object invisibly.
#'
#' @seealso [slate()]
#'
#' @export
summary.slate <- function(object, ...) {
  fit <- object
  
  cat("Coefficients:\n")
  tbl <- rbind(coef=fit$coef, se=sqrt(diag(fit$var_coef)), zstat=fit$zstat, pval=fit$pval)
  print(tbl)
  cat("\n")
  cat("W loadings:\n")
  print(fit$W)
  cat("\n")
  cat("iter = ", fit$iter, ", ")
  cat("converge = ", fit$converge, "\n")
  cat("bic = ", fit$bic, "\n")
  
  invisible(tbl)
}

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
  cat("Computing initial values\n")
  init <- slate_init(Y, Zt, surv_time, surv_event, surv_covar, iter_max=500, tol=1e-6, verbose=FALSE)


  for(i in 1:nalpha){
    alpha <- alpha_grid[i]
    max_lambda <- max_XYZ/alpha

    for(j in 1:nlambda){

      lambda <- lambda_frac[j]*max_lambda

      cat("================ alpha = ", alpha, ", lambda.frac = ", lambda_frac[j], ", lambda = ", lambda, " ==================\n")
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


#' Update Cox regression coefficients for the SLATE model
#'
#' Performs iterative Newton-Raphson updates for Cox proportional hazards
#' regression coefficients using expected latent-variable quantities from the
#' SLATE model.
#'
#' @param coef Numeric vector of initial Cox regression coefficients, including
#'   coefficients for observed survival covariates and latent variables.
#' @param surv_covar Numeric matrix of observed survival covariates.
#' @param surv_event Numeric vector indicating event status, where
#'   `1` denotes an observed event and `0` ` denotes censoring.
#' @param cumhaz_subj Numeric vector of subject-specific cumulative baseline
#'   hazards.
#' @param risk_sets List of risk-set information at each event time.
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param W Numeric factor loading matrix.
#' @param psi Numeric vector of residual variances for observed features.
#' @param Zt Design matrix  for the multivariate features.
#' @param ex.dat List of expected latent-variable quantities returned by
#'   \code{expected_func_x}.
#' @param gh.dat List containing Gauss-Hermite quadrature nodes and weights.
#' @param tol Numeric convergence tolerance for coefficient updates.
#'
#' @return A numeric vector of updated Cox regression coefficients.
#'
#' @details
#' The function iteratively updates the Cox regression coefficients using a
#' Newton-Raphson procedure:
#' \deqn{
#'   \beta^{(k+1)} =
#'   \beta^{(k)} + I(\beta^{(k)})^{-1} U(\beta^{(k)}),
#' }
#' where \eqn{U(\beta)} is the expected score vector and
#' \eqn{I(\beta)} is the observed information matrix returned by
#' \code{expected_cox_score}.
#'
#' At each iteration:
#' \enumerate{
#'   \item The expected Cox score and information matrix are computed.
#'   \item Coefficients are updated using the generalized inverse of the
#'     information matrix.
#'   \item Missing values in the updated coefficients are replaced with zero.
#'   \item Iteration stops when the maximum absolute coefficient change is less
#'     than `tol`, or after 20 iterations.
#' }
#'
#' @seealso
#' \code{\link{expected_cox_score}}, \code{\link[MASS]{ginv}}
update_cox_coef <- function(coef, surv_covar, surv_event, cumhaz_subj, risk_sets,
                            Y, W, psi, Zt, ex.dat, gh.dat, tol){
 coef[is.na(coef)] <- 0

  coef_old <- coef
  for(iter_cox in 1:20){
    score_dat <- expected_cox_score(surv_covar, surv_event, coef_old, cumhaz_subj, risk_sets,
                                    Y, W, psi, Zt, ex.dat, gh.dat)

    coef_new <- coef_old + ginv(score_dat$info) %*% score_dat$score

    coef_new[is.na(coef_new)] <- 0
   
    coef_change <- max(abs(coef_new - coef_old))

    if(coef_change < tol){
      break
    }
    coef_old <- coef_new
  }
  return(coef_new)
}

#' Expand hazard values to subject-level observations
#'
#' Expands hazard or cumulative hazard values defined at event times to
#' subject-specific values aligned with individual survival times.
#'
#' @param hazard Numeric vector of hazard or cumulative hazard values evaluated
#'   at unique event times.
#' @param surv_time Numeric vector of subject-specific follow-up or event times.
#' @param event_time Numeric vector of unique event times corresponding to the
#'   entries in `hazard`.
#'
#' @return A numeric vector of subject-specific hazard values with the same
#' length and ordering as `surv_time`.
#'
#' @details
#' The function matches each subject's survival time to the corresponding event
#' time hazard value using \code{match}. Subjects whose times do not exactly
#' match an event time are assigned the most recent non-missing hazard value
#' carried forward in sorted survival-time order.
#'
#' The procedure is:
#' \enumerate{
#'   \item Match subject survival times to event-time hazard values.
#'   \item Sort subjects by survival time.
#'   \item Carry the most recent observed hazard value forward to fill missing
#'     values.
#'   \item Restore the original subject ordering.
#' }
#'
#' This is commonly used to map baseline hazard or cumulative hazard estimates
#' from Cox-model event times to all study subjects.
expand_hazard <- function(hazard, surv_time, event_time){
  ord <- order(surv_time)
  haz2 <- hazard[match(surv_time, event_time)]
  haz2 <- haz2[ord]
  last <- 0 # cumulative hazard is 0 before the first event
  for(i in seq_along(haz2)){
    if(!is.na(haz2[i])) last <- haz2[i]
    if(is.na(haz2[i]))  haz2[i] <- last
  }
  haz_subj <- haz2[order(ord)]
  return(haz_subj)
}


#' Update baseline hazard estimates
#'
#' Computes baseline hazard increments at each event time using current survival
#' covariate effects and expected latent-risk contributions.
#'
#' @param surv_covar Numeric vector or matrix of observed survival covariates.
#' @param coef_surv Numeric vector of Cox regression coefficients for
#'   `surv_covar`.
#' @param Exb Numeric vector of expected latent-risk contributions, typically
#'   returned by \code{expected_func_x}.
#' @param risk_sets List of risk-set information. Each element should contain
#'   `subj`, the subject indices in the risk set, and `event`, the subject
#'   indices with events at that event time.
#'
#' @return A numeric vector of baseline hazard increments, one value for each
#'   element of `risk_sets`.
#'
#' @details
#' For each event time, the baseline hazard increment is computed as
#' \deqn{
#'   \hat h_0(t_j) =
#'   \frac{d_j}{\sum_{i \in R_j}
#'   \exp(Z_i^\top \gamma) E[\exp(X_i^\top \beta)]},
#' }
#' where `d_j` is the number of observed events at event time `t_j`, `R_j` is
#' the corresponding risk set, `surv_covar` and `coef_surv` define the observed
#' covariate risk score, and `Exb` provides the expected latent-risk term.
#'
#' @seealso
#' \code{\link{expand_hazard}}, \code{\link{expected_func_x}}
update_hazard <- function(surv_covar, coef_surv, Exb, risk_sets){

  rrisk_adj <- exp(surv_covar %*% coef_surv);
  rrisk_tot <- rrisk_adj * Exb
  haz <- rep(0, length(risk_sets))
  for(iset in 1:length(risk_sets)){
    nevents <- length(risk_sets[[iset]]$event)
    ## sums over indices within a risk set
    indices <- risk_sets[[iset]]$subj
    den <- sum(rrisk_tot[indices])
    haz[iset] <- nevents / den
  }
  return(haz)
}

#' Compute the maximum penalty scaling value for SLATE at which all W = 0
#'
#' Computes the maximum value of the penalization term used to determine the
#' largest regularization parameter in the SLATE model fitting procedure.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Design matrix for the multivariate features.
#' @param nlatent Integer number of latent factors.
#'
#' @return A single numeric value giving the maximum scaled absolute value of the
#'   `XYZ` update term:
#'   \deqn{
#'     \max |\mathbf{X}^\top \mathbf{Y} \mathbf{Z}^\top| / n.
#'   }
#'
#' @details
#' The function first estimates an initial latent factor model using
#' \code{mlm_mle_pca}. Residual variances are constructed from the estimated
#' common variance parameter and used to standardize the observed feature matrix
#' and design matrix.
#'
#' Expected latent variables are then computed with \code{mlm_Ex}, and the
#' following matrices are formed:
#' \describe{
#'   \item{`XY`}{Cross-product between expected latent variables and scaled features.}
#'   \item{`ZZ`}{Scaled design covariance matrix.}
#'   \item{`XYZ`}{Cross-product term used in penalized loading updates.}
#' }
#'
#' The returned quantity is the largest absolute entry of `XYZ`, scaled by the
#' number of subjects. This value is commonly used to define the maximum
#' regularization parameter:
#' \deqn{
#'   \lambda_{\max} = \frac{\max |XYZ|}{n \alpha}.
#' }
#'
#' @seealso
#' \code{\link{slate_grid}}, \code{\link{mlm_mle_pca}},
#' \code{\link{mlm_Ex}}
compute_max_XYZ <- function(Y, Zt,  nlatent){
  fit_pca <- mlm_mle_pca(Y,Zt, nlatent)
  W <- fit_pca$W
  sigma <- fit_pca$sigma
  nfeature <- ncol(Y)
  psi <- rep(sigma, nfeature)
  fit_Ex <- mlm_Ex(Y, Zt, W, psi)
  ## min lambda such that all W are not selected

  ## update terms that depend on psi and W
  Ys <- Y /sqrt(psi)
  Zts <- Zt / sqrt(psi)
  Ex <- fit_Ex$Ex
  XY <- t(Ex) %*% Ys
  ZZ <- Zts %*% t(Zts)
  XYZ <- XY %*% t(Zts)

  nsubj <- nrow(Y)

  max_XYZ <- max(abs(XYZ))/(nsubj)
  return(max_XYZ)
}

#' Compute expected latent-variable functions
#'
#' Computes conditional expectations of latent variables and latent-risk
#' functions given observed features, survival data, and current SLATE parameter
#' estimates.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param W Numeric factor loading matrix.
#' @param psi Numeric vector of feature-specific residual variances.
#' @param Zt Design matrix for the multivariate features.
#' @param surv_cov Numeric vector or matrix of observed survival covariates.
#' @param coef_surv Numeric vector of Cox regression coefficients for
#'   `surv_cov`.
#' @param coef_latent Numeric vector of Cox regression coefficients for the
#'   latent variables.
#' @param surv_event Numeric vector indicating event status, where
#'   `1` denotes an observed event and `0` denotes censoring.
#' @param cumhaz_subj Numeric vector of subject-specific cumulative baseline
#'   hazards.
#' @param gh.dat List of Gauss-Hermite quadrature information, typically with
#'   components `node` and `wt`.
#'
#' @return A list with components:
#' \describe{
#'   \item{Ex}{Matrix of conditional expectations `E[X | Y, T, delta]`.}
#'   \item{Exb}{Matrix/vector of conditional expectations
#'     `E[exp(X beta) | Y, T, delta]`.}
#'   \item{Exbx}{Matrix of conditional expectations
#'     `E[exp(X beta) X | Y, T, delta]`.}
#'   \item{Exx}{Array of conditional second moments
#'     `E[X X' | Y, T, delta]`, indexed by subject.}
#'   \item{Exbxx}{Array of conditional weighted second moments
#'     `E[exp(X beta) X X' | Y, T, delta]`, indexed by subject.}
#' }
#'
#' @details
#' The function first computes the conditional Gaussian distribution of latent
#' variables given the observed multivariate features. It then updates this
#' distribution using survival likelihood contributions from the Cox model.
#'
#' The final conditional expectations are evaluated by numerical integration
#' over the Gauss-Hermite quadrature nodes in `gh.dat`.
#'
#' Missing values in `W` are replaced by zero, and missing values in `psi` are
#' replaced by one before the calculations.
#'
#' @seealso
#' \code{\link{slate}}, \code{\link{gauss_hermite_quad}},
#' \code{\link{expected_cox_score}}
expected_func_x <- function(Y, W, psi, Zt, surv_cov, coef_surv,coef_latent,
                            surv_event, cumhaz_subj, gh.dat){
  node <- gh.dat$node
  nnode <- nrow(node)
  wt <- gh.dat$wt
  nsubj <- nrow(Y)
  nlatent <- nrow(W)

  W[is.na(W)] <- 0
  psi[is.na(psi)] <- 1

  Ex <- Exb <- Exbx <-  matrix(0, nrow=nsubj, ncol=nlatent)
  Exx <- Exbxx <- array(0, dim=c(nlatent, nlatent, nsubj))

  # conditional moments of X given Y
  imat <- diag(nlatent)
  psi_inv = diag(1.0 / psi)
  WZt <- W %*% Zt
  inv_lowdim <- ginv(imat + WZt %*% psi_inv %*% t(WZt))
  vinv <-  psi_inv - psi_inv %*% t(WZt) %*% inv_lowdim %*% WZt %*% psi_inv
  Ex_y<-   Y %*% t(WZt %*% vinv)
  var_x_y <-  imat - WZt %*% vinv %*% t(WZt)
  vinv_x_y <- ginv(var_x_y)

  xb <- surv_cov %*% coef_surv
  f1 <- f2 <- rep(0, nrow(node))

  ## conditional moments of X given Y, surv_time, surv_event and current
  ## parameter estimates


  for(isubj in 1:nsubj){

    f1 <-  -cumhaz_subj[isubj] * exp(xb[isubj]  + node %*%coef_latent)

    if (surv_event[isubj] == 1) {
      f1  <- f1 + (node %*%coef_latent);
    }

    f2 <- f2*0
    for(j in 1:nnode){
      delta <- node[j,] - Ex_y[isubj]
      f2[j] = -.5* delta %*% vinv_x_y %*% delta
    }

    ## weight for conditional density of X
    wtf = exp((log(wt) + f1 + f2));
    wtf[is.na(wtf)] <- 0

    ## densum to standardize
    densum = sum(wtf);


    if(is.na(densum) | densum < 1e-8){
      wtf = wtf + 1e-08;
      densum = sum(wtf);
    }


    ## E[X]
    Ex[isubj] <- (t(wtf) %*% node) / densum

    ## E[exp(x*b) * x]
    Exbx[isubj] <- (t(wtf) %*% (exp(node %*%coef_latent) * node))/densum

    ## E[exp(x*b)]
    Exb[isubj] <- (t(wtf) %*% (exp(node %*%coef_latent)))/densum

    ## E[x * x']
    ##  E[exp(x*b) * x * x']
    Exx[,,isubj] = 0 * Exx[,,isubj]

    for (i in 1:nnode) {
      Exx[,,isubj] <-  Exx[,,isubj] +
        wtf[i] * (node[i,] %o% node[i,])
    }
    Exx[,,isubj] <- Exx[,,isubj] / densum;

    ##  E[exp(x*b) * x * x']
    Exbxx[,,isubj] = 0 * Exbxx[,,isubj]
    tmp <- exp(node %*%coef_latent)
    for (i in 1:nnode) {
      Exbxx[,,isubj] <-  Exbxx[,,isubj] +
        wtf[i] * (node[i,] %o% node[i,]) * tmp[i]
    }
    Exbxx[,,isubj] <- Exbxx[,,isubj] / densum;

  }

  return(list(Ex=Ex, Exb=Exb, Exbx=Exbx, Exx=Exx, Exbxx = Exbxx))
}

#' Fit a multivariate latent factor model using EM
#'
#' Fits a multivariate latent factor model of the form
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' using an expectation-maximization (EM) algorithm.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix of feature annotations.
#' @param nlatent Integer number of latent variables (factors) to estimate.
#' @param iter_max Integer maximum number of EM iterations. Default is `500`.
#' @param tol Numeric convergence tolerance for the average log-likelihood.
#'   Default is `1e-6`.
#' @param verbose Logical; if `TRUE`, prints iteration progress information.
#'   Default is `FALSE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{W}{Estimated latent-factor loading matrix.}
#'   \item{psi}{Estimated feature-specific residual variances.}
#'   \item{niter}{Number of EM iterations performed.}
#'   \item{converge}{Logical indicating whether convergence was achieved.}
#' }
#'
#' @details
#' The model assumes:
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where:
#' \describe{
#'   \item{`Y`}{Observed `n x p` feature matrix.}
#'   \item{`X`}{Latent variable matrix with `nlatent` columns.}
#'   \item{`W`}{Latent-factor coefficient matrix.}
#'   \item{`Z`}{Feature annotation/design matrix.}
#'   \item{`E`}{Residual errors with diagonal covariance matrix
#'     `diag(psi)`.}
#' }
#'
#' Initial values are obtained from \code{mlm_mle_pca}. The EM algorithm then
#' alternates between:
#' \enumerate{
#'   \item Computing conditional expectations of latent variables using
#'     \code{mlm_Ex}.
#'   \item Updating the loading matrix `W`.
#'   \item Updating residual variances `psi`.
#'   \item Evaluating the normalized log-likelihood using
#'     \code{mlm_factor_lnlike}.
#' }
#'
#' Convergence is declared when the absolute change in normalized log-likelihood
#' satisfies:
#' \deqn{
#'   |\ell^{(k+1)} - \ell^{(k)}|
#'   < \texttt{tol} \times (|\ell^{(k)}| + 1).
#' }
#'
#' Missing values in `W` are replaced by zero, and missing values in `psi` are
#' replaced by one during updates.
#'
#' @seealso
#' \code{\link{mlm_mle_pca}}, \code{\link{mlm_Ex}},
#' \code{\link{mlm_factor_lnlike}}
#'
#' @export
mlm_factor_em <- function(Y,Zt, nlatent, iter_max = 500, tol=1e-6, verbose=FALSE){
  ## em for fitting x latent matrix as Y = xWz' + E and

  ## Y = XWZ' + diag(sig)
  ## Y = nxp matrix centered and scaled of p features
  ## X = PC's (latent vars), p x l
  ## Z' = mxp matrix of m annotations for the p features
  ## W = l x m matrix of PC's constrained by Y =XWZ' + diag(sig)
  ## nlatent = number of PC's (nlatent = no. col of X = l)

  ## estimate W and psi:

  nfeat <- ncol(Y)
  p <- nfeat
  nsubj <- nrow(Y)
  n <- nsubj
  ## init values by mle of pca
  fit_pca <- mlm_mle_pca(Y,Zt, nlatent)
  W <- fit_pca$W
  sigma <- fit_pca$sigma


  psi <- rep(sigma, nfeat)

  converge <- FALSE

  lnlike_old <- mlm_factor_lnlike(Y, Zt, W, psi)/(n*p)

  for(iter in 1:iter_max){

    ## update terms that depend on psi and W

    exval <- mlm_Ex(Y, Zt, W, psi)
    Ex <- exval$Ex
    Exx_sum <- exval$Exx_sum

    Ys <- Y /sqrt(psi)
    Zts <- Zt / sqrt(psi)

    XY <- t(Ex) %*% Ys
    ZZ <- Zts %*% t(Zts)
    XYZ <- XY %*% t(Zts)

    W <- ginv(Exx_sum) %*% XYZ %*% ginv(ZZ)
    W[is.na(W)] <- 0

    ## update psi (note use Y and Z,not Ys or Zs)
    WZt <- W %*% Zt
    RR <- t(Y) %*% Y - 2*t(Y) %*% Ex %*% WZt + t(WZt) %*% Exx_sum %*% WZt
    psi <- diag(RR)/n
    psi[is.na(psi)] <- 1

    lnlike_new <- mlm_factor_lnlike(Y, Zt, W, psi)/(n*p)

    lnlike_delta <- lnlike_new - lnlike_old
    if(verbose)
    {
      cat("mlm_factor_em iter = ", iter, ", lnlike = ", lnlike_new, ", delta = ", abs(lnlike_delta), ", tol = ", tol * (abs(lnlike_old) + 1) , "\n")
    }

    if(is.na(lnlike_delta)){
      lnlike_delta <- 1e-8
      cat("init lnlike_delta is NA\n")
    }


    if ( abs(lnlike_delta) < (tol * (abs(lnlike_old) + 1)))  {
      converge <- TRUE
      break
    }

    lnlike_old <- lnlike_new
  }
  return(list(W = W, psi=psi, niter=iter, converge=converge))

}

#' Compute maximum likelihood estimates for probabilistic PCA
#'
#' Computes initial parameter estimates for a probabilistic principal component
#' model with annotation-constrained loadings.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix for features.
#' @param nlatent Integer number of latent variables (principal components)
#'   to estimate.
#'
#' @return A list with components:
#' \describe{
#'   \item{W}{Estimated latent-factor loading matrix.}
#'   \item{sigma}{Estimated residual variance parameter.}
#' }
#'
#' @details
#' The function computes the maximum likelihood estimates for probabilistic PCA
#' following Equation 12.45 from Bishop (Pattern Recognition and Machine
#' Learning).
#'
#' The model assumes:
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where:
#' \describe{
#'   \item{`Y`}{Observed `n x p` feature matrix.}
#'   \item{`X`}{Latent variable matrix with `nlatent` columns.}
#'   \item{`W`}{Latent-factor coefficient matrix.}
#'   \item{`Z`}{Feature annotation/design matrix.}
#'   \item{`E`}{Residual Gaussian noise with variance `sigma`.}
#' }
#'
#' The procedure is:
#' \enumerate{
#'   \item Compute the sample covariance matrix of `Y`.
#'   \item Perform eigendecomposition of the covariance matrix.
#'   \item Retain the leading `nlatent` eigenvectors and eigenvalues.
#'   \item Estimate the residual variance from the remaining eigenvalues.
#'   \item Construct the loading matrix `W` using the annotation matrix `Zt`.
#' }
#'
#' Missing values in `W` are replaced by zero, and missing values in `sigma`
#' are replaced by one.
#'
#' @references
#' Bishop, C. M. (2006). \emph{Pattern Recognition and Machine Learning}.
#' Springer. Equation 12.45.
#'
#' @seealso
#' \code{\link{mlm_factor_em}}, \code{\link[MASS]{ginv}}
#'
#' @export
mlm_mle_pca <- function(Y,Zt, nlatent){
  ## MLE of probabilistic principal components (eq 12.45, p 574 of Bishop)
  ## Y = XWZ' + diag(sig)
  ## Y = nxp matrix centered and scaled of p features
  ## X = PC's (latent vars), p x l
  ## Z' = mxp matrix of m annotations for the p features
  ## W = l x m matrix of PC's constrained by Y =XWZ' + diag(sig)
  ## nlatent = number of PC's (nlatent = no. col of X = l)

  cv <- var(Y)
  ev <- eigen(cv)

  U <- ev$vectors[, 1:nlatent, drop=FALSE]

  if(nlatent == 1){
    L <- matrix(ev$values[1], 1, 1)
  } else{
    L <- diag(ev$values[1:nlatent], nrow=nlatent)
  }

  sig <- sum( ev$values[(nlatent+1):length(ev$values)] )
  sig <- sig / (length(ev$values) - nlatent)

  Ident <- diag(rep(1, nlatent))
  zz_inv <- ginv(Zt %*% t(Zt))
  W <- (L - sig*Ident)^.5 %*% t(U) %*% t(Zt) %*% zz_inv

  rownames(W) <- paste("x", 1:nlatent, sep=".")
  colnames(W) <- rownames(Zt)

  W[is.na(W)] <- 0
  sig[is.na(sig)] <- 1

  return(list(W=W, sigma=sig))
}

#' Compute the latent factor model log-likelihood
#'
#' Computes the observed-data log-likelihood for the multivariate latent factor
#' model used in the SLATE framework.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix of feature annotations.
#' @param W Numeric latent-factor loading matrix.
#' @param psi Numeric vector of feature-specific residual variances.
#'
#' @return A numeric scalar giving the log-likelihood value for the latent factor
#'   model.
#'
#' @details
#' The model assumes:
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where:
#' \describe{
#'   \item{`Y`}{Observed feature matrix.}
#'   \item{`X`}{Latent variable matrix.}
#'   \item{`W`}{Latent-factor loading matrix.}
#'   \item{`Z`}{Feature annotation/design matrix.}
#'   \item{`E`}{Residual Gaussian errors with diagonal covariance matrix
#'     `diag(psi)`.}
#' }
#'
#' The function computes the conditional moments of the latent variables and
#' evaluates the Gaussian log-likelihood contribution for each subject.
#'
#' The covariance inverse is computed using the Woodbury matrix identity:
#' \deqn{
#'   V^{-1} =
#'   \Psi^{-1} -
#'   \Psi^{-1} W Z^\top
#'   (I + W Z^\top \Psi^{-1} Z W^\top)^{-1}
#'   W Z^\top \Psi^{-1}.
#' }
#'
#' Missing values in `W` are replaced by zero, and missing values in `psi`
#' are replaced by one before calculations.
#'
#' @seealso
#' \code{\link{mlm_factor_em}}, \code{\link{mlm_Ex}},
#' \code{\link[MASS]{ginv}}
mlm_factor_lnlike <- function(Y, Zt, W, psi){
  W[is.na(W)] <- 0
  psi[is.na(psi)] <- 1

  nsubj <- nrow(Y)
  WZt <- W %*% Zt
  psi_inv <- diag(1/psi)
  nfeat <- ncol(Y)
  nlatent <- nrow(W)
  imat <- diag(1, nrow=nlatent)

  vinv <- psi_inv - psi_inv %*% t(WZt) %*% ginv(imat + WZt %*% psi_inv %*% t(WZt)) %*% WZt %*% psi_inv

  mu_factor <- WZt %*% vinv
  Var <- imat - WZt %*% vinv %*% t(WZt)

  mu_factor[is.na(mu_factor)] <- 0
  Var[is.na(Var)] <- 0

  lnlike <- 0
  for(i in 1:nsubj){

    ez <- mu_factor %*% Y[i,]
    ezz <- Var + ez %*% t(ez)

    ez[is.na(ez)] <- 0
    ezz[is.na(ezz)] <- 0

    lnlike <- lnlike + sum(Y[i,]*Y[i,]/psi) - 2*(Y[i,]/psi) %*% t(WZt) %*% ez +
      sum( diag(t(WZt) %*% ezz %*% WZt)/psi)

  }

  lnlike <- -.5*lnlike - .5*nsubj*sum(log(psi))
  lnlike <- as.vector(lnlike)
  return(lnlike)
}

#' Compute the penalized latent factor objective function
#'
#' Computes the elastic-net penalized negative expected log-likelihood for the
#' latent factor model used in the SLATE framework.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix of feature annotations.
#' @param W Numeric latent-factor loading matrix.
#' @param psi Numeric vector of feature-specific residual variances.
#' @param Ex Matrix of conditional expectations of latent variables.
#' @param Exx_sum Numeric matrix giving the summed conditional second moments of
#'   the latent variables.
#' @param alpha Elastic-net mixing parameter. `alpha = 1` corresponds to a pure
#'   L1 penalty and `alpha = 0` corresponds to a pure L2 penalty.
#' @param lambda Nonnegative regularization parameter controlling the penalty
#'   strength.
#'
#' @return A numeric scalar giving the penalized negative expected
#'   log-likelihood value.
#'
#' @details
#' The latent factor model assumes:
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where:
#' \describe{
#'   \item{`Y`}{Observed `n x p` feature matrix.}
#'   \item{`X`}{Latent variable matrix.}
#'   \item{`W`}{Latent-factor loading matrix.}
#'   \item{`Z`}{Feature annotation/design matrix.}
#'   \item{`E`}{Residual Gaussian noise with diagonal covariance matrix
#'     `diag(psi)`.}
#' }
#'
#' The function computes the normalized expected quadratic loss:
#' \deqn{
#'   \frac{1}{2np}
#'   \left(
#'   \mathrm{tr}(Y^\top Y)
#'   - 2\mathrm{tr}(W^\top E[X]^\top Y)
#'   + \mathrm{tr}(W^\top E[XX^\top] W)
#'   \right),
#' }
#' after scaling by the residual variances.
#'
#' Elastic-net regularization is then added:
#' \deqn{
#'   \lambda \alpha \sum |W|
#'   +
#'   \frac{\lambda (1-\alpha)}{2} \sum W^2.
#' }
#'
#' Missing values in intermediate matrices are replaced with zero before
#' calculations.
#'
#' @seealso
#' \code{\link{mlm_factor_em}}, \code{\link{update_W}},
#' \code{\link{trace_mat}}
mlm_factor_penfunc <- function(Y, Zt, W, psi, Ex, Exx_sum, alpha, lambda){
  ## elastic net penalized neg expected lnlike for
  ## factor model for Y = XWZ' + E
  ## Y is n x p
  ## X is latent n x l
  ## W is loadings l x m
  ## Z' is design matrix for cols of Y, m x p

  n <- nrow(Y)
  p <- ncol(Y)

  Ys <- Y /sqrt(psi)
  Zs <- Zt / sqrt(psi)

  WZ <- W %*% Zs

  WZ[is.na(WZ)] <- 0
  Exx_sum[is.na(Exx_sum)] <- 0

  f <- trace_mat(t(Ys) %*% Ys) -2*trace_mat(t(WZ) %*% t(Ex) %*% Ys) + trace_mat(t(WZ) %*% Exx_sum %*% WZ)
  f <- f / (2*n*p)
  f <- f + lambda*alpha * sum(abs(W))
  f <- f + lambda*(1-alpha) * sum(W^2)/2

  return(f)
}

#' Compute the latent factor model loss function
#'
#' Computes the normalized expected quadratic loss for the latent factor model
#' used in the SLATE framework.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix of feature annotations.
#' @param W Numeric latent-factor loading matrix.
#' @param psi Numeric vector of feature-specific residual variances.
#' @param Ex Matrix of conditional expectations of latent variables.
#' @param Exx_sum Numeric matrix giving the summed conditional second moments of
#'   the latent variables.
#' @param alpha Elastic-net mixing parameter. Included for interface consistency
#'   but not directly used in the loss computation.
#' @param lambda Nonnegative regularization parameter. Included for interface
#'   consistency but not directly used in the loss computation.
#'
#' @return A numeric scalar giving the normalized expected quadratic loss.
#'
#' @details
#' The latent factor model assumes:
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where:
#' \describe{
#'   \item{`Y`}{Observed `n x p` feature matrix.}
#'   \item{`X`}{Latent variable matrix.}
#'   \item{`W`}{Latent-factor loading matrix.}
#'   \item{`Z`}{Feature annotation/design matrix.}
#'   \item{`E`}{Residual Gaussian noise with diagonal covariance matrix
#'     `diag(psi)`.}
#' }
#'
#' The function evaluates the normalized expected quadratic reconstruction loss:
#' \deqn{
#'   \frac{1}{2n}
#'   \left(
#'   \mathrm{tr}(Y^\top Y)
#'   - 2\mathrm{tr}(W^\top E[X]^\top Y)
#'   + \mathrm{tr}(W^\top E[XX^\top] W)
#'   \right),
#' }
#' after scaling by the residual variances.
#'
#' Unlike \code{\link{mlm_factor_penfunc}}, this function does not include the
#' elastic-net penalty term in the returned value.
#'
#' Missing values in `WZ`, `Exx_sum`, and `Ex` are replaced with zero before
#' calculations.
#'
#' @seealso
#' \code{\link{mlm_factor_penfunc}}, \code{\link{mlm_factor_em}},
#' \code{\link{trace_mat}}
mlm_factor_loss <- function(Y, Zt, W, psi, Ex, Exx_sum, alpha, lambda){
  ## elastic net penalized neg expected lnlike for
  ## factor model for Y = XWZ' + E
  ## Y is n x p
  ## X is latent n x l
  ## W is loadings l x m
  ## Z' is design matrix for cols of Y, m x p

  n <- nrow(Y)
  p <- ncol(Y)

  Ys <- Y /sqrt(psi)
  Zs <- Zt / sqrt(psi)

  WZ <- W %*% Zs

  WZ[is.na(WZ)] <- 0
  Exx_sum[is.na(Exx_sum)] <- 0
  Ex[is.na(Ex)] <- 0

  f <- trace_mat(t(Ys) %*% Ys) -2*trace_mat(t(WZ) %*% t(Ex) %*% Ys) + trace_mat(t(WZ) %*% Exx_sum %*% WZ)
  loss <- f / (2*n)


  return(loss)
}

#' Fit a penalized latent factor model by cyclic coordinate descent
#'
#' Fits the latent factor model using cyclic coordinate descent updates for the
#' loading matrix with an elastic-net penalty.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix of feature annotations.
#' @param nlatent Integer number of latent variables to estimate.
#' @param init Optional list of initial values with components `W` and `psi`.
#'   If `NULL`, initial values are obtained from \code{\link{mlm_factor_em}}.
#' @param alpha Elastic-net mixing parameter. `alpha = 1` gives an L1 penalty
#'   and `alpha = 0` gives an L2 penalty.
#' @param lambda Nonnegative regularization parameter.
#' @param maxiter Integer maximum number of coordinate-descent iterations.
#'   Default is `500`.
#' @param tol Numeric convergence tolerance. Default is `1e-5`.
#' @param verbose Logical; if `TRUE`, prints iteration progress. Default is
#'   `FALSE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{W}{Estimated latent-factor loading matrix.}
#'   \item{psi}{Estimated feature-specific residual variances.}
#'   \item{loss}{Final penalized objective value.}
#'   \item{lnlike}{Final latent factor model log-likelihood.}
#'   \item{bic}{Bayesian information criterion.}
#'   \item{df}{Estimated model degrees of freedom.}
#'   \item{iter}{Number of iterations performed.}
#'   \item{converge}{Logical indicating whether convergence was achieved.}
#'   \item{alpha}{Elastic-net mixing parameter used.}
#'   \item{lambda}{Regularization parameter used.}
#' }
#'
#' @details
#' The fitted model is
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where `X` is an unobserved latent matrix, `W` is a loading matrix, and `E`
#' has diagonal covariance `diag(psi)`.
#'
#' Given current conditional moments of the latent variables, each element of
#' `W` is updated using a soft-thresholded coordinate-descent step. After each
#' full sweep through `W`, residual variances `psi`, expected latent moments,
#' the log-likelihood, and the penalized objective are updated.
#'
#' Convergence is declared when either the relative change in the penalized
#' objective or the maximum absolute change in `W` is less than `tol`.
#'
#' @seealso
#' \code{\link{mlm_factor_em}}, \code{\link{mlm_Ex}},
#' \code{\link{mlm_factor_lnlike}}, \code{\link{mlm_factor_df}},
#' \code{\link{soft}}
#'
#' @export
mlm_factor_ccd <- function(Y, Zt, nlatent, init=NULL, alpha, lambda,
                           maxiter = 500, tol = 1e-5, verbose = FALSE) {

  ## elastic net penalized neg expected lnlike for
  ## factor model for Y = XWZ' + E
  ## Y is n x p
  ## X is latent n x nlatent
  ## W is loadings nlatent x m
  ## Z' is design matrix for cols of Y, m x p

  p <- ncol(Y)
  n <- nrow(Y)
  m <- nrow(Zt)

  if(is.null(init)){
    fit_init <- mlm_factor_em(Y, Zt, nlatent)
    W_curr<- fit_init$W
    psi <- fit_init$psi
  } else{
    W_curr  <- init$W
    psi <- init$psi
  }

  W_curr[is.na(W_curr)] <- 0
  psi[is.na(psi)] <- 1

  converge <- FALSE
  W_old <- W_curr

  exval <- mlm_Ex(Y, Zt, W_curr, psi)
  Ex <- exval$Ex
  Exx_sum <- exval$Exx_sum

  lnlike <- mlm_factor_lnlike(Y, Zt, W_curr, psi)

  ##  update penalty function
  penfunc_curr <- -lnlike/(n)  + alpha*lambda*sum(abs(W_curr))  + (1-alpha)*lambda*sum(W_curr^2)/2

  penfunc_old <- penfunc_curr

  for (iter in 1:maxiter) {

    ## update terms that depend on psi and W_curr
    Ys <- Y /sqrt(psi)
    Zts <- Zt / sqrt(psi)


    XY <- t(Ex) %*% Ys
    ZZ <- Zts %*% t(Zts)
    XYZ <- XY %*% t(Zts)
    M <- Exx_sum %*% W_curr %*% ZZ

    ## KKT test for zero entries
    ##grad_mat <- abs(M -XYZ) / (n*p)
    ##active <- abs(grad_mat) > alpha
    ##cat("max abs grad = ", max(abs(grad_mat)), "\n")

    ## update W
    for (irow in 1:nlatent) {
      for (jcol in 1:m) {

        ## skipping inactive doesn't seem to work here so eval all W
        ##if (!active[irow, jcol]) next

        W_ij_old <- W_curr[irow, jcol]

        numer <- (XYZ[irow, jcol] + Exx_sum[irow, irow] * ZZ[jcol,jcol] * W_curr[irow, jcol] - M[irow, jcol])/(n)
        denom <-(Exx_sum[irow, irow] * ZZ[jcol, jcol])/(n) + lambda * (1 - alpha)
        # update using soft-thresholding
        W_ij_new <- soft(numer, alpha * lambda) / denom


        if(is.na(W_ij_new)) W_ij_new <- 0

        W_curr[irow, jcol] <- W_ij_new
        delta <- W_ij_new - W_ij_old
        ## rank-1 update M after update beta
        M <- M + delta*(Exx_sum[,irow] %o% ZZ[jcol,])

      }
    }


    ## update psi (note use Y and Z,not Ys or Zs)
    WZ <- W_curr %*% Zt
    RR <- t(Y) %*% Y - 2*t(Y) %*% Ex %*% WZ + t(WZ) %*% Exx_sum %*% WZ
    psi <- diag(RR)/n

    psi[is.na(psi)] <- 1

    ## update expected values of X
    exval <- mlm_Ex(Y, Zt, W_curr, psi)
    Ex <- exval$Ex
    Exx_sum <- exval$Exx_sum
    ## update loss

    lnlike <- mlm_factor_lnlike(Y, Zt, W_curr, psi)
    penfunc_curr <- -lnlike/(n)  + alpha*lambda*sum(abs(W_curr))  + (1-alpha)*lambda*sum(W_curr^2)/2


    loss_rel_change <- abs(penfunc_curr - penfunc_old) / (abs(penfunc_old) + tol)
    W_change <- max(abs(W_curr - W_old))

    if (verbose) {
      cat("iter =", iter, ", penfunc_curr =", penfunc_curr, ", W range =", range(W_curr),
         ", W_change = ",  W_change, "\n")
    }
    if(is.na(loss_rel_change)){
      loss_rel_change <- 1e-8
      cat("iter loss_rel_change is NA\n")
    }
    if(is.na(W_change)){
      W[is.na(W)] <- 0
      cat("iter W_change is NA\n")
    }


    if (loss_rel_change < tol || W_change < tol) {
      converge <- TRUE
      break
    }


    W_old <- W_curr
    penfunc_old <- penfunc_curr

  }
  W <- W_curr
  lnlike <- mlm_factor_lnlike(Y, Zt, W, psi)
  df_model <- mlm_factor_df(W, psi)
  bic <- -2*lnlike + df_model*log(n)

  return(list(W = W, psi=psi, loss = penfunc_curr, lnlike=lnlike, bic=bic, df=df_model, iter= iter, converge = converge,
              alpha = alpha, lambda = lambda))
}

#' Update latent factor loadings by coordinate descent
#'
#' Updates the latent-factor loading matrix and residual variances for the SLATE
#' latent factor model using cyclic coordinate descent with an elastic-net
#' penalty.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric design matrix of feature annotations.
#' @param W Numeric current latent-factor loading matrix.
#' @param psi Numeric vector of current feature-specific residual variances.
#' @param Ex Matrix of conditional expectations of latent variables.
#' @param Exx_sum Numeric matrix giving the summed conditional second moments of
#'   the latent variables.
#' @param alpha Elastic-net mixing parameter. `alpha = 1` gives an L1 penalty
#'   and `alpha = 0` gives an L2 penalty.
#' @param lambda Nonnegative regularization parameter.
#' @param maxiter Integer maximum number of coordinate-descent iterations.
#'   Default is `100`.
#' @param tol Numeric convergence tolerance. Default is `1e-5`.
#'
#' @return A list with components:
#' \describe{
#'   \item{W}{Updated latent-factor loading matrix.}
#'   \item{psi}{Updated feature-specific residual variances.}
#'   \item{iter}{Number of coordinate-descent iterations performed.}
#'   \item{converge}{Logical indicating whether convergence was achieved.}
#' }
#'
#' @details
#' The function minimizes the elastic-net penalized expected loss for the latent
#' factor model:
#' \deqn{
#'   Y = X W Z^\top + E.
#' }
#'
#' Each element of `W` is updated using a soft-thresholded coordinate-descent
#' step. After each full sweep through the entries of `W`, the residual variances
#' `psi` are updated from the reconstruction residuals, and convergence is
#' assessed using both the relative change in penalized loss and the maximum
#' absolute change in `W`.
#'
#' @seealso
#' \code{\link{soft}}, \code{\link{mlm_factor_loss}},
#' \code{\link{mlm_factor_penfunc}}
update_W <- function(Y, Zt, W, psi, Ex, Exx_sum, alpha, lambda, maxiter = 100, tol = 1e-5) {

  ## elastic net penalized neg expected lnlike for
  ## factor model for Y = XWZ' + E
  ## Y is n x p
  ## X is latent n x nlatent
  ## W is loadings nlatent x m
  ## Z' is design matrix for cols of Y, m x p

  p <- ncol(Y)
  n <- nrow(Y)
  m <- nrow(Zt)

  nlatent <- nrow(W)
  converge <- FALSE
  W_old <- W


  loss <-  mlm_factor_loss(Y, Zt, W, psi, Ex, Exx_sum, alpha, lambda)
  pen_curr <- loss  + lambda*alpha * sum(abs(W)) + lambda*(1-alpha) * sum(W^2)/2
  pen_old <- pen_curr

  for (iter in 1:maxiter) {

    ## update terms that depend on psi and W
    Ys <- Y /sqrt(psi)
    Zts <- Zt / sqrt(psi)


    XY <- t(Ex) %*% Ys
    ZZ <- Zts %*% t(Zts)
    XYZ <- XY %*% t(Zts)
    M <- Exx_sum %*% W %*% ZZ


    ## update W
    for (irow in 1:nlatent) {
      for (jcol in 1:m) {

        W_ij_old <- W[irow, jcol]

        numer <- (XYZ[irow, jcol] + Exx_sum[irow, irow] * ZZ[jcol,jcol] * W[irow, jcol] - M[irow, jcol])/(n)
        denom <-(Exx_sum[irow, irow] * ZZ[jcol, jcol])/(n) + lambda * (1 - alpha)
        # update using soft-thresholding
        W_ij_new <- soft(numer, alpha * lambda) / denom
        W[irow, jcol] <- W_ij_new
        delta <- W_ij_new - W_ij_old
        ## rank-1 update M after update beta
        M <- M + delta*(Exx_sum[,irow] %o% ZZ[jcol,])

      }
    }


    ## update psi (note use Y and Z,not Ys or Zs)
    WZ <- W %*% Zt
    RR <- t(Y) %*% Y - 2*t(Y) %*% Ex %*% WZ + t(WZ) %*% Exx_sum %*% WZ
    psi <- diag(RR)/n


    ## update loss

    loss <-  mlm_factor_loss(Y, Zt, W, psi, Ex, Exx_sum, alpha, lambda)
    pen_curr <- loss  + lambda*alpha * sum(abs(W)) + lambda*(1-alpha) * sum(W^2)/2


    loss_rel_change <- abs(pen_curr - pen_old) / (abs(pen_old) + tol)
    W_change <- max(abs(W - W_old))

    if(is.na(loss_rel_change)){
      loss_rel_change <- 1e-8
      cat("iter loss_rel_change is NA\n")
    }
    if(is.na(W_change)){
      W[is.na(W)] <- 0
      cat("iter W_change is NA\n")
    }


    if (loss_rel_change < tol || W_change < tol) {
      converge <- TRUE
      break
    }


    W_old <- W
    pen_old <- pen_curr

  }

  return(list(W=W, psi=psi, iter=iter, converge=converge))
}

#' Soft-thresholding operator
#'
#' Applies the soft-thresholding operation commonly used in L1-penalized
#' optimization methods such as the lasso and elastic net.
#'
#' @param x Numeric value to be thresholded.
#' @param gamma Nonnegative thresholding parameter.
#'
#' @return A numeric value after soft-thresholding:
#' \deqn{
#'   S(x, \gamma) =
#'   \mathrm{sign}(x)\max(|x| - \gamma, 0).
#' }
#'
#' @details
#' The soft-thresholding operator shrinks values toward zero:
#' \itemize{
#'   \item If `x > gamma`, returns `x - gamma`.
#'   \item If `x < -gamma`, returns `x + gamma`.
#'   \item Otherwise, returns `0`.
#' }
#'
#' Missing values in `x` are replaced by zero before thresholding.
#'
#' This operator is commonly used in coordinate-descent algorithms for
#' L1-penalized regression and latent factor models.
#'
#' @seealso
#' \code{\link{update_W}}, \code{\link{mlm_factor_ccd}}
soft <- function(x, gamma){
  if(is.na(x)) x <- 0
  if(x > gamma){
    s <- x - gamma
  } else if (x < -gamma){
    s <- x + gamma
  }else{
    s <- 0
  }
  return(s)
}

#' Compute the trace of a matrix
#'
#' Computes the trace of a square matrix, defined as the sum of its diagonal
#' elements.
#'
#' @param M Numeric matrix.
#'
#' @return A numeric scalar equal to:
#' \deqn{
#'   \mathrm{tr}(M) = \sum_i M_{ii}.
#' }
#'
#' @details
#' The trace is computed using:
#' \code{sum(diag(M))}.
#'
#' This utility function is used throughout the latent factor model code for
#' compact matrix-expression calculations.
#'
#' @seealso
#' \code{\link[base]{diag}}, \code{\link[base]{sum}}
trace_mat <- function(M){
  return(sum(diag(M)))
}

#' Degrees of Freedom for a Factor-Structured MLM Parameterization
#'
#' Computes the effective number of free parameters (degrees of freedom)
#' for a factor-structured multivariate linear model (MLM) defined by a
#' loading matrix \code{W} and scaling parameters \code{psi}.
#'
#' The parameter count is adjusted to account for the non-identifiability
#' of latent-factor parameterizations. Specifically, the product
#' \eqn{X W Z^T} is invariant under invertible transformations
#' \eqn{X \rightarrow X A},
#' \eqn{W \rightarrow A^{-1} W B^{-1}}, and
#' \eqn{Z \rightarrow Z B},
#' where \eqn{A} and \eqn{B} are invertible matrices. This invariance
#' reduces the number of independent parameters by
#' \eqn{n_{\mathrm{latent}}(n_{\mathrm{latent}} - 1)/2}.
#'
#' The degrees of freedom are computed as
#' \deqn{
#' df = N_W + N_{\psi} -
#' \frac{n_{\mathrm{latent}}(n_{\mathrm{latent}} - 1)}{2},
#' }
#' where \eqn{N_W} is the number of nonzero elements in \code{W} and
#' \eqn{N_{\psi}} is the number of elements of \code{psi} that differ
#' from 1. Values with absolute magnitude less than \code{1e-6} are
#' treated as zero.
#'
#' @param W Numeric matrix of factor loadings with dimensions
#'   \eqn{n_{\mathrm{latent}} \times n_z}.
#' @param psi Numeric vector, matrix, or array of scaling parameters.
#'   Elements equal to 1 are treated as fixed and do not contribute to
#'   the parameter count.
#'
#' @return A non-negative integer giving the effective number of free
#'   parameters.
#'
#' @details
#' The function counts:
#' \itemize{
#'   \item nonzero entries of \code{W}, and
#'   \item entries of \code{psi} that are not equal to 1.
#' }
#' Missing values (\code{NA}) are ignored in both counts.
#'
#' If the adjusted parameter count is negative, the function returns 0.
#'
#' @seealso
#' \code{\link{mlm_factor_ccd}}, \code{\link{slate}}
mlm_factor_df <- function(W, psi){

  ## since XWZt is invariant under transformations
  ##  X -> XA, W -> A^-1 W B^-1, Z -> ZB
  ## for invertible matrices A[ntrait x ntrait]
  ## and B[nz x nz], these constraints allow
  ## subtraction of nlatent*(nlatent-1)/2 from
  ## total number of parameters
  ## total number of parameters = non-zero-W + non-one-psi
  nlatent <- nrow(W)
  nz <- ncol(W)

  eps <- 1e-6
  
  count_W <-   sum(abs(W) >= eps, na.rm=TRUE)
  count_psi <- sum(abs(psi - 1) >= eps, na.rm=TRUE)
  df <- count_W + count_psi - nlatent*(nlatent-1)/2
  
  if(df < 0) df <- 0

  return(df)
}

#' Compute conditional moments of latent variables
#'
#' Computes the conditional expectation of latent variables and the summed
#' conditional second moments for the latent factor model.
#'
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param Zt Numeric matrix of feature annotations or design covariates,
#'   corresponding to \eqn{Z^\top} in the model.
#' @param W Numeric latent-factor loading matrix.
#' @param psi Numeric vector of feature-specific residual variances.
#'
#' @return A list with components:
#' \describe{
#'   \item{Ex}{Matrix of conditional expectations of latent variables.}
#'   \item{Exx_sum}{Summed conditional second-moment matrix
#'     \eqn{\sum_i E[X_i X_i^\top \mid Y_i]}.}
#' }
#'
#' @details
#' The latent factor model is
#' \deqn{
#'   Y = X W Z^\top + E,
#' }
#' where `E` has diagonal covariance matrix `diag(psi)`.
#'
#' The function computes:
#' \deqn{
#'   E[X \mid Y]
#' }
#' and
#' \deqn{
#'   \sum_i E[X_i X_i^\top \mid Y_i].
#' }
#'
#' Missing values in `Ex` and `Exx_sum` are replaced by zero.
#'
#' @seealso
#' \code{\link{mlm_factor_em}}, \code{\link{mlm_factor_lnlike}},
#' \code{\link[MASS]{ginv}}
#' @export
mlm_Ex <- function(Y, Zt, W, psi) {
  ## sum of expected values of X for
  ## factor model for Y = XWZ' + E
  ## Y is n x p
  ## X is latent n x l
  ## W is loadings l x m
  ## Z' is design matrix for cols of Y, m x p

  nlat <- nrow(W)

  imat <- diag(nlat)
  psi_inv = diag(1.0 / psi)
  WZt <- W %*% Zt

  inv_lowdim <- ginv(imat + WZt %*% psi_inv %*% t(WZt))
  vinv <-  psi_inv - psi_inv %*% t(WZt) %*% inv_lowdim %*% WZt %*% psi_inv

  ## conditional moments of X given Y
  Ex <-   Y %*% t(WZt %*% vinv)
  Ex[is.na(Ex)] <- 0

  var_x_y <-  imat - WZt %*% vinv %*% t(WZt)

  Exx_sum <- matrix(0, nlat, nlat)
  nsubj <- nrow(Y)
  for (i in 1:nsubj) {
    mu <- Ex[i, ,drop=TRUE]
    Exx_sum <- Exx_sum + (var_x_y + mu %o% mu)
  }
  Exx_sum[is.na(Exx_sum)] <- 0
  return(list(Ex=Ex, Exx_sum=Exx_sum))
}

#' Generate multivariate Gauss-Hermite quadrature nodes and weights
#'
#' Constructs a multidimensional Gauss-Hermite quadrature grid for numerical
#' integration over Gaussian latent variables.
#'
#' @param dim Integer dimension of the multivariate integration space.
#' @param nnodes Integer number of quadrature nodes per dimension.
#'
#' @return A list with components:
#' \describe{
#'   \item{node}{Matrix of multidimensional quadrature nodes. Each row
#'     corresponds to one quadrature point.}
#'   \item{wt}{Numeric vector of quadrature weights corresponding to the rows
#'     of `node`.}
#' }
#'
#' @details
#' The function first computes one-dimensional Gauss-Hermite quadrature nodes
#' and weights using \code{gauss.quad}. The weights are adjusted to integrate
#' with respect to the standard Gaussian density rather than the Hermite weight
#' function:
#' \deqn{
#'   w_i^\ast = \frac{w_i}{\exp(-x_i^2)}.
#' }
#'
#' A multidimensional integration grid is then formed using the Cartesian
#' product of the one-dimensional nodes and weights:
#' \itemize{
#'   \item `expand.grid` is used to generate all node combinations.
#'   \item Multidimensional weights are computed as products of the
#'     one-dimensional weights.
#' }
#'
#' The resulting quadrature rule is commonly used for approximating expectations
#' involving latent Gaussian variables in the SLATE framework.
#'
#' @seealso
#' \code{\link[statmod]{gauss.quad}}, \code{\link{expected_func_x}}
gauss_hermite_quad <- function(dim, nnodes){
  ## create grid of Gauss Hermite quad points for
  ## multivariate integration of dimension dim

  gh.dat <- gauss.quad(n=nnodes, kind = "hermite", alpha = 0, beta = 0)
  gh.dat$weights <- gh.dat$weights/exp(-gh.dat$nodes^2)

  node.lst <- wt.lst <- list(length=dim)
  for(i in 1:dim){
    node.lst[[i]] <- gh.dat$nodes
    wt.lst[[i]] <- gh.dat$weights
  }

  grid.node <- expand.grid(node.lst)
  grid.node <- as.matrix(data.frame(grid.node))
  grid.wt <-  apply(expand.grid(wt.lst), 1, prod)

  return(list(node=grid.node, wt=grid.wt))
}


#' Compute the expected Cox score and information matrices
#'
#' Computes the expected Cox score vector and information matrices for the SLATE
#' survival model using numerical integration over the latent variables.
#'
#' @param surv_covar Numeric matrix of observed survival covariates.
#' @param surv_event Numeric or logical vector indicating event status, where
#'   `1` or `TRUE` denotes an observed event and `0` or `FALSE` denotes censoring.
#' @param coef Numeric vector of Cox regression coefficients, with observed
#'   covariate coefficients followed by latent-variable coefficients.
#' @param cumhaz_subj Numeric vector of subject-specific cumulative baseline
#'   hazards.
#' @param risk_sets List of risk-set information. Each element should contain
#'   `subj`, the subject indices in the risk set, and `event`, the subject
#'   indices with observed events at that event time.
#' @param Y Numeric matrix of observed features, with subjects in rows and
#'   features in columns.
#' @param W Numeric latent-factor loading matrix.
#' @param psi Numeric vector of feature-specific residual variances.
#' @param Zt Numeric matrix of feature annotations or design covariates,
#'   corresponding to \eqn{Z^\top} in the latent factor model.
#' @param ex.dat List of expected latent-risk quantities, typically returned by
#'   \code{\link{expected_func_x}}.
#' @param gh.dat List of Gauss-Hermite quadrature information, typically with
#'   components `node` and `wt`.
#'
#' @return A list with components:
#' \describe{
#'   \item{score}{Expected Cox score vector.}
#'   \item{info}{Expected Cox information matrix.}
#'   \item{info_louis}{Louis-adjusted observed information matrix, computed as
#'     `info - var_score`.}
#'   \item{var_score}{Estimated variance of the conditional score.}
#' }
#'
#' @details
#' The function separates `coef` into coefficients for observed survival
#' covariates and coefficients for latent variables. For each risk set, it
#' computes risk-set sums by integrating over the latent-variable distribution
#' conditional on the observed features and survival data.
#'
#' Numerical integration is performed using the Gauss-Hermite quadrature nodes
#' and weights in `gh.dat`. The conditional density weights combine:
#' \itemize{
#'   \item the Gaussian conditional distribution of the latent variables given
#'     observed features, and
#'   \item the subject-specific Cox survival likelihood contribution.
#' }
#'
#' The expected information matrix is accumulated from risk-set weighted second
#' moments. The event-specific score and conditional score variance are then used
#' to compute the Louis-adjusted information matrix.
#'
#' @seealso
#' \code{\link{update_cox_coef}}, \code{\link{expected_func_x}},
#' \code{\link{gauss_hermite_quad}}, \code{\link[MASS]{ginv}}
expected_cox_score <- function(surv_covar, surv_event, coef, cumhaz_subj, risk_sets,
                               Y, W, psi, Zt, ex.dat, gh.dat){
  ncov <- ncol(surv_covar)
  ncoef <- length(coef)
  coef_surv <- coef[1:ncov]
  coef_latent <- coef[(ncov+1):ncoef]
  node <- gh.dat$node
  nnode <- nrow(node)
  wt <- gh.dat$wt
  nsubj <- nrow(Y)
  nlatent <- nrow(W)

  xb <- surv_covar %*% coef_surv
  f1 <- f2 <- rep(0, nrow(node))

  ## conditional moments of X given Y, surv_time, surv_event and current
  ## parameter estimates

  # conditional moments of X given Y
  imat <- diag(nlatent)
  psi_inv = diag(1.0 / psi)
  WZt <- W %*% Zt
  inv_lowdim <- ginv(imat + WZt %*% psi_inv %*% t(WZt))
  vinv <-  psi_inv - psi_inv %*% t(WZt) %*% inv_lowdim %*% WZt %*% psi_inv
  Ex_y <-   Y %*% t(WZt %*% vinv)
  var_x_y <-  imat - WZt %*% vinv %*% t(WZt)
  vinv_x_y <- ginv(var_x_y)


  xb <- surv_covar %*% coef_surv
  f1 <- f2 <- rep(0, nrow(node))

  nadj <- length(coef_surv)

  rrisk_adj <- exp(surv_covar %*% coef_surv);
  rrisk_tot <- rrisk_adj * ex.dat$Exb
  ncovar <- ncol(surv_covar)


  coef <- c(coef_surv, coef_latent)
  S0 <- 0
  S1 <- score <- rep(0, length(coef))
  S2 <- info <- var_score <- matrix(0, length(coef), length(coef))

  for(iset in 1:length(risk_sets)){
    S0 <- S0 * 0
    S1 <- S1 * 0
    S2 <- S2 * 0

    indices <- risk_sets[[iset]]$subj

    for(isubj in indices){
      f1 <-  -cumhaz_subj[isubj] * exp(xb[isubj]  + node %*% coef_latent)
      if (surv_event[isubj] == 1) {
        f1  <- f1 + (node %*% coef_latent);
      }
      f2 <- f2*0
      for(j in 1:nnode){
        delta <- node[j,] - Ex_y[isubj]
        f2[j] = -.5* delta %*% vinv_x_y %*% delta
      }

      ## weight for conditional density of X
      wtf = exp((log(wt) + f1 + f2));
      wtf[is.na(wtf)] <- 0
      ## densum to standardize
      densum = sum(wtf);
      if(is.na(densum)){
        densum <- 1e-8
        cat("iter densum is NA\n")
      }


      if (abs(densum) < 1e-08) {
        wtf = wtf + 1e-08;
        densum = sum(wtf);
      }
      wtf <- wtf / densum

      x_surv <- surv_covar[isubj,]
      xmat <- matrix(x_surv, nrow = nnode, ncol = length(x_surv), byrow = TRUE)
      xmat <- cbind(xmat, node)
      eta <- xmat %*% coef
      eta_wt <- as.vector(exp(eta)*wtf)

      S0 <-  S0  + sum(eta_wt)
      S1 <- S1 + apply(xmat*eta_wt, 2, sum)
      S2 <- S2 + t(xmat*eta_wt) %*% xmat
    }

    mu <- (S1/S0)
    nevents <- length(risk_sets[[iset]]$event)
    info <- info + nevents * (S2/S0  - (S1/S0) %o% (S1/S0))

    event_index <- risk_sets[[iset]]$event

    for(isubj in event_index){
      f1 <-  -cumhaz_subj[isubj] * exp(xb[isubj]  + node %*% coef_latent)
      if (surv_event[isubj] == 1) {
        f1  <- f1 + (node %*% coef_latent);
      }
      f2 <- f2*0
      for(j in 1:nnode){
        delta <- node[j,] - Ex_y[isubj]
        f2[j] = -.5* delta %*% vinv_x_y %*% delta
      }

      ## weight for conditional density of X
      wtf = exp((log(wt) + f1 + f2));
      wtf[is.na(wtf)] <- 0

      ## densum to standardize
      densum = sum(wtf);
      if(is.na(densum)){
        densum <- 1e-8
        cat("iter densum is NA\n")
      }
      if (abs(densum) < 1e-08) {
        wtf = wtf + 1e-08;
        densum = sum(wtf);
      }
      wtf <- as.vector(wtf / densum)
      wtf[is.na(wtf)] <- 0

      x_surv <- surv_covar[isubj,]
      xmat <- matrix(x_surv, nrow = nnode, ncol = length(x_surv), byrow = TRUE)
      xmat <- cbind(xmat, node)
      delta <- t(t(xmat) - mu)
      score_subj <- apply(delta * wtf, 2, sum)
      score <- score + score_subj
      cross_prod <-  t(delta*wtf) %*% delta
      var_score <- var_score  +   cross_prod  - score_subj %o% score_subj
    }
  }

  info_louis <- info - var_score

  return(list(score=score, info=info, info_louis=info_louis,
              var_score=var_score))
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
  cat("Computing initial values\n")
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


#' Simulate Small Test Data to Demo slate
#' @export
slate_sim_test_data <- function(){
  
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
  
  
  n <- 25
  p <- 100
  Y <- matrix(rnorm(n*p), nrow=n)
  z_index <- rep(1:10, rep(10,10))
  Z <- matrix(0, nrow=p, ncol=10)
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


