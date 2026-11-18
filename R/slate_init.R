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
