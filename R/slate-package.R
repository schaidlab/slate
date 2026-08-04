## Created: 8/3/2026
## Author: Dan Schaid

#' Sparse Latent factor model for Annotated features and Time-to-Event data (SLATE)
#'
#' The package download, NEWS, and README are available on github and CRAN: \url{https://cran.r-project.org/package=slate}
#'
#' @section Functions:
#'
#' Below are listed some of the functions:
#'
#' \code{\link{slate}}: 
#'
#' \code{\link{slate_sim_test_data}}: simulate test data
#'
#' \code{\link{slate_grid}}:
#' 
#' \code{\link{slate_best_fit}}:
#'
#' @examples
#' library(slate)
#'
#' @importFrom survival coxph coxph.detail
#' @importFrom statmod gauss.quad
#' @importFrom MASS ginv
#' @importFrom stats rnorm runif rexp
#' @importFrom stats pnorm var
"_PACKAGE"


## Commands to build the package using devtools

##install.packages(c("rlang","devtools", "roxygen2", "usethis", "testthat"))
#remove.packages("slate")

#devtools::document()   # creates man/*.Rd and updates NAMESPACE

#devtools::check()      # checks package

##devtools::check_win_devel()

#devtools::build()      # creates slate_*.tar.gz

##usethis::use_vignette("getting-started")

###### steps from arsenal ##########
####################################
# devtools::install()    # installs locally
# devtools::check_man()
# devtools::test()
# devtools::check()
# withr::with_libpaths(c("../testinstalls/", .libPaths()), devtools::install(build_vignettes = TRUE, dependencies = FALSE))
# devtools::build("../slate/")

#### to upload to CRAN
## Update DESCRIPTION, README.md, NEWS.md, and cran-comments.md
# devtools::revdep_check()
# devtools::release()
