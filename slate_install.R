##install.packages(c("rlang","devtools", "roxygen2", "usethis", "testthat"))
remove.packages("slate")
.rs.restartR()
rm(list = ls())
setwd("C:/Users/schaid/Desktop/R_Projects/slate")
devtools::document()   # creates man/*.Rd and updates NAMESPACE

##Consider adding
##importFrom("stats", "rexp", "rnorm", "runif")
##to your NAMESPACE file.

devtools::check()      # checks package

##devtools::check_win_devel()

devtools::build()      # creates slate_*.tar.gz

##usethis::use_vignette("getting-started")

devtools::install()    # installs locally

library(slate)
