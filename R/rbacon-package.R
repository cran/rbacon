#' rbacon
#'
#' Bacon produces Bayesian age-depth models from dated deposits, reconstructing Bayesian 
#' accumulation histories through combining radiocarbon and other dates with prior information (Blaauw and Christen, 2011).
#'
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk> J. Andres Christen <jac@cimat.mx> 
#' @importFrom grDevices dev.copy2pdf dev.cur dev.interactive dev.list dev.off extendrange grey pdf rgb
#' @importFrom graphics abline axis box curve hist image layout legend lines mtext par plot points polygon rect segments text locator
#' @importFrom stats approx coef dbeta density dgamma dnorm dunif lm median quantile rnorm weighted.mean  
#' @importFrom utils packageName read.csv read.table setTxtProgressBar txtProgressBar write.table
#' @importFrom Rcpp evalCpp
#' @importFrom coda as.mcmc gelman.diag mcmc.list
#' @importFrom data.table fread fwrite
#' @importFrom rintcal ccurve draw.dates
#' @useDynLib rbacon, .registration=TRUE
#' @name rbacon
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
