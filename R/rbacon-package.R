#' rbacon
#'
#' Bacon produces Bayesian age-depth models from dated deposits, reconstructing Bayesian 
#' accumulation histories through combining radiocarbon and other dates with prior information (Blaauw and Christen, 2011).
#'
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk> J. Andres Christen <jac@cimat.mx> 
#' @importFrom grDevices dev.copy2pdf dev.cur dev.interactive dev.list dev.copy dev.off extendrange grey pdf cairo_pdf rgb
#' @importFrom graphics abline axis box curve hist image layout legend lines mtext par plot points polygon rect segments text locator
#' @importFrom stats approx coef dbeta density dgamma dnorm dunif lm median quantile rnorm weighted.mean var approxfun setNames
#' @importFrom utils packageName read.csv read.table setTxtProgressBar txtProgressBar write.table
#' @importFrom Rcpp evalCpp
#' @importFrom coda as.mcmc gelman.diag mcmc.list effectiveSize geweke.diag
#' @importFrom data.table fread fwrite
#' @importFrom rice draw.dates F14CtoC14 pMCtoC14 BCADtocalBP calBPtoBCAD
#' @importFrom rintcal ccurve mix.ccurves new.ccdir
#' @useDynLib rbacon, .registration=TRUE
#' @name rbacon
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
