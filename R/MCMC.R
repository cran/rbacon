
#################### functions for post-run checks and adaptations ####################

#' @name scissors
#' @title Remove the first n iterations.
#' @description Removes iterations of the MCMC time series, and then updates the output file.
#' @details Bacon will perform millions of MCMC iterations for each age-model run by default, although only a fraction
#' of these will be stored. In most cases the remaining MCMC iterations will be well mixed (the upper left panel
#' of the fit of the iterations shows no undesirable features such as trends or sudden systematic drops or rises).
#' If the run has a visible remaining burn-in, scissors can be used to remove them.
#' To remove, e.g., the first 300 iterations, type \code{scissors(300)}. To remove the last 300 iterations, type \code{scissors(-300)}. To remove iterations 300 to 600, type \code{scissors(300:600)}.
#'
#' @param burnin Number of iterations to remove  of the iterative time series. If this value is higher than the amount of remaining iterations,
#' a warning is given and the iterations are not removed. If the provided number is negative, the iterations will be removed from the end of the run, not from the start. If a range is given, this range of iterations is removed.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   scissors(100)
#'   agedepth(d.res=50, age.res=50, d.by=10)
#' }
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   scissors(100)
#'   agedepth()
#' }
#'
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#' \url{https://projecteuclid.org/euclid.ba/1339616472}
#' @export
scissors <- function(burnin, set=get('info')) {
  output <- read.table(paste(set$prefix, ".out", sep=""))
  if(length(burnin) > 1) {
    if(length(burnin) >= nrow(output))
      stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
	output <- output[-burnin,]
  } else {
      if(abs(burnin) >= nrow(output))
        stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
	  if(burnin > 0)
	    output <- output[-(1:burnin),] else
          output <- output[-((nrow(output)-abs(burnin)):nrow(output)),]
	  }

  write.table(output, paste(set$prefix, ".out", sep=""), col.names=FALSE, row.names=FALSE)

  info <- get('info')
  info$output <- output
  .assign_to_global ("info", info)
}



#' @name thinner
#' @title Thin iterations.
#' @description Randomly thin iterations by a given proportion, for example if autocorrelation is visible within the MCMC series.
#' @details From all iterations, a proportion is removed with to-be-removed iterations sampled randomly among all iterations.
#' @param proportion Proportion of iterations to remove. Should be between 0 and 1. Default \code{proportion=0.1}.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   thinner(.1)
#'   agedepth(d.res=50, age.res=50, d.by=10)
#' }
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   thinner(.2)
#'   agedepth()
#' }
#'
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#' \url{https://projecteuclid.org/euclid.ba/1339616472}
#' @export
thinner <- function(proportion=0.1, set=get('info')) {
  output <- read.table(paste(set$prefix, ".out", sep=""))
  if(proportion >= 1)
    stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
  proportion <- sample(nrow(output), proportion*nrow(output))
  output <- output[-proportion,]
  write.table(output, paste(set$prefix, ".out", sep=""), col.names=FALSE, row.names=FALSE)

  info <- get('info')
  info$output <- output
  .assign_to_global ("info", info)
}



#' @name Baconvergence
#' @title Test to identify poorly mixed MCMC runs.
#' @description Test how well-mixed and converged the MCMC runs are with the chosen core and settings, by running the core several times and comparing the different runs using the Gelman and Rubin Reduction factor (Brooks and Gelman, 1998).
#' @details Generally Bacon will perform millions of MCMC iterations for each age-model run, although only a fraction
#' of these will be stored. In most cases the remaining MCMC iterations will be well mixed (the upper left panel
#' of the fit of the iterations shows no strange features such as sudden systematic drops or rises).
#'  However if the iterations seem not well mixed, or if too few remain (say less than a few hundred),
#'  then you could check the Gelman and Rubin Reduction Factor. Too high differences (high Factors) between runs
#' indicate poor MCMC mixing. Robust MCMC mixing is indicated by a Gelman and Rubin Reduction factor
#' (Brooks and Gelman, 1998) below the 1.05 safety threshold.
#' @param core Name of the core, given using quotes. Defaults to one of the cores provided with rbacon, \code{core="MSB2K"}.
#' @param runs Amount of runs to test for mixing. Default \code{runs=5}.
#' @param suggest If initial analysis of the data indicates abnormally slow or fast accumulation rates, Bacon will suggest to change the prior.
#' @param verbose Provide feedback on what is happening (default \code{verbose=TRUE}).
#' @param ... additional options that can be given to the Bacon function.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#'   \donttest{
#'     Baconvergence(runs=2, ssize=100, coredir=tempfile()) # a quick-and-dirty toy example
#'   }
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#' \url{https://projecteuclid.org/euclid.ba/1339616472}
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring
#' convergence of iterative simulations.
#' _Journal of Computational and Graphical Statistics, *7*, 434-455.
#' @export
Baconvergence <- function(core="MSB2K", runs=5, suggest=FALSE, verbose=TRUE, ...) {
  MCMC <- list()
  for(i in 1:runs) { # now the other runs
    message("run number", i, "...\n")
    Bacon(core=core, suggest=suggest, run=TRUE, ask=FALSE, ...)
    set <- get('info')
    if(i == 1)
      nm <- set$prefix
    MCMC[[i]] <- read.table(paste0(nm, ".out"))
    Bacon.cleanup()
  }

  lmcmc <- c() # find the shortest run
  for(i in 1:runs)
    lmcmc <- min(lmcmc, nrow(MCMC[[i]]))
  for(i in 1:runs)
    MCMC[[i]] <- MCMC[[i]][1:lmcmc,]

  dims <- ncol(MCMC[[1]])
  rng <- c()
  for(i in 1:runs)
    rng <- range(rng, MCMC[[i]][dims])
  layout(1)
  plot(MCMC[[1]][[dims]], type="l", bty="n", xlab="", ylab="", main="", ylim=rng)
  for(i in 2:runs)
    lines(MCMC[[i]][[dims]], col=i)

  rt <- gelman.diag(mcmc.list(lapply(MCMC, as.mcmc)), autoburnin=FALSE, transform=TRUE, confidence=0.97)
  if(verbose) {
    message("Did ", runs, " Bacon runs.")
    message("Gelman and Rubin Reduction Factor ", rt$mpsrf, " (smaller and closer to 1 is better).")
    if(rt$mpsrf > 1.05)
      message("Probably not a robust MCMC run! Too much difference between runs, above the 1.05 threshold. Increase sample size?\n") else
        message("Robust MCMC mixing, below the 1.05 safety threshold.\n")
  }
}



# calculate the proportion of dates that are within the age-depth model's confidence ranges
overlap <- function(set=get('info'), digits=0, verbose=TRUE) {
  d <- set$dets[,4]
  inside <- rep(1, length(d))
  for(i in 1:length(d)) {
    daterng <- set$calib$probs[[i]]
    daterng <- cbind(cumsum(daterng[,2])/sum(daterng[,2]), daterng[,1])
    daterng <- approx(daterng[,1], daterng[,2], c((1-set$prob)/2, 1-(1-set$prob)/2))$y
    age <- quantile(Bacon.Age.d(d[i], BCAD=FALSE), c((1-set$prob)/2, 1-(1-set$prob)/2))
    daterng <- daterng[!is.na(daterng)]
    if(length(daterng) > 0)
      if(max(daterng) < min(age) || max(age) < min(daterng))
        inside[i] <- 0
  }
  inside <- 100*sum(inside)/length(d)
  if(verbose) 
    message(if(inside < 80) "Warning! Only ", round(inside, digits), "% of the dates overlap with the age-depth model (", 100*set$prob, "% ranges)")
}
