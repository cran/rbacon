# estimate how many MCMC iterations will be ran and returned
bacon.its <- function(ssize, burnin, set=get('info'), ACCEP_EV=20, EVERY_MULT=25, BURN_IN_MULT=3000) {
#   int it = EVERY_MULT * All.Dim() * (BURN_IN_MULT + ssize); bacon.cpp line 109
#   int every=  EVERY_MULT*All.Dim(); line 112
#  Rprintf("bacon: burn in (initial iterations which will be removed): %d\n", All.Dim() * EVERY_MULT * BURN_IN_MULT);  line 174

  dims <- set$K + 2 # accrates, start age, accumulation rate, memory
  store.every <- dims * EVERY_MULT # depends on the amount of parameters
  MCMC.size <- store.every * (ssize + burnin + BURN_IN_MULT) # all iterations
  # MCMC.kept <- MCMC.size - (store.every * BURN_IN_MULT) # removing burnin
  message(" Will run ", prettyNum(MCMC.size, big.mark=","), " iterations and store ", prettyNum(ssize, big.mark=","))
}



#################### functions for post-run checks and adaptations ####################

#' @name scissors
#' @title remove the first n iterations
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
#' @param write Whether or not to write the changes to the output file. Defaults to TRUE.
#' @param save.info By default, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   nrow(info$output)
#'   scissors(100)
#'   nrow(info$output)
#' }
#'
#' @export
scissors <- function(burnin, set=get('info'), write=TRUE, save.info=set$save.info) {
  output <- fastread(paste0(set$prefix, ".out"))
  if(set$isplum)
    plumout <- fastread(paste0(set$prefix, "_plum.out"))
  if(length(burnin) > 1) {
    if(length(burnin) >= nrow(output))
      stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
    output <- output[-burnin,]
    if(set$isplum)
      plumout <- plumout[-burnin,]
  } else {
      if(abs(burnin) >= nrow(output))
        stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
      if(burnin > 0) {
        output <- output[-(1:burnin),]
        if(set$isplum)
          plumout <- plumout[-(1:burnin),]
      } else {
          output <- output[-((nrow(output)-abs(burnin)):nrow(output)),]
          if(set$isplum)
            plumout <- plumout[-((nrow(plumout)-abs(burnin)):nrow(plumout)),]
        }
    }

  if(write)
    fastwrite(output, paste0(set$prefix, ".out"), col.names=FALSE, row.names=FALSE)
  if(set$isplum) {
    if(write)
      fastwrite(plumout, paste0(set$prefix, "_plum.out"), col.names=FALSE, row.names=FALSE)
    set$phi <- plumout[,1]
    set$ps <- plumout[,-1] # can be >1 column
  }
  set$output <- output
  set$Tr <- nrow(output)
  set$Us <- output[,ncol(output)] # is this the correct column?
  if(save.info)
    assign_to_global("info", set)
  invisible(set)
}



#' @name thinner
#' @title thin iterations
#' @description Randomly thin iterations by a given proportion, for example if autocorrelation is visible within the MCMC series.
#' @details From all iterations, a proportion is removed with to-be-removed iterations sampled randomly among all iterations.
#' @param proportion Proportion of iterations to remove. Should be between 0 and 1. Default \code{proportion=0.1}.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param write Whether or not to write the changes to the output file. Defaults to TRUE.
#' @param save.info By default, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   nrow(info$output)
#'   thinner(.2)
#'   nrow(info$output)
#' }
#'
#' @export
thinner <- function(proportion=0.1, set=get('info'), write=TRUE, save.info=set$save.info) {
  output <- fastread(paste0(set$prefix, ".out"))
  if(set$isplum)
    plumout <- fastread(paste0(set$prefix, "_plum.out"))
  if(proportion >= 1)
    stop("cannot remove that many iterations, there would be none left!", call.=FALSE)
  proportion <- sample(nrow(output), proportion*nrow(output))
  output <- output[-proportion,]
  if(write)
    fastwrite(output, paste0(set$prefix, ".out"), col.names=FALSE, row.names=FALSE)

  #info <- get('info')
  set$output <- output
  if(set$isplum) {
    plumout <- plumout[-proportion,]
    if(write)
      fastwrite(plumout, paste0(set$prefix, "_plum.out"), col.names=FALSE, row.names=FALSE)
    set$phi <- plumout[,1]
    set$ps <- plumout[,-1] # could be >1 columns
  }
  if(save.info)
    assign_to_global("info", set)
  invisible(set)
}



#' @name Baconvergence
#' @title test to identify poorly mixed MCMC runs
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
#' @references
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring
#' convergence of iterative simulations.
#' _Journal of Computational and Graphical Statistics_, *7*, 434-455.
#' @export
Baconvergence <- function(core="MSB2K", runs=5, suggest=FALSE, verbose=TRUE, ...) {
  MCMC <- list()
  for(i in 1:runs) { # now the other runs
    message("run number", i, "...\n")
    Bacon(core=core, suggest=suggest, run=TRUE, ask=FALSE, ...)
    set <- get('info')
    if(i == 1)
      nm <- set$prefix
    MCMC[[i]] <- fastread(paste0(nm, ".out"))
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



#' @name MCMC.diagnostics
#' @title test mixing and stationarity of the MCMC run
#' @description Test how well-mixed and stationary the MCMC run is. A good value for the effective sample size ('ess', number of effective independent samples from the MCMC iterations) is >200 (>1000 indicates an excelling mixing). Besides the mixing, stationarity 'z' is also measured (the start of the run is compared with the end). A 'z' below 1.96 (1 standard deviation) indicates no drift, and if it is >2.58 (2 standard deviations) then the MCMC chain is likely drifting.
#' @details Generally Bacon will perform millions of MCMC iterations for each age-model run, although only a fraction
#' of these will be stored. In most cases the remaining MCMC iterations will be well mixed (ess, and also visually check that the upper left panel
#' of the fit of the iterations shows no strange features such as sudden systematic drops or rises).
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param ssize Number of MCMC iterations.
#' @param talk Whether or not to provide feedback. Defaults to TRUE.
#' @return The 'ess' and 'z' scores, together with an evaluation of the values.
#' @examples
#'   \donttest{
#'     Bacon(ssize=100, coredir=tempfile()) # check the reported warnings
#'   }
#' @export
MCMC.diagnostics <- function(set=get("info"), ssize=nrow(set$output), talk=TRUE) {
  energy <- coda::as.mcmc(set$Us)
  ess <- coda::effectiveSize(energy)
  z <- abs(coda::geweke.diag(energy)$z)
  
  ssize.warn <- paste0(" please run using more iterations (ssize >", ssize, ")")
  
  if(length(energy) < 500) {
    message("Warning: very short MCMC chain,", ssize.warn) 
    invisible(NA)
  } else {

    if(talk) {
      if(ess < 10)
        message("Warning, the MCMC run has a very high autocorrelation (effective sample size=", round(ess,2), ", <100. So,", ssize.warn) else
        if(ess < 100)
          message("Warning, poor MCMC mixing (effective sample size=", round(ess,2), ", <100) -" , ssize.warn) else
          if(ess < 200)
            message("MCMC mixing (effective sample size=", round(ess,2),  ", <200) could be better -", ssize.warn) else
            if(ess < 1000)
              message("Good MCMC mixing (effective sample size=", round(ess,2), ", >200)") else
                message("Excellent MCMC mixing (effective sample size=", round(ess,2), ", >1000)")

      if(z < 1.96) # <1 sd
         message("No sign of MCMC drift (z=", round(z,2), ", <1.96), OK") else
         if(z < 2.58) # <2 sd
           message("Warning, there's a hint of MCMC drift (z=", round(z,2), ", >1.96),", ssize.warn) else
             message("Warning, non-stationary MCMC (z=", round(z,2), ", >2.58),", ssize.warn)
    }

    diag <- c(ess, z)
    names(diag) <- c("effective sample size (ess)", "z")
    invisible(diag)
  }
}


# for plum post-run analysis, function no longer used
# model.Pb.overlap <- function(set=get('info'), talk=TRUE, roundby=c()) {
#
#   depths <- set$dets[,4]
#   A_overlap <- c()
#
#   for(i in 1:nrow(set$dets)) {
#     Aseq <- set$Ai$x[[i]]
#     A_modelled_probs <- set$Ai$y[[i]]
#     A_measured_probs <- dnorm(Aseq, set$dets[i,2], set$dets[i,3]) # Plum assumes a normal dist for A, not rbacon's default student-t
#     A_overlap[i] <- rice::overlap(list(cbind(Aseq, A_modelled_probs), cbind(Aseq, A_measured_probs)),
#       talk = FALSE, visualise = FALSE)
#   }
#
#   if(talk) {
#     if(length(roundby) == 0) roundby <- 2
#       min.overlap <- which(A_overlap==min(A_overlap))[1]
#       min.overlap <- cbind(round(A_overlap[min.overlap], roundby), depths[min.overlap])
#       max.overlap <- which(A_overlap==max(A_overlap))[1]
#       max.overlap <- cbind(round(A_overlap[max.overlap], roundby), depths[max.overlap])
#       mean.overlap <- round(mean(A_overlap), roundby)
#
#       message("Average overlap between measured and modelled Pb-210: ", mean.overlap, "%, from ",
#         min.overlap[1], "% at ", min.overlap[2], " ", set$unit, " to ",
#         max.overlap[1], "% at ", max.overlap[2], set$unit)
#   }
#
#   return(cbind(depths, A_overlap))
# }


# # function no longer used
# model.dates.overlap <- function(set=get('info'), talk=TRUE, roundby=c()) {
#   dates <- set$calib$probs
#   depths <- set$dets[,4]
#
#   model.ages <- lapply(depths, function(d) {
#     dens <- density(Bacon.Age.d(d))
#     list(x = dens$x, y = dens$y)
#   })
#
#   overl <- 100*mapply(function(date, modelled) {
#    # rice::coverage(cbind(modelled$x, modelled$y), date,
#    #   visualise = FALSE) # coverage is no longer a rice function
#     rice::overlap(list(cbind(modelled$x, modelled$y), date), visualise=FALSE)
#   }, dates, model.ages)
#
#   if(talk) {
#     if(length(roundby) == 0) roundby <- 2
#       min.overlap <- which(overl==min(overl))[1]
#       min.overlap <- cbind(round(overl[min.overlap], roundby), depths[min.overlap])
#       max.overlap <- which(overl==max(overl))[1]
#       max.overlap <- cbind(round(overl[max.overlap], roundby), depths[max.overlap])
#       mean.overlap <- round(mean(overl), roundby)
#
#       message("Average coverage (% of model covered by each date) ", mean.overlap, "%, from ",
#         min.overlap[1], "% at ", min.overlap[2], " ", set$unit, " to ",
#         max.overlap[1], "% at ", max.overlap[2], set$unit)
#   }
#
#   return(cbind(depths, overl))
# }



model.Pb.hpd <- function(set=get('info'), prob=0.95, decimals=1, verbose=TRUE) {
  depths <- set$dets[,4]
  n <- min(length(depths), length(set$Ai$x))
  dates <- 1:n
  
  A_overlap <- c()
  for(i in dates) {
    Aseq <- set$Ai$x[[i]]
    A_modelled_probs <- set$Ai$y[[i]]
    A_measured_probs <- dnorm(Aseq, set$dets[i,2], set$dets[i,3]) # Plum assumes a normal dist for A, not rbacon's default student-t
    A_overlap[i] <- rice::hpd.overlap(cbind(Aseq, A_modelled_probs), cbind(Aseq, A_measured_probs), prob=prob, add.zeros=TRUE)
  }
  
  mean_A_overlap <- 100*sum(A_overlap)/length(A_overlap)
  
  string <- paste0(round(mean_A_overlap, decimals), 
    "% (", sum(A_overlap), "/", length(A_overlap), 
    ") of the modelled and measured Pb values overlap (",
    100*set$prob, "% hpd ranges)")
  
  if(verbose)
    message(string)
  
  invisible(c(round(mean_A_overlap, decimals), string))
}



model.dates.hpd <- function(set=get('info'), prob=0.95, decimals=1, verbose=TRUE) {
  depths <- set$dets[,4]
  if(set$isplum)
    if(set$hasBaconData)
      depths <- set$detsBacon[,4]
  if(length(set$slump) > 0)
    depths <- toslump(depths, set$slump, remove=FALSE)
  dates <- 1:length(depths)

  get.modelages <- function(i) {
    depth.age <- density(Bacon.Age.d(depths[i], set), na.rm=TRUE)
    cbind(depth.age$x, depth.age$y/sum(depth.age$y))
  }

  # for each dated depth, check if any of its date's hpds fall within any of the model's hpds
  this.overlap <- function(i) 
    rice::hpd.overlap(get.modelages(i), set$calib$probs[[i]], prob=prob, add.zeros=TRUE) 

  # proportion of dates that overlap with the model - at hpd level
  inorout <- sapply(dates, this.overlap)
  frac.in <- length(which(inorout==TRUE)) / length(depths)
  
  string <- paste0(if(frac.in < .80) "Warning! Only ", round(100*frac.in, decimals), "% of the dates (", length(which(inorout==TRUE)), "/", length(depths), ") overlap with the age-depth model (", 100*set$prob, "% hpd ranges)")
  
  if(verbose) 
    message(string)
  
  invisible(c(round(frac.in, decimals+2), string))
}



# calculate the proportion of dates that are within the age-depth model's confidence ranges
overlap.intervals <- function(set=get('info'), digits=0, verbose=TRUE) {
  d <- set$dets[,4]
  top <- ifelse(length(set$d.min) == 0, 1, min(which(d >= set$d.min)))
  bottom <- ifelse(length(set$d.max) == 0, length(d), max(which(d <= set$d.max)))
  these <- top:bottom
  inside <- rep(1, length(these))
  for(i in these) {
    daterng <- set$calib$probs[[i]]
    daterng <- cbind(cumsum(daterng[,2])/sum(daterng[,2]), daterng[,1])
    daterng <- approx(daterng[,1], daterng[,2], c((1-set$prob)/2, 1-(1-set$prob)/2))$y
    age <- quantile(Bacon.Age.d(d[i], set, BCAD=FALSE), c((1-set$prob)/2, 1-(1-set$prob)/2), na.rm=TRUE)
    daterng <- daterng[!is.na(daterng)]
    if(length(daterng) > 0)
      if(max(daterng) < min(age) || max(age) < min(daterng))
        inside[i] <- 0
  }

  inside <- 100*sum(inside)/length(these)
  if(verbose) 
    message(if(inside < 80) "Warning! Only ", round(inside, digits), "% of the dates overlap with the age-depth model (", 100*set$prob, "% ranges)")
  invisible(inside)
}



learning <- function(set=get('info'), decimals=2, talk=TRUE) {
  # accumulation rate
  prioracc.mean <- set$acc.mean # can be multiple entries
  prioracc.shape <- set$acc.shape # can be multiple entries
  prioracc.sd <- sqrt(prioracc.mean^2 / prioracc.shape)

  if(is.na(set$hiatus.depths[1])) { # can also be a boundary
    postacc.mean <- set$post.acc[1]
    postacc.shape <- set$post.acc[2]
    hiatus.text <- ""
  } else {
      postacc.mean <- set$post.acc[,1]
      postacc.shape <- set$post.acc[,2]

      if(is.na(set$boundary[1])) {
        if(length(set$hiatus.depths) == 1) {
          posthiatus.mean <- set$post.hiatus[1]
          posthiatus.shape <- set$post.hiatus[2]
        } else {
            posthiatus.mean <- set$post.hiatus[,1]
            posthiatus.shape <- set$post.hiatus[,2]
          }
        priorhiatus.mean <- set$hiatus.mean
        priorhiatus.shape <- set$hiatus.shape
        priorhiatus.sd <- sqrt(priorhiatus.mean^2 / priorhiatus.shape)

        priorhiatus.precision <- 1 / (priorhiatus.mean^2 / priorhiatus.shape^2)
        posthiatus.precision <- 1 / (posthiatus.mean^2 / posthiatus.shape^2)
        hiatus.learned <- sqrt(posthiatus.precision / priorhiatus.precision)
        hiatus.z <- (posthiatus.mean - set$hiatus.mean) / priorhiatus.sd

        hiatus.text <- paste0("Hiatus learning ratio: ",
          paste(round(hiatus.learned, decimals), collapse = " & "), 
          "; z-difference: ",
          paste(round(hiatus.z, decimals), collapse = " & "))
      }
  }
  
  prioracc.precision <- 1 / (prioracc.mean^2 / prioracc.shape)
  postacc.precision <- 1 / (postacc.mean^2 / postacc.shape)
  acc.learned <- sqrt(postacc.precision / prioracc.precision)
  acc.z <- (postacc.mean - prioracc.mean) / prioracc.sd
  
  acc.text <- paste0("Accumulation learning ratio: ",
    paste(round(acc.learned, decimals), collapse = " & "), 
    "; z-difference: ",
    paste(round(acc.z, decimals), collapse = " & "))
    
  # memory
  priormem.mean <- set$mem.mean
  priormem.strength <- set$mem.strength
  priormem.sd <- sqrt(abs(priormem.mean * (1 - priormem.mean)) / (priormem.strength + 1))
  priormem.precision <- 1 / priormem.sd
  postmem.mean <- set$post.mem[1]
  postmem.strength <- set$post.mem[2]
  postmem.sd <- sqrt(abs(postmem.mean * (1 - postmem.mean)) / (postmem.strength + 1))
  postmem.precision <- 1/postmem.sd
  mem.learned <- sqrt(postmem.precision / priormem.precision)
  mem.z <- (postmem.mean - priormem.mean) / priormem.sd
  
  mem.text <- paste0("Memory learning ratio: ",
    paste(round(mem.learned, decimals), collapse = ", "), 
    "; z-difference: ",
    round(mem.z, decimals), collapse = ", ")
  
  # same for hiatus (if present)? Not with uniform/unbounded prior
  
  if(set$isplum) {
    # influx, phi
    priorphi.mean <- set$phi.mean
    priorphi.shape <- set$phi.shape
    priorphi.sd <- sqrt(priorphi.mean^2 / priorphi.shape) 
    postphi.mean <- set$post.phi[1]
    postphi.shape <- set$post.phi[2]
   
    priorphi.precision <- 1 / (priorphi.mean^2 / priorphi.shape)
    postphi.precision <- 1 / (postphi.mean^2 / postphi.shape)
    phi.learned <- sqrt(postphi.precision / priorphi.precision)
    phi.z <- (postphi.mean - priorphi.mean) / priorphi.sd
  
    phi.text <- paste0("Influx learning ratio: ",
      paste(round(phi.learned, decimals), collapse = ", "), 
      "; z-difference: ",
      round(phi.z, decimals), collapse = ", ")
  
    # supported, sup
    priorsup.mean <- set$s.mean
    priorsup.shape <- set$s.shape
    priorsup.sd <- sqrt(priorsup.mean^2 / priorsup.shape) 
    postsup.mean <- set$post.supp[1]
    postsup.shape <- set$post.supp[2]

    priorsup.precision <- 1 / (priorsup.mean^2 / priorsup.shape)
    postsup.precision <- 1 / (postsup.mean^2 / postsup.shape)
    sup.learned <- sqrt(postsup.precision / priorsup.precision)
    sup.z <- (postsup.mean - priorsup.mean) / priorsup.sd

    sup.text <- paste0("Supported learning ratio: ",
      paste(round(sup.learned, decimals), collapse = ", "), 
      "; z-difference: ",
      round(sup.z, decimals), collapse = ", ")
  }

  if(talk) {
    message(acc.text)
    message(mem.text)
    if(length(set$hiatus.depths) > 0)
      if(is.na(set$boundary[1]))
        message(hiatus.text)
    if(set$isplum) {
      message(phi.text)
      message(sup.text)
    }
  }

  txt <- c(acc.text, mem.text)  
  if(length(set$hiatus.depths[!is.na(set$hiatus.depths)]) > 0)
    if(is.na(set$boundary[1]))
      txt <- c(txt, hiatus.text)
  if(set$isplum)
    txt <- c(txt, phi.text, sup.text)

  invisible(txt)
}



# internal function to write a summary file to the core's directory. Contains information about how the run went, how well the model fits the dates, and how the prior and posterior distributions compare. 
summarise.run <- function(set=get("info"), roundby=roundby, overlap.hpds=overlap.hpds, learn) {
  post.acc <- round(set$post.acc, roundby)
  post.mem <- round(set$post.mem, roundby)
  if(set$isplum) {
    post.phi <- round(set$post.phi, roundby)
    post.supp <- round(set$post.supp, roundby)
  }
  
  string <- paste0("MCMC diagnostics:\n  effective sample size (ESS) ", 
    round(set$MCMCdiagnostics[1], roundby), "\n  integrated autocorrelation time (IAT) ", 
    round(set$ssize/set$MCMCdiagnostics[1], roundby), "\n  MCMC drift (z) ",
    round(set$MCMCdiagnostics[2], roundby), "\n\n",
    overlap.hpds,

    "\n\nAccumulation rate (", 
    set$age.unit, "/", set$depth.unit,  "):\n  prior(s): mean ",
    ifelse(length(post.acc)==1, set$acc.mean, paste(set$acc.mean, collapse = " & ")),
    ", shape ",
    ifelse(length(post.acc)==1, set$acc.shape, paste(set$acc.shape, collapse = " & ")),
    "\n  posterior(s): mean ",
    ifelse(length(post.acc)==1, post.acc[1], paste(post.acc[,1], collapse = " & ")),
    ", shape ",
        ifelse(length(post.acc)==1, post.acc[2], paste(post.acc[,2], collapse = " & ")),
    "\n  ", learn[1],

    "\n\nMemory (between 0 and 1):\n  prior: mean ", set$mem.mean, ", strength ", set$mem.strength,
    "\n  posterior: mean ", post.mem[1], ", strength ", post.mem[2],
    "\n  ", learn[2],

    if(set$isplum) {
      paste("\n\nFlux (phi):\n  prior: mean ", set$phi.mean, ", shape ", set$phi.strength, 
      "\n  posterior: mean ", post.phi[1], ", strength ", post.mem[2],
      "\n  ", learn[3],

      "\n\nSupported:\n  prior: mean ", set$s.mean, ", shape ", set$s.strength, 
      "\n  posterior: mean ", post.supp[1], ", shape ", post.supp[2],
      "\n  ", learn[4])
    },

        "\n\nDetails:\n  better MCMC runs have higher ESS (aim for >200 or >500), lower IAT, and lower z (aim for <1.96=2sd))",
        "\n  learning ratio = precision(posterior) / precision(prior), where precision = 1/sdev. Values >1 indicate that we've learned about the parameter in question\n  z-difference shows by how much the posterior mean has changed relative to the prior's mean and variance (e.g., an absolute z value >1 indicate >1 sdev difference)"
      )
    
      fl <- paste0(set$prefix, "_summary.txt")
      writeLines(string, fl)
    }
