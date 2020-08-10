#' @name Bacon.Age.d
#' @title Output all ages for a single depth.
#' @description Output all MCMC-derived age estimates for a given depth.
#' @details Obtaining an age-depth model is often only a step towards a goal, e.g., plotting a core's
#' fossil series ('proxies') against calendar time. Bacon.Age.d can be used to list all MCMC-derived age estimates for a given (single) depth, for example to calculate mean ages for a depth.
#' @param d The depth of which Bacon age estimates are to be returned. Has to be a single depth.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param its The set of MCMC iterations to be used. Defaults to the entire MCMC output, \code{its=set$output}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param remove Whether or not to remove NA values (ages within slumps)
#' @author Maarten Blaauw, J. Andres Christen
#' @return Outputs all MCMC-derived ages for a given depth.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   ages.d20 = Bacon.Age.d(20)
#'   mean(ages.d20)
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#'  \url{https://projecteuclid.org/euclid.ba/1339616472}
#' @export
Bacon.Age.d <- function(d, set=get('info'), its=set$output, BCAD=set$BCAD, remove=FALSE) {
  if(length(d) > 1)
    stop("Bacon.Age.d can handle one depth at a time only", call.=FALSE)
  if(length(its) == 0)
    stop("core's data not yet in memory. Please run agedepth first\n", call.=FALSE)

  hiatus.depths <- set$hiatus.depths
  if(length(set$slump) > 0) {
   d <- toslump(d, set$slump, remove=remove)
   if(!is.na(hiatus.depths[1]))
     hiatus.depths <- set$slumphiatus
  }

  ages <- numeric(nrow(its)) # to sort R-base problem with c() and loops
  if(!is.na(d))
    if(d >= set$d.min) { # we cannot calculate ages of depths above the top depth
      topages <- as.vector(its[,1]) # ages for the core top
      maxd <- max(which(set$elbows <= d)) # find the relevant sections
      accs <- as.matrix(its[,1+(1:maxd)]) # the accumulation rates xi for each relevant section
      cumaccs <- cbind(0, t(apply(accs, 1, cumsum))) # cumulative accumulation
      ages <- topages + (set$thick * cumaccs[,maxd]) + # topages + xi * dC + ...
	   ((d-set$elbows[maxd]) * accs[,maxd]) # ... remaining bit of lowest section

      # now using a uniform jump, not gamma
      if(!is.na(hiatus.depths[1]))
        for(i in 1:length(hiatus.depths)) {
          above <- max(which(set$elbows < hiatus.depths[i]), 1)[1]
          below <- above + 1
          if(d > set$elbows[above] && d <= set$elbows[below]) { # adapt ages for sections with hiatus
            if(d > hiatus.depths[i]) # April 2020, changed snippets below from [[i]]] to [,i]
              ages <- set$elbow.below[,i] - (set$slope.below[,i] * (set$elbows[below] - d)) else
                ages <- set$elbow.above[,i] + (set$slope.above[,i] * (d - set$elbows[above]))
            }
        }
    }
  if(BCAD)
    ages <- 1950 - ages
  return(c(ages))
}


# First get ages and slopes of c's just above and below hiatus.
# then calculate slope.below and slope.above as extrapolations from the sections below resp. above (option 1, default), no adapting of slopes (option 0),
# or as w-weighted mix of prior and accrates, resp. prior only (option 2).
# then check if these slopes work (no reversals). Those with reversals revert to the original slopes.
# check approach for boundaries. check that elbows in correct order, then set slopes to either ages.below or ages.above???
hiatus.slopes <- function(set=get('info'), hiatus.option=1) {
  elbows <- set$elbows
  hiatus.depths <- set$hiatus.depths
  if(length(set$slump) > 0)
    hiatus.depths <- set$slumphiatus

  its <- cbind(set$output)
  w <- its[,ncol(its)]^(1/set$thick)
  fillvals <- array(NA, dim=c(nrow(its), length(hiatus.depths))) # 17 Apr 2020
#  set$slope.below <- numeric(nrow(its))  # 17 Apr 2020
#  set$slope.above <- numeric(nrow(its))  # 17 Apr 2020
  set$slope.below <- fillvals  # 17 Apr 2020
  set$slope.above <- fillvals  # 17 Apr 2020
  set$elbow.below <- fillvals # 17 Apr 2020
  set$elbow.above <- fillvals # 17 Apr 2020
  set$above <- numeric(1)
  topages <- as.vector(its[,1]) # ages for the core top
  accs <- as.matrix(its[,1+(1:set$K)]) # the accumulation rates xi for each section
  cumaccs <- set$thick * cbind(0, t(apply(accs, 1, cumsum))) # cumulative accumulation
  elbow.ages <- topages + cumaccs

  for(i in 1:length(hiatus.depths)) {
    above <- max(which(elbows < hiatus.depths[i]), 1) ### is this the correct one?
      elbow.below <- elbow.ages[,above+1]
      elbow.above <- elbow.ages[,above]
      orig.slope <- accs[,above]

    if(hiatus.option == 0) { # then do nothing
      slope.below <- orig.slope
      slope.above <- orig.slope
    }
    if(hiatus.option == 1) { # then extrapolate slopes above/below section w hiatus
      slope.below <- its[,above+2]
      slope.above <- its[,above]
    }
    if(hiatus.option == 2) { # then w-weighted for below, and prior-only for above
      slope.below <- w*set$output[,above+2] + (1-w)*set$acc.mean[i+1]
      slope.above <- rep(set$acc.mean[i], nrow(its))
    }

    # now calculate the ages at the hiatus/boundary, coming from below and from above
    ages.above <- elbow.above + (slope.above * (hiatus.depths[i] - elbows[above]))
    ages.below <- elbow.below - (slope.below * (elbows[above+1] - hiatus.depths[i]))

    if(!is.na(set$boundary[1])) { # then set the boundary's elbow at ages.below for both sections
      if(length(set$slumpboundary) > 0)
        boundary <- set$slumpboundary else
          boundary <- set$boundary
      ages.boundary <- elbow.below - (slope.below * (elbows[above+1] - set$boundary[i]))
      slope.above <- (ages.boundary - elbow.above) / (set$boundary[i] - elbows[above])
      slope.below <- (ages.below - ages.boundary) / (elbows[above+1] - set$boundary[i])
    }

    # for sections with reversals, use the original slopes
    reversed <- c(which(elbow.above > ages.above),
      which(ages.above > ages.below),
      which(ages.below > elbow.below),
      which(slope.below < 0), which(slope.above < 0))
   if(length(reversed) > 0) {
     slope.above[reversed] <- orig.slope[reversed]
     slope.below[reversed] <- orig.slope[reversed]
   }

    # store the updated information
    # set$elbow.below[[i]] <- elbow.below # commented Apr 2020
    # set$elbow.above[[i]] <- elbow.above
    # set$slope.below[[i]] <- slope.below
    # set$slope.above[[i]] <- slope.above
    set$elbow.below[,i] <- elbow.below # new April 2020
    set$elbow.above[,i] <- elbow.above # new 
    set$slope.below[,i] <- slope.below # new
    set$slope.above[,i] <- slope.above # new
    
    set$above <- above
  }
  return(set)
}



#' @name Bacon.hist
#' @title Calculate age distributions of depths.
#' @description Calculate the distribution of age estimates of single or multiple depths.
#' @details Age estimates of specific depths can also be plotted.
#' @param d The depth or depths for which a histogram and age ranges should be provided. If multiple depths are given, then just the age ranges, median and means (no graphs) are provided for each depth.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param hist.lab The y-axis is labelled \code{ylab="Frequency"} by default.
#' @param calc.range Calculate ranges? Takes time so can be left out
#' @param hist.lim Limits of the y-axis.
#' @param draw Draw a plot or not. Defaults to \code{draw=TRUE}, however no plots are made if more than one depth \code{d} is provided.
#'  If \code{draw=FALSE}, then the age ranges, median and mean are given for each depth (as four columns).
#' @param prob Age ranges are given as quantiles, e.g., 2.5\% and 97.5\% for the default of 95\% confidence limits (\code{prob=0.95})).
#' @param hist.col Colour of the histogram. Default grey, \code{hist.col=grey(0.5)}.
#' @param hist.border Colour of the histogram's outline. Default dark grey, \code{hist.border=grey(0.2)}.
#' @param range.col Colour of confidence ranges. Defaults to \code{range.col="blue"}.
#' @param med.col Colour of the median. Defaults to \code{med.col="green"}.
#' @param mean.col Colour of the mean. Defaults to \code{mn.col="red"}.
#' @param verbose Provide feedback on what is happening (default \code{verbose=TRUE}).
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot with the histogram and the age ranges, median and mean, or just the age ranges, medians and means if more than one depth \code{d} is given.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   Bacon.hist(20)
#'   Bacon.hist(20:30)
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#' \url{https://projecteuclid.org/euclid.ba/1339616472}
#' @export
Bacon.hist <- function(d, set=get('info'), BCAD=set$BCAD, age.lab=c(), age.lim=c(), hist.lab="Frequency", calc.range=TRUE, hist.lim=c(), draw=TRUE, prob=set$prob, hist.col=grey(0.5), hist.border=grey(.2), range.col="blue", med.col="green", mean.col="red", verbose=TRUE) {
  outfile <- paste(set$prefix, ".out", sep="")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- .Bacon.AnaOut(outfile, set)
    .assign_to_global("set", set)
  }
  hist3 <- function(d, BCAD) {
    hsts <- list(); maxhist <- 0; minhist <- 1
    pb <- txtProgressBar(min=0, max=max(1,length(d)-1), style = 3)
    for(i in 1:length(d)) {
      if(length(d) > 1)
        setTxtProgressBar(pb, i)
      ages <- Bacon.Age.d(d[i], set, BCAD=BCAD)
      if(length(ages) > 0) {
        hst <- density(ages)
        th0 <- min(hst$x)
        th1 <- max(hst$x)
        maxhist <- max(maxhist, hst$y)
        minhist <- min(minhist, max(hst$y))
        n <- length(hst$x)
        counts <- hst$y
        ds <- d[i]
        hsts <- append(hsts, pairlist(list(d=ds, th0=th0, th1=th1, n=n, counts=counts, max=maxhist, min=minhist)))
      } else hsts$d[[i]] <- d[i]
    }
    return(hsts)
  }
  hists <- hist3(d, BCAD)
  .assign_to_global("hists", hists)

  # rng <- c()
  rng <- array(NA, dim=c(length(d), 4)) # to deal with new R which does not like to fill c() using loops
  if(calc.range)
    rng <- Bacon.rng(d, set, BCAD, prob)

  if(length(d)==1)
    if(draw==TRUE) {
      hst <- hists[[1]]
      if(length(age.lab) == 0)
        age.lab <- ifelse(BCAD, "BC/AD", "cal yr BP")
      if(length(age.lim) == 0)
        age.lim <- c(hst$th0, hst$th1)
      if(BCAD)
        age.lim <- rev(age.lim)
      if(length(hist.lim) == 0)
        hist.lim <- c(0, 1.1*hst$max)

      pol <- cbind(c(hst$th0, seq(hst$th0, hst$th1, length=hst$n), hst$th1), c(0, hst$counts, 0))
      plot(0, type="n", xlim=age.lim, ylim=hist.lim, xlab=age.lab, ylab=hist.lab, yaxs="i")
      polygon(pol, col=hist.col, border=hist.border)
      segments(rng[,1], 0, rng[,2], 0, col=range.col, lwd=3)
      points(rng[,3], 0, col=med.col, pch=20)
      points(rng[,4], 0, col=mean.col, pch=20)

      if(verbose) {
        message("mean (", mean.col, "): ", round(rng[4],1), " ", age.lab,
          ", median (", med.col, "): ",  round(rng[3],1), " ", age.lab, "\n")
        message(100*prob, "% range (", range.col, "): ", round(rng[1],1), " to ", round(rng[2],1), " ", age.lab, "\n")
      }
    }
  invisible(rng)
}


# to calculate age ranges
Bacon.rng <- function(d, set=get('info'), BCAD=set$BCAD, prob=set$prob) {
  outfile <- paste(set$prefix, ".out", sep="")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- .Bacon.AnaOut(outfile, set)
    .assign_to_global("set", set)
  }

  if(length(d) > 1)
    pb <- txtProgressBar(min=0, max=max(1, length(d)-1), style=3)
  rng <- array(NA, dim=c(length(d), 4))
  for(i in 1:length(d)) {
    ages <- Bacon.Age.d(d[i], set, BCAD=BCAD)
	  if(length(!is.na(ages)) > 0) {
      rng[i,1:3] <- quantile(ages[!is.na(ages)], c(((1-prob)/2), 1-((1-prob)/2), .5))
      rng[i,4] <- mean(ages[!is.na(ages)])
    }
    if(length(d) > 1)
      setTxtProgressBar(pb, i)
  }
  if(length(d) > 1)
    close(pb)
  return(rng)
}

#' @name agemodel.it
#' @title Extract one age-model iteration
#' @description For one MCMC iteration (it), extract the corresponding age-depth model.
#' @param it The MCMC iteration of which the age-model should be calculated.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A variable with two columns - depth and the age-depth model of a single iteration.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   lines(agemodel.it(5), col="red")
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#' \url{https://projecteuclid.org/euclid.ba/1339616472}
#' @export
agemodel.it <- function(it, set=get('info'), BCAD=set$BCAD) {
  outfile <- paste(set$prefix, ".out", sep="")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- .Bacon.AnaOut(outfile, set)
    .assign_to_global("set", set)
  }
  # does this function work in cores with slumps?
  if(length(set$hiatus.depths) > 0)
    age <- sort(c(set$d, set$hiatus.depths+.001, set$hiatus.depths))
  age <- c()
  for(i in 1:length(set$elbows))
    age[i] <- Bacon.Age.d(set$elbows[i], set, BCAD=BCAD)[it]
  cbind(set$elbows, age)
}


# calculate slumpfree depths
toslump <- function(d, slump, remove=FALSE) {
  d <- sort(d)
  slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
  slices <- c(0, slump[,2] - slump[,1])
  dfree <- d
  for(i in 1:nrow(slump)) {

    inside <- which(d <= slump[i,2]) # find depths within slump, step 1
    inside <- which(d[inside] >= slump[i,1]) # step 2
    below <- which(d >= slump[i,2]) # adapt depths below slumps

	if(length(below) > 0) # depths below slump
      dfree[below] <- dfree[below] - slices[i+1]

    if(length(inside) > 0) # depths within slump
      if(min(d) < max(slump[i,]))
	    if(remove)
          dfree[inside] <- NA else
            dfree[inside] <- slump[i,1] - sum(slices[1:i])
  }
  return(dfree)
}



# calculate original depths. Needed?
fromslump <- function(d, slump) {
  slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
  slices <- slump[,2] - slump[,1]
  dorig <- d # original depths
  for(i in 1:nrow(slump)) {
    below <- which(d > min(slump[i,]))
  if(length(below) > 0)
    dorig[below] <- dorig[below] - slices[i]
  }
  return(dorig)
}
