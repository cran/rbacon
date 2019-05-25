#' @name Bacon.Age.d
#' @title Output all ages for a single depth.
#' @description Output all MCMC-derived age estimates for a given depth.
#' @details Obtaining an age-depth model is often only a step towards a goal, e.g., plotting a core's 
#' fossil series ('proxies') against calendar time. Bacon.Age.d can be used to list all MCMC-derived age estimates for a given (single) depth, for example to calculate mean ages for a depth. 
#' @param d The depth of which Bacon age estimates are to be returned. Has to be a single depth.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param its The set of MCMC iterations to be used. Defaults to the entire MCMC output, \code{its=set$output}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return Outputs all MCMC-derived ages for a given depth. 
#' @examples 
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   ages.d20 = Bacon.Age.d(20)
#'   mean(ages.d20)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#'  \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
Bacon.Age.d <- function(d, set=get('info'), its=set$output, BCAD=set$BCAD) {
  if(length(d) > 1)
    stop("Bacon.Age.d can handle one depth at a time only", call.=FALSE)
  if(length(its) == 0) 
    stop("core's data not yet in memory. Please run agedepth first\n", call.=FALSE)
  
  hiatus.depths <- set$hiatus.depths
  if(length(set$slump) > 0) {
    d <- approx(set$depths, set$slumpfree, d, rule=2)$y
   if(length(set$hiatus.depths) > 0)
     hiatus.depths <- set$slumphiatus
  }   

  ages <- NA
  if(!is.na(d)) 
    if(d >= set$d.min) { # we cannot calculate ages of depths above the top depth
      topages <- cbind(its[,1]) # ages for the core top
      maxd <- max(which(set$d <= d)) # find the relevant sections
      accs <- as.matrix(its[,1+(1:maxd)]) # the accumulation rates xi for each relevant section
      cumaccs <- cbind(0, t(apply(accs, 1, cumsum))) # cumulative accumulation
      ages <- topages + (set$thick * cumaccs[,maxd]) + # topages + xi * dC + ...
	   ((d-set$d[maxd]) * accs[,maxd]) # ... remaining bit of lowest section

      # following Blaauw and Christen 2011, but now using a uniform jump, not gamma.
	  # Also, here we do attempt to model ages for depths in sections with hiatuses, but results might be unexpected.
      # ages BELOW the hiatus use a w of x_k+1 (i.e. the xi of the section below the one containing the hiatus) and the prior, while those above use prior of accrate only (memory lost)
	  # discussion 21 May 2019: Should really plot age-models of sections until the one above the hiatus, then the age-model from the elbow below the hiatus, without drawing the connecting lines between c_k-1 and C_k of a hiatus. The code puts the hiatus at the bottom of each section. Use more&smaller sections. Sections with hiatus are slopes x
      if(!is.na(hiatus.depths[1])) { # adapt ages of sections around hiatuses or boundaries
        for(i in 1:length(hiatus.depths)) {
          above <- max(which(set$d <= hiatus.depths[i]), 1) 
          below <- above + 1
          if(d > set$d[above] && d < set$d[below]) { # then adapt ages

            # ages for the section just below the hiatus, w-weighted mix of xi and prior
            if(below == 2)
              elbow.below <- topages + (set$thick * its[,2]) else
                elbow.below <- topages + (set$thick * apply(its[,2:below], 1, sum))
            w <- set$output[,ncol(set$output)-1]^(1/set$thick) # memory, penultimate column of .out
            acc <- w*its[,below+1] + (1-w)*set$acc.mean[i+1] # weighted mix of xi+1 and prior 
            ages.below <- elbow.below - (acc * (set$d[below] - d))

            # ages above the hiatus - based on accrate prior only since memory is lost
            if(above == 1)
              elbow.above <- topages else
                if(above == 2)
                  elbow.above <- topages + (set$thick * its[,2]) else
                    elbow.above <- topages + (set$thick * apply(its[,2:above], 1, sum))
            ages.above <- elbow.above + (set$acc.mean[i] * (d-set$d[above]) )
              
            # ensure no reversals among the adapted ages 
            ok <- which(ages.below >= ages.above)
            if(d <= hiatus.depths[i]) 
              ages[ok] <- ages.above[ok] else
                ages[ok] <- ages.below[ok]
          }
        }
      }
      if(BCAD)
        ages <- 1950 - ages
      return(c(ages))
    }
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
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot with the histogram and the age ranges, median and mean, or just the age ranges, medians and means if more than one depth \code{d} is given.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   Bacon.hist(20) 
#'   Bacon.hist(20:30) 
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
Bacon.hist <- function(d, set=get('info'), BCAD=set$BCAD, age.lab=c(), age.lim=c(), hist.lab="Frequency", calc.range=TRUE, hist.lim=c(), draw=TRUE, prob=set$prob, hist.col=grey(0.5), hist.border=grey(.2), range.col="blue", med.col="green", mean.col="red") {
  outfile <- paste(set$prefix, ".out", sep="")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- .Bacon.AnaOut(outfile, set)
    .assign_to_global("set", set)
  }
    
  hist3 <- function(d, BCAD) {
    hsts <- c(); maxhist <- 0; minhist <- 1
    pb <- txtProgressBar(min=0, max=max(1,length(d)-1), style = 3)
    for(i in 1:length(d)) {
      if(length(d) > 1)
        setTxtProgressBar(pb, i)
      ages <- Bacon.Age.d(d[i], set, BCAD=BCAD)
      hst <- density(ages)
      th0 <- min(hst$x)
      th1 <- max(hst$x)
      maxhist <- max(maxhist, hst$y)
      minhist <- min(minhist, max(hst$y))
      n <- length(hst$x)
      counts <- hst$y
      ds <- d[i]
      hsts <- append(hsts, pairlist(list(d=ds, th0=th0, th1=th1, n=n, counts=counts, max=maxhist, min=minhist)))
    }  
    return(hsts)
  }
  hists <- hist3(d, BCAD)
  .assign_to_global("hists", hists)
  
  rng <- c()
  if(calc.range)
    rng <- Bacon.rng(d, set, BCAD, prob)
  
  if(length(d)==1 && draw==TRUE) {  
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

    cat("\nmean (", mean.col, "): ", round(rng[4],1), " ", age.lab,
      ", median (", med.col, "): ",  round(rng[3],1), " ", age.lab, "\n", sep="")
    cat(100*prob, "% range (", range.col, "): ", round(rng[1],1), " to ", round(rng[2],1), " ", age.lab, "\n", sep="")
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
    
  pb <- txtProgressBar(min=0, max=max(1, length(d)-1), style=3)
  rng <- array(NA, dim=c(length(d), 4)) 
  for(i in 1:length(d)) {
    ages <- Bacon.Age.d(d[i], set, BCAD=BCAD)
    rng[i,1:3] <- quantile(ages, c(((1-prob)/2), 1-((1-prob)/2), .5))
    rng[i,4] <- mean(ages)
    setTxtProgressBar(pb, i)
  }  
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
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
agemodel.it <- function(it, set=get('info'), BCAD=set$BCAD) {
  outfile <- paste(set$prefix, ".out", sep="")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- .Bacon.AnaOut(outfile, set)
    .assign_to_global("set", set)
  }
  d <- set$d
  if(length(set$hiatus.depths) > 0)
    age <- sort(c(d, set$hiatus.depths+.001, set$hiatus.depths))
  age <- c()
  for(i in 1:length(d))
    age[i] <- Bacon.Age.d(d[i], set, BCAD=BCAD)[it]
  cbind(d,age)
}


excise <- function(d, slump, d.by) {
  dfree <- d
  for(i in 1:nrow(slump)) {
    dup <- which(d <= slump[i,2]) # find depths within slump, part 1
    dup <- which(d[dup] >= slump[i,1]) # part 2
    dfree[dup] <- NA # and set them to NA
    below <- which(d > slump[i,1]) # adapt depths below slumps
    dfree[below] <- dfree[below] - (slump[i,2] - slump[i,1])-d.by
  }
  return(dfree) 
}


