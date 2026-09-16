#' @name Bacon.Age.d
#' @title output all ages for a single depth
#' @description Output all MCMC-derived age estimates for a given depth.
#' @details Obtaining an age-depth model is often only a step towards a goal, e.g., plotting a core's
#' fossil series ('proxies') against calendar time. Bacon.Age.d can be used to list all MCMC-derived age estimates for a given (single) depth, for example to calculate mean ages for a depth. See also Bacon.d.Age which calculates the depths of a single age estimate.
#' @param d The depth of which Bacon age estimates are to be returned. Has to be a single depth.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param its The set of MCMC iterations to be used. Defaults to the entire MCMC output, \code{its=set$output}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param na.rm Whether or not to remove NA values (ages within slumps)
#' @author Maarten Blaauw, J. Andres Christen
#' @return Outputs all MCMC-derived ages for a given depth.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   ages.d20 = Bacon.Age.d(20)
#'   mean(ages.d20)
#' }
#' @export
Bacon.Age.d <- function(d, set=get('info'), its=set$output, BCAD=set$BCAD, na.rm=FALSE) {
  if(length(d) > 1)
    stop("Bacon.Age.d can handle one depth at a time only", call.=FALSE)
  if(length(its) == 0)
    stop("core's data not yet in memory. Please run agedepth first\n", call.=FALSE)

  nr <- nrow(its)
  if(nr==nrow(set$output))
    rows <- seq_len(nr) else
      rows <- match(its[,1], set$output[,1])

  hiatus.depths <- set$hiatus.depths
  elbows <- set$elbows
  if(length(set$slump) > 0) {
    d <- toslump(d, set$slump, remove=na.rm) # commenting this causes ranges to be calculated correctly, but not the greyscales
    #elbows <- toslump(d, set$slump, remove=na.rm)
    if(!is.na(hiatus.depths[1]))
      hiatus.depths <- set$slumphiatus 
  }

  ages <- numeric(nrow(its)) # to sort R-base problem with c() and loops
  if(!is.na(d))
    if(d >= min(elbows)) { # we cannot calculate ages of depths above the top depth; was > d.min before Feb 2021
      maxd <- max(which(set$elbows <= d)) # find the relevant sections
      acc <- as.matrix(its[,1+(1:maxd)]) # the accumulation rates xi for each relevant section
      cumacc <- cbind(0, t(apply(acc, 1, cumsum))) # cumulative accumulation
      ages <- its[,1] + (set$thick * cumacc[,maxd]) + # topages + xi * dC + ...
        ((d-elbows[maxd]) * acc[,maxd]) # ... remaining bit of lowest section

      if(!is.na(hiatus.depths[1]))
        for(i in 1:length(hiatus.depths)) {
          above <- max(which(elbows < hiatus.depths[i]), 1)[1]
          below <- above + 1
          if(d > elbows[above] && d < elbows[below]) { # adapt ages for sections with hiatus; was && d <=
            if(d > hiatus.depths[i]) 
              ages <- set$elbow.below[rows,i] - (set$slope.below[rows,i] * (elbows[below] - d)) else
                ages <- set$elbow.above[rows,i] + (set$slope.above[rows,i] * (d - elbows[above]))
            }
        }
    } else
        ages <- NA # Feb 2021

   if(na.rm)
     ages <- ages[!is.na(ages)]  

   if(BCAD)
     ages <- calBPtoBCAD(ages)

    return(unname(ages)) # unname so that only the values are returned, not also names/entry numbers
}



#' @name Bacon.d.Age
#' @title output all depths for a single age
#' @description Output all MCMC-derived depth estimates for a single given age.
#' @details Obtaining an age-depth model is often only a step towards a goal, e.g., plotting a core's
#' fossil series ('proxies') against calendar time. Bacon.d.Age can be used to list all MCMC-derived depths belonging to a given (single) age, for example to calculate mean depths belonging to a modelled depth. 
#' This function was kindly written and provided by Timon Netzel (Bonn University). See also Bacon.Age.d, which calculates the ages for a single depth.
#' @param age The age estimate for which depths are to be returned. Has to be a single age.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param its The set of MCMC iterations to be used. Defaults to the entire MCMC output, \code{its=set$output}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param na.rm Whether or not to remove NA values (ages within slumps)
#' @author Maarten Blaauw, J. Andres Christen
#' @return Outputs all MCMC-derived depths for a given age.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   ages.d20 = Bacon.Age.d(20)
#'   mean(ages.d20)
#' }
#' @export
Bacon.d.Age <- function(age, set=get("info"), BCAD=set$BCAD, its=set$output, na.rm=FALSE ){

  if(BCAD && !set$BCAD) 
    age <- BCADtocalBP(age)
  if(length(age) > 1) 
    stop("Bacon.d.Age can handle one age at a time only", call. = FALSE)
  if(length(its) == 0) 
    stop("core's data not yet in memory. Please run agedepth first\n", call. = FALSE)

  # create elbow ages
  elbows <- set$elbows
  topages <- as.vector(its[,1]) # ages for the core top
  accs <- as.matrix(its[,1+(1:set$K)]) # the accumulation rates xi for each section
  cumaccs <- set$thick * cbind(0, t(apply(accs, 1, cumsum))) # cumulative accumulation
  elbow.ages <- topages + cumaccs

  if(age <= min(elbow.ages) || age > max(elbow.ages)) 
    message(" Warning, age outside the core's age range!\n")

  # prepare for any slumps
  hiatus.depths <- set$hiatus.depths
  if(length(set$slump) > 0) {
    diff.slumps <- apply(set$slump, 1, diff, na.rm=TRUE) # vector containing the difference of each slump
    if(!is.na(hiatus.depths[1])) 
      hiatus.depths <- set$slumphiatus
  }
    
  # prepare hiatus
  hiatus.check <- rep(0, ncol(elbow.ages)) # vector containing the information for each elbow whether there is a hiatus (1) or no hiatus (0)
  if(!is.na(hiatus.depths[1])) {
    # if there is no slump
    if(length(set$slump) == 0) 
      hiatus.depths <- set$hiatus.depths
    # where are the hiatuses
    for(i in 1:length(hiatus.depths)) 
      hiatus.check[min(which(elbows > hiatus.depths[i]))] <- 1
  }

  # intersections which are below the last elbow
  these.bottom <- which(elbow.ages[, ncol(elbow.ages)] < age) 
    
  # summary of the depths corresponding to the chosen age
  depths <- rep(NA, nrow(its))
  for(i in 2:ncol(elbow.ages)) {
    # consideration of hiatuses
    if(hiatus.check[i] == 1) { 
      for(j in 1:length(hiatus.depths)) {
        above <- max(which(elbows < hiatus.depths[j]), 1)[1]
        below <- above + 1
        # two age ranges for one hiatus (below/above) with the adapted slopes (could be included in the info list)
        ages.above <- set$elbow.above[,j] + (set$slope.above[,j] * (hiatus.depths[j] - elbows[above]))
        ages.below <- set$elbow.below[,j] - (set$slope.below[,j] * (elbows[below] - hiatus.depths[j]))
        # intersections below and above
        these.above <- (elbow.ages[,above] < age) * (ages.above > age) 
        these.below <- (ages.below < age) * (elbow.ages[,below] > age) 
        # find the depths which correspond to the intersections below and above the hiatus
        if(sum(these.above, na.rm=T) > 0) {
          # slopes above hiatus
          accs.above <- (1/set$slope.above[,j])[which(these.above==1)]
          # age difference with the last age range above the hiatus
          age.diff.above <- age - elbow.ages[which(these.above==1), above]
          # depths: last depth above hiatus depth + slopes * age difference
          depths[which(these.above==1)] <- elbows[above] + accs.above * age.diff.above
        }
        if(sum(these.below, na.rm=T) > 0) {
          # slopes below hiatus
          accs.below <- (1/set$slope.below[,j])[which(these.below==1)]
          # age difference with the lower hiatus age range 
          age.diff.below  <- age - ages.below[which(these.below==1)]
          # depths: hiatus depth + slopes * age difference 
          depths[which(these.below==1)] <- hiatus.depths[j] + accs.below * age.diff.below
        }
      }
    } else {
      # consideration of the cases without hiatuses
            
      # find the correct intersections
      these <- (elbow.ages[, i - 1] < age) * (elbow.ages[, i] > age)

      # find the depths which correspond to these intersections
      if(sum(these, na.rm=T) > 0) { 
        # which slopes
        accs <- 1/set$output[which(these == 1), i]
        # difference between the age and the elbow ages found
        age.diff <- age - elbow.ages[which(these == 1), i - 1]
        # depths: elbow depth + slopes *  age difference
        depths[which(these==1)] <- elbows[i-1] + accs * age.diff
      }
            
      # find these intersections which are below the last elbow
      if(length(these.bottom) > 0 & i == ncol(elbow.ages)) {
        # slopes below last elbow
        accs.bottom <- 1/set$output[these.bottom, i]
        # difference between the age and the last elbow ages found
        age.diff.bottom <-  age - elbow.ages[these.bottom, i ]
        # depth: last elbow + slopes *  age difference 
        depths[these.bottom] <- elbows[i-1] + accs.bottom * age.diff.bottom
      }
    }  
  }
    
  # Query: remove NA entries if some entries are not NA (happens for some hiatuses).
  # In these cases the depth vector does not have the same length for different age.
  # In cases where the age entered is younger than min(elbow.ages), each entry in depth is NA and is returned as such.
  if(!all(is.na(depths)) & na.rm) 
    depths <- depths[!is.na(depths)]
  # consideration of slumps
  if(length(set$slump) > 0) {
    for(j in 1:nrow(set$slump)) {
      # which depths are below the slump
      these.below.slump <- (depths >= set$slump[j,1]) 
      below.slump <- sum(these.below.slump, na.rm = T)
      if(below.slump > 0) {
        # shift the found depths with the slump depths above
        these.shift <- which(these.below.slump==1)
        depths[these.shift] <- depths[these.shift] + diff.slumps[j]  
      }
    }
  } 
  return(depths)
}



# Find the ages t of the elbows of the section that contains the hiatus h_i, t_{k} for the elbow c_k just above the hiatus, and t_{k+1} for the elbow c_{k+1} just below it.
# Find the starting age t_s of the hiatus by taking the age of c_{k+1} and the slope x_{k+1} of the section below it, and extrapolating to the depth of hiatus_i. Do the same for the ending age by extrapolating from the slope above the section with the hiatus.
# Now that we know the length of the hiatus (if it's not a boundary), calculate the end age of the hiatus as t_e = t_s + l_i.
# Extrapolate from this end age to the top elbow's age, t_k.
# report the corresponding slopes x.
# for boundaries, the hiatus length is 0, so we find t_s and extrapolate to t_k.
hiatus.slopes <- function(set=get('info')) {
  elbows <- set$elbows
  hiatus.depths <- if(length(set$slump) > 0) 
    set$slumphiatus else set$hiatus.depths
  h <- length(hiatus.depths) # number of hiatuses
  ishiatus <- is.na(set$boundary[1])

  its <- as.matrix(set$output)
  n <- nrow(its) # number of iterations
  
  set$slope.below <- matrix(NA_real_, nrow=n, ncol=h)
  set$slope.above <- matrix(NA_real_, nrow=n, ncol=h)
  set$elbow.below <- matrix(NA_real_, nrow=n, ncol=h)
  set$elbow.above <- matrix(NA_real_, nrow=n, ncol=h)
  if(ishiatus) {
    set$hiatus.start <- matrix(NA_real_, nrow=n, ncol=h)
    set$hiatus.end <- matrix(NA_real_, nrow=n, ncol=h)
  }
  
  accs <- its[,1+(1:set$K)] # the accumulation rates (actually sedimentation times) xi for each section
  elbow.ages <- its[,1] + set$thick * cbind(0, t(apply(accs, 1, cumsum))) # cumulative accumulation at elbows, starting at 0
  hiatus.section <- findInterval(hiatus.depths, elbows) # finds which section(s) overlap(s) the hiatus depth(s)

  for(i in 1:h) {
    k <- hiatus.section[i]
    elbow.below <- elbow.ages[,k+1]
    elbow.above <- elbow.ages[,k]
    slope.below <- accs[,k+1]
    slope.above <- accs[,k-1]
    slope.orig <- accs[,k]
    age.h <- elbow.above + slope.orig * (hiatus.depths[i] - elbows[k])

    if(ishiatus) { # extrapolate using the slopes of the surrounding sections
      hiatus.start <- elbow.below - slope.below*(elbows[k+1] - hiatus.depths[i])
      hiatus.end <- elbow.above + slope.above*(hiatus.depths[i] - elbows[k])
      hiatus.duration <- hiatus.start - hiatus.end
      ok <- !is.na(hiatus.duration) & hiatus.duration >= 0
      slope.below[ok] <- (elbow.below[ok] - hiatus.start[ok]) / (elbows[k+1] - hiatus.depths[i])
      slope.above[ok] <- (hiatus.end[ok] - elbow.above[ok]) / (hiatus.depths[i] - elbows[k])
      slope.above[!ok] <- slope.orig[!ok]
      slope.below[!ok] <- slope.orig[!ok]
      hiatus.start[!ok] <- age.h[!ok] # new Sep '26
      hiatus.end[!ok] <- age.h[!ok] # new Sep '26
      hiatus.duration[!ok] <- 0
    } else { # it's a boundary, and we interpolate to and from hiatus.end
        boundary.age <- elbow.below - slope.below * (elbows[k+1] - set$boundary[i])
        ok <- boundary.age >= elbow.above & boundary.age <= elbow.below
        slope.above[ok] <- (boundary.age[ok] - elbow.above[ok]) /
          (set$boundary[i] - elbows[k])
        slope.below[ok] <- (elbow.below[ok] - boundary.age[ok]) /
          (elbows[k+1] - set$boundary[i])
        slope.above[!ok] <- slope.orig[!ok]
        slope.below[!ok] <- slope.orig[!ok]
      }

    # store the updated information
    set$elbow.below[,i] <- elbow.below
    set$elbow.above[,i] <- elbow.above
    set$slope.below[,i] <- slope.below
    set$slope.above[,i] <- slope.above
    if(ishiatus) {
      set$hiatus.start[,i] <- hiatus.start
      set$hiatus.end[,i] <- hiatus.end
    }
  }
  return(set)
}



#' @name Bacon.hist
#' @title calculate age distributions of depths
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
#' @param progress Show a progress bar (default \code{progress=TRUE}).
#' @param save.info A variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) can be saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A variable called `hists', and a plot with the histogram and the age ranges, median and mean, or just the age ranges, medians and means if more than one depth \code{d} is given.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   Bacon.hist(20)
#'   Bacon.hist(20:30)
#' }
#' @export
Bacon.hist <- function(d, set=get('info'), BCAD=set$BCAD, age.lab=c(), age.lim=c(), hist.lab="Frequency", calc.range=TRUE, hist.lim=c(), draw=TRUE, prob=set$prob, hist.col=grey(0.5), hist.border=grey(.2), range.col="blue", med.col="green", mean.col="red", verbose=TRUE, progress=TRUE, save.info=set$save.info) {
  outfile <- paste0(set$prefix, ".out")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- Bacon.AnaOut(outfile, set, MCMC.resample=FALSE)
    if(save.info)
      assign_to_global("info", set)
  }
  
  if(min(d) < min(set$elbows)) {
    these.d <- d[which(d < min(set$elbows))]
    stop("depth(s) ", paste(these.d, collapse=", "), " lie above the top/minimum depth (", min(set$elbows), " ", set$unit, "). Adjust d.min?")
  }

  hist3 <- function(d, BCAD) {
    hsts <- list(); maxhist <- 0; minhist <- 1
    if(progress && verbose)
      pb <- txtProgressBar(min=0, max=max(1,length(d)-1), style = 3)
    for(i in 1:length(d)) {
      if(progress && verbose)
        if(length(d) > 1) {
          if(length(d) < 500) # progress bar slows things down much when i is large
            setTxtProgressBar(pb, i) else
              if(i %% 10 == 0)
                setTxtProgressBar(pb, i)
        }

      ages <- Bacon.Age.d(d[i], set, BCAD=BCAD)

      if(length(ages[!is.na(ages)]) > 0) { # added !is.na(ages) 21 April 21
        hst <- density(ages)
        th0 <- min(hst$x)
        th1 <- max(hst$x)
        maxhist <- max(maxhist, hst$y)
        minhist <- min(minhist, max(hst$y))
        n <- length(hst$x)
        counts <- hst$y
        ds <- d[i]
        hsts <- append(hsts, pairlist(list(d=ds, th0=th0, th1=th1, 
          n=n, counts=counts, max=maxhist, min=minhist)))
      } else hsts$d[[i]] <- d[i]
    }
    return(hsts)
  }

  hists <- hist3(d, BCAD)
  if(save.info)
    assign_to_global("hists", hists)

  rng <- array(NA, dim=c(length(d), 4)) # R > 4.0 does not like to fill c() using loops
  if(calc.range) {
    if(length(d) == 1)
      rng <- rbind(ageranges(d, set=set, BCAD=BCAD, prob=prob)) else
        rng <- ageranges(d, set=set, BCAD=BCAD, prob=prob)[,-1]
    for(i in 1:length(d))
      hists$rng[[i]] <- rng[i,]
  }

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
      segments(rng[,2], 0, rng[,3], 0, col=range.col, lwd=3)
      points(rng[,4], 0, col=med.col, pch=20)
      points(rng[,5], 0, col=mean.col, pch=20)

      if(verbose) {
        message("mean (", mean.col, "): ", round(rng[4],1), " ", age.lab,
          ", median (", med.col, "): ",  round(rng[3],1), " ", age.lab, "\n")
        message(100*prob, "% range (", range.col, "): ", round(rng[1],1), " to ", round(rng[2],1), " ", age.lab, "\n")
      }
    }
  #invisible(rng)
  invisible(hists)
}



# to calculate age ranges. We're now using (documented) ageranges instead
Bacon.rng <- function(d, set=get('info'), BCAD=set$BCAD, prob=set$prob, verbose=TRUE) {
  outfile <- paste0(set$prefix, ".out")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- Bacon.AnaOut(outfile, set, MCMC.resample=FALSE)
    assign_to_global("set", set) # MB Jan 2024; can this go?
  }

  if(length(d) > 1)
    pb <- txtProgressBar(min=0, max=max(1, length(d)-1), style=3)
  rng <- array(NA, dim=c(length(d), 4))

  for(i in 1:length(d)) {
    ages <- Bacon.Age.d(d[i], set, BCAD=BCAD)
    ages <- ages[!is.na(ages)]
    if(length(ages) > 0) {
      rng[i,1:3] <- quantile(ages, c(((1-prob)/2), 1-((1-prob)/2), .5), na.rm=TRUE)
      rng[i,4] <- mean(ages)
    }
    if(length(d) > 1) {
      if(verbose)
        if(length(d) < 500) # progress bar slows things down when i is large
          setTxtProgressBar(pb, i) else
            if(i %% 10 == 0)
              setTxtProgressBar(pb, i)
      }
  }
  return(rng)
}



#' @name agemodel.it
#' @title extract one age-model iteration
#' @description For one MCMC iteration (it), extract the corresponding age-depth model.
#' @param it The MCMC iteration of which the age-model should be calculated.
#' @param depths The depths for which the age estimates are to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param save.info If TRUE, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A variable with two columns - depth and the age-depth model of a single iteration.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=50, d.res=50, d.by=10)
#'   lines(agemodel.it(5), col="red")
#' }
#' @export
agemodel.it <- function(it, depths=c(), set=get('info'), BCAD=set$BCAD, save.info=set$save.info) {
  if(it <1 || length(it) > 1)
    stop("can only handle one iteration at a time", call.=FALSE)
  outfile <- paste0(set$prefix, ".out")
  if(length(set$output) == 0 || length(set$Tr) == 0) {
    set <- Bacon.AnaOut(outfile, set, MCMC.resample=FALSE)
    if(save.info)
      assign_to_global("info", set) # changed 'set' to 'info'
  }
  out <- set$output
  if(!it %in% 1:nrow(out))
    stop(paste0("You asked for iteration number ", it, ", but the run has only ", nrow(out), " iterations"), call.=FALSE)

  if(length(depths) == 0)
    depths <- set$depths

  ages <- sapply(depths, function(d) {
    a <- Bacon.Age.d(d, set, BCAD=BCAD, its=set$output[it,,drop=FALSE])
    if(length(a) > 1) NA else a # slumps can cause havoc
  })

  return(cbind(depths=depths, ages=ages))
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



#' @name squeeze
#' @title squeeze some depths of a core
#' @description Squeeze or compress depths below a boundary by a certain amount. Accompanies the stretch function; see the stretch function for code on running the accordion
#' @param d The depth(s) to be squeezed
#' @param boundary The depth below which depths should be squeezed
#' @param times The factor by which the depths should be squeezed
#' @author Maarten Blaauw
#' @return The squeezed depth(s)
#' @examples
#'   squeeze(40,25,20)
#' @export
squeeze <- function(d, boundary, times) {
  below <- which(d > boundary)
  if(length(below) > 0)
    d[below] <- boundary + (d[below]-boundary)/times
  return(d)
}



#' @name stretch
#' @title stretch some depths of a core
#' @description Stretch squeezed depths e.g., calculate the original depths of depths that were squeezed. Accompanies the squeeze function.
#' @param d The depth(s) to be stretched
#' @param boundary The depth below which depths should be stretched
#' @param times The factor by which the depths should be stretched
#' @author Maarten Blaauw
#' @return The stretched depth(s)
#' @examples
#'   stretch(25.75,25,20)
#' \dontrun{
#'   # To play the accordion, first squeeze an existing core.
#'   # Let's squeeze the depths below 10 cm core depth 20 times:
#'   Bacon("accordion", 1)
#'   dets <- info$dets
#'   dets[,4] <- squeeze(dets[,4], 10, 20)
#'
#'   # make a new directory for the squeezed core, and place the dets file there:
#'   dir.create("Bacon_runs/squeezed")
#'   write.table(dets, "Bacon_runs/squeezed/squeezed.csv", row.names=FALSE, sep=",")
#'
#'   # now run that squeezed core, adding a boundary (10cm) and adapting the acc.mean prior (20x):
#'   Bacon("squeezed", 1, boundary=10, acc.mean=c(5, 20*5))
#'   # finally, plot while stretching the depths onto the original scale:
#'   agedepth(accordion=c(10,20))
#' }
#' @export
stretch <- function(d, boundary, times) {
  below <- which(d > boundary)
  if(length(below) > 0)
    d[below] <- boundary + (d[below]-boundary)*times
  return(d)
}
