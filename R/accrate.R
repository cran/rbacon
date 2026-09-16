

## Accumulation rate calculations
#' should take into account hiatuses and slumps (I think it does though)
#' @name accrate.depth
#' @title obtain estimated accumulation rates as for any depth of a core
#' @description Obtain accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any depth of a core.
#' @details Considering accumulation rates is crucial for age-depth modelling, and even more so if they are subsequently used for calculating proxy
#' influx values, or interpreted as proxy for environmental change such as carbon accumulation.
#' Bacon deals explicitly with accumulation rate and its variability through defining prior distributions.
#' This function obtains accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations
#' for any depth of a core. Deals with only 1 depth at a time. See also \code{accrate.age}.
#' @param d The depth for which accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
#' @param remove.hiatuses Any hiatuses will affect apparent accumulation rates within sections. Therefore, by default the hiatus jumps will be removed from the accumulation rates within sections containing hiatuses.
#' @param na.rm Remove NA entries. These are NOT removed by default, ensuring that always the same amount of iterations is returned.
#' @param inversion.threshold Very small accumulation rate values will become very large when their inverse is calculated. By default, any accumulation rate smaller than 1e-6 is set to 1e-6.
#' @author Maarten Blaauw, J. Andres Christen
#' @return all MCMC estimates of accumulation rate of the chosen depth.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   d20 <- accrate.depth(20)
#'   hist(d20)
#'   d20 <- accrate.depth(20, cmyr=TRUE) # to calculate accumulation rates in cm/yr
#'   mean(d20)
#' }
#' @export
accrate.depth <- function(d, set=get('info'), cmyr=FALSE, remove.hiatuses=TRUE, na.rm=FALSE, inversion.threshold=1e-6) {
  if(is.na(d))
    return(NA)

  d.hiatus <- set$hiatus.depths
  if(length(set$slump) > 0) {
    slump <- set$slump
    d.hiatus <- set$slumphiatus
    for(i in 1:nrow(slump))
      if(d >= min(slump[i,]) && d < max(slump[i,]))
        return(NA) # don't bother with depths within slumps
    d <- toslump(d, slump) # work with slumpfree depths
  }

  accs.elbows <- set$output[,2:(set$K+1)]
  if(all(!is.na(set$elbows)) && min(set$elbows) <= d && d <= max(set$elbows))
    accs <- unlist(accs.elbows[max(which(set$elbows <= d))]) else
      accs <- NA

  if(!is.na(d.hiatus[1]))
    if(remove.hiatuses)
      for(i in 1:length(d.hiatus)) {
        k <- max(which(set$elbows <= d.hiatus[i])) # elbow just above the hiatus
        if(d >= set$elbows[k] && d < set$elbows[k+1]) # d is within a hiatus
          if(d <= d.hiatus[i])
            accs <- set$slope.above[,i] else
              accs <- set$slope.below[,i]
      }

  accs <- as.numeric(accs)
  if(na.rm)
    accs <- accs[!is.na(accs)]
  if(cmyr) {
    accs[accs < inversion.threshold] <- inversion.threshold
    accs <- 1/accs
  }
  return(accs)
}



# should take into account hiatuses
#' @name accrate.age
#' @title obtain estimated accumulation rates for any age of a core
#' @description Obtain accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any age of a core.
#' @details Considering accumulation rates is crucial for age-depth modelling, and even more so if they are subsequently
#' used for calculating proxy influx values, or interpreted as proxy for environmental change such as carbon accumulation. See also \code{accrate.age.ghost}, \code{accrate.depth} and \code{accrate.depth.ghost}.
#' Bacon deals explicitly with accumulation rate and its variability through defining prior distributions.
#' This function obtains accumulation rates (in years per cm, so actually sedimentation times) as estimated
#' by the MCMC iterations for any age of a core. Deals with only 1 age at a time. See also \code{accrate.depth}.
#' @param age The age for which the accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
#' @param ages The ages of the age-depth model. Not provided by default, but can be provided to speed things up if the function is called repeatedly
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param silent Warn when ages are outside the core's range. Default \code{silent=TRUE}.
#' @param na.rm Remove NA entries. These are NOT removed by default, ensuring that always the same amount of iterations is returned.
#' @author Maarten Blaauw, J. Andres Christen
#' @return all MCMC estimates of accumulation rate of the chosen age.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   accrate.a5000 <- accrate.age(5000)
#'   plot(accrate.a5000, pch='.')
#'   hist(accrate.a5000)
#' }
#' @export
accrate.age <- function(age, set=get('info'), cmyr=FALSE, ages=c(), BCAD=set$BCAD, silent=TRUE, na.rm=FALSE) {
  if(length(ages) == 0)
    ages <- sapply(set$elbows, Bacon.Age.d)
  if(BCAD)
    ages <- BCADtocalBP(ages)

  if(!silent)
    if(age < min(ages) || age > max(ages))
      stop(" Warning, age outside the core's age range!\n")

   hiatus.sections <- c()
   if(!is.na(set$hiatus.depths[1])) # was ... length==0
     for(i in 1:length(set$hiatus.depths))
       hiatus.sections[i] <- max(which(set$elbows <= set$hiatus.depths[i]))

  accs <- rep(NA_real_, nrow(ages)) # suggested by henningte on github
  for(i in 2:ncol(ages)) {
    these <- (ages[,i-1] < age) & (ages[,i] > age)

    # if(length(set$hiatus.depth) >0) then check if the age falls within a hiatus section
    # if it does, then check if it is below or above it, and use set$slope.below or set$slope.above accordingly
    if(any(these)) { # age lies within these age-model iterations
      rows <- which(these)
      accs[rows] <- set$output[rows,i] # assign the accumulation rates

      if(length(hiatus.sections) > 0)
        if((i-1) %in% hiatus.sections) {
          j <- which(hiatus.sections == (i-1))
          below <- age > set$hiatus.start[rows, j] # which of these lie below the hiatus
          above <- age < set$hiatus.end[rows, j] # which lie above it
          inside <- age <= set$hiatus.start[rows, j] & age >= set$hiatus.end[rows, j]

          accs[rows[below]] <- set$slope.below[rows[below],j]
          accs[rows[above]] <- set$slope.above[rows[above],j]
          accs[rows[inside]] <- NA
        }
    }
  }

  if(na.rm)
    accs <- accs[!is.na(accs)]
  if(cmyr)
    accs <- 1/accs

  return(accs)
}



#' @name accrate.depth.summary
#' @title provide a summary of the estimated accumulation rates for any depth of a core
#' @description Obtain a summary (95\% range, 68\% range, 50\%=median, mean) of the accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any depth of a core.
#' @param d The depth for which accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
#' @param remove.hiatuses Hiatuses will affect apparent accumulation rates within sections. Therefore, by default the hiatus jumps will be removed from the accumulation rates within sections containing hiatuses.
#' @param na.rm Remove NA entries. These are NOT removed by default, so that always the same amount of iterations is returned. NAs will however be removed if a core has slumps.
#' @param probs The probability ranges to be returned. Defaults to the minima and maxima of the 95\% and 68\% ranges, as well as the median: \code{probs=c(.025, .16, .84, .975, .5)}.
#' @author Maarten Blaauw
#' @return A summary of the estimated accumulation rate of the chosen depth: minimum of the 95\% interval, minimum of the 68\% interval, maximum of the 68\% interval, maximum of the 95\% interval, median (i.e., 50\%) and mean.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   accrate.depth.summary(20)
#' }
#' @export
accrate.depth.summary <- function(d, set=get('info'), cmyr=FALSE, remove.hiatuses=TRUE, na.rm=FALSE, probs=c(.025, .16, .84, .975, .5)) {
  if(length(d) > 1)
    stop("can handle one depth at a time only")
  accs <- accrate.depth(d, set, cmyr, remove.hiatuses=remove.hiatuses, na.rm=na.rm)
  qu <- quantile(accs, probs, na.rm=na.rm)
  mn <- mean(accs, na.rm=na.rm)
  names(mn) <- "mean"
  return(c(qu, mn))
}



# should take into account slumps and hiatuses
#' @name accrate.age.summary
#' @title provide a summary of the estimated accumulation rates for any age of a core
#' @description Obtain a summary (95\% range, 68\% range, median, mean) of the accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any age of a core.
#' @param age The age for which accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
#' @param na.rm Remove NA entries. These are NOT removed by default, so that always the same amount of iterations is returned.
#' @param probs The probability ranges to be returned. Defaults to the minima and maxima of the 95\% and 68\% ranges, as well as the median: \code{probs=c(.025, .16, .84, .975, .5)}.
#' @author Maarten Blaauw
#' @return A summary of the estimated accumulation rate of the chosen depth: minimum of the 95\% interval, minimum of the 68\% interval, maximum of the 68\% interval, maximum of the 95\% interval, median (i.e., 50\%) and mean.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   accrate.age.summary(5000)
#' }
#' @export
accrate.age.summary <- function(age, set=get('info'), cmyr=FALSE, na.rm=TRUE, probs=c(.025, .16, .84, .975, .5)) {
  if(length(age) > 1)
    stop("can handle one depth at a time only")
  accs <- accrate.age(age, set, cmyr, na.rm=na.rm)
  qu <- quantile(accs, probs, na.rm=na.rm)
  mn <- mean(accs, na.rm=na.rm)
  names(mn) <- "mean"
  return(c(mn, qu))
}



# should take into account slumps and hiatuses
#' @name accrates.core
#' @title provide a summary of the estimated accumulation rates for a range of core depths
#' @description Obtain a summary (95\% range, 68\% range, median, mean) of the accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for a range of depths of a core, and optionally write this as a file to the core directory (ending in '_accrates.txt').
#' @param dseq The sequence of depths for which accumulation rates need to be returned. Defaults to whatever info$dseq is, which most often is a sequence from the top to the bottom of the core at 1 cm increments.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
#' @param remove.hiatuses Any hiatuses will affect apparent accumulation rates within sections. Therefore, by default the hiatus jumps will be removed from the accumulation rates within sections containing hiatuses.
#' @param na.rm Remove NA entries. These are NOT removed by default, so that always the same amount of iterations is returned.
#' @param probs The probability ranges to be returned. Defaults to the minima and maxima of the 95\% and 68\% ranges, as well as the median: \code{probs=c(.025, .16, .84, .975, .5)}.
#' @param round The number of decimals to report. Defaults to \code{round=2}.
#' @param write Whether or not to write the summary to a file, in the core's directory and ending in `_accrates.txt`.
#' @param sep Character to separate the entries within the file. Defaults to a tab, \code{sep="\t"}.
#' @author Maarten Blaauw
#' @return A summary of the estimated accumulation rate for all selected depths: minimum of the 95\% interval, minimum of the 68\% interval, maximum of the 68\% interval, maximum of the 95\% interval, median (i.e., 50\%) and mean. This is optionally written to a file in the core directory.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   myaccrates <- accrates.core()
#' }
#' @export
accrates.core <- function(dseq=c(), set=get('info'), cmyr=FALSE, remove.hiatuses=TRUE, na.rm=TRUE, probs=c(.025, .16, .84, .975, .5), round=2, write=TRUE, sep="\t") {
  if(length(dseq) == 0)
    dseq <- set$depths
  
  if(length(set$slump) > 0)
    na.rm <- TRUE  
  mysummary <- function(dseq)
    accrate.depth.summary(dseq, set, cmyr, remove.hiatuses=remove.hiatuses, na.rm, probs)
  allaccs <- t(sapply(dseq, mysummary))
  allaccs <- cbind(depths=dseq, round(allaccs, round))
  
  if(write) {
    fl <- paste0(set$coredir, set$core, "/", set$core, "_", set$K, "_accrates.txt")
    message("writing the accumulation rate summary to ", fl)
    write.table(allaccs, fl, sep=sep, quote=FALSE, row.names=FALSE)
  } 
  invisible(allaccs)
}




#' @name accrate.depth.ghost
#' @title plot modelled accumulation rates against the depths of a core
#' @description Plot grey-scale representation of modelled accumulation rates over a core's depth. Each section of the core (see Bacon's option \code{"thick"}) will have modelled accumulation rates.
#' @details This plot shows the modelled accumulation rates in grey-scales, where darker grey indicates more likely accumulation rates.
#' Axis limits for accumulation rates are estimated automatically, however upper limits can be very variable (and thus hard to predict)
#' if calculated in cm/yr; therefore you might want to manually adapt the axis limits after plotting with default settings (e.g., \code{acc.lim=c(0,1)}). See also \code{accrate.age.ghost}, \code{accrate.depth} and \code{accrate.age}.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param d The depths for which the accumulation rates are to be calculated. Default to the entire core.
#' @param d.lim Axis limits for the depths.
#' @param acc.lim Axis limits for the accumulation rates.
#' @param d.lab Label for the depth axis.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm. Axis limits are difficult to calculate when \code{cmyr=TRUE}, so a manual adaptation of \code{acc.lim} might be a good idea.
#' @param acc.lab Axis label for the accumulation rate.
#' @param dark The value beyond which any higher values will be clipped. Set to 1 by default. This can be set to lower values if for example only a very small area in the plot obtains the darkest values - then lowering the `dark` value will darker a larger area.
#' @param darkest The darkest grey value is darkest=1 by default; lower values will result in lighter maximum colours; values >1 are not advised.
#' @param cutoff Point below which colours will no longer be printed. Default \code{cutoff=0.001}.
#' @param zero.col The colour where the ghost is 0. Defaults to \code{zero.col="white"}, which together with \code{max.col="black"} results in a greyscale. More creative colour gradients can be implemented by checking the 600+ colours in \code{colours()}.
#' @param max.col The colour where the ghost is at its maximum. Defaults to \code{max.col="black"}, which together with \code{zero.col="white"} results in a greyscale. More creative colour gradients can be implemented by checking the 600+ colours in \code{colours()}.
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param prob Probability ranges. Defaults to \code{prob=0.95}.
#' @param plot.range If \code{plot.range=TRUE}, the confidence ranges (two-tailed; half of the probability at each side) are plotted.
#' @param range.col Colour of the confidence ranges.
#' @param range.lty Line type of the confidence ranges.
#' @param plot.mean If \code{plot.mean=TRUE}, the means are plotted.
#' @param mean.col Colour of the mean accumulation rates.
#' @param mean.lty Type of the mean lines.
#' @param plot.median If \code{plot.mean=TRUE}, the medians are plotted.
#' @param median.col Colour of the median accumulation rates.
#' @param median.lty Type of the median lines.
#' @param rotate.axes The default is to plot the accumulation rates horizontally and the depth vertically (\code{rotate.axes=FALSE}). Change rotate.axes value to rotate axes.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}.
#' @param rev.acc The direction of the accumulation rate axis can be reversed from the default (\code{rev.acc=TRUE}).
#' @param xaxs Extension of x-axis. By default, add some extra white-space at both extremes (\code{xaxs="r"}). See ?par for other options.
#' @param yaxs Extension of y-axis. By default, add no extra white-space at both extremes (\code{yaxs="i"}). See ?par for other options.
#' @param bty Type of box to be drawn around the plot (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param remove.laststep Add a white line to remove spurious lines at the extreme of the graph. Defaults to TRUE.
#' @param use.raster Rasters can be aligned or not in the underlying image function. Setting \code{use.raster=FALSE, default} takes a bit longer to draw and sometimes causes strange lines owing to anti-aliasing. However, the alternative of \code{use.raster=TRUE} causes greyscales on some devices (e.g., OSX quartz) to 'flip'. If this is the case, use 'flip.acc=TRUE'.
#' @param flip.acc When using \code{use.raster=TRUE}, sometimes greyscales are flipped. If this is the case, see if setting \code{flip.acc=TRUE} solves this. 
#' @param remove.hiatuses Any hiatuses will affect apparent accumulation rates within sections. Therefore, by default the hiatus jumps will be removed from the accumulation rates within sections containing hiatuses.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A grey-scale plot of accumulation rate against core depth, and (invisibly) the list of depths and their accumulation rates (ranges, medians, means).
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   layout(1)
#'   tmp <- accrate.depth.ghost()
#'   head(tmp)
#' }
#' @export
accrate.depth.ghost <- function(set=get('info'), d=set$elbows, d.lim=c(), acc.lim=c(), d.lab=c(), cmyr=FALSE, acc.lab=c(), dark=1, darkest=.8, cutoff=0.001, zero.col="white", max.col="black", rgb.scale=c(0,0,0), rgb.res=100, prob=0.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, plot.median=TRUE, median.col="blue", median.lty=2, rotate.axes=FALSE, rev.d=FALSE, rev.acc=FALSE, xaxs="r", yaxs="r", bty="l", remove.laststep=TRUE, use.raster=FALSE, flip.acc=FALSE, remove.hiatuses=TRUE) {

  max.acc <- 0; max.dens <- 0; max.acc2 <- 0
  acc <- list(); min.rng <- numeric(length(d)); max.rng <- numeric(length(d))
  mean.rng <- numeric(length(d)); median.rng <- numeric(length(d))

  slump <- set$slump
  inslump <- integer(0)
  if(length(slump) > 0) # remove depths within slump from the analysis
    for(i in 1:nrow(slump))
      inslump <- c(inslump, which(d >= min(slump[i,]) & d <= max(slump[i,])))
  keep <- setdiff(seq_along(d), inslump) # only work with non-slump depths
  d <- d[keep]

  for(i in 1:length(d)) {
    d.acc <- accrate.depth(d[i], set, cmyr=cmyr, remove.hiatuses=remove.hiatuses)
    if(length(acc.lim) == 0)
      acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr, remove.hiatuses=remove.hiatuses), from=0, na.rm=TRUE) else
        acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr, remove.hiatuses=remove.hiatuses), from=0, to=max(acc.lim, na.rm=TRUE), na.rm=TRUE)
  }

  for(i in 1:length(d)) {
    max.acc <- max(max.acc, acc[[i]]$x, na.rm=TRUE)
    max.acc2 <- max(max.acc2, quantile(acc[[i]]$x, .99, na.rm=TRUE)) # take a value close to the max
    max.dens <- max(max.dens, acc[[i]]$y, na.rm=TRUE)
    accs <- accrate.depth(d[i], set, cmyr=cmyr, remove.hiatuses=remove.hiatuses, na.rm=TRUE)
    quants <- quantile(accs, c((1-prob)/2, 1-((1-prob)/2)), na.rm=TRUE)
    min.rng[i] <- quants[1]
    max.rng[i] <- quants[2]
    mean.rng[i] <- mean(accs)
    median.rng[i] <- median(accs)
  }

#  if(length(inslump) > 0)
#    stored <- cbind(d, min.rng[-inslump], max.rng[-inslump], median.rng[-inslump], mean.rng[-inslump]) else
      stored <- cbind(d, min.rng[keep], max.rng[keep], median.rng[keep], mean.rng[keep])

#  stored <- cbind(d, min.rng[-inslump], max.rng[-inslump], median.rng[-inslump], mean.rng[-inslump])
  colnames(stored) <- c("depth", "min.rng", "max.rng", "median", "mean")

  for(i in 1:length(d)) {  
    acc[[i]]$y <- acc[[i]]$y/(dark*max.dens)
    acc[[i]]$y[acc[[i]]$y > 1] <- 1 # set "dark" to black
    acc[[i]]$y[acc[[i]]$y < cutoff] <- NA # do not plot too light/small values
  }

  if(length(d.lim) == 0)
    d.lim <- range(set$dets[,4], na.rm=TRUE) # avoiding any slump-based definitions of d
  if(length(d.lab) == 0)
    d.lab <- paste0("depth (", set$depth.unit, ")")
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste0("accumulation rate (", set$depth.unit, "/", set$age.unit, ")") else
        acc.lab <- paste0("accumulation rate (", set$age.unit, "/", set$depth.unit, ")")

  if(rev.d)
    d.lim <- rev(d.lim)
  if(length(acc.lim) == 0)
    acc.lim <- c(0, max.acc2)
  if(rev.acc)
    acc.lim <- rev(acc.lim)

  if(is.na(max.col))
    col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(max(accs$y[!is.na(accs$y)]), 0, length=rgb.res)) else
      col <- col.scales(rgb.res, zero.colour=zero.col, max.colour=max.col, dark=dark, darkest=darkest)

  if(rotate.axes) {
    plot(0, type="n", xlab=acc.lab, ylab=d.lab, ylim=d.lim, xlim=acc.lim, bty="n", xaxs=xaxs, yaxs=yaxs)
    for(i in 2:length(d)) {
      accs <- acc[[i-1]]
      #if(is.null(accs))
      #  next
      z <- if (flip.acc) t(rev(accs$y)) else t(accs$y)
      if(deviceIsQuartz()) 
        if(use.raster)
          if(rev.acc)
            z <- t(rev(z))

      image(accs$x, d[c(i - 1, i)], t(z), add=TRUE, col=col, useRaster=use.raster)
    }
    if(plot.range)
      for(i in 2:(length(d))) {
        segments(min.rng[i-1], d[i-1], min.rng[i-1], d[i], col=range.col, lty=range.lty)
        segments(min.rng[i-1], d[i], min.rng[i], d[i], col=range.col, lty=range.lty)
        segments(max.rng[i-1], d[i-1], max.rng[i-1], d[i], col=range.col, lty=range.lty)
        segments(max.rng[i-1], d[i], max.rng[i], d[i], col=range.col, lty=range.lty)
      }
    if(plot.mean)
      for(i in 2:length((d))) {
        segments(mean.rng[i-1], d[i-1], mean.rng[i-1], d[i], col=mean.col, lty=mean.lty)
        segments(mean.rng[i-1], d[i], mean.rng[i], d[i], col=mean.col, lty=mean.lty)
      }
    if(plot.median)
      for(i in 2:length((d))) {
        segments(median.rng[i-1], d[i-1], median.rng[i-1], d[i], col=median.col, lty=median.lty)
        segments(median.rng[i-1], d[i], median.rng[i], d[i], col=median.col, lty=median.lty)
      }
    if(remove.laststep)
      abline(h=min(set$elbows), col="white", lwd=2)
  } else {
      plot(0, type="n", xlab=d.lab, ylab=acc.lab, xlim=d.lim, ylim=acc.lim, bty="n", xaxs=xaxs, yaxs=yaxs)
      for(i in 2:length(d)) {  
        accs <- acc[[i-1]]
        if(is.null(accs))
          next
        z <- if (flip.acc) t(rev(accs$y)) else t(accs$y)
        if(deviceIsQuartz()) 
          if(use.raster)
            z <- t(z[length(z):1])
        image(d[c(i - 1, i)], accs$x, z, add=TRUE, col=col, useRaster=use.raster)
      }

    if(plot.range) {
      lines(d, min.rng[keep], type="s", col=range.col, lty=range.lty, pch=NA)
      lines(d, max.rng[keep], type="s", col=range.col, lty=range.lty, pch=NA)
    }
    if(plot.mean)
      lines(d, mean.rng[keep], type="s", col=mean.col, lty=mean.lty)
    if(plot.median)
      lines(d, median.rng[keep], type="s", col=median.col, lty=median.lty)
    if(remove.laststep)
      abline(v=max(set$elbows), col="white", lwd=1.5)
    }

  box(bty=bty)  
  invisible(stored)
}



#' @name accrate.age.ghost
#' @title plot a core's accumulation rates against calendar time
#' @description Plot a grey-scale representation of a core's estimated accumulation rates against time.
#' @details Calculating accumulation rates against calendar age will take some time to calculate, and might show unexpected
#' rates around the core's maximum ages (only a few of all age-model iterations will reach such ages and they will tend to have
#'  modelled accumulation rates for the lower depths much lower than the other iterations). Axis limits for accumulation rates
#'   are estimated automatically, however upper limits can be very variable (and thus hard to predict) if calculated in \code{cm/yr}.
#'  Therefore you might want to manually adapt the axis limits after plotting with default settings (e.g., \code{acc.lim=c(0,1)}). See also \code{accrate.depth.ghost}, \code{accrate.depth} and \code{accrate.age}.
#' The grey-scale reconstruction around the oldest ages of any reconstruction often indicates very low accumulation rates.
#' This is due to only some MCMC iterations reaching those old ages, and these iterations will have modelled very slow accumulation rates.
#' Currently does not deal well with hiatuses, so do not interpret accumulation rates close to depths with inferred hiatuses.
#' If warning messages such as "In min(x) : no non-missing arguments to min; returning Inf" appear, Bacon likely became confused about whether or not to use BC/AD. Best run your core again with the desired setting for BCAD. 
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param na.rm Remove NA entries. These are NOT removed by default, ensuring that always the same amount of iterations is returned.
#' @param kcal Use kcal BP. Default is \code{kcal=FALSE}.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the plot. Default \code{age.res=400}.
#' @param acc.res Resolution or amount of greyscale pixels to cover the accumulation rate scale plot. Default \code{age.res=400}.
#' @param cutoff Point below which colours will no longer be printed. Default \code{cutoff=0.001}.
#' @param zero.col The colour where the ghost is 0. Defaults to \code{zero.col="white"}, which together with \code{max.col="black"} results in a greyscale. More creative colour gradients can be implemented by checking the 600+ colours in \code{colours()}.
#' @param max.col The colour where the ghost is at its maximum. Defaults to \code{max.col="black"}, which together with \code{zero.col="white"} results in a greyscale. More creative colour gradients can be implemented by checking the 600+ colours in \code{colours()}.
#' @param dark The value beyond which any higher values will be clipped. Set to 1 by default. This can be set to lower values if for example only a very small area in the plot obtains the darkest values - then lowering the `dark` value will darker a larger area.
#' @param darkest The darkest grey value is darkest=1 by default; lower values will result in lighter maximum colours; values >1 are not advised.
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param prob Probability ranges. Defaults to \code{prob=0.95}.
#' @param plot.range If \code{plot.range=TRUE}, the confidence ranges (two-tailed; half of the probability at each side) are plotted.
#' @param range.col Colour of the confidence ranges.
#' @param range.lty Line type of the confidence ranges.
#' @param plot.mean If \code{plot.mean=TRUE}, the means are plotted.
#' @param mean.col Colour of the mean accumulation rates.
#' @param mean.lty Type of the mean lines.
#' @param plot.median If \code{plot.mean=TRUE}, the medians are plotted.
#' @param median.col Colour of the median accumulation rates.
#' @param median.lty Type of the median lines.
#' @param acc.lim Axis limits for the accumulation rates.
#' @param acc.lab Axis label for the accumulation rate.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm. Axis limits are difficult to calculate when \code{cmyr=TRUE}, so a manual adaptation of \code{acc.lim} might be a good idea.
#' @param rotate.axes The default is to plot the calendar age horizontally and accumulation rates vertically. Change to \code{rotate.axes=TRUE} value to rotate axes.
#' @param rev.age The direction of the age axis, which can be reversed using \code{rev.age=TRUE}.
#' @param rev.acc The direction of the accumulation rate axis, which can be reversed (\code{rev.acc=TRUE}.
#' @param xaxs Extension of the x-axis. White space can be added to the vertical axis using \code{xaxs="r"}.
#' @param yaxs Extension of the y-axis. White space can be added to the vertical axis using \code{yaxs="r"}.
#' @param bty Type of box to be drawn around the plot (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param use.raster Rasters can be aligned or not in the underlying image function. Setting \code{use.raster=FALSE} takes a bit longer to draw and sometimes causes strange lines owing to anti-aliasing. Therefore, \code{use.raster=TRUE} is the default, however on some devices (e.g., OSX quartz) this causes greyscales to 'flip'. If this is the case, use 'flip.acc=TRUE'.
#' @param flip.acc When using \code{use.raster=TRUE}, sometimes greyscales are flipped. If this is the case, see if setting \code{flip.acc=TRUE} solves this. 
#' @param flip.age When using \code{use.raster=TRUE}, sometimes greyscales are flipped. If this is the case, see if setting \code{flip.age=TRUE} solves this. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return A greyscale plot of accumulation rate against calendar age, and (invisibly) the list of ages and their accumulation rates (ranges, medians, means).
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=20, d.res=20, d.by=10)
#'   layout(1)
#'   tmp <- accrate.age.ghost(age.res=200, acc.res=100)
#'   head(tmp)
#' }
#' @export
accrate.age.ghost <- function(set=get('info'), age.lim=c(), age.lab=c(), na.rm=TRUE, kcal=FALSE, age.res=400, acc.res=200, cutoff=.001, zero.col="white", max.col="black", dark=1, darkest=1, rgb.scale=c(0,0,0), rgb.res=100, prob=.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, plot.median=TRUE, median.col="blue", median.lty=2, acc.lim=c(), acc.lab=c(), BCAD=set$BCAD, cmyr=FALSE, rotate.axes=FALSE, rev.age=FALSE, rev.acc=FALSE, use.raster=FALSE, flip.acc=FALSE, flip.age=FALSE, xaxs="i", yaxs="i", bty="l") {
  if(length(age.lim) == 0) 
     age.lim <- extendrange(set$ranges[,5]) # just the mean ages, not the extremes
  if(set$BCAD) # was set$BCAD
    age.lim <- rice::BCADtocalBP(age.lim) # work with cal BP internally
  age.seq <- seq(min(age.lim), max(age.lim), length=age.res)
    
  if(length(acc.lim) == 0) {
    acc.lim <- c(0, 1.05*quantile(unlist(set$output[,2:(1+set$K)]), .999, na.rm=TRUE)) # almost-maximum accrate in the output
    #acc.lim <- c(0, 1.05*max(set$output[,2:(1+set$K)], na.rm=TRUE))
    if(cmyr)
      acc.lim <- 1/acc.lim
    acc.lim[is.infinite(acc.lim)] <- 0  
  }  
  acc.seq <- seq(min(acc.lim, na.rm=TRUE), max(acc.lim, na.rm=TRUE), length=acc.res)
  acc_min <- min(acc.lim, na.rm=TRUE)
  acc_max <- max(acc.lim, na.rm=TRUE)
  
  z <- array(0, dim=c(acc.res, age.res)) # accs in rows, ages in columns
  acc.rng <- array(NA, dim=c(age.res, 2))
  acc.mean <- rep(NA, age.res); acc.median <- acc.mean

  # speed things up by not repeatedly calculating ages in accrate.age
  ages <- array(0, dim=c(nrow(set$output), length(set$elbows)))
  for(i in 1:ncol(ages))
    ages[,i] <- Bacon.Age.d(set$elbows[i], set, BCAD=BCAD) # BCAD was F, June '25

  pb <- txtProgressBar(min=0, max=max(1,length(age.seq)-1), style = 3)
  for(i in 1:age.res) {
    setTxtProgressBar(pb, i)
    acc <- accrate.age(age.seq[i], set, cmyr=cmyr, ages=ages, silent=TRUE, BCAD=BCAD, na.rm=na.rm) # BCAD was F, June '25
    acc <- acc[!is.na(acc)]
    if(length(acc) > 1) {
      z[,i] <- density(acc, from=acc_min, to=acc_max, n=acc.res)$y
      acc.rng[i,] <- quantile(acc, c((1-prob)/2, 1-((1-prob)/2)))
      acc.mean[i] <- mean(acc)
      acc.median[i] <- median(acc)
    }
  }
  message("") # print a newline
  
  stored <- cbind(age.seq, acc.rng[,1], acc.rng[,2], acc.median, acc.mean)
  colnames(stored) <- c("ages", "min.rng", "max.rng", "median", "mean")
  z <- t(z) # when using image to draw the greyscales, it will rotate z
  if(flip.acc)
    z <- z[,ncol(z):1] 
  if(flip.age)
    z <- z[nrow(z):1,]  
  #z <- z/(dark*max(z)) # normalise, set dark to black
  z <- z /(dark*quantile(z, .999, na.rm=TRUE))
  z[z>1] <- 1 # avoid values > 1
  z[z<cutoff] <- NA # do not plot very small/light greyscale values

 # if(deviceIsQuartz()) 
 #   if(use.raster)
 #     z <- z[,ncol(z):1]
  if(rev.acc) {
    acc.lim <- rev(acc.lim)  
    if(use.raster)
      if(deviceIsQuartz()) 
         z <- z[,ncol(z):1]
  }
  if(rev.age) {
    age.lim <- rev(age.lim)
    if(use.raster)
      if(deviceIsQuartz()) 
        z <- z[nrow(z):1,]
  }
#  if(BCAD)
#     z <- z[nrow(z):1,]

  if(is.na(max.col))
    cols <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0,1, length=rgb.res)) else
  cols <- col.scales(rgb.res, zero.colour=zero.col, max.colour=max.col, dark=dark, darkest=darkest)

  if(length(age.lab) == 0)
    if(BCAD)
      age.lab <- "BC/AD" else
        age.lab <- ifelse(kcal, "kcal BP", "cal BP")
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste0("accumulation rate (", set$depth.unit, "/", set$age.unit, ")") else
        acc.lab <- paste0("accumulation rate (", set$age.unit, "/", set$depth.unit, ")")

  if(rotate.axes) {
    yaxt <- ifelse(kcal || BCAD, "n", "s")
    plot(0, type="n", ylim=age.lim, ylab=age.lab, xlim=acc.lim, xlab=acc.lab, yaxs=xaxs, xaxs=yaxs, yaxt=yaxt, bty="n")
    if(BCAD)
      axis(2, pretty(age.lim), labels=calBPtoBCAD(pretty(age.lim))) else
        if(kcal)
          axis(2, pretty(age.lim), labels=pretty(age.lim)/1e3)
    #rasterImage(as.raster(img), min(age.seq), min(acc.seq), max(age.seq), max(acc.seq))
    image(acc.seq, age.seq, t(z), col=cols, add=TRUE, useRaster=use.raster)
    if(plot.range) {
      lines(acc.rng[,1], age.seq, pch=".", col=range.col, lty=range.lty)
      lines(acc.rng[,2], age.seq, pch=".", col=range.col, lty=range.lty)
    }
    if(plot.mean) 
      lines(acc.mean, age.seq, col=mean.col, lty=mean.lty)
    if(plot.median) 
      lines(acc.median, age.seq, col=median.col, lty=median.lty)
  } else {
      xaxt <- ifelse(kcal || BCAD, "n", "s")
      plot(0, type="n", xlim=age.lim, xlab=age.lab, ylim=acc.lim, xaxt=xaxt, ylab=acc.lab, xaxs=xaxs, yaxs=yaxs, bty="n")
      if(BCAD)
        axis(1, pretty(age.lim), labels=calBPtoBCAD(pretty(age.lim))) else
        if(kcal)
          axis(1, pretty(age.lim), labels=pretty(age.lim)/1e3)
        image(age.seq, acc.seq, z, add=TRUE, col=cols, useRaster=use.raster)
      # rasterImage(as.raster(img), min(age.seq), min(acc.seq), max(age.seq), max(acc.seq))

      if(plot.range) {
        lines(age.seq, acc.rng[,1], pch=".", col=range.col, lty=range.lty)
        lines(age.seq, acc.rng[,2], pch=".", col=range.col, lty=range.lty)
      }
      if(plot.mean)
        lines(age.seq, acc.mean, col=mean.col, lty=mean.lty)
      if(plot.median)
        lines(age.seq, acc.median, col=median.col, lty=median.lty)
    }

  box(bty=bty)
  invisible(stored)
}



#' @name flux.age.ghost
#' @title plot flux rates for proxies
#' @description Plot grey-scale representation of estimated flux rates for proxies against calendar age.
#' @details To plot flux rates (e.g. pollen grains/cm2/yr) as greyscales,
#' provide a plain text file with headers and the data in columns separated by commas, ending in '_flux.csv'
#' and saved in your core's folder. The first column should contain the depths, and the next columns should contain
#' the proxy concentration values (leaving missing values empty). Then type for example \code{flux.age.ghost(1)} to plot the
#' flux values for the first proxy in the .csv file. Instead of using a _flux.csv file, a flux variable can also be defined
#'  within the R session (consisting of depths and their proxy concentrations in two columns). Then provide the name of this variable, e.g.: \code{flux.age.ghost(flux=flux1)}.
#' See Bacon_runs/MSB2K/MSB2K_flux.csv for an example.
#' @param column Which proxy to use (counting from the column number in the .csv file after the depths column).
#' @param flux Instead of using a file, the data can also be provided as a variable. The first column should be the depths, and the variable 'column' should indicate which column (after the depth column) contains the proxy of interest. For example, if using Plum we could produce a greyscale of the mass accumulation rate: myflux <- info$detsPlum[,c(4,6)];
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param coredir Folder where the core's files \code{core} are and/or will be located. This will be a folder with the core's name, within either the folder \code{coredir='Bacon_runs/'}, or the folder Cores/ if it already exists within R's working directory, or a custom-built folder. For example, use \code{coredir="."} to place the core's folder within the current working directory, or \code{coredir="F:"} if you want to put the core's folder and files on a USB drive loaded under F:.
#' @param remove.hiatuses Hiatuses will affect apparent accumulation rates within sections. Therefore, by default the hiatus jumps will be removed from the accumulation rates within sections containing hiatuses.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param age.rev The direction of the age axis can be reversed using \code{age.rev=TRUE}.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the plot. Default \code{age.res=500}.
#' @param age.compress Since the bottom and top edges of age-models often show weird flux behaviour, the top and bottom half are shaved off by default (1\%, \code{age.compress=-0.01}).
#' @param flux.lim Limits of the flux axes.
#' @param flux.rev The flux axis can be reversed with \code{flux.rev=TRUE}.
#' @param flux.lab Axis labels. Defaults to \code{flux.lab="flux"}.
#' @param flux.res Resolution or amount of greyscale pixels to cover the flux scale of the plot. Default \code{flux.res=500}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param clip.prob Since some iterations show very high flux values, it is useful to clip the top values. Defaults to \code{clip.prob=0.975}.
#' @param plot.range Plot curves that indicate a probability range.
#' @param prob Probability for the flux ranges; defaults to \code{prob=0.95}.
#' @param range.col Grey seems nice (default \code{range.col=grey(0.5)}).
#' @param range.lty Line type of the flux ranges. Defaults to dashed, \code{range.lty=2}.
#' @param plot.mean Plot the mean fluxes. Default TRUE.
#' @param mean.col Red seems nice.
#' @param mean.lty Line type of the means. Defaults to dashed, \code{range.lty=2}.
#' @param plot.median Plot the median fluxes. Default TRUE
#' @param median.col Blue seems nice.
#' @param median.lty Line type of the medians. Defaults to dashed, \code{range.lty=2}.
#' @param flux.cols Colour of the flux ghost graph. Defaults to greyscales, \code{flux.cols=gray.colors(256, start=darkest, end=0)}.
#' @param dark Any (normalised) flux value above this threshold is set to the threshold; default is \code{dark=0.9}.
#' @param darkest The darkest value. Set to 1 for the darkest colour to be black if using a greyscale.
#' @param rotate.axes The default of plotting calendar year on the horizontal axis and fluxes on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param use.raster Rasters can be aligned or not in the underlying image function. By default, we use \code{use.raster=FALSE}. This takes a bit longer to draw and sometimes causes strange lines owing to anti-aliasing.
#' @param xaxs Limits of the horizontal axis. By default does not extend, \code{xaxs="i"}.
#' @param yaxs Limits of the vertical axis. By default does not extend, \code{yaxs="i"}.
#' @param bty Type of box to plot around the graph. Defaults to L-shaped, \code{bty="l"}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of flux rates, and the underlying greyscales, means, medians and ranges (invisibly).
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   flux.age.ghost(1)
#' }
#' @export
flux.age.ghost <- function(column=1, flux=c(), set=get("info"), coredir=set$coredir, remove.hiatuses=TRUE, age.lab=c(), age.lim=c(), age.rev=FALSE, age.res=500, age.compress=-0.01, flux.lim=c(), flux.rev=FALSE, flux.lab="flux (proxy/yr)", flux.res=500, BCAD=set$BCAD, prob=0.95, clip.prob = 0.99, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, plot.median=TRUE, median.col="blue", median.lty=2, flux.cols=col.scales(256, "white", "black", darkest=darkest), dark=.9, darkest=1, rotate.axes=FALSE, use.raster=FALSE, xaxs="i", yaxs="i", bty="l") {
  if(is.null(flux)) {
    if(!dir.exists(file.path(set$coredir, set$core)))
      stop("please provide the folder where the run's _flux.csv file can be found, e.g., coredir='~/Desktop/Bacon'")
    pf <- read.csv(file.path(set$coredir, set$core, paste0(set$core, "_flux.csv")))
    depths <- pf[,1]
    if(!is.numeric(column) || column < 1 || column > ncol(pf) - 1)
      stop("column out of range")
    proxy  <- pf[,column+1]
  } else { # then we assume that flux is provided as columns
     depths <- flux[,1]
     proxy <- flux[,column+1]
     }
  inside <- which(depths >= min(set$d.min) & depths <= max(set$d.max))
  proxy <- proxy[inside] # remove values outside the core's depth range
  depths <- depths[inside]
  D <- length(depths) # number of remaining depth and proxy slices

  if(is.null(age.lim))
    age.lim <- rev(extendrange(set$ranges[,2:3], f=age.compress)) # remove the top/bottom ends, as they often show strange fluxes
  if(age.rev)
    age.lim <- rev(age.lim)
  age.breaks <- seq(min(age.lim), max(age.lim), length.out=age.res+1)
  age.mids <- 0.5 * (age.breaks[-1] + age.breaks[-length(age.breaks)])

  ages.matrix <- sapply(depths, Bacon.Age.d) # ages of depths d
  acc.matrix <- vapply(depths,
    function(d) accrate.depth(d, set=set, remove.hiatuses=remove.hiatuses), numeric(set$Tr)) # accs of d
  flux.matrix <- sweep(acc.matrix, 2, proxy, FUN=function(a, p) p/a) # g / yr/cm

  flux.on.age <- t(apply( # left half of matrix contains the ages, right half fluxes
    cbind(ages.matrix, flux.matrix), 1,
      function(row) approx(row[1:D], row[(D + 1):(2 * D)],
      xout=age.mids, rule=1, yleft=NA, yright=NA)$y))

  flux.upper <- quantile(flux.on.age, clip.prob, na.rm=TRUE) # remove very high flux values
  flux.on.age[] <- pmin(flux.on.age, flux.upper) # retain the values on the limit
  flux.breaks <- seq(0, flux.upper, length.out=flux.res+1) # make a grid
  flux.mids <- 0.5 * (flux.breaks[-1] + flux.breaks[-length(flux.breaks)])
  fluxes <- matrix(0, nrow=flux.res, ncol=age.res) # flux × age

  for(j in seq_len(age.res)) {
    colj <- na.omit(flux.on.age[,j])
    if(length(colj) > 2) {
      d <- density(colj, from=0, to=flux.upper, n=flux.res)
      fluxes[,j] <- d$y / sum(d$y) # normalise each column
    }
  }

  fluxes <- fluxes / max(fluxes) # greyscale between 0 and 1
  fluxes[fluxes>dark] <- dark # if only a very small spot is the darkest, then increase dark
  fluxes <- fluxes / max(fluxes) # normalise again, so that black=black

  mean.flux <- colMeans(flux.on.age, na.rm=TRUE)
  median.flux <- apply(flux.on.age, 2, median, na.rm=TRUE)
  prob.edges <- c(((1-prob)/2), 1 - ((1-prob)/2))
  rng.flux <- t(apply(flux.on.age, 2, quantile, probs=prob.edges, na.rm=TRUE))

  if(length(age.lab) == 0)
    age.lab <- ifelse(BCAD, "BC/AD", "cal BP")
  if(length(flux.lab) == 0)
     flux.lab <- "flux (proxy/yr)"
  if(length(flux.lim) == 0)
    flux.lim <- c(0, clip.prob*flux.upper)
  if(flux.rev)
    flux.lim <- rev(flux.lim)
  
  if(BCAD && !set$BCAD) {
    age.lim <- rev(rice::calBPtoBCAD(age.lim))
    age.mids <- rev(rice::calBPtoBCAD(age.mids))
    fluxes <- fluxes[,ncol(fluxes):1]
    mean.flux <- rev(mean.flux)
    median.flux <- rev(median.flux)
    rng.flux <- rng.flux[nrow(rng.flux):1,]
  }
  if(!BCAD && set$BCAD) {
    age.lim <- rev(rice::BCADtocalBP(age.lim))
    age.mids <- rev(rice::BCADtocalBP(age.mids))
    fluxes <- fluxes[,ncol(fluxes):1]
    mean.flux <- rev(mean.flux)
    median.flux <- rev(median.flux)
    rng.flux <- rng.flux[nrow(rng.flux):1,]
  }

  if(rotate.axes) {
    plot(0, type="n", xaxs=xaxs, yaxs=yaxs, bty=bty, 
      ylim=age.lim, ylab=age.lab, xlim=flux.lim, xlab=flux.lab)
    image(flux.mids, age.mids, t(t(fluxes)), add=TRUE,
      useRaster=use.raster, col=flux.cols)

    if(plot.mean)
      lines(mean.flux, age.mids, col=mean.col, lty=2)
    if(plot.median)
      lines(median.flux, age.mids, col=median.col, lty=2)
    if(plot.range) {
      lines(rng.flux[,1], age.mids, col=range.col, lty=range.lty)
      lines(rng.flux[,2], age.mids, col=range.col, lty=range.lty)
    }
  } else {
    plot(0, type="n", xaxs=xaxs, yaxs=yaxs, bty=bty, 
      xlim=age.lim, xlab=age.lab, ylim=flux.lim, ylab=flux.lab)
    image(age.mids, flux.mids, t(fluxes), add=TRUE,
      useRaster=use.raster, col=flux.cols)

    if(plot.mean)
      lines(age.mids, mean.flux, col=mean.col, lty=2)
    if(plot.median)
      lines(age.mids, median.flux, col=median.col, lty=2)
    if(plot.range) {
      lines(age.mids, rng.flux[,1], col=range.col, lty=range.lty)
      lines(age.mids, rng.flux[,2], col=range.col, lty=range.lty)
    }
  }
  
  invisible(list(age=age.mids, flux=flux.mids, prob=t(fluxes), 
    medians=cbind(age.mids, median.flux), means=cbind(age.mids, mean.flux), 
    ranges=cbind(age.mids, rng.flux)))
}
