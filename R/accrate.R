

## Accumulation rate calculations
#' should take into account hiatuses
#' @name accrate.depth
#' @title Obtain estimated accumulation rates as for any depth of a core.
#' @description Obtain accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any depth of a core.
#' @details Considering accumulation rates is crucial for age-depth modelling, and even more so if they are subsequently used for calculating proxy
#' influx values, or interpreted as proxy for environmental change such as carbon accumulation.
#' Bacon deals explicitly with accumulation rate and its variability through defining prior distributions.
#' This function obtains accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations
#' for any depth of a core. Deals with only 1 depth at a time. See also \code{accrate.age}.
#' @param d The depth for which accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
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
accrate.depth <- function(d, set=get('info'), cmyr=FALSE, na.rm=FALSE, inversion.threshold=1e-6) {
  accs.elbows <- set$output[,2:(set$K+1)]
  if(min(set$elbows) <= d && max(set$elbows) >= d)
    accs <- unlist(accs.elbows[max(which(set$elbows <= d))]) else
      accs <- NA
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
#' @title Obtain estimated accumulation rates for any age of a core.
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

   # these two lines do the same as the loop below, 
   #   but at the same speed and values outside the ages do not get NAs
   # col_indices <- rowSums(ages <= age) + 1 
   # accs <- set$output[cbind(seq_len(nrow(ages)), col_indices)]

  accs <- rep(NA_real_, nrow(ages)) # suggested by henningte on github
  for(i in 2:ncol(ages)) {
    these <- (ages[,i-1] < age) & (ages[,i] > age)
    if(sum(these) > 0) # age lies within these age-model iterations
      accs[which(these>0)] <- set$output[which(these>0),i] # Jan 2023
  }

  if(na.rm)
    accs <- accs[!is.na(accs)]
  if(cmyr)
    accs <- 1/accs

  return(accs)
}



#' @name accrate.depth.summary
#' @title Provide a summary of the estimated accumulation rates for any depth of a core.
#' @description Obtain a summary (95\% range, 68\% range, median, mean) of the accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any depth of a core.
#' @param d The depth for which accumulation rates need to be returned.
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
#'   accrate.depth.summary(20)
#' }
#' @export
accrate.depth.summary <- function(d, set=get('info'), cmyr=FALSE, na.rm=FALSE, probs=c(.025, .16, .84, .975, .5)) {
  if(length(d) > 1)
    stop("can handle one depth at a time only")
  accs <- accrate.depth(d, set, cmyr, na.rm)
  qu <- quantile(accs, probs, na.rm=na.rm)
  mn <- mean(accs, na.rm=na.rm)
  names(mn) <- "mean"
  return(c(qu, mn))
}



#' @name accrate.age.summary
#' @title Provide a summary of the estimated accumulation rates for any age of a core.
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



#' @name accrates.core
#' @title Provide a summary of the estimated accumulation rates for a range of core depths
#' @description Obtain a summary (95\% range, 68\% range, median, mean) of the accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for a range of depths of a core, and optionally write this as a file to the core directory (ending in '_accrates.txt').
#' @param dseq The sequence of depths for which accumulation rates need to be returned. Defaults to whatever info$dseq is, which most often is a sequence from the top to the bottom of the core at 1 cm increments.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm.
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
accrates.core <- function(dseq=c(), set=get('info'), cmyr=FALSE, na.rm=FALSE, probs=c(.025, .16, .84, .975, .5), round=2, write=TRUE, sep="\t") {
  if(length(dseq) == 0)
    dseq <- set$depths
  mysummary <- function(dseq)
    accrate.depth.summary(dseq, set, cmyr, na.rm, probs)
  allaccs <- t(sapply(dseq, mysummary))
  allaccs <- round(allaccs,round)
  
  if(write) {
    fl <- paste0(set$coredir, set$core, "/", set$core, "_", set$K, "_accrates.txt")
    message("writing the accumulation rate summary to ", fl)
	write.table(allaccs, fl, sep=sep, quote=FALSE, row.names=FALSE)
  } 
  invisible(allaccs)
}




#' @name accrate.depth.ghost
#' @title Plot modelled accumulation rates against the depths of a core.
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
#' @param dark The darkest grey value is dark=1 by default; lower values will result in lighter grey but values >1 are not advised.
#' @param cutoff Point below which colours will no longer be printed. Default \code{cutoff=0.001}.
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
accrate.depth.ghost <- function(set=get('info'), d=set$elbows, d.lim=c(), acc.lim=c(), d.lab=c(), cmyr=FALSE, acc.lab=c(), dark=1, cutoff=0.001, rgb.scale=c(0,0,0), rgb.res=100, prob=0.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, plot.median=TRUE, median.col="blue", median.lty=2,  rotate.axes=FALSE, rev.d=FALSE, rev.acc=FALSE, xaxs="r", yaxs="r", bty="l", remove.laststep=TRUE) {
  max.acc <- 0; max.dens <- 0
  acc <- list(); min.rng <- numeric(length(d)); max.rng <- numeric(length(d)); mean.rng <- numeric(length(d)); median.rng <- numeric(length(d))
  for(i in 1:length(d))
    if(length(acc.lim) == 0)
      acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr), from=0) else
        acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr), from=0, to=max(acc.lim))
  for(i in 1:length(d)) {
    max.acc <- max(max.acc, acc[[i]]$x)
    max.dens <- max(max.dens, acc[[i]]$y)
    accs <- accrate.depth(d[i], set, cmyr=cmyr)
    quants <- quantile(accs, c((1-prob)/2, 1-((1-prob)/2)))
    min.rng[i] <- quants[1]
    max.rng[i] <- quants[2]
    mean.rng[i] <- mean(accs)
    median.rng[i] <- median(accs)
   }
  stored <- cbind(d, min.rng, max.rng, median.rng, mean.rng)
  colnames(stored) <- c("depth", "min.rng", "max.rng", "median", "mean")

  for(i in 1:length(d)) {
    acc[[i]]$y <- acc[[i]]$y/(dark*max.dens)
    acc[[i]]$y[acc[[i]]$y > 1] <- 1 # set "dark" to black
    acc[[i]]$y[acc[[i]]$y < cutoff] <- NA # do not plot too light/small values
  }

  if(length(d.lim) == 0)
    d.lim <- range(d)
  if(length(d.lab) == 0)
    d.lab <- paste0("depth (", set$depth.unit, ")")
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste0("accumulation rate (", set$depth.unit, "/", set$age.unit, ")") else
        acc.lab <- paste0("accumulation rate (", set$age.unit, "/", set$depth.unit, ")")

  if(rev.d)
    d.lim <- rev(d.lim)
  if(length(acc.lim) == 0)
    acc.lim <- c(0, max.acc)
  if(rev.acc)
    acc.lim <- rev(acc.lim)

  if(rotate.axes) {
    plot(0, type="n", xlab=acc.lab, ylab=d.lab, ylim=d.lim, xlim=acc.lim, bty="n", xaxs=xaxs, yaxs=yaxs)
    for(i in 2:length(d)) {
      accs <- acc[[i-1]]
      col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(max(accs$y[!is.na(accs$y)]), 0, length=rgb.res)) # was acc[[i]]
      ghost.mirror(accs$x, d[c(i-1, i)], t(1-t(accs$y)), col=col) # was acc[[i]]
    }
    if(plot.range) {
      lines(min.rng, d, type="s", col=range.col, lty=range.lty)
      lines(max.rng, d, type="s", col=range.col, lty=range.lty)
    }
    if(plot.mean)
      lines(mean.rng, d, type="s", col=mean.col, lty=mean.lty)
    if(plot.median)
      lines(median.rng, d, type="s", col=median.col, lty=median.lty)
    if(remove.laststep)
      abline(h=min(set$elbows), col="white", lwd=2)
  } else {
      plot(0, type="n", xlab=d.lab, ylab=acc.lab, xlim=d.lim, ylim=acc.lim, bty="n", xaxs=xaxs, yaxs=yaxs)
      for(i in 2:length(d)) {
        accs <- acc[[i-1]]
        col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(max(accs$y[!is.na(accs$y)]), 0, length=rgb.res)) # was acc[[i]]
        ghost.mirror(d[c(i-1, i)], accs$x, 1-t(accs$y), col=col) # was acc[[i]]
      }
      if(plot.range) {
        lines(d, min.rng, type="s", col=range.col, lty=range.lty, pch=NA)
        lines(d, max.rng, type="s", col=range.col, lty=range.lty, pch=NA)
        }
    if(plot.mean)
      lines(d, mean.rng, type="s", col=mean.col, lty=mean.lty)
    if(plot.median)
      lines(d, median.rng, type="s", col=median.col, lty=median.lty)
    if(remove.laststep)
      abline(v=max(set$elbows), col="white", lwd=1.5)
    }

  box(bty=bty)  
  invisible(stored)
}



#' @name accrate.age.ghost
#' @title Plot a core's accumulation rates against calendar time.
#' @description Plot a grey-scale representation of a core's estimated accumulation rates against time.
#' @details Calculating accumulation rates against calendar age will take some time to calculate, and might show unexpected
#' rates around the core's maximum ages (only a few of all age-model iterations will reach such ages and they will tend to have
#'  modelled accumulation rates for the lower depths much lower than the other iterations). Axis limits for accumulation rates
#'   are estimated automatically, however upper limits can be very variable (and thus hard to predict) if calculated in \code{cm/yr}.
#'  Therefore you might want to manually adapt the axis limits after plotting with default settings (e.g., \code{acc.lim=c(0,1)}). See also \code{accrate.depth.ghost}, \code{accrate.depth} and \code{accrate.age}.
#' The grey-scale reconstruction around the oldest ages of any reconstruction often indicates very low accumulation rates.
#' This is due to only some MCMC iterations reaching those old ages, and these iterations will have modelled very slow accumulation rates.
#' Currently does not deal well with hiatuses, so do not interpret accumulation rates close to depths with inferred hiatuses.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param kcal Use kcal BP. Default is \code{kcal=FALSE}.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the plot. Default \code{age.res=400}.
#' @param acc.res Resolution or amount of greyscale pixels to cover the accumulation rate scale plot. Default \code{age.res=400}.
#' @param cutoff Point below which colours will no longer be printed. Default \code{cutoff=0.001}.
#' @param dark The darkest grey value is dark=1 by default; lower values will result in lighter grey but values >1 are not advised.
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
accrate.age.ghost <- function(set=get('info'), age.lim=c(), age.lab=c(), kcal=FALSE, age.res=400, acc.res=200, cutoff=.001, dark=1, rgb.scale=c(0,0,0), rgb.res=100, prob=.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, plot.median=TRUE, median.col="blue", median.lty=2, acc.lim=c(), acc.lab=c(), BCAD=set$BCAD, cmyr=FALSE, rotate.axes=FALSE, rev.age=FALSE, rev.acc=FALSE, xaxs="i", yaxs="i", bty="l") {
  if(length(age.lim) == 0) 
     age.lim <- extendrange(set$ranges[,5]) # just the mean ages, not the extremes
  if(set$BCAD) # was set$BCAD
    age.lim <- BCADtocalBP(age.lim) # work with cal BP internally
  age.seq <- seq(min(age.lim), max(age.lim), length=age.res)
    
  if(length(acc.lim) == 0) {
    acc.lim <- c(0, 1.05*max(set$output[,2:(1+set$K)])) # maximum accrate in the output
    if(cmyr)
      acc.lim <- 1/acc.lim
    acc.lim[is.infinite(acc.lim)] <- 0  
  }  
  acc.seq <- seq(min(acc.lim, na.rm=TRUE), max(acc.lim, na.rm=TRUE), length=acc.res)
  
  z <- array(0, dim=c(age.res, acc.res))
  acc.rng <- array(NA, dim=c(age.res, 2))
  acc.mean <- rep(NA, age.res); acc.median <- acc.mean

  # speed things up by not repeatedly calculating ages in accrate.age
  ages <- array(0, dim=c(nrow(set$output), length(set$elbows)))
  for(i in 1:ncol(ages))
    ages[,i] <- Bacon.Age.d(set$elbows[i], BCAD=FALSE)

  pb <- txtProgressBar(min=0, max=max(1,length(age.seq)-1), style = 3)
  for(i in 1:age.res) {
    setTxtProgressBar(pb, i)
    acc <- accrate.age(age.seq[i], cmyr=cmyr, ages=ages, silent=TRUE, BCAD=FALSE)
    acc <- acc[!is.na(acc)]
    if(length(acc[!is.na(acc)]) > 1) {
      z[i,] <- density(acc, from=min(acc.lim, na.rm=TRUE), to=max(acc.lim, na.rm=TRUE), n=acc.res)$y
      acc.rng[i,] <- quantile(acc, c((1-prob)/2, 1-((1-prob)/2)))
      acc.mean[i] <- mean(acc)
      acc.median[i] <- median(acc)
    }
  }
  message("\n")
  stored <- cbind(age.seq, acc.rng[,1], acc.rng[,2], acc.median, acc.mean)
  colnames(stored) <- c("ages", "min.rng", "max.rng", "median", "mean")

  z <- z/(dark*max(z)) # normalise, set dark to black
  z[z>1] <- 1 # avoid values > 1
  z[z<cutoff] <- NA # do not plot very small/light greyscale values
  
  if(rev.age)
    age.lim <- rev(age.lim)
  if(rev.acc)
    acc.lim <- rev(acc.lim)
  if(length(age.lab) == 0)
    if(BCAD)
      age.lab <- "BC/AD" else
        age.lab <- ifelse(kcal, "kcal BP", "cal BP")
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste0("accumulation rate (", set$depth.unit, "/", set$age.unit, ")") else
        acc.lab <- paste0("accumulation rate (", set$age.unit, "/", set$depth.unit, ")")

  cols <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0, 1, length=rgb.res))

  if(rotate.axes) {
    yaxt <- ifelse(kcal || BCAD, "n", "s")
    plot(0, type="n", ylim=age.lim, ylab=age.lab, xlim=acc.lim, xlab=acc.lab, yaxs=xaxs, xaxs=yaxs, yaxt=yaxt, bty="n")
    if(BCAD)
      axis(2, pretty(age.lim), labels=calBPtoBCAD(pretty(age.lim))) else
        if(kcal)
          axis(2, pretty(age.lim), labels=pretty(age.lim)/1e3)
    ghost.mirror(acc.seq, age.seq, t(z), col=cols)
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
      ghost.mirror(age.seq, acc.seq, z, col=cols)
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
#' @title Plot flux rates for proxies.
#' @description Plot grey-scale representation of estimated flux rates for proxies against calendar age.
#' @details To plot flux rates (e.g. pollen grains/cm2/yr) as greyscales,
#' provide a plain text file with headers and the data in columns separated by commas, ending in '_flux.csv'
#' and saved in your core's folder. The first column should contain the depths, and the next columns should contain
#' the proxy concentration values (leaving missing values empty). Then type for example \code{flux.age.ghost(1)} to plot the
#' flux values for the first proxy in the .csv file. Instead of using a _flux.csv file, a flux variable can also be defined
#'  within the R session (consisting of depths and their proxy concentrations in two columns). Then provide the name of this variable, e.g.: \code{flux.age.ghost(flux=flux1)}.
#' See Bacon_runs/MSB2K/MSB2K_flux.csv for an example.
#' @param proxy Which proxy to use (counting from the column number in the .csv file after the depths column).
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param yr.lim Deprecated - use age.lim instead
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the plot. Default \code{age.res=200}.
#' @param yr.res Deprecated - use age.res instead
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param flux Define a flux variable within the R session (consisting of depths and their proxy concentrations in two columns) and provide the name of this variable, e.g.:
#' \code{flux.age.ghost(flux=flux1)}. If left empty (\code{flux=c()}), a flux file is expected (see \code{proxy}).
#' @param plot.range Plot curves that indicate a probability range, at resolution of yr.res.
#' @param prob Probability range, defaults to \code{prob=0.8} (10 \% at each side).
#' @param range.col Red seems nice.
#' @param range.lty Line type of the confidence ranges.
#' @param flux.lim Limits of the flux axes.
#' @param flux.lab Axis labels. Defaults to \code{flux.lab="flux"}.
#' @param plot.mean Plot the mean fluxes.
#' @param mean.col Red seems nice.
#' @param mean.lty Line type of the means.
#' @param plot.median Plot the median fluxes.
#' @param median.col Blue seems nice.
#' @param median.lty Line type of the medians.
#' @param upper Maximum flux rates to plot. Defaults to the upper 99\%; \code{upper=0.99}.
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param dark The darkest grey value is \code{dark=1} by default; lower values will result in lighter grey but \code{values >1} are not allowed.
#' @param cutoff Point below which colours will no longer be printed. Default \code{cutoff=0.001}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.lab Deprecated - use age.lab instead
#' @param rotate.axes The default of plotting calendar year on the horizontal axis and fluxes on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.flux The flux axis can be reversed with \code{rev.flux=TRUE}.
#' @param rev.age The direction of the age axis can be reversed using \code{rev.age=TRUE}.
#' @param rev.yr Deprecated - use rev.age instead
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of flux rates.
#' @examples
#' \dontrun{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   flux.age.ghost(1)
#' }
#' @export
flux.age.ghost <- function(proxy=1, age.lim=c(), yr.lim=age.lim, age.res=200, yr.res=age.res, set=get('info'), flux=c(), plot.range=TRUE, prob=.8, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, plot.median=TRUE, median.col="blue", median.lty=2, flux.lim=c(), flux.lab=expression("flux (g cm"^-1*" yr"^-1*")"), upper=.95, rgb.scale=c(0,0,0), rgb.res=100, dark=set$dark, cutoff=0.001, BCAD=set$BCAD, age.lab=c(), yr.lab=age.lab, rotate.axes=FALSE, rev.flux=FALSE, rev.age=FALSE, rev.yr=rev.age) {
  if(length(flux) == 0) { # then read a .csv file, expecting data in columns with headers
    flux <- read.csv(paste0(set$coredir, set$core, "/", set$core, "_flux.csv"))
    flux <- cbind(flux[,1], flux[,1+proxy])
      isNA <- is.na(flux[,2])
      flux <- flux[!isNA,]
  }
  if(length(age.lim) == 0) {
    min.age <- min(set$ranges[,2])
    max.age <- max(set$ranges[,3])
    age.lim <- c(min.age, max.age)
  } else {
      min.age <- min(age.lim)
      max.age <- max(age.lim)
    }

  age.seq <- seq(min(min.age, max.age), max(min.age, max.age), length=age.res)
  fluxes <- array(NA, dim=c(nrow(set$output), length(age.seq)))
  for(i in 1:nrow(set$output)) {
    #setTxtProgressBar(pb, i)
    ages <- as.numeric(set$output[i,1:(ncol(set$output)-1)]) # 1st step to calculate ages for each set$elbows
    ages <- c(ages[1], ages[1]+set$thick * cumsum(ages[2:length(ages)])) # now calculate the ages for each set$elbows
    ages.d <- approx(ages, c(set$elbows, max(set$elbows)+set$thick), age.seq, rule=1)$y # find the depth belonging to each age.seq, NA if none
    ages.i <- floor(approx(ages, (length(set$elbows):0)+1, age.seq, rule=2)$y) # find the column belonging to each age.seq
    flux.d <- approx(flux[,1], flux[,2], ages.d, rule=1)$y # interpolate flux (in depth) to depths belonging to each age.seq
    fluxes[i,] <- flux.d / as.numeric(set$output[i,(1+ages.i)]) # (amount / cm^3) / (yr/cm) = amount * cm-2 * yr-1
    fluxes[is.na(fluxes)] <- 0
  }
  message("\n")
  if(length(flux.lim) == 0)
    flux.lim <- c(0, quantile(fluxes[!is.na(fluxes)], upper))
  max.dens <- 0

  for(i in 1:length(age.seq)) {
    tmp <- fluxes[!is.na(fluxes[,i]),i] # all fluxes that fall at the required age.seq age
    if(length(tmp) > 0)
      max.dens <- max(max.dens, density(tmp, from=0, to=max(flux.lim))$y)
  }

  if(length(age.lim) == 0)
    age.lim <- range(age.seq)
  if(length(age.lab) == 0)
    age.lab <- ifelse(BCAD, "BC/AD", "cal BP")
  if(rotate.axes)
    plot(0, type="n", ylim=age.lim, ylab=age.lab, xlim=flux.lim, xlab=flux.lab, yaxt="n") else
      plot(0, type="n", xlim=age.lim, xlab=age.lab, ylim=flux.lim, ylab=flux.lab, xaxt="n")
  if(BCAD && !set$BCAD) {
    if(rotate.axes)
      axis(2, pretty(age.lim), labels=calBPtoBCAD(pretty(age.lim))) else
        axis(1, pretty(age.lim), labels=calBPtoBCAD(pretty(age.lim)))
  } else
      ifelse(rotate.axes, axis(2), axis(1))

  min.rng <- numeric(length(age.seq)); max.rng <- numeric(length(age.seq)); mean.rng <- numeric(length(age.seq)); median.rng <- numeric(length(age.seq))
  for(i in 2:length(age.seq)) {
    tmp <- fluxes[!is.na(fluxes[,i]),i] # all fluxes that fall at the required age.seq age
    rng <- quantile(tmp, c((1-prob)/2, 1-((1-prob)/2)))
    min.rng[i] <- rng[1]
    max.rng[i] <- rng[2]
    mean.rng[i] <- mean(tmp)
    median.rng[i] <- median(tmp)
    if(length(tmp[tmp>=0]) > 2) {
      flux.hist <- density(tmp, from=0, to=max(flux.lim))
      flux.hist$y <- flux.hist$y - min(flux.hist$y) # no negative fluxes
      flux.hist$y <- flux.hist$y / (dark*max.dens) # normalise
      flux.hist$y[flux.hist$y > 1] <- 1 # no values > 1
      flux.hist$y[flux.hist$y < cutoff] <- NA # do not plot very small/light greyscale values
      col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3],
        seq(0, max(flux.hist$y[!is.na(flux.hist$y)]), length=rgb.res))
      if(rotate.axes)
        ghost.mirror(flux.hist$x, age.seq[c(i-1,i)], matrix(flux.hist$y), col=col) else
          ghost.mirror(age.seq[c(i-1,i)], flux.hist$x, t(matrix(flux.hist$y)), col=col)
    }
  }

  if(plot.range)
    if(rotate.axes) {
      lines(min.rng, age.seq, col=range.col, lty=range.lty)
      lines(max.rng, age.seq, col=range.col, lty=range.lty)
  } else {
      lines(age.seq, min.rng, col=range.col, lty=range.lty)
      lines(age.seq, max.rng, col=range.col, lty=range.lty)
    }
  if(plot.mean)
    if(rotate.axes)
      lines(mean.rng, age.seq, col=mean.col, lty=mean.lty) else
       lines(age.seq, mean.rng, col=mean.col, lty=mean.lty)
  if(plot.median)
    if(rotate.axes)
      lines(median.rng, age.seq, col=median.col, lty=median.lty) else
       lines(age.seq, median.rng, col=median.col, lty=median.lty)
}
