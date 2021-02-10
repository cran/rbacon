

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
#' @author Maarten Blaauw, J. Andres Christen
#' @return all MCMC estimates of accumulation rate of the chosen depth.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   d20 <- accrate.depth(20)
#'   hist(d20)
#'   d20 <- accrate.depth(20, cmyr=TRUE) # to calculate accumulation rates in cm/yr
#'   mean(d20)
#' @export
accrate.depth <- function(d, set=get('info'), cmyr=FALSE) {
  accs.elbows <- set$output[,2:(set$K+1)]
  if(min(set$elbows) <= d && max(set$elbows) >= d)
    accs <- unlist(accs.elbows[max(which(set$elbows <= d))]) else
  #  accs <- set$output[,1+min(which(set$elbows >= d))] else # was 1 + max(which(set$elbows < d))...
      accs <- NA
  if(cmyr) 1/accs else accs
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
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param silent Warn when ages are outside the core's range. Default \code{silent=TRUE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return all MCMC estimates of accumulation rate of the chosen age.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   accrate.a5000 = accrate.age(5000)
#'   plot(accrate.a5000, pch='.')
#'   hist(accrate.a5000)
#' @export
accrate.age <- function(age, set=get('info'), cmyr=FALSE, BCAD=set$BCAD, silent=TRUE) {
   ages <- array(0, dim=c(nrow(set$output), length(set$elbows)))
   for(i in 1:ncol(ages))
     ages[,i] <- Bacon.Age.d(set$elbows[i])

  if(!silent)
    if(age < min(ages) || age > max(ages))
     stop(" Warning, age outside the core's age range!\n")

  accs <- c()
  for(i in 2:ncol(ages)) {
    these <- (ages[,i-1] < age) * (ages[,i] > age)
    if(sum(these) > 0) # age lies within these age-model iterations
      accs <- c(accs, set$output[which(these>0),i]) # was i+1
  }

  if(cmyr)
    accs <- 1/accs
  return(accs)
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
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param prob Probability ranges. Defaults to \code{prob=0.95}.
#' @param plot.range If \code{plot.range=TRUE}, the confidence ranges (two-tailed; half of the probability at each side) are plotted.
#' @param range.col Colour of the confidence ranges.
#' @param range.lty Line type of the confidence ranges.
#' @param plot.mean If \code{plot.mean=TRUE}, the means are plotted.
#' @param mean.col Colour of the mean accumulation rates.
#' @param mean.lty Type of the mean lines.
#' @param rotate.axes The default is to plot the accumulation rates horizontally and the depth vertically (\code{rotate.axes=FALSE}). Change rotate.axes value to rotate axes.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}.
#' @param rev.acc The direction of the accumulation rate axis can be reversed from the default (\code{rev.acc=TRUE}).
#' @author Maarten Blaauw, J. Andres Christen
#' @return A grey-scale plot of accumulation rate against core depth.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#'   layout(1)
#'   accrate.depth.ghost()
#' @export
accrate.depth.ghost <- function(set=get('info'), d=set$elbows, d.lim=c(), acc.lim=c(), d.lab=c(), cmyr=FALSE, acc.lab=c(), dark=1, rgb.scale=c(0,0,0), rgb.res=100, prob=0.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, rotate.axes=FALSE, rev.d=FALSE, rev.acc=FALSE) {
  max.acc <- 0; max.dens <- 0
  acc <- list(); min.rng <- numeric(length(d)); max.rng <- numeric(length(d)); mn.rng <- numeric(length(d))
  for(i in 1:length(d))
    if(length(acc.lim) == 0)
      acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr), from=0) else
        acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr), from=0, to=max(acc.lim))
  for(i in 1:length(d)) {
    max.acc <- max(max.acc, acc[[i]]$x)
    max.dens <- max(max.dens, acc[[i]]$y)
    quants <- quantile(accrate.depth(d[i], set, cmyr=cmyr), c((1-prob)/2, 1-((1-prob)/2)))
    min.rng[i] <- quants[1]
    max.rng[i] <- quants[2]
    mn.rng[i] <- mean(accrate.depth(d[i], set, cmyr=cmyr))
   }

  for(i in 1:length(d)) {
    acc[[i]]$y <- acc[[i]]$y/max.dens
    acc[[i]]$y[acc[[i]]$y > dark] <- dark
  }

  if(length(d.lim) == 0)
    d.lim <- range(d)
  if(length(d.lab) == 0)
    d.lab <- paste("depth (", set$depth.unit, ")", sep="")
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
    plot(0, type="n", xlab=acc.lab, ylab=d.lab, ylim=d.lim, xlim=acc.lim)
    for(i in 2:length(d)) {
      col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0, 1-max(acc[[i-1]]$y), length=rgb.res)) # was acc[[i]]
      image(acc[[i-1]]$x, d[c(i-1, i)], t(1-t(acc[[i-1]]$y)), add=TRUE, col=col) # was acc[[i]]
    }
    if(plot.range) {
      lines(min.rng, d+(set$thick/2), col=range.col, lty=range.lty)
      lines(max.rng, d+(set$thick/2), col=range.col, lty=range.lty)
    }
    if(plot.mean)
      lines(mn.rng, d+(set$thick/2), col=mean.col, lty=mean.lty)
  } else {
      plot(0, type="n", xlab=d.lab, ylab=acc.lab, xlim=d.lim, ylim=acc.lim)
      for(i in 2:length(d)) {
        col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(max(acc[[i-1]]$y), 0, length=rgb.res)) # was acc[[i]]
        image(d[c(i-1, i)], acc[[i-1]]$x, 1-t(acc[[i-1]]$y), add=TRUE, col=col) # was acc[[i]]
        }
      if(plot.range) {
        lines(d+(set$thick/2), min.rng, col=range.col, lty=range.lty)
        lines(d+(set$thick/2), max.rng, col=range.col, lty=range.lty)
        }
    if(plot.mean)
      lines(d+(set$thick/2), mn.rng, col=mean.col, lty=mean.lty)
    }
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
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the plot. Default \code{age.res=400}.
#' @param acc.res Resolution or amount of greyscale pixels to cover the accumulation rate scale plot. Default \code{age.res=400}.
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
#' @return A greyscale plot of accumulation rate against calendar age.
#' @examples
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(age.res=20, d.res=20, d.by=10)
#'   layout(1)
#'   accrate.age.ghost(age.res=20, acc.res=10)
#' }
#' @export
accrate.age.ghost <- function(set=get('info'), age.lim=c(), age.lab=c(), age.res=400, acc.res=200, cutoff=.001, rgb.scale=c(0,0,0), rgb.res=100, prob=.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, acc.lim=c(), acc.lab=c(), BCAD=set$BCAD, cmyr=FALSE, rotate.axes=FALSE, rev.age=FALSE, rev.acc=FALSE, xaxs="i", yaxs="i", bty="l") {

  if(length(age.lim) == 0) 
     age.lim <- extendrange(set$ranges[,5]) # just the mean ages, not the extremes
  if(BCAD)
    age.lim <- 1950 - age.lim
  age.seq <- seq(min(age.lim), max(age.lim), length=age.res)
    
  if(length(acc.lim) == 0)
    acc.lim <-  c(0, 1.05*max(set$output[,2:(1+set$K)])) # maximum accrate in the output
  if(cmyr)
    acc.lim <- 1/acc.lim
  acc.seq <- seq(min(acc.lim), max(acc.lim), length=acc.res)
  breaks <- c(acc.seq, acc.seq[acc.res]+diff(acc.seq[1:2])) # bins of the histogram
    
  z <- array(0, dim=c(age.res, acc.res))
  acc.rng <- array(0, dim=c(age.res, 2))
  acc.mn <- 0
  
  pb <- txtProgressBar(min=0, max=max(1,length(age.seq)-1), style = 3)
  for(i in 1:age.res) {
    setTxtProgressBar(pb, i)
    acc <- accrate.age(age.seq[i], cmyr=cmyr, BCAD=BCAD, silent=TRUE)
    z[i,] <- hist(acc, breaks=breaks, plot=FALSE)$counts
    acc.rng[i,] <- quantile(acc, c((1-prob)/2, 1-((1-prob)/2)))
    acc.mn[i] <- mean(acc)
  }
  message("\n")
  
  z <- z/max(z) # normalise
  z[z<cutoff] <- NA # do not plot very small values
  
  if(BCAD)
    age.lim <- rev(age.lim)
  if(rev.age)
    age.lim <- rev(age.lim)
  if(rev.acc)
    acc.lim <- rev(acc.lim)
  if(length(age.lab) == 0)
    if(BCAD)
      age.lab <- "BC/AD" else
        age.lab <- "cal BP"
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste0("accumulation rate (", set$depth.unit, "/", set$age.unit, ")") else
        acc.lab <- paste0("accumulation rate (", set$age.unit, "/", set$depth.unit, ")")

  cols <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0, 1, length=rgb.res))

  if(rotate.axes) {
    plot(0, type="n", ylim=age.lim, ylab=age.lab, xlim=acc.lim, xlab=acc.lab, yaxs=xaxs, xaxs=yaxs)
    image(acc.seq, age.seq, t(z), col=cols, add=TRUE)
      if(plot.range) {
        lines(acc.rng[,1], age.seq, pch=".", col=range.col, lty=range.lty)
        lines(acc.rng[,2], age.seq, pch=".", col=range.col, lty=range.lty)
      }
    if(plot.mean) 
      lines(acc.mn, age.seq, col=mean.col, lty=mean.lty)
  } else {
       plot(0, type="n", xlim=age.lim, xlab=age.lab, ylim=acc.lim, ylab=acc.lab, xaxs=xaxs, yaxs=yaxs)
       image(age.seq, acc.seq, t(t(z)), col=cols, add=TRUE)
       if(plot.range) {
          lines(age.seq, acc.rng[,1], pch=".", col=range.col, lty=range.lty)
          lines(age.seq, acc.rng[,2], pch=".", col=range.col, lty=range.lty)
        }
      if(plot.mean) 
        lines(age.seq, acc.mn, col=mean.col, lty=mean.lty)
    }

  if(length(bty) > 0)
    box(bty=bty)
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
#' @param upper Maximum flux rates to plot. Defaults to the upper 99\%; \code{upper=0.99}.
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param dark The darkest grey value is \code{dark=1} by default; lower values will result in lighter grey but \code{values >1} are not allowed.
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
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   flux.age.ghost(1)
#' }
#' @export
flux.age.ghost <- function(proxy=1, age.lim=c(), yr.lim=age.lim, age.res=200, yr.res=age.res, set=get('info'), flux=c(), plot.range=TRUE, prob=.8, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, flux.lim=c(), flux.lab="flux", upper=.95, rgb.scale=c(0,0,0), rgb.res=100, dark=set$dark, BCAD=set$BCAD, age.lab=c(), yr.lab=age.lab, rotate.axes=FALSE, rev.flux=FALSE, rev.age=FALSE, rev.yr=rev.age) {
  if(length(flux) == 0) { # then read a .csv file, expecting data in columns with headers
    flux <- read.csv(paste(set$coredir, set$core, "/", set$core, "_flux.csv", sep=""))
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
    ages <- as.numeric(set$output[i,1:(ncol(set$output)-1)]) # 1st step to calculate ages for each set$elbows
    ages <- c(ages[1], ages[1]+set$thick * cumsum(ages[2:length(ages)])) # now calculate the ages for each set$elbows
    ages.d <- approx(ages, c(set$elbows, max(set$elbows)+set$thick), age.seq, rule=1)$y # find the depth belonging to each age.seq, NA if none
    ages.i <- floor(approx(ages, (length(set$elbows):0)+1, age.seq, rule=2)$y) # find the column belonging to each age.seq
    flux.d <- approx(flux[,1], flux[,2], ages.d)$y # interpolate flux (in depth) to depths belonging to each age.seq
    fluxes[i,] <- flux.d / as.numeric(set$output[i,(1+ages.i)]) # (amount / cm^3) / (yr/cm) = amount * cm-2 * yr-1
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
    plot(0, type="n", ylim=age.lim, ylab=age.lab, xlim=flux.lim, xlab=flux.lab) else
      plot(0, type="n", xlim=age.lim, xlab=age.lab, ylim=flux.lim, ylab=flux.lab)
  min.rng <- numeric(length(age.seq)); max.rng <- numeric(length(age.seq)); mn.rng <- numeric(length(age.seq))
  for(i in 2:length(age.seq)) {
    tmp <- fluxes[!is.na(fluxes[,i]),i] # all fluxes that fall at the required age.seq age
    rng <- quantile(tmp, c((1-prob)/2, 1-((1-prob)/2)))
    min.rng[i] <- rng[1]
    max.rng[i] <- rng[2]
    mn.rng[i] <- mean(tmp)
    if(length(tmp[tmp>=0]) > 2) {
      flux.hist <- density(tmp, from=0, to=max(flux.lim))
      flux.hist$y <- flux.hist$y - min(flux.hist$y) # no negative fluxes
      col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0, max(flux.hist$y)/max.dens, length=rgb.res))
      if(rotate.axes)
        image(flux.hist$x, age.seq[c(i-1,i)], matrix(flux.hist$y), add=TRUE, col=col) else
          image(age.seq[c(i-1,i)], flux.hist$x, t(matrix(flux.hist$y)), add=TRUE, col=col)
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
      lines(mn.rng, age.seq, col=mean.col, lty=mean.lty) else
       lines(age.seq, mn.rng, col=mean.col, lty=mean.lty)
}
