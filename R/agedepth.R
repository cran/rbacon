#' @name agedepth 
#' @title Plot an age-depth model
#' @description Plot the age-depth model of a core.
#' @details After loading a previous run, or after running either the \link{scissors} or \link{thinner} command, plot the age-model
#' again using the command \code{agedepth()}. 
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}. See also \code{unit}.
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth).
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param d.by Depth intervals at which ages are calculated. Default 1. Alternative depth intervals can be provided using, e.g., d.\code{by=0.5}.
#' @param depths By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}. Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}.
#' @param depths.file By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.  
#' If \code{depths.file=TRUE}, Bacon will read a file containing the depths for which you require ages. 
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir}, 
#' and its name should start with the core's name and end with '_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param yr.min Minimum calendar age of the age-depth plot.
#' @param yr.max Maximum calendar age of the age-depth plot.
#' @param dark Darkness of the greyscale age-depth model. By default, the darkest grey value is calculated as 10 times the height of the lowest-precision age estimate \code{dark=c()}. Lower values will result in lighter grey but values >1 are not allowed.
#' @param prob Confidence interval to report (between 0 and 1, default 0.95 or 95\%).
#' @param rounded Rounding of years. Default is to round to single years.
#' @param d.res Resolution or amount of greyscale pixels to cover the depth scale of the age-model plot. Default \code{d.res=200}.
#' @param yr.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
#' @param grey.res Grey-scale resolution of the age-depth model. Default \code{grey.res=100}.
#' @param rotate.axes By default, the age-depth model is plotted with the depths on the horizontal axis and ages on the vertical axis. This can be changed with \code{rotate.axes=TRUE}.
#' @param rev.yr The direction of the age axis, which can be reversed using \code{rev.yr=TRUE}.
#' @param rev.d The direction of the depth axis, which can be reversed using \code{rev.d=TRUE}.
#' @param unit Depth units, default \code{unit="cm"}.
#' @param maxcalc Number of depths to calculate ages for. If this is more than \code{maxcalc=500}, a warning will be shown that calculations will take time.
#' @param height The maximum heights of the distributions of the dates on the plot. See also \code{normalise.dists}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
#' @param panels Divide the graph panel. Defaults to 1 graph per panel, \code{panels=layout(1)}. To avoid dividing into panels, use \code{panels=c()}.
#' @param plot.range Whether or not to plot the curves showing the confidence ranges of the age-model. Defaults to (\code{plot.range=TRUE}).
#' @param range.col The colour of the curves showing the confidence ranges of the age-model. Defaults to medium grey (\code{range.col=grey(0.5)}).
#' @param range.lty The line type of the curves showing the confidence ranges of the age-model. Defaults to \code{range.lty=12}.
#' @param mn.col The colour of the mean age-depth model: default \code{mn.col="red"}.
#' @param mn.lty The line type of the mean age-depth model. Default \code{mn.lty=12}.
#' @param med.col The colour of the median age-depth model: not drawn by default \code{med.col=NA}.
#' @param med.lty The line type of the median age-depth model. Default \code{med.lty=12}.
#' @param C14.col The colour of the calibrated ranges of the dates. Default is semi-transparent blue: \code{C14.col=rgb(0,0,1,.35)}.
#' @param C14.border The colours of the borders of calibrated 14C dates. Default is semi-transparent dark blue: \code{C14.border=rgb(0, 0, 1, 0.5)}.
#' @param cal.col The colour of the non-14C dates. Default is semi-transparent blue-green: \code{cal.col=rgb(0,.5,.5,.35)}.
#' @param cal.border The colour of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{cal.border=rgb(0,.5,.5,.5)}. Not used by default.
#' @param dates.col As an alternative to colouring dates based on whether they are 14C or not, sets of dates can be coloured as, e.g., \code{dates.col=colours()[2:100]}.
#' @param hiatus.col The colour of the depths of any hiatuses. Default \code{hiatus.col=grey(0.5)}.
#' @param hiatus.lty The line type of the depths of any hiatuses. Default \code{hiatus.lty=12}.
#' @param greyscale The function to produce a coloured representation of all age-models. Defaults to grey-scales: \code{greyscale=function(x) grey(1-x)}. 
#' @param slump.col Colour of slumps. Defaults to \code{slump.col=grey(0.8)}.
#' @param normalise.dists By default, the distributions of more precise dates will cover less time and will thus peak higher than less precise dates. This can be avoided by specifying \code{normalise.dists=FALSE}.
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), \code{cc=3} for SHCal13 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param title The title of the age-depth model is plotted on the main panel. By default this is the core's name. To leave empty: \code{title=""}. 
#' @param title.location Location of the title. Default \code{title.location='topleft'}.
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param bty Type of box to be drawn around plots (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param mar Plot margins (amount of white space along edges of axes 1-4). Default \code{mar=c(3,3,1,1)}.
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.5, .7, .0)}.
#' @param xaxs Extension of x-axis. By default, add some extra white-space at both extremes (\code{xaxs="r"}). See ?par for other options. 
#' @param yaxs Extension of y-axis. By default, add some extra white-space at both extremes (\code{yaxs="r"}). See ?par for other options. 
#' @param xaxt Whether or not to plot the x-axis. Can be used to adapt axes after a plot. See ?par for other options. 
#' @param yaxt Whether or not to plot the y-axis. Can be used to adapt axes after a plot. See ?par for other options. 
#' @param plot.pdf Produce a pdf file of the age-depth plot.
#' @param model.only By default, panels showing the MCMC iterations and the priors and posteriors for accumulation rate and memory are plotted above the main age-depth model panel. This can be avoided by supplying \code{model.only=TRUE}. Note however that this removes relevant information to evaluate the age-depth model, so we do recommend to present age-models together with these upper panels.
#' @param dates.only By default, the age-depth model is plotted on top of the dates. This can be avoided by supplying \code{dates.only=TRUE}.
#' @param talk Provide a summary of the age ranges after producing the age-depth model graph; default \code{talk=FALSE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of the age-depth model, and estimated ages incl. confidence ranges for each depth.
#' @examples 
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())}
#'   agedepth(yr.res=50, d.res=50, d.by=10)
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   agedepth()
#' }
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
agedepth <- function(set=get('info'), BCAD=set$BCAD, unit="cm", d.lab=c(), yr.lab=c(), d.min=set$d.min, d.max=set$d.max, d.by=set$d.by, depths=set$depths, depths.file=FALSE, yr.min=c(), yr.max=c(), dark=c(), prob=set$prob, rounded=0, d.res=400, yr.res=400, date.res=100, grey.res=100, rotate.axes=FALSE, rev.yr=FALSE, rev.d=FALSE, maxcalc=500, height=15, mirror=TRUE, up=TRUE, cutoff=.001, plot.range=TRUE, panels=layout(1), range.col=grey(.5), range.lty="12", mn.col="red", mn.lty="12", med.col=NA, med.lty="12", C14.col=rgb(0,0,1,.35), C14.border=rgb(0,0,1,.5), cal.col=rgb(0,.5,.5,.35), cal.border=rgb(0,.5,.5,.5), dates.col=c(), hiatus.col=grey(0.5), hiatus.lty="12", greyscale=grey(seq(1, 0, length=grey.res)), slump.col=grey(0.8), normalise.dists=TRUE, cc=set$cc, title=set$core, title.location="topleft", after=set$after, bty="l", mar=c(3,3,1,1), mgp=c(1.5,.7,.0), xaxs="r", yaxs="r", xaxt="s", yaxt="s", plot.pdf=FALSE, dates.only=FALSE, model.only=FALSE, talk=FALSE) {
# Load the output, if it exists
  outp <- paste0(set$prefix, ".out")
  if(file.exists(outp))
    set <- .Bacon.AnaOut(outp, set)

  par(bty=bty, mar=mar, mgp=mgp, yaxs="i")
  if(model.only) 
    panels else { # layout(1) can mess things up if plotting within an existing panel...
      pn <- c(1:3, rep(4,3))
      if(!is.na(set$hiatus.depths)[1])
        if(is.na(set$boundary)[1])
          pn <- c(1:4, rep(5,4))  
      layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7)) 
      .PlotLogPost(set, 0, set$Tr) # convergence information
      .PlotAccPost(set)
      .PlotMemPost(set, set$core, set$K, "", set$mem.strength, set$mem.mean, ds=1, thick=set$thick)
      if(!is.na(set$hiatus.depths[1]))
        if(is.na(set$boundary)[1])
          .PlotHiatusPost(set, set$hiatus.max)
  } 


  # calculate and plot the ranges and 'best' estimates for each required depth
  if(length(depths) > 0)
    d <- sort(depths) else
      d <- seq(d.min, d.max, by=d.by) 
      
  if(length(d) > maxcalc)
    cat("Warning, this will take quite some time to calculate. I suggest increasing d.by to, e.g.", 10*d.by, "\n") # was set$d.by

  # cosmetic; to avoid very steep plotted curves of mean and ranges across a hiatus
  for(i in set$hiatus.depths) 
    d <- sort(unique(c(i+after, i, d)))
  for(i in set$slump)
    d <- sort(unique(c(i+after, i, d)))

  if(talk)
    cat("Calculating age ranges\n")  
  ranges <- Bacon.rng(d, set, BCAD=BCAD, prob=prob) 

  # calculate calendar axis limits

  modelranges <- range(ranges)
  dates <- set$calib$probs
  dateranges <- c()
  for(i in 1:length(dates))
    if(BCAD) 
      dateranges <- range(dateranges, 1950-dates[[i]][,1]) else
        dateranges <- range(dateranges, dates[[i]][,1])
  if(length(yr.min) == 0) 
    yr.min <- min(modelranges, dateranges)
  if(length(yr.max) == 0)
    yr.max <- max(modelranges, dateranges)
  yr.lim <- c(yr.min, yr.max)
  if(BCAD)
    yr.lim <- rev(yr.lim)
  if(rev.yr)
    yr.lim <- rev(yr.lim)
  d.lim <- c(d.max, d.min)
  if(rev.d)
    d.lim <- d.lim[2:1]
    
  if(length(d.lab) == 0)
    d.lab <- paste("Depth (", unit, ")", sep="")
  if(length(yr.lab) == 0)
    yr.lab <- ifelse(BCAD, "BC/AD", "cal yr BP")

  par(xaxs=xaxs, yaxs=yaxs, bty="n")
  if(rotate.axes)
    plot(0, type="n", ylim=d.lim, xlim=yr.lim, ylab=d.lab, xlab=yr.lab, bty="n", xaxt=xaxt, yaxt=yaxt) else
      plot(0, type="n", xlim=d.lim[2:1], ylim=yr.lim, xlab=d.lab, ylab=yr.lab, bty="n", xaxt=xaxt, yaxt=yaxt)

  if(!dates.only) {
    if(talk)
      cat("\nPreparing ghost graph\n")    
    .agedepth.ghost(set, rotate.axes=rotate.axes, BCAD=BCAD, d.res=d.res, yr.res=yr.res, grey.res=grey.res, dark=dark, colours=greyscale, d.min=d.min, d.max=d.max, yr.lim=yr.lim)
  }
  
  if(length(set$slump) > 0 ) 
    for(i in 1:nrow(set$slump))
      if(rotate.axes)
        rect(min(yr.lim)-1e3, set$slump[i,1], max(yr.lim)+1e3, set$slump[i,2], col=slump.col, border=slump.col) else
          rect(set$slump[i,1], min(yr.lim)-1e3, set$slump[i,2], max(yr.lim)+1e3, col=slump.col, border=slump.col)      

  calib.plot(set, BCAD=BCAD, rotate.axes=rotate.axes, height=height, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, dates.col=dates.col, new.plot=FALSE, normalise.dists=normalise.dists)
  legend(title.location, title, bty="n", cex=1.5)
  box(bty=bty)

  hiatus.depths <- set$hiatus.depths
  if(!is.na(set$boundary[1]))
    hiatus.depths <- set$boundary
  if(length(hiatus.depths) > 0) # draw locations of boundaries/hiatuses
    if(rotate.axes)
      abline(h=hiatus.depths, col=hiatus.col, lty=hiatus.lty) else
        abline(v=hiatus.depths, col=hiatus.col, lty=hiatus.lty)

  th <- rbind(1, nrow(ranges))
  if(!is.na(set$hiatus.depths[1])) {
    hi.d <- c()
    for(i in set$hiatus.depths)
      hi.d <- c(hi.d, max(which(d <= i)))
    th <- array(sort(c(1, nrow(ranges), hi.d-1, hi.d)), dim=c(2,length(hi.d)+1))
  }

  if(!dates.only)
    for(i in 1:ncol(th)) {
      h <- th[1,i] : th[2,i]
      if(rotate.axes) {
        lines(ranges[h,1], d[h], col=range.col, lty=range.lty)
        lines(ranges[h,2], d[h], col=range.col, lty=range.lty)
        lines(ranges[h,3], d[h], col=med.col, lty=med.lty) # median
        lines(ranges[h,4], d[h], col=mn.col, lty=mn.lty) # mean
      } else {
          lines(d[h], ranges[h,1], col=range.col, lty=range.lty)
          lines(d[h], ranges[h,2], col=range.col, lty=range.lty)
          lines(d[h], ranges[h,3], col=med.col, lty=med.lty) # median
          lines(d[h], ranges[h,4], col=mn.col, lty=mn.lty) # mean
        }
    }
  
  set$ranges <- cbind(d, round(ranges, rounded))
  colnames(set$ranges) <- c("depth", "min", "max", "median", "mean")
  .assign_to_global("info", set)

  if(plot.pdf)
    if(names(dev.cur()) != "null device")
      dev.copy2pdf(file=paste(set$prefix, ".pdf", sep=""))

  write.table(set$ranges, paste(set$prefix, "_ages.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
  rng <- abs(round(set$ranges[,3]-set$ranges[,2], rounded))
  min.rng <- d[which(rng==min(rng))]
  max.rng <- d[which(rng==max(rng))]
  if(length(min.rng)==1)
    min.rng <- paste(" yr at", min.rng, noquote(set$unit)) else
      min.rng <- paste(" yr between", min(min.rng), "and", max(min.rng), noquote(set$unit))
  if(length(max.rng)==1)
    max.rng <- paste(" yr at", max.rng, set$unit) else
      max.rng <- paste(" yr between", min(max.rng), "and", max(max.rng), noquote(set$unit))
  if(talk) {
    if(!dates.only)
      cat("\nMean ", 100*prob, "% confidence ranges ", round(mean(rng), rounded), " yr, min. ",
        min(rng), min.rng, ", max. ", max(rng), max.rng, "\n", sep="")
    overlap()
  } else
      cat("\n")
}
