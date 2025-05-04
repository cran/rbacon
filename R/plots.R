
#' @name proxy.ghost
#' @title Proxies analysed along the depths of a core can be plotted as 'proxy-ghost' graphs against calendar time while taking into account chronological uncertainties. Here darker grey indicates more likely calendar ages for specific proxy values.
#' @description Proxies analysed along the depths of a core can be plotted as 'proxy-ghost' graphs against calendar time while taking into account chronological uncertainties. Here darker grey indicates more likely calendar ages for specific proxy value.
#' @details Place a csv file with the values of proxies against depth within your core's folder. The values should be in columns separated by commas (default \code{sep=","}), the first column containing the depths and the first line (header) containing the proxy names.
#' The file name should start with the core's name and end with "_proxies.csv". For an example see \code{"Bacon_coredir/MSB2K/MSB2K_proxies.csv"} or \code{"Cores/MSB2K/MSB2K_proxies.csv"}.
#' @param proxy Which proxy to use (counting from the column number in the .csv file after the depths column).
#' @param proxy.lab Label of the proxy axis. Default names are taken from the csv file.
#' @param proxy.res Greyscale pixels are calculated for \code{proxy.res=250} proxy values by default, as a compromise between image quality and calculation speed. If the output looks very pixel-like (e.g., when choosing to plot only part of the record using proxy.lim), set this option to higher values.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{age.res=250} as a compromise between image quality and calculation speed. If the output looks very pixel-like (e.g., when choosing to plot only part of the record using age.lim), set this option to higher values.
#' @param yr.res Deprecated - use age.res instead
#' @param rgb.scale The function to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model. Default \code{rgb.res=100}.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param cutoff Point below which colours will no longer be printed. Default \code{cutoff=0.001}.
#' @param dark By default, the darkest grey value is assigned to the most likely value within the entire core (normalised to 1; \code{dark=1}). By setting dark to, e.g., \code{dark=.8}, all values of and above 0.8 will be darkest (and values below that threshold will be lighter grey the lower their probabilities).
#' @param darkest Darkness of the most likely value. Is black by default (\code{darkest=1}); lower values will result in lighter grey.
#' @param rotate.axes The default is to plot the calendar horizontally, however the plot can be rotated (\code{rotate.axes=TRUE}).
#' @param rev.proxy The proxy axis can be reversed if \code{rev.proxy=TRUE}.
#' @param rev.age The calendar axis can be reversed using \code{rev.age=TRUE}.
#' @param yr.rev Deprecated - use rev.age instead
#' @param plot.mean The mean ages of the proxy values can be added using \code{plot.mean=TRUE}.
#' @param mean.col Colour of the weighted mean ages of the proxy values.
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=NULL}).
#' @param yr.lim Deprecated - use age.lim instead
#' @param proxy.lim Ranges of the proxy axis, calculated automatically by default (\code{proxy.lim=NULL}).
#' @param sep Separator between the fields of the plain text file containing the depth and proxy data.
#' @param xaxs Extension of x-axis. By default, no white-space will be added at the axis extremes (\code{xaxs="i"}). See ?par for other options.
#' @param yaxs Extension of y-axis. By default, no white-space will be added at the axis extremes (\code{xaxs="i"}). See ?par for other options.
#' @param xaxt The x-axis is plotted by default, but this can be switched off using \code{xaxt="n"}.
#' @param yaxt The y-axis is plotted by default, but this can be switched off using \code{yaxt="n"}.
#' @param bty Type of box to be drawn around the plot (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.lab Deprecated - use age.lab instead
#' @param verbose Provide feedback on what is happening (default \code{verbose=TRUE}).
#' @param add Add to an existing graph (default \code{add=FALSE}).
#' @author Maarten Blaauw, J. Andres Christen
#' @return A grey-scale graph of the proxy against calendar age.
#' @examples
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   layout(1)
#'   proxy.ghost()
#' }
#' @export
proxy.ghost <- function(proxy=1, proxy.lab=NULL, proxy.res=250, age.res=200, yr.res=age.res, rgb.scale=c(0,0,0), rgb.res=100, set=get('info'), cutoff=0.001, dark=1, darkest=1, rotate.axes=FALSE, rev.proxy=FALSE, rev.age=FALSE, yr.rev=rev.age, plot.mean=FALSE, mean.col="red", age.lim=NULL, yr.lim=age.lim, proxy.lim=NULL, sep=",", xaxs="i", yaxs="i", xaxt="s", yaxt="s", bty="l", BCAD=set$BCAD, age.lab=ifelse(BCAD, "BC/AD", "cal yr BP"), yr.lab=age.lab, verbose=TRUE, add=FALSE) {
  if(length(set$Tr)==0)
    stop("please first run agedepth()", call.=FALSE)
  proxies <- read.csv(paste0(set$coredir, set$core, "/", set$core, "_proxies.csv"), header=TRUE, sep=sep)
  if(length(proxy.lab)==0)
    proxy.lab <- names(proxies)[proxy+1]
  proxy <- cbind(as.numeric(proxies[,1]), as.numeric(proxies[,proxy+1]))
  proxy <- proxy[!is.na(proxy[,2]),]
  proxy <- proxy[which(proxy[,1] <= set$d.max),]
  proxy <- proxy[which(proxy[,1] >= set$d.min),]
  pr.mn.ages <- approx(set$ranges[,1], set$ranges[,5], proxy[,1], rule=1)$y
  if(length(unique(proxy[,2])) == 1)
    stop("this proxy's values remain constant throughout the core, and cannot be proxy-ghosted!", call.=FALSE)
  proxyseq <- seq(min(proxy[,2]), max(proxy[,2]), length=proxy.res)
#  out <- list(yrseq=c(), binned=c(), maxs=c())
  ds <- NULL
  d.length <- array(1, dim=c(proxy.res, 2))

  for(i in 1:proxy.res) {
    tmp  <- .DepthsOfScore(proxyseq[i], proxy)
    ds <- c(ds, tmp)
    if(i > 1)
      d.length[i,1] <- d.length[(i-1),2]+1
    d.length[i,2] <- d.length[i,1]+length(tmp)-1
    if(length(tmp) == 0)
      d.length[i,] <- d.length[i-1,]
  }
  if(verbose)
    message("Calculating histograms")

  hists <- Bacon.hist(ds, set, calc.range=FALSE) # BCAD always FALSE
  message("\n")

  age.min <- c()
  age.max <- c()
  for(i in 1:length(hists)) {
    age.min <- min(age.min, hists[[i]]$th0)
    age.max <- max(age.max, hists[[i]]$th1)
  }
  age.seq <- seq(age.min, age.max, length=age.res)

  all.counts <- array(0, dim=c(length(hists), length(age.seq)))
  for(i in 1:length(hists))
    all.counts[i,] <- approx(seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n), hists[[i]]$counts, age.seq)$y
  all.counts[is.na(all.counts)] <- 0
  all.counts <- all.counts/max(all.counts)
  all.counts[all.counts > dark] <- dark
  max.counts <- array(0, dim=c(proxy.res, length(age.seq)))
  for(i in 1:proxy.res)
    for(j in 1:length(age.seq))
      max.counts[i,j] <- max(all.counts[d.length[i,1]:d.length[i,2],j])
  if(dark>1)
    stop("dark values larger than 1 are not allowed\n", call.=FALSE) else
      max.counts[max.counts > dark] <- dark
  max.counts[max.counts < cutoff] <- NA # don't plot too small/light values
  if(length(age.lim) == 0)
    if(xaxs=="r")
      age.lim <- extendrange(pretty(age.seq), f=.04) else
        age.lim <- range(age.seq)[2:1]
  max.counts <- max.counts[,ncol(max.counts):1]	# tmp May 2025
  if(rev.proxy)
	  max.counts <- max.counts[nrow(max.counts):1,]
  if(rev.age) {
    age.lim <- age.lim[2:1]
	max.counts <- max.counts[,ncol(max.counts):1]
  }
  if(BCAD) {
    age.lim <- calBPtoBCAD(age.lim)
    max.counts <- max.counts[,ncol(max.counts):1]
    age.seq <- calBPtoBCAD(age.seq)
  }

  if(length(proxy.lim) == 0)
    proxy.lim <- range(proxyseq)
  if(rev.proxy)
    proxy.lim <- rev(proxy.lim)
  col <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0, darkest, length=rgb.res))
  if(rotate.axes) {
    if(!add)
	  plot(0, type="n", xlim=proxy.lim, ylim=age.lim, ylab=age.lab, xlab=proxy.lab, xaxs=xaxs, yaxs=yaxs, xaxt=xaxt, yaxt=yaxt) 
	image(proxyseq, age.seq, max.counts, col=col, add=TRUE, useRaster=TRUE)
    if(plot.mean)
      lines(proxy[,2], pr.mn.ages, col=mean.col)
  } else {
      if(!add)
        plot(0, type="n", ylim=proxy.lim, xlim=age.lim, xlab=age.lab, ylab=proxy.lab, xaxs=xaxs, yaxs=yaxs, xaxt=xaxt, yaxt=yaxt)
	image(age.seq, proxyseq, t(max.counts), col=col, add=TRUE, useRaster=TRUE)  
    if(plot.mean)
      lines(pr.mn.ages, proxy[,2], col=mean.col)
  }
  box(bty=bty)
  
  invisible(list(ages=age.seq, proxy=proxyseq, counts=max.counts, means=cbind(pr.mn.ages, proxy[,2])))
}



# for the proxy.ghost function
.DepthsOfScore <- function(value, dat) {
  d <- c()
  for(i in 1:(nrow(dat)-1)) {
    valueRange <- dat[i:(i+1),2]
    if(min(valueRange) <= value && max(valueRange) >= value) {
      slope <- (dat[i,2] - dat[i+1,2]) / (dat[i,1] - dat[i+1,1])
      intercept <- dat[i,2] - (slope*dat[i,1])
      if(slope == 0 && i > 1)
        d[i-1] <- dat[i,1]
      d <- sort(c(d, (value - intercept) / slope ))
    }
  }
  unique(d)
}



#' @name AgesOfEvents
#' @title Event probabilities against calendar age
#' @description Plot probability curves for events in the core, expressed against calendar age.
#' @details Probabilities of depths with 'events' in an age-modelled core can be plotted against time, taking into account
#' chronological uncertainties (Blaauw et al. 2007). Such events could be for example core depths at which proxies
#' indicate changes toward wetter local conditions. This can be expressed as values between 0 (no event) and 1 (event at 100\% probability)
#' for each depth.
#'
#' Blaauw et al. 2010 propose to estimate probabilities of events by finding specific proxy features such as increasing curves.
#' Probabilities are then estimated through resampling from the proxy values, where low to modest rises of proxy curves result
#' in low event probabilities, and clear proxy rises in high probabilities. A smooth spline can be applied to adapt the balance of short-term vs long-term events. To calculate the event probabilities,
#' produce a file with two columns (depth and corresponding proxy-derived probabilities, separated by white spaces).
#'  Do not provide headers at the file's first line, and save the file with extension "_events.txt" within the core's Bacon folder. See Cores/MSB2K/MSB2K_events.txt (or Bacon_runs/MSB2K/MSB2K_events.txt) for an example. Events are calculated as the probability that an event took place within specific time windows - or more specifically, that the Bacon age-depth model puts depths with assigned event probabilities in that time window.
#'
#' does not yet deal correctly with hiatuses.
#' @param window Width of the window.
#' @param move Step size with which the window moves.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param plot.steps Plot probability values step-wise (defaults to \code{plot.steps=FALSE}, which plots smooth curves instead).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.lab Deprecated - use age.lab instead
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param yr.lim Deprecated - use age.lim instead
#' @param prob.lab Label of the probability axis (default \code{prob.lab="probability"}).
#' @param prob.lim Limits of the probability axis (calculated automatically by default).
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.age The direction of the age axis, which can be reversed using \code{rev.age=TRUE}.
#' @param rev.yr Deprecated - use rev.age instead
#' @param yaxs Extension of the y-axis. Defaults to the exact ranges of the probability values. White space can be added to the vertical axis using \code{yaxs="r"}.
#' @param bty Type of box to be drawn around plots. Draw a box around the graph (\code{"n"} for none, and \code{"l"}, \code{"7"},
#'  \code{"c"}, \code{"u"}, "]" or \code{"o"} for correspondingly shaped boxes).
#' @author Maarten Blaauw, J. Andres Christen
#' @return The resulting probabilities are plotted and saved within the core's folder (file names ending with the window width and "_probs.txt").
#' @examples
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   AgesOfEvents(100, 10)
#' }
#' @references
#' Blaauw, M., Christen, J.A., Mauquoy, D., van der Plicht, J., Bennett, K.D. (2007) Testing the timing of radiocarbon-dated events between proxy archives. _The Holocene_, *17*, 283-288.
#' Blaauw, M., Wohlfarth, B., Christen, J.A., Ampel, L., Veres, D., Hughen, K.A., Preusser, F., Svensson, A. (2010) Were last glacial climate events simultaneous between Greenland and France? A quantitative comparison using non-tuned chronologies. _Journal of Quaternary Science_ *25*, 387-394.
#' @export
AgesOfEvents <- function(window, move, set=get('info'), plot.steps=FALSE, BCAD=set$BCAD, age.lab=c(), yr.lab=age.lab, age.lim=c(), yr.lim=age.lim, prob.lab="probability", prob.lim=c(), rotate.axes=FALSE, rev.age=TRUE, rev.yr=rev.age, yaxs="i", bty="l") {
  if(move == 0)
    stop("I cannot move anywhere if move = 0", call.=FALSE)
  outfile <- paste0(set$prefix, "_", window, "_probs.txt")
  file.create(outfile)
  MCMCname <- paste0(set$prefix, ".out")
  probfile <- paste0(set$coredir, set$core, "/", set$core, "_events.txt")
  if(!file.exists(probfile))
    stop("file with probabilities for events per depth not found! Check the manual", call.=FALSE)
  probs <- fastread(probfile)
  if(!is.numeric(probs[1,1]))
    stop("first line of the _events.txt file should NOT contain titles; please remove them", call.=FALSE)
  if(min(probs[,1]) < min(set$elbows) || max(probs[,1]) > max(set$elbows)) {
    message("some depths in the _events.txt file go beyond the age-model; I will remove them")
    file.rename(probfile, paste0(probfile, "_backup"))
    probs <- probs[which(probs[,1] >= min(set$elbows)),]
    probs <- probs[which(probs[,1] <= max(set$elbows)),]
    fastwrite(probs, probfile, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }

  if(length(age.lim) == 0) {
    min.age <- min(set$ranges[,2])
    max.age <- max(set$ranges[,3])
    age.lim <- c(min.age, max.age)
  } else {
      min.age <- min(age.lim)
      max.age <- max(age.lim)
    }

  events(min.age, max.age, move, window, outfile, MCMCname, nrow(set$output), set$K, set$elbows[1], set$thick, probfile, nrow(probs))
  probs <- fastread(outfile)
  if(BCAD) {
    probs[,1] <- calBPtoBCAD(probs[,1])
    o <- order(probs[,1])
    probs <- probs[o,]
  }
  #close(outfile) # May 2021

  if(plot.steps) {
    d.sort <- sort(rep(1:nrow(probs),2))
    d.sort <- cbind(d.sort[-length(d.sort)], d.sort[-1])
    probs <- cbind(c(min(probs[,1]), probs[d.sort[,1],1], max(probs[,1])), c(0,probs[d.sort[,2],2],0))
  } else
      probs <- cbind(c(min(probs[,1]), probs[,1], max(probs[,1])), c(0,probs[,2],0))
  oldpar <- par(yaxs=yaxs, bty=bty)
  on.exit(par(oldpar))

  if(rev.age)
    age.lim <- rev(age.lim)
  if(length(age.lab) == 0)
    age.lab <- ifelse(BCAD, "BC/AD", "cal BP")
  if(length(prob.lim) == 0)
    prob.lim <- c(0, 1.1*max(probs[,2]))

  if(rotate.axes) {
    plot(probs[,2], probs[,1], type="n", ylab=age.lab, xlab=prob.lab, xlim=prob.lim, ylim=age.lim)
    polygon(probs[,2:1], col="grey")
  } else {
      plot(probs, type="n", xlab=age.lab, ylab=prob.lab, ylim=prob.lim, xlim=age.lim)
      polygon(probs, col="grey")
    }
  if(move > window) 
    message("\nAre you sure you want the window widths to be smaller than the moves?")
  invisible(probs)
}
