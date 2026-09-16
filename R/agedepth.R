#' @title plot an age-depth model
#' @description Plot the age-depth model of a core.
#' @details After loading a previous run, or after running either the \link{scissors} or \link{thinner} command, plot the age-model
#' again using the command \code{agedepth()}.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}. See also \code{depth.unit}.
#' @param age.lab The labels for the calendar axis (default \code{age.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param kcal Use kcal BP. Default is \code{kcal=FALSE}.
#' @param yr.lab Deprecated - use age.lab instead
#' @param acc.lab The labels for the accumulation rate plot (top middle). Default \code{d.lab="Acc. rate (yr/cm)"} (or whatever units you're using).
#' @param mem.lab The labels for the memory plot (top right). Default \code{d.lab="Memory"}.
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth). 
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param d.by Depth intervals at which ages are calculated. Default 1. Alternative depth intervals can be provided using, e.g., d.\code{d.by=0.5}.
#' @param depths By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}. Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}.
#' @param depths.file By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.
#' If \code{depths.file=TRUE}, Bacon will read a file containing the depths for which you require ages.
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir},
#' and its name should start with the core's name and end with '_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param accordion An experimental option to squeeze and stretch depths below a boundary. Needs 2 parameters: boundary depth and compression ratio, e.g., \code{accordion=c(25,20)}. Defaults to inactive, \code{accordion=c()}.
#' @param plotatthesedepths An option to plot ages at different depths (e.g., if a core has been compressed during a run). Use with extreme caution!
#' @param age.min Minimum age of the age-depth plot.
#' @param yr.min Deprecated - use age.min instead.
#' @param age.max Maximum age of the age-depth plot.
#' @param yr.max Deprecated - use age.min instead.
#' @param hiatus.option How to calculate accumulation rates and ages for sections with hiatuses. Either extrapolate from surrounding sections (default, \code{hiatus.option=1}), use a w-weighted mix between the prior and posterior values for depths below the hiatus and prior information only for above the hiatus (\code{hiatus.option=2}), or use the originally calculated slopes (\code{hiatus.option=0}).
#' @param dark The value beyond which any higher values will be clipped. Set to 1 by default. This can be set to lower values if for example only a very small area in the plot obtains the darkest values - then lowering the `dark` value will darker a larger area.
#' @param darkest The darkest grey value is darkest=1 by default; lower values will result in lighter maximum colours; values >1 are not advised.
#' @param prob Confidence interval to report (between 0 and 1, default 0.95 or 95\%).
#' @param rounded Rounding of years. Default is to round to single years (1 digit for plum models).
#' @param d.res Resolution or amount of greyscale pixels to cover the depth scale of the age-model plot. Default \code{d.res=200}.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param yr.res Deprecated - use age.res instead.
#' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
#' @param rotate.axes By default, the age-depth model is plotted with the depths on the horizontal axis and ages on the vertical axis. This can be changed with \code{rotate.axes=TRUE}.
#' @param rev.age The direction of the age axis, which can be reversed using \code{rev.age=TRUE}.
#' @param rev.yr Deprecated - use rev.age instead.
#' @param rev.d The direction of the depth axis, which can be reversed using \code{rev.d=TRUE}.
#' @param use.raster Rasters can be aligned or not in the underlying image function. Setting \code{use.raster=FALSE, default} takes a bit longer to draw and sometimes causes strange lines owing to anti-aliasing. However, the alternative of \code{use.raster=TRUE} causes greyscales on some devices (e.g., OSX quartz) to 'flip'. If this is the case, use 'flip.acc=TRUE'.
#' @param flip.age When using \code{use.raster=TRUE}, sometimes greyscales are flipped. If this is the case, see if setting \code{flip.age=TRUE} can flip the ages back again. 
#' @param flip.d When using \code{use.raster=TRUE}, sometimes greyscales are flipped. If this is the case, see if setting \code{flip.d=TRUE} can flip the depths back again. 
#' @param depth.unit Units of the depths. Defaults to the one provided in the Bacon() command, \code{depth.unit=set$depth.unit}.
#' @param age.unit Units of the ages. Defaults to \code{age.unit="yr"}.
#' @param unit Deprecated and replaced by \code{depth.unit}.
#' @param maxcalc Number of depths to calculate ages for. If this is more than \code{maxcalc=500}, a warning will be shown that calculations will take time.
#' @param height The maximum heights of the distributions of the dates on the plot. See also \code{normalise.dists}.
#' @param calheight Multiplier for the heights of the distributions of dates on the calendar scale. Defaults to \code{calheight=1}.
#' @param ex As an alternative to 'height' and 'calheight', the distribution heights can be set as either a single value (e.g., \code{ex=1}) or for each date (for example, \code{myex=rep(1, nrow(info$dets)); myex[1:5] <- 10; agedepth(ex=myex)}).
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.1}).
#' @param plot.range Whether or not to plot the curves showing the confidence ranges of the age-model. Defaults to (\code{plot.range=TRUE}).
#' @param range.col The colour of the curves showing the confidence ranges of the age-model. Defaults to medium grey (\code{range.col=grey(0.5)}).
#' @param range.lty The line type of the curves showing the confidence ranges of the age-model. Defaults to \code{range.lty=12}.
#' @param range.lwd Widths of the lines of the ranges of the age-depth model. Default \code{range.lwd=1}.
#' @param mn.col The colour of the mean age-depth model: default \code{mn.col="red"}.
#' @param mn.lty The line type of the mean age-depth model. Default \code{mn.lty=12}.
#' @param mn.lwd Width of the line of the mean age-depth model. Default \code{mn.lwd=1}.
#' @param med.col The colour of the median age-depth model: not drawn by default \code{med.col=NA}.
#' @param med.lty The line type of the median age-depth model. Default \code{med.lty=12}.
#' @param med.lwd Width of the line of the median age-depth model. Default \code{med.lwd=1}.
#' @param C14.col The colour of the calibrated ranges of the dates. Default is semi-transparent blue: \code{C14.col=rgb(0,0,1,.35)}.
#' @param C14.border The colours of the borders of calibrated 14C dates. Default is semi-transparent dark blue: \code{C14.border=rgb(0, 0, 1, 0.5)}.
#' @param cal.col The colour of the non-14C dates. Default is semi-transparent blue-green: \code{cal.col=rgb(0,.5,.5,.35)}.
#' @param cal.border The colour of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{cal.border=rgb(0,.5,.5,.5)}. Not used by default.
#' @param dates.col As an alternative to colouring dates based on whether they are 14C or not, sets of dates can be coloured as, e.g., \code{dates.col=colours()[2:100]}.
#' @param pb.background Probability at which total Pb values are considered to have reached background values, or in other words, that their modelled values are at or below supported + detection limit (Al)). Setting this at 0.5 means that any depth with a Pb measurement, where at least half of the iterations model Pb values reaching background values, is flagged as having reached background. The age-model is not extended to any Pb measurements that have reached background.
#' @param pbmodelled.col Colour of the modelled 210Pb values. Defaults to shades of blue: \code{pbmodelled.col=function(x) rgb(0,0,1,x)}.
#' @param pbmeasured.col Colour of the measured 210Pb values (default \code{pbmeasured.col="blue"}). Draws rectangles of the upper and lower depths as well as the Pb values with 95 percent error ranges. 
#' @param pb.lim Axis limits for the Pb-210 data. Calculated automatically by default (\code{pblim=c()}).
#' @param supp.col Colour of supported Pb-210. Defaults to semi-transparent purple, because why not.
#' @param remove.tail Whether or not to remove the tail measurements when plotting. Sometimes automated removal might go wrong, so then this option can be used to avoid removing the tail measurements.
#' @param MCMC.resample After the MCMC run, if there are more MCMC iterations than requested, only the last 'ssize' iterations will be retained. Defaults to TRUE.
#' @param hiatus.col The colour of the depths of any hiatuses. Default \code{hiatus.col=grey(0.5)}.
#' @param hiatus.lty The line type of the depths of any hiatuses. Default \code{hiatus.lty=12}.
#' @param from.col Starting colour of the gradient for the ghost plot of the age-depth model. Defaults to \code{from.col="white"}.
#' @param to.col Colour of the gradient for the ghost plot of the age-depth model. Defaults to a greyscale, \code{to.col="black"}.
#' @param rgb.scale As alternative to 'from.col' and 'to.col', rgb.scale can be used to produce a coloured representation of all age-models. Needs 3 values for the intensity of red, green and blue. Defaults to grey-scales: \code{rgb.scale=c(0,0,0)}, but could also be, say, scales of red (\code{rgb.scale=c(1,0,0)}). 
#' @param rgb.res Resolution of the colour spectrum depicting the age-depth model, if using 'rgb.scale'. Default \code{rgb.res=100}.
#' @param slump.col Colour of slumps. Defaults to \code{slump.col=grey(0.8)}.
#' @param normalise.dists By default, the distributions of more precise dates will cover less time and will thus peak higher than less precise dates. This can be avoided by specifying \code{normalise.dists=FALSE}.
#' @param same.heights Plot the distributions of the dates all at the same maximum height (default \code{same.height=FALSE}).
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine), \code{cc=3} for SHCal20 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param title The title of the age-depth model is plotted on the main panel. By default this is the core's name. To leave empty: \code{title=""}.
#' @param title.location Location of the title. Default \code{title.location='topleft'}.
#' @param title.size Size of the title font. Defaults to \code{title.size=1.5}.
#' @param plot.labels Whether or not to plot labels next to the dated depths. Defaults to \code{FALSE}.
#' @param labels Add labels to the dates (as given by the first column of the .csv file). \code{FALSE} by default.
#' @param label.age Position on the age axis of the date labels. By default draws them before the youngest age (1), but can also draw them after the oldest age (2), or above its mean (3).
#' @param label.offset Offsets of the positions of the labels, giving the x and y offsets. Defaults to c(0,0).
#' @param label.size Size of labels.
#' @param label.col Colour of the labels. Defaults to the colour given to the borders of the dates.
#' @param label.adj  Justification of the labels. Follows R's adj option: A value of '0' produces left-justified text, '0.5' (the default) centered text and '1' right-justified text.
#' @param label.rot Rotation of the label. 0 by default (horizontal).
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param bty Type of box to be drawn around plots (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param mar.left Plot margins for the topleft panel (amount of white space along edges of axes 1-4). Default \code{mar.left=c(3,3,1,1)}.
#' @param mar.middle Plot margins for the middle panel(s) at the top (amount of white space along edges of axes 1-4). Default \code{mar.middle=c(3,3,1,1)}.
#' @param mar.right Plot margins for the topright panel (amount of white space along edges of axes 1-4). Default \code{mar.right=c(3,3,1,1)}.
#' @param mar.main Plot margins for the main panel (amount of white space along edges of axes 1-4). Default \code{mar.main=c(3,3,1,1)}.
#' @param righthand Adapt the righthand margins by a certain amount (default 2) to allow a righthand axis to be plotted (for plum)
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, 0.7, 0.0)}.
#' @param xaxs Extension of x-axis. By default, add some extra white-space at both extremes (\code{xaxs="r"}). See ?par for other options.
#' @param yaxs Extension of y-axis. By default, add no extra white-space at both extremes (\code{yaxs="i"}). See ?par for other options.
#' @param MCMC.col Colour of the MCMC output. Defaults to \code{post.col=grey(0.4))}.
#' @param post.col Colour of the posterior histogram. Defaults to \code{post.col=grey(0.8))}. In case of hiatuses/boundaries, you can provide colours for the individual sections, e.g., post.col=c(2,3,4).
#' @param post.border Colour of the posterior border. Defaults to \code{post.border=grey(0.4))}.  In case of hiatuses/boundaries, you can provide colours for the individual sections, e.g., post.col=c(2,3,4).
#' @param acc.post.col Colour of the posterior of the accumulation rate (yr/cm). Defaults to \code{post.col=grey(0.8))}. Multiple colours can be provided. See post.col.
#' @param acc.post.border Colour of the posterior border of the accumulation rate (yr/cm). Defaults to \code{post.border=grey(0.4))}. Multiple colours can be provided. See post.border.
#' @param mem.post.col Colour of the posterior of the memory. Defaults to \code{post.col=grey(0.8))}. See post.col. 
#' @param mem.post.border Colour of the posterior of the memory. Defaults to \code{post.border=grey(0.4))}. See post.border.
#' @param hiatus.post.col Colour of the posterior of the hiatus length. Defaults to \code{post.col=grey(0.8))}. Multiple colours can be provided. See post.col.
#' @param hiatus.post.border Colour of the border of the hiatus length. Defaults to \code{post.border=grey(0.4))}. Multiple colours can be provided. See post.border.
#' @param phi.post.col Colour of the posterior of the 210Pb influx phi. Defaults to \code{post.col=grey(0.8))}. See post.col.
#' @param phi.post.border Colour of the posterior border of the 210Pb influx phi. Defaults to \code{post.border=grey(0.4))}. See post.border.
#' @param s.post.col Colour of the posterior of the supported 210Pb. Defaults to \code{post.col=grey(0.8))}. See post.col.
#' @param s.post.border Colour of the border of the posterior supported 210Pb. Defaults to \code{post.border=grey(0.4))}. See post.border.
#' @param prior.col Colour of the prior curves. Defaults to light green, \code{prior.col=3)}.
#' @param prior.lwd Line width of the prior curves. Defaults to \code{prior.lwd=2)}.
#' @param prior.fontcol Colour of the font accompanying the posterior histograms. Defaults to red, \code{prior.fontcol=2)}.
#' @param prior.col Colour of the prior curves. Defaults to light green, \code{prior.col=3)}.
#' @param prior.lwd Line width of the prior curves. Defaults to \code{prior.lwd=2)}.
#' @param prior.fontcol Colour of the font accompanying the posterior histograms. Defaults to red, \code{prior.fontcol=2)}.
#' @param prior.ticks Plot tickmarks and values on the vertical axes for the prior and posterior distributions. Defaults to no tick marks (\code{prior.ticks="n"}). Set to \code{prior.ticks="s"} to plot the tick marks. Note that these values are of little practical use, as they correspond poorly to, e.g., the mean and strength values. All that matters is that the areas of both the prior and the posterior distributions sum to 1; wider distributions tend to give lower peaks, and narrower distributions higher peaks. 
#' @param prior.fontsize Font size of the prior, relative to R's standard size. Defaults to \code{prior.fontsize=0.9}.
#' @param toppanel.fontsize Font size of the top panels, relative to R's standard size. Defaults to \code{prior.fontsize=0.9}.
#' @param mainpanel.tickfontsize Font size of values at the tick marks in the main panel, relative to R's standard size. Defaults to \code{mainpanel.tickfontsize=1}.
#' @param mainpanel.labelfontsize Font size of axis labels in the main panel, relative to R's standard size. Defaults to \code{mainpanel.labelsize=1}.
#' @param acc.xlim Horizontal axis limits of the accumulation rate panel. Calculated automatically by default.
#' @param acc.ylim Vertical axis limits of the accumulation rate panel. Calculated automatically by default.
#' @param mem.xlim Horizontal axis limits of the memory panel. Calculated automatically by default.
#' @param mem.ylim Vertical axis limits of the memory panel. Calculated automatically by default.
#' @param hiatus.xlim Horizontal axis limits of the hiatus size panel. Calculated automatically by default.
#' @param hiatus.ylim Vertical axis limits of the hiatus size panel. Calculated automatically by default.
#' @param phi.xlim Horizontal axis limits of the phi panel. Calculated automatically by default.
#' @param phi.ylim Vertical axis limits of the phi panel. Calculated automatically by default.
#' @param supp.xlim Horizontal axis limits of the supported-Pb panel. Calculated automatically by default.
#' @param supp.ylim Vertical axis limits of the supported-Pb panel. Calculated automatically by default.
#' @param xaxt Whether or not to plot the x-axis. Can be used to adapt axes after a plot. See ?par for other options.
#' @param yaxt Whether or not to plot the y-axis. Can be used to adapt axes after a plot. See ?par for other options.
#' @param plot.pb Plot the 210Pb data (if present). Defaults to \code{plot.pb=TRUE}.
#' @param pb.lty Line type of measured Pb-210 data.
#' @param plot.pdf Produce a pdf file of the age-depth plot.
#' @param quartz Use quartz-based pdf plots if available (on Mac OS). Defaults to FALSE.
#' @param cairo Use cairo-based pdf plots if available (on Mac OS). Defaults to FALSE.
#' @param model.only By default, panels showing the MCMC iterations and the priors and posteriors for accumulation rate and memory are plotted above the main age-depth model panel. This can be avoided by supplying \code{model.only=TRUE}. Note however that this removes relevant information to evaluate the age-depth model, so we do recommend to present age-models together with these upper panels.
#' @param dates.only By default, the age-depth model is plotted on top of the dates. This can be avoided by supplying \code{dates.only=TRUE}.
#' @param verbose Provide a summary of the age ranges after producing the age-depth model graph; default \code{verbose=FALSE}.
#' @param roundby Rounding of the values reported at the end of the function. Defaults to 2 decimals.
#' @param save.info By default, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @param write.summary Whether or not to write a file to the core's folder, which summarizes the run (MCMC, model fit, priors, posteriors, learning).
#' @param ssize Number of MCMC iterations to use.
#' @param use.cpp Whether or not to use a cpp function to calculate the ages. Defaults to TRUE, but can be set to FALSE (much slower but less experimental)
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of the age-depth model, and estimated ages incl. confidence ranges for each depth.
#' @examples
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   agedepth()
#' }
#' @export
agedepth <- function(set=get('info'), BCAD=set$BCAD, depth.unit=set$depth.unit, age.unit="yr", unit=depth.unit, d.lab=c(), age.lab=c(), yr.lab=age.lab, kcal=FALSE, acc.lab=c(), mem.lab=c(), d.min=c(), d.max=c(), d.by=c(), depths=set$depths, depths.file=FALSE, accordion=c(), plotatthesedepths=c(), age.min=c(), yr.min=age.min, age.max=c(), yr.max=age.max, hiatus.option=1, dark=1, darkest=1, prob=set$prob, rounded=c(), d.res=500, age.res=500, yr.res=age.res, date.res=100, rotate.axes=FALSE, rev.age=FALSE, rev.yr=rev.age, rev.d=FALSE, use.raster=FALSE, flip.age=FALSE, flip.d=FALSE, maxcalc=500, height=1, calheight=1, ex=1, mirror=TRUE, up=TRUE, cutoff=.1, plot.range=TRUE, range.col=grey(.5), range.lty="12", range.lwd=1, mn.col="red", mn.lty="12", mn.lwd=1, med.col=NA, med.lty="12", med.lwd=1, C14.col=rgb(0,0,1,.35), C14.border=rgb(0,0,1,.5), cal.col=rgb(0,.5,.5,.35), cal.border=rgb(0,.5,.5,.5), dates.col=c(), pb.background=.5, pbmodelled.col=function(x) rgb(0,0,1,.5*x), pbmeasured.col="blue", pb.lim=c(), supp.col=rgb(.5,0,.5,.5), remove.tail=TRUE, MCMC.resample=TRUE, hiatus.col=grey(0.5), hiatus.lty="12", from.col="white", to.col="black", rgb.scale=c(0,0,0), rgb.res=100, slump.col=grey(0.8), normalise.dists=TRUE, same.heights=FALSE, cc=set$cc, title=set$core, title.location="topleft", title.size=1.5, plot.labels=FALSE, labels=c(), label.age=1, label.size=0.8, label.col="black", label.offset=c(0,0), label.adj=c(0.5,0), label.rot=0, after=set$after, bty="l", mar.left=c(3,3,1,.5), mar.middle=c(3,0,1,.5), mar.right=c(3,0,1,.5), mar.main=c(3,3,1,1), righthand=3, mgp=c(1.7,.7,.0), xaxs="r", yaxs="i", MCMC.col=grey(.4), post.col=rgb(0,0,0,.2), post.border=rgb(0,0,0,.4), acc.post.col=c(), acc.post.border=c(), mem.post.col=c(), mem.post.border=c(), hiatus.post.col=c(), hiatus.post.border=c(), phi.post.col=c(), phi.post.border=c(), s.post.col=c(), s.post.border=c(), prior.col=3, prior.lwd=2, prior.fontcol=2, prior.ticks="n", prior.fontsize=0.9, toppanel.fontsize=0.9, mainpanel.tickfontsize=1, mainpanel.labelfontsize=1, acc.xlim=c(), acc.ylim=c(), mem.xlim=c(), mem.ylim=c(), hiatus.xlim=c(), hiatus.ylim=c(), phi.xlim=c(), phi.ylim=c(), supp.xlim=c(), supp.ylim=c(), xaxt="s", yaxt="s", plot.pb=TRUE, pb.lty=1, plot.pdf=FALSE, quartz=FALSE, cairo=FALSE, dates.only=FALSE, model.only=FALSE, verbose=TRUE, roundby=2, save.info=set$save.info, write.summary=TRUE, ssize=4e3, use.cpp=TRUE) {
# Load the output, if it exists
  outp <- paste0(set$prefix, ".out")
  if(file.exists(outp))
    set <- Bacon.AnaOut(outp, set, MCMC.resample)

  set$MCMCdiagnostics <- MCMC.diagnostics(set, ssize, talk=verbose)

  # Plum-specific
  if(set$isplum) {
    outPlum <- paste0(set$prefix, "_plum.out")
    if(file.exists(outPlum))
      set <- Plum.AnaOut(outPlum, set, MCMC.resample)
    if(save.info)
      assign_to_global("info", set)

    # sometimes runs don't go well, with the age-model totally lost. This is indicated by a very peaked posterior for memory, very close to 1
    if(mn <- mean(set$output[,set$K+2]) > 0.95)
      if(mn > set$mem.mean)
        message("\nWarning, this run has a very high posterior memory (", round(mn, 2), ") and probably didn't go very well. Please run again\n")
    bg <- background(set) # calculate which Pb data have likely reached background
    set$background <- bg
  }

  #set$BCAD <- BCAD # tmp May 2021
  # Adapt ages of sections which contain boundaries or hiatuses
  if(!is.na(set$hiatus.depths[1]))
    set <- hiatus.slopes(set)

  if(save.info)
    assign_to_global("info", set)

  oldpar <- par(no.readonly = TRUE)
  #on.exit(par(oldpar))
  par(mar=mar.main, mgp=mgp)

  if(!model.only) {
    par(mar=mar.left, bty=bty, mgp=mgp, xaxs=xaxs, yaxs=yaxs)

    ifelse(set$isplum, pn <- c(1:5, rep(6,5)), pn <- c(1:3, rep(4,3)))
    if(!is.na(set$hiatus.depths[1]))
      if(is.na(set$boundary[1]))
        ifelse(set$isplum, pn <- c(1:6, rep(7,6)), pn <- c(1:4, rep(5,4)))

    layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    PlotLogPost(set, 0, set$Tr, xaxs=xaxs, yaxs=yaxs, panel.size=toppanel.fontsize, col=MCMC.col) # convergence information

    par(mar=mar.middle) # reduce white space
    #on.exit(par(oldpar))

    if(length(acc.post.col) == 0)
      acc.post.col <- post.col
    if(length(acc.post.border) == 0)
      acc.post.border <- post.border
    if(length(mem.post.col) == 0)
      mem.post.col <- post.col
    if(length(mem.post.border) == 0)
      mem.post.border <- post.border
    if(length(hiatus.post.col) == 0)
      hiatus.post.col <- post.col
    if(length(hiatus.post.border) == 0)
      hiatus.post.border <- post.border
    if(length(phi.post.col) == 0)
      phi.post.col <- post.col
    if(length(phi.post.border) == 0)
      phi.post.border <- post.border
    if(length(s.post.col) == 0)
      s.post.col <- post.col
    if(length(s.post.border) == 0)
      s.post.border <- post.border

    set$post.acc <- PlotAccPost(set, depth.unit=depth.unit, age.unit=age.unit, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, acc.xlim=acc.xlim, acc.ylim=acc.ylim, acc.lab=acc.lab, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=acc.post.col, hist.border=acc.post.border)
    set$post.mem <- PlotMemPost(set, set$core, set$K, "", set$mem.strength, set$mem.mean, ds=1, thick=set$thick, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, mem.xlim=mem.xlim, mem.ylim=mem.ylim, mem.lab=mem.lab, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=mem.post.col, hist.border=mem.post.border)

    if(!is.na(set$hiatus.depths[1]))
      if(is.na(set$boundary[1])) {
        gaps <- PlotHiatusPost(set, mn=set$hiatus.mean, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, hiatus.xlim=mem.xlim, hiatus.ylim=mem.ylim, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=hiatus.post.col, hist.border=hiatus.post.border)

        set$post.hiatus <- cbind(gaps$post.mn, gaps$post.shape)
        set$hiatus.sizes <- gaps$gaps
    }

    if(set$isplum) {
       set$post.phi <- PlotPhiPost(set, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, phi.xlim=phi.xlim, phi.ylim=phi.ylim, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=phi.post.col, hist.border=phi.post.border)
       if(set$nPs > 1)
         prior.ticks <- "s" # because with varying supported Pb, the y-axis is important
       par(mar=mar.right)

       set$post.supp <- PlotSuppPost(set, xaxs=xaxs, yaxs=yaxs, yaxt=prior.ticks, prior.size=prior.fontsize, panel.size=toppanel.fontsize, supp.xlim=supp.xlim, supp.ylim=supp.ylim, line.col=prior.col, line.width=prior.lwd, text.col=prior.fontcol, hist.col=s.post.col, hist.border=s.post.border, data.col=supp.col)

       mar.main[4] <- mar.main[4] + righthand # to enable space for righthand axis
    }
    par(mar=mar.main) # new May 2021
  }

  # calculate and plot the ranges and 'best' estimates for each required depth
  depth.limits <- c(d.max, d.min)
  age.limits <- c(age.min, age.max)
  if(length(d.min) == 0)
    d.min <- set$d.min
  if(length(d.max) == 0)
    if(length(accordion) == 2)
      d.max <- stretch(set$d.max, accordion[1], accordion[2]) else
        d.max <- set$d.max
  if(length(d.by) == 0)
    d.by <- set$d.by
  if(depths.file) {
    dfile <- paste0(set$coredir, set$core, "/", set$core, "_depths.txt")
    if(!file.exists(dfile))
      stop("I cannot find the file ", paste0(set$coredir, set$core, "/", set$core, "_depths.txt"), call.=FALSE)
    depths <- fastread(dfile, header=FALSE)[,1]
    if(!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers", call.=FALSE)
  }
  if(length(depths) > 0)
    d <- sort(depths) else
      d <- seq(set$d.min, set$d.max, by=d.by) # not d.min itself as depths < set$d.min cannot be calculated. Same for d.max, best not extrapolate here
  # but if squeezing/stretching has to be done, then:
  if(length(accordion) == 2) {
    d <- seq(set$d.min, stretch(set$d.max, accordion[1], accordion[2]), by=d.by)
    d <- squeeze(d, accordion[1], accordion[2])
  }

  if(length(d) > maxcalc)
    message("Warning, this will take quite some time to calculate. I suggest increasing d.by to, e.g., ", 10*d.by, "\n") # was set$d.by

  for(i in set$hiatus.depths)
    d <- sort(unique(c(i+after, i, d)))
  for(i in set$slump)
    d <- sort(unique(c(i+after, i, d)))
  
  if(set$isplum) # new May 2021
    if(set$ra.case == 0) # renamed from incorrectly named radon.case
      if(!set$hasBaconData) # May 2021
        d <- d[which(d <= max(set$detsOrig[,2]))]

  if(verbose && !use.cpp)
    message("Calculating age ranges...")
  
  modelranges <- c()
  if(length(rounded) == 0)
    rounded <- ifelse(set$isplum, 1, 0) # plum tends to need higher precision
  ranges <- ageranges(d, paste0(set$prefix, "_ages.txt"), verbose=verbose,
    set=set, BCAD=BCAD, prob=prob, roundby=rounded, use.cpp=use.cpp)
  d.rng <- d

  modelranges <- range(ranges[,-1], na.rm=TRUE)
  if(length(set$calib$probs) > 0) {
    dates <- set$calib$probs
  dateranges <- c()
  for(i in 1:length(dates))
    if(BCAD && !set$BCAD) # we want to plot in BCAD but the original is in cal BP
      dateranges <- range(dateranges, 1950-dates[[i]][,1], na.rm=TRUE) else
        dateranges <- range(dateranges, dates[[i]][,1], na.rm=TRUE)
  } else dateranges <- modelranges # assuming plum with no additional dates

 if(length(age.min) == 0)
    age.min <- min(modelranges, dateranges)
  if(length(age.max) == 0)
    age.max <- max(modelranges, dateranges)
  age.lim <- extendrange(c(age.min, age.max), f=0.01)
  # then also adapt depth axis lim?
  
  if(BCAD)
    age.lim <- rev(age.lim)
  if(rev.age)
    age.lim <- rev(age.lim)
  d.lim <- rev(extendrange(c(d.max, d.min), f=0.01))
  if(rev.d)
    d.lim <- d.lim[2:1]

  # recalculate calendar axis limits of the plots
  if(length(depth.limits) > 0) { # depth range has been provided
    if(length(age.limits) == 0) { # no age ranges provided
      minyr <- approx(ranges[,1], ranges[,2], c(d.max, d.min), rule=2)$y
      maxyr <- approx(ranges[,1], ranges[,3], c(d.max, d.min), rule=2)$y
      # has to find age range of the relevant dates
      inside <- which(set$calib$d >= min(d.min, d.max, na.rm=TRUE) & set$calib$d <= max(d.min, d.max, na.rm=TRUE))
      date.ranges <- c()
      if(length(inside)>0)
        for(i in inside)
          date.ranges <- range(date.ranges, set$calib$probs[[i]][,1], na.rm=TRUE)
      age.limits <- range(minyr, maxyr, date.ranges)
      depth.limits <- c(d.max, d.min)
    } else {
        age.limits <- c(age.min, age.max)
        depth.limits <- c(d.max, d.min)
      }
  } else {
      if(length(age.limits) == 0) {
        depth.limits <- c(d.max, d.min)
        if(length(set$calib) == 0) # we have (non-Pb) dates
          dates.limits <- c() else
            dates.limits <- range(unlist(lapply(set$calib$probs, function(x) x[,1])))
        if(BCAD)
          dates.limits <- rice::BCADtocalBP(dates.limits)
        age.limits <- range(ranges[,2:3], dates.limits)
      } else {
          depth.limits <- approx(ranges[,4], ranges[,1], c(age.max, age.min), rule=2)$y
          age.limits <- c(min(age.min, max(ranges[,4]), na.rm=TRUE), max(age.max, min(ranges[,4]), na.rm=TRUE))
          if(BCAD)
           age.limits <- rev(age.limits)
        }
    }
  if(BCAD)
    age.limits <- rev(age.limits)
  if(rev.age)
    age.limits <- rev(age.limits)
  if(rev.d)
    depth.limits <- rev(depth.limits)

  if(length(d.lab) == 0)
    d.lab <- paste0("Depth (", depth.unit, ")")
  if(length(age.lab) == 0)
    age.lab <- ifelse(BCAD, ifelse(min(age.limits >= 0), "cal AD", "cal BC/AD"),
      ifelse(kcal, "kcal BP", paste("cal", age.unit, "BP")))

  if(kcal)
    ifelse(rotate.axes, xaxt <- "n", yaxt <- "n")

  if(set$isplum)
    if(remove.tail) {
      # has to check for non-210Pb dates below the tail
      above <- which(set$background < pb.background)
      if(set$hasBaconData) 
        d.lim <- range(set$dets[above,4], set$detsBacon[,4])[2:1] else
          d.lim[which(d.lim==max(d.lim))] <- max(set$dets[above,4])
    }

  pdf.fl <- paste0(set$prefix, ".pdf")
  if(!dev.interactive() && !isTRUE(getOption("knitr.in.progress"))) {
    if(capabilities("aqua") && quartz) # macOS
      grDevices::quartz(file=pdf.fl, type="pdf") else 
        if(capabilities("cairo") && cairo)
          cairo_pdf(filename=pdf.fl) else 
            pdf(file=pdf.fl)
  }

  if(rotate.axes)
    plot(0, type="n", ylim=depth.limits, xlim=age.limits, ylab=d.lab, xlab=age.lab, bty="n", xaxt=xaxt, yaxt=yaxt, mar=mar.main, cex.axis=mainpanel.tickfontsize, cex.lab=mainpanel.labelfontsize) else
      plot(0, type="n", xlim=depth.limits[2:1], ylim=age.limits, xlab=d.lab, ylab=age.lab, bty="n", xaxt=xaxt, yaxt=yaxt, mar=mar.main, cex.axis=mainpanel.tickfontsize, cex.lab=mainpanel.labelfontsize)
  if(kcal)
    axis(ifelse(rotate.axes, 1, 2), pretty(age.lim), pretty(age.lim/1e3), cex.axis=mainpanel.labelfontsize)

  if(set$isplum) { # new May 2021
    top <- min(set$dets[,4]-set$dets[,5])
    if(d.min > top)
      d.min <- top # because the top 210Pb-dated depth should be at or below d.min
    if(!set$hasBaconData) { # needs more work
      above <- which(set$background < pb.background)
      dmax <- max(set$dets[above,4]) # cut depths that have reached background
      if(remove.tail)
        if(!is.infinite(dmax)) {
          ranges <- ranges[which(d.rng <= dmax),]
          d <- d[which(d <= dmax)]
          d.max <- dmax
        }
    }
  }

  if(!dates.only) {
    if(verbose && !use.cpp)
      message("Preparing ghost graph... ")
    agedepth.ghost(set, rotate.axes=rotate.axes, accordion=accordion, BCAD=BCAD, d.min=d.min, d.max=d.max, d.res=d.res, age.res=age.res, rev.d=rev.d, rev.age=rev.age, rgb.res=rgb.res, dark=dark, from.col=from.col, to.col=to.col, rgb.scale=rgb.scale, age.lim=age.limits, use.raster=use.raster, flip.age=flip.age, flip.d=flip.d, verbose=verbose, use.cpp=use.cpp)
  }

  if(length(set$slump) > 0 )
    for(i in 1:nrow(set$slump))
      if(rotate.axes)
        rect(min(age.lim)-1e3, set$slump[i,1], max(age.lim)+1e3, set$slump[i,2], col=slump.col, border=slump.col) else
          rect(set$slump[i,1], min(age.lim)-1e3, set$slump[i,2], max(age.lim)+1e3, col=slump.col, border=slump.col)

  if(length(plotatthesedepths) == 0) 
    thesedateddepths <- c() else {
      dseq <- seq(d.min, d.max, length=length(plotatthesedepths)) # this OK? May 2023
      if(set$isplum)
        thesedateddepths <- approx(dseq, plotatthesedepths, set$detsBacon[,4])$y else
          thesedateddepths <- approx(dseq, plotatthesedepths, set$dets[,4])$y
    }

  if(set$isplum) {
    if(set$hasBaconData)
      calib.plot(set, dets=set$detsBacon, accordion=accordion, BCAD=BCAD, cc=cc, rotate.axes=rotate.axes, height=height, calheight=calheight, ex=ex, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, dates.col=dates.col, new.plot=FALSE, same.heights=same.heights)

    if(plot.pb)
      set <- draw.pbmodelled(set, BCAD=BCAD, rotate.axes=rotate.axes, age.lim=age.lim, d.lim=c(d.min, d.max), pbmodelled.col=pbmodelled.col, pbmeasured.col=pbmeasured.col, pb.lim=pb.lim, supp.col=supp.col, mgp=mgp, pb.lty=pb.lty, save.info=save.info) # pointing to set since April 2025 
  } else
    calib.plot(set, dets=set$dets, accordion=accordion, BCAD=BCAD, cc=cc, rotate.axes=rotate.axes, height=height, calheight=calheight, ex=ex, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, dates.col=dates.col, new.plot=FALSE, same.heights=same.heights)

  if(plot.labels) {
    if(length(labels) == 0)
      labels <- set$dets[,1]
    for(i in 1:nrow(set$dets)) {
      age <- set$calib$probs[[i]][,1]
      if(BCAD)
        age <- rice::calBPtoBCAD(age)
      if(label.age == 1)
        age.pos <- max(age) else
          age.pos <- min(age)
      d <- set$calib$d[i]
      if(rotate.axes)
         text(age.pos+label.offset[1], d+label.offset[2], labels[i], cex=label.size, col=label.col, adj=label.adj, srt=label.rot) else
           text(d+label.offset[1], age.pos+label.offset[2], labels[i], cex=label.size, col=label.col, adj=label.adj, srt=label.rot)
    }
  }
  legend(title.location, title, bty="n", cex=title.size)
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
    th <- array(sort(c(1, nrow(ranges), hi.d, hi.d)), dim=c(2,length(hi.d)+1))
  }

  if(length(accordion) == 2)
    d <- stretch(d, accordion[1], accordion[2])

  if(!dates.only)
    for(i in 1:ncol(th)) {
      h <- th[1,i] : th[2,i]
      if(rotate.axes) {
        lines(ranges[h,2], d[h], col=range.col, lty=range.lty, lwd=range.lwd)
        lines(ranges[h,3], d[h], col=range.col, lty=range.lty, lwd=range.lwd)
        lines(ranges[h,4], d[h], col=med.col, lty=med.lty, lwd=med.lwd) # median
        lines(ranges[h,5], d[h], col=mn.col, lty=mn.lty, lwd=mn.lwd) # mean
      } else {
          lines(d[h], ranges[h,2], col=range.col, lty=range.lty, lwd=range.lwd)
          lines(d[h], ranges[h,3], col=range.col, lty=range.lty, lwd=range.lwd)
          lines(d[h], ranges[h,4], col=med.col, lty=med.lty, lwd=med.lwd) # median
          lines(d[h], ranges[h,5], col=mn.col, lty=mn.lty, lwd=mn.lwd) # mean
        }
    }

  set$ranges <- ranges

  opened.device <- FALSE
  if(!dev.interactive() && !isTRUE(getOption("knitr.in.progress"))) {
    pdf(file=pdf.fl)
    opened.device <- TRUE
  }
  if(plot.pdf) {
    if(dev.interactive())
      grDevices::dev.copy2pdf(file=pdf.fl) else
        if(opened.device)
         dev.off()
  }

  rng <- abs(ranges[,3]-ranges[,2])
  min.rng <- d[which(rng==min(rng, na.rm=TRUE))]
  max.rng <- d[which(rng==max(rng, na.rm=TRUE))]
  if(length(min.rng)==1)
    min.rng <- paste(age.unit, "at", min.rng, noquote(depth.unit)) else
      min.rng <- paste(age.unit, "between", min(min.rng), "and", max(min.rng), noquote(depth.unit))
  if(length(max.rng)==1)
    max.rng <- paste(age.unit, "at", max.rng, noquote(depth.unit)) else
      max.rng <- paste(age.unit, "between", min(max.rng), "and", max(max.rng), noquote(depth.unit))

  if(verbose) {
    if(!dates.only)
      message("\nMean ", 100*prob, "% confidence ranges ", 
        round(mean(rng), rounded), " ", age.unit, ", min. ",
        min(rng), " ", min.rng, ", max. ", max(rng), " ", max.rng)

    if(set$isplum) {
      set$overlap <- model.Pb.hpd(set, verbose=TRUE, decimals=rounded)
      if(set$hasBaconData) { # then also report the overlap with non-210Pb data
        overlap.hpds <- model.dates.hpd(set, verbose=verbose)
      } else overlap.hpds <- set$overlap
      below <- which(bg > pb.background)
      if(length(below) == 0)
        message("It seems that the Pb-measurements haven't reached background") else {
          message("Pb-210 likely at or below detection limit from ", 
            min(set$dets[below,4]), " ", set$depth.unit, " depth onward: ")
          for(i in seq_along(below))
            cat(set$dets[below[i],4], " ", set$depth.unit,
              " (", round(100*bg[below[i]]), "%)",
              if(i < length(below)) ", " else "\n", sep = "")
        }
    } else
       overlap.hpds <- model.dates.hpd(set, verbose=verbose)

    # report summaries of posteriors
    posteriors <- list(acc=set$post.acc, mem=set$post.mem)
    if(!is.na(hiatus.depths[1]))
      if(is.na(set$boundary[1])) # it's a hiatus, not a boundary
        posteriors$hiatus <- set$post.hiatus # probably also allow for multiple posteriors here
    posteriors <- lapply(posteriors, round, digits=roundby)
    if(!is.na(hiatus.depths[1])) {
      if(is.na(set$boundary[1])) {
        message("Posteriors: accrate mean ", paste0(posteriors$acc[,1], collapse=" & "),
          ", shape ", paste0(posteriors$acc[,2], collapse=" & ",
          ", memory mean ", posteriors$mem[1], ", strength ", posteriors$mem[2],
          ", hiatus mean ", posteriors$hiatus[1], ", shape ", posteriors$hiatus[2]))
        }
    } else
      message("Posteriors: accrate mean ", posteriors$acc[1], ", shape ", posteriors$acc[2],
        ", memory mean ", posteriors$mem[1], ", strength ", posteriors$mem[2])

    if(set$isplum) {
      plumpost <- round(c(set$post.phi, set$post.supp), roundby)
      message("phi mean ", plumpost[1], ", shape ", plumpost[2], 
        ", supported mean ", plumpost[3], ", shape ", plumpost[4])
      overlap.hpds <- model.Pb.hpd(set, verbose=verbose)
    }
  } else {
      overlap.hpds <- model.dates.hpd(set, verbose=verbose)
      #overlap <- as.numeric(overlap.hpds[1])
    } 
  set$overlap.hpds <- overlap.hpds

  # write a summary file to the core's directory
  if(write.summary) {
    learned <- learning(set, talk=verbose)
    summarise.run(set, roundby, overlap.hpds[2], learn=learned)
  }

  if(save.info)
    assign_to_global("info", set)
  
  #par(oldpar) # tmp Jan 2021
  invisible(set)
}
