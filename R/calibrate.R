
#' @name copyCalibrationCurve
#' @title Copy a calibration curve.
#' @description Copy one of the the calibration curves into memory. 
#' @details Copy the radiocarbon calibration curve defined by cc into memory. 
#' @return The calibration curve (invisible).
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), 
#' \code{cc=3} for SHCal13 (southern hemisphere terrestrial). 
#' @param postbomb Use \code{postbomb=TRUE} to get a postbomb calibration curve (default \code{postbbomb=FALSE}).
#' @author Maarten Blaauw, J. Andres Christen
#' @examples 
#' intcal13 <- copyCalibrationCurve(1)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @export
copyCalibrationCurve <- function(cc=1, postbomb=FALSE) {
  if(postbomb) {
    if(cc==1) fl <- "postbomb_NH1.14C" else
      if(cc==2) fl <- "postbomb_NH2.14C" else
        if(cc==3) fl <- "postbomb_NH3.14C" else
          if(cc==4) fl <- "postbomb_SH1-2.14C" else 
            if(cc==5) fl <- "postbomb_SH3.14C" else  
              stop("Calibration curve doesn't exist\n")		       
  } else
  if(cc==1) fl <- "3Col_intcal13.14C" else
    if(cc==2) fl <- "3Col_marine13.14C" else
      if(cc==3) fl <- "3Col_shcal13.14C" else
        stop("Calibration curve doesn't exist\n")  
  cc <- system.file("extdata/Curves", fl, package='rbacon')
  cc <- read.table(cc)
  invisible(cc)
}



#' @name mix.curves
#' @title Build a custom-made, mixed calibration curve. 
#' @description If two curves need to be `mixed' to calibrate, e.g. for dates of mixed terrestrial and marine carbon sources, then this function can be used. 
#' @details The proportional contribution of each of both calibration curves has to be set. 
#'
#' @param proportion Proportion of the first calibration curve required. e.g., change to \code{proportion=0.7} if \code{cc1} should contribute 70\% (and \code{cc2} 30\%) to the mixed curve.
#' @param cc1 The first calibration curve to be mixed. Defaults to the northern hemisphere terrestrial curve IntCal13.
#' @param cc2 The second calibration curve to be mixed. Defaults to the marine curve IntCal13.
#' @param name Name of the new calibration curve.
#' @param dirname Directory where the file will be written. If using the default \code{dirname="."}, 
#' the new curve will be saved in current working directory. 
#' @param offset Any offset and error to be applied to \code{cc2} (default 0 +- 0).
#' @param sep Separator between fields (tab by default, "\\t") 
#' @author Maarten Blaauw, J. Andres Christen
#' @return A file containing the custom-made calibration curve, based on calibration curves \code{cc1} and \code{cc2}.
#' @examples
#' mix.curves()
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
mix.curves <- function(proportion=.5, cc1="3Col_intcal13.14C", cc2="3Col_marine13.14C", name="mixed.14C", dirname=".", offset=c(0,0), sep="\t") {
  ccloc <- paste(system.file("extdata", package='rbacon'), "/Curves/", sep="") 
  dirname <- .validateDirectoryName(dirname)
  
  cc1 <- read.table(paste(ccloc, cc1,  sep=""))
  cc2 <- read.table(paste(ccloc, cc2,  sep=""))
  cc2.mu <- approx(cc2[,1], cc2[,2], cc1[,1], rule=2)$y + offset[1] # interpolate cc2 to the calendar years of cc1
  cc2.error <- approx(cc2[,1], cc2[,3], cc1[,1], rule=2)$y
  cc2.error <- sqrt(cc2.error^2 + offset[2]^2)
  mu <- proportion * cc1[,2] + (1-proportion) * cc2.mu
  error <- proportion * cc1[,3] + (1-proportion) * cc2.error
  write.table(cbind(cc1[,1], mu, error), paste(dirname, name,  sep=""), row.names=FALSE, col.names=FALSE, sep=sep)
}



#' @name pMC.age
#' @title Calculate C14 ages from pMC values.
#' @description Calculate C14 ages from pMC values of radiocarbon dates.
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from pMC values. The reverse function is \link{age.pMC}.
#' @param mn Reported mean of the pMC.
#' @param sdev Reported error of the pMC.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, use \code{1} here.
#' @param decimals Amount of decimals required for the radiocarbon age.
#' @author Maarten Blaauw, J. Andres Christen
#' @return Radiocarbon ages from pMC values. If pMC values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   pMC.age(110, 0.5) # a postbomb date, so with a negative 14C age
#'   pMC.age(80, 0.5) # prebomb dates can also be calculated
#'   pMC.age(.8, 0.005, 1) # pMC expressed against 1 (not against 100\%)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#'  \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
pMC.age <- function(mn, sdev, ratio=100, decimals=0) {
  y <- -8033 * log(mn/ratio)
  sdev <- y - -8033 * log((mn+sdev)/ratio)
  round(c(y, sdev), decimals)
}



#' @name age.pMC
#' @title Calculate pMC values from C14 ages
#' @description Calculate pMC values from radiocarbon ages
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages, 
#' this function can be used to calculate pMC values from radiocarbon ages. The reverse function of \link{pMC.age}.
#' @param mn Reported mean of the 14C age.
#' @param sdev Reported error of the 14C age.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, use \code{1} here.
#' @param decimals Amount of decimals required for the pMC value.
#' @author Maarten Blaauw, J. Andres Christen
#' @return pMC values from C14 ages.
#' @examples
#'   age.pMC(-2000, 20)
#'   age.pMC(-2000, 20, 1)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
age.pMC <- function(mn, sdev, ratio=100, decimals=3) {
  y <- exp(-mn / 8033)
  sdev <- y - exp(-(mn + sdev) / 8033)
  signif(ratio*c(y, sdev), decimals)
}



#' @name add.dates
#' @title Add dates to age-depth plots
#' @description Add dated depths to plots, e.g. to show dates that weren't used in the age-depth model
#' @details Sometimes it is useful to add additional dating information to age-depth plots, e.g., to show outliers or how dates calibrate with different estimated offsets. 
#' @param mn Reported mean of the date. Can be multiple dates. 
#' @param sdev Reported error of the date. Can be multiple dates. 
#' @param depth Depth of the date. 
#' @param cc The calibration curve to use: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), \code{cc=0} for none (dates that are already on the cal BP scale).
#' @param above Treshold for plotting of probability values. Defaults to \code{above=1e-3}.
#' @param ex Exaggeration of probability distribution plots. Defaults to \code{ex=50}.
#' @param normal By default, Bacon uses the student's t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param normalise By default, the date is normalised to an area of 1 (\code{normalise=TRUE}). 
#' @param t.a The dates are treated using the student's t distribution by default (\code{normal=FALSE}). 
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010). 
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file). 
#' For symmetry reasons, t.a must always be equal to t.b-1. 
#' @param t.b The dates are treated using the student's t distribution by default (\code{normal=FALSE}). 
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010). 
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file). 
#' For symmetry reasons, t.a must always be equal to t.b-1. 
#' @param age.res Resolution of the date's distribution. Defaults to \code{date.res=100}.
#' @param times The extent of the range to be calculated for each date. Defaults to \code{times=20}.
#' @param col The colour of the ranges of the date. Default is semi-transparent red: \code{col=rgb(1,0,0,.5)}.
#' @param border The colours of the borders of the date. Default is semi-transparent red: \code{border=rgb(1,0,0,0.5)}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param pch The shape of any marker to be added to the date. Defaults to a cross, \code{pch=4}. To leave empty, use \code{pch=NA}.  
#' @author Maarten Blaauw, J. Andres Christen
#' @return A date's distribution, added to an age-depth plot.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth()
#'   add.dates(5000, 100, 60)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
add.dates <- function(mn, sdev, depth, cc=1, above=1e-3, ex=10, normal=TRUE, normalise=TRUE, t.a=3, t.b=4, age.res=100, times=20, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5), rotate.axes=FALSE, mirror=TRUE, up=TRUE, BCAD=FALSE, pch=4) {
  if(cc > 0)
    cc <- copyCalibrationCurve(cc)
  
  for(i in 1:length(mn)) {
    yrs <- seq(mn[i]-times*sdev[i], mn[i]+times*sdev[i], length=age.res)
    if(length(cc) < 2)
      cc <- cbind(yrs, yrs, rep(0, length(yrs)))
    ages <- approx(cc[,1], cc[,2], yrs)$y
    errors <- approx(cc[,1], cc[,3], yrs)$y

  if(normal)
    probs <- dnorm(ages, mn[i], sqrt(sdev[i] + errors)^2) else
      probs <- (t.b + (mn[i]-ages)^2  / (2*(sdev[i]^2 + errors^2))) ^ (-1*(t.a+0.5))
  if(normalise)
    probs <- probs / sum(probs)
  these <- which(probs >= above)
  if(length(these) > 0) {
    yrs <- yrs[these]
    probs <- probs[these]
  }

  if(!up)
    up <- -1  
  if(BCAD)
    yrs <- 1950 - yrs
  if(mirror)
    pol <- cbind(depth[i] + ex*c(probs, -rev(probs)), c(yrs, rev(yrs))) else
      pol <- cbind(depth[i] - up*ex*c(0, probs,  0), c(min(yrs), yrs, max(yrs)))
  if(rotate.axes)
    pol <- pol[,2:1]
  polygon(pol, col=col, border=border)
  if(length(pch) > 0)
    if(rotate.axes)
      points(mean(yrs), depth[i], col=border, pch=pch) else
        points(depth[i], mean(yrs), col=border, pch=pch)
  }
}



#' @name calib.plot
#' @title Plot the dates
#' @description Produce a plot of the dated depths and their dates
#' @details This function is generally called internally to produce the age-depth graph. 
#' It can be used to produce custom-built graphs. 
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}). 
#' @param rev.yr The direction of the calendar age axis can be reversed from the default (\code{rev.yr=TRUE}) 
#' @param yr.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=c()}).
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param height The heights of the distributions of the dates. See also \code{normalise.dists}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
#' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
#' @param C14.col Colour of the calibrated distributions of the dates. Default is semi-transparent blue: \code{rgb(0,0,1,.35)}.
#' @param C14.border Colours of the borders of calibrated 14C dates. Default is transparent dark blue: cal.col
#' @param cal.col Colour of the non-14C dates in the age-depth plot: default semi-transparent blue-green: \code{rgb(0,.5,.5,.35)}.
#' @param cal.border Colour of the of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{rgb(0,.5,.5,.5)}.
#' @param dates.col As an alternative to colouring dates based on whether they are 14C or not, sets of dates can be coloured as, e.g., \code{dates.col=colours()[2:100]}.
#' @param slump.col Colour of slumps. Defaults to \code{slump.col=grey(0.8)}.
#' @param new.plot Start a new plot (\code{new.plot=TRUE}) or plot over an existing plot (\code{new.plot=FALSE}).
#' @param plot.dists Plot the distributions of the dates (default \code{plot.dists=TRUE}).
#' @param normalise.dists By default, the distributions of more precise dates will cover less time and will thus peak higher than less precise dates. This can be avoided by specifying \code{normalise.dists=FALSE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   calib.plot()
#'
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
### produce plots of the calibrated distributions
calib.plot <- function(set=get('info'), BCAD=set$BCAD, rotate.axes=FALSE, rev.d=FALSE, rev.yr=FALSE, yr.lim=c(), date.res=100, d.lab=c(), yr.lab=c(), height=15, mirror=TRUE, up=TRUE, cutoff=.001, C14.col=rgb(0,0,1,.5), C14.border=rgb(0,0,1,.75), cal.col=rgb(0,.5,.5,.5), cal.border=rgb(0,.5,.5,.75), dates.col=c(), slump.col=grey(0.8), new.plot=TRUE, plot.dists=TRUE, normalise.dists=TRUE) {
  height <- length(set$d.min:set$d.max) * height/50
  if(length(yr.lim) == 0)
    lims <- c()  
  for(i in 1:length(set$calib$probs))
    lims <- c(lims, set$calib$probs[[i]][,1])
  yr.min <- min(lims)
  yr.max <- max(lims)
  if(BCAD) {
    yr.min <- 1950 - yr.min
    yr.max <- 1950 - yr.max
    } 
  if(length(yr.lab) == 0)
    yr.lab <- ifelse(set$BCAD, "BC/AD", "cal yr BP")
  yr.lim <- c(yr.min, yr.max)  
  if(rev.yr)
    yr.lim <- yr.lim[2:1]
  dlim <- range(set$d)
  if(rev.d)
    dlim <- dlim[2:1]
  if(length(d.lab) == 0)
    d.lab <- paste("depth (", set$unit, ")", sep="")  
  if(new.plot)
    if(rotate.axes)
      plot(0, type="n", xlim=yr.lim, ylim=dlim[2:1], xlab=yr.lab, ylab=d.lab, main="") else
        plot(0, type="n", xlim=dlim, ylim=yr.lim, xlab=d.lab, ylab=yr.lab, main="")

  if(length(set$slump) > 0)
    if(rotate.axes)
      abline(h=set$slump, lty=2, col=slump.col) else
        abline(v=set$slump, lty=2, col=slump.col)

  if(plot.dists)
    for(i in 1:length(set$calib$probs)) {
      cal <- set$calib$probs[[i]]
        d <- set$calib$d[[i]]
      if(BCAD)
        cal[,1] <- 1950-cal[,1]
      o <- order(cal[,1])
      if(normalise.dists)
        cal <- cbind(cal[o,1], height*cal[o,2]/sum(cal[,2])) else
          cal <- cbind(cal[o,1], height*cal[o,2]/max(cal[,2]))
      cal <- cal[cal[,2] >= cutoff,]  
      cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), length=200)) # tmp
      if(mirror)
        pol <- cbind(c(d-cal$y, d+rev(cal$y)), c(cal$x, rev(cal$x))) else
         if(up)
           pol <- cbind(d-c(0, cal$y, 0), c(min(cal$x), cal$x, max(cal$x))) else
             pol <- cbind(d+c(0, cal$y, 0), c(min(cal$x), cal$x, max(cal$x)))
      if(rotate.axes)
        pol <- cbind(pol[,2], pol[,1])
      if(ncol(set$dets)==4 || (ncol(set$dets) > 4 && set$dets[i,5] > 0)) {
        col <- C14.col
        border <- C14.border
      } else {
          col <- cal.col
          border <- cal.border
        }
      if(length(dates.col) > 0) {
		col <- dates.col[i]
        border <- dates.col[i]
	  }
      polygon(pol, col=col, border=border)
    }
}


# calibrate C14 dates and calculate distributions for any calendar dates
.bacon.calib <- function(dat, set=get('info'), date.res=100, normal=set$normal, t.a=set$t.a, t.b=set$t.b, delta.R=set$delta.R, delta.STD=set$delta.STD, ccdir="") {
  # read in the curves
  if(set$cc1=="IntCal13" || set$cc1=="\"IntCal13\"")
    cc1 <- read.table(paste(ccdir, "3Col_intcal13.14C",sep="")) else
      cc1 <- read.csv(paste(ccdir, set$cc1, ".14C", sep=""), header=FALSE, skip=11)[,1:3]
  if(set$cc2=="Marine13" || set$cc2=="\"Marine13\"")
    cc2 <- read.table(paste(ccdir, "3Col_marine13.14C",sep="")) else
      cc2 <- read.csv(paste(ccdir, set$cc2, ".14C", sep=""), header=FALSE, skip=11)[,1:3]
  if(set$cc3=="SHCal13" || set$cc3=="\"SHCal13\"")
    cc3 <- read.table(paste(ccdir, "3Col_shcal13.14C",sep="")) else
      cc3 <- read.csv(paste(ccdir, set$cc3, ".14C", sep=""), header=FALSE, skip=11)[,1:3]
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") cc4 <- NA else
    cc4 <- read.table(paste(ccdir, set$cc4, sep=""))[,1:3]

  if(set$postbomb != 0) {
    if(set$postbomb==1) bomb <- read.table(paste(ccdir,"postbomb_NH1.14C", sep=""))[,1:3] else
      if(set$postbomb==2) bomb <- read.table(paste(ccdir,"postbomb_NH2.14C", sep=""))[,1:3] else
        if(set$postbomb==3) bomb <- read.table(paste(ccdir,"postbomb_NH3.14C", sep=""))[,1:3] else
          if(set$postbomb==4) bomb <- read.table(paste(ccdir,"postbomb_SH1-2.14C", sep=""))[,1:3] else
            if(set$postbomb==5) bomb <- read.table(paste(ccdir,"postbomb_SH3.14C", sep=""))[,1:3] else
              stop("Warning, cannot find postbomb curve #", set$postbomb, " (use values of 1 to 5 only)")
      bomb.x <- seq(max(bomb[,1]), min(bomb[,1]), by=-.1) # interpolate
      bomb.y <- approx(bomb[,1], bomb[,2], bomb.x)$y
      bomb.z <- approx(bomb[,1], bomb[,3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level=0)
      if(set$postbomb < 4)
        cc1 <- rbind(bomb, cc1, deparse.level=0) else
          cc3 <- rbind(bomb, cc3, deparse.level=0)
  }

  ## use Gaussian or t (Christen and Perez Radiocarbon 2009) calibration
  if(round(set$t.b-set$t.a) !=1)
    stop("\n Warning! t.b - t.a should always be 1, check the manual")
  d.cal <- function(cc, rcmean, w2, t.a, t.b) {
    if(set$normal)
      cal <- cbind(cc[,1], dnorm(cc[,2], rcmean, sqrt(cc[,3]^2+w2))) else
        cal <- cbind(cc[,1], (t.b+ ((rcmean-cc[,2])^2) / (2*(cc[,3]^2 + w2))) ^ (-1*(t.a+0.5)))
    cal[,2] <- cal[,2]/sum(cal[,2])
    if(length(which(cal[,2]>set$cutoff)) > 5) # ensure that also very precise dates get a range of probabilities
      cal[which(cal[,2]>set$cutoff),] else {
        calx <- seq(min(cal[,1]), max(cal[,1]), length=50)
        caly <- approx(cal[,1], cal[,2], calx)$y
        cbind(calx, caly/sum(caly))
      }
  }

  # now calibrate all dates
  calib <- list(d=dat[,4])
  if(ncol(dat)==4) { # only one type of dates (e.g., calBP, or all IntCal13 C14 dates)
    if(set$cc==0) {
      delta.R <- 0; delta.STD <- 0 # only C14 dates should need correcting for age offsets
      x <- seq(min(dat[,2])-max(100,4*max(dat[,3])), max(dat[,2])+max(100,4*max(dat[,3])), length=date.res)
      ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
    } else
        if(set$cc==1) ccurve <- cc1 else
          if(set$cc==2) ccurve <- cc2 else
            if(set$cc==3) ccurve <- cc3 else
              ccurve <- cc4
    for(i in 1:nrow(dat))
      calib$probs[[i]] <- d.cal(ccurve, dat[i,2]-delta.R, dat[i,3]^2+delta.STD^2, set$t.a, set$t.b)
  } else
      for(i in 1:nrow(dat)) {
        det <- as.numeric(dat[i,])
        if(det[5]==0) {
          x <- seq(det[2]-max(100,4*det[3]), det[2]+max(100,4*det[3]), length=date.res)
          ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
        } else
            if(det[5]==1) ccurve <- cc1 else if(det[5]==2) ccurve <- cc2 else
              if(det[5]==3) ccurve <- cc3 else ccurve <- cc4

        delta.R <- set$delta.R; delta.STD <- set$delta.STD; t.a <- set$t.a; t.b <- set$t.b
        if(length(det) >= 7 && det[5] > 0) { # the user provided age offsets; only for C14 dates
          delta.R <- det[6]
          delta.STD <- det[7]
        }

        if(length(det) >= 9) { # the user provided t.a and t.b values for each date
          t.a <- det[8]
          t.b <- det[9]
          if(round(t.b-t.a) !=1) 
            stop("\n Warning! t.b - t.a should always be 1, check the manual")
        }
        calib$probs[[i]] <- d.cal(ccurve, det[2]-delta.R, det[3]^2+delta.STD^2, t.a, t.b)
      }
  calib
}

