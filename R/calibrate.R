


#' @name add.dates
#' @title Add dates to age-depth plots
#' @description Add dated depths to plots, e.g. to show dates that weren't used in the age-depth model
#' @details Sometimes it is useful to add additional dating information to age-depth plots, e.g., to show outliers or how dates calibrate with different estimated offsets.
#' @param mn Reported mean of the date. Can be multiple dates. Negative numbers indicate postbomb dates (if cc > 0).
#' @param sdev Reported error of the date. Can be multiple dates.
#' @param depth Depth of the date.
#' @param cc The calibration curve to use: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine), \code{cc=0} for none (dates that are already on the cal BP scale).
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param above Treshold for plotting of probability values. Defaults to \code{above=1e-3}.
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}
#' @param normal By default, Bacon uses the student's t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param delta.R Mean of core-wide age offsets (e.g., regional marine offsets).
#' @param delta.STD Error of core-wide age offsets (e.g., regional marine offsets).
#' @param t.a The dates are treated using the student's t distribution by default (\code{normal=FALSE}).
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param t.b The dates are treated using the student's t distribution by default (\code{normal=FALSE}).
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param date.res Resolution of the date's distribution. Defaults to \code{date.res=100}.
#' @param height The heights of the distributions of the dates. See also \code{normalise.dists}.
#' @param calheight Multiplier for the heights of the distributions of dates on the calendar scale. Defaults to \code{calheight=1}.
#' @param agesteps Step size for age units of the distribution. Default \code{agesteps=1}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.005}).
#' @param col The colour of the ranges of the date. Default is semi-transparent red: \code{col=rgb(1,0,0,.5)}.
#' @param border The colours of the borders of the date. Default is semi-transparent red: \code{border=rgb(1,0,0,0.5)}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param pch The shape of any marker to be added to the date. Defaults to a cross, \code{pch=4}. To leave empty, use \code{pch=NA}.
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{ccdir=""}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A date's distribution, added to an age-depth plot.
#' @examples
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth()
#'   add.dates(5000, 100, 60)
#' }
#' @export
add.dates <- function(mn, sdev, depth, cc=1, set=get('info'), above=1e-6, postbomb=0, normal=TRUE, delta.R=set$delta.R, delta.STD=set$delta.STD, t.a=set$t.a, t.b=set$t.b, date.res=100, height=1, calheight=1, agesteps=1, cutoff=0.005, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5), rotate.axes=FALSE, mirror=TRUE, up=TRUE, BCAD=FALSE, pch=4, ccdir="") {
  if(ccdir == "")
    ccdir <- system.file("extdata", package="IntCal")
  ccdir <- validateDirectoryName(ccdir)

  if(mn < 0 && cc > 0)
    if(postbomb == 0)
      stop("Negative C-14 age, please provide a postbomb curve, e.g. postbomb=1")
  
  dat <- cbind(mn, mn, sdev, depth, cc, delta.R, delta.STD, t.a, t.b)
  probs <- bacon.calib(dat, set, date.res, cutoff, postbomb, normal, t.a, t.b, delta.R, delta.STD, ccdir) 

  for(i in 1:nrow(dat)) {
    d <- probs$d[[i]]
    yrs <- probs$probs[[i]][,1]
    cal <- probs$probs[[i]][,2]

    if(BCAD)
      yrs <- 1950-yrs
    cal <- approx(yrs, cal, seq(min(yrs), max(yrs), by=agesteps))
    cal <- cbind(cal$x, cal$y/sum(cal$y))
    cal[,2] <- 3 * height * cal[,2] * (set$d.max - set$d.min) # scale the height of the blobs with the depth axis
    if(dat[i,5] == 0) # cal BP dates could have different relative height:
      cal[,2] <- calheight * cal[,2]

  if(!up)
    up <- -1
  if(mirror)
    pol <- cbind(d + c(cal[,2], -rev(cal[,2])), c(cal[,1], rev(cal[,1]))) else
      pol <- cbind(d - up*c(0, cal[,2],  0), c(min(cal[,1]), cal[,1], max(cal[,1])))
  if(rotate.axes)
    pol <- pol[,2:1]
  polygon(pol, col=col, border=border)
  if(length(pch) > 0)
    if(rotate.axes)
      points(mean(yrs), d, col=border, pch=pch) else
        points(d, mean(yrs), col=border, pch=pch)
    }
  invisible(probs)
}



#' @name calib.plot
#' @title Plot the dates
#' @description Produce a plot of the dated depths and their dates
#' @details This function is generally called internally to produce the age-depth graph.
#' It can be used to produce custom-built graphs.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param cc Calibration curve to be used (defaults to info$cc)
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
#' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
#' @param rev.yr Deprecated - use rev.age instead
#' @param age.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{age.lim=c()}).
#' @param yr.lim Deprecated - use age.lim instead
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
#' @param age.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.lab Deprecated - use age.lab instead
#' @param height The heights of the distributions of the dates. See also \code{normalise.dists}.
#' @param calheight Multiplier for the heights of the distributions of dates on the calendar scale. Defaults to \code{calheight=1}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.1}).
#' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
#' @param C14.col Colour of the calibrated distributions of the dates. Default is semi-transparent blue: \code{rgb(0,0,1,.35)}.
#' @param C14.border Colours of the borders of calibrated 14C dates. Default is transparent dark blue: cal.col
#' @param cal.col Colour of the non-14C dates in the age-depth plot: default semi-transparent blue-green: \code{rgb(0,.5,.5,.35)}.
#' @param cal.border Colour of the of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{rgb(0,.5,.5,.5)}.
#' @param dates.col As an alternative to colouring dates based on whether they are 14C or not, sets of dates can be coloured as, e.g., \code{dates.col=colours()[2:100]}.
#' @param slump.col Colour of slumps. Defaults to \code{slump.col=grey(0.8)}.
#' @param new.plot Start a new plot (\code{new.plot=TRUE}) or plot over an existing plot (\code{new.plot=FALSE}).
#' @param plot.dists Plot the distributions of the dates (default \code{plot.dists=TRUE}).
#' @param same.heights Plot the distributions of the dates all at the same maximum height (default \code{same.height=FALSE}), which instead normalises the distributions (all have an area of 1).
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   calib.plot()
#' @export
### produce plots of the calibrated distributions
calib.plot <- function(set=get('info'), BCAD=set$BCAD, cc=set$cc, rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, rev.yr=rev.age, age.lim=c(), yr.lim=age.lim, date.res=100, d.lab=c(), age.lab=c(), yr.lab=age.lab, height=1, calheight=1, mirror=TRUE, up=TRUE, cutoff=.1, C14.col=rgb(0,0,1,.5), C14.border=rgb(0,0,1,.75), cal.col=rgb(0,.5,.5,.5), cal.border=rgb(0,.5,.5,.75), dates.col=c(), slump.col=grey(0.8), new.plot=TRUE, plot.dists=TRUE, same.heights=FALSE) {
  if(length(age.lim) == 0)
    lims <- c()
  for(i in 1:length(set$calib$probs))
    lims <- c(lims, set$calib$probs[[i]][,1])
  age.min <- min(lims)
  age.max <- max(lims)
  if(BCAD) {
    age.min <- 1950 - age.min
    age.max <- 1950 - age.max
    }
  if(length(age.lab) == 0)
    age.lab <- ifelse(set$BCAD, "BC/AD", paste("cal", set$age.unit, " BP"))
  age.lim <- extendrange(c(age.min, age.max), f=0.01)
  if(rev.age)
    age.lim <- age.lim[2:1]
  dlim <- range(set$elbows)
  if(rev.d)
    dlim <- dlim[2:1]
  if(length(d.lab) == 0)
    d.lab <- paste("depth (", set$depth.unit, ")", sep="")
  if(new.plot)
    if(rotate.axes)
      plot(0, type="n", xlim=age.lim, ylim=dlim[2:1], xlab=age.lab, ylab=d.lab, main="") else
        plot(0, type="n", xlim=dlim, ylim=age.lim, xlab=d.lab, ylab=age.lab, main="")

  if(length(set$slump) > 0)
    if(rotate.axes)
      abline(h=set$slump, lty=2, col=slump.col) else
        abline(v=set$slump, lty=2, col=slump.col)

  maxhght <- 0; agesteps <- c()
  if(plot.dists) 
    for(i in 1:length(set$calib$probs)) {
      maxhght[i] <- max(set$calib$probs[[i]][,2] / sum(set$calib$probs[[i]][,2]))
      agesteps <- min(agesteps, abs(diff(set$calib$probs[[1]][,1])))
    }

  if(plot.dists)
    for(i in 1:length(set$calib$probs)) {
      cal <- cbind(set$calib$probs[[i]])
      dupl <- which(diff(cal[,1]) == 0)
      if(length(dupl) > 0) # avoid warning of collapsing to unique values
        cal <- cal[-dupl,]
      d <- set$calib$d[[i]]
      if(BCAD)
        cal[,1] <- 1950-cal[,1]
      cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), by=agesteps)) # all dists should have the same binsize
      cal <- cbind(cal$x, cal$y/sum(cal$y))
      if(same.heights)
        cal[,2] <- cal[,2]/max(cal[,2])
      cal[,2] <- ((height * cal[,2])/median(maxhght))/25 * (max(dlim) - min(dlim)) # scale the height of the blobs with the core length
      if(ncol(set$dets) > 4 && set$dets[i,5] == 0) # cal BP dates could have different relative height:
        cal[,2] <- calheight * cal[,2]

      x = cal[,1]
      y = cal[,2]
      y = y[!duplicated(x)]
      x = x[!duplicated(x)]
      cal <- approx(x, y, seq(min(x), max(x), length= 100)) # tmp but probably not a bad idea

      if(mirror)
        pol <- cbind(c(d-cal$y, d+rev(cal$y)), c(cal$x, rev(cal$x))) else
         if(up)
           pol <- cbind(d-c(0, cal$y, 0), c(min(cal$x), cal$x, max(cal$x))) else
             pol <- cbind(d+c(0, cal$y, 0), c(min(cal$x), cal$x, max(cal$x)))

      if(rotate.axes)
        pol <- cbind(pol[,2], pol[,1])
      if(ncol(set$dets)==4 && cc > 0 || (ncol(set$dets) > 4 && set$dets[i,5] > 0)) {
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
bacon.calib <- function(dat, set=get('info'), date.res=100, cutoff=0.01, postbomb=set$postbomb, normal=set$normal, t.a=set$t.a, t.b=set$t.b, delta.R=set$delta.R, delta.STD=set$delta.STD, ccdir="") {
  # read in the curves
  if(set$cc1=="IntCal20" || set$cc1=="\"IntCal20\"")
    cc1 <- read.table(paste0(ccdir, "3Col_intcal20.14C")) else
      cc1 <- read.csv(paste0(ccdir, set$cc1, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc2=="Marine20" || set$cc2=="\"Marine20\"")
    cc2 <- read.table(paste0(ccdir, "3Col_marine20.14C")) else
      cc2 <- read.csv(paste0(ccdir, set$cc2, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc3=="SHCal20" || set$cc3=="\"SHCal20\"")
    cc3 <- read.table(paste0(ccdir, "3Col_shcal20.14C")) else
      cc3 <- read.csv(paste0(ccdir, set$cc3, ".14C"), header=FALSE, skip=11)[,1:3]
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") cc4 <- NA else
    cc4 <- read.table(paste0(ccdir, set$cc4))[,1:3]

  if(postbomb != 0) {
    if(postbomb==1) bomb <- read.table(paste0(ccdir,"postbomb_NH1.14C"))[,1:3] else
      if(postbomb==2) bomb <- read.table(paste0(ccdir,"postbomb_NH2.14C"))[,1:3] else
        if(postbomb==3) bomb <- read.table(paste0(ccdir,"postbomb_NH3.14C"))[,1:3] else
          if(postbomb==4) bomb <- read.table(paste0(ccdir,"postbomb_SH1-2.14C"))[,1:3] else
            if(postbomb==5) bomb <- read.table(paste0(ccdir,"postbomb_SH3.14C"))[,1:3] else
              stop("cannot find postbomb curve #", postbomb, " (use values of 1 to 5 only)", call.=FALSE)
      # bomb.x <- seq(max(bomb[,1]), min(bomb[,1]), by=-.1) # interpolate
      bomb <- bomb[order(bomb[,1], decreasing=FALSE),]
      bomb.x <- seq(min(bomb[,1]), max(bomb[,1]), by=.1) # interpolate
      bomb.y <- approx(bomb[,1], bomb[,2], bomb.x)$y
      bomb.z <- approx(bomb[,1], bomb[,3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level=0)
      if(set$postbomb < 4)
        cc1 <- rbind(bomb, cc1, deparse.level=0) else
          cc3 <- rbind(bomb, cc3, deparse.level=0)
  }

  ## use Gaussian or t (Christen and Perez Radiocarbon 2009) calibration
  if(round(set$t.b-set$t.a) != 1)
    stop("t.b - t.a should always be 1, check the manual", call.=FALSE)

  d.cal <- function(cc, rcmean, w2, t.a, t.b) { # formula updated Oct 2020
    if(set$normal)
      cal <- cbind(cc[,1], dnorm(cc[,2], rcmean, sqrt(cc[,3]^2+w2))) else
        cal <- cbind(cc[,1], (t.b + ((rcmean-cc[,2])^2) / (2*(cc[,3]^2 + w2))) ^ (-1*(t.a+0.5))) # student-t
    cal[,2] <- cal[,2] / sum(cal[,2]) # normalise
    
    above <- which(cal[,2]/max(cal[,2]) > cutoff)
    if(length(above) > 5)
      return(cal[min(above):max(above),]) else
        return(cal)
  }

  # now calibrate all dates
  calib <- list(d=dat[,4], cc=dat[,4])
  if(ncol(dat)==4) { # only one type of dates (e.g., calBP, or all IntCal20 C14 dates)
    if(set$cc==0) {
      calib$cc <- rep(0, nrow(dat))
      xsteps <- min(dat[,3])/5 # minimum step size to cover the smallest error
      xseq1 <- seq(min(dat[,2])-(4*max(dat[,3])), max(dat[,2])+(4*max(dat[,3])), length=100*date.res)
      xseq2 <- seq(min(dat[,2])-(4*max(dat[,3])), max(dat[,2])+(4*max(dat[,3])), by=xsteps)
      xlength <- max(length(xseq1), length(xseq2)) # choose the longest vector
      x <- seq(min(dat[,2])-(4*max(dat[,3])), max(dat[,2])+(4*max(dat[,3])), length=xlength)
      ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
    } else {
        calib$cc <- rep(set$cc, nrow(dat))
        if(set$cc==1) ccurve <- cc1 else
          if(set$cc==2) ccurve <- cc2 else
            if(set$cc==3) ccurve <- cc3 else
              ccurve <- cc4
      }
    for(i in 1:nrow(dat))
      calib$probs[[i]] <- d.cal(ccurve, dat[i,2]-delta.R, dat[i,3]^2+delta.STD^2, set$t.a, set$t.b)
  } else {
      calib$cc <- dat[,5]
      for(i in 1:nrow(dat)) {
        dets <- c(NA, as.numeric(dat[i,-1])) # the first column is not numeric
        if(dets[5]==0) {
          x <- seq(dets[2]-(4*dets[3]), dets[2]+(4*dets[3]), length=date.res)
          ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
        } else {
            if(dets[5]==1) ccurve <- cc1 else if(dets[5]==2) ccurve <- cc2 else
              if(dets[5]==3) ccurve <- cc3 else ccurve <- cc4
            }
        delta.R <- set$delta.R; delta.STD <- set$delta.STD; t.a <- set$t.a; t.b <- set$t.b
        if(length(dets) >= 7 && dets[5] > 0) { # the user provided age offsets; only for C14 dates
          delta.R <- dets[6]
          delta.STD <- dets[7]
        }
        if(length(dets) >= 9) { # the user provided t.a and t.b values for each date
          t.a <- dets[8]
          t.b <- dets[9]
          if(round(t.b-t.a) != 1)
            stop("t.b - t.a should always be 1, check the manual", call.=FALSE)
        }
        calib$probs[[i]] <- d.cal(ccurve, dets[2]-delta.R, dets[3]^2+delta.STD^2, t.a, t.b)
      }
    }
  calib
}



#' @name draw.pbmeasured
#' @title Plot the 210Pb data
#' @description Produce a plot of the 210Pb data and their depths
#' @details This function is generally called internally to produce the age-depth graph.
#' It can be used to produce custom-built graphs.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
#' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param pb.lim Minimum and maximum of the 210Pb axis ranges, calculated automatically by default (\code{pb.lim=c()}).
#' @param age.lim Minimum and maximum of the age ranges to be used to plot 210Pb values. Calculated automatically by default (\code{age.lim=c()}).
#' @param d.lim Minimum and maximum depths to plot; calculated automatically by default (\code{d.lim=c()}).
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
#' @param pb.lab The label for the 210Pb axis (default \code{pb.lab="210Pb (Bq/kg)"} or \code{"210Pb (dpm/g)"}).
#' @param pbmeasured.col The label for the measured 210Pb data. \code{pbmeasured.col="blue"}.
#' @param pbmeasured.lty Line type of the measured 210Pb data. Defaults to continuous lines.
#' @param pb.log Use a log scale for the 210Pb-axis (default \code{pb.log=FALSE}).
#' @param supp.col Colour of the supported 210Pb data. Defaults to red: \code{supp.col="red"}.
#' @param newplot make new plot (default TRUE)
#' @param on.agescale Plot the Pb-210 on the cal BP scale. Defaults to FALSE.
#' @author Maarten Blaauw, J. Andres Christen, Marco Aquino-Lopez
#' @return A plot of the measured 210Pb values
#' @export
draw.pbmeasured <- function(set=get('info'), rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, BCAD=set$BCAD, pb.lim=c(), age.lim=c(), d.lim=c(), d.lab=c(), pb.lab=c(), pbmeasured.col="blue", pbmeasured.lty=1, pb.log=FALSE, supp.col="purple", newplot=TRUE, on.agescale=FALSE) {
  depths <- set$detsOrig[,2]
  dns <- set$detsOrig[,3]
  Pb <- set$detsOrig[,4]
  err <- set$detsOrig[,5]
  thickness <- set$detsOrig[,6]
  n <- nrow(set$detsOrig)

  if(length(pb.lim) == 0)
    pb.lim <- extendrange(c(0, Pb+2*err), f=c(0,0.05))

  # translate pb values to cal BP/AD values for plotting on the age axis
  pb2bp <- function(pb, pb.min=pb.lim[1], pb.max=pb.lim[2], agemin=min(age.lim), agemax=max(age.lim), AD=BCAD) {
    if(on.agescale) {
        if(AD) {
          ex <- (agemin - agemax) / (pb.max - pb.min)
          return(agemax + ex*pb)
        } else {
            ex <- (agemax - agemin) / (pb.max - pb.min)
            return(agemin + ex*pb)
        }
      } else
        return(pb)
  }

  if(newplot) {
    if(length(d.lab) == 0)
      d.lab <- paste0("depth (", set$depth.unit, ")")
    if(length(pb.lab) == 0)
      pb.lab <- ifelse(set$Bqkg, "210Pb (Bq/kg)", "210Pb (dpm/g)")

    if(length(d.lim) == 0)
      d.lim <- range(depths, set$supportedData[,3])
    if(rev.d)
      d.lim <- d.lim[2:1]
    if(rotate.axes)
      plot(0, type="n", ylim=d.lim, ylab=d.lab, xlim=pb2bp(pb.lim), xlab=pb.lab) else
        plot(0, type="n", xlim=d.lim, xlab=d.lab, ylim=pb2bp(pb.lim), ylab=pb.lab)
  }

  if(rotate.axes)
    rect(pb2bp(Pb-err), depths-thickness, pb2bp(Pb+err), depths, border=pbmeasured.col, lty=pbmeasured.lty) else
      rect(depths-thickness, pb2bp(Pb-err), depths, pb2bp(Pb+err), lty=pbmeasured.lty, border=pbmeasured.col)

  if(length(set$supportedData) > 0) {
    supp <- set$supportedData[,1]
    supperr <- set$supportedData[,2]
    suppd <- set$supportedData[,3]
    suppthick <- set$supportedData[,4]

    if(rotate.axes)
      rect(pb2bp(supp-supperr), suppd-suppthick, pb2bp(supp+supperr), suppd,
        border=supp.col, lty=pbmeasured.lty) else
        rect(suppd-suppthick, pb2bp(supp-supperr), suppd, pb2bp(supp+supperr),
          border=supp.col, lty=pbmeasured.lty)
  }
}



### for running Plum, but is looked for by generic agedepth() function, so is included in the rbacon code
#' @name draw.pbmodelled
#' @title Plot the 210Pb data
#' @description Produce a plot of the 210Pb data and their depths
#' @details This function is generally called internally to produce the age-depth graph.
#' It can be used to produce custom-built graphs.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}).
#' @param rev.age The direction of the calendar age axis can be reversed from the default (\code{rev.age=TRUE})
#' @param pb.lim Minimum and maximum of the 210Pb axis ranges, calculated automatically by default (\code{pb.lim=c()}).
#' @param d.lim Minimum and maximum depths to plot; calculated automatically by default (\code{d.lim=c()}).
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}.
#' @param pb.lab The label for the 210Pb axis (default \code{pb.lab="210Pb (Bq/kg)"} or \code{"210Pb (dpm/g)"}).
#' @param supp.col Colour of the supported 210Pb data. Defaults to red: \code{supp.col="red"}.
#' @param pbmodelled.col Colour of the modelled 210Pb values. Defaults to scales of blue: \code{pbmodelled.col=function(x) rgb(0,0,1,x)}.
#' @param pbmeasured.col Colour of the measured 210Pb values. Defaults to blue.
#' @param plot.measured Plot the measured 210Pb values (default \code{plot.measured=TRUE}).
#' @param age.lim values of the age axis. Used to calculate where to plot the pb values on the secondary axis
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.7, .7, .0)}.
#' @param pb.lty Line type of measured Pb-210 data.
#' @author Maarten Blaauw, J. Andres Christen, Marco Aquino-Lopez
#' @return A plot of the modelled (and optionally the measured) 210Pb values
#' @export
draw.pbmodelled <- function(set=get('info'), BCAD=set$BCAD, rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, pb.lim=c(), d.lim=c(), d.lab=c(), pb.lab=c(), pbmodelled.col=function(x) rgb(0,0,1,x), pbmeasured.col="blue", supp.col="purple", plot.measured=TRUE, age.lim=c(), mgp=mgp, pb.lty=1) {
  pb <- set$dets[set$dets[,9] == 5,]
  depths <- pb[,4] # set$detsOrig[,2]
  dns <- pb[,6] # set$detsOrig[,3]
  Pb <- pb[,2] # set$detsOrig[,4]
  err <- pb[,3] # set$detsOrig[,5]
  thickness <- pb[,5] # set$detsOrig[,6]
  n <- nrow(pb)

  if(ncol(pb) > 6) {
    supp <- pb[,7] # set$detsOrig[,7]
    supperr <- pb[,8] # set$detsOrig[,8]
  } else {
    supp <- set$supportedData[,1]
    supperr <- set$supportedData[,2]
    suppd <- set$supportedData[,3]
    suppthick <- set$supportedData[,4]
  }

  if(length(d.lab) == 0)
    d.lab <- paste("depth (", set$depth.unit, ")", sep="")
  if(length(pb.lab) == 0)
    pb.lab <- ifelse(set$Bqkg,
      expression(""^210*"Pb (Bq/kg)"),
        expression(""^210*"Pb (dpm/g)"))

  if(length(d.lim) == 0)
    d.lim <- range(depths)
  if(rev.d)
    d.lim <- d.lim[2:1]

  if(length(set$phi) > 0) {
    Ai <- list(x=NULL, y=NULL)
    hght <- 0; pbmin <- c(); pbmax <- 0
    A.rng <- array(0, dim=c(n,2))
    for(i in 1:length(depths)) {
      A <- A.modelled(depths[i]-thickness[i], depths[i], dns[i], set)
      tmp <- density(A)
      Ai$x[[i]] <- tmp$x
      Ai$y[[i]] <- tmp$y
      hght <- max(hght, Ai$y[[i]])
      pbmin <- min(pbmin, Ai$x[[i]])
      pbmax <- max(pbmax, Ai$x[[i]])
      A.rng[i,] <- quantile(A, c((1-set$prob)/2, 1-(1-set$prob)/2))
    }

    if(length(pb.lim) == 0)
      pb.lim <- extendrange(c(0, Pb-2*err, Pb+2*err, pbmax), f=c(0, 0.05))

    pb2bp <- function(pb, pb.min=pb.lim[1], pb.max=pb.lim[2], agemin=min(age.lim), agemax=max(age.lim), AD=BCAD) {
      if(AD) {
        ex <- (agemin - agemax) / (pb.max - pb.min)
        return(agemax + ex*pb)
      } else {
          ex <- (agemax - agemin) / (pb.max - pb.min)
          return(agemin + ex*pb)
      }
    }

    # save the values for later
    set$Ai <- Ai
    set$A.rng <- A.rng
    assign_to_global("info", set, .GlobalEnv) # doesn't work

    this <- ifelse(rotate.axes, 3, 4)
    pretty.pb <- pretty(c(pbmin, pbmax)) # not OK?
    pretty.pb <- pretty(pb.lim)
    onbp <- pb2bp(pretty.pb)
    if(BCAD)
      axis(this, rev(onbp), rev(pretty.pb), col=pbmeasured.col, col.axis=pbmeasured.col, col.lab=pbmeasured.col) else
        axis(this, onbp, pretty.pb, col=pbmeasured.col, col.axis=pbmeasured.col, col.lab=pbmeasured.col)

    mtext(pb.lab, this, 2.5, col=pbmeasured.col, cex=.8)

    for(i in 1:length(depths)) {
      if(BCAD) {
        ages <- pb2bp(rev(Ai$x[[i]]))
        z <- t(rev(Ai$y[[i]]))/hght
      } else {
          ages <- pb2bp(Ai$x[[i]])
          z <- t(Ai$y[[i]])/hght
        }

      if(rotate.axes)
        image(ages, c(depths[i]-thickness[i], depths[i]), z, col=pbmodelled.col(seq(0, 1-max(z),  length=50)), add=TRUE) else
          image(c(depths[i]-thickness[i], depths[i]), ages, z, col=pbmodelled.col(seq(0, 1-max(z), length=50)), add=TRUE)
    }
  }

#   # indicate in redscale which Pb-210 data have most likely reached background
#   if(draw.background) {
#     bg <- background(set)
#     set$background <- bg
#     assign_to_global("info", set, .GlobalEnv)
#     if(rotate.axes)
#       abline(h=set$dets[,4]-(set$dets[,5]/2), col=rgb(bg,0,0,bg), lty=3, lwd=bg) else
#         abline(v=set$dets[,4]-(set$dets[,5]/2), col=rgb(bg,0,0,bg), lty=3, lwd=bg)
#   }

    if(plot.measured)
      draw.pbmeasured(newplot=FALSE, rotate.axes=rotate.axes, BCAD=BCAD, on.agescale=TRUE, pb.lim=pb.lim, age.lim=age.lim, supp.col=supp.col)
}


### for running Plum, but is looked for by generic agedepth() function (through draw.pbmodelled()), so is included in the rbacon code
#' @name A.modelled
#' @title Calculate modelled 210Pb
#' @description Calculate modelled 210Pb values of a sample slice, based on the parameters of the age-model (i.e., time passed since deposition of the bottom and top of the slice), supported and influx
#' @param d.top top depth of the slice
#' @param d.bottom bottom depth of the slice
#' @param dens Density of the slice (in g/cm3)
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param phi The modelled values of the 210Pb influx
#' @param sup The modelled values of the supported 210Pb
#' @author Maarten Blaauw
#' @return a list of modelled values of A
#' @export
A.modelled <- function(d.top, d.bottom, dens, set=get('info'), phi=set$phi, sup=set$ps) {
  if(d.top >= d.bottom)
    stop("\n d.top should be above d.bottom", call.=FALSE)
  dd <- 1
  if(ncol(cbind(sup)) > 1) { # then multiple, varying estimates of supported, find the one belonging to the specified depth interval
    dd <- set$supportedData[,3] # bottom depths
    dd <- max(1, which(dd <= d.bottom))
    sup <- sup[,dd]
  }

  t.top <- Bacon.Age.d(d.top, BCAD=F) - set$theta0
  t.bottom <- Bacon.Age.d(d.bottom, BCAD=F) - set$theta0
  #  multiply <- ifelse(set$Bqkg, 10, 500)
  multiply <- 1 # since Bqkg or dpmg is already set earlier (for set$dets and set$detsPlum)
  return(sup + ((phi / (.03114*multiply*dens) ) * (exp( -.03114*t.top) - exp(-.03114*t.bottom)) ) )
} 



### for running Plum, but is looked for by generic agedepth() function (through draw.pbmodelled()), so is included in the rbacon code
#' @name background
#' @title calculate probabilities that Pb-210 data have reached background levels
#' @description Checks which of the Pb-210 data most likely have reached background levels and thus are below the detection limit Al (probabilities between 0 and 1)
#' @author Maarten Blaauw
#' @return a list of probabilities for each Pb-210 data point
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param Al The detection limit. Default \code{Al=0.1}.
#' @export
background <- function(set=get('info'), Al=set$Al) {
  if(set$isplum) { # works with Pb-210 data only
    pb <- 0
    its <- nrow(set$output)
#    dets <- set$detsOrig[,c(2,6,3)] # we need maxdepth, mindepth, density
    dets <- set$dets[which(set$dets[,9] == 5),4:6] # should leave out any non-Pb data
    ps <- cbind(set$ps)
    for(i in 1:nrow(dets)) {
      As <- A.modelled(dets[i,1]-dets[i,2], dets[i,1], dets[i,3])
      if(set$ra.case == 2)
        ps <- set$ps[,i] else
          ps <- set$ps
      bg <- which((As - ps) <= Al) # which modelled data are at or below the detection limit?
      pb[i] <- length(bg) / its
    }
    return(pb)
  }
}
