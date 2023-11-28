


#' @name add.dates
#' @title Add dates to age-depth plots
#' @description Add dated depths to plots, e.g. to show dates that weren't used in the age-depth model
#' @details Sometimes it is useful to add additional dating information to age-depth plots, e.g., to show outliers or how dates calibrate with different estimated offsets. Calls rintcal's draw.dates function.
#' @param mn Reported mean of the date. Can be multiple dates. Negative numbers indicate postbomb dates (if cc > 0).
#' @param sdev Reported error of the date. Can be multiple dates.
#' @param depth Depth of the date.
#' @param cc The calibration curve to use: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine), \code{cc=0} for none (dates that are already on the cal BP scale).
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param above Threshold for plotting of probability values. Defaults to \code{above=1e-3}.
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}
#' @param normal By default, Bacon uses the t-distribution (Christen and Perez 2009) to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param delta.R Mean of core-wide age offsets (e.g., regional marine offsets).
#' @param delta.STD Error of core-wide age offsets (e.g., regional marine offsets).
#' @param t.a The dates are treated using the t distribution by default (\code{normal=FALSE}).
#' The t model has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param t.b The dates are treated using the t distribution by default (\code{normal=FALSE}).
#' The t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
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
#' @param cc.dir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{cc.dir=c()}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A date's distribution, added to an age-depth plot.
#' @examples
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth()
#'   add.dates(5000, 100, 60)
#' }
#' @export
add.dates <- function(mn, sdev, depth, cc=1, set=get('info'), above=1e-6, postbomb=0, normal=TRUE, delta.R=set$delta.R, delta.STD=set$delta.STD, t.a=set$t.a, t.b=set$t.b, date.res=100, height=.1, calheight=1, agesteps=1, cutoff=0.005, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5), rotate.axes=FALSE, mirror=TRUE, up=TRUE, BCAD=FALSE, pch=4, cc.dir=c()) {

  dists <- draw.dates(mn-delta.R, sqrt(sdev^2+delta.STD^2), depth, cc=cc, postbomb=postbomb, normal=normal, t.a=t.a, t.b=t.b, dist.res=date.res, ex=height, threshold=cutoff, col=col, border=border, draw.hpd=FALSE, rotate.axes=!rotate.axes, mirror=mirror, up=up, cc.dir=cc.dir, add=TRUE, BCAD=BCAD)
  
  if(length(pch) > 0) {
    best <- c()
    for(i in 1:length(mn)) 
      best[i] <- dists[[1]][,i][which(dists[[2]][,i] == max(dists[[2]][,i]))][1] 
    if(rotate.axes)
      points(best, depth, pch=pch, col=col) else
        points(depth, best, pch=pch, col=col)
  }  
}



#' @name calib.plot
#' @title Plot the dates
#' @description Produce a plot of the dated depths and their dates
#' @details This function is generally called internally to produce the age-depth graph.
#' It can be used to produce custom-built graphs.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param dets The set of determinations to be plotted.
#' @param accordion If depths have to be squeezed/stretched, the parameters can be set here. Defaults to being empty, but requires 2 parameters if active, e.g., \code{accordion=c(10,20)}.
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
calib.plot <- function(set=get('info'), dets=set$dets, accordion=c(), BCAD=set$BCAD, cc=set$cc, rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, rev.yr=rev.age, age.lim=c(), yr.lim=age.lim, date.res=100, d.lab=c(), age.lab=c(), yr.lab=age.lab, height=1, calheight=1, mirror=TRUE, up=TRUE, cutoff=.1, C14.col=rgb(0,0,1,.5), C14.border=rgb(0,0,1,.75), cal.col=rgb(0,.5,.5,.5), cal.border=rgb(0,.5,.5,.75), dates.col=c(), slump.col=grey(0.8), new.plot=TRUE, plot.dists=TRUE, same.heights=FALSE) {
	
  # agedepth calls as follows:
  #calib.plot(set, BCAD=BCAD, cc=cc, rotate.axes=rotate.axes, height=height, calheight=calheight, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, dates.col=dates.col, new.plot=FALSE, same.heights=same.heights)
  
  # we have to set ka (kcal?) as an option as well

  #dets <- set$dets # should this be set$detsBacon (if it exists)?
  d.R <- 0; d.STD <- 0
  t.a <- set$t.a; t.b <- set$t.b
  cc <- rep(cc, nrow(dets))
  if(ncol(dets) > 4) {
    cc <- dets[,5]
    if(ncol(dets) >= 7) {
      d.R <- dets[,6]
      d.STD <- dets[,7] 
    }
    if(ncol(dets) > 7) {
      t.a <- dets[,8]
      t.b <- dets[,9]
    }
  }

  #draw.dates(dets[,2]-d.R, sqrt(dets[,3]^2+d.STD^2), dets[,4], cc=cc, BCAD=BCAD, rotate.axes=!rotate.axes, d.rev=rev.d, dist.res=date.res, d.lab=d.lab, age.rev=age.rev, age.lim=age.lim, height=ex, mirror=mirror, up=up, threshold=cutoff, col=C14.col, border=C14.border, cal.col=cal.col, cal.border=cal.border, t.a=t.a, t.b=t.b)
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
    d.lab <- paste0("depth (", set$depth.unit, ")")
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
      if(length(accordion) == 2)
        d <- stretch(d, accordion[1], accordion[2])
      if(BCAD)
        cal[,1] <- 1950-cal[,1]
      if((max(cal[,1]) - min(cal[,1])) > 4*agesteps)
        cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), by=agesteps)) else
          cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), length=100))
      # the above is not ideal because it causes different heights for very precise distributions

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
     # if(ncol(set$dets)==4 && cc > 0 || (ncol(set$dets) > 4 && set$dets[i,5] > 0)) {
      if(cc[i] == 0) {
          col <- cal.col
          border <- cal.border
        } else {
            col <- C14.col
            border <- C14.border
          }

     if(length(dates.col) > 0) {
        col <- dates.col[i]
        border <- dates.col[i]
      }
      polygon(pol, col=col, border=border)
    }
}


# bacon.calib is used by the Bacon function: info$calib <- bacon.calib(dets, info, date.res, cc.dir=cc.dir, cutoff=cutoff)
# it calibrates C14 dates and calculate distributions for any calendar dates
# it then returns d, cc, and probs
bacon.calib <- function(dat, set=get('info'), date.res=100, cutoff=0.01, postbomb=set$postbomb, normal=set$normal, t.a=set$t.a, t.b=set$t.b, delta.R=set$delta.R, delta.STD=set$delta.STD, cc.dir=c()) {
  # read in the curves

  cc1 <- ccurve(set$cc1, cc.dir=cc.dir)
  cc2 <- ccurve(set$cc2, cc.dir=cc.dir)
  cc3 <- ccurve(set$cc3, cc.dir=cc.dir)
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") cc4 <- NA else
     cc4 <- fastread(file.path(cc.dir, set$cc4))[,1:3] # file.path was paste0

  if(postbomb != 0) {
    bomb <- ccurve(postbomb, postbomb=TRUE, glue=FALSE, cc.dir=cc.dir) # glue=FALSE added July 2023
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
        cal <- cbind(cc[,1], (t.b + ((rcmean-cc[,2])^2) / (2*(cc[,3]^2 + w2))) ^ (-1*(t.a+0.5))) # t dist
    cal[,2] <- cal[,2] / sum(cal[,2]) # normalise
    
    above <- which(cal[,2]/max(cal[,2]) > cutoff)
    if(length(above) > 5)
      return(cal[min(above):max(above),]) else
        return(cal)
  }

  calib <- list(d=dat[,4], cc=set$cc)
  if(ncol(dat) == 4) # only one type of dates (e.g., calBP, or all IntCal20 C14 dates)
    dat[,5] <- rep(set$cc, nrow(dat))
  calib$cc <- dat[,5]

  for(i in 1:nrow(dat)) {
    dets <- c(NA, as.numeric(dat[i,-1])) # the first column is not numeric
    if(dets[5]==0) {
      x <- seq(dets[2]-(4*dets[3]), dets[2]+(4*dets[3]), length=date.res)
      calcurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
    } else {
        if(dets[5]==1) calcurve <- cc1 else
          if(dets[5]==2) calcurve <- cc2 else
            if(dets[5]==3) calcurve <- cc3 else
              calcurve <- cc4
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
    calib$probs[[i]] <- d.cal(calcurve, dets[2]-delta.R, dets[3]^2+delta.STD^2, t.a, t.b)
  }

  return(calib)
}
