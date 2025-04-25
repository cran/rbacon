### functions which are for running Plum, but are looked for by generic agedepth() function, so are included in the rbacon code

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
draw.pbmeasured <- function(set=get('info'), rotate.axes=FALSE, rev.d=FALSE, rev.age=FALSE, BCAD=set$BCAD, pb.lim=c(), age.lim=c(), d.lim=c(), d.lab=c(), pb.lab=c(), pbmeasured.col="blue", pbmeasured.lty=2, pb.log=FALSE, supp.col="purple", newplot=TRUE, on.agescale=FALSE) {
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
    d.lab <- paste0("depth (", set$depth.unit, ")")
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
    #assign_to_global("info", set, .GlobalEnv) # doesn't work

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
        ghost.mirror(ages, c(depths[i]-thickness[i], depths[i]), t(z), col=pbmodelled.col(seq(0, 1-max(z),  length=50))) else
          ghost.mirror(c(depths[i]-thickness[i], depths[i]), ages, z, col=pbmodelled.col(seq(0, 1-max(z), length=50)))
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
      draw.pbmeasured(set=set, newplot=FALSE, rotate.axes=rotate.axes, BCAD=BCAD, on.agescale=TRUE, pb.lim=pb.lim, age.lim=age.lim, supp.col=supp.col)	
}



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

  t.top <- Bacon.Age.d(d.top, set=set, BCAD=FALSE) - set$theta0
  t.bottom <- Bacon.Age.d(d.bottom, set=set, BCAD=FALSE) - set$theta0
  #  multiply <- ifelse(set$Bqkg, 10, 500)
  multiply <- 1 # since Bqkg or dpmg is already set earlier (for set$dets and set$detsPlum)
  return(sup + ((phi / (.03114*multiply*dens) ) * (exp( -.03114*t.top) - exp(-.03114*t.bottom)) ) )
}



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
    dets <- set$dets[which(set$dets[,9] == 5),4:6] # only Pb data
    ps <- cbind(set$ps)
    for(i in 1:nrow(dets)) {
      As <- A.modelled(dets[i,1]-dets[i,2], dets[i,1], dets[i,3], set=set)
      if(set$ra.case == 2)
        ps <- set$ps[,i] else
          ps <- set$ps
      bg <- which((As - ps) <= Al) # which modelled data are at or below the detection limit?
      pb[i] <- length(bg) / its
    }
    return(pb)
  }
}
