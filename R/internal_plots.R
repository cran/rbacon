

#################### user-invisible plot functions ####################

# to plot greyscale/ghost graphs of the age-depth model
agedepth.ghost <- function(set=get('info'), dseq=c(), d.min=set$d.min, d.max=set$d.max, accordion=c(), BCAD=set$BCAD, rotate.axes=FALSE, d.res=400, age.res=400, rgb.res=100, dark=c(), rgb.scale=c(0,0,0), cutoff=0.001, age.lim) {
  if(length(dseq) == 0)
    dseq <- seq(d.min, d.max, length=d.res)
  if(set$isplum) # plum has a strange feature with a grey shape appearing
    dseq <- dseq[-1] # at dmin. Thus removing the first depth
  if(length(set$slump) > 0) {
    d.inside <- c()
    for(i in 1:nrow(set$slump)) {
      inside <- which(dseq < max(set$slump[i,]))
      inside <- which(dseq[inside] > min(set$slump[i,]))
      d.inside <- c(d.inside, inside)
    }
  dseq <- dseq[-d.inside]
  }

 if(length(accordion) == 2)
    dseq <- squeeze(dseq, accordion[1], accordion[2])

  Bacon.hist(dseq, set, BCAD=BCAD, calc.range=FALSE, draw=FALSE)
  hists <- get('hists')
  scales <- array(0, dim=c(length(dseq), age.res))
  ageseq <- seq(min(age.lim), max(age.lim), length=age.res)
  for(i in 1:length(hists)) { # was length(dseq)
  if(length(hists[[i]]) < 7)
      ages <- sort(unlist(hists[[i]])) else {
       if(length(hists[[i]]$th0) == 0) # can't calculate ages beyond upper depths
         ages <- NA else
           ages <- seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n)
       if(length(ages[!is.na(ages)]) > 0)
        scales[i,] <- approx(ages, hists[[i]]$counts, ageseq, rule=2)$y
     }
  }
  minmax <- hists[[length(hists)]]$min
  maxmax <- hists[[length(hists)]]$max
  scales <- scales/maxmax # normalise to the height of most precise age estimate
  if(length(dark) == 0)
    dark <- 10 * minmax/maxmax
  scales[scales > dark] <- dark
  scales <- scales/max(scales) # May 2021
  dseq <- sort(dseq)
  if(length(accordion) == 2)
    dseq <- stretch(dseq, accordion[1], accordion[2]) # careful now!
  cols <- rgb(rgb.scale[1], rgb.scale[2], rgb.scale[3], seq(0,1, length=rgb.res))
  
  scales[scales<cutoff] <- NA # so that pixels with probs very close to 0 are not plotted as white but empty
  
  if(rotate.axes)
    image(ageseq, dseq, t(scales), add=TRUE, col=cols, useRaster=FALSE) else
      image(dseq, ageseq, scales, add=TRUE, col=cols, useRaster=FALSE)
}



# Time series of the log of the posterior
PlotLogPost <- function(set, from=0, to=set$Tr, xaxs="i", yaxs="i", panel.size=.9)
  plot(from:(to-1), -set$Us[(from+1):to], type="l",
    ylab="Log of Objective", xlab="Iteration", main="", xaxs=xaxs, yaxs=yaxs, col=grey(0.4), cex.axis=panel.size)



# plot the prior for the accumulation rate
PlotAccPrior <- function(s, mn, set=get('info'), depth.unit=depth.unit, age.unit=age.unit, main="", xlim=c(0, 3*max(mn)), xlab=c(), ylab="Density", add=FALSE, legend=TRUE, csize=.9) {
  o <- order(s, decreasing=TRUE)
#  priors <- unique(cbind(s[o],mn[o])[,1:2])
  priors <- cbind(unique(s[o]), unique(mn[o]))[,1:2]
  x <- 0
  if(length(xlab) == 0)
    xlab <- paste0("Acc. rate (", noquote(age.unit), "/", noquote(depth.unit), ")")

  if(length(priors) == 2) {
    curve(dgamma(x, s, s/mn), col=3, lwd=2, from=0, to=max(xlim), xlim=xlim, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("acc.shape: ", priors[1], "\nacc.mean: ", priors[2])
  } else {
      priors <- priors[order(priors[,1]*priors[,2]),]
      curve(dgamma(x, priors[1,1], priors[1,1]/priors[1,2]), col=3, lwd=2, from=0, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
      for(i in 2:nrow(priors))
        curve(dgamma(x, priors[i,1], priors[i,1]/priors[i,2]), col=3, lwd=2, from=0, xlim=xlim, xlab=xlab, ylab=ylab, add=if(i==1) add else TRUE)
      txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
    }
  if(legend)
    legend("topleft", txt, bty="n", cex=csize, text.col=2, xjust=1)
}



# plot the prior for the memory (= accumulation rate varibility between neighbouring depths)
PlotMemPrior <- function(s, mn, thick, ds=1, set=get('info'), xlab="Memory (ratio)", ylab="Density", main="", add=FALSE, legend=TRUE, csize=.9) {
  o <- order(s, decreasing=TRUE)
#  priors <- unique(cbind(s[o],mn[o])[,1:2])
  priors <- cbind(unique(s[o]), unique(mn[o]))[,1:2]
  x <- 0

  if(length(priors)==2) {
    curve(dbeta(x, s*mn, s*(1-mn)), from=0, to=1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
    txt <- paste0("mem.strength: ", s, "\nmem.mean: ", mn, "\n", set$K, " ", round(thick,3)," ", noquote(set$depth.unit), " sections")
  } else {
	  priors <- priors[order(priors[,1]*priors[,2]),]
      curve(dbeta(x, priors[1,1]*priors[1,2], priors[1,1]*(1-priors[1,2])), from=0, to=1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
      for(i in 2:nrow(priors))
        curve(dbeta(x, priors[i,1]*priors[i,2], priors[i,1]*(1-priors[i,2])), from=0, to=1, col=3, lwd=2, xlab="", ylab="", add=TRUE)
      txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
    }
  if(legend)
    legend("topleft", txt, bty="n", cex=csize, text.col=2, xjust=0)
  warn <- FALSE
  for(i in s)
    for(j in mn)
      if(i*(1-j) <= 1) warn <- 1
  if(warn)
    message("Warning! Chosen memory prior might cause problems.\nmem.strength * (1 - mem.mean) should be smaller than 1 ")
}



# plot the prior for the hiatus length
PlotHiatusPrior <- function(mx=set$hiatus.max, hiatus=set$hiatus.depths, set=get('info'), xlab=paste0("Hiatus size (", set$age.unit, ")"), ylab="Density", main="", xlim=c(0, 1.1*max(mx)), add=FALSE, legend=TRUE, csize=.9) {
  if(add)
    lines(c(0, 0, mx, mx), c(0, 1/mx, 1/mx, 0), col=3, lwd=2) else
      plot(c(0, 0, mx, mx), c(0, 1/mx, 1/mx, 0), xlab=xlab, ylab=ylab, xlim=xlim, type="l", col=3, lwd=2)

  txt <- paste("hiatus.max: ", toString(mx))
  if(legend)
    legend("topleft", txt, bty="n", cex=csize, text.col=2, xjust=0)
}



# plot the Supported prior (for plum)
PlotSuppPrior <- function(set=get('info'), xaxs="i", yaxs="i", legend=TRUE, csize=.9) {
  if(set$Bqkg)
    lab = "s.Bq/Kg" else 
      lab = "dpm/g"

  x <- 0
  s = set$s.shape
  mn = set$s.mean
  xlim = c(0, 3*mn)
  curve(dgamma(x, s, s/mn), col=3, lwd=2, from=0, to=max(xlim), xlim=xlim, xlab=lab, ylab="")
  txt <- paste("Supported", "\nS.shape: ", s, "\nS.mean: ", mn)
  if(legend)
    legend("topleft", txt, bty="n", cex=csize, text.col=2, xjust=0)
}



# to plot the prior for Supply (for plum)
PlotPhiPrior <- function(s, mn, set=get('info'), depth.unit=depth.unit, age.unit=age.unit, main="", xlim=c(0, 3*max(mn)), xlab="Bq/m^2 yr", ylab="", add=FALSE, legend=TRUE, csize=.9) {
  x <- 0
  s = set$phi.shape
  mn = set$phi.mean
  xlim = c(0, 3*mn)
  curve(dgamma(x, s, s/mn), col=3, lwd=2, from=0, to=max(xlim), xlim=xlim, xlab=xlab, ylab=ylab)

  txt <- paste( "Influx", "\nAl: ", toString(round(set$Al,2)), "\nPhi.shape: ", toString(round(s,2)), "\nPhi.mean: ", toString(round(mn,2)) )
  if(legend)
    legend("topleft", txt, bty="n", cex=csize, text.col=2, xjust=0)
}



# plot the posterior (and prior) of the accumulation rate
PlotAccPost <- function(set=get('info'), s=set$acc.shape, mn=set$acc.mean, main="", depth.unit=set$depth.unit, age.unit=set$age.unit, ylab="Frequency", xaxs="i", yaxs="i", yaxt="n", prior.size=.9, panel.size=.9, acc.xlim=c(), acc.ylim=c()) {
  hi <- 2:(set$K-1)
  if(!is.na(set$hiatus.depths)[1])
    for(i in set$hiatus.depths)
      hi <- hi[-max(which(set$elbows < i))]
  post <- c()
  for(i in hi)
    post <- c(post, set$output[[i]])
  post <- density(post, from=0)
  post <- cbind(c(0, post$x, max(post$x)), c(0, post$y, 0))
  maxprior <- dgamma((s-1)/(s/mn), s, s/mn)
  if(is.infinite(max(maxprior)))
    max.y <- max(post[,2]) else
      max.y <- max(maxprior, post[,2])
  if(length(acc.xlim) == 0)
    acc.xlim <- range(0, post[,1], 2*mn)
  if(length(acc.ylim) == 0)
    acc.ylim <- c(0, 1.05*max.y)
  acc.lab <- paste0("Acc. rate (", age.unit, "/", depth.unit, ")")
  plot(0, type="n", xlim=acc.xlim, xlab=acc.lab, ylim=acc.ylim, ylab="", xaxs=xaxs, yaxs=yaxs, yaxt=yaxt, cex.axis=panel.size)
  polygon(post, col=grey(.8), border=grey(.4))
  PlotAccPrior(s, mn, add=TRUE, xlim=acc.xlim, xlab="", ylab=ylab, main=main, csize=prior.size)
}



# plot the posterior (and prior) of the memory
PlotMemPost <- function(set=get('info'), corenam, K, main="", s=set$mem.strength, mn=set$mem.mean, xlab=paste("Memory"), ylab="Density", ds=1, thick, xaxs="i", yaxs="i", yaxt="n", prior.size=.9, panel.size=.9, mem.xlim=c(), mem.ylim=c()) {
  post <- density(set$output[,2+set$K]^(1/set$thick), from=0, to=1) # was output[,set$n] but not ok for rplum
  post <- cbind(c(min(post$x), post$x, max(post$x)), c(0, post$y, 0))
  maxprior <- max(dbeta((0:100)/100, s*mn, s*(1-mn)))
  if(length(mem.xlim) == 0)
    mem.xlim <- range(post[,1])
  if(is.infinite(max(maxprior)))
    max.y <- max(post[,2]) else
      max.y <- max(maxprior, max(post[,2]))
  if(length(mem.ylim) == 0)
    mem.ylim <- c(0, 1.05*max.y)
  plot(0, type="n", xlab=xlab, xlim=mem.xlim, ylim=mem.ylim, ylab="", main="", xaxs=xaxs, yaxs=yaxs, yaxt=yaxt, cex.axis=panel.size)
  polygon(post, col=grey(.8), border=grey(.4))
  PlotMemPrior(s, mn, thick, add=TRUE, xlab="", ylab=ylab, main=main, csize=prior.size)
}



# plot the posterior (and prior) of the hiatus
PlotHiatusPost <- function(set=get('info'), mx=set$hiatus.max, main="", xlab=paste0("Hiatus size (", set$age.unit, ")"), ylab="Frequency", after=set$after, xaxs="i", yaxs="i", yaxt="n", prior.size=.9, panel.size=.9, hiatus.xlim=c(), hiatus.ylim=c()) {
  gaps <- c()
  for(i in set$hiatus.depths) {
    below <- Bacon.Age.d(i+after, set)
    above <- Bacon.Age.d(i-after, set)
    gaps <- c(gaps, below - above)
  }
  if(length(hiatus.xlim) == 0)
    hiatus.xlim <- c(0, 1.1*(max(mx, gaps)))
  max.y <- 1.1/mx
  if(length(gaps) > 1) {
    gaps <- density(gaps, from=0)
    max.y <- max(max.y, gaps$y)
  }
  if(length(hiatus.ylim) == 0)
    hiatus.ylim <- c(0, max.y)
  plot(0, type="n", main="", xlab=xlab, xlim=hiatus.xlim, ylab=ylab, ylim=hiatus.ylim, xaxs=xaxs, yaxs=yaxs, yaxt=yaxt, cex.axis=panel.size)
  if(length(gaps) > 1)
    polygon(cbind(c(min(gaps$x), gaps$x, max(gaps$x)), c(0,gaps$y,0)),
    col=grey(.8), border=grey(.4))
  PlotHiatusPrior(add=TRUE, xlab="", ylab=ylab, main=main, csize=prior.size)
}



# plot the Supported data (for plum)
PlotSuppPost <- function(set=get('info'), xaxs="i", yaxs="i", legend=TRUE, supp.xlim=c(), supp.ylim=c(), yaxt="n", prior.size=.9, panel.size=.9) {
  lab <- ifelse(set$Bqkg, "Bq/kg", "dpm/g")

  if(length(supp.xlim) == 0)
    supp.xlim <- c(min( set$detsPlum[,4]), max(set$detsPlum[,4]))

  if(set$nPs > 1) {
    rng <- array(NA, dim=c(set$nPs, 22)) # 22 is the number of segments to draw. Always???
    for(i in 1:set$nPs) {
      rng[i,1:21] <- quantile( set$ps[,i] , seq(0,2,0.1)/2)
      rng[i,22] <- mean( set$ps[,i] )
    }

    if(length(supp.ylim) == 0)
      supp.ylim <- c(min( rng[,1]), max(rng[,21]))

    plot(0, type="n", ylim=supp.ylim, xlim=supp.xlim, main="", xlab="Depth (cm)", ylab=lab, yaxt=yaxt, cex.axis=panel.size)
    n = 21
    colorby = 1.0 / (n/2)
    for(i in 1:(n/2)) {
      segments(set$detsPlum[,4], rng[,i], set$detsPlum[,4], rng[,(i+1)], grey(1.0-colorby*i), lwd=3)
      segments(set$detsPlum[,4], rng[,n-i], set$detsPlum[,4], rng[,n-(i-1)], grey(1.0-colorby*i), lwd=3)
    }
    lines(set$detsPlum[,4], rng[,22], col="red", lty=12) # mean

  } else {
    post <- density(set$ps)
    suppdata <- set$supportedData[,1:2]
    rng <- range(post$x, suppdata[,1]-suppdata[,2], suppdata[,1]+suppdata[,2])

    if(length(supp.ylim) == 0)
      supp.ylim <- c(0, max(post$y))

    plot(post, type="n", xlab=lab, xlim=rng, main="", ylab="", yaxt=yaxt, cex.axis=panel.size)
    polygon(post, col=grey(.8), border=grey(.4))
    lines(seq(min(set$ps),max(set$ps),.05), dgamma(seq(min(set$ps), max(set$ps), .05), shape=set$s.shape, scale=set$s.mean/set$s.shape), col=3, lwd=2)

    # plot the measurements of supported
    # extend the plot's horizontal axis, range(...)
    y <- seq(0, max(post$y), length=nrow(suppdata)+50)[1+(1:nrow(suppdata))] # at bottom graph
    segments(suppdata[,1]-suppdata[,2], y, suppdata[,1]+suppdata[,2], y, col=rgb(.5,0,.5,.5))
    points(suppdata[,1], y, pch=20, col=rgb(.5,0,.5,.5))
  }
  txt <- paste0("supported", "\ns.shape: ", set$s.shape, "\ns.mean: ", set$s.mean)  

  if(legend)
    legend("topright", txt, bty="n", cex=prior.size, text.col=2, adj=c(0,.2))
}



# plot the Supply data (for plum)
PlotPhiPost <- function(set=get('info'), xlab=paste0("Bq/",expression(m^2)," yr"), ylab="", xaxs="i", yaxs="i", legend=TRUE, yaxt="n", csize=.8, prior.size=.9, panel.size=.9, phi.xlim=c(), phi.ylim=c()) {
  post <- density(set$phi)
  if(length(phi.xlim) == 0)
    phi.xlim <- range(post$x)
  if(length(phi.ylim) == 0)
    phi.ylim <- range(post$y)

  plot(post, type="n", xlim=phi.xlim, ylim=phi.ylim, xlab=xlab, ylab=ylab, main="", yaxt=yaxt, cex.axis=panel.size)
  polygon(post, col=grey(.8), border=grey(.4))

  x <- seq(phi.xlim[1], phi.xlim[2], length=50)
  lines(x, dgamma(x,shape=set$phi.shape,scale=set$phi.mean/set$phi.shape), col=3, lwd=2)

  txt <- paste( "influx", "\nAl: ", toString(round(set$Al,2)), "\nphi.shape: ", toString(round(set$phi.shape,2)), "\nphi.mean: ", toString(round(set$phi.mean,2)) )

  if(legend)
    legend("topright", txt, bty="n", cex=prior.size, text.col=2, adj=c(0,.2))
}
