

#################### user-invisible plot functions ####################

# to plot greyscale/ghost graphs of the age-depth model
.agedepth.ghost <- function(set=get('info'), d.min=set$d.min, d.max=set$d.max, BCAD=set$BCAD, rotate.axes=FALSE, d.res=400, yr.res=400, grey.res=100, dark=c(), colours=grey(seq(1, 0, length=100)), xaxt="s", yaxt="s", yr.lim) {
  dseq <- seq(d.min, d.max, length=d.res)
  Bacon.hist(dseq, set, BCAD=BCAD, calc.range=FALSE)
  hists <- get('hists') 
  
  scales <- array(NA, dim=c(d.res, yr.res))
  yrseq <- seq(min(yr.lim), max(yr.lim), length=yr.res)
  for(i in 1:length(dseq)) {
    yrs <- seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n)
    scales[i,] <- approx(yrs, hists[[i]]$counts, yrseq)$y
  }
  minmax <- hists[[length(hists)]]$min
  maxmax <- hists[[length(hists)]]$max
  scales <- scales/maxmax # normalise to the height of most precise age estimate
  if(length(dark) == 0)
    dark <- 10 * minmax/maxmax
  scales[scales > dark] <- dark 
  
  if(rotate.axes)
    image(yrseq, dseq, t(scales), add=TRUE, col=colours, useRaster=FALSE) else
      image(dseq, yrseq, scales, add=TRUE, col=colours, useRaster=FALSE)
}



# Time series of the log of the posterior
.PlotLogPost <- function(set, from=0, to=set$Tr)
  plot(from:(to-1), -set$Us[(from+1):to], type="l",
    ylab="Log of Objective", xlab="Iteration", main="")



# plot the prior for the accumulation rate
.PlotAccPrior <- function(s, mn, set=get('info'), main="", xlim=c(0, 3*max(mn)), xlab=paste("Acc. rate (yr/", noquote(set$unit), ")", sep=""), ylab="Density", add=FALSE, legend=TRUE, cex=.9) {
  o <- order(s, decreasing=TRUE)
  priors <- unique(cbind(s[o],mn[o])[,1:2])
  x <- 0
  if(length(priors) == 2) {
    curve(dgamma(x, s, s/mn), col=3, lwd=2, from=0, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("acc.shape: ", priors[1], "\nacc.mean: ", priors[2])
  } else {
	  priors <- priors[order(priors[,1]*priors[,2]),]
      curve(dgamma(x, priors[1,1], priors[1,1]/priors[1,2]), col=3, lwd=2, from=0, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
      for(i in 2:nrow(priors))
        curve(dgamma(x, priors[i,1], priors[i,1]/priors[i,2]), col=3, lwd=2, from=0, xlim=xlim, xlab=xlab, ylab=ylab, add=if(i==1) add else TRUE)
      txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
    }
  if(legend)
    legend("topright", txt, bty="n", cex=cex, text.col=2, adj=c(0,.2))
}



# plot the prior for the memory (= accumulation rate varibility between neighbouring depths)
.PlotMemPrior <- function(s, mn, thick, ds=1, set=get('info'), xlab="Memory (ratio)", ylab="Density", main="", add=FALSE, legend=TRUE, cex=.9) {
  o <- order(s, decreasing=TRUE)
  priors <- unique(cbind(s[o],mn[o])[,1:2])
  x <- 0

  if(length(priors)==2) {
    curve(dbeta(x, s*mn, s*(1-mn)), from=0, to=1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("mem.strength: ", s, "\nmem.mean: ", mn, "\n", set$K, " ", round(thick,3), noquote(set$unit), " sections", sep="")
  } else {
	  priors <- priors[order(priors[,1]*priors[,2]),]
      curve(dbeta(x, priors[1,1]*priors[1,2], priors[1,1]*(1-priors[1,2])), from=0, to=1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
      for(i in 2:nrow(priors))
        curve(dbeta(x, priors[i,1]*priors[i,2], priors[i,1]*(1-priors[i,2])), from=0, to=1, col=3, lwd=2, xlab="", ylab="", add=TRUE)
      txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
    }
  if(legend)
    legend("topright", txt, bty="n", cex=cex, text.col=2, adj=c(0,.2))
  warn <- FALSE
  for(i in s)
    for(j in mn)
      if(i*(1-j) <= 1) warn <- 1
  if(warn)
    cat("\nWarning! Chosen memory prior might cause problems.\nmem.strength * (1 - mem.mean) should be smaller than 1\n ")
}



# plot the prior for the hiatus length
.PlotHiatusPrior <- function(mx=set$hiatus.max, hiatus=set$hiatus.depths, set=get('info'), xlab="Hiatus size (yr)", ylab="Density", main="", xlim=c(0, 1.1*max(mx)), add=FALSE, legend=TRUE) {
  if(add)
    lines(c(0, 0, mx, mx), c(0, 1/mx, 1/mx, 0), col=3, lwd=2) else
      plot(c(0, 0, mx, mx), c(0, 1/mx, 1/mx, 0), xlab=xlab, ylab=ylab, xlim=xlim, type="l", col=3, lwd=2)
      
  txt <- paste("hiatus.max: ", toString(mx))
  if(legend)
    legend("topright", txt, bty="n", cex=.7, text.col=2, adj=c(0,.2))
}



# plot the posterior (and prior) of the accumulation rate
.PlotAccPost <- function(set=get('info'), s=set$acc.shape, mn=set$acc.mean, main="", xlab=paste("Acc. rate (yr/", set$unit, ")", sep=""), ylab="Frequency") {
  hi <- 2:(set$K-1)
  if(!is.na(set$hiatus.depths)[1])
    for(i in set$hiatus.depths) 
      hi <- hi[-max(which(set$d < i))]
  post <- c()
  for(i in hi) 
    post <- c(post, set$output[[i]])
  post <- density(post, from=0)
  post <- cbind(c(0, post$x, max(post$x)), c(0, post$y, 0))
  post <<- post
  maxprior <- dgamma((s-1)/(s/mn), s, s/mn)
  if(is.infinite(max(maxprior))) 
    max.y <- max(post[,2]) else
      max.y <- max(maxprior, post[,2])
  lim.x <- range(0, post[,1], 2*mn)
  plot(0, type="n", xlim=lim.x, xlab=xlab, ylim=c(0, 1.05*max.y), ylab="")
  polygon(post, col=grey(.8), border=grey(.4))
  .PlotAccPrior(s, mn, add=TRUE, xlim=range(post[,1]), xlab="", ylab=ylab, main=main)
}



# plot the posterior (and prior) of the memory
.PlotMemPost <- function(set=get('info'), corenam, K, main="", s=set$mem.strength, mn=set$mem.mean, xlab=paste("Memory"), ylab="Density", ds=1, thick) {
  post <- density(set$output[,set$n]^(1/set$thick), from=0, to=1)
  post <- cbind(c(min(post$x), post$x, max(post$x)), c(0, post$y, 0))
  maxprior <- max(dbeta((0:100)/100, s*mn, s*(1-mn)))
  if(is.infinite(max(maxprior))) max.y <- max(post[,2]) else
    max.y <- max(maxprior, max(post[,2]))
  plot(0, type="n", xlab=xlab, xlim=range(post[,1]), ylim=c(0, 1.05*max.y), ylab="", main="")
  polygon(post, col=grey(.8), border=grey(.4))
  .PlotMemPrior(s, mn, thick, add=TRUE, xlab="", ylab=ylab, main=main)
}



# plot the posterior (and prior) of the hiatus
.PlotHiatusPost <- function(set=get('info'), mx=set$hiatus.max, main="", xlim=c(), xlab=paste("Hiatus size (yr)", sep=""), ylab="Frequency", after=set$after) {
  gaps <- c()
  for(i in set$hiatus.depths) {
    below <- Bacon.Age.d(i+after, set)
    above <- Bacon.Age.d(i-after, set)
    gaps <- below - above 
#    reversals <- which(below < above) 
#    gaps <- c(gaps, below[-reversals] - above[-reversals])
  }
  gaps <- density(gaps, from=0)
  max.y <- 1.1*max(1/mx, gaps$y)
  if(length(xlim) == 0)
    xlim <- c(0, 1.1*mx)
  plot(0, type="n", main="", xlab=xlab, xlim=xlim, ylab=ylab, ylim=c(0, max.y))
  polygon(cbind(c(min(gaps$x), gaps$x, max(gaps$x)), c(0,gaps$y,0)),
    col=grey(.8), border=grey(.4))
  .PlotHiatusPrior(add=TRUE, xlab="", ylab=ylab, main=main)
}



