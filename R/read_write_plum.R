# Do a regression to determine the optimal numbers of supported data to use; the minimum is 3
.check.equi <- function(dets) {
  rawdata = dets[,4] #210Pb
  rawsd   = dets[,5] #sd(210pb)
  deps    = dets[,2] #depth

  lendat = length(rawdata)
  numdat = as.integer(.5*length(rawdata))
  usedat = rawdata[(lendat-3):lendat]
  usesd  = rawsd[(lendat-3):lendat]
  usex   = 1:length((lendat-3):lendat)
  usereg = lm(usedat ~ usex, weights=1/(usesd^2))
  reg    = coef(summary(usereg))[2,4]
  est    = coef(summary(usereg))[1,1]
  coe    = 3
  for(i in 1:numdat) {
    usedat = rawdata[(lendat-3-i):lendat]
    usesd  = rawsd[(lendat-3-i):lendat]
    usex   = 1:length((lendat-3-i):lendat)
    usereg = lm(usedat ~ as.numeric(scale(usex)), weights=1/(usesd^2))
    reg1   = coef(summary(usereg))[2,4]
    est1   = mean(usedat) #coef(summary(usereg))[1,1]
    if(reg1 > reg) {
      reg  = reg1
      coe  = (3+i)
      est  = est1
    }
  }

  ans <- readline(message("The regression process proposes using the last", as.integer(coe), "data points as estimates of the supported activity, with a p-value of", round((reg),3), "OK? (y/n) "))
  if(!(ans=="y" || ans=="")) 
    stop("  OK. Please adapt settings.\n\n", call.=FALSE)

  c(coe, reg1)
}

# read the 210Pb dets file
.read.dets.plum <- function(core, coredir, n.supp, date.sample, set=get('info'), sep=",", dec=".", cc=1, Bqkg = TRUE) {
  # a relation between the name of column and its position in HP1C file
  # These are the columns of the plum file
  idColumn       = 1
  plumdataColumn = 4
  stdColumn      = 5
  depthColumn    = 2
  deltaColumn    = 6
  rhoColumn      = 3
  radonColumn    = 7
  sdRadonColumn  = 8
  radonCase      = -1

  csv.file <- paste(coredir,  core, "/", core, ".csv", sep="")
  dat.file <- paste(coredir,  core, "/", core, ".dat", sep="")

  suggested.names <- c("labID","depth(cm)","density(g/cm^3)","210Pb(Bq/kg)","sd(210Pb)","thickness(cm)", "226Ra(Bq/kg)", "sd(226Ra)")
  changed <- 0

  #Read file
  if(file.exists(csv.file)) {
    dets <- read.table(csv.file, header=TRUE, sep=sep)
    if(file.exists(dat.file)) # deal with old .dat files
      if(file.info(csv.file)$mtime < file.info(dat.file)$mtime)
        message("Warning, the .dat file is newer than the .csv file! I will read the .csv file. From now on please modify ", csv.file, ", not ", dat.file, " \n", sep="") else
          message("Reading", csv.file, "\n")
  } else {
    if(file.exists(paste0(csv.file, ".txt"))) {
      file.rename(paste0(csv.file, ".txt"), csv.file)
      message("Removing .txt extension from .csv file")
    } else {
      message("No .csv file found, reading", dat.file, "and converting it to .csv")
      dets <- read.table(dat.file, header=TRUE)
      changed <- 1
    }
  }

  name <- tolower(names(dets))
  commas <- grep(",,", readLines(csv.file)) # check if there are too many commas (e.g., lines with just commas)
  if(length(!is.na(commas)) > 0) # often an artefact of spreadsheet programs
    stop("check the .csv file in a plain-text editor for 'orphan' commas\n", call.=FALSE)

  #check that depths are in ascending order
  if(min(diff(dets[,depthColumn])) < 0) {
    message("Warning, the depths are not in ascending order, I will correct this")
    dets <- dets[ order(dets[,depthColumn]),]
    write.table(dets, csv.file, sep=sep, dec=dec, row.names=FALSE, quote=FALSE)
  }

  # check if the file has some parameters
  if(ncol(dets) == 7) {
    detsOrig <- dets [,-c(7)]

    if( !is.na( dets[1,7] ) ) # param defined on .csv file
      date.sample = dets[1,7]
    if( !is.na(dets[2,7]) )
      n.supp = dets[2,7]
    if( !is.na(dets[3,7]) )
      radonCase = dets[3,7]

    dets <- dets [,-c(7)]
  } else if( ncol(dets) == 9 ) {
    detsOrig <- dets [,-c(9)]

    if( !is.na( dets[1,9] ) ) #param defined on .csv file
      date.sample = dets[1,9]
    if( !is.na(dets[2,9]) )
      n.supp = dets[2,9]
    if( !is.na(dets[3,9]) )
      radonCase = dets[3,9]

    dets <- dets [,-c(9)]
  } else
      detsOrig <- dets

  if( !(is.na(n.supp)&&is.na(radonCase)) ){
    if( n.supp != 0 && radonCase==2 ){
      stop("The radon case can't be 2 when the number of supported measurements is different from 0; check the manual\n ", call.=FALSE )
    }
  }

  dlim = c(0, max(detsOrig[,depthColumn]))

  if(ncol(detsOrig) == 6) {
    age.min <- min( c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5]) )
    age.max <- max( c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5]) )
  } else {
    age.min <- min( c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5],detsOrig[,7]-detsOrig[,8]) )
    age.max <- max( c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5],detsOrig[,7]+detsOrig[,8]) )
  }

  age.lim <- extendrange(c(age.min, age.max), f=0.01)

  layout(matrix(c(1), nrow=1, byrow=TRUE), heights=c(1))
  oldpar <- par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l")
  on.exit(par(oldpar))

  ylab <- ifelse(Bqkg, '210Pb (Bq/kg)', '210Pb (dpm/g)')
  plot(0, type='n', pch=16,col=c(rep('red',nrow(detsOrig)),rep('red',nrow(detsOrig))),
    cex=.3, ylab=ylab, xlab='depth(cm)', xlim = dlim, ylim = age.lim )


  if( ncol(detsOrig) == 6 ) {
    segments(c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5]),
      c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5]),
      lwd=2, col=c(rep(3,nrow(detsOrig)), rep('red',nrow(detsOrig))))
  } else {
    segments(c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5],detsOrig[,7]-detsOrig[,8]),
      c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5],detsOrig[,7]+detsOrig[,8]),
      lwd =2, col=c(rep(3,nrow(detsOrig)), rep('red',nrow(detsOrig))))
  }

  if( core == "HP1C" ) {
    radonColumn = 4
    sdRadonColumn = 5
    supportedData = dets[30:33,c(radonColumn, sdRadonColumn)]
    dets = dets[1:29,]
    if( radonCase < 0 )
      radonCase = 0
    date.sample = 2018.5
  } else if( ncol(dets) == 6 ) {
    radonColumn = 4
    sdRadonColumn = 5

    if( n.supp == 0 ) { # the number of supported data
      if( radonCase < 0 )
        radonCase = 0
      supportedData <- NULL
    } else if( n.supp >0 ){
      if( radonCase < 0 )
        radonCase = 1
      supportedData = dets[(nrow(dets)-n.supp+1):nrow(dets),c(radonColumn, sdRadonColumn)]
      dets = dets[1:(nrow(dets)-n.supp),]
    } else { # do a linear regression to estimate the supported data values
      tmp = .check.equi( dets )
      n.supp = tmp[1]
      supportedData = dets[(nrow(dets)-n.supp+1):nrow(dets),c(radonColumn, sdRadonColumn)]
  #    dets = dets[1:(nrow(dets)-n.supp),]
      if( radonCase < 0 )
        radonCase = 1
    }
  } else if( ncol(dets)  == 8 ){
    radonColumn    = 7
    sdRadonColumn  = 8
    #columns 7 and 8 are supported data
    supportedData = dets[,c(7, 8)]
    dets = dets[,-c(radonColumn,sdRadonColumn)]

    if( max(is.na(supportedData)) > 0 ){
      message("Missing values are detected; the radon case is set to 1\n")

      elim <- NULL
      for(i in 1:nrow(supportedData)){
        if( max(is.na(supportedData[i,])) > 0 ){
          elim <- c( elim, i )
        }
      }

      supportedData = supportedData[ -elim, ]

      radonCase = 1
      ans <- readline(message("Additionaly, you can set a number of supported data to use, do you want to set a number of supported? (y/n)"))
      if(tolower(substr(ans, 1, 1)) == "y") {
        ans <- readline(cat("Ok, n.supp="))
        n.supp = as.integer(ans)
        radonColumn = 4
        sdRadonColumn = 5
        tmp <- dets[(nrow(dets)-n.supp+1):(nrow(dets)),c(radonColumn, sdRadonColumn)]
        names(tmp) <- colnames(supportedData)
        supportedData <- rbind( supportedData,  tmp)
        dets <- dets[1:(nrow(dets)-n.supp),]
      }


    }else if( n.supp > 0 ){
      if( radonCase < 0 )
        radonCase = 1
      radonColumn = 4
      sdRadonColumn = 5
      tmp <- dets[(nrow(dets)-n.supp+1):(nrow(dets)),c(radonColumn, sdRadonColumn)]
      names(tmp) <- colnames(supportedData)
      supportedData <- rbind( supportedData,  tmp)
      dets <- dets[1:(nrow(dets)-n.supp),]
    } else {
      message("\n Plum can assume to have a constant supported 210Pb and use the 226Ra data to infer this one value.")
      message(" Plum can also assume individual supported 210Pb values per measured depth.")
      message(" It is important to consider that this will greatly increase the computing time and it should only be used when clear patterns are observed in the 226Ra data.")
      ans <- readline(cat(" Do you want to use the individual supported 210Pb? (y/n) "))
      if( !tolower(substr(ans, 1, 1)) == "y") {
        message(" OK, using constant supported 201Pb\n")
        radonCase = 1
      } else {
        message(" OK, using individual supported 210Pb per data point.")
        radonCase = 2
      }
    }
  } else {
    stop("Unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.\n", call.=FALSE)
  }

  # more sanity checks for dets values
  if(!is.numeric(dets[,plumdataColumn]) || !is.numeric(dets[,stdColumn]) || !is.numeric(dets[,depthColumn]))
    stop("unexpected values in dets file, I expected numbers. Check the manual.\n", call.=FALSE)

  # more sanity checks for supported values
  if(!is.numeric(dets[,deltaColumn]) || !is.numeric(dets[,rhoColumn]) )
    stop("unexpected values in dets file, I expected numbers. Check the manual.\n", call.=FALSE)


  # if current dets differ from original .csv file, rewrite it
  if(changed > 0){

    write.table(dets, csv.file, sep=paste(sep, "\t", sep=""), dec=dec, row.names=FALSE, col.names=suggested.names[1:ncol(dets)], quote=FALSE)
  }

  dets = dets[,c(idColumn, plumdataColumn, stdColumn, depthColumn, deltaColumn, rhoColumn)]
  list(dets, supportedData, radonCase, date.sample, detsOrig, n.supp)
}


# read in default values, values from previous run, any specified values, and report the desired one. Internal function.
.plum.settings <- function(core, coredir, dets, thick, remember=TRUE, d.min, d.max, d.by, depths.file,
  slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, hiatus.shape,
  BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob,
  defaults, runname, ssize, dark, MinAge, MaxAge, cutoff, age.res, after, age.unit,
  supportedData, date.sample, Al, phi.shape, phi.mean, s.shape, s.mean, radonCase, Bqkg, n.supp) {

  vals <- list(d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob, age.unit) # do these now need the lengths of all the values, e.g. d.min=numeric(1) ?
  valnames <- c("d.min", "d.max", "d.by", "depths.file", "slump", "acc.mean", "acc.shape", "mem.mean", "mem.strength", "boundary", "hiatus.depths", "hiatus.max", "BCAD", "cc", "postbomb", "cc1", "cc2", "cc3", "cc4", "depth.unit", "normal", "t.a", "t.b", "delta.R", "delta.STD", "prob", "age.unit")
  #TODO: modificar para que acepte el vector de soportado y los valores propios de plum, como es el nombre de archivo de soportado y los paremetros como "Al"
  extr <- function(i, def=deffile, pre=prevfile, exists.pre=prevf, rem=remember, sep=" ", isnum=TRUE) {
    if(length(vals[[i]]) > 0) # tmp
      if(any(is.na(vals[[i]]))) {
        ext.def <- strsplit(def[i], sep)[[1]]
        ext.def <- ext.def[-length(ext.def)] # remove description
        if(exists.pre) {
          ext.pre <- strsplit(pre[i], sep)[[1]]
          ext.pre <- ext.pre[-length(ext.pre)] # remove description
          if(def[i] == pre[i]) # values for dev and pre similar, no worries
            ext <- ext.pre else
              if(rem) {
                if(i==13) ifelse(ext.pre, "using BC/AD", "using cal BP") else
                if(i>2) message(" using previous run's value for ", valnames[i], ", ", ext.pre)
                ext <- ext.pre
              } else {
                  if(i==13) ifelse(ext.def, "using BC/AD", "using cal BP") else
                  if(i>2) message(" using default value for ", valnames[i], ", ", ext.def)
                  ext <- ext.def
                }
        } else ext <- ext.def

        if(any(ext=="NA") || any(is.na(ext))) NA else
          if(isnum) as.numeric(ext) else noquote(ext)
      } else
        if(isnum) as.numeric(vals[[i]]) else vals[[i]]
  }

  # read in default values and those of previous run if available
  deffile <- readLines(defaults, n=-1)
  prevfile <- paste(coredir, core, "/", core, "_settings.txt", sep="")
  prevf <- FALSE
  if(file.exists(prevfile)) {
    prevfile <- readLines(prevfile, n=-1)
    if(length(prevfile) > 0) prevf <- TRUE
  }

  #d.min <- extr(1); d.by <- extr(3); depths.file <- extr(4)
  #slump <- extr(5); acc.mean <- extr(6);
  #if(length(acc.shape) == 1)
  #  acc.shape <- extr(7)
  #mem.mean <- extr(8)
  #mem.strength <- extr(9)
  #boundary <- if(is.na(boundary)[1]) NA else sort(extr(10))
  #hiatus.depths <- if(is.na(hiatus.depths)[1]) NA else sort(extr(11))
  #hiatus.max <- extr(12)
  #BCAD <- extr(13); cc <- extr(14); postbomb <- extr(15); cc1 <- extr(16, isnum=FALSE)
  #cc2 <- extr(17, isnum=FALSE); cc3 <- extr(18, isnum=FALSE); cc4 <- extr(19, isnum=FALSE)
  #depth.unit <- extr(20, isnum=FALSE); normal <- extr(21); t.a <- extr(22); t.b <- extr(23)
  #delta.R <- extr(24); delta.STD <- extr(25); prob <- extr(26); age.unit <- extr(27, isnum=FALSE)

  if(is.na(d.min) || d.min=="NA")
    d.min <- min(dets[,4])
  if(is.na(d.max) || d.max=="NA")
    d.max <- max(dets[,4])
  if(length(acc.shape) < length(acc.mean))
    acc.shape <- rep(acc.shape, length(acc.mean)) else
      if(length(acc.shape) > length(acc.mean))
        acc.mean <- rep(acc.mean, length(acc.shape))
  if(length(mem.strength) < length(mem.mean))
    mem.strength <- rep(mem.strength, length(mem.mean)) else
      if(length(mem.strength) > length(mem.mean))
        mem.mean <- rep(mem.mean, length(mem.strength))

  ## produce/update settings file, and return the values
  prevfile <- file(paste(coredir, core, "/", core, "_settings.txt", sep=""), "w")
  scat <- function(m, n="") cat(m, n, sep="", file=prevfile)
  cat(d.min, " #d.min\n", d.max, " #d.max\n", d.by, " #d.by\n",
    depths.file, " #depths.file\n", slump, " #slump\n", sep="", file=prevfile)
  for(i in acc.mean) scat(i, " "); scat("#acc.mean\n")
  for(i in acc.shape) scat(i, " "); scat("#acc.shape\n", "")
  for(i in mem.mean) scat(i, " "); scat("#mem.mean\n", "")
  for(i in mem.strength) scat(i, " "); scat("#mem.strength\n", "")
  for(i in boundary) scat(i, " "); scat("#boundary\n", "")
  for(i in hiatus.depths) scat(i, " "); scat("#hiatus.depths\n", "")
  for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "")
  #for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "") # redundant
  cat(BCAD, " #BCAD\n", cc, " #cc\n", postbomb, " #postbomb\n",
    cc1, " #cc1\n", cc2, " #cc2\n", cc3, " #cc3\n", cc4, " #cc4\n",
    depth.unit, " #depth.unit\n", normal, " #normal\n", t.a, " #t.a\n", t.b, " #t.b\n",
    delta.R, " #delta.R\n", delta.STD, " #d.STD\n", prob, " #prob\n", age.unit, "#age.unit\n", sep="", file=prevfile)

  cat(date.sample, " #date.sample\n", Al, " #Al\n", phi.shape, " #phi.shape\n", phi.mean, " #phi.mean\n",
    s.shape, " #s.shape\n", s.mean, " #s.mean\n", radonCase, " #radonCase\n", Bqkg, " #Bqkg\n", sep="", file=prevfile)

  cat(n.supp, " #n.supp\n", sep="", file=prevfile);

  close(prevfile)

  if(length(MinAge) == 0)
    MinAge <- min(1950 - as.integer(format(Sys.time(), "%Y")), round(dets[,2] - (5*dets[,3])))
  if(length(MaxAge) == 0)
    MaxAge <- max(1e6, round(dets[,2] + (5*dets[,3])))

  theta0 = 1950 - date.sample

  list(core=core, thick=thick, dets=dets, d.min=d.min, d.max=d.max, coredir=core,
    d.by=d.by, depths.file=depths.file, slump=slump,
    acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean,
    mem.strength=mem.strength, boundary=boundary,
    hiatus.depths=hiatus.depths, hiatus.max=hiatus.max,
    BCAD=BCAD, cc=cc, postbomb=postbomb,
    cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=noquote(depth.unit), unit=depth.unit, age.unit=noquote(age.unit), normal=normal,
    t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, date=date(),
    runname=runname, ssize=ssize, dark=dark, MinAge=MinAge, MaxAge=MaxAge,
    cutoff=cutoff, age.res=age.res, after=after,
    supportedData=supportedData, theta0 = theta0, Al=Al, phi.shape=phi.shape, phi.mean=phi.mean, s.shape=s.shape, s.mean=s.mean,
    radonCase=radonCase, Bqkg=Bqkg)
}

#function to merge dets of plum and bacon data
.merge.dets <- function(detsPlum, detsBacon, delta.R, delta.STD, t.a, t.b, cc){
  if( ncol(detsBacon) >= 5 ){
    cc <- detsBacon[,5]
    detsBacon <- detsBacon[,-5]
  }else{
    cc <- array(cc, dim=c(nrow(detsBacon),1))
  }

  if( ncol(detsBacon) < 9 ){

    for(i in (ncol(detsBacon)+1):9){
      if( i == 5){
        col <- array(delta.R, dim=c(nrow(detsBacon),1))
      }else if(i == 6){
        col <- array(delta.STD, dim=c(nrow(detsBacon),1))
      }else if(i == 7){
        col <- array(t.a, dim=c(nrow(detsBacon),1))
      }else if(i == 8){
        col <- array(t.b, dim=c(nrow(detsBacon),1))
      }else if(i==9){
        col <- cc
      }
      detsBacon <- cbind(detsBacon, col)
    }
    colnames(detsBacon) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")
    #print(detsBacon)
  }

  if( ncol(detsPlum) < 9 ){
    for(i in (ncol(detsPlum)+1):9){
      if( i == 5){
        col <- array(delta.R, dim=c(nrow(detsPlum),1))
      }else if(i == 6){
        col <- array(delta.STD, dim=c(nrow(detsPlum),1))
      }else if(i == 7){
        col <- array(t.a, dim=c(nrow(detsPlum),1))
      }else if(i == 8){
        col <- array(t.b, dim=c(nrow(detsPlum),1))
      }else if(i==9){
        col <- array(5, dim=c(nrow(detsPlum),1))
      }
      detsPlum <- cbind(detsPlum, col)
    }
    colnames(detsPlum) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")
    #print(detsPlum)
  }

  dets <- rbind(detsPlum, detsBacon, make.row.names = FALSE)
  dets <- dets[ order(dets[,4]),]

}

# write files to be read by the main Bacon age-depth modelling function
.write.plum.file <- function(set=get('info')) {

  #a relation between the name of column and his position
  #These are the column of the plum file
  idColumn       = 1
  plumdataColumn = 2
  stdColumn      = 3
  depthColumn    = 4
  deltaColumn    = 5
  rhoColumn      = 6

  if(length(set$slump) > 0) {
    dets <- set$slumpdets
    hiatus.depths <- set$slumphiatus
    boundary <- set$slumpboundary
  } else {
    dets <- set$dets
    hiatus.depths <- set$hiatus.depths
    boundary <- set$boundary
  }

  if(set$d.min < min(dets[,depthColumn])) { # repeat relevant row, change error and depth
    # extrap <- c(NA, min(dets[,2]), max(1e5, 1e3*dets[,2], 1e3*dets[,3]), set$d.min, 0)
    dets <- rbind(dets[which(dets[,depthColumn] == min(dets[,depthColumn]))[1],], dets, make.row.names=FALSE)
    dets[1,1] <- NA # calling this "d.min" causes issues
    dets[1,3] <- max(1e5, 1e3*dets[,4], 1e3*dets[,3])
    dets[1,depthColumn] <- set$d.min
  }
  #print(set$d.max)
  #print(max(dets[,depthColumn]))
  if(set$d.max > max(dets[,depthColumn])) { # repeat relevant row, change error and depth
    # extrap <- c(NA, max(dets[,2]), max(1e5, 1e3*dets[,2], 1e3*dets[,3]), set$d.max, 0)
    dets <- rbind(dets, dets[which(dets[,depthColumn] == max(dets[,depthColumn]))[1],], make.row.names=FALSE)
    dets[nrow(dets),1] <- NA # calling this "d.max" causes issues
    dets[nrow(dets),3] <- max(1e5, 1e3*dets[,4], 1e3*dets[,3])
    dets[nrow(dets),depthColumn] <- set$d.max
  }

  #print(dets)

  supportedData <- set$supportedData

  fl <- file(set$bacon.file, "w")
  cat("## Ran on", set$date, "\n\n", file=fl)
  cat("Cal 0 : ConstCal;\nCal 1 : ",
  if(set$cc1=="IntCal13" || set$cc1=="\"IntCal13\"") "IntCal13"
    else noquote(set$cc1), ", ", set$postbomb, ";\nCal 2 : ",
  if(set$cc2=="Marine13" || set$cc2=="\"Marine13\"") "Marine13"
    else noquote(set$cc2), ";\nCal 3 : ",
  if(set$cc3=="SHCal13" || set$cc3=="\"SHCal13\"") "SHCal13"
    else noquote(set$cc3), ", ", set$postbomb, ";",
			if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") set$cc4 <- NULL
    else
      paste("\nCal 4 : GenericCal, ", set$cc4, ";", sep=""), sep="", file=fl)
  cat("\nCal 4 : ConstCal;", sep="", file=fl)

  cat("\n##          alPhi mPhi  alS  mS     Al   theta0  Radon_case  supported_data_file", file=fl)

  cat("\nCal 5 : Plum, ", set$phi.shape, ", ",  set$phi.mean, ", ",  set$s.shape, ", ", set$s.mean, ", ", set$Al, ", ", set$theta0, ", ",
        set$radonCase, ", ", set$plum.file,";", sep="", file=fl)

  #cat("\n##   id.    210Pb   std   depth   delta     rho  t.a t.b cc ... Plum: 210Pb data", file=fl)
  cat("\n##    ", colnames(dets), " ... Plum: 210Pb data",sep=", ", file=fl)

  # we need to send the dets with all columns so pre-processing is needed
  for( i in 1:nrow(dets) ){
    cat( "\nDet ", i-1, " : ", as.character(dets[i,1]),
        " , ", dets[i,2],
        ", ", dets[i,3],
        ", ", dets[i,4],
        ", ", dets[i,5],
        ", ", dets[i,6],
        ", ", dets[i,7],
        ", ", dets[i,8],
        ", ", dets[i,9],
        ";", sep="", file=fl)
  }

  if(!is.na(hiatus.depths[1])) {
    if(is.null(boundary[1]))
      message("\n  Hiatus set at depth(s) ", hiatus.depths) else
        message("  Boundary set at depth(s) ", boundary)
    if(length(set$acc.shape)==1)
      set$acc.shape <- rep(set$acc.shape, length(hiatus.depths)+1)
    if(length(set$acc.mean)==1)
      set$acc.mean <- rep(set$acc.mean, length(hiatus.depths)+1)
    if(length(set$hiatus.max)==1)
      set$hiatus.max <- rep(set$hiatus.max, length(hiatus.depths))
#      if(length(set$hiatus.shape)==1)
#        set$hiatus.shape <- rep(set$hiatus.shape, length(set$hiatus.depths))
    .assign_to_global("info", set)
    cat("\n\n### Depths and priors for fixed hiatuses, in descending order",
      "\n##### cm  alpha beta      ha     hb", file=fl)
    for(i in length(hiatus.depths):1)
      cat("\nHiatus ", i-1, ":  ", hiatus.depths[i], ",  ", set$acc.shape[i+1],
        ",  ", set$acc.shape[i+1]/set$acc.mean[i+1], ",  ", .1, # last value (h.a) was NA but this conflicts with setting initial values for hiatus length
        ",  ", set$hiatus.max[i], ";", sep="", file=fl)
  }

  cK <- set$d.min+(set$thick*set$K)
  ### final parameters - dmax now calculated as dmin+(dC*K)
  if( is.na(set$seed) ){
  wrapup <- paste("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax",
    "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
    ",  ", set$theta0-.02, ",  ", 26500, ",  ", set$theta0-0.01, ",  ", set$theta0+0.01,
    ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
    ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
    ", ", cK, ";\n", sep="")
  }else{
    wrapup <- paste("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax  seed",
      "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
      ",  ", set$theta0-.02, ",  ", 26500, ",  ", set$theta0-0.01, ",  ", set$theta0+0.01,
      ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
      ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
      ", ", cK, ", ", set$seed, ";\n", sep="")
  }
  cat(wrapup, file=fl)
  close(fl)

  fl <- file(set$plum.file, "w")
  for (i in 1:nrow(supportedData)){
    for (j in 1:ncol(supportedData)){
      cat( supportedData[i,j], " ", sep="", file= fl )
    }
    cat("\n", file=fl)
  }
  close(fl)
}

# function to read output files into memory
.Plum.AnaOut <- function(fnam, set=get('info')) {
  out <- read.table(fnam)
  n <- ncol(out)-1
  set$nPs  <- n
  set$TrPs <- nrow(out)
  set$phi  <- out[,1]
  set$ps   <- out[,2:(n+1)]
  set
}

# read the dets file, converting old formats to new ones if so required
.read.dets.plumbacon <- function(core, otherdates, coredir, set=get('info'), sep=",", dec=".", cc=1) {
  # if a .csv file exists, read it (checking that it is more recent than any .dat file in the folder). Otherwise, read the .dat file, check the columns, report back if >4 (>5?) columns, and convert to .csv (report this also)
  csv.file <- paste(coredir,  core, "/", otherdates, ".csv", sep="")
  dat.file <- paste(coredir,  core, "/", otherdates, ".dat", sep="")

  dR.names <- c("r", "d", "d.r", "dr", "deltar", "r.mn", "rm", "rmn", "res.mean", "res.mn", "delta.r")
  dSTD.names <- c("d.std", "std", "std.1", "dstd", "r.std", "rstd", "res.sd", "delta.std", "deltastd")
  ta.names <- c("t", "t.a", "ta", "sta")
  tb.names <- c("t", "t.b", "tb", "stb")
  cc.names <- c("c", "cc")
  suggested.names <- c("labID", "age", "error", "depth", "cc", "dR", "dSTD", "ta", "tb")
  changed <- 0

  if(file.exists(csv.file)) {
    dets <- read.table(csv.file, header=TRUE, sep=sep)
    if(file.exists(dat.file)) # deal with old .dat files
      if(file.info(csv.file)$mtime < file.info(dat.file)$mtime)
        message("Warning, the .dat file is newer than the .csv file! I will read the .csv file. From now on please modify ", csv.file, ", not ", dat.file) else
          message("Reading", csv.file)
    } else {
      if(file.exists(paste0(csv.file, ".txt"))) {
        file.rename(paste0(csv.file, ".txt"), csv.file)
        message("Removing .txt extension from .csv file")
      } else {
        message("No .csv file found, reading", dat.file, "and converting it to .csv")
        dets <- read.table(dat.file, header=TRUE)
        changed <- 1
        }
    }
  name <- tolower(names(dets))
  commas <- grep(",,", readLines(csv.file)) # check if there are too many commas (e.g., lines with just commas)
  if(length(!is.na(commas)) > 0) # often an artifact of spreadsheet programs
    stop("check the .csv file in a plain-text editor for 'orphan' commas\n", call.=FALSE)

  # check if 'classic' dets file, which has a different column order from the current default
  if(ncol(dets) > 4)
    if(ncol(dets) == 5) { # then probably a 'new' dets file
      if((name[5] %in% cc.names) && min(dets[,5]) >= 0 && max(dets[,5]) <= 4) {} else # extra check for correct values
        stop("unexpected name or values in fifth column (cc, should be between 0 and 4). Please check the manual for guidelines in producing a correct .csv file.\n", call.=FALSE)
    } else
      if(ncol(dets) == 6) { # probably an 'old' file: dR, dSTD, but could also be cc and delta.R (so no column for delta.STD)
        if(name[5] %in% dR.names && name[6] %in% dSTD.names) {
			cat("\nHELP!!! 6!!!\n")
          dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
          message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file)
          changed <- 1
        } else
	      stop("unexpected names for columns 5/6. If you want to include delta.R, also add a column for delta.STD. Check the manual for guidelines to producing a correct .csv file.\n", call.=FALSE)
      } else
        if(ncol(dets) == 7) { # probably a 'new' file: cc, dR, dSTD
          if(name[5] %in% cc.names && min(dets[,5]) >= 0 && max(dets[,5]) <= 4 &&
            name[6] %in% dR.names && name[7] %in% dSTD.names)
              {} else
                 stop("unexpected column names, order or values in dets file. \nPlease check the manual for correct dets file formats.\n", call.=FALSE)
        } else
          if(ncol(dets) == 8) { # probably an 'old' file: dR, dSTD, ta, tb
            if(name[5] %in% dR.names && name[6] %in% dSTD.names)
            if(name[7] %in% ta.names && name[8] %in% tb.names)
            if(range(dets[,8] - dets[,7]) == c(1,1)) { # check that these set expected student-t values
              dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
              message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file)
              changed <- 1
            } else
              stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file", call.=FALSE)
          } else
            if(ncol(dets) == 9) { # most complex case, many checks needed
              if(name[9] %in% cc.names && # we're almost sure that this is a 'classic' dets file
                min(dets[,9]) >= 0 && max(dets[,9]) <= 4 && # check that this sets calibration curves
                  range(dets[,8] - dets[,7]) == c(1,1) && # check that these set expected student-t values
                    name[5] %in% dR.names && name[6] %in% dSTD.names && # column names as expected?
                      name[7] %in% ta.names && name[8] %in% tb.names) { # column names as expected?
                        dets <- dets[,c(1:4,9,5:8)] # shuffle colums around
                        message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD, t.a, t.b, cc. \nAdapting column order and saving as", csv.file)
                        changed <- 1
                      } else
                        if(name[5] %in% cc.names && # oh, probably a 'new' file from more recent Bacon
                          min(dets[,5]) >= 0 && max(dets[,5]) <= 4 && # check that this sets cal.curves
                            range(dets[,9] - dets[,8]) == c(1,1) && # columns 8-9 set student-t correctly
                              name[8] %in% ta.names && name[9] %in% tb.names && # and are correctly named
                                name[6] %in% dR.names && name[7] %in% dSTD.names) # all lights are green
                                  {} else
                                     stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file", call.=FALSE)
            } else
              stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.\n", call.=FALSE)

  # more sanity checks
  if(!is.numeric(dets[,2]) || !is.numeric(dets[,3]) || !is.numeric(dets[,4]))
    stop("unexpected values in dets file, I expected numbers. Check the manual.\n", call.=FALSE)
  if(min(dets[,3]) <= 0) {
    message("Warning, zero year errors don't exist in Bacon's world. I will increase them to 1 ", set$age.unit, " yr")
    dets[dets[,3] <= 0,3] <- 1
    changed <- 1
  }
  if( nrow(dets) > 1 && min(diff(dets[,4])) < 0) {
    message("Warning, the depths are not in ascending order, I will correct this")
    dets <- dets[ order(dets[,4]), ]
    changed <- 1
  }

  # if current dets differ from original .csv file, rewrite it
  if(changed > 0)
    write.table(dets, csv.file, sep=paste(sep, "\t", sep=""), dec=dec, row.names=FALSE, col.names=suggested.names[1:ncol(dets)], quote=FALSE)
  dets
}
