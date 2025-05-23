validateDirectoryName <- function(dir) {
  if(!dir.exists(dir))
    dir.create(dir, recursive=TRUE)
  dir <- suppressWarnings(normalizePath(dir))
  lastchar <- substr(dir, nchar(dir), nchar(dir))
  if(lastchar != "/" & lastchar != "\\" & lastchar != "" & lastchar != "." )
    dir <- paste0(dir, "/") # does this work in Windows?
  return(dir)
}



#' @name set.initvals
#' @title Set initial values for the Bacon MCMC run.
#' @description Select initial values th0 and th1 for a Bacon MCMC run and write them into a file that can be read by Bacon.
#' @details By default, the initial MCMC values th0 and th1 of the Bacon age-depth model (upper ages and accumulation rate
#' for each model section) are estimated randomly. Since version 3.1.0, these starting values can
#' also be provided in a file with extension _bacon.init, placed
#' within the core's folder. This file will need to have two rows,
#' each for one of the two initial sets of parameters required (the
#' t-walk requires two starting estimates for all MCMC parameters).
#' If such a file is found (and correctly formatted), Bacon will use
#' the values within this file as starting points for the MCMC run.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A .bacon.init file
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param core The name of the core for which a bacon.init file needs to be made
#' @param values use this if you wish to provide the values (2 rows with starting age, accumulation rates for each model section, and memory parameter w).
#' @param click use this if you wish to use the cursor to manually select age-depth points from the current graphic device. This is the default option. Right click once you have selected all datapoints, or provide the number of expected datapoints as a value (e.g., \code{click=5}).
#' @export
set.initvals <- function(set=get('info'), core=set$core, values=c(), click=1) {
  elbows <- set$elbows

  if(length(values) == 0) { # then the user will click on an existing graph to select age-depth points
    message("Please select your initval age-depth points on the graph by clicking with your left mouse button. When you are done selecting, use right-click.")

    if(set$isplum)
      if(set$hasBaconData) {
        par(mfg=c(2,1)) # bottom-left panel
        message("Please select your data from the bottom-left panel with the cal BP/C14 dates")
      }

    if(click == 1) # user right-clicks once done
      draft.agedepth <- locator() else
        draft.agedepth <- locator(click)
    points(draft.agedepth, col=2, pch=4, cex=5, lwd=2)
    message("Lovely.")

    if(set$isplum)
      if(set$hasBaconData)
        par(mfg=c(2,2)) # back to bottom-right panel

    ages <- approx(draft.agedepth$x, draft.agedepth$y, elbows, rule=2)$y
    agedepth <- cbind(elbows, ages)
    accs <- diff(ages) / diff(elbows)
    accs <- c(accs[1], accs)
    accs1 <- abs(jitter(accs)) # add some scatter to ensure that th
    accs2 <- abs(jitter(accs))
    ages1 <- jitter(ages[1])
    ages2 <- jitter(ages[2])
    w1 <- jitter(0.5) # memory
    w2 <- jitter(0.5)
    init1 <- c(ages1, accs1, w1)
    init2 <- c(ages2, accs2, w2)
    } else { # user provides 2 sets of starting points
      init1 <- values[1,]
      init2 <- values[2,]
    }
  write.table(rbind(init1, init2), paste0(set$bacon.file, ".init"), sep="\t", col.names=FALSE, row.names=FALSE)
}



#' @name clam2bacon
#' @title Translate clam .csv files to Bacon .csv files.
#' @description Reads a clam .csv file containing the dates, and transforms it into a Bacon .csv file.
#' @details Please ensure that if the clam file has offset (d.R) estimates, that errors (d.STD) are provided manually, since these values cannot be determined automatically from the clam .csv file.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A Bacon .csv file
#' @param core The name of the core for which a clam .csv.file needs to be translated into a Bacon .csv file
#' @param clamdir The directory where the clam runs reside. Defaults to \code{coredir="clam_runs"}.
#' @param bacondir The directory where the Bacon runs reside. Defaults to \code{coredir="Bacon_runs"}.
#' @param sep The separator for the .csv files. Defaults to \code{sep=","}.
#' @param cc Calibration curve for C-14 dates: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine),
#' @export
clam2bacon <- function(core, clamdir="clam_runs", bacondir="Bacon_runs", sep=",", cc=1) {
  clamfl <- read.csv(paste0(clamdir, "/", core, "/", core, ".csv"), sep=sep)
  ID <- as.character(clamfl[,1])
  C14 <- clamfl[,2]
  calBP <- clamfl[,3]
  error <- clamfl[,4]
  res <- clamfl[,5]
  d <- clamfl[,6]

  cc <- rep(cc, nrow(clamfl))
  cc[which(!is.na(calBP))] <- 0 # cal BP ages get cc=0
  ages <- C14
  ages[which(is.na(ages))] <- calBP[which(!is.na(calBP))]

  if(length(res[!is.na(res)]) > 0)
    message("Warning, please add reservoir effects manually to the bacon .csv file")

  baconfl <- cbind(ID, ages, error, d, cc)
  bacondir <- paste0(bacondir, "/", core)
  if(!dir.exists(bacondir))
    dir.create(bacondir)
  fastwrite(baconfl, paste0(bacondir, "/", core, ".csv"), sep=sep, row.names=FALSE, quote=FALSE)
}



#' @name bacon2clam
#' @title Translate Bacon .csv files to clam .csv files.
#' @description Reads a Bacon .csv file containing the dates, and transforms it into a clam .csv file.
#' @details Assumes that Bacon .csv files with 4 columns indicate 14C dates. Please make sure this is correct.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A clam .csv file
#' @param core The name of the core for which a Bacon .csv.file needs to be translated into a clam .csv file
#' @param bacondir The directory where the Bacon runs reside. Defaults to \code{coredir="Bacon_runs"}.
#' @param clamdir The directory where the clam runs reside. Defaults to \code{coredir="clam_runs"}.
#' @param sep The separator for the .csv files. Defaults to \code{sep=","}.
#' @examples{
#' \donttest{
#'  tmpfl <- tempfile()
#'   Bacon(run=FALSE, ask=FALSE, coredir=tmpfl)
#'   bacon2clam("MSB2K", bacondir=tmpfl, clamdir=tmpfl)
#'  }
#' }
#' @export
bacon2clam <- function(core, bacondir="Bacon_runs", clamdir="clam_runs", sep=",") {
  baconfl <- read.csv(paste0(bacondir, "/", core, "/", core, ".csv"), sep=sep)
  ID <- as.character(baconfl[,1])
  ages <- baconfl[,2]
  error <- baconfl[,3]
  d <- baconfl[,4]

  C14 <- rep(NA, length(ages))
  calBP <- C14
  res <- C14

  if(ncol(baconfl) == 4) {
    C14 <- ages
    message("Warning, I am assuming that these are radiocarbon dates, not cal BP")
  }
  if(ncol(baconfl) >= 5) {# then cc specified
    calBP[which(baconfl[,5] == 0)] <- ages[which(baconfl[,5] == 0)]
    C14[which(baconfl[,5] > 0)] <- ages[which(baconfl[,5] > 0)]
  }
  if(ncol(baconfl) >= 7) { # then also offsets specified
    res <- baconfl[,6]
    error <- sqrt(error^2 + baconfl[,7]^2)
  }

  clamfl <- cbind(ID, C14, calBP, error, res, d)
  clamfl[is.na(clamfl)] <- ""
  clamdir <- paste0(clamdir, "/", core)
  if(!dir.exists(clamdir))
    dir.create(clamdir)
  fastwrite(clamfl, paste0(clamdir, "/", core, ".csv"), sep=sep, row.names=FALSE, quote=FALSE)
}



#' @name Bacon.cleanup
#' @title Remove files made to produce the current core's age-depth model.
#' @description Remove files ending in .bacon, .plum (if it exists), .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @details If cores behave badly, you can try cleaning up previous runs and settings, by
#' removing files *.bacon, *.plum, *.out, *.pdf, *_ages.txt, and *_settings.txt of current core.
#' @return A message stating that the files and settings of this run have been deleted.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   Bacon.cleanup()
#' @export
Bacon.cleanup <- function(set=get('info')) {
  files <- c(paste0(set$prefix, ".bacon"), paste0(set$prefix, ".plum"), paste0(set$prefix, ".out"),
    paste0(set$prefix, ".pdf"), paste0(set$prefix, "_ages.txt"),
    paste0(set$coredir,set$core, "/", set$core, "_settings.txt"))
  for(i in files)
    if(file.exists(i))
      tmp <- file.remove(i)
  if(exists("tmp"))
    rm(tmp)
#  if(exists('info')) # new Oct 2020
#    rm(info)
  message("Previous runs of core ", set$core, " with thick=", set$thick, " ", set$depth.unit, " deleted. Now try running the core again\n")
}



# If coredir is left empty, check for a folder named Cores in the current working directory, and if this doesn't exist, for a folder called Bacon_runs (make this folder if it doesn't exist yet and if the user agrees).
# Check if we have write access. If not, tell the user to provide a different, writeable location for coredir.
assign_coredir <- function(coredir, core, ask=TRUE, isPlum=FALSE) {
  ifelse(isPlum, runs <- "Plum_runs", runs <- "Bacon_runs")
  if(coredir == "") {
    if(dir.exists("Cores"))
      coredir <- "Cores" else
        if(dir.exists(runs))
          coredir <- runs else {
            coredir <- runs
            if(!ask)
              ans <- "y" else
                ans <- readline(paste0("I will create a folder called ", coredir, ", is that OK? (y/n)  "))

              if(tolower(substr(ans,1,1)) == "y")
                wdir <- dir.create(coredir, FALSE) else
                  stop("No problem. Please provide an alternative folder location using coredir\n", call.=FALSE)
            if(!wdir)
              stop("cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.", call.=FALSE)
        }
  } else {
    if(!dir.exists(coredir))
        wdir <- dir.create(coredir, FALSE)
      if(!dir.exists(coredir)) # if it still doesn't exist, we probably don't have enough permissions
        stop("cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.", call.=FALSE)
  }
  coredir <- validateDirectoryName(coredir)
  cat("The run's files will be put in this folder: ", coredir, core, "\n", sep="")
  return(coredir)
}



# read the dets file, converting old formats to new ones if so required
read.dets <- function(core, coredir, othername=c(), set=get('info'), sep=",", dec=".", cc=1) {
  # if a .csv file exists, read it (checking that it is more recent than any .dat file in the folder). Otherwise, read the .dat file, check the columns, report back if >4 (>5?) columns, and convert to .csv (report this also)
  if(length(othername)> 0) {
    csv.file <- paste0(coredir, core, "/", othername)
    dat.file <- csv.file
  } else {
      csv.file <- paste0(coredir, core, "/", core, ".csv")
      dat.file <- paste0(coredir, core, "/", core, ".dat")
    }

  dR.names <- c("r", "d", "d.r", "dr", "deltar", "r.mn", "rm", "rmn", "res.mean", "res.mn", "delta.r")
  dSTD.names <- c("d.std", "std", "std.1", "dstd", "r.std", "rstd", "res.sd", "delta.std", "deltastd")
  ta.names <- c("t", "t.a", "ta", "sta")
  tb.names <- c("t", "t.b", "tb", "stb")
  cc.names <- c("c", "cc")
  suggested.names <- c("labID", "age", "error", "depth", "cc", "dR", "dSTD", "ta", "tb")
  changed <- 0

  if(nchar(sep) != 1L) # fwrite wants sep to be of 1 character length only
    message("Warning, sep should be 1 character only, please adapt")

  if(file.exists(csv.file)) {
    dets <- fastread(csv.file, header=TRUE, sep=sep)
    if(file.exists(dat.file)) # deal with old .dat files
      if(file.mtime(csv.file) < file.mtime(dat.file))
        message("Warning, the .dat file is newer than the .csv file! I will read the .csv file. From now on please modify ", csv.file, ", not ", dat.file) else
          message("Reading", csv.file)
    } else {
      if(file.exists(paste0(csv.file, ".txt"))) {
        file.rename(paste0(csv.file, ".txt"), csv.file)
        message("Removing .txt extension from .csv file")
      } else {
        message("No .csv file found, reading", dat.file, " and converting it to .csv")
        dets <- fastread(dat.file, header=TRUE)
        changed <- 1
        }
    }
  name <- tolower(names(dets))
  commas <- grep(",,", readLines(csv.file)) # check if there are too many commas (e.g., lines with just commas)
  if(length(!is.na(commas)) > 0) # often an artefact of spreadsheet programs
    stop("check the .csv file in a plain-text editor for 'orphan' commas\n", call.=FALSE)

  # check if 'classic' dets file, which has a different column order from the current default
  if(ncol(dets) > 4)
    if(ncol(dets) == 5) { # then probably a 'new' dets file
      ok <- 0
      # if((name[5] %in% cc.names) && (min(dets[,5])[1] >= 0) && (max(dets[,5])[1] <= 4)) {} else
      if(name[5] %in% cc.names) ok <- ok+1
      if(min(dets[,5])[1] >= 0) ok <- ok+1
      if(max(dets[,5])[1] <= 4) ok <- ok+1
      if(ok < 3)
        stop("unexpected name or values in fifth column (cc, should be between 0 and 4). Please check the manual for guidelines in producing a correct .csv file.\n", call.=FALSE)
    } else
      if(ncol(dets) == 6) { # probably an 'old' file: dR, dSTD, but could also be cc and delta.R (so no column for delta.STD)
        ok <- 0
        if(name[5] %in% dR.names) ok <- ok+1
        if(name[6] %in% dSTD.names) ok <- ok+1
        if(ok == 2) {
          dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
          message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file)
          changed <- 1
        } else
          stop("unexpected names for columns 5/6. If you want to include delta.R, also add a column for delta.STD. Check the manual for guidelines to producing a correct .csv file.\n", call.=FALSE)
      } else
        if(ncol(dets) == 7) { # probably a 'new' file: cc, dR, dSTD
          ok <- 0
          if(name[5] %in% cc.names) ok <- ok+1
          if(min(dets[,5])[1] >= 0) ok <- ok+1
          if(max(dets[,5])[1] <= 4) ok <- ok+1
          if(name[6] %in% dR.names) ok <- ok+1
          if(name[7] %in% dSTD.names) ok <- ok+1
          if(ok < 5)
            stop("unexpected column names, order or values in dets file. \nPlease check the manual for correct dets file formats.\n", call.=FALSE)
        } else
          if(ncol(dets) == 8) { # probably an 'old' file: dR, dSTD, ta, tb
            ok <- 0
            if(name[5] %in% dR.names) ok <- ok+1
            if(name[6] %in% dSTD.names) ok <- ok+1
            if(name[7] %in% ta.names) ok <- ok+1
            if(name[8] %in% tb.names) ok <- ok+1
            #if(all(range(dets[,8] - dets[,7]) == c(1,1)) ok <- ok+1
            if(all(dets[,8] == dets[,7]+1)) ok <- ok + 1
            if(ok == 5) { # check that these set expected t distribution values
              dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
              message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file)
              changed <- 1
            } else
              stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file", call.=FALSE)
          } else
            if(ncol(dets) == 9) { # most complex case, many checks needed
              ok <- 0
              if(name[9] %in% cc.names) ok <- ok+1 # almost sure this is a 'classic' dets file
              if(min(dets[,9])[1] >= 0) ok <- ok+1
              if(max(dets[,9])[1] <= 4) ok <- ok+1
              tab <- range(dets[,8] - dets[,7])
              if(tab[1] == 1) ok <- ok+1
              if(tab[2] == 1) ok <- ok+1
              if(name[5] %in% dR.names) ok <- ok+1
              if(name[6] %in% dSTD.names) ok <- ok+1
              if(name[7] %in% ta.names) ok <- ok+1
              if(name[8] %in% tb.names) ok <- ok+1
              if(ok == 9) { # column names as expected?
                dets <- dets[,c(1:4,9,5:8)] # shuffle columns around
                message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD, t.a, t.b, cc. \nAdapting column order and saving as", csv.file)
                changed <- 1
              } else { # probably a 'new' file from more recent Bacon
                ok <- 0
                if(name[5] %in% cc.names) ok <- ok+1
                if(min(dets[,5])[1] >= 0) ok <- ok+1
                if(max(dets[,5])[1] <= 4) ok <- ok+1
                tab <- range(dets[,9] - dets[,8])
                if(tab[1] == 1) ok <- ok+1
                if(tab[2] == 1) ok <- ok+1
                if(name[8] %in% ta.names) ok <- ok+1
                if(name[9] %in% tb.names) ok <- ok+1
                if(name[6] %in% dR.names) ok <- ok+1
                if(name[7] %in% dSTD.names) ok <- ok+1
                if(ok < 9)
                  stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file", call.=FALSE)
               }
            } else
              stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.\n", call.=FALSE)

  # more sanity checks
  ok <- 0
  if(is.numeric(dets[,2])) ok <- ok+1
  if(is.numeric(dets[,3])) ok <- ok+1
  if(is.numeric(dets[,4])) ok <- ok+1
  if(ok < 3)
    stop("unexpected values in dets file, I expected numbers. Check the manual.\n", call.=FALSE)
  if(min(dets[,3]) <= 0) {
    message("Warning, zero year errors don't exist in Bacon's world. I will increase them to 1 ", set$age.unit, " yr")
    dets[dets[,3] <= 0,3] <- 1
    changed <- 1
  }
  if(min(0, diff(dets[,4])) < 0) { # added 0 (for if just 1 row of dates)
    stop("The depths in your .csv file are not in ascending order. Please correct this.")
    # dets <- dets[ order(dets[,4]), ] #CHANGED: se elimina "set" antes de dets, por un error en uso del objeto
    # changed <- 1
  }
  if(any(grepl(",", dets[,1])))  {# commas in the lab IDs, which can confuse bacon (c++)
    dets[,1] <- gsub(",", "", dets[,1])
    message("Warning, I removed commas from the lab ID fields")
  }
  if(any(grepl('"', dets[,1]))) {# commas in the lab IDs, which can confuse bacon (c++)
    dets[,1] <- gsub('"', '', dets[,1])
    message("Warning, I removed quotation marks")
  }

  # if current dets differ from original .csv file, rewrite it
  if(changed > 0)
    fastwrite(as.data.frame(dets), csv.file, sep=sep, dec=dec, row.names=FALSE, col.names=suggested.names[1:ncol(dets)], quote=FALSE) # sep was paste0(sep, "\t") but fwrite needs sep to be exactly 1 char
  dets
}



# read in default values, values from previous run, any specified values, and report the desired one. Internal function.
Bacon.settings <- function(core, coredir, dets, thick, remember=TRUE, d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, hiatus.shape, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob, defaults, runname, ssize, dark, youngest.age, oldest.age, cutoff, age.res, after, age.unit) {

  vals <- list(d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob, age.unit) # do these now need the length for each parameter where available?
  valnames <- c("d.min", "d.max", "d.by", "depths.file", "slump", "acc.mean", "acc.shape", "mem.mean", "mem.strength", "boundary", "hiatus.depths", "hiatus.max", "BCAD", "cc", "postbomb", "cc1", "cc2", "cc3", "cc4", "depth.unit", "normal", "t.a", "t.b", "delta.R", "delta.STD", "prob", "age.unit")

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
                if(i>2) cat(" using previous run's value for ", valnames[i], ", ", ext.pre, "\n", sep="")
                ext <- ext.pre
              } else {
                  if(i==13) ifelse(ext.def, "using BC/AD", "using cal BP") else
                  if(i>2) cat(" using default value for ", valnames[i], ", ", ext.def, "\n", sep="")
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
  prevfile <- paste0(coredir, core, "/", core, "_settings.txt")
  prevf <- FALSE
  if(file.exists(prevfile)) {
    prevfile <- readLines(prevfile, n=-1)
    if(length(prevfile) > 0) prevf <- TRUE
  }

  d.min <- extr(1); d.max <- extr(2); d.by <- extr(3); depths.file <- extr(4)
  slump <- extr(5); acc.mean <- extr(6);
  if(length(acc.shape) == 1)
    acc.shape <- extr(7)
  mem.mean <- extr(8)
  mem.strength <- extr(9)
  boundary <- if(is.na(boundary)[1]) NA else sort(extr(10))
  hiatus.depths <- if(is.na(hiatus.depths)[1]) NA else sort(extr(11))
  hiatus.max <- extr(12)
  BCAD <- extr(13); cc <- extr(14); postbomb <- extr(15); cc1 <- extr(16, isnum=FALSE)
  cc2 <- extr(17, isnum=FALSE); cc3 <- extr(18, isnum=FALSE); cc4 <- extr(19, isnum=FALSE)
  depth.unit <- extr(20, isnum=FALSE); normal <- extr(21); t.a <- extr(22); t.b <- extr(23)
  delta.R <- extr(24); delta.STD <- extr(25); prob <- extr(26); age.unit <- extr(27, isnum=FALSE)

  if(is.na(d.min)) # removed || d.min == "NA" 10 April 2020
    d.min <- min(dets[,4])
  if(is.na(d.max)) # removed || d.max == "NA" 10 April 2020
    d.max <- max(dets[,4])
  if(length(acc.shape) < length(acc.mean))
    acc.shape <- rep(acc.shape, length(acc.mean)) else
      if(length(acc.shape) > length(acc.mean))
        acc.mean <- rep(acc.mean, length(acc.shape))
  if(length(mem.strength) < length(mem.mean))
    mem.strength <- rep(mem.strength, length(mem.mean)) else
      if(length(mem.strength) > length(mem.mean))
        mem.mean <- rep(mem.mean, length(mem.strength))

  ### produce/update settings file, and return the values
  prevfile <- file(paste0(coredir, core, "/", core, "_settings.txt"), "w")
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
  for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "") # redundant
  cat(BCAD, " #BCAD\n", cc, " #cc\n", postbomb, " #postbomb\n",
    cc1, " #cc1\n", cc2, " #cc2\n", cc3, " #cc3\n", cc4, " #cc4\n",
    depth.unit, " #depth.unit\n", normal, " #normal\n", t.a, " #t.a\n", t.b, " #t.b\n",
    delta.R, " #delta.R\n", delta.STD, " #d.STD\n", prob, " #prob\n", age.unit, "#age.unit\n", sep="", file=prevfile)
  close(prevfile)

  if(length(youngest.age) == 0)
    youngest.age <- min(1950 - as.integer(format(Sys.time(), "%Y")))#, round(dets[,2] - (5*dets[,3])))
  if(length(oldest.age) == 0)
    oldest.age <- max(1e6, round(dets[,2] + (5*dets[,3])))

  list(core=core, thick=thick, dets=dets, d.min=d.min, d.max=d.max,
    d.by=d.by, depths.file=depths.file, slump=slump,
    acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean,
    mem.strength=mem.strength, boundary=boundary,
    hiatus.depths=hiatus.depths, hiatus.max=hiatus.max,
    BCAD=BCAD, cc=cc, postbomb=postbomb,
    cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=noquote(depth.unit), unit=depth.unit, age.unit=noquote(age.unit), normal=normal,
    t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, date=date(),
    runname=runname, ssize=ssize, dark=dark, youngest.age=youngest.age, oldest.age=oldest.age,
    cutoff=cutoff, age.res=age.res, after=after)
}



# write files to be read by the main Bacon age-depth modelling function
write.Bacon.file <- function(set=get('info'), younger.than=c(), older.than=c(), save.info=set$save.info) {

  if(length(set$slump) > 0) {
    dets <- set$slumpdets
    hiatus.depths <- set$slumphiatus
    boundary <- set$slumpboundary
  } else {
      dets <- set$dets
      hiatus.depths <- set$hiatus.depths
      boundary <- set$boundary
    }
  if(set$d.min < min(dets[,4])) { # repeat relevant row, change error and depth
    # extrap <- c(NA, min(dets[,2]), max(1e5, 1e3*dets[,2], 1e3*dets[,3]), set$d.min, 0)
    dets <- rbind(dets[which(dets[,4] == min(dets[,4]))[1],], dets, make.row.names=FALSE)
    dets[1,1] <- NA # calling this "d.min" causes issues
    dets[1,3] <- max(1e5, 1e3*dets[,2], 1e3*dets[,3])
    dets[1,4] <- set$d.min
  }
  if(set$d.max > max(dets[,4])) { # repeat relevant row, change error and depth
    # extrap <- c(NA, max(dets[,2]), max(1e5, 1e3*dets[,2], 1e3*dets[,3]), set$d.max, 0)
    dets <- rbind(dets, dets[which(dets[,4] == max(dets[,4]))[1],], make.row.names=FALSE)
    dets[nrow(dets),1] <- NA # calling this "d.max" causes issues
    dets[nrow(dets),3] <- max(1e5, 1e3*dets[,2], 1e3*dets[,3])
    dets[nrow(dets),4] <- set$d.max
  }

  fl <- file(set$bacon.file, "w")
  cat("## Ran on", set$date, "\n\n", file=fl)
  cat("Cal 0 : ConstCal;\nCal 1 : ",
  if(set$cc1=="IntCal20" || set$cc1=="\"IntCal20\"") "IntCal20"
    else noquote(set$cc1), ", ", set$postbomb, ";\nCal 2 : ",
  if(set$cc2=="Marine20" || set$cc2=="\"Marine20\"") "Marine20"
    else noquote(set$cc2), ";\nCal 3 : ",
  if(set$cc3=="SHCal20" || set$cc3=="\"SHCal20\"") "SHCal20"
    else noquote(set$cc3), ", ", set$postbomb, ";",
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") set$cc4 <- c()
    else
      paste0("\nCal 4 : GenericCal, ", set$cc4, ";", sep=""), file=fl)
  cat("\n\n##   id.   age    std   depth  delta.R  delta.STD     t.a   t.b   cc", file=fl)

  if(ncol(dets) == 4) { # then we need to provide some constants once only
    cat("\nDet 0 : ", as.character(dets[1,1]), " ,  ", dets[1,2], ",  ",
      dets[1,3], ",  ", dets[1,4], ",  ", set$delta.R, ",  ", set$delta.STD,
      ",  ", set$t.a, ",  ", set$t.b, ",  ", set$cc, ";", sep="", file=fl)
      if(nrow(dets)>1)
        for(i in 2:nrow(dets)) {
          if(i %in% older.than) # is an "older-than" date
            cat("\nDetCensor ", i-1, " : ",  as.character(dets[i,1]),
            " , ", dets[i,2], ", ", dets[i,3], ", ", dets[i,4],
            ";", sep="", file=fl) else
              if(i %in% younger.than)
                cat("\nDetCensorE ", i-1, " : ",  as.character(dets[i,1]),
                " , ", dets[i,2], ", ", dets[i,3], ", ", dets[i,4],
                ";", sep="", file=fl) else
                cat("\nDet ", i-1, " : ",  as.character(dets[i,1]),
                " , ", dets[i,2], ", ", dets[i,3], ", ", dets[i,4],
                ";", sep="", file=fl)
        }
  } else { # use additional columns provided within the dates file
      cc <- dets[,5]
      delta.R <- rep(set$delta.R, nrow(dets))
      delta.R[cc==0] <- 0 # only apply dR to C14 dates
      delta.STD <- rep(set$delta.STD, nrow(dets))
      delta.STD[cc==0] <- 0 # only apply dR to C14 dates
      t.a <- rep(set$t.a, nrow(dets))
      t.b <- rep(set$t.b, nrow(dets))

      if(ncol(dets) >= 7) {
        delta.R <- dets[,6]
        delta.STD <- dets[,7]
      }
      if(ncol(dets) >= 9) {
        t.a <- dets[,8]
        t.b <- dets[,9]
      }

      for(i in 1:nrow(dets)) {
        firstentry <- "\nDet "
        if(i %in% older.than)
          firstentry <- "\nDetCensor "
        if(i %in% younger.than)
          firstentry <- "\nDetCensorE "  
        cat(firstentry, i-1, " : ",  as.character(dets[i,1]), " , ",
          dets[i,2], ", ", dets[i,3], ", ", dets[i,4],  ",  ",
          delta.R[i], ",  ", delta.STD[i], ",  ", t.a[i], ",  ", t.b[i], ",  ",
          cc[i], ";", sep="", file=fl) 
      }
    }

  if(!is.na(hiatus.depths[1])) {
    if(is.null(boundary[1]))
      message("  Hiatus set at depth(s)", paste("", hiatus.depths)) else
        message("  Boundary set at depth(s) ",  paste("", boundary))
    if(length(set$acc.shape) == 1)
      set$acc.shape <- rep(set$acc.shape, length(hiatus.depths)+1)
    if(length(set$acc.mean) == 1)
      set$acc.mean <- rep(set$acc.mean, length(hiatus.depths)+1)
    if(length(set$hiatus.max) == 1)
      set$hiatus.max <- rep(set$hiatus.max, length(hiatus.depths))
#      if(length(set$hiatus.shape)==1)
#        set$hiatus.shape <- rep(set$hiatus.shape, length(set$hiatus.depths))
    if(save.info)
      assign_to_global("info", set)

    cat("\n\n### Depths and priors for fixed hiatuses, in descending order",
      "\n##### cm  alpha beta      ha     hb", file=fl)
    for(i in length(hiatus.depths):1)
      cat("\nHiatus ", i-1, ":  ", hiatus.depths[i], ",  ", set$acc.shape[i+1],
        ",  ", set$acc.shape[i+1]/set$acc.mean[i+1], ",  ", .9, # last value (h.a) was NA but this conflicts with setting initial values for hiatus length # .9 was 0.1, should be <1
        ",  ", set$hiatus.max[i], ";", sep="", file=fl)
  }

  cK <- set$d.min+(set$thick*set$K)
  ### final parameters - dmax now calculated as dmin+(dC*K)
  if(is.na(set$seed)) {
    wrapup <- paste("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax",
      "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
      ",  ", set$youngest.age, ",  ", set$oldest.age, ",  ", set$th0[1], ",  ", set$th0[2],
      ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
      ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
      ", ", cK, ";\n", sep="")
    cat(wrapup, file=fl)
  } else {
    wrapup <- paste("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax seed",
      "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
      ",  ", set$youngest.age, ",  ", set$oldest.age, ",  ", set$th0[1], ",  ", set$th0[2],
      ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
      ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
      ", ", cK, ", ", as.integer(set$seed), ";\n", sep="")
    cat(wrapup, file=fl)
  }
  close(fl)
}



# internal functions to speed up reading and writing files, using the data.table R package if present
fastread <- function(fnam, ...)
  if("data.table" %in% (.packages()))
    as.data.frame(fread(fnam, ...)) else
      read.table(fnam, ...)



fastwrite <- function(out, fnam, ...)
  if("data.table" %in% (.packages()))
    fwrite(as.data.frame(out), fnam, ...) else
      write.table(out, fnam, ...)



# function to read output files into memory
Bacon.AnaOut <- function(fnam, set=get('info'), MCMC.resample=TRUE) {
  out <- fastread(fnam) # was read.table
  if(MCMC.resample)
    if(set$ssize < nrow(out)) { # MB Aug 2022
      ss <- (nrow(out) - set$ssize + 1):nrow(out) # select the last ssize its only
      out <- out[ss,] # MB Aug 2022
      fastwrite(out, fnam, col.names=FALSE, row.names=FALSE) # MB Dec 2022
    }
  n <- ncol(out)-1
  set$n <- n
  set$Tr <- nrow(out)
  set$Us <- out[,n+1]
  set$output <- out[,1:n]
  return(set)
}



# function to read plum output files into memory, updated May 2021
Plum.AnaOut <- function(fnam, set=get('info'), MCMC.resample=TRUE) {
  out <- fastread(fnam) # was read.table
  if(MCMC.resample)
    if(set$ssize < nrow(out)) { # MB Aug 2022
      ss <- (nrow(out) - set$ssize + 1):nrow(out) # select the last ssize its only
      out <- out[ss,] # MB Aug 2022
      fastwrite(out, fnam, col.names=FALSE, row.names=FALSE) # MB Dec 2022
    }
  n <- ncol(out)-1
  set$nPs <- n
  set$TrPs <- nrow(out)
  set$phi  <- out[,1]
  set$ps   <- out[,2:(n+1)]
  return(set)
}



# function to load results in global environment
# parameter position defaults to 1, which equals an assignment to the global environment
assign_to_global <- function(key, val, pos=1) {
  assign(key, val, envir=as.environment(pos) )
}



#' @name ageranges
#' @title save age ranges of depth(s) 
#' @description Calculates an age summary (min/max probability range, median and mean) of a depth or a series of depths. The depths are returned and can be saved to a file.  
#' @return The minimum and maximum of the probability range (default 2.5% of both sides), median and mean ages for a depth or series of depths.
#' @param d Depth or depths for which age summaries are to be calculated. If left empty, a sequence of depths is calculated from the top to the bottom of the core, see 'd.by'.
#' @param file The name of the file where the age summary has to be stored. If empty (default), the ages are returned.
#' @param sep Separator for the fields, if saving to a file (defaults to a tab, "\\t").
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param na.rm Whether or not NAs are to be removed. Defaults to \code{na.rm=FALSE}. NAs will still be returned for depths outside of the core's depth range.
#' @param prob Probability range. Half of the range is taken away from both sides of the distribution (e.g., 2.5\% for the default of \code{prob=0.95}).
#' @param d.by Steps for calculation of depths, if `d` is left empty. Defaults to steps of 1. 
#' @param roundby Rounding for the age estimates. Defaults to 1 decimal.
#' @param show.progress Show a progress bar. Defaults to TRUE.
#' @export
ageranges <- function(d=c(), file=c(), sep="\t", set=get("info"), BCAD=set$BCAD, na.rm=FALSE, prob=0.95, d.by=1, roundby=1, show.progress=TRUE) {
  if(length(d) == 0)
    d <- seq(set$d.min, set$d.max, by=d.by)
  
  n <- length(d)
  mins <- numeric(length(d))
  maxs <- numeric(length(d))
  medians <- numeric(length(d))
  means <- numeric(length(d))

  if(show.progress)
    if(n > 50) 
      pb <- txtProgressBar(min=0, max=max(1, length(d)-1), style=3)
	   
  for(i in seq_len(n)) {
	if(show.progress)
      if(n >= 50)
        if(n < 500 || i %% 10 == 0) 
          setTxtProgressBar(pb, i)
	  
    ages <- Bacon.Age.d(d[i], set, BCAD=BCAD, na.rm=na.rm)
    quan <- quantile(ages, c((1-prob)/2, 1-((1-prob)/2), .5), na.rm=na.rm)
    mins[i] <- quan[1]
    maxs[i] <- quan[2]
    medians[i] <- quan[3]
    means[i] <- mean(ages) 	
  }
  if(show.progress && n>50)
    message("")
  
  mins <- round(mins, roundby)
  maxs <- round(maxs, roundby)
  medians <- round(medians, roundby)
  means <- round(means, roundby)
  
  if(length(d) == 1) {
    names <- c(paste0("min.", 100*prob), paste0("max.", 100 * prob),
      "median", "mean")	  
	summ <- setNames(c(mins, maxs, medians, means), names)
  } else {
      summ <- data.frame(depth = d, mins, maxs, medians, means)
      colnames(summ) <- c("depth",
        paste0("min.", 100 * prob), paste0("max.", 100*prob),
        "median", "mean")
    }
  if(length(file) > 0)
    fastwrite(summ, file, col.names=TRUE, row.names=FALSE, sep=sep, quote=FALSE)
  return(summ)	  
}
