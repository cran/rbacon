#' rbacon
#'
#' Bacon produces Bayesian age-depth models from dated deposits, reconstructing Bayesian 
#' accumulation histories through combining radiocarbon and other dates with prior information (Blaauw and Christen, 2011).
#'
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk> J. Andres Christen <jac@cimat.mx> 
#' @importFrom grDevices dev.cur dev.off pdf dev.copy2pdf grey rgb dev.list extendrange
#' @importFrom graphics abline box curve hist image layout legend lines par plot points polygon segments rect axis mtext
#' @importFrom stats approx dbeta density dgamma dnorm dunif lm quantile rnorm weighted.mean coef
#' @importFrom utils read.csv read.table write.table packageName txtProgressBar setTxtProgressBar
#' @importFrom Rcpp evalCpp
#' @importFrom coda gelman.diag mcmc.list as.mcmc

#' @useDynLib rbacon
#' @name rbacon
NULL

# do: ensure that agedepth uses info$depth.unit and same for age, adapt lowest section so that there is always a section below the lowermost dated depth, check fs::path(dir, data_name) as cross-platform alternative to specifying paths, can ssize be predicted more accurately?,  why do we warn that "acc.shape cannot be equal to acc.mean"?

# done: 

# for future versions: check updated src code (which has the rplum code) to see where the models get too fat, add vignette(s). produce greyscale proxy graph with proxy uncertainties?, smooth bacon, check/adapt behaviour of AgesOfEvents around hiatuses, add function to estimate best thickness, F14C, if hiatus or boundary plot acc.posts of the individual sections?, allow for asymmetric cal BP errors (e.g. read from files), make more consistent use of dark for all functions (incl. flux and accrate.age.ghost), remove darkest?, introduce write.Bacon function to write files only once user agrees with the model, proxy.ghost very slow with long/detailed cores - optimization possible?, check again if/how/when Bacon gets confused by Windows usernames with non-ascii characters (works fine on Mac)

#' @name Bacon
#' @title Main age-depth modelling function
#' @description This is the main age-depth modelling function of the rbacon package.
#' @details Bacon is an approach to age-depth modelling that uses Bayesian statistics in order to reconstruct Bayesian
#' accumulation histories for deposits, through combining radiocarbon and other dates with prior information ('Blaauw' and 'Christen', 2011).
#'
#' Bacon divides a core into many thin vertical sections (by default of \code{thick=5} cm thickness),
#' and through millions of Markov Chain Monte Carlo (MCMC) iterations estimates
#' the accumulation rate (in years/cm; so more correctly, sedimentation times) for each of these sections.
#' Combined with an estimated starting date for the first section, these accumulation rates then form the age-depth model.
#' The accumulation rates are constrained by prior information on the accumulation rate (\code{acc.mean, acc.shape)} and its
#' variability between neighbouring depths, or "memory" (\code{mem.mean, mem.strength}). Hiatuses can be introduced as well, also constrained by prior information (\code{hiatus.max}).
#'
#' Although Bacon works with any kind of absolute dates (e.g., OSL, tephra or other dates on a calendar scale),
#' it is often used to age-model 14C-dated sequences. Radiocarbon dates should be calibrated using either IntCal20
#' (for terrestrial northern hemisphere material; Reimer et al., 2020), Marine20 (for marine dates; Hughen et al., 2020),
#' SHCal20 (for southern hemisphere dates; Hogg et al., 2020) or any other calibration curve (see below), while modern 14C
#' dates are calibrated using one of the post-bomb calibration curves (NH1, NH2 or NH3 for the northern hemisphere,
#' SH1-2 or SH3 for the southern hemisphere; Hua et al., 2013). See \url{http://calib.org/CALIBomb/} if you are unsure which
#' postbomb curve you need. If Bacon finds postbomb dates (negative 14C ages) and you haven't specified a postbomb curve,
#' you will be prompted. Provide postbomb curves as, e.g., \code{postbomb=1} for the NH1 postbomb curve (2 for NH2, 3 for NH3, 4 for SH1-2, 5 for SH3).
#'
#' For calendar dates, i.e. dates that are already on the calendar scale and thus should not be calibrated, set\code{cc=0}.
#'
#' @param core Name of the core, given using quotes. Defaults to one of the cores provided with rbacon, \code{core="MSB2K"}.
#' An alternative core provided with this package is RLGH3 (Jones et al., 1989).
#' To run your own core, produce a .csv file with the dates as outlined in the manual, add a folder with the core's name to the default directory for cores (see \code{coredir}), and save the .csv file there. For example, the file's location and name could be \code{Bacon_runs/MyCore/MyCore.csv}. Then run Bacon as follows: \code{Bacon("MyCore")}
#' @param thick Bacon will divide the core into sections of equal thickness specified by thick (default \code{thick=5}).
#' @param coredir Folder where the core's files \code{core} are and/or will be located. This will be a folder with the core's name, within either the folder \code{coredir='Bacon_runs/'}, or the folder Cores/ if it already exists within R's working directory, or a custom-built folder. For example, use \code{coredir="."} to place the core's folder within the current working directory, or \code{coredir="F:"} if you want to put the core's folder and files on a USB drive loaded under F:.
#' Thinner (and thus more) sections will result in smoother age-models, but too many sections can cause `run-away' models.
#' @param prob Confidence interval to report. This should lie between 0 and 1, default 0.95 (95 \%).
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth).
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param d.by Depth intervals at which ages are calculated. Defaults to \code{d.by=1}.
#' @param seed Seed used for C++ executions, if it is not assigned then the seed is set by system.
#' @param depth.unit Units of the depths. Defaults to \code{depth.unit="cm"}.
#' @param age.unit Units of the ages. Defaults to \code{age.unit="yr"}.
#' @param unit Deprecated and replaced by \code{depth.unit}.
#' @param depths By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.
#' Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}.
#' @param depths.file By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.
#' If \code{depths.file=TRUE}, Bacon will read a file containing the depths for which you require ages.
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir},
#' and its name should start with the core's name and end with `_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param acc.shape The prior for the accumulation rate consists of a gamma distribution with two parameters.
#' Its shape is set by acc.shape (default \code{acc.shape=1.5}; higher values result in more peaked shapes).
#' @param acc.mean The accumulation rate prior consists of a gamma distribution with two parameters. Its mean is set by acc.mean (default \code{acc.mean=20} yr/cm (or whatever age or depth units are chosen),
#' which can be changed to, e.g., 5, 10 or 50 for different kinds of deposits). Multiple values can be given in case of hiatuses or boundaries, e.g., Bacon(hiatus.depths=23, acc.mean=c(5,20))
#' @param mem.strength The prior for the memory (dependence of accumulation rate between neighbouring depths) is a beta distribution, which looks much like the gamma distribution.
#'  but its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.strength=4}
#'  (higher values result in more peaked shapes) allow for a large range of posterior memory values.
#' @param mem.mean The prior for the memory is a beta distribution, which looks much like the gamma distribution but
#' its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.mean=0.7}
#' allow for a large range of posterior memory values.
#' @param boundary The assumed depths of any boundary, which divides sections of different accumulation rate regimes (e.g., as indicated by major change in the stratigraphy). No hiatus is assumed between these sections, and memory is reset crossing the boundary. Different accumulation priors can be set for the sections above and below the boundary, e.g., \code{acc.mean=c(5, 20)}. See also \code{hiatus.depths}, \code{mem.mean}, \code{acc.mean} and \code{acc.shape}. Setting many boundaries might not work, and having more than one boundary per model section (see \code{'thick'}) might not work either.
#' @param hiatus.depths The assumed depths for any hiatus should be provided as, e.g.,
#' \code{hiatus.depths=20} for one at 20cm depth, and \code{hiatus.depths=c(20,40)} for two hiatuses at 20 and 40 cm depth.
#' @param hiatus.max The prior for the maximum length of the hiatus. Hiatus length is a uniform distribution, with equal probabilities between 0 and \code{hiatus.max} yr (or whatever other \code{age.unit} is chosen).
#' @param add Add a value to the maximum hiatus length if a boundary is chosen. Defaults to 100 yr (or whatever other age unit is chosen). Can be adapted if Bacon complains that the parameters are out of support.
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param cc Calibration curve for C-14 dates: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine),
#' \code{cc=3} for SHCal20 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param cc1 For northern hemisphere terrestrial 14C dates (IntCal20).
#' @param cc2 For marine 14C dates (Marine20).
#' @param cc3 For southern hemisphere 14C dates (SHCal20).
#' @param cc4 Use an alternative curve (3 columns: cal BP, 14C age, error, separated by white spaces and saved as a plain-text file). See \code{ccdir}.
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{ccdir=""} since they are loaded into R's memory.
#' For example, use \code{ccdir="."} to choose current working directory, or \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}. Note that all calibration curves should reside in the same directory. If you want to add a custom-built curve, put it in the directory where the default calibration curves are (probably \code{list.files(paste0(.libPaths(), "/rbacon/extdata/Curves/"))}).
#' Alternatively produce a new folder, and add your curve as well as the default calibration curves there (cc1, cc2 and cc3; e.g., \code{write.table(copyCalibrationCurve(1), "./3Col_intcal20.14C", sep="\t")}.)
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}
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
#' @param normal By default, Bacon uses the student's t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param suggest If initial analysis of the data indicates abnormally slow or fast accumulation rates, Bacon will suggest to change the prior.
#'  Also, if the length of the core would cause too few or too many sections with the default settings, Bacon will suggest an alternative section thickness \code{thick}.
#'  Accept these suggested alternative settings by typing "y" (or "yes please" if you prefer to be polite), or leave as is by typing "n" (or anything else, really). To get rid of these suggestions, use \code{suggest=FALSE}.
#' @param reswarn Bacon will warn you if the number of sections lies outside the safe range (default between 10 and 200 sections;
#' \code{reswarn=c(10,200)}). Too few sections could lead to an `elbowy' model while with too many sections the modelling process can get lost,
#'  resulting in age-models far away from the dated depths.
#' @param remember Bacon will try to remember which settings you have applied to your cores (default \code{remember=TRUE}). If you run into inconsistencies or other problems,
#' try running your core again with \code{remember=FALSE}, or, start cleanly by typing \code{Bacon.cleanup()}.
#' @param ask By default Bacon will ask you to confirm that you want to run the core with the provided settings. Disable this using \code{ask=FALSE} (e.g., for batch runs).
#' @param run In order to load an existing Bacon run instead of producing a new one, you can use \code{run=FALSE}.
#' @param defaults Name of the file containing settings for the core. For internal use only - do not change.
#' @param sep Separator between the fields of the plain text file containing the dating information. Default \code{sep=","}.
#' @param dec Character for decimal points. Default to \code{dec="."}.
#' @param runname Text to add to the corename for specific runs, e.g., \code{runname="MyCore_Test1"}.
#' @param slump Upper and lower depths of any sections of assumed abrupt accumulation, that require excising before age-modelling (and adding after age-modelling). Requires pairs of depths, e.g., \code{slump=c(10,15,60,67)} for slumps at 67-60 and 15-10 cm core depth.
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param ssize The approximate amount of iterations to store at the end of the MCMC run. Default 2000; decrease for faster (but less reliable) runs or increase for cores where the MCMC mixing (panel at upper-left corner of age-model graph) appears problematic.
#' @param th0 Starting years for the MCMC iterations.
#' @param burnin Amount of initial, likely sub-optimal MCMC iterations that will be removed.
#' @param MinAge Minimum age limit for Bacon runs, default at current year in cal BP. To set plot limits, use \code{yr.min} instead.
#' @param MaxAge Maximum age limit for Bacon runs, default at 1,000,000 cal BP. To set plot limits, use \code{yr.max} instead.
#' @param MinYr Deprecated - use MinAge instead.
#' @param MaxYr Deprecated - use MaxAge instead.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
#' @param plot.pdf Produce a pdf file of the age-depth plot. Defaults to \code{plot.pdf=TRUE} after a Bacon run.
#' @param dark Darkness of the greyscale age-depth model. The darkest grey value is \code{dark=1} by default.
#' Lower values will result in lighter grey but values >1 are not allowed.
#' @param date.res Date distributions are plotted using \code{date.res=100} segments by default.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param yr.res Deprecated - use age.res instead
#' @param close.connections Internal option to close connections after a run. Default \code{close.connections=TRUE}.
#' @param verbose Provide feedback on what is happening (default \code{verbose=TRUE}).
#' @param ... options for the age-depth graph. See \link{agedepth} and \link{calib.plot}
#' @author Maarten Blaauw, J. Andres Christen
#' @return An age-depth model graph, its age estimates, and a summary.
#' @examples
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())
#' }
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   Bacon(cc=2, delta.R=80, delta.STD=40, coredir=tempfile())
#' }
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474.
#' \url{https://projecteuclid.org/euclid.ba/1339616472}
#'
#' Christen, J.A., Perez E., S., 2010. A new robust statistical model for radiocarbon data. Radiocarbon 51, 1047-1059. 
#'
#' Reimer et al., 2020. The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0–55 cal kBP). Radiocarbon 62. doi: 10.1017/RDC.2020.41
#'
#' Hogg et al. 2020 SHCal20 Southern Hemisphere calibration, 0-55,000 years cal BP. Radiocarbon 62. doi: 10.1017/RDC.2020.59
#'
#' Hughen et al. 2020 Marine20-the marine radiocarbon age calibration curve (0-55,000 cal BP). Radiocarbon 62. doi: 10.1017/RDC.2020.68.
#'
#' Hua, Q., Barbetti, M., Rakowski, A.Z., 2013. Atmospheric radiocarbon for the period 1950-2010.
#' Radiocarbon 55(4), <doi:10.2458/azu_js_rc.v55i2.16177>.
#'
#' Jones, V.J., Stevenson, A.C., Battarbee, R.W., 1989. Acidification of lakes in Galloway, south west Scotland
#' - a diatom and pollen study of the post-glacial history of the Round Loch of Glenhead.
#' Journal of Ecology 77: 1-23.
#'
#' @export
Bacon <- function(core="MSB2K", thick=5, coredir="", prob=0.95, d.min=NA, d.max=NA, d.by=1, seed=NA, depths.file=FALSE, depths=c(), depth.unit="cm", age.unit="yr", unit=depth.unit, acc.shape=1.5, acc.mean=20, mem.strength=4, mem.mean=0.7, boundary=NA, hiatus.depths=NA, hiatus.max=10000, add=c(), after=.0001/thick, cc=1, cc1="IntCal20", cc2="Marine20", cc3="SHCal20", cc4="ConstCal", ccdir="", postbomb=0, delta.R=0, delta.STD=0, t.a=3, t.b=4, normal=FALSE, suggest=TRUE, reswarn=c(10,200), remember=TRUE, ask=TRUE, run=TRUE, defaults="defaultBacon_settings.txt", sep=",", dec=".", runname="", slump=c(), BCAD=FALSE, ssize=2000, th0=c(), burnin=min(500, ssize), MinAge=c(), MaxAge=c(), MinYr=MinAge, MaxYr=MaxAge, cutoff=.001, plot.pdf=TRUE, dark=1, date.res=100, age.res=200, yr.res=age.res, close.connections=TRUE, verbose=TRUE, ...) {
  # Check coredir and if required, copy example file in core directory
  coredir <- assign_coredir(coredir, core, ask)
  if(core == "MSB2K" || core == "RLGH3") {
    dir.create(paste(coredir, core, "/", sep=""), showWarnings = FALSE, recursive = TRUE)
    fileCopy <- system.file(paste("extdata/Cores/", core, sep=""), package="rbacon") # change to package rbacon
    file.copy(fileCopy, coredir, recursive = TRUE, overwrite=FALSE)
  }

  # set the calibration curve
  if(ccdir == "")
    ccdir <- paste(system.file("extdata", package=packageName()), "/Curves/", sep="")
  ccdir <- .validateDirectoryName(ccdir)

  # default_settings.txt is located within system.file
  defaults <- system.file("extdata", defaults, package=packageName())
  # read in the data, adapt settings from defaults if needed
  dets <- .read.dets(core, coredir, sep=sep, dec=dec, cc=cc)
  # give feedback about calibration curves used
  if(ncol(dets) > 4 && length(cc) > 0) {
    cc.csv <- unique(dets[,5])
	if(verbose) {
      if(length(cc.csv) == 1) {
        if(cc.csv != cc)
          message(" Using calibration curve specified within the .csv file,", cc[cc.csv], "\n")
      } else
        if(min(cc.csv) == 0)
          message(" Using a mix of cal BP and calibrated C-14 dates\n")
        else
          message(" Using several C-14 calibration curves\n")
    }
  }

  if(suggest) { # adapt prior for mean accumulation rate?
    sugg <- sapply(c(1,2,5), function(x) x*10^(-1:2)) # some suggested "round" values
    ballpacc <- lm(dets[,2]*1.1 ~ dets[,4])$coefficients[2] # very rough acc.rate estimates, uncalibrated dates
    ballpacc <- abs(sugg - ballpacc) # get absolute differences between given acc.mean and suggested ones
    ballpacc <- ballpacc[ballpacc > 0] # do not suggest 0
    sugg <- sugg[order(ballpacc)[1]] # suggest rounded acc.rate with lowest absolute difference
    if(!sugg %in% acc.mean) {
      ans <- readline(message(" Ballpark estimates suggest changing the prior for acc.mean to ", sugg, " ", age.unit, "/", depth.unit, ". OK? (y/N) "))
      if(tolower(substr(ans,1,1)) == "y")
        acc.mean <- sugg else
          message(" No problem, using the provided prior")
    }
  }

  if(!is.na(boundary[1]))
    boundary <- sort(unique(boundary))
  if(!is.na(hiatus.depths[1])) {
    hiatus.depths <- sort(unique(hiatus.depths))
    if(length(acc.mean) == 1)
      acc.mean <- rep(acc.mean, length(hiatus.depths)+1)
  }

  info <- .Bacon.settings(core=core, coredir=coredir, dets=dets, thick=thick, remember=remember, d.min=d.min, d.max=d.max, d.by=d.by, depths.file=depths.file, slump=slump, acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean, mem.strength=mem.strength, boundary=boundary, hiatus.depths=hiatus.depths, hiatus.max=hiatus.max, BCAD=BCAD, cc=cc, postbomb=postbomb, cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=depth.unit, normal=normal, t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, defaults=defaults, runname=runname, ssize=ssize, dark=dark, MinAge=MinAge, MaxAge=MaxAge, cutoff=cutoff, age.res=age.res, after=after, age.unit=age.unit)
  .assign_to_global("info", info)
  info$coredir <- coredir
  info$seed <- seed
  info$isplum <- FALSE

  ### check for initial mistakes
  if(any(info$acc.shape == info$acc.mean))
    stop("acc.shape cannot be equal to acc.mean", call.=FALSE)
  if(info$t.b - info$t.a != 1)
    stop("t.b - t.a should always be 1, check the manual", call.=FALSE)
  if(min(acc.shape) < 1)
    warning("\nWarning, using values <1 for acc.shape might cause unexpected results\n", call.=TRUE)

  ### calibrate dates
  if(info$cc > 0) # confirm we are using radiocarbon dates
    if(info$postbomb == 0 && ((ncol(info$dets) == 4 && min(info$dets[,2]) < 0) ||
      ncol(info$dets)>4 && max(info$dets[,5]) > 0 && min(info$dets[info$dets[,5] > 0,2]) < 0))
        stop("you have negative C14 ages so should select a postbomb curve", call.=FALSE)
  info$calib <- .bacon.calib(dets, info, date.res, ccdir=ccdir)

  ### find some relevant values
  info$rng <- c()
  for(i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[which(tmp[,2]>cutoff),1])
  }
  if(length(th0)==0) # provide two ball-park/initial age estimates
    info$th0 <- round(rnorm(2, max(MinAge, dets[1,2]), dets[1,3]))
  info$th0[info$th0 < info$MinAge] <- info$MinAge # otherwise twalk will not start

  ### assign depths
  if(length(depths) == 0)
    depths <- seq(info$d.min, info$d.max, by=d.by) # was info$d.by
  if(depths.file) {
    dfile <- paste0(info$coredir, info$core, "/", info$core, "_depths.txt")
    if(!file.exists(dfile))
      stop("I cannot find the file ", paste0(info$coredir, info$core, "/", info$core, "_depths.txt"), call.=FALSE)
    depths <- read.table(dfile, header=FALSE)[,1]
    if(!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers", call.=FALSE)
  }
  info$depths <- depths
  if(min(depths) < info$d.min)
    info$d.min <- min(depths)
  if(max(depths) > info$d.max)
    info$d.max <- max(depths)

  info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
  info$K <- length(info$elbows)
  info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model

  # stop and warn if hiatus.depths conflict with other parameters
  if(!is.na(info$hiatus.depths[1]) || !is.na(info$boundary[1])) {
    ifelse(is.na(info$boundary[1]), hd <- info$hiatus.depths, hd <- info$boundary)
    if(min(hd) < info$d.min) # hiatus above core top
      stop("cannot have hiatus above the core's top depth. Adapt hiatus.depths or d.min.", call.=FALSE)
    if(max(hd)+info$thick > info$d.max)
      stop("the age-depth model should have at least one section below the one containing the deepest hiatus. Adapt thick or d.max?", call.=FALSE)
    if(length(hd) > 1) { # then check for how far separated hiatuses are
      above <- c()
      for(i in hd)
        above <- c(above, max(which(info$elbows <= i))) # find the section top of each hiatus
        if(any(diff(above) < 2)) # stop if fewer than 2 section elbows separating hiatuses
          stop("we need at least 2 section elbows between hiatuses. Choose fewer hiatuses, different depths, more sections (decrease thick) or a different d.min.\n ", call.=FALSE)
    }
  }

   ans <- "n"
    if(suggest)
      if(length(reswarn) == 2)
        if(info$K < min(reswarn)) {
          sugg <- pretty(thick*(info$K/min(reswarn)), 10)
          sugg <- min(sugg[sugg>0])
          ans <- readline(message(" Warning, the current value for thick, ", thick, ", will result in very few age-model sections (", info$K, ", not very flexible). Suggested maximum value for thick: ", sugg, " OK? (y/n) "))
        } else
          if(info$K > max(reswarn)) {
            sugg <- max(pretty(thick*(info$K/max(reswarn))))
            ans <- readline(message(" Warning, the current value for thick, ", thick, ", will result in very many age-model sections (", info$K, ", possibly hard to run). Suggested minimum value for thick: ", sugg, " OK? (y/n) "))
          }
    if(tolower(substr(ans, 1, 1)) == "y") {
      message(" OK, setting thick to ", sugg, "\n")
      thick <- sugg
      info$thick = thick #CHANGED: if the answer is "yes", the global thick value is not updated
      info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
      if(length(info$slump) > 0) # why here, and not a few lines later?
        info$elbows <- seq(floor(info$d.min), toslump(ceiling(info$d.max), info$slump), by=thick)
      info$K <- length(info$elbows)
      info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model
    }

  ### prepare for any slumps
  if(length(slump) > 0) {
    if(length(slump) %% 2 == 1)
      stop("slumps need both upper and lower depths. Please check the manual", call.=FALSE)
    slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
    info$slump <- slump

    slumpdmax <- toslump(ceiling(info$d.max), slump)
	info$elbows <- seq(floor(info$d.min), slumpdmax, by=thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model

	info$slumpfree <- toslump(depths, slump)
    info$slumphiatus <- toslump(info$hiatus.depths, slump) # check
	if(!is.na(info$boundary[1])) {
      info$slumpboundary <- toslump(info$boundary, slump) # check
      info$slumphiatus <- info$slumpboundary
    }
    slumpdets <- info$dets
    slumpdets[,4] <- toslump(slumpdets[,4], slump, remove=FALSE) # dates within slumps are not removed
    info$slumpdets <- slumpdets[!is.na(slumpdets[,4]),]
  }

  ### produce files
  info$prefix <- paste(coredir, core, "/", core, runname, "_", info$K, sep="")
  info$coredir <- coredir
  info$bacon.file <- paste(info$prefix, ".bacon", sep="")
  if(!file.exists(outfile <- paste(info$prefix, ".out", sep="")))
    file.create(outfile)

  ### store values (again) for future manipulations
  if(BCAD)
    info$BCAD <- TRUE
  if(!is.na(boundary[1])) {
    if(length(slump) > 0)
      boundary <- info$slumpboundary
    info$hiatus.depths <- boundary
    if(length(add) == 0)
      add <- info$acc.mean # then add a short (max)hiatus, large enough not to crash Bacon but not affect the chronology much. Needs more work
    info$hiatus.max <- add
  }
  .assign_to_global("info", info)

  prepare <- function() {
    ### plot initial data and priors
    pn <- c(1,2,3,3)
    if(!is.na(info$hiatus.depths[1])) # was ...hiatus.depths)[1])
      if(is.na(info$boundary[1]))
        pn <- c(1,2,3,4,4,4)
    layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    oldpar <- par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l")
	on.exit(par(oldpar))         	
    .PlotAccPrior(info$acc.shape, info$acc.mean, depth.unit=depth.unit, age.unit=age.unit)
    .PlotMemPrior(info$mem.strength, info$mem.mean, thick)
    if(!is.na(info$hiatus.depths)[1])
      if(is.na(info$boundary)[1])
        .PlotHiatusPrior(info$hiatus.max, info$hiatus.depths)
    calib.plot(info, BCAD=BCAD)
    legend("topleft", core, bty="n", cex=1.5)
  }

  cook <- function() {
    txt <- paste(info$prefix, ".bacon", sep="")
    bacon(txt, outfile, ssize, ccdir)
    scissors(burnin, info)
    agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=TRUE, age.unit=age.unit, depth.unit=depth.unit, ...)
    if(plot.pdf)
      if(interactive())
        if(length(dev.list()) > 0)
          dev.copy2pdf(file=paste0(info$prefix, ".pdf")) else {
            pdf(file=paste0(info$prefix, ".pdf"))
            agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=FALSE, age.unit=age.unit, depth.unit=depth.unit, ...)
            dev.off()
          }
  }

### run bacon if initial graphs seem OK; run automatically, not at all, or only plot the age-depth model
  .write.Bacon.file(info)
  if(!run)
    prepare() else
      if(!ask)
        cook() else {
          prepare()
          ans <- readline(message(" Run ", core, " with ", info$K, " sections? (Y/n) "))
          ans <- tolower(substr(ans,1,1))[1]
          if(ans=="y" || ans=="")
            cook() else
              message("  OK. Please adapt settings")
        }
  if(close.connections)
    closeAllConnections()
}
