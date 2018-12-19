#' rbacon
#' 
#' Bacon produces Bayesian age-depth models from dated deposits, reconstructing Bayesian 
#' accumulation histories through combining radiocarbon and other dates with prior information (Blaauw and Christen, 2011).
#' 
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk> J. Andres Christen <jac@cimat.mx>
#' @importFrom grDevices dev.cur dev.off pdf dev.copy2pdf grey rgb dev.list
#' @importFrom graphics abline box curve hist image layout legend lines par plot points polygon segments rect
#' @importFrom stats approx dbeta density dgamma dnorm dunif lm quantile rnorm weighted.mean
#' @importFrom utils read.csv read.table write.table packageName txtProgressBar setTxtProgressBar 
#' @importFrom Rcpp evalCpp
#' @importFrom coda gelman.diag mcmc.list as.mcmc 
#' @useDynLib rbacon
#' @name rbacon
NULL

# done: 

# do: check how to compile on Win-systems the updated twalk.h with the enhanced buffer size (decides when to keep MCMC calculations in memory and when to write them to a .out file) 


# for future versions: check behaviour of AgesOfEvents around hiatuses, add function to estimate best thickness (once published), check w Andres if correction to remove -5.0 cal BP in IntCal13 and SHCal13 curves was done correctly, F14C, if hiatus or boundary plot acc.posts of the individual sections?, allow for asymmetric cal BP errors (e.g. read from files), make more consistent use of dark for all functions (incl. flux and accrate.age.ghost), remove darkest?, introduce write.Bacon function to write files only once user agrees with the model, proxy.ghost very slow with long/detailed cores - optimization possible?, check again if/how Bacon gets confused by Windows usernames with non-ascii characters (works fine on Mac)

#' @name Bacon 
#' @title Main age-depth modelling function
#' @description This is the main age-depth modelling function of the rbacon package.
#' @details Bacon is an approach to age-depth modelling that uses Bayesian statistics in order to reconstruct Bayesian 
#' accumulation histories for deposits, through combining radiocarbon and other dates with prior information (Blaauw and Christen, 2011).
#'
#' Bacon divides a core into many thin vertical sections (by default of \code{thick=5} cm thickness), 
#' and through millions of Markov Chain Monte Carlo (MCMC) iterations estimates 
#' the accumulation rate (in years/cm; so more correctly, sedimentation times) for each of these sections. 
#' Combined with an estimated starting date for the first section, these accumulation rates then form the age-depth model. 
#' The accumulation rates are constrained by prior information on the accumulation rate (\code{acc.mean, acc.shape)} and its 
#' variability between neighbouring depths, or "memory" (\code{mem.mean, mem.strength}). Hiatuses can be introduced as well, also constrained by prior information (\code{hiatus.max}).
#'
#' Although Bacon works with any kind of absolute dates (e.g., OSL, tephra or other dates on a calendar scale), 
#' it is often used to age-model 14C-dated sequences. Radiocarbon dates should be calibrated using either IntCal13 
#' (for terrestrial northern hemisphere material; Reimer et al., 2013), Marine13 (for marine dates; Reimer et al., 2013), 
#' SHCal13 (for southern hemisphere dates; Hogg et al., 2013) or any other calibration curve (see below), while modern 14C 
#' dates are calibrated using one of the post-bomb calibration curves (NH1, NH2 or NH3 for the northern hemisphere, 
#' SH1-2 or SH3 for the southern hemisphere; Hua et al., 2013). See \url{http://calib.org/CALIBomb} if you are unsure which 
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
#' @param unit Units of the depths. Note that the default prior for accumulation rate assumes the default \code{unit="cm"}.
#' @param depths By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.  
#' Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}.
#' @param depths.file By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.  
#' If \code{depths.file=TRUE}, Bacon will read a file containing the depths for which you require ages. 
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir}, 
#' and its name should start with the core's name and end with `_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param acc.shape The prior for the accumulation rate consists of a gamma distribution with two parameters. 
#' Its shape is set by acc.shape (default \code{acc.shape=1.5}; higher values result in more peaked shapes).
#' @param acc.mean The accumulation rate prior consists of a gamma distribution with two parameters. Its mean is set by acc.mean (default \code{acc.mean=20} yr/cm, 
#' which can be changed to, e.g., 5, 10 or 50 for different kinds of deposits).
#' @param mem.strength The prior for the memory (dependence of accumulation rate between neighbouring depths) is a beta distribution, which looks much like the gamma distribution.
#'  but its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.strength=4}
#'  (higher values result in more peaked shapes) allow for a large range of posterior memory values.
#' @param mem.mean The prior for the memory is a beta distribution, which looks much like the gamma distribution but 
#' its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.mean=0.7}
#' allow for a large range of posterior memory values. 
#' @param boundary The assumed depths of any boundary, which divides sections of different accumulation rate regimes (e.g., as indicated by major change in the stratigraphy). No hiatus is assumed between these sections, and memory is reset crossing the boundary. Different accumulation priors can be set for the sections above and below the boundary, e.g., \code{acc.mean=c(5, 20)}. See also \code{hiatus.depths}, \code{mem.mean}, \code{acc.mean} and \code{acc.shape}.
#' @param hiatus.depths The assumed depths for any hiatus should be provided as, e.g., 
#' \code{hiatus.depths=20} for one at 20cm depth, and \code{hiatus.depths=c(20,40)} for two hiatuses at 20 and 40 cm depth.
#' @param hiatus.max The prior for the maximum length of the hiatus. Hiatus length is a uniform distribution, with equal probabilities between 0 and \code{hiatus.max} yr.
#' @param add Add a value to the maximum hiatus length if a boundary is chosen. Defaults to 100 yr. Can be adapted if Bacon complains that the parameters are out of support. 
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param cc Calibration curve for C-14 dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), 
#' \code{cc=3} for SHCal13 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param cc1 For northern hemisphere terrestrial 14C dates (IntCal13).
#' @param cc2 For marine 14C dates (Marine13).
#' @param cc3 For southern hemisphere 14C dates (SHCal13).
#' @param cc4 Use an alternative curve (3 columns: cal BP, 14C age, error, separated by white spaces and saved as a plain-text file). See \code{ccdir}.
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{ccdir=""} since they are loaded into R's memory. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}. Note that all calibration curves should reside in the same directory. If you want to add a custom-built curve, put it in the directory where the default calibration curves are, or alternatively produce a new folder, and add your curve as well as the default calibration curves there (e.g., \code{write.table(copyCalibrationCurve(1), "./3Col_intcal13.14C", sep="\t")}.)
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
#' @param MinYr Minimum age limit for Bacon runs, default at current year in cal BP. To set plot limits, use \code{yr.min} instead.
#' @param MaxYr Maximum age limit for Bacon runs, default at 1,000,000 cal BP. To set plot limits, use \code{yr.max} instead.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}). 
#' @param plot.pdf Produce a pdf file of the age-depth plot. Defaults to \code{plot.pdf=TRUE} after a Bacon run. 
#' @param dark Darkness of the greyscale age-depth model. The darkest grey value is \code{dark=1} by default.
#' Lower values will result in lighter grey but values >1 are not allowed.
#' @param date.res Date distributions are plotted using \code{date.res=100} segments by default.
#' @param yr.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
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
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' 
#' Christen, J.A., Perez E., S., 2010. A new robust statistical model for radiocarbon data. Radiocarbon 51, 1047-1059.
#' 
#' Reimer, P.J., Bard, E., Bayliss, A., Beck, J.W., Blackwell, P.G., Bronk Ramsey, C., Buck, C.E., Edwards,
#' R.L., Friedrich, M., Grootes, P.M., Guilderson, T.P., Haflidason, H., Hajdas, I., Hatte, C., Heaton, T.J.,
#' Hoffmann, D.L., Hogg, A.G., Hughen, K.A., Kaiser, K.F., Kromer, B., Manning, S.W., Niu, M., Reimer, R.W.,
#' Richards, D.A., Scott, M.E., Southon, J.R., Turney, C.S.M., van der Plicht, J., 2013. IntCal13 and
#' Marine13 radiocarbon age calibration curves 0-50,000 yr cal BP. Radiocarbon 55(4), 1869-1887
#'
#' Hogg, A.G., Hua, Q., Blackwell, P.G., Buck, C.E., Guilderson, T.P., Heaton, T.J., Niu, M., Palmer, J.,
#' Reimer, P.J., Reimer, R., Turney, C.S.M., Zimmerman, S.R.H., 2013. ShCal13 Southern Hemisphere
#' calibration, 0-50,000 cal yr BP. Radiocarbon 55(4), doi:10.2458/azu_js_rc.55.16783.
#'
#' Hua, Q., Barbetti, M., Rakowski, A.Z., 2013. Atmospheric radiocarbon for the period 1950-2010. 
#' Radiocarbon 55(4), doi:10.2458/azu_js_rc.v55i2.16177.
#'
#' Jones, V.J., Stevenson, A.C., Battarbee, R.W., 1989. Acidification of lakes in Galloway, south west Scotland
#' - a diatom and pollen study of the post-glacial history of the Round Loch of Glenhead. 
#' Journal of Ecology 77: 1-23.
#'
#' @export
Bacon <- function(core="MSB2K", thick=5, coredir="", prob=0.95, d.min=NA, d.max=NA, d.by=1, depths.file=FALSE, depths=c(), unit="cm", acc.shape=1.5, acc.mean=20, mem.strength=4, mem.mean=0.7, boundary=NA, hiatus.depths=NA, hiatus.max=10000, add=100, after=.000001, cc=1, cc1="IntCal13", cc2="Marine13", cc3="SHCal13", cc4="ConstCal", ccdir="", postbomb=0, delta.R=0, delta.STD=0, t.a=3, t.b=4, normal=FALSE, suggest=TRUE, reswarn=c(10,200), remember=TRUE, ask=TRUE, run=TRUE, defaults="default_settings.txt", sep=",", dec=".", runname="", slump=c(), BCAD=FALSE, ssize=2000, th0=c(), burnin=min(500, ssize), MinYr=c(), MaxYr=c(), cutoff=.001, plot.pdf=TRUE, dark=1, date.res=100, yr.res=200, ...) {
  # Check coredir and if required, copy example file in core directory
  coredir <- assign_coredir(coredir, core, ask)
  if(core == "MSB2K" || core == "RLGH3") {
    dir.create(paste(coredir, core, "/", sep=""), showWarnings = FALSE, recursive = TRUE)
    fileCopy <- system.file(paste("extdata/Cores/", core, sep=""), package="rbacon")
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
    if(length(cc.csv) == 1) {
      if(cc.csv != cc)
        cat(" Using calibration curve specified within the .csv file,", cc[cc.csv], "\n")
      } else
        if(min(cc.csv) == 0)
          cat(" Using a mix of cal BP and calibrated C-14 dates\n") else
            cat(" Using several C-14 calibration curves\n")
    }

    if(suggest) { # adapt prior for mean accumulation rate?
      sugg <- sapply(c(1,2,5), function(x) x*10^(-1:2)) # some suggested "round" values
      ballpacc <- lm(dets[,2]*1.1 ~ dets[,4])$coefficients[2] # very rough acc.rate estimates, uncalibrated dates
      ballpacc <- abs(sugg - ballpacc) # get absolute differences between given acc.mean and suggested ones
      ballpacc <- ballpacc[ballpacc > 0] # do not suggest 0
      sugg <- sugg[order(ballpacc)[1]] # suggest rounded acc.rate with lowest absolute difference
      if(!sugg %in% acc.mean) {
        ans <- readline(cat(" Ballpark estimates suggest changing the prior for acc.mean to ", sugg, " yr/", unit, ". OK? (y/n)  ", sep=""))
        if(tolower(substr(ans,1,1)) == "y")
          acc.mean <- sugg else
            cat(" No problem, using prior acc.mean=", acc.mean, " yr/", unit, "\n", sep="")
      }
    }

  info <- .Bacon.settings(core=core, coredir=coredir, dets=dets, thick=thick, remember=remember, d.min=d.min, d.max=d.max, d.by=d.by, depths.file=depths.file, slump=slump, acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean, mem.strength=mem.strength, boundary=boundary, hiatus.depths=hiatus.depths, hiatus.max=hiatus.max, BCAD=BCAD, cc=cc, postbomb=postbomb, cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, unit=unit, normal=normal, t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, defaults=defaults, runname=runname, ssize=ssize, dark=dark, MinYr=MinYr, MaxYr=MaxYr, cutoff=cutoff, yr.res=yr.res, after=after)
  .assign_to_global("info", info)
  info$coredir <- coredir
 
  ### check for initial mistakes
  if(any(info$acc.shape == info$acc.mean))
    stop("\n Warning! acc.shape cannot be equal to acc.mean", call.=FALSE)  
  if(info$t.b - info$t.a != 1)
    stop("\n Warning! t.b - t.a should always be 1, check the manual")

  ### calibrate dates
  if(info$cc > 0) # confirm we are using radiocarbon dates
    if(info$postbomb == 0 && ((ncol(info$dets) == 4 && min(info$dets[,2]) < 0) ||
      ncol(info$dets)>4 && max(info$dets[,5]) > 0 && min(info$dets[info$dets[,5] > 0,2]) < 0))
        stop("\nWarning, you have negative C14 ages so should select a postbomb curve")
  info$calib <- .bacon.calib(dets, info, date.res, ccdir=ccdir)

  ### find some relevant values
  info$rng <- c()
  for(i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[which(tmp[,2]>cutoff),1])
  }
  if(length(th0)==0) # provide two ball-park/initial age estimates
    info$th0 <- round(rnorm(2, max(MinYr, dets[1,2]), dets[1,3]))
  info$th0[info$th0 < info$MinYr] <- info$MinYr # otherwise twalk will not start

  ### assign depths, possibly suggest alternative value for thick
  if(length(depths) == 0)
    depths <- seq(info$d.min, info$d.max, by=info$d.by) 
  if(depths.file) { 
    dfile <- paste0(info$coredir, info$core, "/", info$core, "_depths.txt")
    if(!file.exists(dfile))
      stop(" Warning! I cannot find the file ", paste0(info$coredir, info$core, "/", info$core, "_depths.txt"), "\n")
    depths <- read.table(dfile, header=FALSE)[,1]
    if(!is.numeric(depths[1]))
      stop(" Warning! File should contain numbers only, no headers\n")
  }
  info$depths <- depths	
  if(min(depths) < info$d.min)
    info$d.min <- min(depths)
  if(max(depths) > info$d.max)
    info$d.max <- max(depths)

  info$d <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
  info$K <- length(info$d)
  ans <- "n"
  if(suggest)
    if(length(reswarn) == 2)
      if(info$K < min(reswarn)) {
        sugg <- pretty(thick*(info$K/min(reswarn)), 10)
        sugg <- min(sugg[sugg>0])
        ans <- readline(cat(" Warning, the current value for thick, ", thick, ", will result in very few age-model sections (", info$K, ", not very flexible). Suggested maximum value for thick: ", sugg, " OK? (y/n) ", sep=""))
      } else
        if(info$K > max(reswarn)) {
          sugg <- max(pretty(thick*(info$K/max(reswarn))))
          ans <- readline(cat(" Warning, the current value for thick, ", thick, ", will result in very many age-model sections (", info$K, ", possibly hard to run). Suggested minimum value for thick: ", sugg, " OK? (y/n) ", sep=""))
        }
  if(tolower(substr(ans, 1, 1)) == "y") {
    cat(" OK, setting thick to ", sugg, "\n")
    thick <- sugg
    info$d <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
    info$K <- length(info$d)
  }

  ### prepare for any slumps
  if(length(slump) > 0) {
    if(length(slump) %% 2 == 1)
      stop("\n Warning, slumps need both upper and lower depths. Please check the manual", call.=FALSE)
    slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
      info$slump <- slump
    info$slumpfree <- excise(depths, slump, info$d.by)
    info$slumpboundary <- excise(info$boundary, slump, info$d.by) # check
    info$slumphiatus <- excise(info$hiatus.depths, slump, info$d.by) # check
    if(!is.na(info$slumpboundary[1]))
      info$slumphiatus <- info$slumpboundary
    slumpdets <- info$dets
    slumpdets[,4] <- excise(slumpdets[,4], slump, info$d.by)
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
  if(!is.na(boundary)[1]) {
    if(length(slump) > 0)
      boundary <- info$slumpboundary  
    info$hiatus.depths <- boundary
    info$hiatus.max <- add  
  }
  .assign_to_global("info", info)  

  prepare <- function() {
    ### plot initial data and priors
    pn <- c(1,2,3,3)
    if(!is.na(info$hiatus.depths)[1])
      if(is.na(info$boundary)[1])
        pn <- c(1,2,3,4,4,4)
    layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l")
    .PlotAccPrior(info$acc.shape, info$acc.mean)
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
    agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, talk=TRUE, ...)
    if(plot.pdf) 
      if(interactive())
        if(length(dev.list()) > 0 ) # then try to save plotting time 
          dev.copy2pdf(file=paste0(info$prefix, ".pdf")) else {
            pdf(file=paste(info$prefix, ".pdf", sep=""))
            agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, talk=FALSE, ...)
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
          ans <- readline(cat("  Run", core, "with", info$K, "sections? (y/n) "))
          if(tolower(substr(ans,1,1))[1]=="y")
            cook() else cat("  OK. Please adapt settings.\n\n")
          }
  # closeAllConnections()
}

