#' rbacon
#' 
#' Bacon is an approach to age-depth modelling that uses Bayesian statistics to reconstruct Bayesian 
#' accumulation histories for deposits, through combining radiocarbon and other dates with prior information Blaauw and Christen, 2011).
#' 
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk> J. Andres Christen <jac@cimat.mx>
#' @importFrom grDevices dev.cur dev.off pdf dev.copy2pdf grey rgb
#' @importFrom graphics abline box curve hist image layout legend lines par plot points polygon segments
#' @importFrom stats approx dbeta density dgamma dnorm lm quantile rnorm weighted.mean
#' @importFrom utils read.csv read.table write.table packageName
#' @importFrom Rcpp evalCpp
#' @importFrom coda gelman.diag mcmc.list as.mcmc 
#' @useDynLib rbacon
#' @name rbacon
NULL  

# for future versions: slump, F14C, 0-yr hiatuses w certain hiatus parameter combinations, enhanced age calculations around hiatuses, if hiatus plot acc.posts of the individual sections?, allow for asymmetric cal BP errors (e.g. read from files), check if postbomb dates really taken into account (if MinYr=-1e3), check more consistent use of dark and darkest for all functions (incl. flux and accrate.age.ghost)

# Done version 2.3.1.2: ensured more predictable behaviour if R is started in a non-writable directory (e.g. default R on Windows). Now uses tempfile() correctly in Bacon.hist. Added confidence ranges to accrate.age.ghost and flux.age.ghost. Calculates mean ages better (based on all its and not on hists, which are sensitive to choice of bins)

# Done version 2.3.1.1: Official R package, new options: depths (in addition to depths.file), default core directory now called 'Bacon_runs' (if the Cores folder doesn't already exist within the directoy where R is working), added option to not plot x or y axis (xaxt, yaxt), added option to not plot the date distributions mirrored, renamed weighted means of age estimates to means (means of histograms), updated documentation, renamed function flux.age, plot.accrate.age and plot.accrate.depth to flux.age.ghost, accrate.age.ghost and accrate.depth.ghost, respectively. Bugs squashed: BCAD did not cause correct plots, many sundry bugs

# Done version 2.2: changed .hpd to _ages.txt since many users get tricked by the extension, Bacon.hist gives 95\% ranges, mid and wmean, 
# BCAD (sort of...), settings file, removed calc.every (gave problems with long cores), warn if any dateF errors 0 change to 1 optional depths.file for age calculations, language check cpp files, killed hist bug that assumed integers for res,
# changed to .csv files as default with option to use .dat files, suggest to adapt prior for acc.mean if initial estimate
# much different from default (20), allowing for different separator (e.g. French/Scandinavian use ';' not ',', 
# also check use decimal points), renamed 'res' to hopefully more intuitive 'thick', d.R/d.STD 
# for mixed dates where cc=0, enhanced flexibility MaxYr/MinYr (e.g. now depends on dates if cc=0), 
# updated bacon's hist function (now bin/hist2), now reads depths from a file, added option to change axes 
# order in age-depth graphs, added cleanup function to remove prior files etc., updated calibration curves to IntCal13, 
# option in agedepth to only plot the age-model (so not the upper panels)


#################### workhorse functions ####################

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
#' variability between neighbouring depths, or 'memory' (\code{mem.mean, mem.strength}). Hiatuses can be introduced as well, also constrained by prior information (\code{hiatus.mean, hiatus.strength}).
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
#' @param coredir Subfolder where the core's files \code{core} are and will be located. This will be a folder with the core's name, within either the folder \code{coredir='Bacon_runs/'}, or the folder Cores/ if it already exists within R's working directory, or a custom-built folder. For example, use \code{coredir="."} to place the core's folder within the current working directory, or \code{coredir="F:"} if you want to put the core's folder and files on a USB drive loaded under F:.
#' Thinner (and thus more) sections will result in smoother age-models, but too many sections can cause 'run-away' models.
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
#' and its name should start with the core's name and end with '_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param acc.shape The prior for the accumulation rate consists of a gamma distribution with two parameters. 
#' Its shape is set by acc.shape (default \code{acc.shape=1.5}; higher values result in more peaked shapes).
#' @param acc.mean The accumulation rate prior consists of a gamma distribution with two parameters. Its mean is set by acc.mean (default \code{acc.mean=20} yr/cm, 
#' which can be changed to, e.g., 5, 10 or 50 for different kinds of deposits).
#' @param mem.strength The prior for the memory is a beta distribution, which looks much like the gamma distribution
#'  but its values are always between 0 (no assumed memory) and 1 (100 \% memory). Its default settings of \code{mem.strength=4}
#'  (higher values result in more peaked shapes) allow for a large range of posterior memory values.
#' @param mem.mean The prior for the memory is a beta distribution, which looks much like the gamma distribution but 
#' its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.mean=0.7}
#' allow for a large range of posterior memory values. 
#' @param hiatus.depths The assumed depths for any hiatus should be provided as, e.g., 
#' \code{hiatus.depths=20} for one at 20cm depth, and \code{hiatus.depths=c(20,40)} for two hiatuses at 20 and 40 cm depth.
#' @param hiatus.shape The prior for the length of the hiatus, which is a gamma distribution with two parameters. 
#' Its shape is set by hiatus.shape (default \code{acc.shape=1}).
#' @param hiatus.mean The prior for the length of the hiatus, which is a gamma distribution with two parameters. 
#' Its mean is set by hiatus.mean (default \code{acc.mean=1000}).
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param cc Calibration curve for C-14 dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), 
#' \code{cc=3} for SHCal13 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param cc1 For northern hemisphere terrestrial 14C dates (IntCal13).
#' @param cc2 For marine 14C dates (Marine13).
#' @param cc3 For southern hemisphere 14C dates (SHCal13).
#' @param cc4 Use an alternative curve (3 columns: cal BP, 14C age, error, separated by white spaces and saved as a plain-text file).
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{ccdir=""} since they are loaded into R's memory. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}
#' @param deltaR Mean of core-wide age offsets (e.g., regional marine offsets).
#' @param deltaSTD Error of core-wide age offsets (e.g., regional marine offsets). 
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
#'  Accept these suggested alternative settings by typing 'y' (or 'yes please' if you prefer to be polite), or leave as is by typing 'n' (or anything else, really). To get rid of these suggestions, use \code{suggest=FALSE}. 
#' @param reswarn Bacon will warn you if the number of sections lies outside the safe range (default between 10 and 200 sections; 
#' \code{reswarn=c(10,200)}). Too few sections could lead to an 'elbowy' model while with too many sections the modelling process can get lost,
#'  resulting in age-models far away from the dated depths. 
#' @param remember Bacon will try to remember which settings you've applied to your cores (default \code{remember=TRUE}). If you run into inconsistencies or other problems, 
#' try running your core again with \code{remember=FALSE}, or, start cleanly by typing \code{Bacon.cleanup()}. 
#' @param ask By default Bacon will ask you to confirm that you want to run the core with the provided settings. Disable this using \code{ask=FALSE} (e.g., for batch runs).
#' @param run In order to load an existing Bacon run instead of producing a new one, you can use \code{run=FALSE}.
#' @param defaults Name of the file containing settings for the core. For internal use only - do not change. 
#' @param sep Separator between the fields of the plain text file containing the dating information. Default \code{sep=","}.
#' @param dec Character for decimal points. Default to \code{dec="."}.
#' @param runname Text to add to the corename for specific runs, e.g., \code{runname="MyCore_Test1"}.
#' @param slump (Not implemented yet).
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param ssize The approximate amount of iterations to store at the end of the MCMC run. Default 2000; decrease for faster (but less reliable) runs or increase for cores where the MCMC mixing (panel at upper-left corner of age-model graph) appears problematic.
#' @param th0 Starting years for the MCMC iterations.
#' @param burnin Amount of initial, likely sub-optimal MCMC iterations that will be removed.
#' @param MinYr Minimum age limit for Bacon runs, default at -1000 cal BP. To set plot limits, use \code{yr.min} instead.
#' @param MaxYr Maximum age limit for Bacon runs, default at 1,000,000 cal BP. To set plot limits, use \code{yr.max} instead.
#' @param find.round Temporary internal measure to deal with rounding problems at hiatuses; default 4. 
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
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   Bacon(cc=2, deltaR=80, deltaSTD=40, coredir=tempfile())
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
Bacon <- function(core="MSB2K", thick=5, coredir="", prob=0.95, d.min=NA, d.max=NA, d.by=1, depths.file=FALSE, depths=c(), unit="cm", acc.shape=1.5, acc.mean=20, mem.strength=4, mem.mean=0.7, hiatus.depths=NA, hiatus.shape=1, hiatus.mean=1000, after=.0001, cc=1, cc1="IntCal13", cc2="Marine13", cc3="SHCal13", cc4="ConstCal", ccdir="", postbomb=0, deltaR=0, deltaSTD=0, t.a=3, t.b=4, normal=FALSE, suggest=TRUE, reswarn=c(10,200), remember=TRUE, ask=TRUE, run=TRUE, defaults="default_settings.txt", sep=",", dec=".", runname="", slump=NA, BCAD=FALSE, ssize=2000, th0=c(), burnin=min(200, ssize), MinYr=c(), MaxYr=c(), find.round=4, cutoff=.001, plot.pdf=TRUE, dark=1, date.res=100, yr.res=200, ...) {
  # If coredir is left empty, check for a folder named Cores in the current working directory, and if this doesn't exist, for a folder called Bacon_runs (make this folder if it doesn't exist yet).
  # Check if we have write access. If not, tell the user to provide a different, writeable location for coredir. 
  if(coredir == "") {
    if(dir.exists("Cores")) 
      coredir <- "Cores" else
        if(dir.exists("Bacon_runs"))
	      coredir <- "Bacon_runs" else {
			coredir <- "Bacon_runs"			  
            wdir <- dir.create(coredir, FALSE)
			if(!wdir)
			  stop("Cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.")
	      }    
    } else {
	    if(!dir.exists(coredir))
          wdir <- dir.create(coredir, FALSE)
        if(!dir.exists(coredir)) # if it still doesn't exist, we probably don't have enough permissions
          stop("Cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.")
	}
  coredir <- .validateDirectoryName(coredir)
  cat("The run's files will be put in this folder: ", coredir, core, "\n", sep="")
	  
  # Copy example file in core directory
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

    if(suggest) {
      # adapt prior for mean accumulation rate?
      sugg <- sapply(c(1,2,5), function(x) x*10^(-1:2)) # some suggested 'round' values
      ballpacc <- lm(dets[,2]*1.1 ~ dets[,4])$coefficients[2] # very rough acc.rate estimates, uncalibrated dates
      ballpacc <- abs(sugg - ballpacc) # get absolute differences between given acc.mean and suggested ones
      ballpacc <- ballpacc[ballpacc > 0] # don't suggest 0; tmp MB
      sugg <- sugg[order(ballpacc)[1]] # suggest rounded acc.rate with lowest absolute difference
      if(sugg != acc.mean) {
        ans <- readline(cat(" Ballpark estimates suggest changing the prior for acc.mean to ", sugg, " yr/", unit, ". OK? (y/n)  ", sep=""))
        if(tolower(substr(ans,1,1)) == "y")
          acc.mean <- sugg else
            cat(" No problem, using prior acc.mean=", acc.mean, " yr/", unit, "\n", sep="")
      }
    }

  info <- .Bacon.settings(core=core, coredir=coredir, dets=dets, thick=thick, remember=remember, d.min=d.min, d.max=d.max, d.by=d.by, depths.file=depths.file, slump=slump, acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean, mem.strength=mem.strength, hiatus.depths=hiatus.depths, hiatus.mean=hiatus.mean, hiatus.shape=hiatus.shape, BCAD=BCAD, cc=cc, postbomb=postbomb, cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, unit=unit, normal=normal, t.a=t.a, t.b=t.b, deltaR=deltaR, deltaSTD=deltaSTD, prob=prob, defaults=defaults, runname=runname, ssize=ssize, dark=dark, MinYr=MinYr, MaxYr=MaxYr, cutoff=cutoff, yr.res=yr.res, after=after, find.round=find.round)
  .assign_to_global("info", info)
  info$coredir <- coredir
 
  ### check for initial mistakes
  if(any(info$acc.shape == info$acc.mean))
    stop("\n Warning! acc.shape cannot be equal to acc.mean", call.=FALSE)
  if(info$t.b - info$t.a != 1)
    stop("\n Warning! t.b - t.a should always be 1, check the manual")

  ### calibrate dates
  if(info$cc > 0) # confirm we're using radiocarbon dates
    if(info$postbomb == 0 && ((ncol(info$dets)==4 && min(info$dets[,2]) < 0) ||
      ncol(info$dets)>4 && max(info$dets[,5]) > 0 && min(info$dets[info$dets[,5] > 0,2]) < 0))
        stop("\nWarning, you have negative C14 ages so should select a postbomb curve")
  info$calib <- .bacon.calib(dets, info, date.res, ccdir=ccdir)

  ### find some relevant values
  info$rng <- c()
  for(i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[which(tmp[,2]>cutoff),1])
  }
  if(length(th0)==0) # provide two ball-park initial age estimates
    info$th0 <- round(rnorm(2, max(MinYr, dets[1,2]), dets[1,3]))
  info$th0[info$th0 < info$MinYr] <- info$MinYr # otherwise twalk will not start

  ### assign depths, possibly suggest alternative value for thick
  info$d <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
  info$K <- length(info$d)
  ans <- "n"
  if(suggest)
    if(length(reswarn)==2)
      if(info$K < min(reswarn)) {
        sugg <- min(pretty(thick*(info$K/min(reswarn))))
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

  ### produce files
  info$prefix <- paste(coredir, core, "/", core, runname, "_", info$K, sep="")
  info$coredir <- coredir
  
  info$bacon.file <- paste(info$prefix, ".bacon", sep="")
  if(!file.exists(outfile <- paste(info$prefix, ".out", sep="")))
    file.create(outfile)

  ### store values (again) for future manipulations
  .assign_to_global("info", info)

  prepare <- function() {
    ### plot initial data and priors
    layout(matrix(if(is.na(info$hiatus.depths)[1]) c(1,2,3,3) else c(1,2,3,4,4,4),
      nrow=2, byrow=TRUE), heights=c(.3,.7))
    par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0))
    .PlotAccPrior(info$acc.shape, info$acc.mean)
    .PlotMemPrior(info$mem.strength, info$mem.mean, thick)
    if(!is.na(info$hiatus.depths)[1])
      .PlotHiatusPrior(info$hiatus.shape, info$hiatus.mean, info$hiatus.depths)
    calib.plot(info, ...)
    legend("topleft", core, bty="n", cex=1.5)
  }
   
  cook <- function() {
    txt <- paste(info$prefix, ".bacon", sep="")
    bacon(txt, outfile, ssize, ccdir)
    scissors(burnin, info)
    agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, ...)
    if(plot.pdf) {
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
  closeAllConnections()
}



#' @name agedepth 
#' @title Plot an age-depth model
#' @description Plot the age-depth model of a core.
#' @details After loading a previous run, or after running either the \link{scissors} or \link{thinner} command, plot the age-model
#' again using the command \code{agedepth()}. 
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param d.lab The labels for the depth axis. Default \code{d.lab="Depth (cm)"}. See also \code{unit}.
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth).
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param d.by Depth intervals at which ages are calculated. Default 1. Alternative depth intervals can be provided using, e.g., d.\code{by=0.5}.
#' @param depths By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}. Alternative depths can be provided as, e.g., \code{depths=seq(0, 100, length=500)} or as a file, e.g., \code{depths=read.table("CoreDepths.txt"}. See also \code{depths.file}.
#' @param depths.file By default, Bacon will calculate the ages for the depths \code{d.min} to \code{d.max} in steps of \code{d.by}.  
#' If \code{depths.file=TRUE}, Bacon will read a file containing the depths for which you require ages. 
#' This file, containing the depths in a single column without a header, should be stored within \code{coredir}, 
#' and its name should start with the core's name and end with '_depths.txt'. Then specify \code{depths.file=TRUE} (default \code{FALSE}). See also \code{depths}.
#' @param yr.min Minimum calendar age of the age-depth plot.
#' @param yr.max Maximum calendar age of the age-depth plot.
#' @param dark Darkness of the greyscale age-depth model. The darkest grey value is \code{dark=1} by default; lower values will result in lighter grey but values >1 are not allowed.
#' @param prob Confidence interval to report (between 0 and 1, default 0.95 or 95\%).
#' @param rounded Rounding of years. Default is to round to single years.
#' @param d.res Resolution or amount of greyscale pixels to cover the depth scale of the age-model plot. Default \code{d.res=200}.
#' @param yr.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param date.res Date distributions are plotted using \code{date.res=100} points by default.
#' @param grey.res Grey-scale resolution of the age-depth model. Default \code{grey.res=100}.
#' @param rotate.axes By default, the age-depth model is plotted with the depths on the horizontal axis and ages on the vertical axis. This can be changed with \code{rotate.axes=TRUE}.
#' @param rev.yr The direction of the age axis, which can be reversed using \code{rev.yr=TRUE}.
#' @param rev.d The direction of the depth axis, which can be reversed using \code{rev.d=TRUE}.
#' @param unit Depth units, default \code{unit="cm"}.
#' @param maxcalc Number of depths to calculate ages for. If this is more than \code{maxcalc=500}, a warning will be shown that calculations will take time.
#' @param height The maximum heights of the distributions of the dates on the plot. See also \code{normalise.dists}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
#' @param panels Divide the graph panel. Defaults to 1 graph per panel, \code{panels=layout(1)}. To avoid dividing into panels, use \code{panels=c()}.
#' @param plot.range Whether or not to plot the curves showing the confidence ranges of the age-model. Defaults to (\code{plot.range=TRUE}).
#' @param range.col The colour of the curves showing the confidence ranges of the age-model. Defaults to medium grey (\code{range.col=grey(0.5)}).
#' @param range.lty The line type of the curves showing the confidence ranges of the age-model. Defaults to \code{range.lty=12}.
#' @param mn.col The colour of the mean age-depth model: default \code{mn.col="red"}.
#' @param mn.lty The line type of the mean age-depth model. Default \code{mn.lty=12}.
#' @param med.col The colour of the median age-depth model: not drawn by default \code{med.col=NA}.
#' @param med.lty The line type of the median age-depth model. Default \code{med.lty=12}.
#' @param C14.col The colour of the calibrated ranges of the dates. Default is semi-transparent blue: \code{C14.col=rgb(0,0,1,.35)}.
#' @param C14.border The colours of the borders of calibrated 14C dates. Default is semi-transparent dark blue: \code{C14.border=rgb(0, 0, 1, 0.5)}.
#' @param cal.col The colour of the non-14C dates. Default is semi-transparent blue-green: \code{cal.col=rgb(0,.5,.5,.35)}.
#' @param cal.border The colour of the border of non-14C dates in the age-depth plot: default semi-transparent dark blue-green: \code{cal.border=rgb(0,.5,.5,.5)}.
#' @param hiatus.col The colour of the depths of any hiatuses. Default \code{hiatus.col=grey(0.5)}.
#' @param hiatus.lty The line type of the depths of any hiatuses. Default \code{hiatus.lty=12}.
#' @param greyscale The function to produce a grey-scale representation of all age-models. Default \code{greyscale=function(x) grey(1-x)}. 
#' @param bins The amount of bins in histograms. Calculated automatically by default: \code{bins=c()}. See also \code{yr.res}.
#' @param normalise.dists By default, the distributions of more precise dates will cover less time and will thus peak higher than less precise dates. This can be avoided by specifying \code{normalise.dists=FALSE}.
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), \code{cc=3} for SHCal13 (southern hemisphere terrestrial). For dates that are already on the cal BP scale use \code{cc=0}.
#' @param title The title of the age-depth model is plotted on the main panel. By default this is the core's name. To leave empty: \code{title=""}. 
#' @param title.location Location of the title. Default \code{title.location='topleft'}.
#' @param after Sets a short section above and below hiatus.depths within which to calculate ages. For internal calculations - do not change.
#' @param bty Type of box to be drawn around plots (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param mar Plot margins (amount of white space along edges of axes 1-4). Default \code{mar=c(3,3,1,1)}.
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted). Defaults to \code{mgp=c(1.5, .7, .0)}.
#' @param xaxs Extension of x-axis. By default, add some extra white-space at both extremes (\code{xaxs="r"}). See ?par for other options. 
#' @param yaxs Extension of y-axis. By default, add some extra white-space at both extremes (\code{yaxs="r"}). See ?par for other options. 
#' @param xaxt Whether or not to plot the x-axis. Can be used to adapt axes after a plot. See ?par for other options. 
#' @param yaxt Whether or not to plot the y-axis. Can be used to adapt axes after a plot. See ?par for other options. 
#' @param plot.pdf Produce a pdf file of the age-depth plot.
#' @param model.only By default, panels showing the MCMC iterations and the priors and posteriors for accumulation rate and memory are plotted above the main age-depth model panel. This can be avoided by supplying \code{model.only=TRUE}. Note however that this removes relevant information to evaluate the age-depth model, so we do recommend to present age-models together with these upper panels.
#' @param dates.only By default, the age-depth model is plotted on top of the dates. This can be avoided by supplying \code{dates.only=TRUE}.
#' @param cleanup Temporary files used to produce the plot are removed automatically after the plot has been made by default; \code{cleanup=TRUE}. 
#' @param talk Provide a summary of the age ranges after producing the age-depth model graph; default \code{talk=TRUE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of the age-depth model, and estimated ages incl. confidence ranges for each depth.
#' @examples 
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())}
#'   agedepth(yr.res=50)
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   agedepth()
#' }
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
agedepth <- function(set=get('info'), BCAD=FALSE, unit="cm", d.lab=c(), yr.lab=c(), d.min=set$d.min, d.max=set$d.max, d.by=set$d.by, depths=c(), depths.file=FALSE, yr.min=c(), yr.max=c(), bins=c(), dark=set$dark, prob=set$prob, rounded=0, d.res=200, yr.res=200, date.res=100, grey.res=100, rotate.axes=FALSE, rev.yr=FALSE, rev.d=FALSE, maxcalc=500, height=15, mirror=TRUE, up=TRUE, cutoff=.001, plot.range=TRUE, panels=layout(1), range.col=grey(0.5), range.lty="12", mn.col="red", mn.lty="12", med.col=NA, med.lty="12", C14.col=rgb(0,0,1,.35), C14.border=rgb(0,0,1,.5), cal.col=rgb(0,.5,.5,.35), cal.border=rgb(0,.5,.5,.5), hiatus.col=grey(0.5), hiatus.lty="12", greyscale=function(x) grey(1-x), normalise.dists=TRUE, cc=set$cc, title=set$core, title.location="topleft", after=set$after, bty="l", mar=c(3,3,1,1), mgp=c(1.5,.7,.0), xaxs="r", yaxs="r", xaxt="s", yaxt="s", plot.pdf=FALSE, dates.only=FALSE, model.only=FALSE, cleanup=TRUE, talk=TRUE) {
    # Load the output, if it exists
    outp <- paste(set$prefix, ".out", sep="")
    if(file.exists(outp))
      lngth <- length(readLines(outp))
    set <- .Bacon.AnaOut(paste(set$prefix, ".out", sep=""), set)

    par(bty=bty, mar=mar, mgp=mgp, yaxs="i")
    if(model.only) 
      panels else { # layout(1) can mess things up if plotting within an existing panel...
        layout(matrix(if(is.na(set$hiatus.depths)[1]) c(1:3, rep(4,3)) else c(1:4, rep(5,4)), nrow=2, byrow=TRUE), heights=c(.3,.7)) 
        .PlotLogPost(set, 0, set$Tr) # convergence information
        .PlotAccPost(set)
        .PlotMemPost(set, set$core, set$K, "", set$mem.strength, set$mem.mean, ds=1, thick=set$thick)
        if(!is.na(set$hiatus.depths[1]))
          .PlotHiatusPost(set, set$hiatus.shape, set$hiatus.mean)
    } 

    # calculate calendar axis limits
    if(length(BCAD) > 0)
      if(BCAD != set$BCAD) 
        set$BCAD <- BCAD
    modelranges <- c()
    modelranges <- Bacon.hist(range(set$d), set, bins=bins, Plot=FALSE, prob=prob, yr.res=yr.res)
    dates <- set$calib$probs
    dateranges <- c()
    for(i in 1:length(dates))
      dateranges <- range(dateranges, dates[[i]][,1])
    if(set$BCAD)
      dateranges <- 1950 - dateranges	 
    
    if(length(yr.min) == 0) 
      yr.min <- min(modelranges, dateranges)
    if(length(yr.max) == 0)
      yr.max <- max(modelranges, dateranges)
    yr.lim <- c(yr.min, yr.max)
    if(set$BCAD) yr.lim <- rev(yr.lim)
    if(rev.yr) yr.lim <- rev(yr.lim)
    d.lim=c(d.max, d.min)
    if(rev.d)
      dlim <- dlim[2:1]
    
    if(length(d.lab) == 0)
      d.lab <- paste("Depth (", unit, ")", sep="")
        if(length(yr.lab) == 0)
      yr.lab <- ifelse(set$BCAD, "BC/AD", "cal yr BP")

    par(xaxs=xaxs, yaxs=yaxs, bty="n")
    if(rotate.axes)
      plot(0, type="n", ylim=d.lim, xlim=yr.lim, ylab=d.lab, xlab=yr.lab, bty="n", xaxt=xaxt, yaxt=yaxt) else
        plot(0, type="n", xlim=d.lim[2:1], ylim=yr.lim, xlab=d.lab, ylab=yr.lab, bty="n", xaxt=xaxt, yaxt=yaxt)

    if(!dates.only)
      .depth.ghost(set, rotate.axes=rotate.axes, d.res=d.res, yr.res=yr.res, grey.res=grey.res, dark=dark, greyscale=greyscale, d.min=d.min, d.max=d.max, cleanup=cleanup)

    calib.plot(set, rotate.axes=rotate.axes, height=height, mirror=mirror, up=up, date.res=date.res, cutoff=cutoff, C14.col=C14.col, C14.border=C14.border, cal.col=cal.col, cal.border=cal.border, new.plot=FALSE, normalise.dists=normalise.dists)
    legend(title.location, title, bty="n", cex=1.5)
    box(bty=bty)

    # now calculate and plot the ranges and 'best' estimates for each required depth
    if(depths.file) {
      dfile <- paste(set$coredir, set$core, "/", set$core, "_depths.txt", sep="")
      if(!file.exists(dfile))
        stop(" Warning! I cannot find the file ", paste(set$coredir, set$core, "/", set$core, "_depths.txt", sep=""), "\n")
      d <- read.table(dfile, header=FALSE)[,1]
      if(!is.numeric(d[1]))
        stop(" Warning! File should contain numbers only, no headers\n")
      d <- sort(as.numeric(d))
    } else
      if(length(depths) > 0)
        d <- sort(depths) else
          d <- seq(d.min, d.max, by=d.by) 
      
    if(length(d) > maxcalc)
      cat(" Warning, this will take quite some time to calculate. I suggest increasing d.by to, e.g.", 10*set$d.by, "\n")

    for(i in set$hiatus.depths)
      if(i %in% d) d <- sort(c(i-after, d)) else 
        d <- sort(c(i-after, i, d)) ### tmp

    abline(v=set$hiatus.depths, col=hiatus.col, lty=hiatus.lty)

    ranges <- Bacon.hist(d, set, prob=prob, yr.res=yr.res, cleanup=cleanup)

    th <- rbind(1, nrow(ranges))
    if(!is.na(set$hiatus.depths[1])) {
      hi.d <- c()
    for(i in set$hiatus.depths)
      hi.d <- c(hi.d, max(which(d<=i)))
    th <- array(sort(c(1, nrow(ranges), hi.d-1, hi.d)), dim=c(2,length(hi.d)+1))
    }

    if(!dates.only)
      for(i in 1:ncol(th)) {
        h <- th[1,i] : th[2,i]
        if(rotate.axes) {
          lines(ranges[h,1], d[h], col=range.col, lty=range.lty)
          lines(ranges[h,2], d[h], col=range.col, lty=range.lty)
          lines(ranges[h,3], d[h], col=med.col, lty=med.lty) # median
          lines(ranges[h,4], d[h], col=mn.col, lty=mn.lty) # mean
        } else {
            lines(d[h], ranges[h,1], col=range.col, lty=range.lty)
            lines(d[h], ranges[h,2], col=range.col, lty=range.lty)
            lines(d[h], ranges[h,3], col=med.col, lty=med.lty) # median
            lines(d[h], ranges[h,4], col=mn.col, lty=mn.lty) # mean
          }
      }
 
    set$ranges <- cbind(d, round(ranges, rounded))
    colnames(set$ranges) <- c("depth", "min", "max", "median", "wmean")
    
    .assign_to_global("info", set)

    if(plot.pdf)
      if(names(dev.cur()) != "null device")
        dev.copy2pdf(file=paste(set$prefix, ".pdf", sep=""))

    write.table(set$ranges, paste(set$prefix, "_ages.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
    rng <- abs(round(set$ranges[,3]-set$ranges[,2], rounded))
    min.rng <- d[which(rng==min(rng))]
    max.rng <- d[which(rng==max(rng))]
    if(length(min.rng)==1) min.rng <- paste(" yr at", min.rng, noquote(set$unit)) else
      min.rng <- paste(" yr between", min(min.rng), "and", max(min.rng), noquote(set$unit))
    if(length(max.rng)==1) max.rng <- paste(" yr at", max.rng, set$unit) else
      max.rng <- paste(" yr between", min(max.rng), "and", max(max.rng), noquote(set$unit))
    if(talk)
  	  if(!dates.only)
        cat("Mean ", 100*prob, "% confidence ranges ", round(mean(rng), rounded), " yr, min. ",
          min(rng), min.rng, ", max. ", max(rng), max.rng, "\n\n", sep="")
}



#' @name Bacon_runs
#' @title List the folders present in the current core directory.
#' @description Lists all folders located within the core's directory.
#' @details The directory is either "Bacon_runs", "Cores" or a custom-named one. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return A list of folders
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   Bacon_runs()
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @export
Bacon_runs <- function()
  list.files(get('info')$coredir)



#' @name Bacon.cleanup 
#' @title Remove files made to produce the current core's age-depth model.
#' @description Remove files .bacon, .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @details If cores behave badly, you can try cleaning up previous runs and settings, by
#' removing files .bacon, .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @return A message stating that the files and settings of this run have been deleted.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}. 
#' @author Maarten Blaauw, J. Andres Christen
#' @examples 
#'   Bacon(run=FALSE, coredir=tempfile())
#'   Bacon.cleanup()
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @export
Bacon.cleanup <- function(set=get('info')) { 
  files <- c(paste(set$prefix, ".bacon", sep=""), paste(set$prefix, ".out", sep=""),
    paste(set$prefix, ".pdf", sep=""), paste(set$prefix, "_ages.txt", sep=""),
    paste(set$coredir,set$core, "/", set$core, "_settings.txt", sep=""))
  for(i in files)
    if(file.exists(i))
      tmp <- file.remove(i)
  if(exists("tmp")) rm(tmp)
    cat("Previous Bacon runs of core", set$core, "with thick =", set$thick, "deleted. Now try running the core again\n")
}




#################### functions for post-run checks and adaptations ####################

#' @name scissors
#' @title Remove the first n iterations.
#' @description Removes the first n iterations of the MCMC time series, and then updates the output file.
#' @details Bacon will perform millions of MCMC iterations for each age-model run by default, although only a fraction 
#' of these will be stored. In most cases the remaining MCMC iterations will be well mixed (the upper left panel 
#' of the fit of the iterations shows no undesirable features such as trends or sudden systematic drops or rises).
#' If the run has a visible remaining burn-in, scissors can be used to remove them. 
#' To remove, e.g., the first 300 iterations, type \code{scissors(300)}. 
#' This will remove these iterations in the run's output file.
#'
#' @param burnin Number of iterations to remove. If this value is higher than the amount of remaining iterations, 
#' a warning is given and the iterations are not removed.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples 
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   scissors(100)
#'   agedepth(d.res=50)
#' }
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   scissors(100)
#'   agedepth()
#' }
#'
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
scissors <- function(burnin, set=get('info')) {
  output <- read.table(paste(set$prefix, ".out", sep=""))
  if(burnin >= nrow(output))
    stop("\nCannot remove that many iterations, there would be none left!\n")
  output <- output[burnin:nrow(output),]
  write.table(output, paste(set$prefix, ".out", sep=""), col.names=FALSE, row.names=FALSE)
    
  info <- get('info')
  info$output <- output
  .assign_to_global ("info", info)
}



#' @name thinner
#' @title Thin iterations by given proportion, for example if autocorrelation is visible within the MCMC series, and then update the output file.
#' @description Thin iterations by given proportion, for example if autocorrelation is visible within the MCMC series.
#' @details From all iterations, a proportion is removed with to-be-removed iterations sampled randomly among all iterations.
#' @param proportion Proportion of iterations to remove. Should be between 0 and 1. Default \code{proportion=0.1}.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples 
#' \dontshow{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   thinner(.1)
#'   agedepth(d.res=50)
#' }
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   thinner(.2)
#'   agedepth()
#' }
#'
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
thinner <- function(proportion=0.1, set=get('info')) {
  output <- read.table(paste(set$prefix, ".out", sep=""))
  if(proportion >= 1)
    stop("\nCannot remove that many iterations, there would be none left!\n\n", call.=FALSE)
  proportion <- sample(nrow(output), proportion*nrow(output))
  output <- output[-proportion,]
  write.table(output, paste(set$prefix, ".out", sep=""), col.names=FALSE, row.names=FALSE)
    
  info <- get('info')
  info$output <- output
  .assign_to_global ("info", info)
}



#' @name Baconvergence 
#' @title Test to identify poorly mixed MCMC runs.
#' @description Test how well-mixed and converged the MCMC runs are with the chosen core and settings, by running the core several times and comparing the different runs using the Gelman and Rubin Reduction factor (Brooks and Gelman, 1998).
#' @details Generally Bacon will perform millions of MCMC iterations for each age-model run, although only a fraction 
#' of these will be stored. In most cases the remaining MCMC iterations will be well mixed (the upper left panel 
#' of the fit of the iterations shows no strange features such as sudden systematic drops or rises).
#'  However if the iterations seem not well mixed, or if too few remain (say less than a few hundred), 
#'  then you could check the Gelman and Rubin Reduction Factor. Too high differences (high Factors) between runs 
#' indicate poor MCMC mixing. Robust MCMC mixing is indicated by a Gelman and Rubin Reduction factor 
#' (Brooks and Gelman, 1998) below the 1.05 safety threshold.
#' @param core Name of the core, given using quotes. Defaults to one of the cores provided with rbacon, \code{core="MSB2K"}. 
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param runs Amount of runs to test for mixing. Default \code{runs=5}.
#' @param ... additional options that can be given to the Bacon function. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return NA
#' @examples
#'   \donttest{
#'     Baconvergence(runs=2, ssize=100, coredir=tempfile()) # a quick-and-dirty toy example
#'   }
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring
#' convergence of iterative simulations. 
#' _Journal of Computational and Graphical Statistics_, *7*, 434-455.
#' @export
Baconvergence <- function(core="MSB2K", set=get('info'), runs=5, ...) {
  MCMC <- list()
  for(i in 1:runs) {
    cat("run number", i, "...\n")
    Bacon(core=core, ask=FALSE, suggest=FALSE, run=TRUE, ...)
    MCMC[[i]] <- read.table(paste(set$prefix, ".out", sep=""))
  }
  
  lmcmc <- c() # find the shortest run
  for(i in 1:runs)
    lmcmc <- min(lmcmc, nrow(MCMC[[i]]))
  for(i in 1:runs)
    MCMC[[i]] <- MCMC[[i]][1:lmcmc,]

  dims <- ncol(MCMC[[1]])
  rng <- c()
  for(i in 1:runs)
    rng <- range(rng, MCMC[[i]][dims])
  layout(1)
  plot(MCMC[[1]][[dims]], type="l", bty="n", xlab="", ylab="", main="", ylim=rng)
  for(i in 2:runs)
    lines(MCMC[[i]][[dims]], col=i)

  rt <- gelman.diag(mcmc.list(lapply(MCMC, as.mcmc)), autoburnin=FALSE, transform=TRUE, confidence=0.97)
  cat("Did", runs, "Bacon runs.\n")
  cat("Gelman and Rubin Reduction Factor", rt$mpsrf, " (smaller and closer to 1 is better)\n")
  if(rt$mpsrf > 1.05)
    cat("Probably not a robust MCMC run! Too much difference between runs, above the 1.05 threshold. Increase sample size?\n", sep="") else
      cat("Robust MCMC mixing, below the 1.05 safety threshold.\n", sep="")
}



#' @name Bacon.Age.d
#' @title Output all ages for a single depth.
#' @description Output all MCMC-derived age estimates for a given depth.
#' @details Obtaining an age-depth model is often only a step towards a goal, e.g., plotting a core's 
#' fossil series ('proxies') against calendar time. Bacon.Age.d can be used to list all MCMC-derived age estimates for a given (single) depth, for example to calculate mean ages for a depth. 
#' @param d The depth of which Bacon age estimates are to be returned. Has to be a single depth.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param its The set of MCMC iterations to be used. Defaults to the entire MCMC output, \code{its=set$output}.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return Outputs all MCMC-derived ages for a given depth. 
#' @examples 
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   ages.d20 = Bacon.Age.d(20)
#'   mean(ages.d20)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#'  \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
Bacon.Age.d <- function(d, set=get('info'), its=set$output, BCAD=set$BCAD) {
  if(length(d) > 1)
	stop("Bacon.Age.d can handle one depth at a time only.\n")
  if(length(its) == 0) {
	#stop("Please run agedepth() first to re-load an existing Bacon run")
    its <- read.table(paste(set$prefix, ".out", sep="") )
  }	
  elbows <- cbind(its[,1])
  accs <- its[,2:(ncol(its)-1)]
  for(i in 2:ncol(accs))
    elbows <- cbind(elbows, elbows[,ncol(elbows)] + (set$thick * accs[,i-1]))

  if(d %in% set$d)
    ages <- elbows[, which(set$d == d)] else {
      maxd <- max(which(set$d < d))
      ages <- elbows[,maxd] + ((d-set$d[maxd]) * accs[,maxd])
    }
   if(!is.na(set$hiatus.depths)[1])
     for(hi in set$hiatus.depths) {
       below <- min(which(set$d > hi), set$K-1)+1
       above <- max(which(set$d < hi))
       if(d > set$d[above] && d < set$d[below]) {
         start.hiatus <- elbows[,below] - (its[,1+below] * (set$d[below] - hi))
         end.hiatus <- elbows[,above] + (its[,above] * (hi - set$d[above]))
         ok <- which(end.hiatus < start.hiatus)
         if(d < hi)
         ages[ok] <- elbows[ok,above] + (its[ok,above] * (d - set$d[above])) else
           ages[ok] <- elbows[ok,below] - (its[ok,1+below] * (set$d[below] - d))
       }
     }
  if(BCAD)
    ages <- 1950 - ages
  ages
}



#################### Graphs for post-run analysis ####################

#' @name Bacon.hist
#' @title Calculate age distributions of depths.
#' @description Calculate the distribution of age estimates of single or multiple depths.
#' @details Age estimates of specific depths can also be plotted.
#' @param d The depth or depths for which a histogram and age ranges should be provided. If multiple depths are given, then just the age ranges, median and means (no graphs) are provided for each depth.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param bins Bacon will try to automatically estimate reasonable bin sizes. Changing the histograms bin size could alter the shape of the histogram considerably. Values below c. 50 are generally not recommended.
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=c()}).
#' @param hist.lab The y-axis is labelled \code{ylab="Frequency"} by default.
#' @param hist.lim Limits of the y-axis.
#' @param Plot Make a plot or not. Defaults to \code{Plot=TRUE}, however no plots are made if more than one depth \code{d} is provided.
#'  If \code{Plot=FALSE}, then the age ranges, median and mean are given for each depth (as four columns).
#' @param yr.res Bacon will try to automatically estimate reasonable bin sizes. Changing the histograms bin size could alter the shape of the histogram considerably. Values below c. 50 are generally not recommended. See also \code{bins}.
#' @param prob Age ranges are given as quantiles, e.g., 2.5\% and 97.5\% for the default of 95\% confidence limits (\code{prob=0.95})). 
#' @param hist.col Colour of the histogram. Default grey, \code{hist.col=grey(0.5)}.
#' @param hist.border Colour of the histogram's outline. Default dark grey, \code{hist.border=grey(0.2)}.
#' @param range.col Colour of confidence ranges. Defaults to \code{range.col="blue"}.
#' @param med.col Colour of the median. Defaults to \code{med.col="green"}.
#' @param mean.col Colour of the mean. Defaults to \code{mn.col="red"}.
#' @param cleanup The drawing of the plot uses temporary files, which are removed afterward by default (\code{cleanup=TRUE}). 
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot with the histogram and the age ranges, median and mean, or just the age ranges, medians and means if more than one depth \code{d} is given.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   Bacon.hist(20) 
#'   Bacon.hist(20:30) 
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
Bacon.hist <- function(d, set=get('info'), bins=c(), yr.lab=c(), yr.lim=c(), hist.lab="Frequency", hist.lim=c(), Plot=TRUE, yr.res=200, prob=set$prob, hist.col=grey(0.5), hist.border=grey(.2), range.col="blue", med.col="green", mean.col="red", cleanup=TRUE) {
  tmpfile <- tempfile()
  depthfile <- tempfile()    
  write(d, depthfile, 1)
  outfile <- paste(set$prefix, ".out", sep="")
  if(!exists("set$output") || !exists("set$Tr")) {
    set <- .Bacon.AnaOut(outfile, set)		
    .assign_to_global("set", set)
  }
    
  if(length(bins) == 0)
    bins <- max(50,if(exists("set$output"))2*ceiling(sqrt(nrow(set$output))))

  hist2(outfile, set$Tr, set$d.min, set$thick, set$K, bins, tmpfile, length(d), depthfile)
  source(tmpfile)
  if(cleanup) 
    file.remove(tmpfile, depthfile)

  # here check for hiatuses and re-calculate any affected depths d. 
  if(!is.na(set$hiatus.depths)[1]) {
    d.hists <- c()
    for(j in 1:length(hists))
      d.hists[j] <- hists[[j]]$d
    for(hi in set$hiatus.depths)
      for(i in d) {
        below <- min(which(set$d > hi), set$K-1)+1
        above <- max(which(set$d < hi))
        if(i < set$d[below] && i > set$d[above]) {
          this <- c()
          if(any( round(d.hists, set$find.round) == round(i,set$find.round)) )
            this <- max(which(round(d.hists,set$find.round) == round(i,set$find.round)))[1] # ugly
          if(length(this) > 0) {
            newdist <- hist(Bacon.Age.d(i, set), bins, plot=FALSE)
            hists[[this]]$th0 <- round(min(newdist$mids))  # th0
            hists[[this]]$th1 <- round(max(newdist$mids))  # th1
            hists[[this]]$n <- length(newdist$mids)        # n
            hists[[this]]$counts <- round(newdist$counts)
            .assign_to_global("hists", hists)
          } 
        }
      }
  }

  ranges.hst <- function(hst, d) {
	mn <- weighted.mean(hst[,1], hst[,2]) # this needs to be a real mean, based on Bacon.Age.d or so
    hst <- cbind(hst[,1], cumsum(hst[,2])/sum(hst[,2]))
    qu <- approx(hst[,2], hst[,1], c((1-prob)/2, 1-((1-prob)/2), .5), rule=2)
    c(qu$y, mn)
  }

  rng <- array(NA, dim=c(length(d), 4))
  for(i in 1:length(d)) {
    yrs <- seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n)
    if(set$BCAD) yrs <- 1950 - yrs
    hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), length=yr.res))
    if(length(hst$y[hst$y>0]) < 2)
    hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), length=50)) # for very narrow ranges
    rng[i,] <- ranges.hst(cbind(hst$x, hst$y), d[i])
    if(exists("set$Tr"))
      if(exists("set$output")) # then get a better estimate of the mean (not disturbed by bin sizes)
  	    rng[i,4] <- mean(Bacon.Age.d(d[i]))
  }

  rm(tmpfile, depthfile)
  if(length(d)==1 && Plot==TRUE) {
    if(length(yr.lab) == 0)
      yr.lab <- ifelse(set$BCAD, "BC/AD", "cal yr BP")
    if(length(yr.lim) == 0)
      yr.lim <- range(hst$x)
    if(length(hist.lim) == 0)
      hist.lim <- c(0, max(hst$y))
	# mn <- mean(Bacon.Age.d(d))
    pol <- cbind(c(min(hst$x), hst$x, max(hst$x)), c(0, hst$y, 0))
    plot(0, type="n", xlim=yr.lim, ylim=hist.lim, xlab=yr.lab, ylab=hist.lab, yaxs="i")
    polygon(pol, col=hist.col, border=hist.border)
    segments(rng[1], 0, rng[2], 0, col=range.col, lwd=3)
	points(rng[3], 0, col=med.col, pch=20)
    points(rng[4], 0, col=mean.col, pch=20)

    cat("  mean (", mean.col, "): ", round(rng[4],1), " ", yr.lab,
      ", median (", med.col, "): ",  round(rng[3],1), " ", yr.lab, "\n", sep="")
    cat(100*prob, "% range (", range.col, "): ", round(rng[1],1), " to ", round(rng[2],1), yr.lab, "\n", sep="")
  } else
  return(rng)
}



#' @name proxy.ghost
#' @title Proxies analysed along the depths of a core can be plotted as 'proxy-ghost' graphs against calendar time while taking into account chronological uncertainties. Here darker grey indicates more likely calendar ages for specific proxy values.
#' @description Proxies analysed along the depths of a core can be plotted as 'proxy-ghost' graphs against calendar time while taking into account chronological uncertainties. Here darker grey indicates more likely calendar ages for specific proxy value.
#' @details Place a csv file with the values of proxies against depth within your core's folder. The values should be in columns separated by commas (default \code{sep=","}), the first column containing the depths and the first line (header) containing the proxy names. 
#' The file name should start with the core's name and end with "_proxies.csv". For an example see \code{"Bacon_coredir/MSB2K/MSB2K_proxies.csv"} or \code{"Cores/MSB2K/MSB2K_proxies.csv"}.
#' @param proxy Which proxy to use (counting from the column number in the .csv file after the depths column). 
#' @param proxy.lab Label of the proxy axis. Default names are taken from the csv file.
#' @param proxy.res Greyscale pixels are calculated for \code{proxy.res=200} proxy values by default. 
#' @param yr.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param grey.res Grey-scale resolution of the proxy graph. Default \code{grey.res=100}. See also \code{bins}.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param bins Bin sizes to calculate the grey-scales are calculated automatically. See also \code{grey.res}. 
#' @param dark By default, the darkest grey value is assigned to the most likely value within the entire core (normalised to 1; \code{dark=1}). By setting dark to, e.g., \code{dark=.8}, all values of and above 0.8 will be darkest (and values below that threshold will be lighter grey the lower their probabilities).
#' @param darkest Darkness of the most likely value. Is black by default (\code{darkest=1}); lower values will result in lighter grey.  
#' @param rotate.axes The default is to plot the calendar horizontally, however the plot can be rotated (\code{rotate.axes=TRUE}). 
#' @param proxy.rev The proxy axis can be reversed if \code{proxy.rev=TRUE}.
#' @param yr.rev The calendar axis can be reversed using \code{yr.rev=TRUE}.
#' @param plot.mean The mean ages of the proxy values can be added using \code{plot.mean=TRUE}.
#' @param mean.col Colour of the weighted mean ages of the proxy values.
#' @param yr.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=c()}).
#' @param proxy.lim Ranges of the proxy axis, calculated automatically by default (\code{proxy.lim=c()}).
#' @param sep Separator between the fields of the plain text file containing the depth and proxy data.
#' @param xaxs Extension of x-axis. By default, no white-space will be added at the axis extremes (\code{xaxs="i"}). See ?par for other options. 
#' @param yaxs Extension of y-axis. By default, no white-space will be added at the axis extremes (\code{xaxs="i"}). See ?par for other options. 
#' @param xaxt The x-axis is plotted by default, but this can be switched off using \code{xaxt="n"}.
#' @param yaxt The y-axis is plotted by default, but this can be switched off using \code{yaxt="n"}.
#' @param bty Type of box to be drawn around the plot (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param cleanup Temporary files to produce the greyscales are removed automatically after the plot has been produced (\code{cleanup=TRUE}). 
#' @author Maarten Blaauw, J. Andres Christen
#' @return A grey-scale graph of the proxy against calendar age. 
#' @examples
#' \donttest{
#'   Bacon(ask=FALSE, coredir=tempfile())
#'   layout(1)
#'   proxy.ghost()
#' }
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
proxy.ghost <- function(proxy=1, proxy.lab=c(), proxy.res=200, yr.res=200, grey.res=100, set=get('info'), bins=c(), dark=1, darkest=1, rotate.axes=FALSE, proxy.rev=FALSE, yr.rev=FALSE, plot.mean=FALSE, mean.col="red", yr.lim=c(), proxy.lim=c(), sep=",", xaxs="i", yaxs="i", xaxt="s", yaxt="s", bty="l", BCAD=set$BCAD, yr.lab=ifelse(BCAD, "BC/AD", "cal yr BP"), cleanup=TRUE) {
  if(length(set$Tr)==0)
    stop("\nPlease first run agedepth()\n\n")
  proxies <- read.csv(paste(set$coredir, set$core, "/", set$core, "_proxies.csv", sep=""), header=TRUE, sep=sep)
  if(length(proxy.lab)==0)
    proxy.lab <- names(proxies)[proxy+1]
  proxy <- cbind(as.numeric(proxies[,1]), as.numeric(proxies[,proxy+1]))
  proxy <- proxy[!is.na(proxy[,2]),]
  proxy <- proxy[which(proxy[,1] <= max(set$d)),]
  proxy <- proxy[which(proxy[,1] >= min(set$d)),]
  pr.mn.ages <- approx(set$ranges[,1], set$ranges[,5], proxy[,1], rule=1)$y
  if(length(unique(proxy[,2]))==1)
    stop("\nThis proxy's values remain constant throughout the core, and cannot be proxy-ghosted!\n\n")
  proxyseq <- seq(min(proxy[,2]), max(proxy[,2]), length=proxy.res)
  out <- list(yrseq=c(), binned=c(), maxs=c())
  ds <- c()
  d.length <- array(1, dim=c(proxy.res, 2))

  for(i in 1:proxy.res) {
    tmp  <- .DepthsOfScore(proxyseq[i], proxy)
    ds <- c(ds, tmp)
    if(i > 1)
      d.length[i,1] <- d.length[(i-1),2]+1
    d.length[i,2] <- d.length[i,1]+length(tmp)-1
    if(length(tmp) == 0)
      d.length[i,] <- d.length[i-1,]
  }
  cat("\nCalculating histograms... \n")

  Bacon.hist(ds, set, bins, Plot=FALSE, cleanup=cleanup)
  yr.min <- c()
  yr.max <- c()
  min.rng <- c(); max.rng <- c()
  hists <- get('hists') 
    
  for(i in 1:length(hists)) {
    yr.min <- min(yr.min, hists[[i]]$th0)
    yr.max <- max(yr.max, hists[[i]]$th1)
  }
  yr.seq <- seq(yr.min, yr.max, length=yr.res)

  all.counts <- array(0, dim=c(length(hists), length(yr.seq)))
  for(i in 1:length(hists))
    all.counts[i,] <- approx(seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n), hists[[i]]$counts, yr.seq)$y
  all.counts[is.na(all.counts)] <- 0
  all.counts <- all.counts/max(all.counts)
  all.counts[all.counts > dark] <- dark
  max.counts <- array(0, dim=c(proxy.res, length(yr.seq)))
  for(i in 1:proxy.res)
    for(j in 1:length(yr.seq))
      max.counts[i,j] <- max(all.counts[d.length[i,1]:d.length[i,2],j])
  if(dark>1)
    stop("Warning, dark values larger than 1 are not allowed\n") else  
      max.counts[max.counts > dark] <- dark
  if(length(yr.lim)==0)
    if(xaxs=="r")
      yr.lim <- range(pretty(c(1.04*max(yr.seq), .96*min(yr.seq)))) else
        yr.lim <- range(yr.seq)[2:1]
  if(yr.rev)
    yr.lim <- yr.lim[2:1]
  if(BCAD) {
    yr.lim <- 1950-yr.lim
    max.counts <- max.counts[,ncol(max.counts):1]
    yr.seq <- 1950-rev(yr.seq)
  }

  if(length(proxy.lim)==0)
    proxy.lim <- range(proxyseq)
  if(proxy.rev)
    proxy.lim <- proxy.lim[2:1]
  if(rotate.axes) {
    image(proxyseq, yr.seq, max.counts, xlim=proxy.lim, ylim=yr.lim, col=grey(seq(1, 1-darkest, length=grey.res)), ylab=yr.lab, xlab=proxy.lab, xaxs=xaxs, yaxs=yaxs, xaxt=xaxt, yaxt=yaxt)
    if(plot.mean) 
      lines(proxy[,2], pr.mn.ages, col=mean.col)
  } else {
    image(yr.seq, proxyseq, t(max.counts), xlim=yr.lim, ylim=proxy.lim, col=grey(seq(1, 1-darkest, length=grey.res)), xlab=yr.lab, ylab=proxy.lab, xaxs=xaxs, yaxs=yaxs, xaxt=xaxt, yaxt=yaxt)
    if(plot.mean)
      lines(pr.mn.ages, proxy[,2], col=mean.col)
}
  box(bty=bty)
}



## Accumulation rate calculations
#' should take into account hiatuses
#' @name accrate.depth
#' @title Obtain estimated accumulation rates as for any depth of a core.
#' @description Obtain accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any depth of a core.
#' @details Considering accumulation rates is crucial for age-depth modelling, and even more so if they are subsequently used for calculating proxy
#' influx values, or interpreted as proxy for environmental change such as carbon accumulation. 
#' Bacon deals explicitly with accumulation rate and its variability through defining prior distributions.
#' This function obtains accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations 
#' for any depth of a core. Deals with only 1 depth at a time. See also \code{accrate.age}.
#' @param d The depth for which accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return all MCMC estimates of accumulation rate of the chosen depth.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile()) 
#'   agedepth(yr.res=50)
#'   d20 <- accrate.depth(20)
#'   hist(d20)
#'   d20 <- accrate.depth(20, cmyr=TRUE) # to calculate accumulation rates in cm/yr
#'   mean(d20)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
accrate.depth <- function(d, set=get('info'), cmyr=FALSE) {
  if(is.na(set$hiatus.depths)) {
    if(min(set$d) <= d && max(set$d) >= d)
      accs <- set$output[,1+max(which(set$d <= d))] else accs <- NA
  } else
    for(hi in set$hiatus.depths) {
      whichbelow <- min(which(set$d > d))
      whichabove <- max(which(set$d < d))
      if(set$d[whichbelow] > hi && set$d[whichabove] < hi)
        if(d < hi)
          accs <- set$output[,1+1+whichbelow] else
            accs <- set$output[,0+whichabove]
    }
  if(cmyr) 1/accs else accs
}



# should take into account hiatuses
#' @name accrate.age
#' @title Obtain estimated accumulation rates for any age of a core.
#' @description Obtain accumulation rates (in years per cm, so actually sedimentation times) as estimated by the MCMC iterations for any age of a core.
#' @details Considering accumulation rates is crucial for age-depth modelling, and even more so if they are subsequently 
#' used for calculating proxy influx values, or interpreted as proxy for environmental change such as carbon accumulation. See also \code{accrate.age.ghost}, \code{accrate.depth} and \code{accrate.depth.ghost}.
#' Bacon deals explicitly with accumulation rate and its variability through defining prior distributions. 
#' This function obtains accumulation rates (in years per cm, so actually sedimentation times) as estimated 
#' by the MCMC iterations for any age of a core. Deals with only 1 age at a time. See also \code{accrate.depth}.
#' @param age The age for which the accumulation rates need to be returned.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return all MCMC estimates of accumulation rate of the chosen age.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   accrate.a5000 = accrate.age(5000)
#'   plot(accrate.a5000, pch='.')   
#'   hist(accrate.a5000)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#'  \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
accrate.age <- function(age, set=get('info'), cmyr=FALSE) {
 ages <- cbind(set$output[,1])
 for(i in 1:set$K)
   ages <- cbind(ages, ages[,i] + (set$thick * (set$output[,i+1])))

 if(age < min(ages) || age > max(ages))
   cat(" Warning, age outside the core's age range!\n")
 accs <- c()
 for(i in 2:ncol(ages)) {
   these <- (ages[,i-1] < age) * (ages[,i] > age)
   if(sum(these) > 0) # age lies within these age-model iterations
     accs <- c(accs, set$output[which(these>0),i+1])
  }
  if(cmyr) 
    accs <- 1/accs 
  return(accs)
}



#' @name accrate.depth.ghost
#' @title Plot modelled accumulation rates against the depths of a core.
#' @description Plot grey-scale representation of modelled accumulation rates over a core's depth. Each section of the core (see Bacon's option \code{"thick"}) will have modelled accumulation rates. 
#' @details This plot shows the modelled accumulation rates in grey-scales, where darker grey indicates more likely accumulation rates.
#' Axis limits for accumulation rates are estimated automatically, however upper limits can be very variable (and thus hard to predict) 
#' if calculated in cm/yr; therefore you might want to manually adapt the axis limits after plotting with default settings (e.g., \code{acc.lim=c(0,1)}). See also \code{accrate.age.ghost}, \code{accrate.depth} and \code{accrate.age}.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param d The depths for which the accumulation rates are to be calculated. Default to the entire core. 
#' @param d.lim Axis limits for the depths.
#' @param acc.lim Axis limits for the accumulation rates.
#' @param d.lab Label for the depth axis. 
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm. Axis limits are difficult to calculate when \code{cmyr=TRUE}, so a manual adaptation of \code{acc.lim} might be a good idea.
#' @param acc.lab Axis label for the accumulation rate. 
#' @param dark The darkest grey value is dark=1 by default; lower values will result in lighter grey but values >1 are not advised. 
#' @param grey.res Grey-scale resolution. Default \code{grey.res=100}. See also \code{bins}.
#' @param prob Probability ranges. Defaults to \code{prob=0.95}.
#' @param plot.range If \code{plot.range=TRUE}, the confidence ranges (two-tailed; half of the probability at each side) are plotted. 
#' @param range.col Colour of the confidence ranges.
#' @param range.lty Line type of the confidence ranges.
#' @param plot.mean If \code{plot.mean=TRUE}, the means are plotted. 
#' @param mean.col Colour of the mean accumulation rates.
#' @param mean.lty Type of the mean lines.
#' @param rotate.axes The default is to plot the accumulation rates horizontally and the depth vertically (\code{rotate.axes=FALSE}). Change rotate.axes value to rotate axes.
#' @param rev.d The direction of the depth axis can be reversed from the default (\code{rev.d=TRUE}.
#' @param rev.acc The direction of the accumulation rate axis can be reversed from the default (\code{rev.acc=TRUE}).
#' @author Maarten Blaauw, J. Andres Christen
#' @return A grey-scale plot of accumulation rate against core depth.
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   layout(1)
#'   accrate.depth.ghost()
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
accrate.depth.ghost <- function(set=get('info'), d=set$d, d.lim=c(), acc.lim=c(), d.lab=c(), cmyr=FALSE, acc.lab=c(), dark=1, grey.res=100, prob=0.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, rotate.axes=FALSE, rev.d=FALSE, rev.acc=FALSE) {
  max.acc <- 0; max.dens <- 0
  acc <- list(); min.rng <- c(); max.rng <- c(); mn.rng <- c()
  for(i in 1:length(d))
    if(length(acc.lim) == 0)
      acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr), from=0) else
        acc[[i]] <- density(accrate.depth(d[i], set, cmyr=cmyr), from=0, to=max(acc.lim))
  for(i in 1:length(d)) {
    max.acc <- max(max.acc, acc[[i]]$x)
    max.dens <- max(max.dens, acc[[i]]$y)
	quants <- quantile(accrate.depth(d[i], set, cmyr=cmyr), c((1-prob)/2, 1-((1-prob)/2)))
	min.rng[i] <- quants[1]
	max.rng[i] <- quants[2]
	mn.rng[i] <- mean(accrate.depth(d[i], set, cmyr=cmyr))
   }
    
  for(i in 1:length(d)) {
    acc[[i]]$y <- acc[[i]]$y/max.dens
    acc[[i]]$y[acc[[i]]$y > dark] <- dark
  }

  if(length(d.lim) == 0)
    d.lim <- range(set$d)  
  if(length(d.lab) == 0)
    d.lab <- paste("depth (", set$unit, ")", sep="")
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste("accumulation rate (", set$unit, "/yr)", sep="") else
        acc.lab <- paste("accumulation rate (yr/", set$unit, ")", sep="") 
  
  if(rev.d) 
    d.lim <- rev(d.lim)
  if(length(acc.lim) == 0)
    acc.lim <- c(0, max.acc)
  if(rev.acc)
    acc.lim <- rev(acc.lim)

  if(rotate.axes) {
    plot(0, type="n", xlab=acc.lab, ylab=d.lab, ylim=d.lim, xlim=acc.lim)
    for(i in 2:length(d))
      image(acc[[i]]$x, d[c(i-1, i)], t(1-t(acc[[i]]$y)), add=TRUE, col=grey(seq(1-max(acc[[i]]$y), 1, length=grey.res)))
    if(plot.range) {
      lines(min.rng, d-(set$thick/2), col=range.col, lty=range.lty)  	
      lines(max.rng, d-(set$thick/2), col=range.col, lty=range.lty)
    }
	if(plot.mean)    	
      lines(mn.rng, d-(set$thick/2), col=mean.col, lty=mean.lty)  		 
  } else {
      plot(0, type="n", xlab=d.lab, ylab=acc.lab, xlim=d.lim, ylim=acc.lim)
      for(i in 2:length(d))
        image(d[c(i-1, i)], acc[[i]]$x, 1-t(acc[[i]]$y), add=TRUE, col=grey(seq(1-max(acc[[i]]$y), 1, length=grey.res)))
      if(plot.range) {
        lines(d-(set$thick/2), min.rng, col=range.col, lty=range.lty)  		 
        lines(d-(set$thick/2), max.rng, col=range.col, lty=range.lty)
        }
    if(plot.mean)
      lines(d-(set$thick/2), mn.rng, col=mean.col, lty=mean.lty)  		 
    }
}



#' @name accrate.age.ghost
#' @title Plot a core's accumulation rates against calendar time. 
#' @description Plot a grey-scale representation of a core's estimated accumulation rates against time. 
#' @details Calculating accumulation rates against calendar age will take some time to calculate, and might show unexpected 
#' rates around the core's maximum ages (only a few of all age-model iterations will reach such ages and they will tend to have
#'  modelled accumulation rates for the lower depths much lower than the other iterations). Axis limits for accumulation rates
#'   are estimated automatically, however upper limits can be very variable (and thus hard to predict) if calculated in \code{cm/yr}.
#'  Therefore you might want to manually adapt the axis limits after plotting with default settings (e.g., \code{acc.lim=c(0,1)}). See also \code{accrate.depth.ghost}, \code{accrate.depth} and \code{accrate.age}.
#' The grey-scale reconstruction around the oldest ages of any reconstruction often indicates very low accumulation rates. 
#' This is due to only some MCMC iterations reaching those old ages, and these iterations will have modelled very slow accumulation rates. 
#' Currently does not deal well with hiatuses, so do not interpret accumulation rates close to depths with inferred hiatuses.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param yr.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=c()}).
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param grey.res Resolution of greyscales. Default \code{grey.res=50}, which does not aim to poke fun at a famous novel.
#' @param prob Probability ranges. Defaults to \code{prob=0.95}.
#' @param plot.range If \code{plot.range=TRUE}, the confidence ranges (two-tailed; half of the probability at each side) are plotted. 
#' @param range.col Colour of the confidence ranges.
#' @param range.lty Line type of the confidence ranges.
#' @param plot.mean If \code{plot.mean=TRUE}, the means are plotted. 
#' @param mean.col Colour of the mean accumulation rates.
#' @param mean.lty Type of the mean lines.
#' @param acc.lim Axis limits for the accumulation rates.
#' @param acc.lab Axis label for the accumulation rate. 
#' @param upper Maximum accumulation rates to plot. Defaults to the upper 99\%; \code{upper=0.99}.
#' @param dark The darkest grey value is dark=1 by default; lower values will result in lighter grey but values >1 are not advised. 
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param cmyr Accumulation rates can be calculated in cm/year or year/cm. By default \code{cmyr=FALSE} and accumulation rates are calculated in year per cm. Axis limits are difficult to calculate when \code{cmyr=TRUE}, so a manual adaptation of \code{acc.lim} might be a good idea.
#' @param rotate.axes The default is to plot the calendar age horizontally and accumulation rates vertically. Change to \code{rotate.axes=TRUE} value to rotate axes.
#' @param rev.yr The direction of the age axis, which can be reversed using \code{rev.yr=TRUE}.
#' @param rev.acc The direction of the accumulation rate axis, which can be reversed (\code{rev.acc=TRUE}.
#' @param xaxs Extension of the x-axis. White space can be added to the vertical axis using \code{xaxs="r"}.	
#' @param yaxs Extension of the y-axis. White space can be added to the vertical axis using \code{yaxs="r"}.
#' @param bty Type of box to be drawn around the plot (\code{"n"} for none, and \code{"l"} (default), \code{"7"}, \code{"c"}, \code{"u"}, or \code{"o"} for correspondingly shaped boxes).
#' @author Maarten Blaauw, J. Andres Christen
#' @return A greyscale plot of accumulation rate against calendar age. 
#' @examples
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   layout(1)
#'   accrate.age.ghost()
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
accrate.age.ghost <- function(set=get('info'), yr.lim=c(), yr.lab=c(), yr.res=200, grey.res=50, prob=.95, plot.range=TRUE, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, acc.lim=c(), acc.lab=c(), upper=0.99, dark=50, BCAD=set$BCAD, cmyr=FALSE, rotate.axes=FALSE, rev.yr=FALSE, rev.acc=FALSE, xaxs="i", yaxs="i", bty="o") {
  if(length(yr.lim) == 0) {
    min.age <- min(set$ranges[,2])
    max.age <- max(set$ranges[,3])
  } else {
      min.age <- min(yr.lim)
      max.age <- max(yr.lim)
    }
   
  yr.seq <- seq(min.age, max.age, length=yr.res)
  max.y <- 0; all.x <- c()
  hist.list <- list(x=c(), y=c(), min.rng=c(), max.rng=c(), mn.rng=c())
  for(i in 1:length(yr.seq)) {
    if(!(i %% 50))
      cat(".")
      acc <- accrate.age(yr.seq[i], set, cmyr=cmyr)
    if(cmyr) acc <- rev(acc)
	accs <- acc	
    if(length(acc) > 2)
    if(length(acc.lim) == 0)
      acc <- density(acc, from=0) else
        acc <- density(acc, from=0, to=max(acc.lim))
	hist.list$yr[[i]] <- yr.seq[i]		
    hist.list$x[[i]] <- acc$x
    hist.list$y[[i]] <- acc$y/sum(acc$y)
    rng <- quantile(accs, c((1-prob)/2, 1-((1-prob)/2)))
	hist.list$mn.rng[[i]] <- mean(accs)	
	hist.list$min.rng[[i]] <- rng[1]
	hist.list$max.rng[[i]] <- rng[2]	
    max.y <- max(max.y, hist.list$y[[i]])
    all.x <- c(all.x, acc$x)
  }
  cat("\n")
  
  if(BCAD) 
    yr.seq <- 1950 - yr.seq
  if(length(yr.lim) == 0)
    yr.lim <- range(yr.seq)
  if(length(acc.lim) == 0)
    acc.lim <- c(0, 1.1*quantile(all.x, upper))
  if(rev.yr)
    yr.lim <- rev(yr.lim)
  if(rev.acc)
    acc.lim <- rev(acc.lim)
  if(length(yr.lab) == 0)
	  if(BCAD)
        yr.lab <- "BC/AD" else
          yr.lab <- "cal BP"    
  if(length(acc.lab) == 0)
    if(cmyr)
      acc.lab <- paste("accumulation rate (", set$unit, "/yr)", sep="") else
        acc.lab <- paste("accumulation rate (yr/", set$unit, ")", sep="") 
  if(rotate.axes) {
    plot(0, type="n", xlim=acc.lim, xlab=acc.lab, ylim=yr.lim, ylab=yr.lab, xaxs=xaxs, yaxs=yaxs)
    for(i in 2:length(yr.seq))
      image(sort(hist.list$x[[i]]), hist.list$yr[c(i-1,i)], t(t(matrix(hist.list$y[[i]]))),
        col=grey(seq(1, 1-min(1,max(hist.list$y[[i]])*dark/max.y), length=grey.res)), add=TRUE)
    } else {
      plot(0, type="n", xlim=yr.lim, xlab=yr.lab, ylim=acc.lim, ylab=acc.lab, xaxs=xaxs, yaxs=yaxs)
      for(i in 2:length(yr.seq))
        image(sort(hist.list$yr[c(i-1,i)]), hist.list$x[[i]], t(matrix(hist.list$y[[i]])),
          col=grey(seq(1, 1-min(1,max(hist.list$y[[i]])*dark/max.y), length=grey.res)), add=TRUE)
      }
   if(plot.range)
	 if(rotate.axes) {
       lines(hist.list$min.rng, yr.seq, pch=".", col=range.col, lty=range.lty)
       lines(hist.list$max.rng, yr.seq, pch=".", col=range.col, lty=range.lty)
	 } else {
		 lines(yr.seq, hist.list$min.rng, pch=".", col=range.col, lty=range.lty)
		 lines(yr.seq, hist.list$max.rng, pch=".", col=range.col, lty=range.lty)
	   } 	  
	 if(plot.mean)
       if(rotate.axes)
         lines(hist.list$mn.rng, yr.seq, pch=".", col=mean.col, lty=mean.lty) else
           lines(yr.seq, hist.list$mn.rng, pch=".", col=mean.col, lty=mean.lty)
      
  if(length(bty) > 0)
	box(bty=bty)  	  
}



#' @name flux.age.ghost
#' @title Plot flux rates for proxies.
#' @description Plot grey-scale representation of estimated flux rates for proxies against calendar age.
#' @details To plot flux rates (e.g. pollen grains/cm2/yr) as greyscales, 
#' provide a plain text file with headers and the data in columns separated by commas, ending in '_flux.csv'
#' and saved in your core's folder. The first column should contain the depths, and the next columns should contain 
#' the proxy concentration values (leaving missing values empty). Then type for example \code{flux.age.ghost(1)} to plot the 
#' flux values for the first proxy in the .csv file. Instead of using a _flux.csv file, a flux variable can also be defined
#'  within the R session (consisting of depths and their proxy concentrations in two columns). Then provide the name of this variable, e.g.: \code{flux.age.ghost(flux=flux1)}. 
#' See Bacon_runs/MSB2K/MSB2K_flux.csv for an example.
#' @param proxy Which proxy to use (counting from the column number in the .csv file after the depths column). 
#' @param yr.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=c()}).
#' @param yr.res Resolution or amount of greyscale pixels to cover the age scale of the plot. Default \code{yr.res=200}.
#' @param set Detailed information of the current run, stored within this session's memory as variable info.
#' @param flux Define a flux variable within the R session (consisting of depths and their proxy concentrations in two columns) and provide the name of this variable, e.g.: 
#' \code{flux.age.ghost(flux=flux1)}. If left empty (\code{flux=c()}), a flux file is expected (see \code{proxy}).
#' @param plot.range Plot curves that indicate a probability range, at resolution of yr.res.
#' @param prob	Probability range, defaults to \code{prob=0.8} (10 \% at each side). 
#' @param range.col Red seems nice. 
#' @param range.lty Line type of the confidence ranges.
#' @param flux.lim Limits of the flux axes.
#' @param flux.lab Axis labels. Defaults to \code{flux.lab="flux"}.
#' @param plot.mean Plot the mean fluxes.
#' @param mean.col Red seems nice.
#' @param mean.lty Line type of the means.
#' @param upper Maximum flux rates to plot. Defaults to the upper 99\%; \code{upper=0.99}.
#' @param dark The darkest grey value is \code{dark=1} by default; lower values will result in lighter grey but \code{values >1} are not allowed.
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param rotate.axes The default of plotting calendar year on the horizontal axis and fluxes on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.flux The flux axis can be reversed with \code{rev.flux=TRUE}.
#' @param rev.yr The direction of the age axis can be reversed using \code{rev.yr=TRUE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A plot of flux rates.
#' @examples
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   flux.age.ghost(1)
#' }
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
flux.age.ghost <- function(proxy=1, yr.lim=c(), yr.res=200, set=get('info'), flux=c(), plot.range=TRUE, prob=.8, range.col=grey(0.5), range.lty=2, plot.mean=TRUE, mean.col="red", mean.lty=2, flux.lim=c(), flux.lab="flux", upper=.95, dark=set$dark, BCAD=set$BCAD, yr.lab=c(), rotate.axes=FALSE, rev.flux=FALSE, rev.yr=FALSE) {
  if(length(flux) == 0) { # then read a .csv file, expecting data in columns with headers
    flux <- read.csv(paste(set$coredir, set$core, "/", set$core, "_flux.csv", sep=""))
    flux <- cbind(flux[,1], flux[,1+proxy])
      isNA <- is.na(flux[,2])
      flux <- flux[!isNA,]
  }
  if(length(yr.lim) == 0) {
    min.age <- min(set$ranges[,2])
    max.age <- max(set$ranges[,3])
    yr.lim <- c(min.age, max.age)
  } else {
      min.age <- min(yr.lim)
      max.age <- max(yr.lim)
    }
  
  age.seq <- seq(min(min.age, max.age), max(min.age, max.age), length=yr.res)
  fluxes <- array(NA, dim=c(nrow(set$output), length(age.seq)))
  for(i in 1:nrow(set$output)) {
    if(!(i %% 50)) cat(".")
    ages <- as.numeric(set$output[i,1:(ncol(set$output)-1)]) # 1st step to calculate ages for each set$d
    ages <- c(ages[1], ages[1]+set$thick * cumsum(ages[2:length(ages)])) # now calculate the ages for each set$d
    ages.d <- approx(ages, c(set$d, max(set$d)+set$thick), age.seq, rule=1)$y # find the depth belonging to each age.seq, NA if none
    ages.i <- floor(approx(ages, (length(set$d):0)+1, age.seq, rule=2)$y) # find the column belonging to each age.seq
    flux.d <- approx(flux[,1], flux[,2], ages.d)$y # interpolate flux (in depth) to depths belonging to each age.seq
    fluxes[i,] <- flux.d / as.numeric(set$output[i,(1+ages.i)]) # (amount / cm^3) / (yr/cm) = amount * cm-2 * yr-1
  }
  cat("\n")
  if(length(flux.lim) == 0)
    flux.lim <- c(0, quantile(fluxes[!is.na(fluxes)], upper))
  max.dens <- 0
  for(i in 1:length(age.seq)) {
    tmp <- fluxes[!is.na(fluxes[,i]),i] # all fluxes that fall at the required age.seq age
    if(length(tmp) > 0)
      max.dens <- max(max.dens, density(tmp, from=0, to=max(flux.lim))$y)
  }

  if(length(yr.lim) == 0)
    yr.lim <- range(age.seq)
  if(length(yr.lab) == 0)
    yr.lab <- ifelse(BCAD, "BC/AD", "cal BP")
  if(rotate.axes)
    plot(0, type="n", ylim=yr.lim, ylab=yr.lab, xlim=flux.lim, xlab=flux.lab) else
      plot(0, type="n", xlim=yr.lim, xlab=yr.lab, ylim=flux.lim, ylab=flux.lab)
  min.rng <- c(); max.rng <- c(); mn.rng <- c() 	  	  
  for(i in 2:length(age.seq)) {
    tmp <- fluxes[!is.na(fluxes[,i]),i] # all fluxes that fall at the required age.seq age
	rng <- quantile(tmp, c((1-prob)/2, 1-((1-prob)/2)))
	min.rng[i] <- rng[1]
	max.rng[i] <- rng[2]
	mn.rng[i] <- mean(tmp)
    if(length(tmp[tmp>=0]) > 2) {
      flux.hist <- density(tmp, from=0, to=max(flux.lim))
      flux.hist$y <- flux.hist$y - min(flux.hist$y) # no negative fluxes
      if(rotate.axes) 
        image(flux.hist$x, age.seq[c(i-1,i)], matrix(flux.hist$y/max.dens), add=TRUE,
          col=grey(seq(1, max(0, 1-dark*(max(flux.hist$y)/max.dens)), length=100))) else 
          image(age.seq[c(i-1,i)], flux.hist$x, t(matrix(flux.hist$y/max.dens)), add=TRUE,
            col=grey(seq(1, max(0, 1-dark*(max(flux.hist$y)/max.dens)), length=100)))	
    }
  } 
 if(plot.range)
   if(rotate.axes) {
     lines(min.rng, age.seq, col=range.col, lty=range.lty)
     lines(max.rng, age.seq, col=range.col, lty=range.lty)
 } else {
     lines(age.seq, min.rng, col=range.col, lty=range.lty)
     lines(age.seq, max.rng, col=range.col, lty=range.lty)	
   }
  if(plot.mean)
    if(rotate.axes)
   	 lines(mn.rng, age.seq, col=mean.col, lty=mean.lty) else
       lines(age.seq, mn.rng, col=mean.col, lty=mean.lty)  
}



#' @name AgesOfEvents
#' @title Event probabilities against calendar age
#' @description Plot probability curves for events in the core, expressed against calendar age.
#' @details Probabilities of depths with 'events' in an age-modelled core can be plotted against time, taking into account 
#' chronological uncertainties (Blaauw et al. 2007). Such events could be for example core depths at which proxies 
#' indicate changes toward wetter local conditions. This can be expressed as values between 0 (no event) and 1 (event at 100\% probability) 
#' for each depth. 
#' 
#' Blaauw et al. 2010 propose to estimate probabilities of events by finding specific proxy features such as increasing curves. 
#' Probabilities are then estimated through resampling from the proxy values, where low to modest rises of proxy curves result 
#' in low event probabilities, and clear proxy rises in high probabilities. A smooth spline can be applied to adapt the balance of short-term vs long-term events. To calculate the event probabilities, 
#' produce a file with two columns (depth and corresponding proxy-derived probabilities, separated by white spaces).
#'  Do not provide headers at the file's first line, and save the file with extension "_events.txt" within the core's Bacon folder. See Cores/MSB2K/MSB2K_events.txt (or Bacon_runs/MSB2K/MSB2K_events.txt) for an example. Events are calculated as the probability that an event took place within specific time windows - or more specifically, that the Bacon age-depth model puts depths with assigned event probabilities in that time window.  
#'
#' does not yet deal correctly with hiatuses.
#' @param window Width of the window.
#' @param move Step size with which the window moves.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @param plot.steps Plot probability values step-wise (defaults to \code{plot.steps=FALSE}, which plots smooth curves instead).
#' @param BCAD The calendar scale of graphs and age output-files is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @param yr.lab The labels for the calendar axis (default \code{yr.lab="cal BP"} or \code{"BC/AD"} if \code{BCAD=TRUE}).
#' @param yr.lim Minimum and maximum calendar age ranges, calculated automatically by default (\code{yr.lim=c()}).
#' @param prob.lab Label of the probability axis (default \code{prob.lab="probability"}).
#' @param prob.lim Limits of the probability axis (calculated automatically by default).
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param rev.yr The direction of the age axis, which can be reversed using \code{rev.yr=TRUE}.
#' @param yaxs Extension of the y-axis. Defaults to the exact ranges of the probability values. White space can be added to the vertical axis using \code{yaxs="r"}.
#' @param bty Type of box to be drawn around plots. Draw a box around the graph (\code{"n"} for none, and \code{"l"}, \code{"7"},
#'  \code{"c"}, \code{"u"}, "]" or \code{"o"} for correspondingly shaped boxes).
#' @author Maarten Blaauw, J. Andres Christen
#' @return The resulting probabilities are plotted and saved within the core's folder (file names ending with the window width and "_probs.txt").
#' @examples
#' \donttest{
#'   Bacon(run=FALSE, coredir=tempfile())
#'   agedepth(yr.res=50)
#'   AgesOfEvents(100, 10)
#' }
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M., Christen, J.A., Mauquoy, D., van der Plicht, J., Bennett, K.D. (2007) Testing the timing of radiocarbon-dated events between proxy archives. _The Holocene_, *17*, 283-288.
#' Blaauw, M., Wohlfarth, B., Christen, J.A., Ampel, L., Veres, D., Hughen, K.A., Preusser, F., Svensson, A. (2010) Were last glacial climate events simultaneous between Greenland and France? A quantitative comparison using non-tuned chronologies. _Journal of Quaternary Science_ *25*, 387-394.
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
AgesOfEvents <- function(window, move, set=get('info'), plot.steps=FALSE, BCAD=set$BCAD, yr.lab=c(), yr.lim=c(), prob.lab="probability", prob.lim=c(), rotate.axes=FALSE, rev.yr=TRUE, yaxs="i", bty="l") {
  if(move == 0)
    stop("I cannot move anywhere if move = 0\n")
  outfile <- paste(set$prefix, "_", window, "_probs.txt", sep="")
  file.create(outfile)
  MCMCname <- paste(set$prefix, ".out", sep="")
  probfile <- paste(set$coredir, set$core, "/", set$core, "_events.txt", sep="")
  if(!file.exists(probfile))
    stop("\nFile with probabilities for events per depth not found! Check the manual\n\n")
  probs <- read.table(probfile)
  if(!is.numeric(probs[1,1]))
    stop("\nFirst line of the _events.txt file should NOT contain titles; please remove them\n\n")
  if(min(probs[,1]) < min(set$d) || max(probs[,1]) > max(set$d)) {
    cat("\nSome depths in the _events.txt file go beyond the age-model; I will remove them\n\n")
    file.rename(probfile, paste(probfile, "_backup", sep=""))
    probs <- probs[which(probs[,1] >= min(set$d)),]
    probs <- probs[which(probs[,1] <= max(set$d)),]
    write.table(probs, probfile, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  
  if(length(yr.lim) == 0) {
    min.age=min(set$ranges[,2])
    max.age=max(set$ranges[,3])
    yr.lim=c(min.age, max.age)
  } else {
      min.age <- min(yr.lim)
      max.age <- max(yr.lim)
    }
  
  events(min.age, max.age, move, window, outfile, MCMCname, nrow(set$output), set$K, set$d[1], set$thick, probfile, nrow(probs))
  probs <- read.table(outfile)
  if(BCAD) {
    probs[,1] <- 1950 - probs[,1]
    o <- order(probs[,1])
    probs <- probs[o,]
  }

  if(plot.steps) {
    d.sort <- sort(rep(1:nrow(probs),2))
    d.sort <- cbind(d.sort[-length(d.sort)], d.sort[-1])
    probs <- cbind(c(min(probs[,1]), probs[d.sort[,1],1], max(probs[,1])), c(0,probs[d.sort[,2],2],0))
  } else
      probs <- cbind(c(min(probs[,1]), probs[,1], max(probs[,1])), c(0,probs[,2],0))
  par(yaxs=yaxs, bty=bty)

  if(rev.yr)
    yr.lim <- rev(yr.lim)
  if(length(yr.lab) == 0)
    yr.lab <- ifelse(BCAD, "BC/AD", "cal BP")
  if(length(prob.lim) == 0)
    prob.lim <- c(0, 1.1*max(probs[,2]))

  if(rotate.axes) {
    plot(probs[,2], probs[,1], type="n", ylab=yr.lab, xlab=prob.lab, xlim=prob.lim, ylim=yr.lim)
    polygon(probs[,2:1], col="grey")
  } else {
      plot(probs, type="n", xlab=yr.lab, ylab=prob.lab, ylim=prob.lim, xlim=yr.lim)
      polygon(probs, col="grey")
    }
  if(move > window) cat("\nAre you sure you want the window widths to be smaller than the moves?\n")
  invisible(probs)
}



#################### user-invisible plot functions ####################

# to plot greyscale/ghost graphs of age-depth model
.depth.ghost <- function(set=get('info'), bins=c(), d.min=set$d.min, d.max=set$d.max, rotate.axes=FALSE, d.res=200, yr.res=200, grey.res=100, dark=set$dark, greyscale=function(x) grey(1-x), cleanup=TRUE, xaxt="s", yaxt="s") {
  d <- seq(d.min, d.max, length=d.res) 
  Bacon.hist(d, set, bins, yr.res=yr.res, cleanup=cleanup)
  d.jumps <- diff(d)[1]
  hists <- get('hists') 
    
  peak <- 0
  for(i in 1:(length(d))) {
    yrs <- seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n)
    hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), length=yr.res))
    if(length(hst$y>0) < 5) # for narrow distributions, go narrow
      hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), length=20))
    peak <- max(peak, max(hst$y))
  }

  for(i in 1:length(d)) {
    yrs <- seq(hists[[i]]$th0, hists[[i]]$th1, length=hists[[i]]$n)
    hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), length=yr.res))
    if(length(hst$y>0) < 5) # for narrow distributions, go narrow
      hst <- approx(yrs, hists[[i]]$counts, seq(min(yrs), max(yrs), length=20))
    if(set$BCAD) {
      hst$x <- rev(1950-hst$x)
      hst$y <- rev(hst$y)
    }
    hst$y <- hst$y/peak
    hst$y[hst$y>dark] <- max(hst$y)
    if(dark>1) stop(" Warning, the dark value should be <1")
    gr <- seq(0, dark, length=grey.res)
    if(rotate.axes)
      image(hst$x, c(d[i]-(d.jumps/2), d[i]+(d.jumps/2)), t(t(hst$y)), add=TRUE, col=greyscale(gr)) else
        image(c(d[i]-(d.jumps/2), d[i]+(d.jumps/2)), hst$x, t(hst$y), add=TRUE, col=greyscale(gr))
  }
}


# Time series of the log of the posterior
.PlotLogPost <- function(set, from=0, to=set$Tr)
  plot(from:(to-1), -set$Us[(from+1):to], type="l",
    ylab="Log of Objective", xlab="Iteration", main="")



# calibrate C14 dates and calculate distributions for any calendar dates
.bacon.calib <- function(dat, set=get('info'), date.res=100, normal=set$normal, t.a=set$t.a, t.b=set$t.b, deltaR=set$deltaR, deltaSTD=set$deltaSTD, ccdir="") {
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
      deltaR <- 0; deltaSTD <- 0 # only C14 dates should need correcting for age offsets
      x <- seq(min(dat[,2])-max(100,4*max(dat[,3])), max(dat[,2])+max(100,4*max(dat[,3])), length=date.res)
      ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
    } else
        if(set$cc==1) ccurve <- cc1 else
          if(set$cc==2) ccurve <- cc2 else
            if(set$cc==3) ccurve <- cc3 else
              ccurve <- cc4
    for(i in 1:nrow(dat))
      calib$probs[[i]] <- d.cal(ccurve, dat[i,2]-deltaR, dat[i,3]^2+deltaSTD^2, set$t.a, set$t.b)
  } else
      for(i in 1:nrow(dat)) {
        det <- as.numeric(dat[i,])
        if(det[5]==0) {
          x <- seq(det[2]-max(100,4*det[3]), det[2]+max(100,4*det[3]), length=date.res)
          ccurve <- cbind(x, x, rep(0,length(x))) # dummy 1:1 curve
        } else
            if(det[5]==1) ccurve <- cc1 else if(det[5]==2) ccurve <- cc2 else
              if(det[5]==3) ccurve <- cc3 else ccurve <- cc4

        deltaR <- set$deltaR; deltaSTD <- set$deltaSTD; t.a <- set$t.a; t.b <- set$t.b
        if(length(det) >= 7 && det[5] > 0) { # the user provided age offsets; only for C14 dates
          deltaR <- det[6]
          deltaSTD <- det[7]
        }

        if(length(det) >= 9) { # the user provided t.a and t.b values for each date
          t.a <- det[8]
          t.b <- det[9]
          if(round(t.b-t.a) !=1) 
            stop("\n Warning! t.b - t.a should always be 1, check the manual")
        }
        calib$probs[[i]] <- d.cal(ccurve, det[2]-deltaR, det[3]^2+deltaSTD^2, t.a, t.b)
      }
  calib
}



# plot the prior for the accumulation rate
.PlotAccPrior <- function(s, mn, set=get('info'), main="", xlim=c(0, 3*max(mn)), xlab=paste("Acc. rate (yr/", noquote(set$unit), ")", sep=""), ylab="Density", add=FALSE, legend=TRUE, cex=.9) {
  o <- order(s, decreasing=TRUE)
  priors <- unique(cbind(s[o],mn[o])[,1:2])
  x <- 0
  if(length(priors) == 2) {
    curve(dgamma(x, s, s/mn), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("acc.shape: ", priors[1], "\nacc.mean: ", priors[2])
  } else {
      curve(dgamma(x, priors[1,1], priors[1,1]/priors[1,2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
      for(i in 2:nrow(priors))
        curve(dgamma(x, priors[i,1], priors[i,1]/priors[i,2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=if(i==1) add else TRUE)
      txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
    }
  if(legend)
    legend("topright", txt, bty="n", cex=cex, text.col=2, adj=c(0,0))
}



# plot the prior for the memory (= accumulation rate varibility between neighbouring depths)
.PlotMemPrior <- function(s, mn, thick, ds=1, set=get('info'), xlab="Memory (ratio)", ylab="Density", main="", add=FALSE, legend=TRUE, cex=.9) {
  o <- order(s, decreasing=TRUE)
  priors <- unique(cbind(s[o],mn[o])[,1:2])
  x <- 0

  if(length(priors)==2) {
    curve(dbeta(x, s*mn, s*(1-mn)), 0, 1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("mem.strength: ", s, "\nmem.mean: ", mn, "\n", set$K, " ", thick, noquote(set$unit), " sections", sep="")
  } else {
      curve(dbeta(x, priors[1,1]*priors[1,2], priors[1,1]*(1-priors[1,2])), 0, 1, col=3, lwd=2, xlab=xlab, ylab=ylab, add=add)
      for(i in 2:nrow(priors))
        curve(dbeta(x, priors[i,1]*priors[i,2], priors[i,1]*(1-priors[i,2])), 0, 1, col=3, lwd=2, xlab="", ylab="", add=TRUE)
      txt <- paste("acc.shape: ", toString(priors[,1]), "\nacc.mean: ", toString(priors[,2]))
    }
  if(legend)
    legend("topright", txt, bty="n", cex=cex, text.col=2, adj=c(0,0))
  warn <- FALSE
  for(i in s)
    for(j in mn)
      if(i*(1-j) <= 1) warn <- 1
  if(warn) cat("\nWarning! Chosen memory prior might cause problems.\nmem.strength * (1 - mem.mean) should be smaller than 1\n ")
}



# plot the prior for the hiatus length
.PlotHiatusPrior <- function(s=set$hiatus.shape, mn=set$hiatus.mean, hiatus=set$hiatus.depths, set=get('info'), xlab="Hiatus size (yr)", ylab="Density", main="", xlim=c(0, 3*max(mn)), add=FALSE, legend=TRUE) {
  o <- order(s, decreasing=TRUE)
  priors <- unique(cbind(s[o],mn[o])[,1:2])
  x <- 0

  if(length(priors) == 2) {
    curve(dgamma(x, priors[1], priors[1]/priors[2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=add)
    txt <- paste("h.shape: ", toString(priors[1]), "\nh.mean: ", toString(priors[2]), "\nd: ", toString(hiatus))
  } else
    for(i in 2:nrow(priors)) {
      curve(dgamma(x, priors[i,1], priors[i,1]/priors[i,2]), col=3, lwd=2, xlim=xlim, xlab=xlab, ylab=ylab, add=if(i==1) add else TRUE)
      txt <- paste("h.shape: ", toString(priors[,1]), "\nh.mean: ", toString(priors[,2]), "\nd: ", toString(hiatus))
    }
  if(legend)
    legend("topright", txt, bty="n", cex=.7, text.col=2, adj=c(0,0))
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
  post <- density(post)
  maxprior <- dgamma((s-1)/(s/mn), s, s/mn)
  if(is.infinite(max(maxprior))) 
    max.y <- max(post$y) else
      max.y <- max(maxprior, post$y)
  lim.x <- range(0, post$x, 2*mn)
  plot(0, type="n", xlim=lim.x, xlab=xlab, ylim=c(0, 1.1*max.y), ylab="")
  polygon(post, col=grey(.8), border=grey(.4))
  .PlotAccPrior(s, mn, add=TRUE, xlim=range(post$x), xlab="", ylab=ylab, main=main)
}



# plot the posterior (and prior) of the memory
.PlotMemPost <- function(set=get('info'), corenam, K, main="", s=set$mem.strength, mn=set$mem.mean, xlab=paste("Memory"), ylab="Density", ds=1, thick) {
  post <- density(set$output[,set$n]^(1/set$thick), from=0, to=1)
  post <- cbind(c(min(post$x), post$x, max(post$x)), c(0, post$y, 0))
  maxprior <- max(dbeta((0:100)/100, s*mn, s*(1-mn)))
  if(is.infinite(max(maxprior))) max.y <- max(post[,2]) else
    max.y <- max(maxprior, max(post[,2]))
  plot(0, type="n", xlab=xlab, xlim=c(0,1), ylim=c(0, 1.1*max.y), ylab="", main="")
  polygon(post, col=grey(.8), border=grey(.4))
  .PlotMemPrior(s, mn, thick, add=TRUE, xlab="", ylab=ylab, main=main)
}



# plot the posterior (and prior) of the hiatus
.PlotHiatusPost <- function(set=get('info'), shape=set$hiatus.shape, mn=set$hiatus.mean, main="", xlim=c(0, 3*max(set$acc.mean)), xlab=paste("Hiatus size (yr)", sep=""), ylab="Frequency", minbreaks=10, after=set$after) {
  gaps <- c()
  for(i in set$hiatus.depths) {
    below <- Bacon.Age.d(i+after, set)
    above <- Bacon.Age.d(i-after, set)
    gaps <- c(gaps, below - above)
  }
  gaps <- density(gaps, from=0)
  max.x <- max(gaps$x, xlim)
  max.prior <- dgamma((shape-1)/(shape/mn), shape, shape/mn)
  if(is.infinite(max(max.prior)))
    max.y <- max(gaps$y) else
      max.y <- max(max.prior, gaps$y)
  plot(0, type="n", main="", xlab=xlab, xlim=c(0, max.x), ylab=ylab, ylim=c(0, 1.1*max.y))
  polygon(cbind(c(min(gaps$x), gaps$x, max(gaps$x)), c(0,gaps$y,0)),
    col=grey(.8), border=grey(.4))
  .PlotHiatusPrior(add=TRUE, xlim=c(0, max.x), xlab="", ylab=ylab, main=main)
}



#################### Functions for dates incl. calibration of C-14 dates ####################

#' @name copyCalibrationCurve
#' @title Copy a calibration curve.
#' @description Copy one of the the calibration curves into memory. 
#' @details Copy the radiocarbon calibration curve defined by cc into memory. 
#' @return The calibration curve (invisible).
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), 
#' \code{cc=3} for SHCal13 (southern hemisphere terrestrial). 
#' @author Maarten Blaauw, J. Andres Christen
#' @examples 
#' intcal13 <- copyCalibrationCurve(1)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @export
copyCalibrationCurve <- function(cc=1) {
  nameFile="3Col_intcal13.14C" 
  if(cc>=1 && cc <=3) {
    if(cc==1) nameFile="3Col_intcal13.14C" else
      if(cc==2) nameFile="3Col_marine13.14C" else
        if(cc==3) nameFile="3Col_shcal13.14C"
  } else
      stop("Calibration curve doesn't exist\n") 		 
  cc <- system.file("extdata", paste("Curves/", nameFile, sep=""), package='rbacon')
  cc <- read.table(cc)
  invisible(cc)
}



#' @name mix.curves
#' @title Build a custom-made, mixed calibration curve. 
#' @description If two curves need to be 'mixed' to calibrate, e.g. for dates of mixed terrestrial and marine carbon sources, then this function can be used. 
#' @details The proportional contribution of each of both calibration curves has to be set. 
#'
#' @param proportion Proportion of the first calibration curve required. e.g., change to \code{proportion=0.7} if \code{cc1} should contribute 70\% (and \code{cc2} 30\%) to the mixed curve.
#' @param cc1 The first calibration curve to be mixed. Defaults to the northern hemisphere terrestrial curve IntCal13.
#' @param cc2 The second calibration curve to be mixed. Defaults to the marine curve IntCal13.
#' @param name Name of the new calibration curve.
#' @param dirname Directory where the file will be written. If using the default \code{dirname="."}, 
#' the new curve will be saved in current working directory. 
#' @param offset Any offset and error to be applied to \code{cc2} (default 0 +- 0).
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
mix.curves <- function(proportion=.5, cc1="3Col_intcal13.14C", cc2="3Col_marine13.14C", name="mixed.14C", dirname=".", offset=c(0,0)) {
  ccloc <- paste(system.file("extdata", package='rbacon'), "/Curves/", sep="") 
  dirname <- .validateDirectoryName(dirname)
  
  cc1 <- read.table(paste(ccloc, cc1,  sep=""))
  cc2 <- read.table(paste(ccloc, cc2,  sep=""))
  cc2.mu <- approx(cc2[,1], cc2[,2], cc1[,1], rule=2)$y + offset[1] # interpolate cc2 to the calendar years of cc1
  cc2.error <- approx(cc2[,1], cc2[,3], cc1[,1], rule=2)$y
  cc2.error <- sqrt(cc2.error^2 + offset[2]^2)
  mu <- proportion * cc1[,2] + (1-proportion) * cc2.mu
  error <- proportion * cc1[,3] + (1-proportion) * cc2.error
  write.table(cbind(cc1[,1], mu, error), paste(dirname, name,  sep="") , row.names=FALSE, col.names=FALSE, sep="\t")
}



#' @name pMC.age
#' @title Calculate C14 ages from pmC values.
#' @description Calculate C14 ages from pmC values of radiocarbon dates.
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



#' @name add.date
#' @title Add dates to age-depth plots
#' @description Add dated depths to plots, e.g. to show dates that weren't used in the age-depth model
#' @details Sometimes it is useful to add additional dating information to age-depth plots, e.g., to show outliers or how dates calibrate with different estimated offsets. 
#' @param mn Reported mean of the date. Can be multiple dates. 
#' @param sdev Reported error of the date. Can be multiple dates. 
#' @param depth Depth of the date. 
#' @param cc The calibration curve to use: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), \code{cc=0} for none (dates that are already on the cal BP scale).
#' @param above Treshold for plotting of probability values. Defaults to \code{above=1e-3}.
#' @param exx Exxagaration of probability distribution plots. Defaults to \code{exx=50}.
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
#'   add.date(5000, 100, 60)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
add.date <- function(mn, sdev, depth, cc=1, above=1e-3, exx=50, normal=TRUE, normalise=TRUE, t.a=3, t.b=4, age.res=100, times=20, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5), rotate.axes=FALSE, mirror=TRUE, up=TRUE, BCAD=FALSE, pch=4) {
  if(cc > 0)
    cc = copyCalibrationCurve(cc)
  
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
  	  pol <- cbind(depth[i] + exx*c(probs, -rev(probs)), c(yrs, rev(yrs))) else
        pol <- cbind(depth[i] - up*exx*c(0, probs,  0), c(min(yrs), yrs, max(yrs)))
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
calib.plot <- function(set=get('info'), rotate.axes=FALSE, rev.d=FALSE, rev.yr=FALSE, yr.lim=c(), date.res=100, d.lab=c(), yr.lab=c(), height=15, mirror=TRUE, up=TRUE, cutoff=.001, C14.col=rgb(0,0,1,.5), C14.border=rgb(0,0,1,.75), cal.col=rgb(0,.5,.5,.5), cal.border=rgb(0,.5,.5,.75), new.plot=TRUE, plot.dists=TRUE, normalise.dists=TRUE) {
  height <- length(set$d.min:set$d.max) * height/50
  if(length(yr.lim) == 0)
    lims <- c()  
  for(i in 1:length(set$calib$probs))
    lims <- c(lims, set$calib$probs[[i]][,1])
  yr.min <- min(lims)
  yr.max <- max(lims)
  if(set$BCAD) {
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

  if(plot.dists)
    for(i in 1:length(set$calib$probs)) {
      cal <- set$calib$probs[[i]]
        d <- set$calib$d[[i]]
      if(set$BCAD)
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
      polygon(pol, col=col, border=border)
    }
}



#################### internal, user-invisible functions ####################

.validateDirectoryName <- function(dir) {
  if(!dir.exists(dir))
    dir.create(dir, recursive=TRUE)
  dir <- suppressWarnings(normalizePath(dir))
  lastchar <- substr(dir, nchar(dir), nchar(dir))
  if(lastchar != "/" & lastchar != "\\" & lastchar != "" & lastchar != "." )
    dir <- paste(dir, "/", sep="") # does this work in Windows?
  return(dir)
}



# for the proxy.ghost function
.DepthsOfScore <- function(value, dat) {
  d <- c()
  for(i in 1:(nrow(dat)-1)) {
    valueRange <- dat[i:(i+1),2]
    if(min(valueRange) <= value && max(valueRange) >= value) {
      slope <- (dat[i,2] - dat[i+1,2]) / (dat[i,1] - dat[i+1,1])
      intercept <- dat[i,2] - (slope*dat[i,1])
      if(slope==0) d[i-1] <- dat[i,1]
      d <- sort(c(d, (value - intercept) / slope ))
    }
  }
  unique(d)
}



# read the dets file, converting old formats to new ones if so required
.read.dets <- function(core, coredir, set=get('info'), sep=",", dec=".", cc=1) {
  # if a .csv file exists, read it (checking that it is more recent than any .dat file in the folder). Otherwise, read the .dat file, check the columns, report back if >4 (>5?) columns, and convert to .csv (report this also)
  csv.file <- paste(coredir,  core, "/", core, ".csv", sep="")
  dat.file <- paste(coredir,  core, "/", core, ".dat", sep="")

  dR.names <- c("r", "d", "d.r", "dr", "deltar", "r.mn", "rm", "rmn", "res.mean", "res.mn")
  dSTD.names <- c("r", "d", "d.std", "std", "std.1", "dstd", "r.std", "rstd", "res.sd")
  ta.names <- c("t", "t.a", "ta", "sta")
  tb.names <- c("t", "t.b", "tb", "stb")
  cc.names <- c("c", "cc")
  suggested.names <- c("labID", "age", "error", "depth", "cc", "dR", "dSTD", "ta", "tb")
  changed <- 0

  if(file.exists(csv.file)) {
    dets <- read.table(csv.file, header=TRUE, sep=sep)
    if(file.exists(dat.file)) # deal with old .dat files
      if(file.info(csv.file)$mtime < file.info(dat.file)$mtime)
        cat("Warning, the .dat file is newer than the .csv file! I will read the .csv file. From now on please modify ", csv.file, ", not ", dat.file, " \n", sep="") else
          cat("Reading", csv.file, "\n")
    } else {
      cat("No .csv file found, reading", dat.file, "and converting it to .csv\n")
      dets <- read.table(dat.file, header=TRUE)
      changed <- 1
    }
  name <- tolower(names(dets))

  # check if 'classic' dets file, which has a different column order from the current default
  if(ncol(dets) > 4)
    if(ncol(dets) == 5) { # then probably a 'new' dets file
      if((name[5] %in% cc.names) && min(dets[,5]) >= 0 && max(dets[,5]) <= 4) {} else# extra check for correct values
        stop("Error! Unexpected name or values in fifth column (cc, should be between 0 and 4). Please check the manual for guidelines in producing a correct .csv file.\n")
    } else
      if(ncol(dets) == 6) { # probably an 'old' file: dR, dSTD
        if(name[5] %in% dR.names && name[6] %in% dSTD.names) {
          dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
          cat(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file, "\n")
          changed <- 1
        }
      } else
        if(ncol(dets) == 7) { # probably a 'new' file: cc, dR, dSTD
          if(name[5] %in% cc.names && min(dets[,5]) >= 0 && max(dets[,5]) <= 4 &&
            name[6] %in% dR.names && name[7] %in% dSTD.names) 
              {} else
                stop("Error! Unexpected column names, order or values in dets file. \nPlease check the manual for correct dets file formats.\n")
        } else
          if(ncol(dets) == 8) { # probably an 'old' file: dR, dSTD, ta, tb
            if(name[5] %in% dR.names && name[6] %in% dSTD.names)
            if(name[7] %in% ta.names && name[8] %in% tb.names)
            if(range(dets[,8] - dets[,7]) == c(1,1)) { # check that these set expected student-t values
              dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
              cat(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file, "\n")
              changed <- 1
            } else
              stop("Error! Unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file")
          } else
            if(ncol(dets) == 9) { # most complex case, many checks needed
              if(name[9] %in% cc.names && # we're almost sure that this is a 'classic' dets file
                min(dets[,9]) >= 0 && max(dets[,9]) <= 4 && # check that this sets calibration curves
                  range(dets[,8] - dets[,7]) == c(1,1) && # check that these set expected student-t values
                    name[5] %in% dR.names && name[6] %in% dSTD.names && # column names as expected?
                      name[7] %in% ta.names && name[8] %in% tb.names) { # column names as expected?
                        dets <- dets[,c(1:4,9,5:8)] # shuffle colums around
                        cat(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD, t.a, t.b, cc. \nAdapting column order and saving as", csv.file, "\n")
                        changed <- 1
                      } else
                        if(name[5] %in% cc.names && # oh, probably a 'new' file from more recent Bacon
                          min(dets[,5]) >= 0 && max(dets[,5]) <= 4 && # check that this sets cal.curves
                            range(dets[,9] - dets[,8]) == c(1,1) && # columns 8-9 set student-t correctly
                              name[8] %in% ta.names && name[9] %in% tb.names && # and are correctly named
                                name[6] %in% dR.names && name[7] %in% dSTD.names) # all lights are green
                                  {} else
                                     stop("Error! Unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file")
            } else
              stop("Error! Unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.\n")

  # more sanity checks
  if(!is.numeric(dets[,2]) || !is.numeric(dets[,3]) || !is.numeric(dets[,4]))
    stop("Error, unexpected values in dets file, I expected numbers. Check the manual.\n", call.=FALSE)
  if(min(dets[,3]) <= 0) {
    cat("Warning, zero year errors don't exist in Bacon's world. I will increase them to 1yr.\n")
    dets[dets[,3] <= 0,3] <- 1
    changed <- 1
  }
  if(min(diff(dets[,4]) < 0)) {
    cat("\nWarning, the depths are not in ascending order, I will correct this")
    dets <- dets[order(set$dets[,4]),]
    changed <- 1
  }

  # if current dets differ from original .csv file, rewrite it
  if(changed > 0)
    write.table(dets, csv.file, sep=paste(sep, "\t", sep=""), dec=dec, row.names=FALSE, col.names=suggested.names[1:ncol(dets)], quote=FALSE)
  dets
}



# read in default values, values from previous run, any specified values, and report the desired one. Internal function. 
.Bacon.settings <- function(core, coredir, dets, thick, remember=TRUE, d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, hiatus.depths, hiatus.mean, hiatus.shape, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, unit, normal, t.a, t.b, deltaR, deltaSTD, prob, defaults, runname, ssize, dark, MinYr, MaxYr, cutoff, yr.res, after, find.round) {
  
  vals <- list(d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, hiatus.depths, hiatus.mean, hiatus.shape, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, unit, normal, t.a, t.b, deltaR, deltaSTD, prob)
  valnames <- c("d.min", "d.max", "d.by", "depths.file", "slump", "acc.mean", "acc.shape", "mem.mean", "mem.strength", "hiatus.depths", "hiatus.mean", "hiatus.shape", "BCAD", "cc", "postbomb", "cc1", "cc2", "cc3", "cc4", "unit", "normal", "t.a", "t.b", "deltaR", "deltaSTD", "prob")

  extr <- function(i, def=deffile, pre=prevfile, exists.pre=prevf, rem=remember, sep=" ", isnum=TRUE) {
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
  prevfile <- paste(coredir, core, "/", core, "_settings.txt", sep="")
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
  hiatus.depths <- if(is.na(hiatus.depths)[1]) NA else extr(10)
  hiatus.mean <- extr(11); hiatus.shape <- extr(12)
  BCAD <- extr(13); cc <- extr(14); postbomb <- extr(15); cc1 <- extr(16, isnum=FALSE)
  cc2 <- extr(17, isnum=FALSE); cc3 <- extr(18, isnum=FALSE); cc4 <- extr(19, isnum=FALSE)
  unit <- extr(20, isnum=FALSE); normal <- extr(21); t.a <- extr(22); t.b <- extr(23)
  deltaR <- extr(24); deltaSTD <- extr(25); prob <- extr(26)

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

  ### produce/update settings file, and return the values
  prevfile <- file(paste(coredir, core, "/", core, "_settings.txt", sep=""), "w")
  scat <- function(m, n="") cat(m, n, sep="", file=prevfile)
  cat(d.min, " #d.min\n", d.max, " #d.max\n", d.by, " #d.by\n",
    depths.file, " #depths.file\n", slump, " #slump\n", sep="", file=prevfile)
  for(i in acc.mean) scat(i, " "); scat("#acc.mean\n")
  for(i in acc.shape) scat(i, " "); scat("#acc.shape\n", "")
  for(i in mem.mean) scat(i, " "); scat("#mem.mean\n", "")
  for(i in mem.strength) scat(i, " "); scat("#mem.strength\n", "")
  for(i in hiatus.depths) scat(i, " "); scat("#hiatus.depths\n", "")
  for(i in hiatus.mean) scat(i, " "); scat("#hiatus.mean\n", "")
  for(i in hiatus.shape) scat(i, " "); scat("#hiatus.shape\n", "")
  cat(BCAD, " #BCAD\n", cc, " #cc\n", postbomb, " #postbomb\n",
    cc1, " #cc1\n", cc2, " #cc2\n", cc3, " #cc3\n", cc4, " #cc4\n",
    unit, " #unit\n", normal, " #normal\n", t.a, " #t.a\n", t.b, " #t.b\n",
    deltaR, " #deltaR\n", deltaSTD, " #d.STD\n", prob, " #prob\n", sep="", file=prevfile)
  close(prevfile)

  if(length(MinYr) == 0)
    MinYr <- min(-1e3, round(dets[,2] - (5*dets[,3])))
  if(length(MaxYr) == 0)
    MaxYr <- max(1e6, round(dets[,2] + (5*dets[,3])))

  list(core=core, thick=thick, dets=dets, d.min=d.min, d.max=d.max,
    d.by=d.by, depths.file=depths.file, slump=slump,
    acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean,
    mem.strength=mem.strength, hiatus.depths=hiatus.depths, hiatus.mean=hiatus.mean,
    hiatus.shape=hiatus.shape, BCAD=BCAD, cc=cc, postbomb=postbomb,
    cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, unit=noquote(unit), normal=normal,
    t.a=t.a, t.b=t.b, deltaR=deltaR, deltaSTD=deltaSTD, prob=prob, date=date(),
    runname=runname, ssize=ssize, dark=dark, MinYr=MinYr, MaxYr=MaxYr, 
    cutoff=cutoff, yr.res=yr.res, after=after, find.round=find.round)
}



# write files to be read by the main Bacon age-depth modelling function
.write.Bacon.file <- function(set=get('info')) {
  if(set$d.min < min(set$dets[,4])) {
    extrap <- c(NA, min(set$dets[,2]), max(1e4, 100*set$dets[,3]), set$d.min, 0)
    dets <- rbind(extrap, set$dets, deparse.level=0)
  }
  if(set$d.max > max(set$dets[,4])) {
    extrap <- c(NA, max(set$dets[,2]), max(1e4, 100*set$dets[,3]), set$d.max, 0)
    dets <- rbind(set$dets, extrap, deparse.level=1)
  }

  fl <- file(set$bacon.file, "w")
  cat("## Ran on", set$date, "\n\n", file=fl)
  cat("Cal 0 : ConstCal;\nCal 1 : ",
  if(set$cc1=="IntCal13" || set$cc1=="\"IntCal13\"") "IntCal13" else noquote(set$cc1),
  ", ", set$postbomb, ";\nCal 2 : ",
  if(set$cc2=="Marine13" || set$cc2=="\"Marine13\"") "Marine13" else noquote(set$cc2),
  ";\nCal 3 : ",
  if(set$cc3=="SHCal13" || set$cc3=="\"SHCal13\"") "SHCal13" else noquote(set$cc3), ", ", set$postbomb, ";",
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") set$cc4 <- c() else
  # if(.Platform$OS == "windows")
  #  paste("\nCal 4 : GenericCal, Curves\\", set$cc4, ";", sep="") else
  paste("\nCal 4 : GenericCal, ", set$cc4, ";", sep=""), sep="", file=fl)
  cat("\n\n##   id.   yr    std   depth  deltaR  deltaSTD     t.a   t.b   cc", file=fl)

  if(ncol(set$dets) == 4) { # then we need to provide some constants once only
    cat("\nDet 0 : ", as.character(set$dets[1,1]), " ,  ", set$dets[1,2], ",  ",
      set$dets[1,3], ",  ", set$dets[1,4], ",  ", set$deltaR, ",  ", set$deltaSTD,
      ",  ", set$t.a, ",  ", set$t.b, ",  ", set$cc, ";", sep="", file=fl)
      if(nrow(set$dets)>1)
        for(i in 2:nrow(set$dets))
          cat("\nDet ", i-1, " : ",  as.character(set$dets[i,1]),
          " , ", set$dets[i,2], ", ", set$dets[i,3], ", ", set$dets[i,4],
          ";", sep="", file=fl)
  } else { # use additional columns provided within the .dat file
      cc <- set$dets[,5]
      deltaR <- rep(set$deltaR, nrow(set$dets))
      deltaR[cc==0] <- 0 # only apply dR to C14 dates
      deltaSTD <- rep(set$deltaSTD, nrow(set$dets))
      deltaSTD[cc==0] <- 0 # only apply dR to C14 dates
      t.a <- rep(set$t.a, nrow(set$dets))
      t.b <- rep(set$t.b, nrow(set$dets))

      if(ncol(set$dets) >= 7) {
        deltaR <- set$dets[,6]
        deltaSTD <- set$dets[,7]
      }
      if(ncol(set$dets) >= 9) {
        t.a <- set$dets[,8]
        t.b <- set$dets[,9]
      }

      for(i in 1:nrow(set$dets))
        cat("\nDet ", i-1, " : ",  as.character(set$dets[i,1]), " , ",
          set$dets[i,2], ", ", set$dets[i,3], ", ", set$dets[i,4],  ",  ",
          deltaR[i], ",  ", deltaSTD[i], ",  ", t.a[i], ",  ", t.b[i], ",  ",
          cc[i], ";", sep="", file=fl)
    }

  if(!is.na(set$hiatus.depths)[1]) { ### hiatus(es) inferred. Prior values should be provided starting from top section
    cat("\n  Hiatus set at depth(s)", set$hiatus.depths, "\n")
      if(length(set$acc.shape)==1)
        set$acc.shape <- rep(set$acc.shape, length(set$hiatus.depths)+1)
      if(length(set$acc.mean)==1)
        set$acc.mean <- rep(set$acc.mean, length(set$hiatus.depths)+1)
      if(length(set$hiatus.mean)==1)
        set$hiatus.mean <- rep(set$hiatus.mean, length(set$hiatus.depths))
      if(length(set$hiatus.shape)==1)
        set$hiatus.shape <- rep(set$hiatus.shape, length(set$hiatus.depths))
      .assign_to_global ("info", set)
      cat("\n\n### Depths and priors for fixed hiatuses, in descending order",
        "\n##### cm  alpha beta      ha     hb", file=fl)
      for(i in length(set$hiatus.depths):1)
        cat("\nHiatus ", i-1, ":  ", set$hiatus.depth[i], ",  ", set$acc.shape[i+1],
          ",  ", set$acc.shape[i+1]/set$acc.mean[i+1], ",  ", set$hiatus.shape[i],
          ",  ", set$hiatus.shape[i]/set$hiatus.mean[i], ";", sep="", file=fl)
  }

  ### final parameters
  wrapup <- paste("\n\n##\t\t K   MinYr   MaxYr   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax",
    "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
    ",  ", set$MinYr, ",  ", set$MaxYr, ",  ", set$th0[1], ",  ", set$th0[2],
    ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
    ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
    ", ", set$d.max, ";\n", sep="")
  cat(wrapup, file=fl)
  close(fl)
}



# function to read output files into memory
.Bacon.AnaOut <- function(fnam, set=get('info')) {
  out <- read.table(fnam)
  n <- ncol(out)-1
  set$n <- n
  set$Tr <- nrow(out)
  set$Us <- out[,n+1]
  set$output <- out[,1:n]
  set
}



# function to load results in global environment
# parameter position defaults to 1, which equals an assignment to the global environment
.assign_to_global <- function(key, val, pos=1) {
  assign(key, val, envir=as.environment(pos) )
}
