# in the inst/dev/ folder, there is now a testBaconplots.Rmd function which automates plotting and checking many functions. There is also a file render-plots.yml which can be used to test many plots on a range of github systems (ubuntu, fedora and windows). Produced html files can be downloaded and checked locally. To do this, the file has to be placed in .github/workflows/.

# Check if we can/should return to using a gamma distribution instead of a uniform one for the hiatus

# make a function to include e.g. cumulative weight/pollen instead of depths - 'fake' depths

# do: check that overlap function continues to function (sometimes reports 0% overlap when the dates fit well), check rplum bugs w youngest.age (is the bug in rbacon or in rplum?) and w larger-than-previous error sizes

# replacing the plotting of the calibrated distributions by rice's functions doesn't seem to speed up anything, so keeping the original method in place for now.

# for future versions: add function to estimate best thick value, check if a less ugly solution can be found to internal_plots.R at line 26 (hists length < 7). This happens when there are some very precise dates causing non-creation of th0/th1, investigate the slowness of plotting after the Bacon run (not only dates, also the model's 95% ranges etc.), produce proxy.ghost graph with proxy uncertainties?, check/adapt behaviour of AgesOfEvents around hiatuses, if hiatus or boundary plot acc.posts of the individual sections?, allow for asymmetric cal BP errors (e.g. read from files), proxy.ghost very slow with long/detailed cores - optimization possible?, check again if/how/when Bacon gets confused by Windows usernames with non-ascii characters (works fine on Mac; use normalizePath or other R-based solutions)

# read https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Registering-native-routines for linking between rbacon and rplum. Currently done using utils::getFromNamespace which is basically a way to allow :::

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
#' The accumulation rates are constrained by prior information on the accumulation rate (\code{acc.mean, acc.shape}) and its
#' variability between neighbouring depths, or "memory" (\code{mem.mean, mem.strength}). Hiatuses can be introduced as well, also constrained by prior information (\code{hiatus.max}).
#'
#' Although Bacon works with any kind of absolute dates (e.g., OSL, tephra or other dates on a calendar scale),
#' it is often used to age-model 14C-dated sequences. Radiocarbon dates should be calibrated using either IntCal20
#' (for terrestrial northern hemisphere material; Reimer et al., 2020), Marine20 (for marine dates; Hughen et al., 2020),
#' SHCal20 (for southern hemisphere dates; Hogg et al., 2020) or any other calibration curve (see below), while modern 14C
#' dates are calibrated using one of the post-bomb calibration curves (NH1, NH2 or NH3 for the northern hemisphere,
#' SH1-2 or SH3 for the southern hemisphere; Hua et al., 2022). See \url{http://calib.org/CALIBomb/} if you are unsure which
#' postbomb curve you need. If Bacon finds postbomb dates (negative 14C ages) and you haven't specified a postbomb curve,
#' you will be prompted. Provide postbomb curves as, e.g., \code{postbomb=1} for the NH1 postbomb curve (2 for NH2, 3 for NH3, 4 for SH1-2, 5 for SH3).
#' For calendar dates, i.e. dates that are already on the calendar scale and thus should not be calibrated, set\code{cc=0}.
#'
#' Since version 3.1.0, rbacon can also handle younger-than and older-than ages, with the model aiming to either go 'above'
#' or 'below' such dates as requested. If the resulting combination of parameters becomes problematic (e.g., no initial
#' combination of parameters can be found that obeys the priors or is in chronological order), then the output will often be wrong.
#' If so, using the function \link{set.initvals} could help.
#'
#' By default, the initial MCMC values of the Bacon age-depth model (upper ages and accumulation rate for each model section)
#' are estimated randomly. Since version 3.1.0, these starting values can also be provided in a file with extension _bacon.init,
#' placed within the core's folder. This file will need to have two rows, each for one of the two initial sets of parameters required
#' (the t-walk requires two starting estimates for all MCMC parameters).
#' If such a file is found (and correctly formatted), Bacon will use the values within this file
#' as starting points for the MCMC run. See function \link{set.initvals} for more information.
#'
#' From version 2.5.1 on (i.e., since February 2021), the default memory prior has changed to \code{mem.mean=0.5}
#' and \code{mem.strength=10}. Previously used c++ code contained a bug which caused the prior information for the memory not to be
#' taken into account correctly. Now that this bug has been repaired, the default memory prior has been updated such that it should work
#' for most types of cores, and should result in similar output to previous versions of Bacon. There is no need to re-do previous runs.
#' However, it is considered good practice to test the impact of different settings on a site's age-depth model
#' (e.g., thick, acc.mean, acc.shape, mem.mean, acc.strength).

#' @param core Name of the core, given using quotes. Defaults to one of the cores provided with rbacon, \code{core="MSB2K"}.
#' An alternative core provided with this package is RLGH3 (Jones et al., 1989).
#' To run your own core, produce a .csv file with the dates as outlined in the manual, add a folder with the core's name to the default directory for cores (see \code{coredir}), and save the .csv file there. For example, the file's location and name could be \code{Bacon_runs/MyCore/MyCore.csv}. Then run Bacon as follows: \code{Bacon("MyCore")}
#' @param thick Bacon will divide the core into sections of equal thickness specified by thick (default \code{thick=5}).
#' @param coredir Folder where the core's files \code{core} are and/or will be located. This will be a folder with the core's name, within either the folder \code{coredir='Bacon_runs/'}, or the folder Cores/ if it already exists within R's working directory, or a custom-built folder. For example, use \code{coredir="."} to place the core's folder within the current working directory, or \code{coredir="F:"} if you want to put the core's folder and files on a USB drive loaded under F:.
#' Thinner (and thus more) sections will result in smoother age-models, but too many sections can cause `run-away' models.
#' @param prob Confidence interval to report. This should lie between 0 and 1, default 0.95 (95 \%).
#' @param d.min Minimum depth of age-depth model (use this to extrapolate to depths higher than the top dated depth).
#' @param d.max Maximum depth of age-depth model (use this to extrapolate to depths below the bottom dated depth).
#' @param add.bottom Add a model section at the bottom of the core, in order to ensure the bottommost date is taken into account. Default \code{add.bottom=TRUE}. This is a new option and can cause age-models to differ from previous version. Please re-run the model if in doubt.
#' @param d.by Depth intervals at which ages are calculated. Defaults to \code{d.by=1}. Please ensure that the value of d.by is smaller than that of 'thick', otherwise plots might turn out wrong.
#' @param seed Seed used for C++ executions. If it is not assigned (\code{seed=NA}; default) then the seed is set by system.
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
#'  but its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.strength=10}
#'  (higher values result in more peaked shapes) allow for a large range of posterior memory values. Please note that the default memory prior has been updated from rbacon version 2.5.1 on, to repair a bug. 
#' @param mem.mean The prior for the memory is a beta distribution, which looks much like the gamma distribution but
#' its values are always between 0 (no assumed memory) and 1 (100\% memory). Its default settings of \code{mem.mean=0.5}
#' allow for a large range of posterior memory values. Please note that the default memory prior has been updated from rbacon version 2.5.1. on, to repair a bug. 
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
#' @param cc4 Provide the name of an alternative curve (3 columns: cal BP, 14C age, error, separated by white spaces and saved as a plain-text file). It is important here to first produce a tailor-made folder for your and the default calibration curves to live in. See \code{cc.dir}. Defaults to \code{cc4="mixed.14C"}. 
#' @param cc.dir Directory where the calibration curves for C14 dates \code{cc} are located. By default uses the location of the rintcal package which provides the calibration curves. If you want to use custom-made calibration curves, first set up a new folder using the function new.ccdir() in the rintcal package, e.g., \code{new.ccdir="MyCurves"}, then place the custom curve in that folder using \code{rintcal::mix.ccurves(, cc.dir="MyCurves", save=TRUE)}.
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}
#' @param F14C Radiocarbon ages can be provided as F14C values. If doing so, please indicate here which dates were entered as F14C (e.g., if the first 4 dates are in F14C, write \code{F14C=1:4}). The F14C values in your .csv file will then be replaced by their corresponding C14 ages.
#' @param pMC Radiocarbon ages can be provided as pMC values. If doing so, please indicate here which dates were entered as pMC (e.g., if the first 4 dates are in pMC, write \code{pMC=1:4}). The pMC values in your .csv file will then be replaced by their corresponding C14 ages.
#' @param delta.R Mean of core-wide age offsets (e.g., regional marine offsets).
#' @param delta.STD Error of core-wide age offsets (e.g., regional marine offsets).
#' @param t.a The dates are treated using the t distribution (Christen and Perez 2009) by default (\code{normal=FALSE}).
#' This t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param t.b The dates are treated using t distribution by default (\code{normal=FALSE}).
#' The t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2009).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param normal By default, Bacon uses the t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param suggest If initial analysis of the data indicates abnormally slow or fast accumulation rates, Bacon will suggest to change the prior.
#' @param accept.suggestions Automatically accept the suggested values. Use with care. Default \code{accept.suggestions=FALSE}.
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
#' @param remove Whether or not to remove depths within slumps. Defaults to \code{remove=FALSE}.
#' @param BCAD The calendar scale of graphs and age output-files is in cal BP (calendar or calibrated years before the present, where the present is AD 1950) by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @param ssize The amount of iterations to store at the end of the MCMC run. Default 4000; decrease for faster (but less reliable) runs or increase for cores where the MCMC mixing (panel at upper-left corner of age-model graph) appears problematic.
#' @param th0 Starting years for the MCMC iterations. These are randomly chosen by default.
#' @param burnin Amount of initial, likely sub-optimal MCMC iterations that will be removed.
#' @param MinAge Deprecated - use youngest.age instead.
#' @param youngest.age Minimum age limit for Bacon runs, default at current year in cal BP. To set plot limits, use \code{age.min} instead.
#' @param MaxAge Deprecated - use oldest.age instead.
#' @param oldest.age Maximum age limit for Bacon runs, default at 1,000,000 cal BP. To set plot limits, use \code{age.max} instead.
#' @param cutoff Avoid plotting very low probabilities of date distributions (default \code{cutoff=0.001}).
#' @param plot.pdf Produce a pdf file of the age-depth plot. Defaults to \code{plot.pdf=TRUE} after a Bacon run.
#' @param dark Darkness of the greyscale age-depth model. The darkest grey value is \code{dark=1} by default.
#' Lower values will result in lighter grey but values >1 are not allowed.
#' @param date.res Date distributions are plotted using \code{date.res=100} segments by default.
#' @param age.res Resolution or amount of greyscale pixels to cover the age scale of the age-model plot. Default \code{yr.res=200}.
#' @param yr.res Deprecated - use age.res instead
#' @param close.connections Internal option to close connections after a run. Default \code{close.connections=TRUE}.
#' @param save.info By default, a variable called `info' with relevant information about the run (e.g., core name, priors, settings, ages, output) is saved into the working directory. Note that this will overwrite any existing variable with the same name - as an alternative, one could run, e.g., \code{myvar <- Bacon()}, followed by supplying the variable \code{myvar} in any subsequent commands.
#' @param older.than an option to enable dates at the limit of C-14 dating. If there are older.than dates, they tell us that the core should be older than a certain age at that depth. For example, if the 7th and 8th dates in the core's .csv file are older-than dates, use as \code{older.than=c(7,8)}. The MCMC run could be problematic if the older-than ages do not fit with the other information.
#' @param younger.than an option to provide younger-than ages, for example a historical pollen marker. If there are younger-than dates, they tell us that the core should be younger than a certain age at that depth. For example, if the 7th and 8th dates in the core's .csv file are younger.than dates, use as \code{younger.than=c(7,8)}. The MCMC run could be problematic if the younger.than ages do not fit with the other information.
#' @param save.elbowages If you want to have a file with the MCMC-derived ages for all the age-depth model's elbows, set \code{save.elbowages=TRUE} and a file with the ages will be saved in the core's folder, starting with the core name, followed by its number of sections, d.min, and section thickness, and ending in "_elbowages.txt".
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
#' @references
#' Blaauw, M. and Christen, J.A., 2011. Flexible paleoclimate age-depth models using an autoregressive gamma process. Bayesian Anal. 6(3): 457-474.
#'
#' Christen, J.A., Perez E., S., 2010. A new robust statistical model for radiocarbon data. Radiocarbon 51: 1047-1059.
#'
#' Reimer et al., 2020. The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0–55 cal kBP). Radiocarbon 62(4): 725-757. \doi{10.1017/RDC.2020.41}
#'
#' Hogg et al. 2020 SHCal20 Southern Hemisphere calibration, 0-55,000 years cal BP. Radiocarbon 62(4): 759-778. \doi{10.1017/RDC.2020.59}
#'
#' Hughen et al. 2020 Marine20-the marine radiocarbon age calibration curve (0-55,000 cal BP). Radiocarbon 62(4): 779-820.  \doi{10.1017/RDC.2020.68}
#'
#' Hua et al. 2022 Atmospheric radiocarbon for the period 1950-2019. Radiocarbon 64(4), 723-745, \doi{10.1017/RDC.2021.95}
#'
#' Jones, V.J., Stevenson, A.C., Battarbee, R.W., 1989. Acidification of lakes in Galloway, south west Scotland
#' - a diatom and pollen study of the post-glacial history of the Round Loch of Glenhead.
#' Journal of Ecology 77: 1-23.
#'
#' @export
Bacon <- function(core="MSB2K", thick=5, coredir="", prob=0.95, d.min=NA, d.max=NA, add.bottom=TRUE, d.by=1, seed=NA, depths.file=FALSE, depths=c(), depth.unit="cm", age.unit="yr", unit=depth.unit, acc.shape=1.5, acc.mean=20, mem.strength=10, mem.mean=0.5, boundary=NA, hiatus.depths=NA, hiatus.max=10000, add=c(), after=.0001/thick, cc=1, cc1="IntCal20", cc2="Marine20", cc3="SHCal20", cc4="ConstCal", cc.dir=c(), postbomb=0, F14C=c(), pMC=c(), delta.R=0, delta.STD=0, t.a=3, t.b=4, normal=FALSE, suggest=TRUE, accept.suggestions=FALSE, reswarn=c(10,200), remember=TRUE, ask=TRUE, run=TRUE, defaults="defaultBacon_settings.txt", sep=",", dec=".", runname="", slump=c(), remove=FALSE, BCAD=FALSE, ssize=4000, th0=c(), burnin=min(500, ssize), youngest.age=c(), oldest.age=c(), MinAge=c(), MaxAge=c(), cutoff=.01, plot.pdf=TRUE, dark=1, date.res=100, age.res=200, yr.res=age.res, close.connections=TRUE, save.info=TRUE, older.than=c(), younger.than=c(), save.elbowages=FALSE, verbose=TRUE, ...) {
  # Check coredir and if required, copy example files into core directory
  coredir <- assign_coredir(coredir, core, ask, isPlum=FALSE)
  if(core == "MSB2K" || core == "RLGH3") {
    if(!dir.exists(file.path(coredir, core))) {
      dir.create(file.path(coredir, core), showWarnings = FALSE, recursive = TRUE)
      fileCopy <- system.file(paste0("extdata/Cores/", core), package="rbacon")
      file.copy(fileCopy, coredir, recursive = TRUE, overwrite=FALSE)
      }
    }

  # set the calibration curve
  if(length(cc.dir) == 0)
    cc.dir <- system.file("extdata", package="rintcal")
  cc.dir <- validateDirectoryName(cc.dir)

  # default_settings.txt is located within system.file
  defaults <- system.file("extdata", defaults, package=packageName())
  # read in the data, adapt settings from defaults if needed
  dets <- read.dets(core, coredir, sep=sep, dec=dec, cc=cc)
  # give feedback about calibration curves used
  if(ncol(dets) > 4 && length(cc) > 0) {
    cc.csv <- unique(dets[,5])
	if(verbose) {
      if(length(cc.csv) == 1) {
        if(cc.csv != cc)
          message(" Using calibration curve specified within the .csv file,", cc.csv, "\n")
      } else
        if(min(cc.csv) == 0)
          message(" Using a mix of cal BP and calibrated C-14 dates\n")
        else
          message(" Using several C-14 calibration curves\n")
    }
  }

  # Oct 2024
  if(length(F14C) > 0) {
	if(min(dets[F14C,2]) < 0 || max(dets[F14C,2]) > 3) 
      stop("The F14C values cannot be negative and are unlikely to be >3. Are you sure these values are in F14C?")		
    asC14 <- F14CtoC14(dets[F14C,2], dets[F14C,3])
	dets[F14C,2] <- asC14[,1]
	dets[F14C,3] <- asC14[,2]
    csv.file <- paste0(coredir, core, "/", core, ".csv")
	fastwrite(as.data.frame(dets), csv.file, sep=sep, dec=dec, row.names=FALSE, quote=FALSE) 
	message(paste("replaced F14C values with C14 ages in", csv.file))  
  }
  if(length(pMC) > 0) {
	if(min(dets[pMC,2]) < 0 || max(dets[pMC,2]) > 300) 
      stop("The pMC values cannot be negative and are unlikely to be >300. Are you sure these values are in pMC?")		
    asC14 <- pMCtoC14(dets[pMC,2], dets[pMC,3])
	dets[pMC,2] <- asC14[,1]
	dets[pMC,3] <- asC14[,2]
    csv.file <- paste0(coredir, core, "/", core, ".csv")
	fastwrite(as.data.frame(dets), csv.file, sep=sep, dec=dec, row.names=FALSE, quote=FALSE) 
	message(paste("replaced pMC values with C14 ages in", csv.file))  
  }

  if(suggest) { # adapt prior for mean accumulation rate?
    sugg <- sapply(c(1,2,5), function(x) x*10^(-1:2)) # some suggested "round" values
    ballpacc <- lm(dets[,2]*1.1 ~ dets[,4])$coefficients[2] # very rough acc.rate estimates, uncalibrated dates
    ballpacc <- abs(sugg - ballpacc) # get absolute differences between given acc.mean and suggested ones
    ballpacc <- ballpacc[ballpacc > 0] # do not suggest 0
    sugg <- sugg[order(ballpacc)[1]] # suggest rounded acc.rate with lowest absolute difference
    if(!sugg %in% acc.mean) 
      if(accept.suggestions) { # new Oct '20
        acc.mean <- sugg
        message("Adapting acc.mean to ", sugg, " ", age.unit, "/", depth.unit)
    } else {
        ans <- readline(message(" Ballpark estimates suggest changing the prior for acc.mean to ", sugg, " ", age.unit, "/", depth.unit, ". OK? (y/N) "))
        if(tolower(substr(ans,1,1)) == "y")
          acc.mean <- sugg else
            message(" No problem, using the provided prior")
      }
    }

  if(thick < d.by)
    warning("Please set d.by to a value smaller than that of thick", .call=TRUE)

  # check values for the prior's mean, Jan 2021
  if(mem.mean < 0 || mem.mean >1)
    stop("The prior for the mean of the memory should be between 0 and 1", call.=FALSE)
  if(length(mem.mean) > 1)
    stop("Can only use one value for mem.mean across a core", call.=FALSE)
  if(length(mem.strength) > 1)
    stop("Can only use one value for mem.strength across a core", call.=FALSE)
    
  if(!is.na(boundary[1])) {
    boundary <- sort(unique(boundary)) 
    if(length(acc.mean) == 1) # August 2024
      acc.mean <- rep(acc.mean, length(boundary)+1)	
  }
  if(!is.na(hiatus.depths[1])) {
    hiatus.depths <- sort(unique(hiatus.depths))
    if(length(acc.mean) == 1) # why not for boundary?
      acc.mean <- rep(acc.mean, length(hiatus.depths)+1)
  }

  # set reasonable boundaries for the ages
#  if(length(MinAge) == 0)
#    MinAge <- min(1950 - as.integer(format(Sys.time(), "%Y")))#, round(dets[,2] - (5*dets[,3])))
#  if(length(MaxAge) == 0)
#    MaxAge <- max(1e6, round(dets[,2] + (5*dets[,3])))

  info <- Bacon.settings(core=core, coredir=coredir, dets=dets, thick=thick, remember=remember, d.min=d.min, d.max=d.max, d.by=d.by, depths.file=depths.file, slump=slump, acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean, mem.strength=mem.strength, boundary=boundary, hiatus.depths=hiatus.depths, hiatus.max=hiatus.max, BCAD=BCAD, cc=cc, postbomb=postbomb, cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=depth.unit, normal=normal, t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, defaults=defaults, runname=runname, ssize=ssize, dark=dark, youngest.age=youngest.age, oldest.age=oldest.age, cutoff=cutoff, age.res=age.res, after=after, age.unit=age.unit)
  
  # optionally, make the info variable available in the working environment (default, but will overwrite any existing variable with the name 'info')
  info$save.info <- save.info
  if(save.info) 
    assign_to_global("info", info)

  info$coredir <- coredir
  if(is.na(seed))
    seed <-sample(1:1e6, 1) # sample an integer
  set.seed(seed) # Nov 2020
  info$seed <- seed
  info$isplum <- FALSE

  ### check for initial mistakes
  if(length(MinAge) > 0)
    warning("please do not use MinAge, instead use the option youngest.age", call.=FALSE)
  if(length(MaxAge) > 0)
    warning("please do not use MaxAge, instead use the option oldest.age", call.=FALSE)
  #if(any(info$acc.shape == info$acc.mean))
  #  stop("acc.shape cannot be equal to acc.mean", call.=FALSE)
  if(info$t.b - info$t.a != 1)
    warning("t.b - t.a should always be 1, check the manual", call.=FALSE)
  if(min(acc.shape) < 1)
    warning("\nWarning, using values <1 for acc.shape might cause unexpected results\n", call.=TRUE)

  ### calibrate dates
  negativeages <- FALSE
  if(info$postbomb == 0) # no postbomb curve selected
    if(ncol(info$dets) == 4) { # not much info on which calcurve to use
      if(info$cc > 0) # we are using C-14 dates
        if(min(info$dets[,2]) < 0) # negative C-14 ages
          negativeages <- TRUE
    } else # then we have a cc column, which is helpful
        if(max(info$dets[,5]) > 0) # using radiocarbon dates
          if(min(info$dets[which(info$dets[,5] > 0),2]) < 0) # negative C14 ages
            negativeages <- TRUE
  if(negativeages)
    stop("you have negative C14 ages so should select a postbomb curve", call.=FALSE)
  info$calib <- bacon.calib(dets, info, date.res, cc.dir=cc.dir, cutoff=cutoff)
  
  ### find some relevant values
  info$rng <- c()
  for(i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[,1]) # removed cutoff selection from tmp
  }

  if(length(th0) == 0) # provide two ball-park/initial age estimates
    info$th0 <- round(rnorm(2, max(youngest.age, dets[1,2]), dets[1,3]))
  info$th0[info$th0 < info$youngest.age] <- info$youngest.age # otherwise twalk will not start

  ### assign depths
  if(length(depths) == 0)
    depths <- seq(info$d.min, info$d.max, by=d.by) # was info$d.by
  if(depths.file) {
    dfile <- paste0(info$coredir, info$core, "/", info$core, "_depths.txt")
    if(!file.exists(dfile))
      stop("I cannot find the file ", paste0(info$coredir, info$core, "/", info$core, "_depths.txt"), call.=FALSE)
    depths <- fastread(dfile, header=FALSE)[,1]
    if(!is.numeric(depths[1]))
      stop("File should contain numbers only, no headers", call.=FALSE)
  }
  info$depths <- depths
  if(min(depths) < info$d.min)
    info$d.min <- min(depths)
  if(max(depths) > info$d.max)
    info$d.max <- max(depths)

  info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)
  if(add.bottom)  # new October 2020
    info$elbows <- c(info$elbows, max(info$elbows)+thick) # new October 2020
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
          if(accept.suggestions) 
            ans <- "y" else 
              ans <- readline(message(" Warning, the current value for thick, ", thick, ", will result in very few age-model sections (", info$K, ", not very flexible). Suggested maximum value for thick: ", sugg, " OK? (y/n) "))
        } else
          if(info$K > max(reswarn)) {
            sugg <- max(pretty(thick*(info$K/max(reswarn))))
            if(accept.suggestions) 
              ans <- "y" else
                ans <- readline(message(" Warning, the current value for thick, ", thick, ", will result in very many age-model sections (", info$K, ", possibly hard to run). Suggested minimum value for thick: ", sugg, " OK? (y/n) "))
          }
    if(tolower(substr(ans, 1, 1)) == "y") {
      message(" Setting thick to ", sugg, "\n")
      thick <- sugg
      info$thick <- thick #CHANGED: if the answer is "yes", the global thick value is not updated
      info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), by=thick)

      if(length(info$slump) > 0) # why here, and not a few lines later?
        info$elbows <- seq(floor(info$d.min), toslump(ceiling(info$d.max), info$slump, remove=remove), by=thick)
      info$K <- length(info$elbows)
      info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model
    }

  ### prepare for any slumps
  if(length(slump) > 0) {
    if(length(slump) %% 2 == 1)
      stop("slumps need both upper and lower depths. Please check the manual", call.=FALSE)
    slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
    info$slump <- slump

    slumpdmax <- toslump(ceiling(info$d.max), slump, remove=remove)
    info$elbows <- seq(floor(info$d.min), slumpdmax, by=thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min+(info$thick*info$K) # the maximum depth to be used by the bacon model

    info$slumpfree <- toslump(depths, slump, remove=remove)
    info$slumphiatus <- toslump(info$hiatus.depths, slump, remove=remove) # check
    if(!is.na(info$boundary[1])) {
      info$slumpboundary <- toslump(info$boundary, slump, remove=remove) # check
      info$slumphiatus <- info$slumpboundary
    }
    slumpdets <- info$dets
    slumpdets[,4] <- toslump(slumpdets[,4], slump, remove=remove) # by default, dates within slumps are not removed
    info$slumpdets <- slumpdets[!is.na(slumpdets[,4]),]
  }

  ### produce files
  info$prefix <- paste0(coredir, core, "/", core, runname, "_", info$K)
  info$coredir <- coredir
  info$bacon.file <- paste0(info$prefix, ".bacon")
  if(!file.exists(outfile <- paste0(info$prefix, ".out")))
    file.create(outfile)
  
  ### if the dates file has been modified after the outfile, suggest to clean up 
  if(file.mtime(outfile) < file.mtime(paste0(info$coredir, core, "/", core, ".csv")))
    message("Warning! The file with the dates seems newer than the run you are loading. If any dates have been added/changed/removed?, then please run Bacon.cleanup()")

  ### store values (again) for future manipulations
  if(BCAD)
    info$BCAD <- TRUE
  if(!is.na(boundary[1])) {
    if(length(slump) > 0)
      boundary <- info$slumpboundary
    info$hiatus.depths <- boundary
    if(length(add) == 0)
      add <- max(1, 1.5*max(info$acc.mean)) # then add a short (max)hiatus, large enough not to crash Bacon but not affect the chronology much. Needs more work
    info$hiatus.max <- add
  }
  if(save.info)
    assign_to_global("info", info)

  prepare <- function() {
    ### plot initial data and priors
    pn <- c(1,2,3,3)
    if(!is.na(info$hiatus.depths[1])) # was ...hiatus.depths)[1])
      if(is.na(info$boundary[1]))
        pn <- c(1,2,3,4,4,4)
    layout(matrix(pn, nrow=2, byrow=TRUE), heights=c(.3,.7))
    oldpar <- par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l", yaxs="i")
    on.exit(par(oldpar))
    PlotAccPrior(info$acc.shape, info$acc.mean, depth.unit=depth.unit, age.unit=age.unit)
    PlotMemPrior(info$mem.strength, info$mem.mean, thick, info)

    if(!is.na(info$hiatus.depths)[1])
      if(is.na(info$boundary)[1])
        PlotHiatusPrior(info$hiatus.max, info$hiatus.depths)
    calib.plot(info, BCAD=BCAD)
    legend("topleft", core, bty="n", cex=1.5)
  }

  cook <- function() {
    bacon.its(ssize, burnin, info) # information on amounts of iterations
    txt <- paste0(info$prefix, ".bacon")
    #cat("this is the bacon file: ", txt)
    bacon(txt, as.character(outfile), ssize+burnin, cc.dir)
    info <- scissors(burnin, info, save.info=save.info)
	output <- info$output # tmp
    info <- agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=TRUE, age.unit=age.unit, depth.unit=depth.unit, save.info=save.info, ...)
	info$output <- output
    #    cat(mean(info$Tr)) # this is to check how hists and info get saved

    if(plot.pdf)
      if(dev.interactive())
        export.pdf(paste0(info$prefix, ".pdf")) else {
          if(capabilities("cairo"))
            cairo_pdf(filename=paste0(info$prefix, ".pdf")) else 
              pdf(file=paste0(info$prefix, ".pdf"))
            agedepth(info, BCAD=BCAD, depths.file=depths.file, depths=depths, verbose=FALSE, age.unit=age.unit, depth.unit=depth.unit, save.info=FALSE, ...)
            dev.off()
        }		
  }

### run bacon if initial graphs seem OK; run automatically, not at all, or only plot the age-depth model
  write.Bacon.file(info, younger.than=younger.than, older.than=older.than, save.info=save.info)
  if(!run)
    prepare() else
      if(!ask)
        cook() else {
          prepare()
          if(accept.suggestions)
            ans <- "y" else
              ans <- readline(message(" Run ", core, " with ", info$K, " sections? (Y/n) "))
          ans <- tolower(substr(ans,1,1))[1]
          if(ans=="y" || ans=="")
            cook() else
              message("  OK. Please adapt settings")
        }
 # if(close.connections)
 #   close(outfile)

  if(save.elbowages) {
    saved <- sapply(info$elbows, Bacon.Age.d)
    write.table(saved, paste0(info$prefix, "_",info$d.min, "_", thick, "_elbowages.txt"), row.names=FALSE, col.names=FALSE)
  }

  invisible(info) # MB Jan 2024
}



#' @name tofu
#' @title Bacon for vegans
#' @details A vegan wrapper for Bacon - does everything Bacon does, but without the meat.
#' @param ... options for the Bacon command. See \link{Bacon}
#' @return A tofu age-model
#' @export
tofu <- function(...)
  Bacon(...)
