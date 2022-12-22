# rbacon 3.0.0
* accrates.depth.ghost() and accrate.age.ghost() now invisibly return the ranges, medians and means for each depth resp. age for subsequent use, e.g., tmp <- accrates.depth.ghost(); head(tmp)
* MCMC iterations are now stored in the .out files irrespective of whether they were accepted or rejected
* ssize is now much more predictable (if there are more rows in the output file than set by ssize, rbacon keeps only the last n=ssize rows)
* files are read and written faster (assuming that the data.table R package is installed and loaded, which can be a problem on Macs)
* corrected the help description for Bacon.d.Age() with thanks to Henningte
* the length of the ages output of accrate.age() is now the same even if there are NAs in the output (with thanks to Henningte)
* corrected (hopefully) a bug in read.dets related to logical comparisons with variable lengths causing errors in R>=4.2 (with thanks to Nick McKay for reporting)
* accrate.age.ghost gains a kcal option
* repaired warning message about lengths of logical tests
* now links to rintcal package (renamed from IntCal)
* better Plum plots
* runs can now be interrupted by pressing ctrl+c (experimental)

# rbacon 2.5.8
* some minor updates to the vignettes
* corrected a bug in input.cpp which caused a gcc warning in debian

# rbacon 2.5.7
* added an option to agedepth to plot date labels (plot.labels)
* added vignettes
* if ask=FALSE, Bacon now does not ask before writing a new folder (if so required)

# rbacon 2.5.6
* removed closeAllConnections() as requested by Kurt Hornik (CRAN)
* adapted agedepth() for further Plum corrections
* corrected the behaviour of the dark option in agedepth(), accrate.depth.ghost(), accrate.age.ghost(), flux.age.ghost() and proxy.ghost()
* renamed rplum's option radon.case to ra.case
* added plot.median to flux.age.ghost() and accrate.age.ghost()
* repaired BCAD in flux.age.ghost()
* removed the error message that acc.shape cannot be equal to acc.mean (shouldn't be a problem any more)

# rbacon 2.5.5
* removed revdep folder which caused issues when submitting to CRAN

# rbacon 2.5.4
* further separation of rplum and rbacon
* the functions thinner and scissors now deal with Plum runs
* repaired bugs in accrate.age and accrate.age.ghost related to the option BCAD

# rbacon 2.5.3
* corrected a bug where postbomb dates could not be plotted owing to wrong by sign
* further improvements to how agedepth deals with plot margins
* corrected error when d.min was set

# rbacon 2.5.2
* optimised accrate.age.ghost()
* added options to modify the margins of the individual panels in the agedepth plot: mar.left, mar.middle, mar.right, mar.main
* added an option to plot the tickmarks and labels on the vertical axes of the prior panels: prior.ticks="s" (default "n")
* removed the panels option in agedepth() as it didn't work as expected and is better done outside rbacon functions (e.g., layout(1); agedepth(model.only=T))
* added median curves for accrate.age.ghost and accrate.depth.ghost (as means can be influenced by extreme values)

# rbacon 2.5.1
* Adapted the default prior for memory to 0.5 (mean) and 10 (strength), to repair a bug with the original bacon c++ code. This default should work with most cores and give similar results to the previous settings for the memory prior
* Updated c/c++ code (as used in version 2.4.1 with some minor additional updates) 
* pMC.age and other IntCal functions should now load as expected (without having to specify, e.g., IntCal::pMC.age)
* the heights of the calibrated distributions should now scale better according to how precise they are
* the add.dates function now handles postbomb dates. It can also store the calibrated information using, e.g., tmp <- add.dates(2450,30,20); tmp
* The greyscale age-depth graph is now more easily exported to external graphics editors, because areas with very low probabilities are now left empty (instead of plotted as white)
* corrected a bug where ages above d.min received incorrect ages
* Added an option prior.ticks to show tickmarks and values on the vertical axes of the panels that show the prior distributions. These are not drawn by default, as they don't provide much information and clutter the graphs
* Added new options title.size and prior.fontsize for the size of the fonts of the core's title and the red information on settings in the top panels, respectively
* Repaired the functions accrate.depth(), accrate.age(), accrate.depth.ghost() and accrate.age.ghost()
* agemodel.it now treats the upper depth of a core as expected

# rbacon 2.5.0
* updated src/kernel.cpp and src/twalk.h, to repair a bug in one of the moves ('hop'). This means we can now add the updated MCMC code of version 2.4.0 again and accommodate code to run 210Pb-dated cores (via the package rplum) 
* Radiocarbon calibration curves are now loaded from the imported IntCal R package, and have been removed from the rbacon package to save space and remove duplication
* Added option rgb.scale to draw shades of other colours than grey, e.g. red: rgb.scales(1,0,0), for the functions agedepth, accrate.depth.ghost, accrate.age.ghost, proxy.ghost and flux.age.ghost (based on an idea kindly provided by Oliver Wilson).
* Related to rgb.scale, the resolution of the colours has been renamed from grey.res to rgb.res
* added option 'add' to add proxy.ghost graphs to existing plots (based on an idea by Oliver Wilson)
* Repaired bug with greyscales accrate.age.ghost
* if the file with the dates has been modified more recently than a loaded run (e.g., dates could have been added, removed or changed), then a warning is now given that Bacon.cleanup should be ran
* Added a new function Bacon.d.Age to provide the depths belonging to a specific modelled age (kindly contributed by Timon Netzel)
* Depth units are now handled better by the agedepth function
* New option accept.suggestions, which automatically accepts suggestions regarding acc.rate and thick. Use with caution (this option was kindly suggested by Quinn Asena)
* By default, a section is now added below the bottom-most dated depth, in order to ensure that the this depth is always taken into account. Defaults to add.bottom=TRUE. 
* The calibrated distributions should now be of the same size again, so that more precise dates peaks more than less precise ones (suggested by Tiffany Napier). 

# rbacon 2.4.3
* replaced 'cat' with 'message' or 'warning' where possible
* updated to the IntCal20 calibration curves (Reimer et al., 2020)

# rbacon 2.4.2
* Reverted the c/c++ code back to that of version 2.3.9.1, owing to problems with the code introduced in version 2.4.0 (posteriors of memory and age-model are apparently too wide)
* Repaired bug that caused an error when using slumps or boundaries
* Enhanced behaviour of rotate.axes option

# rbacon 2.4.1
* Updated the code to deal with changes in how base-R deals with c() in loops, as suggested by Martin Maechler's e-mail 29 February 2020

# rbacon 2.4.0
* The MCMC code has been updated to remove bugs and to accommodate runs with the upcoming 'rplum' package for 210Pb dating
* Added functions which are required to run the 'rplum' package (although 'rbacon' does not require 'rplum' to be installed)

# rbacon 2.3.9.1
* Added a new option calheight, which acts as a multiplier for the relative height of non-14C dates
* Set default for y-axis to have no space added after the extreme values (yaxs="i"); x-axis has some space added by default (xaxs="r")
* Added a bit of space to d.max and d.min in the main age-depth graph, to accommodate age blobs
* New option kcal, which gives tick marks every 1,000 cal years (default kcal=FALSE)
* Corrected an error when running a core with 4 columns in the .csv file and cc=0
* depth.unit and age.unit now work correctly when provided as options in Bacon or agedepth
* Corrected a bug where thickness (dC) was sometimes internally set to wrong values
* Redid hiatuses: If a core has one or more hiatuses, then variables slopes.above and slopes.below are made for each hiatus, and used internally to adapt ages and accumulation rates for each depth below and above a hiatus within a section containing a hiatus. 
* Slumps, hiatuses and boundaries have gone through a thorough check and should now work better than they did before. Reports of weird behaviour welcome!
* Renamed info\$d to info\$elbows (internal; for better consistency with the naming of parameters within the Bacon paper)

# rbacon 2.3.8
* repaired a bug in cal.h which prevented the postbomb curve postbomb_SH3 from being used
* repaired bug where the prior for the accumulation rate would not always be drawn entirely
* Bacon.hist now takes alternative values for prob into account (e.g., prob=.68)
* The agedepth function now deals better with d.min and d.max values
* Colours of cal BP dates now as expected when cc=0 is provided as Bacon option 
* The fit of the dates to the age-model is now reported correctly also when BCAD=TRUE
* Date distributions should now plot as expected over a wider range of values
* New option acc.lab to provide alternative label for the accumulation rate axis (top-middle panel of the main agedepth graph)
* When provided, d.max or d.min are now dealt with better if extra columns are provided for dR/dSTD and/or t.a/t.b in the core's .csv file
* New options depth.unit (default 'cm') and age.unit (default 'yr'), deprecating the previous poorly named option 'unit' which defaulted to 'cm'. So can now also deal with, e.g., 'Ma' and 'km'
* Replaced occurrences of yr with the more generic unit of age (deprecate yr.min, yr.max, MinYr, MaxYr)
* When Bacon asks for confirmation to run a core (Y/n), the user can now simply press Enter instead of having to type y first. Similarly, by default suggestions to adapt the prior accumulation rate are not accepted (y/N)
* Enhanced drawing of very precise ages (e.g., 1 yr)
* Bacon now stops if there are less than 2 sections between neighbouring hiatuses
* A warning is now given if acc.shape <1 (since this results in weirdly shaped gamma prior distributions)
* An error is thrown when the core's .csv file has 'orphan' commas (can happen if the file was made in a spreadsheet program - check in a plain-text editor)
* add.dates now plots better when mirror=FALSE
* More consistent error messages

# rbacon 2.3.7
* adapted cpp code to allow for more than 10 hiatuses/boundaries (now limited to 50)
* corrected bug causing a warning when a hiatus was set with multiple acc.mean priors provided
* now ensures that hiatus or boundary depths are in the correct order (ascending in depth)

# rbacon 2.3.6
* further enhancements to memory usage in MCMC calculations (bacon.h)

# rbacon 2.3.5
* added function agemodel.it to extract single iterations of a Bacon age-depth model
* added functions clam2bacon and bacon2clam to translate Bacon dates files into clam files et vice versa (inspired by a suggestion from Dewey Dunnington)
* corrected behaviour of boundary and hiatus (especially if together with slumps)
* iterations with age reversals across a hiatus are now removed
* removed closeAllConnections (suggested by Dewey Dunnington)
* Added option to change the field separator to mix.curves (suggested by Thomas Dye)
* MinYr now defaults to the current year (1950 - as.integer(format(Sys.time(), "\%Y")))
* added option in the scissors function to remove a specific range of iterations (e.g., iterations 400 to 800, or the first/last 300)
* produced separate R files for groups of functions
* Bacon now stops if it finds 6 columns with unexpected names in the .csv file. If provided with a delta.R column, Bacon expects a delta.STD column as well. 
* Added an option dates.col to colour sets of dates (suggestion by Greg Cooper)
* enhancements in bacon.h of MCMC calculations and memory usage 

# rbacon 2.3.4
* faster drawing of greyscale plots (though still slower yet better than in version 2.3.1.1 and before)
* added progress bar to functions that can be slow
* repaired a bug in calculating how many dates fall within the model range
* delta.R is now accepted as a header for the dates file

# rbacon 2.3.3
* added an option to include slumps (sort of - more testing still welcome). Example: Bacon(slump=c(50, 52, 60, 70)) for two slumps between 50-52 and 60-70 cm depth
* date-files with .csv.txt extensions are now renamed to .csv (and informing us that it did so)
* default darkness of age-depth greyscale now adapts to a ratio between most and least precise sections (so that very imprecise sections still show some grey)
* repaired option depths (e.g., Bacon(depths=0:100))
* repaired height of prior distribution axes
* repaired Baconvergence()
* added a commentary after each run, mentioning the proportion of dates that lie within the age-depth model's range (some sort of 'agreement')

# rbacon 2.3.2
* Added option boundary, which sets hiatus length to (close to) 0. This leaves the hiatus functionality more or less unchanged, and should cause less confusion with setting hiatus.depths even if no hiatus is desired.
* Enhanced plotting and age calculation of depths close to hiatuses or boundaries.
* Ensured more predictable behaviour if R is started in a non-writable directory (e.g. plain, non-Rstudio R on Windows). 
* Added confidence ranges to accrate.age.ghost and accrate.depth.ghost.
* Enhanced calculation of mean and median (now based on age distribution, not on a derived histogram).
* Corrected behaviour of title.location.
* Corrected many sundry bugs related to plotting, especially with hiatuses or with BCAD=TRUE.
* Added a `NEWS.md` file to track changes to the package.

# rbacon 2.3.1.1
* Now a CRAN R package (not called bacon since that name was already taken).
* Default core directory now Bacon_runs. Other directories can be given, for more flexibility in workflows of users. 
* Calibration curves can be put in a user-specified directory ccdir (hidden by default).
* New function copyCalibrationCurve() to copy calibration curves into an R's session.
* Renamed several options to be more consistent, d.R and d.STD now named delta.R and delta.STD.
* Can now provide depths to be calculated as a variable, as alternative to using a file with depths.
* Added option to not plot x or y axis (xaxt, yaxt).
* Added option to not plot the date distributions mirrored.
* New function Baconvergence() to test for good mixing of MCMC runs. 
* Renamed weighted means of age estimates to means.
* Updated documentation.
* Renamed functions flux.age, plot.accrate.age and plot.accrate.depth to flux.age.ghost, accrate.age.ghost and accrate.depth.ghost, respectively. 
* BCAD dealt with more correctly.
* Repaired many sundry bugs.

# Bacon 2.2
* Updated to 14C calibration curves IntCal13, Marine13 and SHCal13.
* Changed .hpd to _ages.txt since many users get tricked by the extension.
* Changed from .dat files to .csv files as these are more documented and easier to open and edit by users.
* Separator for .csv file can be adapted.
* Renamed res to hopefully more intuitive thick (thickness of sections)
* Added d.R and d.STD
* Bacon.hist gives 95% ranges, mid and wmean, and reads from a file instead of from the command line. 
* Added options to change axis limits, orientation and rotation. 
* BCAD introduced, though not yet working entirely as expected.
* Different prior for acc.mean suggested if initial estimates indicate that this would be beneficial. 
* Introduced a settings file.
* removed calc.every (gave problems with long cores).
* Killed hist bug that assumed integers.
* Language cleanup of cpp files.
* Added option to remove unnecessary files after a run.
* Added option in agedepth to only plot the age-model (so not the upper panels).
* Many bug fixes in the Bacon.R and underlying C/C++ codes
