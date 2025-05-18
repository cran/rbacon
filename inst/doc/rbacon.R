## ----eval=FALSE---------------------------------------------------------------
# install.packages('rbacon')

## -----------------------------------------------------------------------------
library(rbacon)

## ----eval=FALSE---------------------------------------------------------------
# Bacon()

## ----eval=FALSE---------------------------------------------------------------
# ?Bacon

## ----eval=FALSE---------------------------------------------------------------
# ?agedepth

## ----eval=FALSE---------------------------------------------------------------
# Bacon('RLGH3')

## ----eval=FALSE---------------------------------------------------------------
# Bacon('RLGH3')

## ----eval=FALSE---------------------------------------------------------------
# Bacon('RLGH3', acc.mean=50, acc.shape=100)

## ----eval=FALSE---------------------------------------------------------------
# Bacon("RLGH3", acc.mean=50, hiatus.depths=125, hiatus.max=1000)

## ----eval=FALSE---------------------------------------------------------------
# Bacon("RLGH3", acc.mean=50, slump=c(180, 120, 40, 30))

## ----eval=FALSE---------------------------------------------------------------
# mydir <- tempdir()
# Bacon(coredir=mydir)

## ----eval=FALSE---------------------------------------------------------------
# Bacon(suggest=FALSE)
# Bacon(accept.suggestions=TRUE)

## ----eval=FALSE---------------------------------------------------------------
# allcores <- list.files("Bacon_runs")
# for(i in allcores)
#   try(Bacon(i, accept.suggestions=TRUE))

## ----eval=FALSE---------------------------------------------------------------
# Baconvergence("MSB2K", thick=5, runs=5, ssize=100, coredir=tempfile())

## ----eval=FALSE---------------------------------------------------------------
# proxy.ghost(7)

## ----eval=FALSE---------------------------------------------------------------
# accrate.depth.ghost()
# accrate.age.ghost()

## ----echo=FALSE, include=FALSE------------------------------------------------
require(rbacon)
Bacon("MSB2K", ask=FALSE, coredir=tempdir(), suggest=FALSE)
agedepth()

## ----fig.width=5, fig.height=4------------------------------------------------
Bacon.hist(20)

## ----fig.width=5, fig.height=4------------------------------------------------
a.d20 <- Bacon.Age.d(20)
summary(a.d20)
hist(a.d20)

## ----fig.width=5, fig.height=4------------------------------------------------
a.d30 <- Bacon.Age.d(30)
a.d20 <- Bacon.Age.d(20)
summary(a.d30-a.d20)
hist(a.d30-a.d20)

## -----------------------------------------------------------------------------
acc.d20 <- accrate.depth(20)
summary(acc.d20)
acc.a4500 <- accrate.age(4500)
summary(acc.a4500)

## -----------------------------------------------------------------------------
accrate.depth.summary(20)
accrate.age.summary(4500)

## -----------------------------------------------------------------------------
accrates.core()

## ----eval=FALSE---------------------------------------------------------------
# Bacon.hist(23.45)
# depth23.45 <- Bacon.Age.d(23.45)
# hist(depth23.45)
# mean(depth23.45)

## -----------------------------------------------------------------------------
ageranges(10:20)

## ----eval=FALSE---------------------------------------------------------------
# Bacon("MyCore", run=FALSE)
# agedepth()
# # or if you've set the accumulation rate prior to something different than the default:
# Bacon("MyCore", run=FALSE, acc.mean=c(20,5))
# agedepth()

## ----eval=FALSE---------------------------------------------------------------
# Bacon("MyCore", thick=1)
# Bacon("MyCore", 1)

## ----eval=FALSE---------------------------------------------------------------
# Bacon("MyCore", boundary=80, acc.mean=c(10,50))

## ----eval=FALSE---------------------------------------------------------------
# Bacon.hist(23.45)
# depth23.45 <- Bacon.Age.d(23.45)
# hist(depth23.45)
# mean(depth23.45)
# quantile(depth23.45, probs=c(.025, .975))

## ----eval=FALSE---------------------------------------------------------------
# depth23.45 <- Bacon.Age.d(23.45)
# depth12.34 <- Bacon.Age.d(12.34)
# passed <- depth23.45 - depth12.34
# hist(passed)
# mean(passed)
# quantile(passed, probs=c(.025, .975))

## ----eval=FALSE---------------------------------------------------------------
# my_fancy_colours <- sample(colours(), nrow(info$dets))
# agedepth(dates.col=my_fancy_colours)

