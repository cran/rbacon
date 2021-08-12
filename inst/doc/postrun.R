## ---- eval=FALSE--------------------------------------------------------------
#  Baconvergence("MSB2K", thick=5, runs=5, ssize=100, coredir=tempfile())

## ---- eval=FALSE--------------------------------------------------------------
#  proxy.ghost(7)

## ---- eval=FALSE--------------------------------------------------------------
#  accrate.depth.ghost()
#  accrate.age.ghost()

## ---- echo=FALSE, include=FALSE-----------------------------------------------
require(rbacon)
Bacon("MSB2K", ask=FALSE, coredir=tempdir(), suggest=FALSE)
agedepth()

## -----------------------------------------------------------------------------
Bacon.hist(20)

## -----------------------------------------------------------------------------
a.d20 <- Bacon.Age.d(20)
summary(a.d20)
hist(a.d20)

## -----------------------------------------------------------------------------
a.d30 <- Bacon.Age.d(30)
a.d20 <- Bacon.Age.d(20)
summary(a.d30-a.d20)
hist(a.d30-a.d20)

## -----------------------------------------------------------------------------
acc.d20 <- accrate.depth(20)
summary(acc.d20)
acc.a4500 <- accrate.age(4500)
summary(acc.a4500)

