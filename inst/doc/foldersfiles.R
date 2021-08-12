## ---- eval=FALSE--------------------------------------------------------------
#  mydir <- tempdir()
#  Bacon(coredir=mydir)

## ---- eval=FALSE--------------------------------------------------------------
#  Bacon(suggest=FALSE)
#  Bacon(accept.suggestions=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  allcores <- list.files("Bacon_runs")
#  for(i in allcores)
#    try(Bacon(i, accept.suggestions=TRUE))

