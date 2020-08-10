
################################################################
### Program to translate the IntCal20 curves to the standard format:
### cal. BP   rc BP    std. error,
### in cal. BP increasing order, with no commas.
#################################################################


### Read the terrestrial and Marine curves files as distributed in
### http://intcal.org/curves:

cc.terr <- read.csv( "intcal20.14c", header=FALSE, skip=11, sep=",")
cc.marine <- read.csv( "marine20.14c", header=FALSE, skip=11)
cc.south <- read.csv( "shcal20.14c", header=FALSE, skip=11, sep=",")

### Write them in the desired format:

write.table( cc.terr[nrow(cc.terr):1,1:3], file="3Col_intcal20.14C", row.names=FALSE, col.names=FALSE)
write.table( cc.marine[nrow(cc.marine):1,1:3], file="3Col_marine20.14C", row.names=FALSE, col.names=FALSE)
write.table( cc.south[nrow(cc.south):1,1:3], file="3Col_shcal20.14C", row.names=FALSE, col.names=FALSE)
