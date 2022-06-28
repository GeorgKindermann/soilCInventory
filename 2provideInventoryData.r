pathInvDat <- "./dat/"

invPos <- read.table(header = TRUE, text="
x      y      nn  plotId
413573 274537 456 504100
413773 274537 478 504108")
write.table(invPos, paste0(pathInvDat, "posInv.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

Dd <- read.table(header = TRUE, text="
plotId       date peri share
504100 2005-04-27    1     1
504100 2009-06-12    2     1
504100 2015-04-10    3   0.7
504100 2020-08-05    4   0.7
504108 2005-01-27    1     1
504108 2009-02-12    2     1
504108 2015-03-10    3     1
504108 2020-04-05    4     1")
write.table(Dd, paste0(pathInvDat, "dateInv.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)

Ds <- read.table(header = TRUE, text="
plotId treeId species lastAlivePeri whatOut
504100      1       2             3 n
504100      2      10            NA 0.9st.5ba
504100      3       1            NA 1st1ba
504108      1       1            NA n
504108      2       1            NA stba")
write.table(Ds, paste0(pathInvDat, "speciesInv.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
#lastAlivePeri .. needs to be set for standing dead trees
#whatOut .. what will be harvested and to what share:
# n .. nothing, st .. stem, ba .. bark, br .. branch, le .. leaf, su .. stump,
# cr .. coarseRoot, fr .. fineRoot, se .. seed

Dt <- read.table(header = TRUE, text="
plotId treeId peri nrep grep    d    h  hka
504100      1    1    0   16  13.1 12.7  4.9
504100      1    2    0   16  16.6 16.4  8.1
504100      1    3  600    0  19.9 19.6 10.2
504100      1    4  400    0  23.3 22.3 11.2
504100      2    1 1100    0 12.1 11.7  3.9
504100      2    2  800    0 15.6 17.4  7.1
504100      2    3  650    0 18.9 18.6  9.2
504100      3    4 8000    0  1.0  2.0  0.0
504108      1    1  600    0 32.9 28.1 15.0
504108      2    4 9000    0  1.9  1.5  0.3")
write.table(Dt, paste0(pathInvDat, "treeInv.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
#nrep .. Sampling proportional to stemmnumber
#grep .. Sampling proportional to basal area

