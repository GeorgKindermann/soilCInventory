OUTfileYassoCC <- "./dat/yassoCCC.txt"
fileBm2Soil <- "./dat/bm2soilC.txt"

Cb2c <- read.csv("./bmExpansion/coefBm2C.csv")
Dbms <- read.table(fileBm2Soil, header=TRUE)

FRAK <- read.table(header=TRUE, text="
  bm           b2c        cc
  bmLeaf       leaf       leaf
  bmBranch     branch     branch
  bmStem       stem       stem
  bmCoarseRoot coarseRoot stem
  bmFineRoot   fineRoot   fineRoot
  bmStump      stem       stem
  bmSeed       seed       leaf
  bmBark       bark       branch
")

. <- Dbms[FRAK$bm] / 1000 *  #from kg to t
     Cb2c[match(Dbms$species, Cb2c$species), FRAK$b2c] #from bm dry matter to C

Ccc <- read.csv("./bmExpansion/coefChemicalComposition.csv", comment.char = "#")
tt <- unique(Ccc$species)
Ccc <- simplify2array(lapply(split(Ccc[c("species", "a", "w", "e", "n")],
                            Ccc$LitterType)
     , \(x) as.matrix(`row.names<-`(x[match(tt, x$species), -1], x$species))))
rm(tt)

. <- Map(`*`, .[FRAK$bm], asplit(Ccc[Dbms$species, , FRAK$cc], 3))

. <- rbind(
  cbind(d=0, Reduce('+', .[c("bmLeaf", "bmFineRoot", "bmSeed", "bmBark")]))
, cbind(d=round(Dbms$d * 0.3), Reduce('+', .[c("bmBranch", "bmCoarseRoot")]))
, cbind(d=round(Dbms$d * 0.8), .[["bmStem"]])
, cbind(d=round(Dbms$d), .[["bmStump"]]))

. <- aggregate(.~ plotId + year + d, cbind(Dbms[c("plotId", "year")], .), sum)

. <- do.call(rbind, lapply(split(., paste(.$plotId, .$year)), function(x) {
  i <- which(rowSums(x[c("a", "w", "e", "n")]) > 0)
  if(length(i) > 0) x[i,]
  else transform(x[1,], d=0)
}))

#. <- .[order(.$plotId, .$d, .$year),]

write.table(., OUTfileYassoCC, quote = FALSE, row.names = FALSE, col.names = TRUE)
