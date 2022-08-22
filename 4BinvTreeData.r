OUTfileBm <- "./dat/bm2soilB.txt"
#pathInvDat <- "./dat/"
pathInvDat <- "./2invData/"
fileTemperature <- "./dat/tmInvPkt.txt"
filePrecipitation <- "./dat/rrInvPkt.txt"

library(Rcpp)
sourceCpp("./bmExpansion/bmExpansion.cc")

Dwt <- read.table(fileTemperature)
Dwp <- read.table(filePrecipitation)
.Id <- sort(unique(Dwt[,3]))
.Y <- sort(unique(Dwt[,1] %/% 100))
Dw <- array(c(Dwt[order(Dwt[,3], Dwt[,1]), 2], Dwp[order(Dwp[,3], Dwp[,1]), 2])
           , c(12, length(.Y), length(.Id), 2)
         , list(month=1:12, year=.Y, plotId=.Id, what=c("t", "p")))
rm(Dwt, Dwp, .Id, .Y, fileTemperature, filePrecipitation)

Ds <- read.table(paste0(pathInvDat, "speciesInv.txt.xz"), header=TRUE)
Dt <- read.table(paste0(pathInvDat, "treeInv.txt.xz"), header=TRUE)
Dt$alive <- 1
. <- c("plotId","treeId")
Dt$alive[Dt$peri > Ds$lastAlivePeri[match(do.call(paste, Dt[.]), do.call(paste, Ds[.]))]] <- 0
Dt <- Dt[Dt$alive == 1,]
Dt <- Dt[Dt$peri == 9,]
Dt$peri <- NULL

. <- merge(Ds[c("plotId", "treeId", "species")], Dt, all.y=TRUE)
.$whatOut <- "n"
.$year <- 0
.$nout <- 0
.$gout <- 0
. <- .[c("plotId", "treeId", "species", "whatOut", "year", "nrep", "grep", "d", "h", "hka", "alive", "nout", "gout")]

Fmeasurement <- tempfile()
Fbiomass <- tempfile()
write.table(., Fmeasurement, , quote = FALSE, na = "0", row.names = FALSE,
            col.names = FALSE)
rm(.)
calcBm(Fmeasurement, Fbiomass)

Db <- read.table(Fbiomass, header=TRUE)
Db <- Db[do.call(order, Db[c("plotId", "treeId", "year")]),]

Ctr <- read.csv("./bmExpansion/coefTrunoverRates.csv")

FRAK <- c("leaf", "branch", "stem", "coarseRoot", "fineRoot", "stump", "seed", "bark")
i <- match(Db$species, Ctr$species)
j <- match(tolower(FRAK), tolower(substring(names(Db), 3)))
k <- match(tolower(sub("stump", "stem", FRAK, fixed=TRUE)), tolower(names(Ctr)))
nrepG <- ifelse(Db$d > 0, Db$Grepjeha / (Db$d/200)^2 / pi, 0)
#Input from living trees
res <- Db[,j] * Ctr[i, k] *  (Db$Nrepjeha + nrepG) * Db$alive

#Input from dead trees in forest
#res <- res + sweep(Db[,j], 2, c(0, 0, 1, 1, 0, 1, 0, 1), `*`) *
#  (Db$Nrepjeha + nrepG) * (1 - Db$alive) / 35 #years as soil inflow

write.table(cbind(Db[c("plotId", "treeId", "species", "year", "d")], res),
            OUTfileBm, quote = FALSE, row.names = FALSE, col.names = TRUE)
