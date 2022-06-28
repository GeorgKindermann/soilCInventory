OUTfileBm <- "./dat/bm2soil.txt"
OUTfileYassoWetter <- "./dat/yassoWetter.txt"
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

Dd <- read.table(paste0(pathInvDat, "dateInv.txt.xz"), header=TRUE)
Dd$date <- as.Date(Dd$date)
Ds <- read.table(paste0(pathInvDat, "speciesInv.txt.xz"), header=TRUE)
Dt <- read.table(paste0(pathInvDat, "treeInv.txt.xz"), header=TRUE)
Dt$alive <- 1
. <- c("plotId","treeId")
Dt$alive[Dt$peri > Ds$lastAlivePeri[match(do.call(paste, Dt[.]), do.call(paste, Ds[.]))]] <- 0

Dd$year <- sapply(seq_len(nrow(Dd)), function(line) {
y <- as.character(Dd$plotId[line])
x <- Dd$date[line]
year <- format(x, "%Y")
dates <- seq(as.Date(paste0(year, "-1-1")), length = 13, by = "1 months")
dates <- as.POSIXct(paste(dates[1:12], "00:00:00")) + (diff(dates) / 2)
. <- Dw[,year,y,"t"]
i <- which.max(.)
gpStart <- as.POSIXct(approxfun(.[1:i], dates[1:i])(5), origin=dates[1] - as.numeric(dates[1]))
gpEnd <- as.POSIXct(approxfun(.[i:12], dates[i:12])(5), origin=dates[1] - as.numeric(dates[1]))
gpStart <- pmin(gpStart, as.POSIXct(paste0(year, "-6-1")), na.rm=TRUE)
gpEnd <- pmax(gpEnd, as.POSIXct(paste0(year, "-9-1")), na.rm=TRUE)
as.numeric(year) + pmin(1 ,pmax(0, as.numeric(difftime(x, gpStart, units="days")) / as.numeric(difftime(gpEnd, gpStart, units="days"))))
})

Dd$yearC <- sapply(Dd$date, \(x) {
  as.numeric(format(x,"%Y")) + #Get Year an add
    (as.numeric(format(x,"%j")) - 0.5) / #Day of the year divided by
    as.numeric(format(as.Date(paste0(format(x,"%Y"), "-12-31")),"%j")) #days of the year
})
#Weight with days in month?
. <- list(t = round(rowMeans(aperm(Dw[,,,"t"], 3:1), dims=2), 1)
         , p = round(rowSums(aperm(Dw[,,,"p"], 3:1), dims=2))
         , amp = round(apply(Dw[,,,"t"], 3:2, \(x) diff(range(x))/2), 1) )
. <- do.call(cbind
  , lapply(names(.), \(x) `colnames<-`(.[[x]], paste0(x, colnames(.[[x]])))) )
write.table(., OUTfileYassoWetter, quote = FALSE, row.names = TRUE, col.names = TRUE)
rm(Dw)

obsE <- list2env(split(Dd[c("peri", "year", "yearC")], Dd$plotId))
treeE <- list2env(split(Dt, Dt$plotId))
rm(Dd, Dt)

. <- do.call(rbind, lapply(ls(obsE), function(plotId) {
  obs <- get(plotId, obsE)
  obs <- obs[order(obs$year),]
  trees <- get0(plotId, treeE)
  yearStart <- ceiling(max(min(obs$year), min(obs$yearC)))
  yearEnd <- floor(min(max(obs$year), max(obs$yearC)))
  if(yearEnd < yearStart) yearEnd <- yearStart <- round(mean(obs$year))
  years <- yearStart:yearEnd
  if(nrow(obs) > 1 & !is.null(trees)) {
    . <- do.call(rbind, lapply(split(trees, trees$treeId), function(tree) {
      . <- merge(obs, tree[c("peri", "nrep", "grep", "d", "h", "hka", "alive")], all.x=TRUE)
      .$REFn <- .$nrep
      .$REFg <- .$grep
      .$nrep[is.na(.$nrep)] <- 0
      .$grep[is.na(.$grep)] <- 0
      for(COL in c("d", "h", "hka", "alive", "REFn", "REFg")) {
        i <- diff(is.na(.[[COL]]))
        j <- which(i == -1)
        .[[COL]][j] <- .[[COL]][j+1]
        j <- which(i == 1)
        .[[COL]][j+1] <- .[[COL]][j]
      }
      data.frame(plotId
               , treeId = tree$treeId[1]
               , year = years
               , nrep = approxfun(.$yearC, .$nrep)(years)
               , grep = approxfun(.$yearC, .$grep)(years)
               , d = approxfun(.$year, .$d, rule=2)(years)
               , h = approxfun(.$year, .$h, rule=2)(years)
               , hka = approxfun(.$year, .$hka, rule=2)(years)
               , alive = approxfun(.$yearC, .$alive, rule=2)(years)
               , nout = pmax(0, diff(approxfun(.$yearC, .$REFn, rule=2)(c(years[1], years+1)) - approxfun(.$yearC, .$nrep, rule=2)(c(years[1], years+1))))
               , gout = pmax(0, diff(approxfun(.$yearC, .$REFg, rule=2)(c(years[1], years+1)) - approxfun(.$yearC, .$grep, rule=2)(c(years[1], years+1))))
                )
    }))
  } else {
    if(is.null(trees)) {assign(".", NULL)
    } else {
      . <- data.frame(plotId
                    , treeId = trees$treeId
                    , year = years
                    , trees[c("nrep", "grep", "d", "h", "hka", "alive")]
                    , nout = 0, gout = 0)
      }
  }
  . <- .[.$nrep > 0 | .$grep > 0 | .$nout > 0 | .$gout > 0,]
  i <- years[!years %in% unique(.$year)]
  if(length(i) > 0) {
    . <- rbind(., data.frame(plotId, treeId = NA, year = i, nrep = 0, grep = 0
                           , d = 0, h = 0, hka = 0, alive=1, nout = 0, gout=0))}
  .
}))
rm(obsE, treeE)

. <- merge(Ds[c("plotId", "treeId", "species", "whatOut")], ., all.y=TRUE)
rm(Ds)
.$treeId[is.na(.$treeId)] <- 0
.$species[is.na(.$species)] <- 1
.$whatOut[is.na(.$whatOut)] <- "n"

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
nOut <- Db$nout + ifelse(Db$Grepjeha > 0, nrepG * Db$gout / Db$Grepjeha, 0)
#Input from living trees
res <- Db[,j] * Ctr[i, k] *  (Db$Nrepjeha + nrepG - nOut) * Db$alive

nt <- which(c(TRUE, Db$plotId[-nrow(Db)] != Db$plotId[-1] | Db$treeId[-nrow(Db)] != Db$treeId[-1]))
mort <- c(0, diff((Db$Nrepjeha + nrepG) * (1 - Db$alive))) + nOut * (1 - Db$alive)
mort[nt] <- 0
#Input from trees which have died and stay in forest
res <- res + sweep(Db[,j], 2, c(1, 1, 0, 0, 1, 0, .5, 0), `*`) * mort

#Input from dead trees in forest
res <- res + apply(sweep(rbind(Db[1,j], Db[,j]), 2, c(0, 0, 1, 1, 0, 1, 0, 1),
                         `*`), 2, function(x) {-pmin(0, diff(x))}) *
  (Db$Nrepjeha + nrepG - nOut) * (1 - Db$alive)

#Input from removed trees
res <- res + Db[,j] * nOut *
  sapply(c("le", "br", "st", "cr", "fr", "su", "se", "ba"), \(x) {
    m <- numeric(length(Db$whatOut))
    . <- regexpr(paste0("[0-9.]*(?=", x, ")"), Db$whatOut, perl=TRUE)
    m[.>0] <- 1
    m[attr(., "match.length") > 0] <- as.numeric(regmatches(Db$whatOut, .)[attr(., "match.length")[.>0] > 0])
    if(x=="se") pmax(0, 0.5 - m) else 1 - m
  })

write.table(cbind(Db[c("plotId", "treeId", "species", "year", "d")], res),
            OUTfileBm, quote = FALSE, row.names = FALSE, col.names = TRUE)
