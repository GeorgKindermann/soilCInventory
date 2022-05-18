pathInvDat <- "./dat/"
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

Dd <- read.table(paste0(pathInvDat, "dateInv.txt"), header=TRUE)
Dd$date <- as.Date(Dd$date)
Ds <- read.table(paste0(pathInvDat, "speciesInv.txt"), header=TRUE)
Dt <- read.table(paste0(pathInvDat, "treeInv.txt"), header=TRUE)
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
### HIER NOCh DIE WERTE FUER YASSO BESTIMMEN ###
rm(Dw)

obsE <- list2env(split(Dd[c("peri", "year")], Dd$plotId))
treeE <- list2env(split(Dt, Dt$plotId))
rm(Dd, Dt)

. <- do.call(rbind, lapply(ls(obsE), function(plotId) {
  obs <- get(plotId, obsE)
  obs <- obs[order(obs$year),]
  trees <- get0(plotId, treeE)
  yearStart <- ceiling(min(obs$year))
  yearEnd <- floor(max(obs$year))
  if(yearEnd < yearStart) yearEnd <- yearStart <- round(mean(obs$year))
  years <- yearStart:yearEnd
  if(nrow(obs) > 1 & !is.null(trees)) {
    . <- do.call(rbind, lapply(split(trees, trees$treeId), function(tree) {
      . <- merge(obs, tree[c("peri", "nrep", "d", "h", "hka", "alive")], all.x=TRUE)
      OUT <- .$nrep[tail(which(diff(is.na(.$nrep)) == 1), 1)]
      .$nrep[is.na(.$nrep)] <- 0
      for(COL in c("d", "h", "hka", "alive")) {
        i <- diff(is.na(.[[COL]]))
        j <- which(i == -1)
        .[[COL]][j] <- .[[COL]][j+1]
        j <- which(i == 1)
        .[[COL]][j+1] <- .[[COL]][j]
      }
      .f <- approxfun(.$year, .$d, rule=2)
      transform( data.frame(plotId
                          , treeId = tree$treeId[1]
                          , year = years
                          , nrep = approxfun(.$year, .f(.$year)^2 * .$nrep)(years) / .f(years)^2
                          , d = .f(years)
                          , h = approxfun(.$year, .$h, rule=2)(years)
                          , hka = approxfun(.$year, .$hka, rule=2)(years)
                          , alive = approxfun(.$year, .$alive, rule=2)(years))
              , nout = if(length(OUT) > 0) rev(diff(c(0, pmin(OUT, cumsum(diff(c(0, rev(nrep)))))))) else 0
                )
    }))
  } else {
    if(is.null(trees)) {assign(".", NULL)
    } else {
      . <- data.frame(plotId
                    , treeId = trees$treeId
                    , year = years
                    , trees[c("nrep", "d", "h", "hka", "alive")]
                    , nout = 0)
      }
  }
  . <- .[.$nrep > 0,]
  i <- years[!years %in% unique(.$year)]
  if(length(i) > 0) {
    . <- rbind(., data.frame(plotId, treeId = NA, year = i, nrep = 0
                           , d = 0, h = 0, hka = 0, alive=1, nout = 0))}
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

Ctr <- read.csv("./bmExpansion/coefTrunoverRates.csv")
Cb2c <- read.csv("./bmExpansion/coefBm2C.csv")
Ccc <- read.csv("./bmExpansion/coefChemicalComposition.csv")
