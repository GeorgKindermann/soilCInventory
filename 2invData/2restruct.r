#COL <- c("peri", "RW", "HW", "PB", "T", "BAUMNR", "BNR_NEU", "AZI", "DIST", "AZI_ALT", "DIST_ALT", "BA", "BART", "BHD", "HOEHE", "KRONHO", "N", "D")
#I have used in a first run BAUMNR, BNR_NEU, AZI_ALT and DIST_ALT to test if they have relevant info
COL <- c("peri", "RW", "HW", "PB", "T", "AZI", "DIST", "BA", "BART", "BHD", "HOEHE", "KRONHO", "N", "D")

#N 'Nutzungsart: 0..Natürlicher Abgang, 1..Kahlhieb, Rodung und Trassenaufhieb, 2..Standraumerweiterung und Läuterung, 3..positive Auslese, 4..negative Auslese (Niederdurchforstung), 5..Verjüngungshieb, 6..Räumung, 7..Entrümpelung, 8..Einzelstammentnahme und Kleinflächennutzung (bis 500 m2), 9..Zufallsnutzung' ?B, F, K, N
#D 'Dürrling: 0..kein Dürrling, 1..unbedeutender Dürrling, 2..bedeutender Dürrling'
#BRUCH is currently not used

. <- do.call(rbind, lapply(list.files(pattern="^bam\\d+\\.csv\\.bz2$"), function(FNAME) {
  x <- read.table(FNAME, sep="|", header = TRUE)
  names(x)[1] <- substring(names(x)[1], 3)
  x <- x[-1,]
  x <- type.convert(x, as.is=TRUE)
  x$peri <- as.integer(sub("^bam(\\d+).*", "\\1", FNAME))
  x[setdiff(COL, names(x))] <- NA
  x[COL]
}))

#In peri 3 and 4 column BART is called BA
i <- which(.$peri %in% 3:4)
.$BART[i] <- .$BA[i]
.$BA <- NULL

#Remove intermediate points
. <- .[(.$RW + .$HW) %% 2 == 0,]
. <- .[.$PB %in% c(0,8,16,24),]

#Tree ID is in this data: RW HW PBFL AZI DIST
. <- .[do.call(order, .[c("RW", "HW", "PB", "AZI", "DIST", "peri")]),]

#Not needed. Assumption that as long as tree is measured it is not removed
#i <- do.call(paste, .[c("RW", "HW", "PB", "AZI", "DIST")])
#i <- c(i[-length(i)] != i[-1], FALSE)
#.$removed <- c((grepl("[0-9BFKN]", .$N) & !grepl("FALSE", .$N))[-1], FALSE)
#is.na(.$removed) <- i
##Observations of trees which have been declared to be removed
#i <- ave(.$removed, do.call(paste, .[c("RW", "HW", "PB", "AZI", "DIST")])
#       , FUN = function(x) {seq_along(x) <= min(length(x), which(x))} )

.$living <- is.na(.$D) | .$D == 0
#Live dead tree again?
i <- ave(.$living, do.call(paste, .[c("RW", "HW", "PB", "AZI", "DIST")])
       , FUN = function(x) {sum(x) >= min(length(x)+1, which(!x))} )

.$N <- NULL
.$D <- NULL
. <- .[!is.na(.$BHD),]

#Two rows for one period?
i <- ave(.$peri, do.call(paste, .[c("RW", "HW", "PB", "AZI", "DIST")]), FUN = anyDuplicated)
.[i>0,]
#Remove the duplicated rows
. <- .[!duplicated(do.call(paste, .[c("RW", "HW", "PB", "AZI", "DIST", "peri")])),]

#Missing periods?
i <- ave(.$peri, do.call(paste, .[c("RW", "HW", "PB", "AZI", "DIST")])
  , FUN = function(x) {
    y <- range(x)
    (1 + diff(y)) != (length(x) + (y[1] < 8 & y[2] > 8))
  } )
.[i==1,]
#Yes
#Currently it looks like that no new tree is entering at same AZI and DIST
#Observationgaps will be filled up automaticaly in the following process

#set crown heigth of 0 to NA
is.na(.$KRONHO) <- .$KRONHO == 0

#Fill missing heigt to crown base
i <- !.$living & is.na(.$KRONHO)
.$KRONHO[i] <- .$HOEHE[i]
#
.$kh1 <- .$KRONHO
for(GRP in list(c("RW", "HW", "PB", "AZI", "DIST")
          , c("RW", "HW", "PB", "BART", "peri")
          , c("RW", "HW", "PB", "BART")
          , c("RW", "HW", "PB", "peri")
          , c("RW", "HW", "PB")
          , "BART"
            ) ) {
  i <- as.factor(do.call(paste, .[GRP]))
  .$kh2 <- unsplit(lapply(split(.[c("HOEHE", "KRONHO")], i), function(x) {
    if(length(unique(x$KRONHO[!is.na(x$KRONHO)])) > 1 & length(unique(x$HOEHE[!is.na(x$KRONHO)])) > 1)
      as.integer(pmin(x$HOEHE, round(pmax(0, approxfun(x$HOEHE, x$KRONHO / x$HOEHE, rule=2)(x$HOEHE)) * x$HOEHE)))
    else x$KRONHO
  }), i)
  i <- which(is.na(.$kh1))
  .$kh1[i] <- .$kh2[i]
  .$kh2 <- NULL
}
#
i <- which(is.na(.$kh1))
.$kh1[i] <- as.integer(pmin(.$HOEHE[i], pmax(0, round(0.5 * .$HOEHE[i]))))
#
.$KRONHO <- .$kh1
.$kh1 <- NULL

#Add unique tree identifying numbers
#i <- as.factor(do.call(paste, .[c("RW", "HW", "PB")]))
#.$treeId <- unsplit(lapply(split(.[c("AZI", "DIST")], i), function(x) {
#  cumsum(!duplicated(x))
#}), i)
.$treeId <- do.call(paste, c(.[c("AZI", "DIST")], sep="_"))
.[c("AZI", "DIST")] <- NULL

.$plotId <- do.call(paste, c(.[c("RW", "HW", "PB")], sep="_"))
.[c("RW", "HW", "PB", "T")] <- NULL

#Change units
i <- c("BHD", "HOEHE", "KRONHO")
.[i] <- .[i]/10

#grep and nrep
.$grep <- 0
.$nrep <- 0
.$nrep[.$BHD < 10.4] <- 10000/2.6^2/pi
.$nrep[.$BHD > 39.08 & .$peri > 8] <- 10000/9.77^2/pi
.$grep[.$nrep == 0] <- 4


Ds <- do.call(rbind, lapply(unname(split(.[c("plotId", "treeId", "peri", "BART", "living")], paste(.$plotId, .$treeId))), function(x) {
  data.frame(x[c("plotId", "treeId")][1,], species = tail(x$BART, 1)
           , lastAlivePeri = x$peri[match(-1L, diff(x$living))]
           , whatOut = ".95st.8ba")
  }))
write.table(Ds, xzfile("speciesInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.frame(.[c("plotId", "treeId", "peri", "nrep", "grep")], d=.$BHD, h=.$HOEHE, hka=.$KRONHO), xzfile("treeInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)


COL <- c("peri", "RW", "HW", "PB", "DATUM", "KG_WALD")

. <- do.call(rbind, lapply(list.files(pattern="^prf\\d+\\.csv\\.bz2$"), function(FNAME) {
  x <- read.table(FNAME, sep="|", header = TRUE)
  names(x)[1] <- substring(names(x)[1], 3)
  x <- x[-1,]
  x <- type.convert(x, as.is=TRUE)
  x$peri <- as.integer(sub("^prf(\\d+).*", "\\1", FNAME))
  x[setdiff(COL, names(x))] <- NA
  x[COL]
}))

. <- .[grep("[0-9]", .$DATUM),]
i <- regexpr("\\d{2}$", .$DATUM)
j <- regmatches(.$DATUM, i)
regmatches(.$DATUM, i) <- paste0(19 + (as.integer(j) < 60), j)
.$DATUM <- as.Date(.$DATUM, "%d.%m.%Y")

.$plotId <- do.call(paste, c(.[c("RW", "HW", "PB")], sep="_"))
.[c("RW", "HW", "PB")] <- NULL

. <- .[ave(.$KG_WALD>0, .$plotId, FUN=any),]

write.table(data.frame(.["plotId"], date=.$DATUM, .["peri"], share=.$KG_WALD/10), xzfile("dateInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)




