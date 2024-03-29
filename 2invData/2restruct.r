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
           , lastAlivePeri = if(x$living[1]) x$peri[match(-1L, diff(x$living))] else x$peri[1] - 1
           , whatOut = ".95st.8ba")
  }))

Dt <- data.frame(.[c("plotId", "treeId", "peri", "nrep", "grep")], d=.$BHD, h=.$HOEHE, hka=.$KRONHO)


COL <- c("peri", "RW", "HW", "PB", "T", "B", "BA", "KG_WALD")

#B: (peri 3-4) Betriebsart: Hochwald - 1..Wirtschaftswald, 3..Schutzwald im Ertrag, 4..Schutzwald ausser Ertrag, 5..Holzboden ausser Ertrag; Ausschlagwald - 6..Land-Ausschlagwald, 7..Auen-Ausschlagwald, 8..Holzboden ausser Ertrag

#BA: (peri > 4) Betriebsart: ABCD; A=1..Hochwald, A=2..Ausschlagwald, A=3..(peri>8) Nichtwald (Codierungen für BCD haben andere Bedeutung), B=0..Keine Zuordnung möglich, B=1..Land-Ausschlagwald, B=2..Auen-Ausschlagwald, C=1..Wirtschaftswald, C=2..Schutzwald im Ertrag, C=3..Schutzwald außer Ertrag, D=1..Holzboden, D=2..Holzboden außer Ertrag, D=3..Strauchfläche, D=4..unbegehbarer Schutzwald außer Ertrag

. <- do.call(rbind, lapply(list.files(pattern="^tfl\\d+\\.csv\\.bz2$"), function(FNAME) {
  x <- read.table(FNAME, sep="|", header = TRUE)
  names(x)[1] <- substring(names(x)[1], 3)
  x <- x[-1,]
  x <- type.convert(x, as.is=TRUE)
  x$peri <- as.integer(sub("^tfl(\\d+).*", "\\1", FNAME))
  x[setdiff(COL, names(x))] <- NA
  x[COL]
}))
.$B[.$peri > 4] <- NA
.$sur <- .$B %in% c(1, 3, 6, 7) |
.$BA %in% c(1011, 1021, 2111, 2211) | 
(.$peri > 8 & .$BA %in% c(1031, 1013, 1023, 1033, 2113, 2213))
.$plotId <- do.call(paste, c(.[c("RW", "HW", "PB")], sep="_"))
. <- aggregate(KG_WALD ~ peri + plotId + T, .[.$sur,], mean)
Da <- aggregate(cbind(share = KG_WALD/10) ~ peri + plotId, .[.$KG_WALD > 0,], sum)


COL <- c("peri", "RW", "HW", "PB", "DATUM")

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

#Remove intermediate points
. <- .[(.$RW + .$HW) %% 2 == 0,]
. <- .[.$PB %in% c(0,8,16,24),]

.$plotId <- do.call(paste, c(.[c("RW", "HW", "PB")], sep="_"))
.[c("RW", "HW", "PB")] <- NULL

s <- paste(.$peri, .$plotId)
i <- which(duplicated(s))
.[s %in% s[i],]
. <- .[-i,]

. <- merge(., Da)

tt <- table(.$plotId)
. <- .[.$plotId %in% names(tt)[tt > 1],]

Dd <- data.frame(.["plotId"], .["peri"], date=.$DATUM, .["share"])

i <- match(paste(Dt$plotId, Dt$peri), paste(Dd$plotId, Dd$peri))
j <- which(!is.na(i))
Dt <- Dt[j,]
Dt[, c("nrep", "grep")] <- Dt[, c("nrep", "grep")] / Dd$share[i[j]]

Ds <- Ds[paste(Ds$plotId, Ds$treeId) %in% unique(paste(Dt$plotId, Dt$treeId)),]

write.table(Ds, xzfile("speciesInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(Dt, xzfile("treeInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(Dd, xzfile("dateInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)


x <- read.table("tfl9.csv.bz2", sep="|", header = TRUE)
x <- x[-1,]
x <- type.convert(x, as.is=TRUE)
x <- x[c("RW", "HW", "PB", "T", "BA", "KG_WALD")]

x <- x[x$BA %in% c(1011, 1013,1021, 1023, 1031, 1033) & x$KG_WALD > 0, ]
x$plotId <- do.call(paste, c(x[c("RW", "HW", "PB")], sep="_"))
x <- aggregate(KG_WALD ~ plotId + BA, x, sum)
write.table(x, xzfile("managentformInv.txt.xz"), quote = FALSE, row.names = FALSE, col.names = TRUE)

table(x$BA)
. <- split(x$BA, x$plotId)
table(lengths(.))
.[lengths(.) == 2]
