#SoilC when harvested trees stay in forest
t1 <- as.matrix(read.table("soilCOnPktC.txt.bz2", check.names = FALSE))
t2 <- as.matrix(read.table("share.txt.bz2", check.names = FALSE))
#SoilC equilibrium with infall from living trees of period 9
t3 <- as.matrix(read.table("soilCOnPktB.txt.bz2", check.names = FALSE))
t4 <- read.table("shareB.txt.bz2")[,1]

t5 <- sort(unique(c(rownames(t1), rownames(t3))))
t6 <- sort(unique(c(colnames(t1), colnames(t3))))
SC <- array(0, c(length(t5), length(t6), 4), list(plotId=t5, year=t6, what=c("C0", "W0", "Cl", "Wl")))
i <- match(t5, rownames(t1))
j <- match(t6, colnames(t1))
SC[,,1] <- t1[i,j]
SC[,,2] <- t2[i,j]
i <- match(t5, rownames(t3))
j <- match(t6, colnames(t3))
SC[,,3] <- t3[i,j]
SC[,,4] <- t4[i]
rm(t1, t2, t3, t4, t5, t6, i, j)


MF <- read.table("./2invData/managentformInv.txt.xz", header=TRUE)
table(MF$BA)
#1    .. Hochwald
# 0
#  1  .. Wirtschaftswald
#  2  .. Schutzwald mit Holznutzung
#  3  .. Schutzwald ohne Holznutzung
#   1 .. Holzboden
#   3 .. Strauchflaeche
. <- split(MF$BA, MF$plotId)
table(lengths(.))
.[lengths(.) > 1]
MFP <- MF[MF$plotId %in% names(.)[lengths(.) == 1],] #Only with one type


SaE <- read.table("./datThomas/SaE.txt.bz2")[,1]
SiEe <- read.table("./datThomas/SiEextensiv.txt.bz2")[,1]
intersect(SaE, SiEe)

i <- MFP$plotId %in% SaE
table(i, MFP$BA)
i <- MFP$plotId %in% SiEe
table(i, MFP$BA)
MFP$BA <- paste0(MFP$BA, c("","e")[i+1])
MFM <- MFP[match(dimnames(SC)$plotId, MFP$plotId),]


x <- rowsum(SC[,,3] * SC[,,4], MFM$BA, na.rm=TRUE) /
  rowsum(+!is.na(SC[,,3]) * SC[,,4], MFM$BA, na.rm=TRUE)
#
. <- x["1031",] / x["1021e",]
c1 <- mean(.)
plot(as.integer(names(.)), ., xlab="Year")
abline(h=c1)
#
. <- x["1033",] / x["1023",]
c3 <- mean(.)
plot(as.integer(names(.)), ., xlab="Year")
abline(h=c3)
#
colSums(!is.na(SC[,,1])) # -> 1986 - 2016
j <- dimnames(SC)$year %in% 1986:2016
i <- which(apply(!is.na(SC[,j,1]), 1, all))
x <- rowsum(SC[i,j,1] * SC[i,j,2], MFM$BA[i]) / rowsum(SC[i,j,2], MFM$BA[i])
x["1021e",] * c1
x["1023",] * c3
## No differenciation between Strauchflaeche
i <- pmin(2, match(MFM$BA, c("1021e", "1031", "1033")))
j <- which(!is.na(i))
x <- rowsum(SC[j,,3] * SC[j,,4], i[j], na.rm=TRUE) /
  rowsum(+!is.na(SC[j,,3]) * SC[j,,4], i[j], na.rm=TRUE)
#
. <- x["2",] / x["1",]
c1 <- mean(.)
plot(as.integer(names(.)), ., xlab="Year")
abline(h=c1)
#
k <- dimnames(SC)$year %in% 1986:2016
l <- apply(!is.na(SC[,k,1]), 1, all)
l <- whihc(l[j])
x <- rowsum(SC[l,k,1] * SC[l,k,2], i[l]) / rowsum(SC[l,k,2], i[l])
x["1",] * c1

