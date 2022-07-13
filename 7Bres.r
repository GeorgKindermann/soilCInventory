x <- as.matrix(read.table("soilCOnPktB.txt.bz2", check.names = FALSE))
y <- read.table("shareB.txt.bz2")[,1]

. <- colSums(x * y) / sum(y)
. <- . * 120 / .["1990"]
plot(names(.), ., type="l", xlab="Jahr", ylab="Boden C [tC/ha]", main="Mit statischem Input von ÖWI-Periode 9", ylim=c(0,max(.)))

