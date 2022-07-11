x <- as.matrix(read.table("soilCOnPkt.txt.bz2", check.names = FALSE))
y <- as.matrix(read.table("share.txt.bz2", check.names = FALSE))

#Only for points which have observations in all years
i <- which(apply(!is.na(x), 1, all))
. <- colSums(x[i,] * y[i,]) / colSums(y[i,])
. <- . * 120 / .["1990"]
plot(names(.), ., type="l", xlab="Jahr", ylab="Boden C [tC/ha]", main="Pkt mit durchgehender Beobachtung", ylim=c(0,max(.)))

#All points with an observations in the specific year
. <- colSums(x * y, TRUE) / colSums(!is.na(x) * y)
. <- . * 120 / .["1990"]
plot(names(.), ., type="l", xlab="Jahr", ylab="Boden C [tC/ha]", main="Pkt mit Beobachtung im Jahr", ylim=c(0,max(.)))

