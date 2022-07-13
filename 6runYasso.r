fileYassoCC <- "./dat/yassoCC.txt"
fileYassoWetter <- "./dat/yassoWetter.txt"

leach <- 0

Dw <- read.table(fileYassoWetter, header=TRUE)
. <- sapply(c("t", "p", "amp"), \(s) {
  . <- regexpr(paste0("(?<=^", s, ")\\d+$"), names(Dw), perl=TRUE)
  y <- as.integer(regmatches(names(Dw), .))
  setNames(which(.>0)[order(y)], y)
})
Dw <- array(unlist(Dw[as.vector(.)], FALSE, FALSE)
          , c(nrow(Dw), nrow(.), ncol(.))
    , list(plotId = rownames(Dw), year = rownames(.), what = colnames(.)))


Dcc <- read.table(fileYassoCC, header=TRUE)


#Average dimension
. <- rowSums(Dcc[c("a", "w", "e", "n")])
#t1 <- tapply(Dcc$d * ., Dcc$plotId, sum) / tapply(., Dcc$plotId, sum)
DM <- (sum(Dcc$d^.5 * .) / sum(.))^2
#Average inflow
. <- aggregate(cbind(cIn = rowSums(cbind(a, w, e, n))) ~ plotId + year, Dcc, sum)
t1 <- aggregate(year ~ plotId, ., min)
t2 <- .[!paste(.$plotId, .$year) %in% paste(t1$plotId, t1$year),] #Skip first year
#Mittleres Wetter
t3 <- apply(Dw, c(1,3), median)
i <- match(t2$plotId, dimnames(t3)$plotId)
#a <- loess(t2$cIn ~ t2$year + t3[i,])
. <- data.frame(t2[c("cIn", "year")], t3[i,])
a <- glm(cIn ~ year + t + p + amp, ., family=quasipoisson(link = "log"))
M <- proportions(colSums(Dcc[c("a", "w", "e", "n")]))
i <- match(t1$plotId, dimnames(t3)$plotId)
spInC <- data.frame(t1["plotId"], outer(predict(a, data.frame(t1["year"], t3[i,]), type="response"), M))
spInW <- apply(Dw[,dimnames(Dw)$year < 1980,], c(1,3), median)


library(Rcpp)
sourceCpp("y15c_subroutine.cc")
res <- lapply(split(Dcc, Dcc$plotId), \(x) {
  years <- (min(x$year):max(x$year))[-1]
  DT <- diff(years)
  DT <- c(DT, DT[length(DT)]) #Assuming last observation is as long as the previous
  i <- match(x$plotId[1], dimnames(spInW)$plotId)
  avgCIn <- c(unlist(spInC[match(x$plotId[1], spInC$plotId),c("a", "w", "e", "n")]), 0)
  yassoSet(spInW[i, "t"], spInW[i, "p"], spInW[i, "amp"], DM, leach)
  yassoTimespan(1)
  cc <- yassoSpin(avgCIn)
  #Run using climate of previous years with avgCIn
  .w <- Dw[as.character(x$plotId[1]), dimnames(Dw)$year < years[1],]
  DTP <- diff(as.numeric(dimnames(Dw)$year[dimnames(Dw)$year <= years[1]]))
  cc <- Reduce(\(cc, i) {
    yassoSet(.w[i, "t"], .w[i, "p"], .w[i, "amp"], DM, leach)
    yassoTimespan(DTP[i])
    yassoNext(cc, avgCIn)
  }, seq_along(DTP), cc)
  #Fade out of this spin up c
  .w <- Dw[as.character(x$plotId[1]), as.character(years),]
  cc <- do.call(rbind, Reduce(\(cc, i) {
    yassoSet(.w[i, "t"], .w[i, "p"], .w[i, "amp"], DM, leach)
    yassoTimespan(DT[i])
    yassoNext(cc, numeric(5))
  }, seq_along(years), cc, accumulate = TRUE)) |> rowSums()
  CSoil <- rowSums(cbind(cc, sapply(unique(x$d), \(d) {
    y <- x[x$d==d,]
    . <- as.matrix(y[match(years, y$year), c("a", "w", "e", "n")])
    .[is.na(.)] <- 0
    do.call(rbind, Reduce(\(cc, i) {
      yassoSet(.w[i, "t"], .w[i, "p"], .w[i, "amp"], d, leach)
      yassoTimespan(DT[i])
      yassoNext(cc, c(.[i,], 0))
    }, seq_along(years), numeric(5), accumulate = TRUE)) |>
    rowSums()
  })))
  cbind(year=c(years[1], years[1] + cumsum(DT))-1, CSoil)
})

YR <- range(sapply(res, function(x) range(x[,"year"])))
YR <- YR[1]:YR[2]

. <- do.call(rbind, lapply(res, function(x) x[match(YR, x[,"year"]), "CSoil"]))
colnames(.) <- YR
. <- round(., 2)

Dd <- read.table("./2invData/dateInv.txt.xz", header=TRUE)
Dd$date <- as.Date(Dd$date)
tt <- split(Dd[c("date", "share")], Dd$plotId)

tt <- do.call(rbind, lapply(rownames(.), function(i) {
  x <- tt[[i]]
  approxfun(x$date, x$share, rule=2)(as.Date(paste0(YR, "-1-1")))
}))

write.table(., bzfile("soilCOnPkt.txt.bz2"), quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(tt, bzfile("share.txt.bz2"), quote = FALSE, row.names = FALSE, col.names = FALSE)
