fileYassoCC <- "./dat/yassoCC.txt"
fileYassoWetter <- "./dat/yassoWetter.txt"

source("y15_subroutine.r")
theta <- c(0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26)
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
DM <- sum(Dcc$d * .) / sum(.)
#Average inflow
. <- aggregate(cbind(a, w, e, n) ~ plotId + year, Dcc, sum)
t2 <- aggregate(cbind(a, w, e, n) ~ plotId, ., mean)
#Mittleres Wetter
t3 <- apply(Dw, c(1,3), median)
#
M <- proportions(colSums(t2[c("a", "w", "e", "n")]))
. <- rowSums(t2[c("a", "w", "e", "n")])
i <- match(t2$plotId, dimnames(t3)$plotId)
#a <- loess(log(.) ~ scale(t3[i,]))
a <- loess(log(.) ~ t3[i,])
BC <- mean(.) / mean(exp(predict(a)))
spInC <- data.frame(t2["plotId"], outer(exp(predict(a)) * BC, M))
spInW <- apply(Dw[,dimnames(Dw)$year < 1980,], c(1,3), median)


res <- lapply(split(Dcc, Dcc$plotId), \(x) {
  years <- min(x$year):max(x$year)
  DT <- diff(years)
  DT <- c(DT, DT[length(DT)]) #Assuming last observation is as long as the previous
  i <- match(x$plotId[1], dimnames(spInW)$plotId)
  avgCIn <- c(unlist(spInC[match(x$plotId[1], spInC$plotId),c("a", "w", "e", "n")]), 0)
  cc <- yasso.getSpin(theta, spInW[i, "t"], spInW[i, "p"], spInW[i, "amp"], DM, leach, avgCIn)
  #Run to climate of previous years with avgCIn
  .w <- Dw[as.character(x$plotId[1]), dimnames(Dw)$year < years[1],]
  DTP <- diff(as.numeric(dimnames(Dw)$year[dimnames(Dw)$year <= years[1]]))
  cc <- Reduce(\(cc, i) {
    yasso.getNextTimestep(theta, .w[i, "t"], .w[i, "p"], .w[i, "amp"], DM, leach, cc, avgCIn, DTP[i])
  }, seq_along(DTP), cc)
  #Fade out of this spin up c
  .w <- Dw[as.character(x$plotId[1]), as.character(years),]
  cc <- do.call(rbind, Reduce(\(cc, i) {
    yasso.getNextTimestep(theta, .w[i, "t"], .w[i, "p"], .w[i, "amp"], DM, leach, cc, numeric(5), DT[i])
  }, seq_along(years), cc, accumulate = TRUE)) |> rowSums()
  CSoil <- rowSums(cbind(cc, sapply(unique(x$d), \(d) {
    y <- x[x$d==d,]
    . <- as.matrix(y[match(years, y$year), c("a", "w", "e", "n")])
    .[is.na(.)] <- 0
    tt <- colMeans(.w)
    cc <- numeric(5)
    do.call(rbind, Reduce(\(cc, i) {
      yasso.getNextTimestep(theta, .w[i, "t"], .w[i, "p"], .w[i, "amp"], d, leach, cc, c(.[i,], 0), DT[i])
    }, seq_along(years), cc, accumulate = TRUE)) |>
    rowSums()
  })))
  cbind(year=c(years[1], years[1] + cumsum(DT))-1, CSoil)
})
