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

res <- lapply(split(Dcc, Dcc$plotId), \(x) {
  years <- min(x$year):max(x$year)
  DT <- diff(years)
  DT <- c(DT, DT[length(DT)]) #Assuming last observation is as long as the previous
  .w <- Dw[as.character(x$plotId[1]), as.character(years),]
  CSoil <- rowSums(sapply(unique(x$d), \(d) {
    y <- x[x$d==d,]
    . <- as.matrix(y[match(years, y$year), c("a", "w", "e", "n")])
    .[is.na(.)] <- 0
    tt <- colMeans(.w)
    cc <- yasso.getSpin(theta, tt["t"], tt["p"], tt["amp"], d, leach, c(colMeans(.), 0)) #Spinnup with which data?
    do.call(rbind, Reduce(\(cc, i) {
      yasso.getNextTimestep(theta, .w[i, "t"], .w[i, "p"], .w[i, "amp"], d, leach, cc, c(.[i,], 0), DT[i])
    }, seq_along(years), cc, accumulate = TRUE)) |>
    rowSums()
  }))
  cbind(year=c(years[1], years[1] + cumsum(DT))-1, CSoil)
})
