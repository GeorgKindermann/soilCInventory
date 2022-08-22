fileYassoCC <- "./dat/yassoCCB.txt"
fileYassoWetter <- "./dat/yassoWetter.txt"

STARYEAR <- 1980

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

library(Rcpp)
sourceCpp("y15c_subroutine.cc")
yassoTimespan(1)

res <- lapply(split(Dcc, Dcc$plotId), \(x) {
#x <- split(Dcc, Dcc$plotId)[[1]]
  i <- match(x$plotId[1], dimnames(Dw)$plotId)
  .w <- Dw[i, dimnames(Dw)$year >= STARYEAR,]
  rowSums(sapply(split(x, x$d), \(y) {
    #y <- split(x, x$d)[[1]]
    CIn <- c(unlist(y[1,c("a", "w", "e", "n")]), h=0)
    sapply(seq_len(nrow(.w)), function(i) {
      yassoSet(.w[i, "t"], .w[i, "p"], .w[i, "amp"], y$d, leach)
      sum(yassoSpin(CIn))
    })
  }))
})

YR <- dimnames(Dw)$year[dimnames(Dw)$year >= STARYEAR]

. <- do.call(rbind, res)
colnames(.) <- YR
. <- round(., 2)

Dd <- read.table("./2invData/dateInv.txt.xz", header=TRUE)
Dd <- Dd[Dd$peri==9, c("plotId", "share")]
tt <- Dd$share[match(rownames(.), Dd$plotId)]
#tt <- Dd[match(rownames(.), Dd$plotId), "share", drop=FALSE]

write.table(., bzfile("soilCOnPktB.txt.bz2"), quote = FALSE, row.names = TRUE, col.names = TRUE)

write(tt, bzfile("shareB.txt.bz2"), 1)
