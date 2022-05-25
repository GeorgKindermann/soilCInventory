dirOut <- "./dat/"
dirMeteo <- "/dat/meteo/zamg/spartacus/month/"
fileAltiMeteo <- "/dat/meteo/zamg/spartacus/SPARTACUS-MONTHLY_DEM_MASK.nc"
fileInvPoint <- "./dat/posInv.txt"
timeRange <- c(1981, 2020)

library(Rcpp)
library(ncdf4)

sourceCpp("kinterFix.cc")

unlink(paste0(dirOut, "tmInvPkt.txt"))
unlink(paste0(dirOut, "rrInvPkt.txt"))
tf <- tempfile()
L <- data.frame(what = c("Tm", "RR")
              , kiter = c("./kinterFix -o DIROUTtmInvPkt.txt -m WEATHERDATAIN -p fileInvPoint -s 24 -r 500 -fx 0.01 -fy 0.01 -fz 1 -gd -0.0047 -gl -0.01 -gh 0.002 -rl -25 -rh 30"
                        , "./kinterFix -o DIROUTrrInvPkt.txt -m WEATHERDATAIN -p fileInvPoint -s 24 -r 500 -fx 0.01 -fy 0.01 -fz 1 -gd 0.028 -gl 0. -gh 0.08 -rl 0 -rh 600"))
L$kiter <- sub("DIROUT", dirOut, L$kiter, fixed=TRUE)
L$kiter <- sub("WEATHERDATAIN", tf, L$kiter, fixed=TRUE)
L$kiter <- sub("fileInvPoint", fileInvPoint, L$kiter, fixed=TRUE)

nn <- ncvar_get(nc_open(fileAltiMeteo), "dem")
for(i in seq_len(nrow(L))) {
  FNS <- list.files(dirMeteo, paste0("SPARTACUS-MONTHLY_", L$what[i], "_\\d+.nc"))
  . <- as.integer(regmatches(FNS, regexpr("\\d{4}", FNS)))
  FNS <- FNS[. >= timeRange[1] & . <= timeRange[2]]
  for(FN in FNS) {
    nc <- nc_open(paste0(dirMeteo, FN))
    .  <- ncvar_get(nc, L$what[i])
    .x <- ncvar_get(nc, "x")
    .y <- ncvar_get(nc, "y")
    dat <- data.frame(x=.x, y=rep(.y, each=length(.x)), nn=as.vector(nn)
                    , date = rep(paste0(gsub("\\D+", "", FN), sprintf("%02d", 1:12)), each = length(.x) * length(.y))
                    , val=as.vector(.))
    dat <- dat[complete.cases(dat),]
    write.table(dat, tf, quote = FALSE, row.names = FALSE, col.names = FALSE)
    kinterFix(L$kiter[i])
  }
}




