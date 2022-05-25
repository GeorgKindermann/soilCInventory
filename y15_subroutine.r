library(Matrix)

yasso15SetA <- function(theta, avgT, sumP, ampT, diam, leach) {
  tabs <- c(1:4, 32, 35)
  theta[tabs] = -abs(theta[tabs])

  A <- matrix(0, 5, 5)

  m3 <- sumP/1000
#temperature annual cycle approximation
  m0 <- (1./sqrt(2.)-1.)/pi
  m1 <- 1./sqrt(2.)/pi
  m2 <- (1.-1./sqrt(2.))/pi
  clim <- 4. * ampT
  te <- 1:4
  te[1] = avgT + clim * m0
  te[2] = avgT - clim * m1
  te[3] = avgT + clim * m2
  te[4] = avgT + clim * m1

#Average temperature dependence
  tem = sum(exp(theta[22]*te+theta[23]*te^2)) / 4
  temN = sum(exp(theta[24]*te+theta[25]*te^2)) / 4
  temH = sum(exp(theta[26]*te+theta[27]*te^2)) / 4

#Precipitation dependence
  tem = tem * (1.-exp(theta[28] * m3));
  temN = temN * (1.-exp(theta[29] * m3));
  temH = temH * (1.-exp(theta[30] * m3));

#Size class dependence -- no effect if d == 0.0
  size_dep <- 1
  if(diam > 0.) {size_dep <- min(1., (1. + theta[33]*diam + theta[34]*diam^2)^theta[35])}
#Calculating matrix a (will work ok despite the sign of alphas)
  A[1+(0:2*6)] = theta[1:3]*tem*size_dep
  A[4,4] = theta[4]*temN*size_dep
  A[5,5] = theta[32]*temH #no size effect in humus
  dAbs <- abs(A[1+(0:3*6)])
  idx <- 5
  for(i in 0:3) {
    for(j in 0:3) {
      if(i!=j) {
        A[1+i*5+j] = theta[idx] * dAbs[1+j];
        idx <- idx+1
      }
    }
  }
#mass flows AWEN -> H (size effect is present here)
  A[20+1:4] = theta[31] * dAbs
#Leaching (no leaching for humus) 
  if(leach < 0.) {A[1+6*0:3] = A[1+6*0:3] + leach * m3}
  t(A)
}

yasso.getSpin <- function(theta, avgT, sumP, ampT, diam, leach, infall) {
  A <- yasso15SetA(theta, avgT, sumP, ampT, diam, leach)
  solve(A, infall) * -1
}

yasso.getNextTimestep <- function(theta, avgT, sumP, ampT, diam, leach, init, infall, time) {
  A <- yasso15SetA(theta, avgT, sumP, ampT, diam, leach)
  z1 <- A %*% init
  z1 <- z1 + infall
  At <- A * time
  z2 <- Matrix::expm(At) %*% z1
#Alternatives 
#library(expm)
#z2 <- expAtv(At,z1)$eAtv
#
#library(expoRkit)
#z2 <- expv(At, z1)
  z2 <- z2 - infall
  as.vector(solve(A, z2))
}
