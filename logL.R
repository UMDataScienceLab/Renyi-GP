#~~~~~~~~~~ covariance functions and likelihood functions ~~~~~~~~~~#

require(plgp) # construct covariance matrix with D dimension inputs


### covariance functions ###

cuu <- function(a,b,L){
  d <- plgp::distance(a,b) #plgp::distance is a square!
  res=L[1]^2; (L[2]^2)*exp(-0.5*d/res)
}

cfu <- function(a,b,L){
  d <- plgp::distance(a,b) #plgp::distance is a square!
  res=L[1]^2; (L[2]^2)*exp(-0.5*d/res)
}


cff <- function(a,b,L){
  d <- plgp::distance(a,b) #plgp::distance is a square!
  res=L[1]^2; (L[2]^2)*exp(-0.5*d/res)
}


# cfu=function(a,b,L) {d <- plgp::distance(a,b); res=L[1]^2+L[3]^2; L[2]*sqrt(L[1]^2/res)*exp(-0.5*d/res)}
# 
# cff=function(a,b,L) {
#   d=plgp::distance(a,b)
#   res=L[1]^2+2*L[3]^2
#   L[2]^2 * sqrt(L[1]^2/res) * exp(-0.5*d/res)
# }


### log-marginal likelihood ###

logL=function(H,fn)
{
  
  r_cuu <- cuu(z, z, H); r_cfu <- cfu(x, z, H); r_cff <- cff(x, x, H)
  fuf <- r_cfu%*%mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)
  B <- H[6] * fuf + diag(H[2*n+2]^2, dim(fuf)[1]) + (1-H[6]) * r_cff
  
  RHS <- cbind(y, r_cfu)
  
  B_inv_y <- mBCG_noT(B, RHS, maxiter = MAT_ITER)
  
  A <- diag(H[2*n+2]^2, dim(fuf)[1]) + (1-H[6]) * r_cff + H[6] * fuf
  ch <- chol(A + diag(0, ncol(A) ), pivot=TRUE, tol=1e-16)
  logdeter_A <- 2*(sum(log(diag(ch))))
  
  # logdeter <- log(det(  diag(1, ncol(r_cuu)) + H[6]* mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%RHS[,2:(ncol(RHS))]  )) + logdeter_A
  logdeter <- 0 + logdeter_A
  
  B2 <-  diag(1, ncol(r_cuu)) + (1-2*H[6]) * mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%B_inv_y[,2:(ncol(RHS))]

  logdeter2 <- log(det(B2))
  
  # likelihood
  a <- -0.5*log(2*pi)*length(y) + 
    (-1/(2-2*H[6]))*logdeter - 
    0.5*t(y)%*%B_inv_y[,1] + 
    (-H[6]/(2*(1-H[6]) ) )*( log(1/H[4]^2)*length(y) + logdeter2) 
  
  
  return(as.numeric(-a)) # neg-lilelihood
  
}

### gradient ###

logL_grad=function(H,fn) {return(nl.grad(H,fn))} ## gradient
