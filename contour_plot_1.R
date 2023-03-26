## This file is used to generate Figure 1 in our main paper
## Please adjust the values of alpha to get different contour plots

## global parameters ##

require(MASS)
require(nloptr) ## optimizer
require(pracma) ## numerical linear algebra
source("mBCG_noT.R")
source("logL.R")

MAT_ITER = 5000 # max number of iterations for mBCG
n = 1 # number of individual
# opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000, print_level = 3) # optimizer
# x0 <- c(2, rep(1,2*n), 1, 0, alpha) # initial parameters
iter = 1 # number of experiments
input_size = 50 # number of input data points
test_size = 200 # number of testing points

#~~~~~ simulation: data generation ~~~~~#

set.seed(2021)
x <- seq(0,5, length.out = input_size); x_test <- sort(runif(test_size, 0, 5)) # input data
sigma_f = 1.5 # GP variance scale parameter
sigma_e = 0.1 # GP additive noise variance parameter
length = 0.1 # GP length parameter
cov_matrix <- cff(x, x, c(length, sigma_f)) + diag(sigma_e^2, input_size) # covariance matrix
y <- output <- mvrnorm(1, rep(0, input_size), cov_matrix) # output data

### plot training data ###

plot(y~x, type="l")
trains <- list(x); tests <- list(x_test); 
z=seq( min(x) , max(x), length.out=10 ) ## pseudo-inputs or inducing points


### covariance functions ###

### u is latent variable ###

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


### log-marginal likelihood ###

alpha = 0 ## adjust alpha here and make new plot

logL=function(H, fn)
{
  # H[1] length parameter
  # H[2] variance scale parameter
  # H[3] additive noise variance parameter
  
  alpha = 0 ## adjust alpha here and make new plot
  
  r_cuu <- cuu(z, z, H); r_cfu <- cfu(x, z, H); r_cff <- cff(x, x, H) # calculate covariance matrices
  fuf <- r_cfu%*%mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER) # calculate low-rank approximation matrix Q
  B <- alpha * fuf + diag(H[3]^2, dim(fuf)[1]) + (1-alpha) * r_cff # calculate the covariance of alpha-ELBO objective
  
  RHS <- cbind(y, r_cfu)
  
  B_inv_y <- mBCG_noT(B, RHS, maxiter = MAT_ITER) # run conjugate method for efficient matrix product
  
  A <- diag(H[3]^2, dim(fuf)[1]) + (1-alpha) * r_cff + alpha * fuf # calculate matrix \xi in Eq. (6) defined in Appendix C
  ch <- chol(A + diag(0.1^2, ncol(A) ), pivot=TRUE, tol=1e-32) # cholesky decomposition 
  logdeter_A <- 2*(sum(log(diag(ch)))) # calculate log-determinant of \xi
  
  # logdeter <- log(det(  diag(1, ncol(r_cuu)) + H[6]* mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%RHS[,2:(ncol(RHS))]  )) + logdeter_A
  logdeter <- logdeter_A
  
  B2 <-  diag(1, ncol(r_cuu)) + (1-2*alpha) * mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%B_inv_y[,2:(ncol(RHS))] # calculate C_x in Appendix C
  
  logdeter2 <- log(det(B2)) # calculate the log determinant
  
  # return alpha-ELBO objective
  a <- -0.5*log(2*pi)*length(y) + 
    (-1/(2-2*alpha))*logdeter - 
    0.5*t(y)%*%B_inv_y[,1] + 
    (-alpha/(2*(1-alpha) ) )*( log(1/H[3]^2)*length(y) + logdeter2) 
  
  
  return(as.numeric(-a)) # negative alpha-ELBO
  
}

### with the first set of initialization, the optimizer converges to

set.seed(2021)
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000, print_level = 3) # optimizer
x0 <- c( 2*rand() , 2*rand(), 2*rand(), 0, 0, alpha) # change the fourth argument to smaller number
logL_grad=function(H,fn) {return(nl.grad(H,fn))} ## gradient
one <- nloptr(x0=x0, eval_f = logL, eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 5), alpha), ub = c(rep(Inf, 5), alpha))
H0 <- one$solution ## optimal solution

### with the second set of initialization, the optimizer converges to

set.seed(2023)
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000, print_level = 3) # optimizer
x0 <- c( 2*rand() , 2*rand(), 20*rand(), 0, 0, alpha) # change the fourth argument to smaller number
logL_grad=function(H,fn) {return(nl.grad(H,fn))} ## gradient
one <- nloptr(x0=x0, eval_f = logL, eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 5), alpha), ub = c(rep(Inf, 5), alpha))
H0 <- one$solution ## optimal solution

##################################
##### create a contour plot #####
##################################

# length_1 = 20
# length_2 = 20
# length_3 = 30
# H_1 <- seq(-1, 20, length.out = length_1)
# H_2 <- seq(-1, 3, length.out = length_2)
# H_3 <- seq(0.5, 10, length.out = length_3)
# H_input <- expand.grid(H_1, H_2, H_3)
# H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3, 3)
# H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3))
# H_outer_temp = H_outer
# # H_outer <- outer(X = H_1, Y = H_2, FUN = logL)
# 
# library(plotly)
# 
# H_outer = H_outer_temp
# 
# cols <- hcl.colors(25, "YlOrRd")
# cols <- hcl.colors(25, hcl.pals()[35])
# library(gplots)
# filled.contour(x = H_2, y = H_3, z = -(H_outer[1,,]), col=rev(cols)) # this plot log likelihood !

##################################

# the first set of initialization

length_1 = 1
length_2 = 70
length_3 = 50
H_1 <- seq(-0.102128, -0.102128, length.out = length_1)
H_2 <- seq(1.4, 2, length.out = length_2)
H_3 <- seq(-0.400000, 0.1, length.out = length_3)
H_input <- expand.grid(H_1, H_2, H_3)
H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3, 3)
H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3))
H_outer_temp = H_outer
# H_outer <- outer(X = H_1, Y = H_2, FUN = logL)

library(plotly)

H_outer = H_outer_temp

cols <- hcl.colors(25, "YlOrRd")
cols <- hcl.colors(25, hcl.pals()[35])
library(gplots)
filled.contour(x = H_2, y = H_3, z = -(H_outer[1,,]), col=rev(cols)) # this plot log likelihood !


##################################

# the second set of initialization

length_1 = 1
length_2 = 50
length_3 = 50
H_1 <- seq(0.883749, 0.883749, length.out = length_1)
H_2 <- seq(-0.1, 0.1, length.out = length_2)
H_3 <- seq(1.5, 1.8, length.out = length_3)
H_input <- expand.grid(H_1, H_2, H_3)
H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3, 3)
H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3))
H_outer_temp = H_outer
# H_outer <- outer(X = H_1, Y = H_2, FUN = logL)

library(plotly)

H_outer = H_outer_temp

cols <- hcl.colors(25, "YlOrRd")
cols <- hcl.colors(25, hcl.pals()[35])
library(gplots)
filled.contour(x = H_2, y = H_3, z = -(H_outer[1,,]), col=rev(cols)) # this plot log likelihood !


####################################

