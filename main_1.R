#~~~~~~~~~~ code for submission ~~~~~~~~~~#

## required packages and source code ##

source("logL.R") # ccovariance functions and likelihood functions
source("mBCG_noT.R") # mBCG algorithm
source("grlee12.r") # simulation data

require(Matrix) ## create diagonal matrix
require(nloptr) ## optimizer
require(numDeriv) ## gradient
require(MASS) ## basic package
library(survival) ## survival model
library(mnormt) ## multivaraite normal
library(nlme) ## mixed effect model
library(glmmML) ## Gauss-Hermite quadrature
library(AlgDesign)  ## factorial design
library(statmod)  ## Gauss-Legendre quadrature
require(ggplot2) # plots
library(pracma) ## Trace function
library(rootSolve) ## root solver
library(pROC) ## roc curve
require(boot) ## glm
library(e1071) ## svm
library(sendmailR) ## send email to you when code is finished
require(dplyr) ## data manipulation
require(doBy) ## collapse data

## global parameters ##

MAT_ITER = 5000 # max number of iterations for mBCG
alpha = 0.5 # value of alpha
n = 1 # number of individual
K = 1 # number of latent variable
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000, print_level = 3) # optimizer
x0 <- c(2, rep(1,2*n), 1, 0, alpha) # initial parameters
iter = 1 # number of experiments

#~~~~~ simulation 1: 1 dimension ~~~~~#

RMSE <- NULL

for(i in 1:iter){
  
  MAT_ITER = 5000
  x <- seq(0.5,2.5, length.out = 100)
  y <- grlee12(x)
  #y <- (y-mean(y))/std(y)
  plot(y~x)
  trains <- list(x)
  z=seq( min(x) , max(x), length.out=10 ) ## pseudo-input
  
  now <- proc.time()
  x0 <- c( 2*rand() , rep( rand() ,2*n), rand(), 0, alpha)
  one <- nloptr(x0=x0,eval_f = logL,eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 5), alpha), ub = c(rep(Inf, 5), alpha) )
  H0 <- one$solution ## optimal solution 
  later <- proc.time()-now
  
  r_cuu=cuu(z,z,H0);r_cfu=cfu(x, z, H0);r_cff=cff(x, x, H0)
  fuf <- r_cfu%*%mBCG_noT(r_cuu,t(r_cfu))
  BB <- H0[6] * fuf + diag(H0[2*n+2]^2, dim(fuf)[1]) + (1-H0[6]) * r_cff
  r1_main=solve(BB,y); r2_main=mBCG_noT(r_cuu,t(r_cfu))
  ypred = cfu(trains[[1]],z,H0[c(1,1+1,1+n+1)])%*%r2_main%*%r1_main 
  lines(ypred~x, col="red")
  RMSE <- c(RMSE, sqrt( sum((ypred-y)^2) ) )
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

