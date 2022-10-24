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
alpha = 0.99 # initial value of alpha
n = 1 # number of individual
K = 1 # number of latent variable
opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000, print_level = 3) # optimizer
x0 <- c(2, rep(1,2*n), 1, 0, alpha) # initial parameters
# iter = 1 # number of experiments
input_size = 20
test_size = 100
#~~~~~ simulation 1: 1 dimension ~~~~~#

RMSE <- NULL

for(i in 1:1){
  
  MAT_ITER = 5000
  x <- seq(0,10, length.out = input_size); x_test <- sort(runif(test_size, 0, 10)) # replace it by input data
  y <- sin(x)  + rnorm(input_size, 0, 0.05) # replace it by output data
  #y <- (y-mean(y))/std(y)
  # plot(y~x)
  trains <- list(x); tests <- list(x_test); 
  z=seq( min(x) , max(x), length.out=10 ) ## pseudo-input
  
  now <- proc.time()
  x0 <- c( rand() , rep( 1*rand() ,2*n), 0.1, 0, alpha) # change the fourth argument to smaller number
  
  for(j in 1:10){
    print(j)
    print(alpha)

    one <- nloptr(x0=x0, eval_f = logL,eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 5), alpha), ub = c(rep(Inf, 5), alpha))
    alpha = alpha - 0.99/10
    H0 <- one$solution ## optimal solution
    H0[6] = alpha
    x0 = H0
    
  }
  
  print(alpha)

  later <- proc.time()-now
  
  r_cuu=cuu(z,z,H0);r_cfu=cfu(x, z, H0);r_cff=cff(x, x, H0)
  fuf <- r_cfu%*%mBCG_noT(r_cuu,t(r_cfu))
  BB <- H0[6] * fuf + diag(H0[2*n+2]^2, dim(fuf)[1]) + (1-H0[6]) * r_cff
  
  # xxx <- cfu(tests[[1]],z,H0[c(1,1+1,1+n+1)])
  
  r1_main=solve(BB,y); r2_main=mBCG_noT(r_cuu,t(r_cfu)) #r_cfu
  ypred = cfu(tests[[1]],z,H0[c(1,1+1,1+n+1)])%*%r2_main%*%r1_main 
  lines(ypred~x_test, col="red")
  # RMSE <- c(RMSE, sqrt( sum((ypred-sin(x_test))^2) ) )
  
}

# RMSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

