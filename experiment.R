
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
source("mBCG_noT.R") # mBCG algorithm

#~~~~~~~~~ GD-4, 4 dimension example  ~~~~~~~~~~#

MAT_ITER = 5000 # max number of iterations for mBCG
alpha = 0.99 # initial value of alpha
n = 1 # number of individual
K = 1 # number of latent variable

require(geozoo) # generate high dimensional data
source("griewank.R") # simulation data

dim = 4 # this is a 4 dimensional example
input_size = 300 # you can change input size here
test_size = 1000 # you can change test size here

x <- cube.solid.random(dim, input_size)  # generate inputs
z <- cube.solid.random(dim, 10); z <- z$points # generate psedo-inputs, change number based on sample size
x <- x$points



x_test <- cube.solid.random(dim, test_size)
x_test <- x_test$points

#
# for(i in 1:4){
#   x[,i] <- (x[,i]-mean(x[,i]))/std(x[,i])
#   z[,i] <- (z[,i]-mean(z[,i]))/std(z[,i])
#   # x_test[,i] <- (x_test[,i]-mean(x_test[,i]))/std(x_test[,i])
# }
#

y <- apply(x, 1, griewank) # generate outputs

## y <- (y-mean(y))/std(y) ##

input_size = dim(x)[1]
test_size = dim(x_test)[1]

### covariance functions and log-marginal likelihood ###

### here, u is the latent variable ###

### !!! Please refer to "different_length.R" if you want to use different length parameter for each dimension ###

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


x0 <- c( rand() , rand(), 0, rand(), 0, alpha) 

x_pool <- x
y_pool <- y

alpha = 0.99

for(j in 1:20){
  print(j)
  print(alpha)
  
  
  for(k in 1:30){
    
    set.seed(k+j+2020)
    
    ## take batched data at each iteration and optimize
    
    batch_index <- sample(c(21:(input_size-20)), 1)
    batch_index <- c( (batch_index-20) : (batch_index+20) )
    
    x <- x_pool[batch_index,]
    y <- y_pool[batch_index]
    opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2, print_level = 3) # optimizer
    try(one <- nloptr(x0=x0, eval_f = logL,eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 5), alpha), ub = c(rep(Inf, 5), alpha)), silent=TRUE)
    H0 <- one$solution ## optimal solution
    H0[6] = alpha
    x0 = H0
  }
  
  alpha = alpha - 0.99/20
  H0[6] = alpha
  x0 = H0
  
  
}

# print(alpha) 

y_test <- apply(x_test, 1, griewank) # generate test output

# for(i in 1:4) x_test[,i] <- (x_test[,i]-mean(x_test[,i]))/std(x_test[,i]) ##
# y_test <- (y_test-mean(y_test))/std(y_test) ## 

# regular prediction #
# ypred = cff(x_test, x_pool, H0)%*%solve(cff(x_pool, x_pool, H0)+ diag(H0[2*n+2]^2, input_size))%*%y_pool 
########################


### BBMM prediction ###
K_plus_sigma = cff(x_pool, x_pool, H0)+ diag(H0[2*n+2]^2, input_size) 
K_plus_sigma_inv_y <- mBCG_noT(K_plus_sigma , as.matrix(y_pool), maxiter = MAT_ITER)
ypred = cff(x_test, x_pool, H0)%*%K_plus_sigma_inv_y
#######################

sqrt(sum((y_test-ypred)^2))/sqrt(test_size) # RMSE

















#~~~~~~~~~~ GRAMACY & LEE  ~~~~~~~~~~#


## required packages and source code ##

# source("logL.R") # ccovariance functions and likelihood functions

source("grlee12.r") # simulation data
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





## global parameters ##

MAT_ITER = 5000 # max number of iterations for mBCG
alpha = 0.99 # initial value of alpha
n = 1 # number of individual
K = 1 # number of latent variable
# opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 100, print_level = 3) # optimizer
# x0 <- c(2, rep(1,2*n), 1, 0, alpha) # initial parameters
# iter = 1 # number of experiments
input_size = 300
test_size = 2000
#~~~~~ simulation 1: 1 dimension ~~~~~#

RMSE <- NULL


  
alpha = 0.99 # initial value of alpha
MAT_ITER = 5000
set.seed(2020)
x <- seq(0.5,1.5, length.out = input_size); x_test <- sort(runif(test_size, 0.5, 1.5)) # replace it by input data
y <- grlee12(x) + rnorm(input_size, 0, 0.05) # replace it by output data
y_test <- grlee12(x_test); #y_test <- (y_test-mean(y_test))/std(y_test)
#x_test <- (x_test-mean(x_test))/std(x_test)
#x <- (x-mean(x))/std(x)
#y <- (y-mean(y))/std(y)
plot(y~x)
trains <- list(x); tests <- list(x_test); 
z=seq( min(x) , max(x), length.out=10 ) ## pseudo-input
#z <- (z-mean(z))/std(z)

now <- proc.time()
set.seed(2000)
x0 <- c( rand() , rand(), 0, rand(), 0, alpha) # change the fourth argument to smaller number

x_pool <- x
y_pool <- y

for(j in 1:20){
  print(j)
  print(alpha)
  
  
  for(k in 1:30){
    
    # batch_index = c( (1+50*(k-1)):(50*k)  ) 
    set.seed(k+j+2020)
    ## take batched data at each iteration and optimize
    batch_index <- sample(c(21:(input_size-20)), 1)
    batch_index <- c( (batch_index-20) : (batch_index+20) )
    # batch_index[(batch_index<1)] <- 1
    # batch_index[(batch_index>input_size)] <- input_size
    x <- x_pool[batch_index]
    y <- y_pool[batch_index]
    opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2, print_level = 3) # optimizer
    try(one <- nloptr(x0=x0, eval_f = logL,eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 5), alpha), ub = c(rep(Inf, 5), alpha)), silent=TRUE)
    H0 <- one$solution ## optimal solution
    H0[6] = alpha
    x0 = H0
  }
  
  alpha = alpha - 0.99/20
  H0[6] = alpha
  x0 = H0

  
}

print(alpha)

ypred = cff(x_test, x_pool, H0)%*%solve(cff(x_pool, x_pool, H0)+ diag(H0[2*n+2]^2, input_size))%*%y_pool








