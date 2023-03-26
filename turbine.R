## This file tests our algorithm using the real-life dataset. Please download the CMAPSSData file from our GitHub


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

#~~~~~~~~~ turbine data  ~~~~~~~~~~#

MAT_ITER = 5000 # max number of iterations for mBCG
alpha = 0.99 # initial value of alpha
n = 1 # number of individual
K = 1 # number of latent variable

setwd("CMAPSSData/CMAPSSData")
FD001 <- read.table("train_FD001.txt", head=FALSE)
Data <- FD001[(FD001[, 1]==1),]
range01 <- function(x){(x-min(x))/(max(x)-min(x)) + 0.0} 

x_total <- range01(Data[,2]) # input

y_total <- Data[,9] # output
mean_y = mean(y_total)
std_y = sqrt(var(y_total))
y_total <- (y_total-mean_y)/std_y

n_total = length(y_total) # total number of samples
input_size = round(n_total*0.6) # training size
test_size = n_total-input_size # testing size

set.seed(2022)
train_index <- sort(sample(c(1:n_total), size = input_size, replace = F))
test_index <- sort(setdiff( c(1:n_total),  train_index))

x <- x_total[train_index] # training inputs
y <- y_total[train_index] # training outputs
z <- seq(0, 1, length.out = 10)  # generate psedo-inputs, change number based on sample size



x_test <- x_total[test_index]
y_test <- y_total[test_index]

### covariance functions and log-marginal likelihood ###

### here, u is the latent variable ###


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
  
  ## alpha-ELBO objective
  
  r_cuu <- cuu(z, z, H); r_cfu <- cfu(x, z, H); r_cff <- cff(x, x, H)
  fuf <- r_cfu%*%mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)
  B <- H[6] * fuf + diag(H[2*n+2]^2, dim(fuf)[1]) + (1-H[6]) * r_cff
  
  RHS <- cbind(y, r_cfu)
  
  B_inv_y <- mBCG_noT(B, RHS, maxiter = MAT_ITER)
  
  A <- diag(H[2*n+2]^2, dim(fuf)[1]) + (1-H[6]) * r_cff + H[6] * fuf
  ch <- chol(A + diag(0.001, ncol(A) ), pivot=TRUE, tol=1e-32)
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

### BBMM prediction ###
K_plus_sigma = cff(x_pool, x_pool, H0)+ diag(H0[2*n+2]^2, input_size) 
K_plus_sigma_inv_y <- mBCG_noT(K_plus_sigma , as.matrix(y_pool), maxiter = MAT_ITER)
ypred = cff(x_test, x_pool, H0)%*%K_plus_sigma_inv_y
#######################

sqrt(sum((y_test-ypred)^2))/sqrt(test_size) # RMSE







