
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

dim = 4
input_size = 300
test_size = 1000

x <- cube.solid.random(dim, input_size)  # generate inputs
z <- cube.solid.random(dim, 10); z <- z$points # generate psedo-inputs, change number based on sample size
x <- x$points


x_test <- cube.solid.random(dim, test_size)
x_test <- x_test$points

y <- apply(x, 1, griewank) # generate outputs

input_size = dim(x)[1]
test_size = dim(x_test)[1]

### covariance functions and log-marginal likelihood ###

### here, u is the latent variable ###

### different length parameter for each dimension

cuu <- function(a,b,L){
  
  a = t(apply(a, 1, "/", L[3:(3+(dim-1))]^2))
  b = t(apply(b, 1, "/", L[3:(3+(dim-1))]^2))
  
  d <- plgp::distance(a,b) #plgp::distance is a square!
  # d <- t(apply(d, 1, "/", L[3:(3+(dim-1))]^2))
  (L[2]^2)*exp(-0.5*d)
}

cfu <- function(a,b,L){
  a = t(apply(a, 1, "/", L[3:(3+(dim-1))]^2))
  b = t(apply(b, 1, "/", L[3:(3+(dim-1))]^2))
  
  d <- plgp::distance(a,b) #plgp::distance is a square!
  #d <- t(apply(d, 1, "/", L[3:(3+(dim-1))]^2))
  (L[2]^2)*exp(-0.5*d)
}

cff <- function(a,b,L){
  a = t(apply(a, 1, "/", L[3:(3+(dim-1))]^2))
  b = t(apply(b, 1, "/", L[3:(3+(dim-1))]^2))
  d <- plgp::distance(a,b) #plgp::distance is a square!
  #d <- t(apply(d, 1, "/", L[3:(3+(dim-1))]^2))
  (L[2]^2)*exp(-0.5*d)
}

logL=function(H,fn)
{
  
  r_cuu <- cuu(z, z, H); r_cfu <- cfu(x, z, H); r_cff <- cff(x, x, H)
  fuf <- r_cfu%*%mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)
  B <- H[2+dim+1] * fuf + diag(H[1]^2, dim(fuf)[1]) + (1-H[2+dim+1]) * r_cff
  
  RHS <- cbind(y, r_cfu)
  
  B_inv_y <- mBCG_noT(B, RHS, maxiter = MAT_ITER)
  
  A <- diag(H[1]^2, dim(fuf)[1]) + (1-H[2+dim+1]) * r_cff + H[2+dim+1] * fuf
  ch <- chol(A + diag(0.05, ncol(A) ), pivot=TRUE, tol=1e-16)
  logdeter_A <- 2*(sum(log(diag(ch))))
  
  # logdeter <- log(det(  diag(1, ncol(r_cuu)) + H[6]* mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%RHS[,2:(ncol(RHS))]  )) + logdeter_A
  logdeter <- 0 + logdeter_A
  
  B2 <-  diag(1, ncol(r_cuu)) + (1-2*H[2+dim+1]) * mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%B_inv_y[,2:(ncol(RHS))]
  
  logdeter2 <- log(det(B2))
  
  # likelihood
  a <- -0.5*log(2*pi)*length(y) + 
    (-1/(2-2*H[2+dim+1]))*logdeter - 
    0.5*t(y)%*%B_inv_y[,1] + 
    (-H[2+dim+1]/(2*(1-H[2+dim+1]) ) )*( log(1/H[1]^2)*length(y) + logdeter2) 
  
  
  return(as.numeric(-a)) # neg-lilelihood
  
}

### gradient ###

logL_grad=function(H,fn) {return(nl.grad(H,fn))} ## gradient


x0 <- c( rand() , rand(), rep(rand(), dim), alpha) 

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
    try(one <- nloptr(x0=x0, eval_f = logL,eval_grad_f = logL_grad,opts = opts,fn = logL, lb = c(rep(-Inf, 6), alpha), ub = c(rep(Inf, 6), alpha)), silent=FALSE)
    H0 <- one$solution ## optimal solution
    H0[2+dim+1] = alpha
    x0 = H0
  }
  
  alpha = alpha - 0.99/20
  H0[2+dim+1] = alpha
  x0 = H0
  
  
}

# print(alpha) 

y_test <- apply(x_test, 1, griewank) # generate test output

ypred = cff(x_test, x_pool, H0)%*%solve(cff(x_pool, x_pool, H0)+ diag(H0[1]^2, input_size))%*%y_pool # prediction

sqrt(sum((y_test-ypred)^2))/sqrt(test_size) # RMSE
