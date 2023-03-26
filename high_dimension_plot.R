## This file is used to generate Figures 1 and 2 in Appendix

require(MASS)
require(nloptr) ## optimizer
require(pracma) ##
source("mBCG_noT.R")

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
# alpha = 0.99 # initial value of alpha
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

logL=function(H)
{
  alpha = 0
  
  r_cuu <- cuu(z, z, H); r_cfu <- cfu(x, z, H); r_cff <- cff(x, x, H)
  fuf <- r_cfu%*%mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)
  B <- alpha * fuf + diag(H[1]^2, dim(fuf)[1]) + (1-alpha) * r_cff
  
  RHS <- cbind(y, r_cfu)
  
  B_inv_y <- mBCG_noT(B, RHS, maxiter = MAT_ITER)
  
  A <- diag(H[1]^2, dim(fuf)[1]) + (1-alpha) * r_cff + alpha * fuf
  ch <- chol(A + diag(0.05, ncol(A) ), pivot=TRUE, tol=1e-32)
  logdeter_A <- 2*(sum(log(diag(ch))))
  
  # logdeter <- log(det(  diag(1, ncol(r_cuu)) + H[6]* mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%RHS[,2:(ncol(RHS))]  )) + logdeter_A
  logdeter <- 0 + logdeter_A
  
  B2 <-  diag(1, ncol(r_cuu)) + (1-2*alpha) * mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%B_inv_y[,2:(ncol(RHS))]
  
  logdeter2 <- log(det(B2))
  
  # likelihood
  a <- -0.5*log(2*pi)*length(y) + 
    (-1/(2-2*alpha))*logdeter - 
    0.5*t(y)%*%B_inv_y[,1] + 
    (-alpha/(2*(1-alpha) ) )*( log(1/H[1]^2)*length(y) + logdeter2) 
  
  
  return(as.numeric(-a)) # neg-lilelihood
  
}

### gradient ###

# logL_grad=function(H,fn) {return(nl.grad(H,fn))} ## gradient


###################
# length_1 = 1
# length_2 = 1
# length_3 = 1
# length_4 = 1
# length_5 = 5
# length_6 = 5
# H_1 <- seq(-0.510411, length.out = length_1)
# H_2 <- seq(0.788280, length.out = length_2)
# H_3 <- seq(0.796116, length.out = length_3)
# H_4 <- seq(1.999872, length.out = length_4)
# H_5 <- seq(1, 3, length.out = length_5)
# H_6 <- seq(1, 3, length.out = length_6)
# H_input <- expand.grid(H_1, H_2, H_3, H_4, H_5, H_6)
# H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3*length_4*length_5*length_6, 6)
# H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3,length_4,length_5, length_6))
# H_outer_temp = H_outer
# # H_outer <- outer(X = H_1, Y = H_2, FUN = logL)
# 
# library(plotly)
# # surface contour
# 
# 
# # fig <- plot_ly(x = H_1, y = H_2, z = -H_outer[,,3], type = "contour")
# # fig #%>% layout(xaxis= list(showticklabels = FALSE))
# 
# H_outer = H_outer_temp
# 
# 
# cols <- hcl.colors(25, "YlOrRd")
# cols <- hcl.colors(25, hcl.pals()[35])
# library(gplots)
# filled.contour(x = H_5, y = H_6, z = -(H_outer[1,,,,,]), col=rev(cols)) # this plot log likelihood !
# 
# # , col = colorpanel(30, "white", "grey10")

###################

# sol x = ( -0.050914, -0.063265, 2.881757, 2.996898, 3.070571, 3.143977, 0.049500 )


length_1 = 1
length_2 = 1
length_3 = 1
length_4 = 1
length_5 = 5
length_6 = 5
H_1 <- seq(-0.050914, length.out = length_1)
H_2 <- seq(-0.063265, length.out = length_2)
H_3 <- seq(2.881757, length.out = length_3)
H_4 <- seq(2.996898, length.out = length_4)
H_5 <- seq(0.5, 4, length.out = length_5)
H_6 <- seq(0.5, 4, length.out = length_6)
H_input <- expand.grid(H_1, H_2, H_3, H_4, H_5, H_6)
H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3*length_4*length_5*length_6, 6)
H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3,length_4,length_5, length_6))
H_outer_temp = H_outer
# H_outer <- outer(X = H_1, Y = H_2, FUN = logL)

library(plotly)
# surface contour


# fig <- plot_ly(x = H_1, y = H_2, z = -H_outer[,,3], type = "contour")
# fig #%>% layout(xaxis= list(showticklabels = FALSE))

H_outer = H_outer_temp


cols <- hcl.colors(30, "YlOrRd")
cols <- hcl.colors(30, hcl.pals()[35])
library(gplots)
filled.contour(x = H_5, y = H_6, z = -(H_outer[1,,,,,]), col=rev(cols)) # this plot log likelihood !


########################################

# sol x = ( -0.940909, -0.690089, 7.613586, 7.631193, 7.641433, 7.581369, 0.049500 )


length_1 = 1
length_2 = 1
length_3 = 1
length_4 = 1
length_5 = 5
length_6 = 5
H_1 <- seq(-0.940909, length.out = length_1)
H_2 <- seq(-0.690089, length.out = length_2)
H_3 <- seq( 7.613586, length.out = length_3)
H_4 <- seq(7.631193, length.out = length_4)
H_5 <- seq(5, 11, length.out = length_5)
H_6 <- seq(5, 11, length.out = length_6)
H_input <- expand.grid(H_1, H_2, H_3, H_4, H_5, H_6)
H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3*length_4*length_5*length_6, 6)
H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3,length_4,length_5, length_6))
H_outer_temp = H_outer
# H_outer <- outer(X = H_1, Y = H_2, FUN = logL)

library(plotly)
# surface contour


# fig <- plot_ly(x = H_1, y = H_2, z = -H_outer[,,3], type = "contour")
# fig #%>% layout(xaxis= list(showticklabels = FALSE))

H_outer = H_outer_temp


cols <- hcl.colors(30, "YlOrRd")
cols <- hcl.colors(30, hcl.pals()[35])
library(gplots)
filled.contour(x = H_5, y = H_6, z = -(H_outer[1,,,,,]), col=rev(cols)) # this plot log likelihood !
