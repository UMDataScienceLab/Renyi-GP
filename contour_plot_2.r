## global parameters ##

require(MASS)
source("mBCG_noT.R")
source("logL.R")

MAT_ITER = 5000 # max number of iterations for mBCG
n = 1 # number of individual
K = 1 # number of latent variable
# opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 10000, print_level = 3) # optimizer
# x0 <- c(2, rep(1,2*n), 1, 0, alpha) # initial parameters
iter = 1 # number of experiments
input_size = 50
test_size = 200
#~~~~~ simulation: data generation ~~~~~#

set.seed(2021)
# MAT_ITER = 5000
x <- seq(0,5, length.out = input_size); x_test <- sort(runif(test_size, 0, 5)) # replace it by input data
sigma_f = 1.5 # GP parameter
sigma_e = 0.01 # GP parameter
length = 0.01 # GP parameter
cov_matrix <- cff(x, x, c(length, sigma_f)) + diag(sigma_e^2, input_size) # covariance matrix ef
y <- output <- mvrnorm(1, rep(0, input_size), cov_matrix) # output data

### plot training data ###

plot(y~x, type="l")
trains <- list(x); tests <- list(x_test); 
z=seq( min(x) , max(x), length.out=10 ) ## pseudo-input


### covariance functions ###
### u is latent function ###

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

logL=function(H)
{
  alpha = 0.99
  
  r_cuu <- cuu(z, z, H); r_cfu <- cfu(x, z, H); r_cff <- cff(x, x, H)
  fuf <- r_cfu%*%mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)
  B <- alpha * fuf + diag(H[3]^2, dim(fuf)[1]) + (1-alpha) * r_cff
  
  RHS <- cbind(y, r_cfu)
  
  B_inv_y <- mBCG_noT(B, RHS, maxiter = MAT_ITER)
  
  A <- diag(H[3]^2, dim(fuf)[1]) + (1-alpha) * r_cff + alpha * fuf
  ch <- chol(A + diag(0, ncol(A) ), pivot=TRUE, tol=1e-16)
  logdeter_A <- 2*(sum(log(diag(ch))))
  
  # logdeter <- log(det(  diag(1, ncol(r_cuu)) + H[6]* mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%RHS[,2:(ncol(RHS))]  )) + logdeter_A
  logdeter <- 0 + logdeter_A
  
  B2 <-  diag(1, ncol(r_cuu)) + (1-2*alpha) * mBCG_noT(r_cuu, t(r_cfu), maxiter = MAT_ITER)%*%B_inv_y[,2:(ncol(RHS))]
  
  logdeter2 <- log(det(B2))
  
  # likelihood
  a <- -0.5*log(2*pi)*length(y) + 
    (-1/(2-2*alpha))*logdeter - 
    0.5*t(y)%*%B_inv_y[,1] + 
    (-alpha/(2*(1-alpha) ) )*( log(1/H[3]^2)*length(y) + logdeter2) 
  
  
  return(as.numeric(-a)) # neg-lilelihood
  
}

length_1 = 45
length_2 = 45
length_3 = 10
H_1 <- seq(-1, 15, length.out = length_1)
H_2 <- seq(0.01, 3, length.out = length_2)
H_3 <- seq(0.6, 3, length.out = length_3)
H_input <- expand.grid(H_1, H_2, H_3)
H_input <- matrix(as.matrix(H_input), length_1*length_2*length_3, 3)
H_outer <- array(apply(H_input, 1, logL), dim=c(length_1, length_2, length_3))
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
filled.contour(x = H_1, y = H_2, z = -(H_outer[,,2]), col=rev(cols)) # this plot log likelihood !

# , col = colorpanel(30, "white", "grey10")
