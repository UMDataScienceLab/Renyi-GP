#~~~~~~~~~~ mBCG algorithm without T matrices ~~~~~~~~~~#

mBCG_noT <- function(A, B, P, maxiter = 10000, tol = 1e-05){
  
  # mBCG algorithm without returning tridiagonal matrices
  
  # A - covariance matrix (want to calculate inverse)
  # B - matrix 
  # P - preconditioning matrix, default is diagonal scaling matrix
  # maxiter - should be large enough to guarantee convergence
  
  # return - A^{-1}B
  
  if (missing(P)) {
    dA <- diag(A)
    dA[which(dA == 0)] = 1e-04
    Pinv = diag(1/dA, nrow = nrow(A))
  }
  
  
  U0 = matrix(0, nrow(B), ncol(B)) 
  R = B - A %*% U0
  Z = Pinv %*% (R)
  D = Z
  
  iter = 0
  sumr2 = norm(R)
  
  n = dim(A)[1]
  mm = dim(A)[2]
  one = matrix(1, n, 1)
  
  while (sumr2 > (tol) & iter <= (maxiter) ) {
    
    iter = iter + 1
    Ap = A %*% D
    a = as.numeric( (t(R*Z)%*%one)/(t(D*Ap)%*%one) )
    U0 = U0 + t(diag(a) %*% t(D))
    R1 = R - t(diag(a) %*% t(Ap))
    Z1 = Pinv %*% R1
    bet = as.numeric( (t(Z1*R1)%*%one)/(t(Z*R)%*%one) )
    D = Z1 + t(diag(bet) %*% t(D))
    Z = Z1
    R = R1
    sumr2 = norm(R)
    
  }
  
  return(U0)
  
}
