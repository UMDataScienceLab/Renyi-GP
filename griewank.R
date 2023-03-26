## The griewank function from https://www.sfu.ca/~ssurjano/griewank.html

griewank <- function(xx)
{
  
  ii <- c(1:length(xx))
  sum <- sum(xx^2/4000)
  prod <- prod(cos(xx/sqrt(ii)))
  
  y <- sum - prod + 1
  return(y)
}
