## Gramacy & Lee (2012) Function from https://www.sfu.ca/~ssurjano/grlee12.html

grlee12 <- function(x)
{
  
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4
  
  y <- term1 + term2
  return(y)
}
