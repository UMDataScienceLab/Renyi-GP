## Gramacy & Lee (2012) Function

grlee12 <- function(x)
{
  
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4
  
  y <- term1 + term2
  return(y)
}
