BoxCoxTransformY <- function(y, lambda) {

  n <- length(y)
  geom.mean <- exp(sum(log(y))/n)

  if (lambda == 0) {

    return(geom.mean * log(y))

  } else {

    return((y ^ lambda - 1)/(lambda * geom.mean ^ (lambda - 1)))

  }
}