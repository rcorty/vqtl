LOD <- function(alt, null) {

  if (any(identical(alt, NA), identical(null, NA))) {
    return(NA)
  }

  if (!identical(class(alt), class(null))) {
    stop('Can only calculate LOD on models of the same class.')
  }

  if (!inherits(x = alt, what = c('dglm', 'hglm'))) {
    stop('Can only calcualte LOD on models of class dglm or hglm.')
  }

  if (inherits(x = alt, what = 'dglm')) {
    return(0.5*(null$m2loglik - alt$m2loglik )/log(10))
  }

  if (inherits(x = alt, what = 'hglm')) {
    return((alt$likelihood$hlik - null$likelihood$hlik)/log(10))
  }

}