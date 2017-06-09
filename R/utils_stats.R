LOD <- function(alt, null) {

  if(any(identical(x = alt, NA), identical(x = null, NA)))
    return(NA)

  stopifnot('dglm' %in% class(alt))
  stopifnot('dglm' %in% class(null))

  return(0.5*(null$m2loglik - alt$m2loglik )/log(10))
}