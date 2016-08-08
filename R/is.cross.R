#' @title is.cross
#' @name is.cross
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @param x The object being tested for whether or not is is a cross.
#'
#' @return TRUE if x is a cross object, FALSE otherwise.
#' @export
#'
#' @examples
#' is.cross(3)
#' is.cross(qtl::sim.cross(map = qtl::sim.map()))
#'
is.cross <- function(x) {

  # class of x
  if(!('cross' %in% class(x)))
    return(FALSE)

  # pheno, an element of x
  if (!('pheno' %in% names(x)))
    return(FALSE)

  if (!('data.frame' %in% class(x[['pheno']])))
    return(FALSE)


  # geno, an element of x
  if (!('geno' %in% names(x)))
    return(FALSE)

  if (!('list' %in% class(x[['geno']])))
    return(FALSE)


  # chromosomes, elements of geno
  if (!all(sapply(X = x[['geno']], FUN = class) %in% c('A', 'X')))
    return(FALSE)

  if (!all(sapply(X = x[['geno']], FUN = function(x) ('data' %in% names(x)))))
      return(FALSE)


  # consistent number of individuals across pheno and all genos
  pheno.n <- nrow(x[['pheno']])
  geno.ns <- sapply(X = x[['geno']], FUN = function(x) nrow(x[['data']]))
  if (any(pheno.n != geno.ns))
    return(FALSE)

  # length of map matches number of genotypes in each chromosome
  if (!all(sapply(X = x[['geno']], FUN = function(x) ncol(x[['data']]) == length(x[['map']]))))
    return(FALSE)


  return(TRUE)
}



#' @title is.f2.cross
#' @name is.f2.cross
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams is.cross
#'
#' @return TRUE if x is a cross object of type F2, FALSE otherwise
#' @export
#'
#' @examples
#' is.cross(3)
#' is.cross(qtl::sim.cross(map = qtl::sim.map()))
#'
is.f2.cross <- function(x) {

  if (!('f2' %in% class(x)))
    return(FALSE)

  return(is.cross(x))
}




#' @title is.cross.w.genoprobs
#' @name is.cross.w.genoprobs
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams is.cross
#'
#' @return TRUE if x is a cross object with valid genoprobs for each chromosome,
#' FALSE otherwise
#'
#' @export
#'
#' @examples
#' a <- qtl::sim.cross(map = qtl::sim.map())
#' is.cross.w.genoprobs(x = a)
#' b <- qtl::calc.genoprob(cross = a)
#' is.cross.w.genoprobs(x = b)
#'
is.cross.w.genoprobs <- function(x) {

  # check that there is a 'prob' element in each chromosome
  if (!all(sapply(X = x[['geno']], FUN = function(x) 'prob' %in% names(x))))
    return(FALSE)

  # check that the dimension of each prob element matches the number of individuals
  is.valid.prob.array <- function(a, cross, chr) {
    array.dim <- dim(a)
    if (array.dim[1] != qtl::nind(cross))
      return(FALSE)
    if (array.dim[2] < ncol(chr[['data']]))
      return(FALSE)
    if (array.dim[3] != 3)
      return(FALSE)
    return(TRUE)
  }
  if (!all(sapply(X = x[['geno']], FUN = function(chr) is.valid.prob.array(a = chr[['prob']], cross = x, chr = chr))))
    return(FALSE)

  # all values in the prob array must be between 0 and 1
  sapply(X = x[['geno']], FUN = function(x) all(x[['prob']] > 0 & x[['prob']] < 1))

  return(is.cross(x))
}
