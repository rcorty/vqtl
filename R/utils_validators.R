formulae.is.valid.for.cross_ <- function(cross,
                                         formulae) {

  mean.formula <- formulae[['mean.formula']]
  var.formula <- formulae[['var.formula']]

  # the response in 'mean.formula' must be a phenotype in the cross
  phen.names <- names(cross$pheno)
  if (!(all.vars(mean.formula[[2]]) %in% phen.names)) {
    return(FALSE)
  }
  # see http://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object/13218055#13218055
  # for a long discussion about how to pull response from a formula...


  # build up list of allowable covar names for mean and variance sub-models
  marker.names <- colnames(qtl::pull.geno(cross = cross))
  allowable.covar.names <- c(phen.names, marker.names, paste0(marker.names, '_add'), paste0(marker.names, '_dom'))
  allowable.mean.covar.names <- c(allowable.covar.names, 'mean.QTL.add', 'mean.QTL.dom')
  allowable.var.covar.names <- c(allowable.covar.names, 'var.QTL.add', 'var.QTL.dom')

  # extract covariate names from mean and variance sub-models
  mean.covars <- all.vars(mean.formula[[3]])
  var.covars <- all.vars(var.formula)

  # all covariates must be allowable
  if (!all(mean.covars %in% allowable.mean.covar.names)) {
    return(FALSE)
  }

  if (!all(var.covars %in% allowable.var.covar.names)) {
    return(FALSE)
  }

  return(TRUE)
}






#' @title is.scanonevar
#' @rdname utils
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @description utilities for working with scanonevar objects
#'
#' @return TRUE is X is a scanonevar object, FALSE otherwise.
#' @export
#'
#' @examples
#' is.scanonevar(x = 3)
#'
#' test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 4), n.mar = 5))
#' test.cross <- qtl::calc.genoprob(cross = test.cross, step = 2)
#'
#' x <- scanonevar(cross = test.cross)
#' is.scanonevar(x)
#'
is.scanonevar <- function(x) {

  if (!('scanonevar' %in% class(x))) {
    return(FALSE)
  }

  if (!all(c('meta', 'result') %in% names(x))) {
    return(FALSE)
  }

  meta <- x[['meta']]
  result <- x[['result']]

  # valid meta
  # if (!any(identical(names(meta), c('cross', 'modeling.df', 'formulae', 'scan.types', 'model' ,'chrs')),
  #          identical(names(meta), c('cross', 'modeling.df', 'formulae', 'scan.types','chrs')))) {
  #   return(FALSE)
  # }

  if (!(is.cross(meta[['cross']]))) {
    return(FALSE)
  }

  # more possible


  # valid result
  if (!all(c('loc.name', 'chr', 'pos') %in% names(result)))
    return(FALSE)

  if (all(!all(c('mean.lod', 'mean.asymp.p') %in% names(result)),
          !all(c('var.lod', 'var.asymp.p') %in% names(result))))
    return(FALSE)

  # check that LOD scores are greater than 0

  # check that p-values are between 0 and 1

  return(TRUE)
}


#' @title is.scanonevar.w.perms
#' @rdname utils
#'
#' @param x object being tested
#'
#' @return TRUE if x is a scanone var with perms (typically,
#' outputted from scanonevar.perm), and FALSE otherwise.
#'
is.scanonevar.w.perms <- function(x) {

  if (!('perms' %in% names(x))) {
    return(FALSE)
  }

  if (!is.scanonevar(x)) {
    return(FALSE)
  }

  if (!any(grep(pattern = 'empir.p', names(x[['result']])))) {
    return(FALSE)
  }

  return(TRUE)
}

#' @title is.cross
#' @rdname utils
#' @author Robert W. Corty \email{rcorty@@gmail.com}
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
#' @rdname utils
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


#' @title is.bc.cross
#' @rdname utils
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams is.cross
#'
#' @return TRUE if x is a cross object of type 'bc' (backcross), FALSE otherwise
#' @export
#'
#' @examples
#' is.cross(3)
#' is.cross(qtl::sim.cross(map = qtl::sim.map()))
#'
is.f2.cross <- function(x) {

  if (!('bc' %in% class(x)))
    return(FALSE)

  return(is.cross(x))
}



#' @title is.cross.w.genoprobs
#' @rdname utils
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
    if (!(array.dim[3] %in% c(2, 3)))
      return(FALSE)
    return(TRUE)
  }
  if (!all(sapply(X = x[['geno']], FUN = function(chr) is.valid.prob.array(a = chr[['prob']], cross = x, chr = chr))))
    return(FALSE)

  # all values in the prob array must be between 0 and 1
  sapply(X = x[['geno']], FUN = function(x) all(x[['prob']] > 0 & x[['prob']] < 1))

  return(is.cross(x))
}
