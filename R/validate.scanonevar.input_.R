#' @title validate.scanonevar.input_
#' @name validate.scanonevar.input_
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams scanonevar
#'
#' @return TRUE is parameters are valid for scanonevar, errors with a
#' helpful error message if parameters are invalid.
#'
#' @export
#' @examples
#' x <- 27599
#'
validate.scanonevar.input_ <- function(cross,
                                       mean.formula,
                                       var.formula) {

  # todo rewrite this function using if () { stop(message) }
  # rather than stopifnot to improve error messages

  # 'cross' must be a valid cross object
  stopifnot(is.cross(cross))


  # the response in 'mean.formula' must be a phenotype in the cross
  stopifnot(as.character(mean.formula[[2]]) %in% names(cross[['pheno']]))
  # would like to have a less fragile way to extract the LHS of mean.formula
  # formula.tools::lhs(mean.formula) gave an error I couldn't figure out


  # build up list of allowable covar names for mean and variance sub-models
  phen.names <- names(cross[['pheno']])
  marker.names <- unlist(lapply(X = cross[['geno']],
                                FUN = function(chr) names(chr[['map']])))
  allowable.covar.names <- c(phen.names, marker.names)
  allowable.mean.covar.names <- c(allowable.covar.names, 'mean.QTL.add', 'mean.QTL.dom')
  allowable.var.covar.names <- c(allowable.covar.names, 'var.QTL.add', 'var.QTL.dom')

  # extract covariate names from mean and variance sub-models
  mean.covars <- labels(object = terms(x = mean.formula))
  var.covars <- labels(object = terms(x = var.formula))
  # would like to have a less fragile way to extract this information from the formulae

  # all non-keyword covariates must be allowable
  stopifnot(all(mean.covars %in% allowable.mean.covar.names))
  stopifnot(all(var.covars %in% allowable.var.covar.names))


  # make sure at least one keyword appears in the right formula
  stopifnot(any(c('mean.QTL.add', 'mean.QTL.dom') %in% mean.covars,
                c('var.QTL.add', 'var.QTL.dom') %in% var.covars))


  return(TRUE)
}