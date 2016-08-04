#' @title validate.scanonevar.input_
#' @name validate.scanonevar.input_
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams scanonevar
#'
#' @return TRUE is parameters are valid for scanonevar, FALSE otherwise
#'
#' @export
#' @examples
#' x <- 27599
validate.scanonevar.input_ <- function(cross,
                                       mean.formula,
                                       var.formula) {

  stopifnot(is.cross(cross))

  # would like to have a less fragile way to get the LHS of mean.formula
  # formula.tools::lhs(mean.formula) gave an error
  stopifnot(as.character(mean.formula[[2]]) %in% names(cross[['pheno']]))

  phen.names <- names(cross[['pheno']])
  marker.names <- unlist(lapply(X = cross[['geno']],
                                FUN = function(chr) names(chr[['map']])))
  allowable.covar.names <- c(phen.names, marker.names)

  mean.covars <- labels(object = terms(x = mean.formula))
  non.keyword.mean.covars <- mean.covars[!(mean.covars %in% c('mean.QTL.add', 'mean.QTL.dom'))]

  var.covars <- labels(object = terms(x = var.formula))
  non.keywords.var.covars <- var.covars[!(var.covars %in% c('var.QTL.add', 'var.QTL.dom'))]

  stopifnot(!any(is.na(match(x = non.keyword.mean.covars,
                             table = allowable.covar.names))))
  stopifnot(!any(is.na(match(x = non.keywords.var.covars,
                             table = allowable.covar.names))))

  return(TRUE)
}