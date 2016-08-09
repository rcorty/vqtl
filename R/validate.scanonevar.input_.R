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
                                       var.formula,
                                       chrs = qtl::chrnames(cross)) {

  # each argument must be individually valid
  stopifnot(is.cross(cross))
  stopifnot(is.mean.formula(mean.formula))
  stopifnot(is.var.formula(var.formula))
  chrs <- match.arg(arg = chrs, several.ok = TRUE)
  allow.no.qtl <- match.arg(allow.no.qtl)

  # formulae must be valid for use with cross
  stopifnot(formulae.are.valid.for.cross_(cross = cross,
                                          mean.formula = mean.formula,
                                          var.formula = var.formula))

  # arguments must be valid for use in scanonevar
  stopifnot(formulae.are.valid.for.scanonevar_(mean.formula = mean.formula,
                                               var.formula = var.formula))


  return(TRUE)
}