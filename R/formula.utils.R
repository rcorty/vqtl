#' @title replace.markers.with.add.dom_
#' @rdname formula.utils
#'
#' @inheritParams scanonevar
#'
#' @return
#' @export
#'
#' @examples
replace.markers.with.add.dom_ <- function(cross,
                                          mean.formula,
                                          var.formula) {

  marker.names <- colnames(qtl::pull.geno(cross = cross))

  mean.covar.names <- labels(terms(mean.formula))
  var.covar.names <- labels(terms(var.formula))

  mean.marker.covars <- mean.covar.names[mean.covar.names %in% marker.names]
  var.marker.covars <- var.covar.names[var.covar.names %in% marker.names]

  for (mean.marker.covar in mean.marker.covars) {
    new.terms <- paste0('(',
                        paste0(mean.marker.covar,
                               c('_add', '_dom'),
                               collapse = '+'),
                        ')')

    mean.formula <- reformulate(termlabels = gsub(pattern = mean.marker.covar,
                                                  replacement = new.terms,
                                                  x = labels(terms(mean.formula))),
                                response = mean.formula[[2]])
  }

  for (var.marker.covar in var.marker.covars) {
    new.terms <- paste0('(',
                        paste0(var.marker.covar,
                               c('_add', '_dom'),
                               collapse = '+'),
                        ')')
    var.formula <- reformulate(termlabels = gsub(pattern = var.marker.covar,
                                                 replacement = new.terms,
                                                 x = labels(terms(var.formula))),
                               response = var.formula[[2]])
  }

  return(list(mean.formula = mean.formula,
              var.formula = var.formula))
}



is.mean.formula <- function(x) {

  # mean.formula must have a LHS, operator, and RHS
  if (length(x) != 3)
    return(FALSE)

  # first element must be a squiggle
  if (x[[1]] != '~')
    return(FALSE)

  # LHS must have exactly one variable
  if (length(all.vars(x[[2]])) != 1)
    return(FALSE)

  return(TRUE)
}


is.var.formula <- function(x) {

  # var.formula must have an operator and a RHS
  if (length(x) != 2)
    return(FALSE)

  if (x[[1]] == '~')
    return(FALSE)

  return(TRUE)
}

















