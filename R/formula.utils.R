#' @title is.mean.formula
#' @rdname internals
#'
#' @param x object to be tested whether it is a valid mean.formula
#'
#' @return TRUE if x is a valid mean.formula, FALSE otherwise
#'
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

  # consider checking whether RHS has any variables
  # i.e. should y ~ 1 be a valid mean.formula
  # I omitted this check because I think in some cases,
  # maybe related to plotting, this is a valid mean.formula
  # condition is: length(all.vars(x[[3]])) > 0

  return(TRUE)
}



#' @title is.var.formula
#' @rdname internals
#'
#' @param x object to be tested whether it is a valid var.formula
#'
#' @return TRUE if X is a valid var.formula, FALSE otherwise
#'
is.var.formula <- function(x) {

  # var.formula must have an operator and a RHS
  if (length(x) != 2)
    return(FALSE)

  if (x[[1]] != '~')
    return(FALSE)

  return(TRUE)
}




#' @title is.formulae
#' @rdname internals
#'
#' @param x object to be tested whether it is a valid formulae
#'
#' @return TRUE if x is a valid formulae, FALSE otherwise
#'
is.formulae <- function(x) {

  if (!identical(names(x), c('mean.formula', 'var.formula')))
    return(FALSE)

  if (any(!is.mean.formula(x[['mean.formula']]), !is.var.formula(x[['var.formula']])))
    return(FALSE)

  return(TRUE)
}



#' @title make.formulae
#' @rdname internals
#'
#' @param mean.formula a mean.formula
#' @param var.forula a var.formula
#'
#' @return a formulae, a list of length two where the first element is named
#' 'mean.formula' and is a mean.formula and the second element is named
#' 'var.formula' and is a var.formula.
#'
make.formulae <- function(mean.formula, var.formula) {
  stopifnot(is.mean.formula(mean.formula), is.var.formula(var.formula))
  return(list(mean.formula = mean.formula,
              var.formula = var.formula))
}



#' @title replace.markers.with.add.dom_
#' @rdname internals
#'
#' @inheritParams scanonevar
#'
#' @return A a list with two elements, mean.formula and var.formula,
#' where any markers from \code{cross} in the input mean.formula or var.formula
#' have been replaced with (marker.name_add + marker.name_dom).
#'
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
                                                 x = labels(terms(var.formula))))
  }

  return(list(mean.formula = mean.formula,
              var.formula = var.formula))
}





#' @title remove.qtl.terms_
#' @rdname internals
#'
#' @param formulae The formulae from which QTL terms will be removed
#'
#' @return a formulae object similar to the input formulae, but with QTL terms removed
#'
remove.qtl.terms_ <- function(formulae) {

  stopifnot(is.formulae(formulae))

  mean.formula <- formulae[['mean.formula']]
  var.formula <- formulae[['var.formula']]

  mean.covar.name <- all.vars(mean.formula[[3]])
  var.covar.names <- all.vars(formulae[['var.formula']])

  # identify terms that are QTL keywords
  mean.qtl.idxs <- grep(pattern = 'mean.QTL', x = mean.covar.name)
  var.qtl.idxs <- grep(pattern = 'var.QTL', x = var.covar.names)

  # if no qtl terms, 'mean.null' is NULL and no mean testing will be done
  # if no non-qtl terms, rhs is just 1
  if (length(mean.qtl.idxs) == 0) {
    mean.null.formula <- mean.formula
  } else if (length(mean.qtl.idxs) == length(labels(mean.covar.name))) {
    mean.null.formula <- reformulate(termlabels = '1', response = mean.formula[[2]])
  } else {
    mean.null.formula <- formula(drop.terms(termobj = terms(mean.formula),
                                            dropx = mean.qtl.idxs,
                                            keep.response = TRUE))
  }

  # if no qtl terms, 'var.null' is NULL and no var testing will be done
  # if no non-qtl terms, rhs is just 1
  if (length(var.qtl.idxs) == 0) {
    var.null.formula <- formulae[['var.formula']]
  } else if (length(var.qtl.idxs) == length(labels(var.covar.names))) {
    var.null.formula <- reformulate(termlabels = '1', response = NULL)
  } else {
    var.null.formula <- formula(drop.terms(termobj = terms(var.formula),
                                           dropx = var.qtl.idxs,
                                           keep.response = FALSE))
  }

  return(list(mean.formula = mean.null.formula,
              var.formula = var.null.formula))
}




