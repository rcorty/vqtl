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



is.var.formula <- function(x) {

  # var.formula must have an operator and a RHS
  if (length(x) != 2)
    return(FALSE)

  if (x[[1]] != '~')
    return(FALSE)

  return(TRUE)
}




is.formulae <- function(x) {

  if (!(all(c('mean.formula', 'var.formula') %in% names(x))))
    return(FALSE)

  if (any(!is.mean.formula(x[['mean.formula']]), !is.var.formula(x[['var.formula']])))
    return(FALSE)

  return(TRUE)
}



make.formulae_ <- function(mean.formula, var.formula) {
  stopifnot(is.mean.formula(mean.formula), is.var.formula(var.formula))
  return(list(mean.formula = mean.formula,
              var.formula = var.formula))
}



replace.markers.with.add.dom_ <- function(cross,
                                          mean.formula,
                                          var.formula) {

  marker.names <- colnames(qtl::pull.geno(cross = cross))

  mean.covar.names <- labels(stats::terms(mean.formula))
  var.covar.names <- labels(stats::terms(var.formula))

  mean.marker.covars <- mean.covar.names[mean.covar.names %in% marker.names]
  var.marker.covars <- var.covar.names[var.covar.names %in% marker.names]

  for (mean.marker.covar in mean.marker.covars) {
    new.terms <- paste0('(', paste0(mean.marker.covar,
                                    c('_add', '_dom'),
                                    collapse = '+'), ')')

    mean.formula <- stats::reformulate(termlabels = gsub(pattern = mean.marker.covar,
                                                         replacement = new.terms,
                                                         x = labels(stats::terms(mean.formula))),
                                       response = mean.formula[[2]])
  }

  for (var.marker.covar in var.marker.covars) {
    new.terms <- paste0('(', paste0(var.marker.covar,
                                    c('_add', '_dom'),
                                    collapse = '+'), ')')
    var.formula <- stats::reformulate(termlabels = gsub(pattern = var.marker.covar,
                                                        replacement = new.terms,
                                                        x = labels(stats::terms(var.formula))))
  }

  return(list(mean.formula = mean.formula,
              var.formula = var.formula))
}



remove.qtl.terms_ <- function(formulae) {

  stopifnot(is.formulae(formulae))

  mean.formula <- formulae[['mean.formula']]
  var.formula <- formulae[['var.formula']]

  mean.covar.name <- all.vars(mean.formula[[3]])
  var.covar.names <- all.vars(var.formula)

  # identify terms that are QTL keywords
  mean.qtl.idxs <- grep(pattern = 'mean.QTL', x = mean.covar.name)
  var.qtl.idxs <- grep(pattern = 'var.QTL', x = var.covar.names)

  # if no qtl terms, 'mean.null' is NULL and no mean testing will be done
  # if no non-qtl terms, rhs is just 1
  if (length(mean.qtl.idxs) == 0) {
    mean.null.formula <- mean.formula
  } else if (length(mean.qtl.idxs) == length(labels(mean.covar.name))) {
    mean.null.formula <- stats::reformulate(termlabels = '1', response = mean.formula[[2]])
  } else {
    mean.null.formula <- stats::formula(stats::drop.terms(termobj = stats::terms(mean.formula),
                                                          dropx = mean.qtl.idxs,
                                                          keep.response = TRUE))
  }

  # if no qtl terms, 'var.null' is NULL and no var testing will be done
  # if no non-qtl terms, rhs is just 1
  if (length(var.qtl.idxs) == 0) {
    var.null.formula <- var.formula
  } else if (length(var.qtl.idxs) == length(labels(var.covar.names))) {
    var.null.formula <- stats::reformulate(termlabels = '1', response = NULL)
  } else {
    var.null.formula <- stats::formula(stats::drop.terms(termobj = stats::terms(var.formula),
                                                         dropx = var.qtl.idxs,
                                                         keep.response = FALSE))
  }

  return(list(mean.formula = mean.null.formula,
              var.formula = var.null.formula))
}


has_a_random_term <- function(f) {
  any(grepl(pattern = '\\|', x = labels(terms(f))))
}

