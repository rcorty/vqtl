formulae_is_valid_ <- function(formulae) {

  if (!is.formulae(formulae)) {
    return(FALSE)
  }

  # must have at least one QTL term used appropriately
  mean.covars <- all.vars(formulae[['mean.formula']][[3]])
  var.covars <- all.vars(formulae[['var.formula']])
  if (all(!any(c('mean.QTL.add', 'mean.QTL.dom') %in% mean.covars),
          !any(c('var.QTL.add', 'var.QTL.dom') %in% var.covars))) {
    return(FALSE)
  }

  return(TRUE)
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

  mean.formula <- stats::update(old = mean.formula, new = ~ . -mean.QTL.add)
  mean.formula <- stats::update(old = mean.formula, new = ~ . -mean.QTL.dom)

  var.formula <- stats::update(old = var.formula, new = ~ . -var.QTL.add)
  var.formula <- stats::update(old = var.formula, new = ~ . -var.QTL.dom)

  return(list(mean.formula = mean.formula,
              var.formula = var.formula))
}


has_a_random_term <- function(f) {
  any(grepl(pattern = '\\|', x = labels(stats::terms(f))))
}

