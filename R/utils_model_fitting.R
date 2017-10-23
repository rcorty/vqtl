tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

fit_dglm <- function(mf, df, data, family) {
  dglm::dglm(formula = mf, dformula = df, data = data, method = 'ml', family = family)
}

fit_hglm <- function(mf, df, data, family) {
  hglm::hglm2(meanmodel = mf, disp = df, data = data, calc.like = TRUE, family = family)
}

fit_dhglm <- function(mf, df, data) {

}

fit_model <- function(formulae,
                      data,
                      model = c('dglm', 'hglm', 'dhglm'),
                      mean = c('alt', 'null'),
                      var = c('alt', 'null'),
                      permute_what = c('none', 'mean', 'var', 'both'),
                      the.perm,
                      family) {

  model <- match.arg(arg = model)
  mean <- match.arg(arg = mean)
  var <- match.arg(arg = var)
  permute_what <- match.arg(arg = permute_what)

  fitter <- switch(EXPR = model,
                   dglm = fit_dglm,
                   hglm = fit_hglm,
                   dhglm = fit_dhglm)

  mf <- switch(mean, alt = formulae[['mean.alt.formula']], null = formulae[['mean.null.formula']])
  vf <- switch(var, alt = formulae[['var.alt.formula']], null = formulae[['var.null.formula']])

  data <- switch(EXPR = permute_what,
                 none = data,
                 mean = permute.mean.QTL.terms_(df = data, the.perm = the.perm),
                 var = permute.var.QTL.terms_(df = data, the.perm = the.perm),
                 both = permute.QTL.terms_(df = data, the.perm = the.perm))

  tryNA(do.call(what = fitter,
                args = list(mf = mf, df =  vf, data = data, family = family)))
}



log_lik <- function(f) {

  if (inherits(x = f, what = 'dglm')) {
    if (abs(f$m2loglik) > 1e8) { return(NA) }
    return(-0.5*f$m2loglik)
  }
  if (inherits(x = f, what = 'hglm')) {
    if (abs(f$likelihood$hlik) > 1e8) { return(NA) }
    return(f$likelihood$hlik)
  }
  return(stats::logLik(object = f))
}

LRT <- function(alt, null) {

  if (any(identical(alt, NA), identical(null, NA))) {
    return(NA)
  }

  if (!identical(class(alt), class(null))) {
    stop('Can only calculate LOD on models of the same class.')
  }

  if (!inherits(x = alt, what = c('dglm', 'hglm'))) {
    stop('Can only calcualte LOD on models of class dglm or hglm.')
  }

  LRT <- 2*(log_lik(alt) - log_lik(null))

}

LOD <- function(alt, null) {

  return(0.5*LRT(alt = alt, null = null)/log(10))

}

LOD_from_LLs <- function(null_ll, alt_ll) {
  return((alt_ll - null_ll)/log(10))
}

LRT_from_LLs <- function(null_ll, alt_ll) {
  return(2*(alt_ll - null_ll))
}

#' percent variance explained
#'
#' @param LOD the log odds between the null and alternative model
#' @param n the number of observations
#'
#' @return pve
#' @export
#'
pve <- function(LOD, n) {
  1 - 10^(-0.5*LOD/n)
}