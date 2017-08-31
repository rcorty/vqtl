tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

dglm.ml <- function(mf, df, data) {
  dglm::dglm(formula = mf, dformula = df, data = data, method = 'ml')
}

hglm_w_ll <- function(mf, df, data) {
  hglm::hglm2(meanmodel = mf, disp = df, data = data, calc.like = TRUE)
}


fit_model <- function(formulae,
                      data,
                      model = c('dglm', 'hglm'),
                      mean = c('alt', 'null'),
                      var = c('alt', 'null'),
                      permute_what = c('none', 'mean', 'var', 'both'),
                      the.perm) {

  model <- match.arg(arg = model)
  mean <- match.arg(arg = mean)
  var <- match.arg(arg = var)
  permute_what <- match.arg(arg = permute_what)

  fitter <- switch(EXPR = model, dglm = dglm.ml, hglm = hglm_w_ll)

  mf <- switch(mean, alt = formulae[['mean.alt.formula']], null = formulae[['mean.null.formula']])
  vf <- switch(var, alt = formulae[['var.alt.formula']], null = formulae[['var.null.formula']])

  data <- switch(EXPR = permute_what,
                 none = data,
                 mean = permute.mean.QTL.terms_(df = data, the.perm = the.perm),
                 var = permute.var.QTL.terms_(df = data, the.perm = the.perm),
                 both = permute.QTL.terms_(df = data, the.perm = the.perm))

  tryNA(do.call(what = fitter,
                args = list(mf = mf, df =  vf, data = data)))
}



log_lik <- function(f) {

  if (inherits(x = f, what = 'dglm')) {
    return(-0.5*f$m2loglik)
  }
  if (inherits(x = f, what = 'hglm')) {
    return(f$likelihood$hlik)
  }
  return(stats::logLik(object = f))
}


LOD <- function(alt, null) {

  if (any(identical(alt, NA), identical(null, NA))) {
    return(NA)
  }

  if (!identical(class(alt), class(null))) {
    stop('Can only calculate LOD on models of the same class.')
  }

  if (!inherits(x = alt, what = c('dglm', 'hglm'))) {
    stop('Can only calcualte LOD on models of class dglm or hglm.')
  }

  if (inherits(x = alt, what = 'dglm')) {
    return(0.5*(null$m2loglik - alt$m2loglik )/log(10))
  }

  if (inherits(x = alt, what = 'hglm')) {
    return((alt$likelihood$hlik - null$likelihood$hlik)/log(10))
  }

}