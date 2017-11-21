tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

fit_dglm <- function(mf, df, data, glm_family) {#, obs_weights) {
  dglm::dglm(formula = mf, dformula = df, data = data, method = 'ml', family = glm_family)#, weights = obs_weights)
}

fit_hglm <- function(mf, df, data, glm_family) {#, obs_weights) {
  hglm::hglm2(meanmodel = mf, disp = df, data = data, calc.like = TRUE, family = glm_family)#, weights = obs_weights)
}

fit_dhglm <- function(mf, df, data) {
  stop('dhglm not yet implemented.')
}

fit_model <- function(formulae,
                      data,
                      mean = c('alt', 'null'),
                      var = c('alt', 'null'),
                      model = c('dglm', 'hglm', 'dhglm'),
                      glm_family = c('gaussian', 'poisson'),
                      permute_what = c('none', 'mean', 'var', 'both'),
                      the.perm = seq(from = 1, to = nrow(data)),
                      obs_weights = rep(1, nrow(data))) {

  mean <- match.arg(arg = mean)
  var <- match.arg(arg = var)
  model <- match.arg(arg = model)
  glm_family <- match.arg(arg = glm_family)
  permute_what <- match.arg(arg = permute_what)

  mf <- switch(mean, alt = formulae[['mean.alt.formula']], null = formulae[['mean.null.formula']])
  vf <- switch(var, alt = formulae[['var.alt.formula']], null = formulae[['var.null.formula']])

  fit_model <- switch(EXPR = model,
                      dglm = fit_dglm,
                      hglm = fit_hglm,
                      dhglm = fit_dhglm)

  glm_family <- switch(EXPR = glm_family,
                       gaussian = stats::gaussian,
                       poisson = stats::poisson)

  data <- switch(EXPR = permute_what,
                 none = data,
                 mean = permute.mean.QTL.terms_(df = data, the.perm = the.perm),
                 var = permute.var.QTL.terms_(df = data, the.perm = the.perm),
                 both = permute.QTL.terms_(df = data, the.perm = the.perm))

  tryNA(fit_model(mf = mf, df =  vf, data = data, glm_family = glm_family))#, obs_weights = obs_weights)
  # tryNA(do.call(what = fit_model,
                # args = list(mf = mf, df =  vf, data = data, glm_family = glm_family, weights = obs_weights)))
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