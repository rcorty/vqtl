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
                      var = c('alt', 'null')) {

  model <- match.arg(arg = model)
  mean <- match.arg(arg = mean)
  var <- match.arg(arg = var)

  fitter <- switch(EXPR = model, dglm = dglm.ml, hglm = hglm_w_ll)

  mf <- switch(mean, alt = formulae[['mean.alt.formula']], null = formulae[['mean.null.formula']])
  vf <- switch(var, alt = formulae[['var.alt.formula']], null = formulae[['var.null.formula']])

  tryNA(do.call(what = fitter,
                args = list(mf = mf, df =  vf, data = data)))
}

