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

  mf <- swtich(mean, alt = formulae[['mean.alt.formula']], null = formulae[['mean.null.formula']])
  vf <- swtich(var, alt = formulae[['var.alt.formula']], null = formulae[['var.null.formula']])

  do.call(what = purrr::compose(tryNA, fitter),
          args = list(mf = mf, df =  vf, data = data))
}

#
# # the alternative model
# fit.model.mv_ <- function(formulae, df, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.alt.formula']],
#                                      df =  formulae[['var.alt.formula']],
#                                      data = df))
# }
#
#
# # the three null models
# fit.model.m0_ <- function(formulae, df, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.alt.formula']],
#                                      df = formulae[['var.null.formula']],
#                                      data = df))
# }
#
# fit.model.0v_ <- function(formulae, df, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.null.formula']],
#                                      df = formulae[['var.alt.formula']],
#                                      data = df))
# }
#
# fit.model.00_ <- function(formulae, df, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.null.formula']],
#                                      df = formulae[['var.null.formula']],
#                                      data = df))
# }
#
#
# # the three permuted alternative models
# fit.model.m.star.v.star_ <- function(formulae, df, the.perm, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.alt.formula']],
#                                      df = formulae[['var.alt.formula']],
#                                      data = permute.QTL.terms_(df = df,
#                                                                the.perm = the.perm)))
#
# }
#
# fit.model.m.star.v_ <- function(formulae, df, the.perm, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.alt.formula']],
#                                      df = formulae[['var.alt.formula']],
#                                      data = permute.mean.QTL.terms_(df = df,
#                                                                     the.perm = the.perm)))
#
# }
#
# fit.model.m.v.star_ <- function(formulae, df, the.perm, model = c('dglm', 'hglm')) {
#
#   model <- match.arg(arg = model)
#
#   fitter <- switch(EXPR = model,
#                    dglm = dglm.ml,
#                    hglm = hglm_w_ll)
#
#   do.call(what = fitter, args = list(mf = formulae[['mean.alt.formula']],
#                                      df = formulae[['var.alt.formula']],
#                                      data = permute.var.QTL.terms_(df = df,
#                                                                    the.perm = the.perm)))
#
# }
#
#
#
#
#
#
#
#
