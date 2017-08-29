# hglm_w_ll <- function(mf, df, data) {
#   hglm::hglm2(meanmodel = mf,
#               disp = df,
#               data = data,
#               calc.like = TRUE)
# }
#
# # the alternative model
# fit.model.mv_ <- function(formulae, df) {
#
#   return(tryNA(hglm_w_ll(mf = formulae[['mean.alt.formula']],
#                          df =  formulae[['var.alt.formula']],
#                          data = df)))
# }
#
#
# # the three null models
# fit.model.m0_ <- function(formulae, df) {
#
#   return(tryNA(hglm_w_ll(mf = formulae[['mean.alt.formula']],
#                          df = formulae[['var.null.formula']],
#                          data = df)))
# }
#
# fit.model.0v_ <- function(formulae, df) {
#
#   return(tryNA(hglm_w_ll(mf = formulae[['mean.null.formula']],
#                          df = formulae[['var.alt.formula']],
#                          data = df)))
# }
#
# fit.model.00_ <- function(formulae, df) {
#
#   return(tryNA(hglm_w_ll(mf = formulae[['mean.null.formula']],
#                          df = formulae[['var.null.formula']],
#                          data = df)))
# }
#
#
# # the three permuted alternative models
# fit.model.m.star.v.star_ <- function(formulae, df, the.perm) {
#
#   return(tryNA(hglm_w_ll(mf = formulae[['mean.alt.formula']],
#                          df = formulae[['var.alt.formula']],
#                          data = permute.QTL.terms_(df = df,
#                                                    the.perm = the.perm))))
# }
#
# fit.model.m.star.v_ <- function(formulae, df, the.perm) {
#
#   return(tryNA(hglm_w_ll(mf = formulae[['mean.alt.formula']],
#                          df = formulae[['var.alt.formula']],
#                          data = permute.mean.QTL.terms_(df = df,
#                                                         the.perm = the.perm))))
# }
#
# fit.model.m.v.star_ <- function(formulae, df, the.perm) {
#
#   return(tryNA(hglm_w_lls(mf = formulae[['mean.alt.formula']],
#                           df = formulae[['var.alt.formula']],
#                           data = permute.var.QTL.terms_(df = df,
#                                                         the.perm = the.perm))))
# }
#
#
