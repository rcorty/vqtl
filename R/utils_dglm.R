tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

dglm.ml <- function(mf, df, data) {
  dglm::dglm(formula = mf, dformula = df, data = data, method = 'ml')
}

# the alternative model
fit.model.mv_ <- function(formulae, df) {

  return(tryNA(dglm.ml(mf = formulae[['mean.alt.formula']],
                       df =  formulae[['var.alt.formula']],
                       data = df)))
}


# the three null models
fit.model.m0_ <- function(formulae, df) {

  return(tryNA(dglm.ml(mf = formulae[['mean.alt.formula']],
                       df = formulae[['var.null.formula']],
                       data = df)))
}

fit.model.0v_ <- function(formulae, df) {

  return(tryNA(dglm.ml(mf = formulae[['mean.null.formula']],
                       df = formulae[['var.alt.formula']],
                       data = df)))
}

fit.model.00_ <- function(formulae, df) {

  return(tryNA(dglm.ml(mf = formulae[['mean.null.formula']],
                       df = formulae[['var.null.formula']],
                       data = df)))
}


# the three permuted alternative models
fit.model.m.star.v.star_ <- function(formulae, df, the.perm) {

  return(tryNA(dglm.ml(mf = formulae[['mean.alt.formula']],
                       df = formulae[['var.alt.formula']],
                       data = permute.QTL.terms_(df = df,
                                                 the.perm = the.perm))))
}

fit.model.m.star.v_ <- function(formulae, df, the.perm) {

  return(tryNA(dglm.ml(mf = formulae[['mean.alt.formula']],
                       df = formulae[['var.alt.formula']],
                       data = permute.mean.QTL.terms_(df = df,
                                                      the.perm = the.perm))))
}

fit.model.m.v.star_ <- function(formulae, df, the.perm) {

  return(tryNA(dglm.ml(mf = formulae[['mean.alt.formula']],
                       df = formulae[['var.alt.formula']],
                       data = permute.var.QTL.terms_(df = df,
                                                     the.perm = the.perm))))
}



permute.mean.QTL.terms_ <- function(df, the.perm = sample(x = nrow(df))) {

  mean.qtl.col.names <- grep(pattern = 'mean.QTL', x = names(df), value = TRUE)
  for (mean.qtl.col.name in mean.qtl.col.names) {
    df[[mean.qtl.col.name]] <- df[[mean.qtl.col.name]][the.perm]
  }
  return(df)
}


permute.var.QTL.terms_ <- function(df, the.perm = sample(x = nrow(df))) {

  var.qtl.col.names <- grep(pattern = 'var.QTL', x = names(df), value = TRUE)
  for (var.qtl.col.name in var.qtl.col.names) {
    df[[var.qtl.col.name]] <- df[[var.qtl.col.name]][the.perm]
  }
  return(df)
}

permute.QTL.terms_ <- function(df, the.perm = sample(x = nrow(df))) {

  df2 <- permute.var.QTL.terms_(df = df, the.perm = the.perm)
  df3 <- permute.mean.QTL.terms_(df = df2, the.perm = the.perm)
  return(df3)
}




