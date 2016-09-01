tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}

# the alternative model
fit.model.mv_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.alt.formula']],
                          dformula = formulae[['var.alt.formula']],
                          data = df)))
}


# the three null models
fit.model.m0_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.alt.formula']],
                          dformula = formulae[['var.null.formula']],
                          data = df)))
}

fit.model.0v_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.null.formula']],
                          dformula = formulae[['var.alt.formula']],
                          data = df)))
}

fit.model.00_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.null.formula']],
                          dformula = formulae[['var.null.formula']],
                          data = df)))
}


# the three permuted alternative models
fit.model.m.star.v.star_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.alt.formula']],
                          dformula = formulae[['var.alt.formula']],
                          data = permute.QTL.terms_(df))))
}

fit.model.m.star.v_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.alt.formula']],
                          dformula = formulae[['var.alt.formula']],
                          data = permute.mean.QTL.terms_(df))))
}

fit.model.m.v.star_ <- function(formulae, df) {

  return(tryNA(dglm::dglm(formula = formulae[['mean.alt.formula']],
                          dformula = formulae[['var.alt.formula']],
                          data = permute.var.QTL.terms_(df))))
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

permute.QTL.terms_ <- function(df) {

  the.perm <- sample(x = nrow(df))
  df2 <- permute.var.QTL.terms_(df = df, the.perm = the.perm)
  df3 <- permute.mean.QTL.terms_(df = df2, the.perm = the.perm)
  return(df3)
}




