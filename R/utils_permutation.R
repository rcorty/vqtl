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