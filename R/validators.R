formulae.are.valid.for.cross_ <- function(cross = cross,
                                          mean.formula = mean.formula,
                                          var.formula = var.formula) {

  # the response in 'mean.formula' must be a phenotype in the cross
  phen.names <- names(qtl::pull.pheno(cross = cross))
  stopifnot(all.vars(mean.formula[[2]]) %in% phen.names)
  # see http://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object/13218055#13218055
  # for a long discussion about how to pull response from a formula...


  # build up list of allowable covar names for mean and variance sub-models
  marker.names <- colnames(qtl::pull.geno(cross = cross))
  allowable.covar.names <- c(phen.names, marker.names, paste0(marker.names, '_add'), paste0(marker.names, '_dom'))
  allowable.mean.covar.names <- c(allowable.covar.names, 'mean.QTL.add', 'mean.QTL.dom')
  allowable.var.covar.names <- c(allowable.covar.names, 'var.QTL.add', 'var.QTL.dom')

  # extract covariate names from mean and variance sub-models
  mean.covars <- all.vars(mean.formula[[3]])
  var.covars <- all.vars(var.formula)

  # all covariates must be allowable
  stopifnot(all(mean.covars %in% allowable.mean.covar.names))
  stopifnot(all(var.covars %in% allowable.var.covar.names))
}





formulae.are.valid.for.scanonevar_ <- function(mean.formula,
                                               var.formula) {

  # make sure a 'QTL' term is used correctly somewhere
  mean.covars <- all.vars(mean.formula[[3]])
  var.covars <- all.vars(var.formula)
  stopifnot(any(c('mean.QTL.add', 'mean.QTL.dom') %in% mean.covars,
                c('var.QTL.add', 'var.QTL.dom') %in% var.covars))

}