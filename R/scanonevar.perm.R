#' @title scanonevar
#' @name scanonevar
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @description \code{scanonevar} conducts a genome scan in an experimental
#' cross, accommodating covariate effects in residual variance and identifying
#' genetic effects on residual variance.
#'
#' @param cross The \code{cross}, built by \pkg{qtl} to be used in mapping
#' @param mean.formula The formula to describe the mean of the phenotype.
#' Keywords are mean.QTL.add and mean.QTL.dom for the additive and dominance
#' components of the QTL effect on the mean.
#' @param var.formula The formula to describe the residual variance of the
#' phenotype.  Keywords are var.QTL.add and var.QTL.dom for the additive and
#' dominance components of the QTL effect on residual phenotype variance.
#'
#' @return 27599
#' @export
#'
#' @examples
#' set.seed(27599)
#' test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), n.mar = 5), n.ind = 50)
#' scanonevar(cross = test.cross)
#'
scanonevar.perm <- function(sov,
                            n.perms) {

  stopifnot(is.scanonevar(sov))

  # get inputs into a format that is easy for scanonevar_ to use
  wrangled.inputs <- wrangle.scanonevar.input_(cross = sov[['meta']][['cross']],
                                               mean.formula = sov[['meta']][['formulae']][['mean.alt.formula']],
                                               var.formula = sov[['meta']][['formulae']][['var.alt.formula']],
                                               chrs = sov[['meta']][['chrs']])

  # execute the scan
  result <- scanonevar.perm_(modeling.df = wrangled.inputs$modeling.df,
                             loc.info.df = wrangled.inputs$loc.info.df,
                             genoprob.df = wrangled.inputs$genoprob.df,
                             scan.types = wrangled.inputs$scan.types,
                             scan.formulae = wrangled.inputs$scan.formulae,
                             n.perms = n.perms)

  sov <- list(meta = sov[['meta']],
              result = result)

  class(sov) <- c('scanonevar.perm', class(sov))
  return(sov)
}



scanonevar.perm_ <- function(modeling.df,
                             loc.info.df,
                             genoprob.df,
                             scan.types,
                             scan.formulae,
                             n.perms) {

  if ('mean' %in% scan.types) {
    scanonevar.meanperm_()
  }

  if ('var' %in% scan.types) {
    scanonevar.varperm()
  }

  if ('joint' %in% scan.types) {
    scanonevar.jointperm()
  }

  return(result)
}


scanonevar.meanperm_ <- function(modeling.df,
                                 loc.info.df,
                                 genoprob.df,
                                 scan.types,
                                 scan.formulae,
                                 n.perms) {

  result <- initialize.scanonevar.result_(loc.info.df = loc.info.df,
                                          scan.types = 'mean',
                                          scan.formulae = scan.formulae)

  df <- sum(grepl(pattern = 'mean.QTL', x = labels(terms(scan.formulae[['mean.alt.formula']]))))

  # fit the mean null model across the genome
  for (loc.idx in 1:nrow(result)) {

    # fill modeling.df with the genoprobs at the focal loc
    this.loc.name <- result[['loc.name']][loc.idx]
    loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                   loc.name == this.loc.name)

    this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                          loc.genoprobs = loc.genoprobs,
                                                          model.formulae = scan.formulae)

    # hacky way to accomodate x chromosome...
    # I'm OK with having 0 for the dominance components.  That feels very R-ish to me
    # and dglm handles it naturally, ignoring them in the model fitting and
    # giving an effect estimate of NA.  I don't love adjusting the df this way...
    if (any(this.loc.modeling.df[['mean.QTL.dom']] != 0)) {
      this.loc.mean.df <- mean.df
    } else {
      this.loc.mean.df <- mean.df - 1
    }
    if (any(this.loc.modeling.df[['var.QTL.dom']] != 0)) {
      this.loc.var.df <- var.df
    } else {
      this.loc.var.df <- var.df - 1
    }

    mean.null.fit <- tryCatch(expr = dglm::dglm(formula = scan.formulae[['mean.null.formula']],
                                                dformula = scan.formulae[['var.alt.formula']],
                                                data = this.loc.modeling.df),
                              warning = function(w) NA,
                              error = function(e) NA,
                              finally = NA)
    result[['mean.lod']][loc.idx] <- LOD(alt = alternative.fit, null = mean.null.fit)
    result[['mean.asymp.p']][loc.idx] <- pchisq(q = result[['mean.lod']][loc.idx], df = this.loc.mean.df, lower.tail = FALSE)

  }

  # n. perms times, fit an alternative model with mean permuted


}