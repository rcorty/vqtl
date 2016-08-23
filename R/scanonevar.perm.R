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
  result <- scanonevar.perm_(sov = sov[['result']],
                             modeling.df = wrangled.inputs$modeling.df,
                             loc.info.df = wrangled.inputs$loc.info.df,
                             genoprob.df = wrangled.inputs$genoprob.df,
                             scan.types = wrangled.inputs$scan.types,
                             scan.formulae = wrangled.inputs$scan.formulae,
                             n.perms = n.perms)

  sov <- list(meta = sov[['meta']],
              result = result)

  class(sov) <- c('scanonevar', class(sov))
  return(sov)
}



scanonevar.perm_ <- function(sov,
                             modeling.df,
                             loc.info.df,
                             genoprob.df,
                             scan.types,
                             scan.formulae,
                             n.perms) {

  this.context.permutation.max.finder <- function(alt.fitter, null.fitter) {
    permutation.max.finder(alt.fitter = alt.fitter,
                           null.fitter = null.fitter,
                           modeling.df = modeling.df,
                           loc.info.df = loc.info.df,
                           genoprob.df = genoprob.df,
                           scan.formulae = scan.formulae,
                           n.perms = n.perms)
  }

  if ('mean' %in% scan.types) {
    message('Starting mean permutations...')
    mean.lod.maxes <- this.context.permutation.max.finder(alt.fitter = fit.model.m.star.v_,
                                                          null.fitter = fit.model.0v_)
    mean.evd <- evd::fgev(x = mean.lod.maxes)
    sov[['mean.empir.p']] <- pgev(q = sov[['mean.lod']],
                                  loc = fitted(mean.evd)[1],
                                  scale = fitted(mean.evd)[2],
                                  shape = fitted(mean.evd)[3],
                                  lower.tail = FALSE)
  }

  if ('var' %in% scan.types) {
    message('Starting variance permutations...')
    var.lod.maxes <- this.context.permutation.max.finder(alt.fitter = fit.model.m.v.star_,
                                                         null.fitter = fit.model.m0_)
    var.evd <- evd::fgev(x = var.lod.maxes)
    sov[['var.empir.p']] <- pgev(q = sov[['var.lod']],
                                 loc = fitted(var.evd)[1],
                                 scale = fitted(var.evd)[2],
                                 shape = fitted(var.evd)[3],
                                 lower.tail = FALSE)
  }

  if ('joint' %in% scan.types) {
    message('Starting joint mean-variance permutations...')
    joint.lod.maxes <- this.context.permutation.max.finder(alt.fitter = fit.model.m.star.v.star_,
                                                           null.fitter = fit.model.00_)
    joint.evd <- evd::fgev(x = joint.lod.maxes)
    sov[['joint.empir.p']] <- pgev(q = sov[['joint.lod']],
                                   loc = fitted(joint.evd)[1],
                                   scale = fitted(joint.evd)[2],
                                   shape = fitted(joint.evd)[3],
                                   lower.tail = FALSE)
  }

  return(sov)
}


permutation.max.finder <- function(alt.fitter,
                                   null.fitter,
                                   modeling.df,
                                   loc.info.df,
                                   genoprob.df,
                                   scan.formulae,
                                   n.perms) {

  result <- initialize.scanonevar.result_(loc.info.df = loc.info.df,
                                          scan.types = NA,
                                          scan.formulae = scan.formulae)
  result[['null.ll']] <- result[['alt.ll']] <- NA


  # fit the mean null model across the genome
  for (loc.idx in 1:nrow(result)) {

    # fill modeling.df with the genoprobs at the focal loc
    this.loc.name <- result[['loc.name']][loc.idx]
    loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                   loc.name == this.loc.name)

    this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                          loc.genoprobs = loc.genoprobs,
                                                          model.formulae = scan.formulae)

    null.fit <- null.fitter(formulae = scan.formulae, df = this.loc.modeling.df)

    if (!identical(NA, null.fit)) {
      result[['null.ll']][loc.idx] <- null.fit$m2loglik
    }
  }

  # n. perms times, fit an alternative model with mean permuted
  # subtract null ll to compute LOD score at each locus
  # take the max of the LOD scores
  max.lods <- rep(NA, n.perms)
  pb <- utils::txtProgressBar(min = 1, max = n.perms, style = 3)
  for (perm.idx in 1:n.perms) {

    utils::setTxtProgressBar(pb = pb, value = perm.idx)

    for (loc.idx in 1:nrow(result)) {

      # fill modeling.df with the genoprobs at the focal loc
      this.loc.name <- result[['loc.name']][loc.idx]
      loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                     loc.name == this.loc.name)

      this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                            loc.genoprobs = loc.genoprobs,
                                                            model.formulae = scan.formulae)

      alt.fit <- alt.fitter(formulae = scan.formulae, df = this.loc.modeling.df)

      if (!identical(NA, alt.fit)) {
        result[['alt.ll']][loc.idx] <- alt.fit$m2loglik
      }
    }

    max.lods[perm.idx] <- max(0.5*(result[['null.ll']] - result[['alt.ll']]), na.rm = TRUE)
  }

  return(max.lods)
}
