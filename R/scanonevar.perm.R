#' @title scanonevar.perm
#' @name scanonevar.perm
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @description \code{scanonevar.perm} conducts many permuted forms of
#' the \code{scanonevar} inputted, to assess the statistical significance of the results
#' in the inputted scanonevar in a FWER-controlling manner.
#'
#' @param sov the scanonevar whose significance should be assessed empirically in an FWER-controlling method
#' @param n.perms the number of permutations to do
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
              result = result[['sov']],
              perms = result[['perms']])

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

  perms <- list()

  if ('mean' %in% scan.types) {
    message('Starting mean permutations...')
    mean.lod.maxes <- this.context.permutation.max.finder(alt.fitter = fit.model.m.star.v_,
                                                          null.fitter = fit.model.0v_)
    perms[['mean']] <- dplyr::bind_cols(list(test = rep('mean', nrow(mean.lod.maxes))),
                                             mean.lod.maxes)
  }

  if ('var' %in% scan.types) {
    message('Starting variance permutations...')
    var.lod.maxes <- this.context.permutation.max.finder(alt.fitter = fit.model.m.v.star_,
                                                         null.fitter = fit.model.m0_)
    perms[['var']] <- dplyr::bind_cols(list(test = rep('var', nrow(var.lod.maxes))),
                                       var.lod.maxes)
  }

  if ('joint' %in% scan.types) {
    message('Starting joint mean-variance permutations...')
    joint.lod.maxes <- this.context.permutation.max.finder(alt.fitter = fit.model.m.star.v.star_,
                                                           null.fitter = fit.model.00_)
    perms[['joint']] <- dplyr::bind_cols(list(test = rep('joint', nrow(joint.lod.maxes))),
                                         joint.lod.maxes)
  }

  # calc empir ps here rather than in the mean, var, and joint 'ifs'
  perms <- dplyr::bind_rows(perms)
  sov <- calc.empir.ps(sov, perms)

  return(list(sov = sov,
              perms = perms))
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
  max.lods <- list()
  pb <- utils::txtProgressBar(min = 0, max = n.perms, style = 3)  # +1 here makes n.perms = 10 give, 10%, 20%, etc as progress updates
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

    maxes <- result %>%
      dplyr::mutate(LOD.score = 0.5*(null.ll - alt.ll)) %>%
      dplyr::group_by(chr.type) %>%
      dplyr::summarise(max.lod = max(LOD.score, na.rm = TRUE))

    max.lods[[perm.idx]] <- maxes
  }

  return(dplyr::bind_rows(max.lods))
}



calc.empir.ps <- function(sov, perms) {

  tests <- unique(perms[['test']])
  chr.types <- unique(perms[['chr.type']])

  for (this.test in tests) {

    for (this.chr.type in chr.types) {
      the.evd <- evd::fgev(x = perms %>% dplyr::filter(test == this.test, chr.type == this.chr.type) %>% .[['max.lod']])
      idxs <- sov[['chr.type']] == this.chr.type

      sov[[paste0(this.test, '.empir.p')]][idxs] <- evd::pgev(q = sov[[paste0(this.test, '.lod')]][idxs],
                                                              loc = fitted(the.evd)[1],
                                                              scale = fitted(the.evd)[2],
                                                              shape = fitted(the.evd)[3],
                                                              lower.tail = FALSE)
    }
  }

  return(sov)
}

