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
#' @param random.seed value to start the random number generator at, for reproducibility
#' @param n.cores number of cores to use for the permutations
#' @param silent Should all messaging be suppressed?
#'
#' @return 27599
#' @export
#'
#' @importFrom foreach %dopar%
#'
#' @examples
#' set.seed(27599)
#' test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), n.mar = 5), n.ind = 50)
#' scanonevar(cross = test.cross)
#'
scanonevar.perm <- function(sov,
                            n.perms,
                            random.seed = 27599,
                            n.cores = 1,
                            silent = TRUE) {

  stopifnot(is.scanonevar(sov))

  # no way to split perms over cores...so never use more cores than there are perms
  n.cores <- min(n.cores, n.perms)

  # get inputs into a format that is easy for scanonevar_ to use
  wrangled.inputs <- wrangle.scanonevar.input_(cross = sov[['meta']][['cross']],
                                               mean.formula = sov[['meta']][['formulae']][['mean.alt.formula']],
                                               var.formula = sov[['meta']][['formulae']][['var.alt.formula']],
                                               chrs = sov[['meta']][['chrs']])

  # execute the scan
  result <- scanonevar.perm_(sov = sov[['result']],
                             modeling.df = wrangled.inputs$modeling.df,
                             genoprob.df = wrangled.inputs$genoprob.df,
                             loc.info.df = wrangled.inputs$loc.info.df,
                             scan.types = wrangled.inputs$scan.types,
                             scan.formulae = wrangled.inputs$scan.formulae,
                             model = wrangled.inputs$model,
                             cross_type = class(sov[['meta']][['cross']])[1],
                             n.perms = n.perms,
                             seed = random.seed,
                             n.cores = n.cores,
                             silent = silent)

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
                             model,
                             cross_type,
                             n.perms,
                             seed,
                             n.cores,
                             silent) {

  this.context.permutation.max.finder <- function(alt.fitter, null.fitter) {
    permutation.max.finder(alt.fitter = alt.fitter,
                           null.fitter = null.fitter,
                           modeling.df = modeling.df,
                           loc.info.df = loc.info.df,
                           genoprob.df = genoprob.df,
                           scan.formulae = scan.formulae,
                           model = model,
                           cross_type = cross_type,
                           n.perms = n.perms,
                           seed = seed,
                           n.cores = n.cores)
  }

  if (!is.numeric(n.cores)) {
    stop('n.cores must be numeric')
  }
  if (n.cores != 1) {
    cl <- parallel::makeCluster(spec = n.cores)
    doParallel::registerDoParallel(cl = cl)
  }

  perms <- list()

  if ('mean' %in% scan.types) {
    if (!silent) { message('Starting mean permutations...') }

    mean.lod.maxes <- this.context.permutation.max.finder(
      alt.fitter = function(formulae, df, the.perm) {
        fit_model(formulae = formulae,
                  data = df,
                  model = model,
                  mean = 'alt',
                  var = 'alt',
                  permute_what = 'mean',
                  the.perm = the.perm)
      },
      null.fitter = function(formulae, df, the.perm) {
        fit_model(formulae = formulae,
                  data = df,
                  model = model,
                  mean = 'null',
                  var = 'alt')
      })

    if (!silent) { message('Finished mean permutations...') }
    perms[['mean']] <- dplyr::bind_cols(list(test = rep('mean', nrow(mean.lod.maxes))),
                                        mean.lod.maxes)
  }

  if ('var' %in% scan.types) {
    if (!silent) {  message('Starting variance permutations...') }

    var.lod.maxes <- this.context.permutation.max.finder(
      alt.fitter = function(formulae, df, the.perm) {
        fit_model(formulae = formulae,
                  data = df,
                  model = model,
                  mean = 'alt',
                  var = 'alt',
                  permute_what = 'var',
                  the.perm = the.perm)
      },
      null.fitter = function(formulae, df, the.perm) {
        fit_model(formulae = formulae,
                  data = df,
                  model = model,
                  mean = 'alt',
                  var = 'null')
      })

    if (!silent) { message('Finished variance permutations...') }
    perms[['var']] <- dplyr::bind_cols(list(test = rep('var', nrow(var.lod.maxes))),
                                       var.lod.maxes)
  }

  if ('joint' %in% scan.types) {
    if (!silent) { message('Starting joint mean-variance permutations...') }
    joint.lod.maxes <- var.lod.maxes <- this.context.permutation.max.finder(
      alt.fitter = function(formulae, df, the.perm) {
        fit_model(formulae = formulae,
                  data = df,
                  model = model,
                  mean = 'alt',
                  var = 'alt',
                  permute_what = 'both',
                  the.perm = the.perm)
      },
      null.fitter = function(formulae, df, the.perm) {
        fit_model(formulae = formulae,
                  data = df,
                  model = model,
                  mean = 'null',
                  var = 'null')
      })

    if (!silent) {  message('Finished joint mean-variance permutations...') }
    perms[['joint']] <- dplyr::bind_cols(list(test = rep('joint', nrow(joint.lod.maxes))),
                                         joint.lod.maxes)
  }

  if (n.cores != 1) {
    parallel::stopCluster(cl)
  }

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
                                   model,
                                   cross_type,
                                   n.perms,
                                   seed,
                                   n.cores) {

  loc.name <- null.ll <- alt.ll <- chr.type <- LOD.score <- 'fake_global_for_CRAN'

  result <- initialize.scanonevar.result_(loc.info.df = loc.info.df,
                                          scan.types = NA,
                                          scan.formulae = scan.formulae)
  result[['null.ll']] <- result[['alt.ll']] <- NA


  # fit the null model across the genome
  for (loc.idx in 1:nrow(result)) {

    # fill modeling.df with the genoprobs at the focal loc
    this.loc.name <- result[['loc.name']][loc.idx]
    loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                   loc.name == this.loc.name)

    this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                          loc.genoprobs = loc.genoprobs,
                                                          model.formulae = scan.formulae,
                                                          cross_type = cross_type)

    null.fit <- null.fitter(formulae = scan.formulae, df = this.loc.modeling.df)

    result[['null.ll']][loc.idx] <- tryNA(log_lik(null.fit))
  }

  # n.perms times, fit the alternative model with the focal locus permuted
  # subtract null ll to compute LOD score at each locus
  # take the max of the LOD scores

  # single core version
  if (n.cores == 1) {

    max.lods <- list()
    pb <- utils::txtProgressBar(min = 0, max = n.perms, style = 3)

    for (perm.idx in 1:n.perms) {

      utils::setTxtProgressBar(pb = pb, value = perm.idx)
      set.seed(seed = seed + perm.idx)
      the.perm <- sample(x = nrow(this.loc.modeling.df))

      for (loc.idx in 1:nrow(result)) {

        # fill modeling.df with the genoprobs at the focal loc
        this.loc.name <- result[['loc.name']][loc.idx]
        loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                       loc.name == this.loc.name)

        this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                              loc.genoprobs = loc.genoprobs,
                                                              model.formulae = scan.formulae,
                                                              cross_type = cross_type)

        alt.fit <- alt.fitter(formulae = scan.formulae, df = this.loc.modeling.df, the.perm = the.perm)

        result[['alt.ll']][loc.idx] <- tryNA(log_lik(alt.fit))
      }

      maxes <- result %>%
        dplyr::mutate(LOD.score = LOD_from_LLs(null_ll = null.ll, alt_ll = alt.ll)) %>%
        dplyr::group_by(chr.type) %>%
        dplyr::summarise(max.lod = max(LOD.score, na.rm = TRUE))

      max.lods[[perm.idx]] <- maxes
    }
  }

  # multi-core version
  if (n.cores != 1) {

    max.lods <- foreach::foreach(perm.idx = 1:n.perms) %dopar% {

      set.seed(seed = seed + perm.idx)
      the.perm <- sample(x = nrow(this.loc.modeling.df))

      for (loc.idx in 1:nrow(result)) {

        # fill modeling.df with the genoprobs at the focal loc
        this.loc.name <- result[['loc.name']][loc.idx]
        loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                       loc.name == this.loc.name)

        this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                              loc.genoprobs = loc.genoprobs,
                                                              model.formulae = scan.formulae,
                                                              cross_type = cross_type)

        alt.fit <- alt.fitter(formulae = scan.formulae, df = this.loc.modeling.df, the.perm = the.perm)

        result[['alt.ll']][loc.idx] <- tryNA(log_lik(alt.fit))
      }

      maxes <- result %>%
        dplyr::mutate(LOD.score = LOD_from_LLs(null_ll = null.ll, alt_ll = alt.ll)) %>%
        dplyr::group_by(chr.type) %>%
        dplyr::summarise(max.lod = max(LOD.score, na.rm = TRUE))

      maxes
    }
  }

  return(dplyr::bind_rows(max.lods))
}



calc.empir.ps <- function(sov, perms) {

  test <- chr.type <- fitted <- max.lod <- 'fake_global_for_CRAN'

  tests <- unique(perms[['test']])
  chr.types <- unique(perms[['chr.type']])

  for (this.test in tests) {

    for (this.chr.type in chr.types) {
      the.evd <- perms %>% dplyr::filter(test == this.test, chr.type == this.chr.type) %>% dplyr::pull(max.lod) %>% stats::na.exclude() %>% evd::fgev(std.err = FALSE)
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

