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
#' dglm model will be fit if mean.formula has only fixed effects.
#' hglm model will be fit if mean.formula has one or more random effects.
#' @param var.formula The formula to describe the residual variance of the
#' phenotype.  Keywords are var.QTL.add and var.QTL.dom for the additive and
#' dominance components of the QTL effect on residual phenotype variance.
#' var.formula must have only fixed effects.
#' @param chrs chromosomes to scan
#' @param return.covar.effects Should covariate effects estimated at each locus be returned?
#'
#' @return 27599
#' @export
#'
#' @examples
#' set.seed(27599)
#' test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), n.mar = 5), n.ind = 50)
#' scanonevar(cross = test.cross)
#'
scanonevar <- function(cross,
                       mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
                       var.formula = ~ var.QTL.add + var.QTL.dom,
                       chrs = qtl::chrnames(cross = cross),
                       return.covar.effects = FALSE) {

  # give an informative error message if input is invalid
  validate.scanonevar.input_(cross = cross,
                             mean.formula = mean.formula,
                             var.formula = var.formula,
                             chrs = chrs)

  # get inputs into a format that is easy for scanonevar_ to use
  wrangled.inputs <- wrangle.scanonevar.input_(cross = cross,
                                               mean.formula = mean.formula,
                                               var.formula = var.formula,
                                               chrs = chrs)

  # save meta-data on this scan
  meta <- pull.scanonevar.meta_(wrangled.inputs = wrangled.inputs,
                                chrs = chrs)

  # execute the scan
  result <- scanonevar_(modeling.df = wrangled.inputs$modeling.df,
                        loc.info.df = wrangled.inputs$loc.info.df,
                        genoprob.df = wrangled.inputs$genoprob.df,
                        scan.types = wrangled.inputs$scan.types,
                        scan.formulae = wrangled.inputs$scan.formulae,
                        model = wrangled.inputs$model,
                        return.covar.effects = return.covar.effects,
                        cross_type = class(meta$cross)[1])

  sov <- list(meta = meta,
              result = result)

  class(sov) <- c('scanonevar', class(sov))
  return(sov)
}




validate.scanonevar.input_ <- function(cross,
                                       mean.formula,
                                       var.formula,
                                       chrs) {

  # check validity of inputs
  stopifnot(is.cross(cross))
  stopifnot(all(chrs %in% qtl::chrnames(cross)))

  # make and then check formulae
  formulae <- make.formulae_(mean.formula, var.formula)
  stopifnot(formulae_is_valid_(formulae = formulae))

  # formulae must be valid for use with cross
  stopifnot(formulae.is.valid.for.cross_(cross = cross,
                                         formulae = formulae))

  return(TRUE)
}




wrangle.scanonevar.input_ <- function(cross,
                                      mean.formula,
                                      var.formula,
                                      chrs) {

  if (!is.cross.w.genoprobs(x = cross)) {
    message("calculating genoprobs with stepwidth = 5, off.end = 0, error.prob = 1e-4, map.function = 'haldane'")
    cross <- qtl::calc.genoprob(cross = cross, step = 5)
  }

  loc.info.df <- wrangle.loc.info.df_(cross = cross,
                                      chrs = chrs)

  genoprob.df.long <- wrangle.genoprob.df_(cross = cross)

  scan.types <- wrangle.scanonevar.types_(mean.formula, var.formula)

  model <- ifelse(test = has_a_random_term(mean.formula),
                  yes = 'hglm',
                  no = 'dglm')

  scan.formulae <- wrangle.scanonevar.formulae_(cross, mean.formula, var.formula)

  modeling.df <- wrangle.scanonevar.modeling.df_(cross = cross,
                                                 genoprobs = genoprob.df.long,
                                                 scan.formulae = scan.formulae)

  return(list(cross = cross,
              scan.types = scan.types,
              scan.formulae = scan.formulae,
              modeling.df = modeling.df,
              loc.info.df = loc.info.df,
              genoprob.df = genoprob.df.long,
              model = model))
}




wrangle.scanonevar.types_ <- function(mean.formula,
                                      var.formula) {

  mean.qtl.idxs <- grep(pattern = 'mean.QTL', x = labels(stats::terms(mean.formula)))
  var.qtl.idxs <- grep(pattern = 'var.QTL', x = labels(stats::terms(var.formula)))

  if (all(mean.qtl.idxs, !var.qtl.idxs)) {
    return('mean')
  }
  if (all(!mean.qtl.idxs, var.qtl.idxs)) {
    return('var')
  }
  if (all(mean.qtl.idxs, var.qtl.idxs)) {
    return(c('mean', 'var', 'joint'))
  }

  stop('Should never get here.')
}



wrangle.scanonevar.formulae_ <- function(cross,
                                         mean.formula,
                                         var.formula) {

  # first, replace terms that are simply marker names with marker.name_add + marker.name_dom
  alt.formulae <- replace.markers.with.add.dom_(cross,
                                                mean.formula,
                                                var.formula)

  null.formulae <- remove.qtl.terms_(formulae = alt.formulae)

  # this way of adding to a list doesn't add anything when RHS is NULL
  scan.formulae <- list()
  scan.formulae[['mean.alt.formula']] <- alt.formulae[['mean.formula']]
  scan.formulae[['mean.null.formula']] <-  null.formulae[['mean.formula']]
  scan.formulae[['var.alt.formula']] <- alt.formulae[['var.formula']]
  scan.formulae[['var.null.formula']] <-  null.formulae[['var.formula']]

  return(scan.formulae)
}




wrangle.scanonevar.modeling.df_ <- function(cross,
                                            scan.formulae,
                                            genoprobs) {

  response.model.df <- make.response.model.df_(cross = cross,
                                               formulae = scan.formulae)

  qtl.covar.model.df <- make.qtl.covar.model.df_(cross = cross,
                                                 formulae = scan.formulae)

  phen.covar.model.df <- make.phen.covar.model.df_(cross = cross,
                                                   formulae = scan.formulae)

  genet.covar.model.df <- make.genet.covar.add.dom.model.df_(cross = cross,
                                                             formulae = scan.formulae,
                                                             genoprobs = genoprobs)

  return(dplyr::bind_cols(response.model.df,
                          qtl.covar.model.df,
                          phen.covar.model.df,
                          genet.covar.model.df))
}






pull.scanonevar.meta_ <- function(wrangled.inputs,
                                  chrs) {

  meta <- list(cross = wrangled.inputs$cross,
               modeling.df = wrangled.inputs$modeling.df,
               formulae = wrangled.inputs$scan.formulae,
               scan.types = wrangled.inputs$scan.types,
               model = wrangled.inputs$model,
               chrs = chrs)

  class(meta) <- c('scanonevar.meta', class(meta))

  return(meta)
}



scanonevar_ <- function(modeling.df,
                        genoprob.df,
                        loc.info.df,
                        scan.types,
                        scan.formulae,
                        model,
                        return.covar.effects,
                        cross_type) {

  loc.name <- 'fake_global_for_CRAN'

  result <- initialize.scanonevar.result_(loc.info.df = loc.info.df,
                                          scan.types = scan.types,
                                          scan.formulae = scan.formulae,
                                          return.covar.effects = return.covar.effects)

  mean.df <- sum(grepl(pattern = 'mean.QTL', x = labels(stats::terms(scan.formulae[['mean.alt.formula']]))))
  var.df <- sum(grepl(pattern = 'var.QTL', x = labels(stats::terms(scan.formulae[['var.alt.formula']]))))

  # loop over loci
  for (loc.idx in 1:nrow(result)) {

    # fill modeling.df with the genoprobs at the focal loc
    this.loc.name <- result[['loc.name']][loc.idx]
    loc.genoprobs <- dplyr::filter(.data = genoprob.df,
                                   loc.name == this.loc.name)

    # puts a column of 0's for any dominance components if we are on X chromosome
    # dglm handles this naturally, ignoring it in model fitting and giving effect estimate of NA
    this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df,
                                                          loc.genoprobs = loc.genoprobs,
                                                          model.formulae = scan.formulae,
                                                          cross_type = cross_type)

    #  hacky way to adjust the df on the X chr
    if (all(this.loc.modeling.df[['mean.QTL.dom']] == 0)) {
      this.loc.mean.df <- mean.df - 1
    } else {
      this.loc.mean.df <- mean.df
    }
    if (all(this.loc.modeling.df[['var.QTL.dom']] == 0)) {
      this.loc.var.df <- var.df - 1
    } else {
      this.loc.var.df <- var.df
    }

    # Fit the alternative model at this loc
    alternative.fit <- fit_model(formulae = scan.formulae,
                                 data = this.loc.modeling.df,
                                 model = model,
                                 mean = 'alt',
                                 var = 'alt')

    # if requested, save effect estimates
    # may be safer to do with with name-matching, but this seems to work for now
    # all this hullabaloo is because coef(summary()) drops NAs, but coef() doesn't
    if (all(return.covar.effects)) {

      if (identical(alternative.fit, NA)) {

        if (loc.idx == 1) { stop('Cant fit model on locus 1.  Due to programming weirdness, cant return effect estimates.') }
        coef_mtx <- rbind(coef_mtx, rep(NA, ncol(coef_mtx)))

      } else {

        mean_coef_mtx <- stats::coef(summary(alternative.fit))
        mean_ests <- mean_ses <-  rep(NA, length(stats::coef(alternative.fit)))

        mean_ests[!is.na(stats::coef(alternative.fit))] <- mean_coef_mtx[,'Estimate']
        names(mean_ests) <- paste0('m__est__', names(stats::coef(alternative.fit)))

        mean_ses[!is.na(stats::coef(alternative.fit))] <- mean_coef_mtx[,'Std. Error']
        names(mean_ses) <- paste0('m__se__', names(stats::coef(alternative.fit)))


        disp_fit <- alternative.fit$dispersion.fit
        disp_coef_mtx <- stats::coef(summary(disp_fit))
        disp_ests <- disp_ses <-  rep(NA, length(stats::coef(disp_fit)))

        disp_ests[!is.na(stats::coef(disp_fit))] <- disp_coef_mtx[,'Estimate']
        names(disp_ests) <- paste0('v__est__', names(stats::coef(disp_fit)))

        disp_ses[!is.na(stats::coef(disp_fit))] <- disp_coef_mtx[,'Std. Error']
        names(disp_ses) <- paste0('v__se__', names(stats::coef(disp_fit)))

        # collate ests and ses
        if (loc.idx == 1) {
          coef_mtx <- matrix(data = c(mean_ests, mean_ses, disp_ests, disp_ses), nrow = 1)
        } else {
          coef_mtx <- rbind(coef_mtx, c(mean_ests, mean_ses, disp_ests, disp_ses))
        }
      }

    }

    # Fit the appropriate null models at this loc
    # and save LOD score and asymptotic p-value
    if ('mean' %in% scan.types) {
      mean.null.fit <- fit_model(formulae = scan.formulae,
                                 data = this.loc.modeling.df,
                                 model = model,
                                 mean = 'null',
                                 var = 'alt')
      result[['mean.lod']][loc.idx] <- LOD(alt = alternative.fit, null = mean.null.fit)
      result[['mean.asymp.p']][loc.idx] <- stats::pchisq(q = result[['mean.lod']][loc.idx], df = this.loc.mean.df, lower.tail = FALSE)
    }

    if ('var' %in% scan.types) {
      var.null.fit <- fit_model(formulae = scan.formulae,
                                data = this.loc.modeling.df,
                                model = model,
                                mean = 'alt',
                                var = 'null')
      result[['var.lod']][loc.idx] <- LOD(alt = alternative.fit, null = var.null.fit)
      result[['var.asymp.p']][loc.idx] <- stats::pchisq(q = result[['var.lod']][loc.idx], df = this.loc.var.df, lower.tail = FALSE)
    }

    if ('joint' %in% scan.types) {
      joint.null.fit <- fit_model(formulae = scan.formulae,
                                  data = this.loc.modeling.df,
                                  model = model,
                                  mean = 'null',
                                  var = 'null')
      result[['joint.lod']][loc.idx] <- LOD(alt = alternative.fit, null = joint.null.fit)
      result[['joint.asymp.p']][loc.idx] <- stats::pchisq(q = result[['joint.lod']][loc.idx], df = this.loc.mean.df + this.loc.var.df, lower.tail = FALSE)
    }

  }

  if (return.covar.effects) {
    result <- cbind(result, coef_mtx)
  }

  class(result) <- c('scanonevar_result', class(result))

  return(result)
}




initialize.scanonevar.result_ <- function(loc.info.df,
                                          scan.types,
                                          scan.formulae,
                                          return.covar.effects = FALSE) {

  chr.type <- chr <- loc.name <- pos <- 'fake_global_for_CRAN'

  result <- dplyr::select(.data = loc.info.df,
                          chr.type,
                          chr,
                          loc.name,
                          pos)

  if ('mean' %in% scan.types)
    result[['mean.asymp.p']] <- result[['mean.lod']] <- NA
  if ('var' %in% scan.types)
    result[['var.asymp.p']] <- result[['var.lod']] <- NA
  if ('joint' %in% scan.types)
    result[['joint.asymp.p']] <- result[['joint.lod']] <- NA

  # if (return.covar.effects) {
  #   result[['(Intercept)_mef']] <- result[['(Intercept)_mse']] <- NA
  #   for (mean.covar.name in labels(terms(scan.formulae[['mean.alt.formula']]))) {
  #     result[[paste0(mean.covar.name, '_mef')]] <- NA
  #     result[[paste0(mean.covar.name, '_mse')]] <- NA
  #   }
  #   result[['(Intercept)_vef']] <- result[['(Intercept)_vse']] <- NA
  #   for (var.covar.name in labels(terms(scan.formulae[['var.alt.formula']]))) {
  #     result[[paste0(var.covar.name, '_vef')]] <- NA
  #     result[[paste0(var.covar.name, '_vse')]] <- NA
  #   }
  # }

  return(result)
}




#' @title c.scanonevar
#' @name c.scanonevar
#'
#' @param ... the scanonevar objects with permutations to be combined
#'
#' @description combines scanonevar objects that have permutations to improve the precision of the p-value estimates.
#'
#' @return a scanonevar object that is the concatenation of the inputted
#' scanonevars
#'
#' @export
#'
c.scanonevar <- function(...) {
  sovs <- list(...)

  validate.c.scanonevar.input_(sovs)

  first.sov <- sovs[[1]]
  new.perms <- first.sov[['perms']]
  for (sov.idx in 2:(length(sovs))) {
    new.perms <- dplyr::bind_rows(new.perms,
                                  sovs[[sov.idx]][['perms']])
  }
  first.sov[['perms']] <- new.perms

  first.sov[['result']] <- calc.empir.ps(sov = first.sov[['result']],
                                         perms = new.perms)

  return(first.sov)

}


validate.c.scanonevar.input_ <- function(sovs) {

  stopifnot(all(sapply(X = sovs, FUN = is.scanonevar.w.perms)))

  first.sov <- sovs[[1]]
  for (sov.idx in 2:(length(sovs))) {
    # meta has to be the same
    stopifnot(identical(x = first.sov[['meta']], y = sovs[[sov.idx]][['meta']],
                        ignore.environment = TRUE))

    # result has to be the same
    # todo

    # perms have to have the same names...but shouldn't be the same
    stopifnot(all.equal(names(first.sov[['perms']]),
                        names(sovs[[sov.idx]][['perms']])))
  }


}


#' @title summary.scanonevar
#'
#' @author Robert Corty \email{rcorty@@gmail.com}
#'
#' @description \code{summary.scanonevar} prints out the loci in a scanonevar object
#'   that exceed \code{thresh}.  It is an S3 generic for summary().  It handles scanonevar
#'   objects in both LOD units and empirical p value units.
#'
#' @param object the scanonevar object to be summarized
#' @param thresh the threshold over which (for LODs) or under which (for emprirical p values)
#'   a locus will be printed.
#' @param units Which units should be used to summarise?  'lod', 'asymp.p', or 'empir.p'
#' @param ... additional arguments controlling the summary
#'
#' @return None.  Only prints results to screen.
#'
#' @details none
#'
#' @method summary scanonevar
#' @export
#'
summary.scanonevar <- function(object, units = c('lod', 'asymp.p', 'empir.p'), thresh, ...) {

  # hack to get R CMD CHECK to run without NOTEs that these globals are undefined
  mean.lod <- var.lod <- joint.lod <- 'fake_global_for_CRAN'
  mean.asymp.p <- var.asymp.p <- joint.asymp.p <- 'fake_global_for_CRAN'
  mean.empir.p <- var.empir.p <- joint.empir.p <- 'fake_global_for_CRAN'
  chr <- pos <- loc.name <- 'fake_global_for_CRAN'

  units <- match.arg(arg = units)

  if (!any(class(object) == "scanonevar")) {
    stop("Input should have class \"scanonevar\".")
  }

  if (missing(thresh)) {
    thresh <- switch(EXPR = units,
                     lod = 3,
                     asymp.p = 0.05,
                     empir.p = 0.05)
  }

  return <- list()

  if (units == 'lod') {

    return[['mQTL']] <- object$result %>%
      dplyr::filter(mean.lod > thresh) %>%
      dplyr::select(chr, pos, loc.name, mean.lod)

    return[['vQTL']] <- object$result %>%
      dplyr::filter(var.lod > thresh) %>%
      dplyr::select(chr, pos, loc.name, var.lod)

    return[['mvQTL']] <- object$result %>%
      dplyr::filter(joint.lod > thresh) %>%
      dplyr::select(chr, pos, loc.name, joint.lod)
  }

  if (units == 'asymp.p') {

    return[['mQTL']] <- object$result %>%
      dplyr::filter(mean.asymp.p > thresh) %>%
      dplyr::select(chr, pos, loc.name, mean.asymp.p)

    return[['vQTL']] <- object$result %>%
      dplyr::filter(var.asymp.p > thresh) %>%
      dplyr::select(chr, pos, loc.name, var.asymp.p)

    return[['mvQTL']] <- object$result %>%
      dplyr::filter(joint.asymp.p > thresh) %>%
      dplyr::select(chr, pos, loc.name, joint.asymp.p)
  }

  if (units == 'empir.p') {

    return[['mQTL']] <- object$result %>%
      dplyr::filter(mean.empir.p > thresh) %>%
      dplyr::select(chr, pos, loc.name, mean.empir.p)

    return[['vQTL']] <- object$result %>%
      dplyr::filter(var.empir.p > thresh) %>%
      dplyr::select(chr, pos, loc.name, var.empir.p)

    return[['mvQTL']] <- object$result %>%
      dplyr::filter(joint.empir.p > thresh) %>%
      dplyr::select(chr, pos, loc.name, joint.empir.p)
  }

  return(return)
}


