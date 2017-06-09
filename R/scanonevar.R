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
  meta <- pull.scanonevar.meta_(wrangled.inputs,
                                chrs)

  # execute the scan
  result <- scanonevar_(modeling.df = wrangled.inputs$modeling.df,
                        loc.info.df = wrangled.inputs$loc.info.df,
                        genoprob.df = wrangled.inputs$genoprob.df,
                        scan.types = wrangled.inputs$scan.types,
                        scan.formulae = wrangled.inputs$scan.formulae,
                        return.covar.effects = return.covar.effects)

  sov <- list(meta = meta,
              result = result)

  class(sov) <- c('scanonevar', class(sov))
  return(sov)
}




validate.scanonevar.input_ <- function(cross,
                                       mean.formula,
                                       var.formula,
                                       chrs) {

  # each argument must be individually valid
  stopifnot(is.cross(cross))
  stopifnot(is.mean.formula(mean.formula))
  stopifnot(is.var.formula(var.formula))
  stopifnot(all(chrs %in% qtl::chrnames(cross)))

  formulae <- make.formulae_(mean.formula, var.formula)

  # formulae must be valid for use in scanonevar
  stopifnot(formulae_is_valid_(formulae = formulae))

  # formulae must be valid for use with cross
  stopifnot(formulae.is.valid.for.cross_(cross = cross,
                                         formulae = formulae))

  return(TRUE)
}




formulae_is_valid_ <- function(formulae) {

  if (!is.formulae(formulae)) {
    return(FALSE)
  }

  # must have at least one QTL term used appropriately
  mean.covars <- all.vars(formulae[['mean.formula']][[3]])
  var.covars <- all.vars(formulae[['var.formula']])
  if (all(!any(c('mean.QTL.add', 'mean.QTL.dom') %in% mean.covars),
          !any(c('var.QTL.add', 'var.QTL.dom') %in% var.covars))) {
    return(FALSE)
  }

  return(TRUE)
}




wrangle.scanonevar.input_ <- function(cross,
                                      mean.formula,
                                      var.formula,
                                      chrs) {

  if (!is.cross.w.genoprobs(x = cross)) {
    message("calculating genoprobs with stepwidth = 2, off.end = 0, error.prob = 1e-4, map.function = 'haldane'")
    cross <- qtl::calc.genoprob(cross = cross, step = 5)
  }

  loc.info.df <- wrangle.loc.info.df_(cross = cross,
                                      chrs = chrs)

  genoprob.df.long <- wrangle.genoprob.df_(cross = cross)

  scan.types <- wrangle.scanonevar.types_(mean.formula, var.formula)

  scan.formulae <- wrangle.scanonevar.formulae_(cross, mean.formula, var.formula)

  modeling.df <- wrangle.scanonevar.modeling.df_(cross = cross,
                                                 genoprobs = genoprob.df.long,
                                                 scan.formulae = scan.formulae)

  return(list(cross = cross,
              scan.types = scan.types,
              scan.formulae = scan.formulae,
              modeling.df = modeling.df,
              loc.info.df = loc.info.df,
              genoprob.df = genoprob.df.long))
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
               chrs = chrs)
  class(meta) <- c('scanonevar.meta', class(meta))

  return(meta)
}



scanonevar_ <- function(modeling.df,
                        genoprob.df,
                        loc.info.df,
                        scan.types,
                        scan.formulae,
                        return.covar.effects) {

  loc.name <- 'fake_global_for_CRAN'

  result <- initialize.scanonevar.result_(loc.info.df = loc.info.df,
                                          scan.types = scan.types,
                                          scan.formulae = scan.formulae,
                                          return.covar.effects = return.covar.effects)

  mean.df <- sum(grepl(pattern = 'mean.QTL', x = labels(stats::terms(scan.formulae[['mean.alt.formula']]))))
  var.df <- sum(grepl(pattern = 'var.QTL', x = labels(stats::terms(scan.formulae[['var.alt.formula']]))))

  # do this outside the loop because it doesn't change with locus, so only needs to be done once
  if ('joint' %in% scan.types) {
    joint.null.fit <- fit.model.00_(formulae = scan.formulae, df = modeling.df)
  }

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
                                                          model.formulae = scan.formulae)

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
    alternative.fit <- fit.model.mv_(formulae = scan.formulae, df = this.loc.modeling.df)

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
        names(mean_ests) <- paste0('mef_', names(stats::coef(alternative.fit)))

        mean_ses[!is.na(stats::coef(alternative.fit))] <- mean_coef_mtx[,'Std. Error']
        names(mean_ses) <- paste0('mse_', names(stats::coef(alternative.fit)))


        disp_fit <- alternative.fit$dispersion.fit
        disp_coef_mtx <- stats::coef(summary(disp_fit))
        disp_ests <- disp_ses <-  rep(NA, length(stats::coef(disp_fit)))

        disp_ests[!is.na(stats::coef(disp_fit))] <- disp_coef_mtx[,'Estimate']
        names(disp_ests) <- paste0('vef_', names(stats::coef(disp_fit)))

        disp_ses[!is.na(stats::coef(disp_fit))] <- disp_coef_mtx[,'Std. Error']
        names(disp_ses) <- paste0('vse_', names(stats::coef(disp_fit)))

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
      mean.null.fit <- fit.model.0v_(formulae = scan.formulae, df = this.loc.modeling.df)
      result[['mean.lod']][loc.idx] <- LOD(alt = alternative.fit, null = mean.null.fit)
      result[['mean.asymp.p']][loc.idx] <- stats::pchisq(q = result[['mean.lod']][loc.idx], df = this.loc.mean.df, lower.tail = FALSE)
    }

    if ('var' %in% scan.types) {
      var.null.fit <- fit.model.m0_(formulae = scan.formulae, df = this.loc.modeling.df)
      result[['var.lod']][loc.idx] <- LOD(alt = alternative.fit, null = var.null.fit)
      result[['var.asymp.p']][loc.idx] <- stats::pchisq(q = result[['var.lod']][loc.idx], df = this.loc.var.df, lower.tail = FALSE)
    }

    if ('joint' %in% scan.types) {
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



make.loc.specific.modeling.df <- function(general.modeling.df,
                                          loc.genoprobs,
                                          model.formulae) {

  modeling.df <- general.modeling.df

  if ('mean.QTL.add' %in% labels(stats::terms(model.formulae[['mean.alt.formula']]))) {
    modeling.df[['mean.QTL.add']] <- additive.component_(genoprobs.long = loc.genoprobs)
  }
  if ('mean.QTL.dom' %in% labels(stats::terms(model.formulae[['mean.alt.formula']]))) {
    modeling.df[['mean.QTL.dom']] <- dominance.component_(genoprobs.long = loc.genoprobs)
  }
  if ('var.QTL.add' %in% labels(stats::terms(model.formulae[['var.alt.formula']]))) {
    modeling.df[['var.QTL.add']] <- additive.component_(genoprobs.long = loc.genoprobs)
  }
  if ('var.QTL.dom' %in% labels(stats::terms(model.formulae[['var.alt.formula']]))) {
    modeling.df[['var.QTL.dom']] <- dominance.component_(genoprobs.long = loc.genoprobs)
  }

  return(modeling.df)
}





