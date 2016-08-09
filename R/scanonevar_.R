#' @title scanonevar_
#' @name scanonevar_
#' @author Robert W. Corty
#'
#' @param modeling.df the tbl_df that contains the response, all covariates,
#' and NA columns appropriate for marker effects
#' @param genoprob.df the tbl_df that contains the probability of each genotype
#' for each individual at each loc, where a 'loc' is either a marker or a pseudomarker
#' @param loc.info.df the tbl_df that contains information about the
#' locs that are being scanned, minimally a 'loc.name', 'chr' (short for
#' 'chromosome') and 'pos' (short for 'position') for each loc.
#' @param scan.types Type(s) of scan(s) to be conducted, one of: 'mean',
#' 'var', or c('mean', 'var', 'joint')
#'
#' @inheritParams scanonevar
#'
#' @return a scanonevar object
#' @export
#'
#' @examples
#' x <- 27599
#'
scanonevar_ <- function(modeling.df,
                        genoprob.df,
                        loc.info.df,
                        scan.types,
                        scan.formulae,
                        return.covar.effects) {

  result <- initialize.scanonevar.result_(loc.info.df,
                                          scan.types,
                                          scan.formulae,
                                          return.covar.effects)

  mean.df <- sum(grepl(pattern = 'mean.QTL', x = labels(terms(scan.formulae[['mean.alt.formula']]))))
  var.df <- sum(grepl(pattern = 'var.QTL', x = labels(terms(scan.formulae[['var.alt.formula']]))))

  if ('joint' %in% scan.types) {
    joint.null.fit <- tryCatch(expr = dglm::dglm(formula = scan.formulae[['mean.null.formula']],
                                                 dformula = scan.formulae[['var.null.formula']],
                                                 data = modeling.df),
                               warning = function(w) NA,
                               error = function(e) NA,
                               finally = NA)
  }


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

    # Fit the alternative model at this loc
    alternative.fit <- tryCatch(expr = dglm::dglm(formula = scan.formulae[['mean.alt.formula']],
                                                  dformula = scan.formulae[['var.alt.formula']],
                                                  data = this.loc.modeling.df),
                                warning = function(w) NA,
                                error = function(e) NA,
                                finally = NA)

    # if requested, save effect estimates
    # may be safer to do with with name-matching, but this seems to work for now
    if (all(return.covar.effects, !identical(alternative.fit, NA))) {
      result[loc.idx, paste0(names(coef(alternative.fit)), '_mef')] <- coef(alternative.fit)
      result[loc.idx, paste0(names(coef(alternative.fit$dispersion.fit)), '_vef')] <- coef(alternative.fit$dispersion.fit)
    }

    # Fit the appropriate null models at this loc
    # and save LOD score and asymptotic p-value
    if ('mean' %in% scan.types) {
      mean.null.fit <- tryCatch(expr = dglm::dglm(formula = scan.formulae[['mean.null.formula']],
                                                  dformula = scan.formulae[['var.alt.formula']],
                                                  data = this.loc.modeling.df),
                                warning = function(w) NA,
                                error = function(e) NA,
                                finally = NA)
      result[['mean.lod']][loc.idx] <- LOD(alt = alternative.fit, null = mean.null.fit)
      result[['mean.asymp.p']][loc.idx] <- pchisq(q = result[['mean.lod']][loc.idx], df = this.loc.mean.df, lower.tail = FALSE)
    }

    if ('var' %in% scan.types) {
      var.null.fit <- tryCatch(expr = dglm::dglm(formula = scan.formulae[['mean.alt.formula']],
                                                 dformula = scan.formulae[['var.null.formula']],
                                                 data = this.loc.modeling.df),
                               warning = function(w) NA,
                               error = function(e) NA,
                               finally = NA)
      result[['var.lod']][loc.idx] <- LOD(alt = alternative.fit, null = var.null.fit)
      result[['var.asymp.p']][loc.idx] <- pchisq(q = result[['var.lod']][loc.idx], df = this.loc.var.df, lower.tail = FALSE)
    }

    if ('joint' %in% scan.types) {
      result[['joint.lod']][loc.idx] <- LOD(alt = alternative.fit, null = joint.null.fit)
      result[['joint.asymp.p']][loc.idx] <- pchisq(q = result[['joint.lod']][loc.idx], df = this.loc.mean.df + this.loc.var.df, lower.tail = FALSE)
    }

  }


  attr(x = result, which = 'scan.types') <- scan.types
  attr(x = result, which = 'mean.formula') <- scan.formulae[['mean.alt.formula']]
  attr(x = result, which = 'var.formula') <- scan.formulae[['var.alt.formula']]

  return(result)
}






#' @title initialize.scanonevar.result_
#' @rdname scanonevar_
#'
#' @inheritParams scanonevar_
#'
#' @return A tbl_df with the metadata for a scanonevar, once filled in by
#' the scanonevar function, it will be a scanonevar object.
#'
#' @details This internal function should not typically be called by a user.
#'
#' @examples
#' x <- dplyr::data_frame(loc.name = c('m1', 'm2'),
#'                        chr = c(1, 1),
#'                        pos = c(5, 10))
#' initialize.scanonevar.result_(loc.info.df = x,
#'                               scan.types = 'mean')
#'
initialize.scanonevar.result_ <- function(loc.info.df,
                                          scan.types,
                                          scan.formulae,
                                          return.covar.effects) {

  result <- dplyr::select(.data = loc.info.df,
                          loc.name,
                          chr,
                          pos)

  if ('mean' %in% scan.types)
    result[['mean.asymp.p']] <- result[['mean.lod']] <- NA
  if ('var' %in% scan.types)
    result[['var.asymp.p']] <- result[['var.lod']] <- NA
  if ('joint' %in% scan.types)
    result[['joint.asymp.p']] <- result[['joint.lod']] <- NA

  if (return.covar.effects) {
    result[['(Intercept)_mef']] <- NA
    for (mean.covar.name in labels(terms(scan.formulae[['mean.alt.formula']]))) {
      result[[paste0(mean.covar.name, '_mef')]] <- NA
    }
    result[['(Intercept)_vef']] <- NA
    for (var.covar.name in labels(terms(scan.formulae[['var.alt.formula']]))) {
      result[[paste0(var.covar.name, '_vef')]] <- NA
    }
  }

  return(result)
}



#' @title make.loc.specific.modeling.df
#' @rdname scanonevar_
#'
#' @param general.modeling.df The modeling df that will be used across the scan,
#' should have NA for all QTL columns.  These columns will be filled in in the result.
#' @param loc.genoprobs The genoprobs of all individuals of all allels at this loc.
#' @param model.formulae The formulae used for the scan.
#'
#' @return a modeling df appropriate for modeling at one specific loc
#' @export
#'
make.loc.specific.modeling.df <- function(general.modeling.df,
                                          loc.genoprobs,
                                          model.formulae) {

  modeling.df <- general.modeling.df

  if ('mean.QTL.add' %in% labels(terms(model.formulae[['mean.alt.formula']]))) {
    modeling.df[['mean.QTL.add']] <- additive.component(genoprobs.long = loc.genoprobs)
  }
  if ('mean.QTL.dom' %in% labels(terms(model.formulae[['mean.alt.formula']]))) {
    modeling.df[['mean.QTL.dom']] <- dominance.component(genoprobs.long = loc.genoprobs)
  }
  if ('var.QTL.add' %in% labels(terms(model.formulae[['var.alt.formula']]))) {
    modeling.df[['var.QTL.add']] <- additive.component(genoprobs.long = loc.genoprobs)
  }
  if ('var.QTL.dom' %in% labels(terms(model.formulae[['var.alt.formula']]))) {
    modeling.df[['var.QTL.dom']] <- dominance.component(genoprobs.long = loc.genoprobs)
  }

  return(modeling.df)
}


