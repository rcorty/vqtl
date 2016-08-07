#' @title scanonevar_
#' @name scanonevar_
#' @author Robert W. Corty
#'
#' @param modeling.df the tbl_df that contains the response, all covariates,
#' and NA columns appropriate for marker effects
#' @param genoprob.df the tbl_df that contains the probability of each genotype
#' for each individual at each marker or pseudomarker
#' @param marker.info.df the tbl_df that contains information about the
#' markers that are being scanned, minimally a 'marker.name', 'chr' (short for
#' 'chromosome') and 'pos' (short for 'position') for each marker.
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
                        marker.info.df,
                        scan.types,
                        scan.formulae) {

  # todo: get rid of this is marker.info.df ends up having just these three cols
  result <- initialize.scanonevar.result_(marker.info.df,
                                          scan.types)

  mean.df <- sum(grepl(pattern = 'mean.QTL', x = labels(terms(scan.formulae[['mean.alt.formula']]))))
  var.df <- sum(grepl(pattern = 'var.QTL', x = labels(terms(scan.formulae[['var.alt.formula']]))))

  if ('joint' %in% scan.types) {
    joint.null.fit <- dglm::dglm(formula = scan.formulae[['mean.null.formulae']],
                                 dformula = scan.formulae[['var.null.formulae']],
                                 data = modeling.df)
  }


  for (marker.idx in 1:nrow(result)) {

    # fill modeling.df with the genoprobs at the focal marker or pseudomarker
    marker.name <- result[['marker.name']][marker.idx]
    marker.genoprobs <- dplyr::filter(.data = genoprob.df,
                                      marker.name == marker.name)
    # need to widen this since it's coming in long format

    this.marker.modeling.df <- make.marker.specific.modeling.df(general.modeling.df = modeling.df,
                                                                marker.genoprobs = marker.genoprobs,
                                                                model.formulae = scan.formulae)

    # Fit the alternative model at this marker
    alternative.fit <- dglm::dglm(formula = scan.formulae[['mean.alt.formula']],
                                  dformula = scan.formulae[['var.alt.formula']],
                                  data = this.marker.modeling.df)


    # Fit the appropriate null models at this marker
    # and save LOD score and asymptotic p-value
    if ('mean' %in% scan.types) {
      mean.null.fit <- dglm::dglm(formula = scan.formulae[['mean.null.formula']],
                                  dformula = scan.formulae[['var.alt.formula']],
                                  data = this.marker.modeling.df)
      result[['mean.lod ']] <- LOD(alt = alternative.fit, null = mean.null.fit)
      result[['mean.asymp.p']] <- pchisq(q = result[['mean.lod ']], df = mean.df, lower.tail = FALSE)
    }

    if ('var' %in% scan.types) {
      var.null.fit <- dglm::dglm(formula = scan.formulae[['mean.alt.formula']],
                                 dformula = scan.formulae[['var.null.formula']],
                                 data = this.marker.modeling.df)
      result[['mean.lod ']] <- LOD(alt = alternative.fit, null = var.null.fit)
      result[['var.asymp.p']] <- pchisq(q = result[['var.lod ']], df = var.df, lower.tail = FALSE)
    }

    if ('joint' %in% scan.types) {
      result[['joint.p']] <- LRS(alt = alternative.fit,
                                 null = joint.null.fit)
      result[['joint.lod ']] <- LOD(alt = alternative.fit, null = joint.null.fit)
      result[['joint.asymp.p']] <- pchisq(q = result[['joint.lod ']], df = mean.df + var.df, lower.tail = FALSE)
    }

  }


  attr(x = result, which = 'scan.types') <- scan.types
  attr(x = result, which = 'mean.formula') <- mean.formula
  attr(x = result, which = 'var.formula') <- var.formual

  return(result)
}






#' @title initialize.scanonevar.result_
#' @name initizlize.scanonevar.result_
#' @author Robert W. Corty
#'
#' @inheritParams scanonevar_
#'
#' @return A tbl_df with the metadata for a scanonevar, once filled in by
#' the scanonevar function, it will be a scanonevar object.
#'
#' @details This internal function should not typically be called by a user.
#'
#' @examples
#' x <- dplyr::data_frame(marker.name = c('m1', 'm2'),
#'                        chr = c(1, 1),
#'                        pos = c(5, 10))
#' initialize.scanonevar.result_(marker.info.df = x,
#'                               scan.types = 'mean')
#'
initialize.scanonevar.result_ <- function(marker.info.df,
                                          scan.types) {

  result <- dplyr::select(.data = marker.info.df,
                          marker.name,
                          chr,
                          pos)

  if ('mean' %in% scan.types)
    result[['mean.lod']] <- result[['mean.asymp.p']] <- NA
  if ('var' %in% scan.types)
    result[['var.lod']] <- result[['var.asymp']] <- NA
  if ('joint' %in% scan.types)
    result[['joint.lod']] <- result[['joint.asymp.p']] <- NA

  return(result)
}



#' @title LOD
#' @name LOD
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @param alt alternative DGLM
#' @param null null DGLM
#'
#' @return log( likelihood of alternative model / likelihood of null model)
#' @export
#'
#' @examples
#' x <- 27599
LOD <- function(alt, null) {

  stopifnot('dglm' %in% class(alt))
  stopifnot('dglm' %in% class(null))

  return(0.5*(null$m2loglik - alt$m2loglik ))
}





make.marker.specific.modeling.df <- function(general.modeling.df,
                                             marker.genoprobs,
                                             model.formula) {

  if ('mean.QTL.add' %in% labels(terms(model.formula[['mean.alt.formula']]))) {
    general.modeling.df[['mean.QTL.add']] <- additive.component(marker.genoprobs)
  }
  if ('mean.QTL.dom' %in% labels(terms(model.formula[['mean.alt.formula']]))) {
    general.modeling.df[['mean.QTL.dom']] <- dominance.component(marker.genoprobs)
  }
  if ('var.QTL.add' %in% labels(terms(model.formula[['var.alt.formula']]))) {
    general.modeling.df[['var.QTL.add']] <- additive.component(marker.genoprobs)
  }
  if ('var.QTL.dom' %in% labels(terms(model.formula[['var.alt.formula']]))) {
    general.modeling.df[['var.QTL.dom']] <- dominance.component(marker.genoprobs)
  }

  return(modeling.df)
}