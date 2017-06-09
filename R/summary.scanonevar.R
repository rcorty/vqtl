#' @title summary.scannevar
#' @rdname scanonevar
#'
#'
#'
#' @param object the scanonevar to summarize
#' @param units the units
#' @param thresh the threshold
#' @param ... ignored
#'
#' @return a summary of the scanonevar object
#' @export
#'
summary.scanonevar <- function(object, units = c('LOD', 'empir.p'), thresh, ...) {

  chr.type <- chr <- loc.name <- pos <- 'fake_global_for_CRAN'
  mean.empir.p <- var.empir.p <- joint.empir.p <- 'fake_global_for_CRAN'

  stopifnot(is.scanonevar(object))
  units <- match.arg(units)
  if (missing(thresh)) {
    thresh <- switch(EXPR = units,
                     'LOD' = 3,
                     'empir.p' = -log10(0.05))
  }

  output <- list()

  phen_name <- object$meta$formulae$mean.alt.formula[[2]]

  if ('perms' %in% names(object)) {
    output[['intro']] <- paste('a scanonevar of phenotype', phen_name, 'with permutations')
  } else {
    output[['intro']] <- paste('a scanonevar of phenotype', phen_name, 'without permutations')
  }

  result <- object[['result']]
  if (units == 'LOD') {
    output[['mean.peaks']] <- result %>%
      dplyr::slice(get.peak.idxs(v = result[['mean.lod']], thresh = thresh)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('mean'))

    output[['var.peaks']] <- result %>%
      dplyr::slice(get.peak.idxs(v = result[['var.lod']], thresh = thresh)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('var'))

    output[['joint.peaks']] <- result %>%
      dplyr::slice(get.peak.idxs(v = result[['joint.lod']], thresh = thresh)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('joint'))
  }

  if (units == 'empir.p') {
    # output[['mean.peaks']] <- result %>%
    #   dplyr::slice(get.peak.idxs(v = -log10(result[['mean.empir.p']]), thresh = thresh)) %>%
    #   dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('mean'))

    output[['mean_peaks']] <- result %>% dplyr::filter(-log10(mean.empir.p) > thresh) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('mean'))

    # output[['var.peaks']] <- result %>%
    #   dplyr::slice(get.peak.idxs(v = -log10(result[['var.empir.p']]), thresh = thresh)) %>%
    #   dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('var'))

    output[['var_peaks']] <- result %>% dplyr::filter(-log10(var.empir.p) > thresh) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('var'))

    # output[['joint.peaks']] <- result %>%
    #   dplyr::slice(get.peak.idxs(v = -log10(result[['joint.empir.p']]), thresh = thresh)) %>%
    #   dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('joint'))

    output[['joint_peaks']] <- result %>% dplyr::filter(-log10(joint.empir.p) > thresh) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('joint'))
  }

  return(output)
}




get.peak.idxs <- function(v, thresh) {

  above.thresh <- prev.above.thresh <- next.above.thresh <- peak.start <- peak.end <- cum.starts <- cum.ends <- idx <- 'fake_global_for_CRAN'
  dplyr::data_frame(idx = 1:length(v), v = v) %>%
    dplyr::mutate(above.thresh = v > thresh,
           prev.above.thresh = dplyr::lag(above.thresh, default = FALSE),
           next.above.thresh = dplyr::lead(above.thresh, default = FALSE),
           peak.start = above.thresh & !prev.above.thresh,
           peak.end = above.thresh & !next.above.thresh,
           cum.starts = cumsum(peak.start),
           cum.ends = cumsum(peak.end)) %>%
    dplyr::filter(cum.starts != cum.ends) %>%
    dplyr::pull(idx)

  # changed interface to simply return idxs
    # group_by(cum.starts) %>%
    # summarise(peak.start = min(idx),
    #           peak.end = max(idx),
    #           peak.height = max(v)) %>%
    # select(-cum.starts)

}