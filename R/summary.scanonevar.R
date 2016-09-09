#' @title summary.scannevar
#' @rdname summary.scanonevar
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

  stopifnot(is.scanonevar(object))
  units <- match.arg(units)
  if (missing(thresh)) {
    thresh <- ifelse(test = units == 'LOD',
                     yes = 4,
                     no = 3)
  }

  output <- list()

  if ('perms' %in% names(object)) {
    output[['intro']] <- 'a scanonevar object with permutations'
  } else {
    output[['intro']] <- 'a scanonevar object with no permutations'
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
    output[['mean.peaks']] <- result %>%
      dplyr::slice(get.peak.idxs(v = -log10(result[['mean.empir.p']]), thresh = thresh)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('mean'))

    output[['var.peaks']] <- result %>%
      dplyr::slice(get.peak.idxs(v = -log10(result[['var.empir.p']]), thresh = thresh)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('var'))

    output[['joint.peaks']] <- result %>%
      dplyr::slice(get.peak.idxs(v = -log10(result[['joint.empir.p']]), thresh = thresh)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, dplyr::matches('joint'))
  }

  return(output)
}




get.peak.idxs <- function(v, thresh) {

  dplyr::data_frame(idx = 1:length(v), v = v) %>%
    dplyr::mutate(above.thresh = v > thresh,
           prev.above.thresh = dplyr::lag(above.thresh, default = FALSE),
           next.above.thresh = dplyr::lead(above.thresh, default = FALSE),
           peak.start = above.thresh & !prev.above.thresh,
           peak.end = above.thresh & !next.above.thresh,
           cum.starts = cumsum(peak.start),
           cum.ends = cumsum(peak.end)) %>%
    dplyr::filter(cum.starts != cum.ends) %>%
    .[['idx']]
    # changing interface to simply return idxs
    # group_by(cum.starts) %>%
    # summarise(peak.start = min(idx),
    #           peak.end = max(idx),
    #           peak.height = max(v)) %>%
    # select(-cum.starts)

}