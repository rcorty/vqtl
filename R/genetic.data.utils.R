#' @title additive.component
#' @rdname scanonevar_
#'
#' @param genoprobs.long The genoprobs from which the additive genetic
#' component should be extracted.  The 'long' format implies that each row
#' has information on the genoprob of one individual having one allele.
#'
#' @return A vector of the additive genetic component at the locus
#' @export
#'
additive.component <- function(genoprobs.long) {

  alleles <- unique(genoprobs.long[['allele']])

  genoprobs.wide <- tidyr::spread(data = genoprobs.long,
                                  key = allele,
                                  value = genoprob)

  # genoprobs.wide2 <- dplyr::select(.data = genoprobs.wide, one_of(alleles))

  if (all(alleles %in% c('AA', 'AB', 'BB'))) {
    return(genoprobs.wide[['AA']] - genoprobs.wide[['BB']])
  } else if (all(alleles %in% c('g1', 'g2'))) {
    return(genoprobs.wide[['g2']])
  } else {
    stop(paste("Can't determine additive component of loc with alleles:", alleles))
  }
}


#' @title dominance.component
#' @rdname scanonevar_
#'
#' @param genoprobs.long The genoprobs from which the additive genetic
#' component should be extracted.  The 'long' format implies that each row
#' has information on the genoprob of one individual having one allele.
#'
#' @return A vector of the dominance genetic component at the locus
#' @export
#'
dominance.component <- function(genoprobs.long) {

  alleles <- unique(genoprobs.long[['allele']])

  genoprobs.wide <- tidyr::spread(data = genoprobs.long,
                                  key = allele,
                                  value = genoprob)

  # genoprobs.wide2 <- dplyr::select(.data = genoprobs.wide, one_of(alleles))

  if (all(alleles %in% c('AA', 'AB', 'BB'))) {
    return(genoprobs.wide[['AB']])
  } else if (all(alleles %in% c('g1', 'g2'))) {
    return(0)
  } else {
    stop(paste("Can't determine additive component of loc with alleles:", alleles))
  }
}