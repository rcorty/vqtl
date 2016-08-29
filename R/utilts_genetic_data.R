#' @title wrangle.loc.info.df_
#' @rdname internals
#'
#' @inheritParams scanonevar
#'
wrangle.loc.info.df_ <- function (cross, chrs = qtl::chrnames(cross)) {

  # NB: names(class(x)) is the notation to get the name of a chromosome
  loc.info.from.chr <- function(x) {
    chr.type <- class(x)
    chr.name <- names(class(x))
    prob.map <- attr(x = x[['prob']], which = 'map')
    names.starting.with.loc.idxs <- grep(pattern = '^loc', names(prob.map))
    names(prob.map)[names.starting.with.loc.idxs] <-
      paste0('chr', chr.name, '_', names(prob.map)[names.starting.with.loc.idxs])

    return(dplyr::data_frame(chr.type = chr.type,
                             chr = chr.name,
                             loc.name = names(prob.map),
                             pos = prob.map))
  }

  cross[['geno']] <- cross[['geno']][qtl::chrnames(cross = cross) %in% chrs]
  return(dplyr::bind_rows(lapply(X = cross[['geno']],
                                 FUN = loc.info.from.chr)))
}



#' @title wrangle.genoprob.df_
#' @rdname internals
#'
#' @inheritParams scanonevar
#'
wrangle.genoprob.df_ <- function(cross) {

  genoprobs.from.chr <- function(x) {
    prob.tbl <- x[['prob']]
    num.width <- max(nchar(as.character(1:dim(prob.tbl)[1])))
    if (is.null(dimnames(prob.tbl)[[1]])) {
      dimnames(prob.tbl)[[1]] <- paste0('org',
                                        stringr::str_pad(string = 1:dim(prob.tbl)[1],
                                                         width = num.width,
                                                         pad = '0'))
    }
    names.starting.with.loc.idxs <- grep(pattern = '^loc', dimnames(prob.tbl)[[2]])
    dimnames(prob.tbl)[[2]][names.starting.with.loc.idxs] <-
      paste0('chr', names(class(x)), '_', dimnames(prob.tbl)[[2]][names.starting.with.loc.idxs])
    return(as.data.frame.table(prob.tbl, stringsAsFactors = FALSE))
  }

  genoprob.df <- dplyr::tbl_df(dplyr::bind_rows(lapply(X = cross[['geno']],
                                                       FUN = genoprobs.from.chr)))
  names(genoprob.df) <- c('iid', 'loc.name', 'allele', 'genoprob')
  return(genoprob.df)
}



#' @title make.response.model.df_
#' @rdname internals
#'
#' @inheritParams scanonevar
#'
#' @return a tibble of the response in mean.formua
#' @export
#'
make.response.model.df_ <- function(cross,
                                    formulae) {

  stopifnot(is.cross(cross))

  response.name <- all.vars(formulae[['mean.alt.formula']][[2]])

  # not sure if there's a better way to do this
  df <- dplyr::as_data_frame(stats::setNames(list(qtl::pull.pheno(cross = cross,
                                                                  pheno.col = response.name)),
                                             response.name))

  return(df)
}



#' @title make.qtl.covar.model.df
#' @rdname internals
#'
#' @inheritParams scanonevar
#'
#' @return a tibble of the qtl keyword in formulae
#' @export
make.qtl.covar.model.df_ <- function(cross,
                                     formulae) {

  stopifnot(is.cross(cross))

  # get all covariate names
  mean.covar.names <- labels(terms(formulae[['mean.alt.formula']]))
  var.covar.names <- labels(terms(formulae[['var.alt.formula']]))

  # get the covariate names that are keywords and add them as NA to l
  mean.keywords <- c('mean.QTL.add', 'mean.QTL.dom')
  mean.keyword.covar.names <- mean.covar.names[mean.covar.names %in% mean.keywords]
  var.keywords <- c('var.QTL.add', 'var.QTL.dom')
  var.keyword.covar.names <- var.covar.names[var.covar.names %in% var.keywords]

  df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  for (keyword in c(mean.keyword.covar.names, var.keyword.covar.names)) {
    df[[keyword]] <- rep(NA, qtl::nind(cross))
  }
  df[['placeholder']] <- NULL

  return(df)
}



make.phen.covar.model.df_ <- function(cross,
                                      formulae) {

  stopifnot(is.cross(cross))

  # get all covariate names
  mean.covar.names <- labels(terms(formulae[['mean.alt.formula']]))
  var.covar.names <- labels(terms(formulae[['var.alt.formula']]))

  # get the phenotype names
  mean.phen.covar.names <- mean.covar.names[mean.covar.names %in% names(cross[['pheno']])]
  var.phen.covar.names <- var.covar.names[var.covar.names %in% names(cross[['pheno']])]

  df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  for (phen.covar.name in unique(c(mean.phen.covar.names, var.phen.covar.names))) {
    df[[phen.covar.name]] <- cross[['pheno']][[phen.covar.name]]
  }
  df[['placeholder']] <- NULL

  return(df)
}






make.genet.covar.add.dom.model.df_ <- function(cross,
                                               formulae,
                                               genoprobs) {

  stopifnot(is.cross(cross))

  # get all covariate names
  mean.covar.names <- labels(terms(formulae[['mean.alt.formula']]))
  var.covar.names <- labels(terms(formulae[['var.alt.formula']]))

  # set up the names we're looking for
  marker.names <- colnames(qtl::pull.geno(cross))
  add.marker.names <- paste0(marker.names, '_add')
  dom.marker.names <- paste0(marker.names, '_dom')

  # get the covariate names that are marker.name_add or marker.name_dom and add them to modeling.df
  # There should be no covariates that are just marker.name at this point
  mean.add.marker.covar.names <- mean.covar.names[mean.covar.names %in% add.marker.names]
  mean.dom.marker.covar.names <- mean.covar.names[mean.covar.names %in% dom.marker.names]
  var.add.marker.covar.names <- var.covar.names[var.covar.names %in% add.marker.names]
  var.dom.marker.covar.names <- var.covar.names[var.covar.names %in% dom.marker.names]


  # turn it into a tibble
  df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  for (add.marker.covar.name in unique(c(mean.add.marker.covar.names, var.add.marker.covar.names))) {
    this.marker.genoprobs <- dplyr::filter(.data = genoprobs,
                                           loc.name == substr(x = add.marker.covar.name,
                                                              start = 1,
                                                              stop = nchar(add.marker.covar.name) - 4))
    df[[add.marker.covar.name]] <- additive.component(this.marker.genoprobs)
  }


  for (dom.marker.covar.name in unique(c(mean.dom.marker.covar.names, var.dom.marker.covar.names))) {
    this.marker.genoprobs <- dplyr::filter(.data = genoprobs,
                                           loc.name == substr(x = dom.marker.covar.name,
                                                              start = 1,
                                                              stop = nchar(dom.marker.covar.name) - 4))
    df[[dom.marker.covar.name]] <- dominance.component(this.marker.genoprobs)
  }
  df[['placeholder']] <- NULL

  return(df)
}






#' @title additive.component
#' @rdname internals
#'
#' @param genoprobs.long The genoprobs from which the additive genetic
#' component should be extracted.  The 'long' format implies that each row
#' has information on the genoprob of one individual having one allele.
#'
#' @return A vector of the additive genetic component at the locus
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
#' @rdname internals
#'
#' @param genoprobs.long The genoprobs from which the additive genetic
#' component should be extracted.  The 'long' format implies that each row
#' has information on the genoprob of one individual having one allele.
#'
#' @return A vector of the dominance genetic component at the locus
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