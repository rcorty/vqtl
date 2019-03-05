wrangle.loc.info.df_ <- function (cross, chrs = qtl::chrnames(cross)) {

  # NB: names(class(x)) is the notation to get the name of a chromosome
  loc.info.from.chr.name <- function(chr.name) {
    x <- cross[['geno']][[chr.name]]
    chr.type <- class(x)
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
  return(dplyr::bind_rows(lapply(X = chrs,
                                 FUN = loc.info.from.chr.name)))
}



wrangle.genoprob.df_ <- function(cross, chrs = qtl::chrnames(cross)) {

  genoprobs.from.chr.name <- function(chr.name) {
    x <- cross[['geno']][[chr.name]]
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
      paste0('chr', chr.name, '_', dimnames(prob.tbl)[[2]][names.starting.with.loc.idxs])
    return(as.data.frame.table(prob.tbl, stringsAsFactors = FALSE))
  }

  genoprob.df <- dplyr::tbl_df(dplyr::bind_rows(lapply(X = chrs,
                                                       FUN = genoprobs.from.chr.name)))
  names(genoprob.df) <- c('iid', 'loc.name', 'allele', 'genoprob')
  return(genoprob.df)
}



make.response.model.df_ <- function(cross,
                                    formulae = NULL,
                                    response.name = all.vars(formulae[['mean.alt.formula']][[2]])) {

  stopifnot(is.cross(cross))

  # not sure if there's a better way to do this
  # df <- dplyr::as_data_frame(stats::setNames(list(qtl::pull.pheno(cross = cross,
  #                                                                 pheno.col = response.name)),
  #                                            response.name))
  #
  df <- stats::setNames(object = dplyr::tbl_df(qtl::pull.pheno(cross = cross,
                                                               pheno.col = response.name)),
                        nm = response.name)

  return(df)
}



make.qtl.covar.model.df_ <- function(cross,
                                     formulae) {

  stopifnot(is.cross(cross))

  # get all covariate names
  mean.covar.names <- labels(stats::terms(formulae[['mean.alt.formula']]))
  var.covar.names <- labels(stats::terms(formulae[['var.alt.formula']]))

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
                                      formulae,
                                      phen.names = NULL) {

  stopifnot(is.cross(cross))

  if (is.null(phen.names)) {
    # get all covariate names
    mean.covar.names <- all.vars(formulae[['mean.alt.formula']][[3]])
    var.covar.names <- all.vars(formulae[['var.alt.formula']])

    # get the phenotype names
    mean.phen.covar.names <- mean.covar.names[mean.covar.names %in% names(cross[['pheno']])]
    var.phen.covar.names <- var.covar.names[var.covar.names %in% names(cross[['pheno']])]

    phen.names <- unique(c(mean.phen.covar.names, var.phen.covar.names))
  }

  df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  for (phen.covar.name in phen.names) {
    # df <- cbind(df, model.matrix(object = formula(paste('~', phen.covar.name)),
    #                              data = cross[['pheno']]))
    df[[phen.covar.name]] <- cross[['pheno']][[phen.covar.name]]
  }
  df[['placeholder']] <- NULL

  return(df)
}






make.genet.covar.add.dom.model.df_ <- function(cross,
                                               formulae,
                                               genoprobs) {

  loc.name <- 'fake_global_for_CRAN'

  stopifnot(is.cross(cross))

  # get all covariate names
  mean.covar.names <- labels(stats::terms(formulae[['mean.alt.formula']]))
  var.covar.names <- labels(stats::terms(formulae[['var.alt.formula']]))

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
    df[[add.marker.covar.name]] <- additive.component_(genoprobs.long = this.marker.genoprobs,
                                                       cross_type = class(cross)[1])
  }
  df[['placeholder']] <- NULL

  if (class(cross)[1] %in% 'f2') {

    for (dom.marker.covar.name in unique(c(mean.dom.marker.covar.names, var.dom.marker.covar.names))) {
      this.marker.genoprobs <- dplyr::filter(.data = genoprobs,
                                             loc.name == substr(x = dom.marker.covar.name,
                                                                start = 1,
                                                                stop = nchar(dom.marker.covar.name) - 4))
      df[[dom.marker.covar.name]] <- dominance.component_(genoprobs.long = this.marker.genoprobs,
                                                          cross_type = class(cross)[1])
    }

  }

  return(df)
}






additive.component_ <- function(genoprobs.long,
                                cross_type = c('f2', 'bc')) {

  allele <- genoprob <- 'fake_global_for_CRAN'
  cross_type <- match.arg(arg = cross_type)

  alleles <- unique(genoprobs.long[['allele']])
  #
  # genoprobs.wide <- tidyr::spread(data = genoprobs.long,
  #                                 key = allele,
  #                                 value = genoprob)

  if (cross_type == 'f2') {

    # if we are on an autosome
    if (all(alleles %in% c('AA', 'AB', 'BB'))) {
      return(2*genoprobs.long$genoprob[genoprobs.long$allele == 'AA'] + genoprobs.long$genoprob[genoprobs.long$allele == 'AB'])
    } else if (all(alleles %in% c('g1', 'g2'))) {
      return(genoprobs.long$genoprob[genoprobs.long$allele == 'g2'])
    } else {
      stop(paste("Can't determine additive component of loc with alleles:", alleles))
    }
  }

  if (cross_type == 'bc') {
    if (all(alleles %in% c('AA', 'AB'))) {
      return(genoprobs.long$genoprob[genoprobs.long$allele == 'AB'])
    } else if (all(alleles %in% c('g1', 'g2'))) {
      return(genoprobs.long$genoprob[genoprobs.long$allele == 'g2'])
    } else {
      stop(paste("Can't determine additive component of loc with alleles:", alleles))
    }
  }

}

dominance.component_ <- function(genoprobs.long,
                                 cross_type = c('f2', 'bc')) {

  allele <- genoprob <- 'fake_global_for_CRAN'
  cross_type <- match.arg(arg = cross_type)

  alleles <- unique(genoprobs.long[['allele']])
  #
  # genoprobs.wide <- tidyr::spread(data = genoprobs.long,
  #                                 key = allele,
  #                                 value = genoprob)

  if (cross_type == 'f2') {
    if (all(alleles %in% c('AA', 'AB', 'BB'))) {
      return(genoprobs.long$genoprob[genoprobs.long$allele == 'AB'])
    } else if (all(alleles %in% c('g1', 'g2'))) {
      return(0)
    } else {
      stop(paste("Can't determine dominance component of loc with alleles:", alleles))
    }
  }

  if (cross_type == 'bc') {
    stop(paste("Can't use dominance component in a backcross."))
  }
}



make.genet.marker.model.df_ <- function(cross,
                                        marker.names) {

  stopifnot(is.cross(cross))

  # turn it into a tibble
  df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  for (marker.name in marker.names) {
    df[[marker.name]] <- qtl::pull.geno(cross = cross)[,marker.name]
  }

  df[['placeholder']] <- NULL

  return(df)
}


make.loc.specific.modeling.df <- function(general.modeling.df,
                                          loc.genoprobs,
                                          model.formulae,
                                          cross_type) {

  modeling.df <- general.modeling.df


  # note that we're not actually checking whether these components are used
  # it's OK if they get put in modeling.df and don't get used
  modeling.df[['mean.QTL.add']] <- modeling.df[['var.QTL.add']] <- additive.component_(genoprobs.long = loc.genoprobs, cross_type = cross_type)

  modeling.df[['mean.QTL.dom']] <- modeling.df[['var.QTL.dom']] <-  dominance.component_(genoprobs.long = loc.genoprobs, cross_type = cross_type)

  return(modeling.df)
}


