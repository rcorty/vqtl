#' @title wrangle.scanonevar.input_
#' @name wrangle.scanonevar.input_
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams scanonevar
#'
#' @return a list with three elements: 'model.df', the data.frame to be used
#' at every genetic locus in the upcoming scanonevar_, 'scan.df', the
#' data.frame that will be picked from, one genetic at a time to fill the
#' 'QTL' variable in the model and 'scan.types', which indicates which types
#' of scan to do: mqtl, vqtl, and/or mvqtl.
#'
#' @examples
#' x <- 27599
#' @export
wrangle.scanonevar.input_ <- function(cross,
                                      mean.formula,
                                      var.formula,
                                      chrs) {

  if (!is.cross.w.genoprobs(x = cross)) {
    message("calculating genoprobs with stepwidth = 2, off.end = 0, error.prob = 1e-4, map.function = 'haldane'")
    cross <- qtl::calc.genoprob(cross = cross, step = 2)
  }

  # qtl::pull.geno and qtl::pull.markers don't provide information on pseudomarkers
  # qtl::pull.genoprob doesn't give the position of the pseudomarkers
  loc.info.df <- wrangle.loc.info.df_(cross = cross,
                                      chrs = chrs)

  # could probably reshape the output from qtl::pull.genoprob to accomplish the same thing
  genoprob.df.long <- wrangle.genoprob.df_(cross = cross)

  scan.types <- wrangle.scan.types_(mean.formula, var.formula)

  scan.formulae <- wrangle.scan.formulae_(mean.formula, var.formula)

  modeling.df <- wrangle.modeling.df_(cross = cross,
                                      genoprobs = genoprob.df.long,
                                      scan.formulae = scan.formulae)

  return(list(scan.types = scan.types,
              scan.formulae = scan.formulae,
              modeling.df = modeling.df,
              loc.info.df = loc.info.df,
              genoprob.df = genoprob.df.long))
}





#' @title wrangle.loc.info.df_
#' @rdname wrangle.scanonevar.input_
#'
#' @inheritParams wrangle.scanonevar.input_
#' @export
wrangle.loc.info.df_ <- function (cross, chrs) {

  loc.info.from.chr <- function(x) {
    chr.name <- names(class(x))
    prob.map <- attr(x = x[['prob']], which = 'map')
    names.starting.with.loc.idxs <- grep(pattern = '^loc', names(prob.map))
    names(prob.map)[names.starting.with.loc.idxs] <-
      paste0('chr', chr.name, '_', names(prob.map)[names.starting.with.loc.idxs])

    return(dplyr::data_frame(loc.name = names(prob.map),
                             chr = chr.name,
                             pos = prob.map))
  }
  # NB: names(class(x)) is the notation to get the name of a chromosome
  cross[['geno']] <- cross[['geno']][qtl::chrnames(cross = cross) %in% chrs]
  return(dplyr::bind_rows(lapply(X = cross[['geno']],
                                 FUN = loc.info.from.chr)))
}


#' @title wrangle.genoprob.df_
#' @rdname wrangle.scanonevar.input_
#'
#' @inheritParams wrangle.scanonevar.input_
#' @export
wrangle.genoprob.df_ <- function(cross) {

  genoprobs.from.chr <- function(x) {
    prob.tbl <- x[['prob']]
    num.width <- max(nchar(as.character(1:dim(prob.tbl)[1])))
    if (is.null(dimnames(prob.tbl)[[1]])) {
      dimnames(prob.tbl)[[1]] <- paste0('org',
                                        stringr::str_pad(string = 1:dim(prob.tbl)[1],
                                                         width = num.width))
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


#' @title wrangle.scan.types_
#' @rdname wrangle.scanonevar.input_
#'
#' @inheritParams wrangle.scanonevar.input_
#' @export
wrangle.scan.types_ <- function(mean.formula,
                                var.formula) {

  mean.qtl.idxs <- grep(pattern = 'mean.QTL', x = labels(terms(mean.formula)))
  var.qtl.idxs <- grep(pattern = 'var.QTL', x = labels(terms(var.formula)))

  if (all(mean.qtl.idxs, !var.qtl.idxs))
    return('mean')
  if (all(!mean.qtl.idxs, var.qtl.idxs))
    return('var')
  if (all(mean.qtl.idxs, var.qtl.idxs))
    return(c('mean', 'var', 'joint'))

  stop('Should never get here.')
}


#' @title wrangle.scan.formulae_
#' @rdname wrangle.scanonevar.input_
#'
#' @inheritParams wrangle.scanonevar.input_
#' @export
wrangle.scan.formulae_ <- function(mean.formula,
                                   var.formula) {

  mean.terms <- labels(terms(mean.formula))
  var.terms <- labels(terms(var.formula))

  # second, identify terms that are 'keywords' and remove them to make the null formulae
  mean.qtl.idxs <- grep(pattern = 'mean.QTL', x = mean.terms)
  var.qtl.idxs <- grep(pattern = 'var.QTL', x = var.terms)

  # if no qtl terms, 'mean.null' is NULL and no mean testing will be done
  # if no non-qtl terms, rhs is just 1
  if (length(mean.qtl.idxs) == 0) {
    mean.null.formula <- NULL
  } else if (length(mean.qtl.idxs) == length(labels(terms(mean.formula)))) {
    mean.null.formula <- reformulate(termlabels = '1', response = mean.formula[[2]])
  } else {
    mean.null.formula <- formula(drop.terms(termobj = terms(mean.formula),
                                            dropx = mean.qtl.idxs,
                                            keep.response = TRUE))
  }

  # if no qtl terms, 'var.null' is NULL and no var testing will be done
  # if no non-qtl terms, rhs is just 1
  if (length(var.qtl.idxs) == 0) {
    var.null.formula <- NULL
  } else if (length(var.qtl.idxs) == length(labels(terms(var.formula)))) {
    var.null.formula <- reformulate(termlabels = '1', response = NULL)
  } else {
    var.null.formula <- formula(drop.terms(termobj = terms(var.formula),
                                            dropx = var.qtl.idxs,
                                            keep.response = FALSE))
  }

  # this way of adding to a list doesn't add anything when RHS is NULL
  scan.formulae <- list()
  scan.formulae[['mean.alt.formula']] <- mean.formula
  scan.formulae[['mean.null.formula']] <-  mean.null.formula
  scan.formulae[['var.alt.formula']] <- var.formula
  scan.formulae[['var.null.formula']] <-  var.null.formula

  return(scan.formulae)
}


#' @title wrangle.modeling.df_
#' @rdname wrangle.scanonevar.input_
#'
#' @inheritParams wrangle.scanonevar.input_
#' @export
wrangle.modeling.df_ <- function(cross,
                                 scan.formulae,
                                 genoprobs) {

  # start making modeling.df with response
  response.name <- as.character(scan.formulae[['mean.alt.formula']][[2]])
  modeling.df <- dplyr::select_(.data = dplyr::tbl_df(cross[['pheno']]),
                                response.name)

  # get all covariate names
  mean.covar.names <- labels(object = terms(x = scan.formulae[['mean.alt.formula']]))
  var.covar.names <- labels(object = terms(x = scan.formulae[['var.alt.formula']]))


  # get the covariate names that are keywords and add them as NA columns to modeling.df
  mean.keywords <- c('mean.QTL.add', 'mean.QTL.dom')
  mean.keyword.covar.names <- mean.covar.names[mean.covar.names %in% mean.keywords]
  var.keywords <- c('var.QTL.add', 'var.QTL.dom')
  var.keyword.covar.names <- var.covar.names[var.covar.names %in% var.keywords]
  for (keyword in c(mean.keyword.covar.names, var.keyword.covar.names)) {
    modeling.df[[keyword]] <- NA
  }


  # get the covariate names that are phenotypes and add them to modeling.df
  mean.phen.covar.names <- mean.covar.names[mean.covar.names %in% names(cross[['pheno']])]
  var.phen.covar.names <- var.covar.names[var.covar.names %in% names(cross[['pheno']])]
  for (phen.covar.name in unique(c(mean.phen.covar.names, var.phen.covar.names))) {
    modeling.df[[phen.covar.name]] <- cross[['pheno']][[phen.covar.name]]
  }


  # get the covariate names that are markers and add them to modeling.df
  marker.names <- colnames(qtl::pull.geno(cross))
  add.marker.names <- paste0(marker.names, '_add')
  dom.marker.names <- paste0(marker.names, '_dom')
  mean.add.marker.covar.names <- mean.covar.names[mean.covar.names %in% add.marker.names]
  mean.dom.marker.covar.names <- mean.covar.names[mean.covar.names %in% dom.marker.names]
  var.add.marker.covar.names <- var.covar.names[var.covar.names %in% add.marker.names]
  var.dom.marker.covar.names <- var.covar.names[var.covar.names %in% dom.marker.names]

  for (add.marker.covar.name in unique(c(mean.add.marker.covar.names, var.add.marker.covar.names))) {
    this.marker.genoprobs <- dplyr::filter(.data = genoprobs,
                                           loc.name == substr(x = add.marker.covar.name,
                                                              start = 1,
                                                              stop = nchar(add.marker.covar.name) - 4))
    modeling.df[[add.marker.covar.name]] <- additive.component(this.marker.genoprobs)
  }


  for (dom.marker.covar.name in unique(c(mean.dom.marker.covar.names, var.dom.marker.covar.names))) {
    this.marker.genoprobs <- dplyr::filter(.data = genoprobs,
                                           loc.name == substr(x = dom.marker.covar.name,
                                                              start = 1,
                                                              stop = nchar(dom.marker.covar.name) - 4))
    modeling.df[[dom.marker.covar.name]] <- dominance.component(this.marker.genoprobs)
  }

  return(modeling.df)
}