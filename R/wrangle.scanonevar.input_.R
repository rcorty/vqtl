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
#'
wrangle.scanonevar.input_ <- function(cross,
                                      mean.formula,
                                      var.formula) {


  if (!is.cross.w.genoprobs(x = cross)) {
    message('calculating genoprobs with stepwidth = 2')
    cross <- qtl::calc.genoprob(cross = cross, step = 2)
  }

  # build scan.df
  scan.df <- dplyr::bind_rows(lapply(X = cross[['geno']],
                                     FUN = function(x) as.data.frame.table(x[['prob']], stringsAsFactors = FALSE)))
  names(scan.df) <- c('iid', 'marker.name', 'allele', 'genoprob')
  iids <- unique(scan.df[['iid']])
  marker.names <- unique(scan.df[['marker.name']])

  # todo: return this wider -- AA, AB, and BB should be columns


  response.name <- as.character(mean.formula[[2]])
  mean.covar.names <- labels(object = terms(x = mean.formula))
  var.covar.names <- labels(object = terms(x = var.formula))

  # start building mapping.df by pulling all necessary columns from cross[['pheno']]
  mapping.df <- dplyr::select(.data = dplyr::tbl_df(cross[['pheno']]),
                              dplyr::one_of(c(response.name,
                                              mean.covar.names,
                                              var.covar.names)))


  # continue building mapping.df by adding keyword columns as NA
  mean.keywords <- c('mean.QTL.add', 'mean.QTL.dom')
  for (mean.keyword in mean.keywords) {
    if (mean.keyword %in% mean.covar.names)
      mapping.df[[mean.keyword]] <- NA
  }

  var.keywords <- c('var.QTL.add', 'var.QTL.dom')
  for (var.keyword in var.keywords) {
    if (var.keyword %in% var.covar.names)
      mapping.df[[var.keyword]] <- NA
  }

  # continue building mapping.df by adding genetic markers
  marker.names <- unlist(lapply(X = cross[['geno']],
                                FUN = function(chr) names(chr[['map']])))



  return(list(modeling.df = modeling.df,
              marker.info.df = marker.info.df,
              genoprob.df = genoprob.df,
              scan.types = scan.types,
              scan.formulae = scan.formulae))


}