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
#'
#' @return 27599
#' @export
#'
#' @examples
#' x <- 27599
#'
#'
scanonevar <- function(cross,
                       mean.formula = phenotype ~ mean.QTL.add + mean.QTL.dom,
                       var.formula = ~ var.QTL.add + var.QTL.dom,
                       chrs = qtl::chrnames(cross = cross)) {

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

  # execute the scan(s)
  scanonevar_(modeling.df = wrangled.inputs$modeling.df,
              loc.info.df = wrangled.inputs$loc.info.df,
              genoprob.df = wrangled.inputs$genoprob.df,
              scan.types = wrangled.inputs$scan.types,
              scan.formulae = wrangled.inputs$scan.formulae)
}
