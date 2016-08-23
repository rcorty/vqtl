#' #' @title Plot Predictive Interval for Categorical Genotype/Phenotype Groups
#' #' @name predictive.plot
#' #' @author Robert Corty \email{rcorty@@gmail.com}
#' #'
#' #' @description \code{predictive.plot} should be used to visually investigate loci identified
#' #'   with plot.scanonevar or summary.scanonevar.  The user can specify the same mean and variance
#' #'   formulae that were used in the scan, or specify new formulae to investigate interactions.
#' #'
#' #' @param cross The cross object to be plotted
#' #' @param mean.formula the mean formula
#' #' @param var.formula the var formula
#' #' @param title Optionally, title for the plot.  Defaults to 'Predictive of [response phenotype] from
#' #'   [predictive phenotype (e.g. sex)] and [marker name]
#' #' @param xlim Optionally specify x-axis limits.  Defaults to data-dependent.
#' #' @param ylim Optionally specify y-axis limits.  Defaults to data.dependent.
#' #'
#' #' @return The plot
#' #'
#' #' @export
#' #'
#' #' @examples
#' #'   set.seed(27599)
#' #'   my.cross <- sim.cross(map = sim.map(), type = 'f2')
#' #'   my.cross <- calc.genoprob(my.cross)
#' #'
#' prediction.ci.plot <- function(cross,
#'                                mean.formula,
#'                                var.formula,
#'                                group.by,
#'                                title = paste('Predictive Mean and Variance Confidence Intervals of', response.phen),
#'                                xlim = NA,
#'                                ylim = NA) {
#'
#'   # resuse validation from scanonevar, with 'allow.no.qtl' = TRUE (default is FALSE)
#'   validate.scanonevar.input_(cross = cross,
#'                              mean.formula = mean.formula,
#'                              var.formula = var.formula,
#'                              chrs = qtl::chrnames(cross),
#'                              allow.no.qtl = TRUE)
#'
#'   # reuse data wrangling from scanonevar, with chrs = NULL (suppresses df that forms skeleton of mapping result)
#'   wrangled.inputs <- wrangle.prediction.ci.plot.input_(cross = cross,
#'                                                        mean.formula = mean.formula,
#'                                                        var.formula = var.formula,
#'                                                        chrs = NULL)
#'
#'   dglm.fit <- dglm::dglm(formula = wrangled.inputs$scan.formulae[['mean.alt.formula']],
#'                          dformula = wrangled.inputs$scan.formulae[['var.alt.formula']],
#'                          data = wrangled.inputs$modeling.df)
#'
#'   mean.predictions <- predict(dglm.fit, se.fit = TRUE)
#'   mean.estimates <- mean.predictions$fit
#'   mean.standard.errors <- mean.predictions$se.fit
#'
#'   var.predictions <- predict(dglm.fit$dispersion.fit, se.fit = TRUE)
#'   var.estimates <- var.predictions$fit
#'   var.standard.errors <- var.predictions$se.fit
#'
#'   indiv.predictions <- dplyr::bind_cols(wrangled.inputs$modeling.df,
#'                                         dplyr::data_frame(indiv.mean.estim = mean.estimates,
#'                                                           indiv.mean.lb = mean.estimates - mean.standard.errors,
#'                                                           indiv.mean.ub = mean.estimates + mean.standard.errors,
#'                                                           indiv.var.estim = exp(var.estimates),
#'                                                           indiv.var.lb = exp(var.estimates - var.standard.errors),
#'                                                           indiv.var.ub = exp(var.estimates + var.standard.errors)))
#'
#'
#'   groupl.predictions <- indiv.predictions %>%
#'     group_by_(group.by) %>%
#'     summarise(group.mean.estim = mean(indiv.mean.estim),
#'               group.mean.lb = mean(indiv.mean.lb),
#'               group.mean.ub = mean(indiv.mean.ub),
#'               group.var.estim = mean(indiv.var.estim),
#'               group.var.lb = mean(indiv.var.lb),
#'               group.var.ub = mean(indiv.var.ub)) %>%
#'     arrange(phen, genotype)
#'
#'   x <- 3
#'
#' }
