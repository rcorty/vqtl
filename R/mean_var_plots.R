#' @title sample_stats_mean_var_plot
#' @rdname mean_var_plots
#'
#' @param cross the cross object from which phenotypes and genotypes are drawn
#' @param phenotype.name the name of the phenotype
#' @param grouping.factor.names the factors to group individuals by
#'
#' @return Nothing, just plot.
#' @export
#'
sample_stats_mean_var_plot <- function(cross,
                                       phenotype.name,
                                       grouping.factor.names) {

  validate.mean.var.sample.plot.input(cross, phenotype.name, grouping.factor.names)

  marker.names <- grouping.factor.names[grouping.factor.names %in% colnames(qtl::pull.geno(cross))]
  phen.names <- grouping.factor.names[grouping.factor.names %in% names(qtl::pull.pheno(cross))]

  plotting.df <- dplyr::bind_cols(make.response.model.df_(cross = cross,
                                                          response.name = phenotype.name),
                                  make.genet.marker.model.df_(cross = cross,
                                                              marker.names = marker.names),
                                  make.phen.covar.model.df_(cross = cross,
                                                            phen.names = phen.names))

  for (gf.name in grouping.factor.names) {
    plotting.df[[gf.name]] <- factor(plotting.df[[gf.name]])
  }

  p <- plotting.df %>%
    dplyr::group_by_(.dots = grouping.factor.names) %>%
    dplyr::summarise_(mean = lazyeval::interp(~mean(var, na.rm = TRUE), var = as.name(phenotype.name)),
                      sd = lazyeval::interp(~sd(var, na.rm = TRUE), var = as.name(phenotype.name)),
                      mean.se = quote(sd/sqrt(n())),
                      sd.se = quote(sqrt(2)*sd^2/sqrt(n() - 1))) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = mean, y = sd)) +
    ggplot2::geom_segment(mapping = ggplot2::aes_string(x = 'mean - mean.se', xend = 'mean + mean.se', yend = 'sd', color = grouping.factor.names[1])) +
    ggplot2::geom_segment(mapping = ggplot2::aes_string(y = 'sd - sd.se', yend = 'sd + sd.se', xend = 'mean', color = grouping.factor.names[1])) +
    ggplot2::geom_path(size = 3, alpha = 0.5, mapping = ggplot2::aes_string(color = grouping.factor.names[1]))

  if (length(grouping.factor.names) == 1) {
    p <- p + ggplot2::geom_point(size = 3, ggplot2::aes_string(color = grouping.factor.names[1]))
  } else {
    p <- p + ggplot2::geom_point(size = 3, ggplot2::aes_string(color = grouping.factor.names[1], shape = grouping.factor.names[2]))
  }

  return(p)
}



validate.sample_stats_mean_var_plot.input <- function(cross,
                                                      phenotype.name,
                                                      grouping.factor.names) {

  stopifnot(is.cross(cross))

  stopifnot(phenotype.name %in% names(qtl::pull.pheno(cross)))

  for (gf.name in grouping.factor.names) {

    stopifnot(xor((gf.name %in% colnames(qtl::pull.geno(cross))),
                  (gf.name %in% names(qtl::pull.pheno(cross)))))
  }
}



#' @title modeled_effects_mean_var_plot
#' @rdname mean_var_plots
#'
#' @param cross the cross
#' @param phenotype.name the name of the phenotype of interest
#' @param focal.covariate.names the focal covariates, whose effects will be plotted.  Markers or phenotypes.
#' @param nuisance.covariate.names the nuisance covariates, whose effects will be modeled, then marginalized over.  Markers or phenotypes.
#'
#' @return nothing, just the plot.
#' @export
#'
modeled_effects_mean_var_plot <- function(cross,
                                          phenotype.name,
                                          focal.covariate.names,
                                          nuisance.covariate.names) {

  validate.modeled_effects_mean_var_plot.input(cross = cross,
                                               phenotype.name = phenotype.name,
                                               focal.covariate.names = focal.covariate.names,
                                               nuisance.covariate.names = nuisance.covariate.names)

  # use modeling.df from sov, and add in phen names and marker names
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))

  marker.names <- c(focal.covariate.names[focal.covariate.names %in% colnames(qtl::pull.geno(cross = cross))],
                    nuisance.covariate.names[nuisance.covariate.names %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.covariate.names[focal.covariate.names %in% colnames(qtl::pull.pheno(cross = cross))],
                  nuisance.covariate.names[nuisance.covariate.names %in% colnames(qtl::pull.pheno(cross = cross))])

  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(qtl::pull.geno(cross = cross)[,marker.name])
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[,phen.name])
  }

  modeling.df[[placeholder]] <- NULL

  modeling.df <- dplyr::bind_cols(response.df,
                                  marker.df,
                                  phen.df)



  # use null formulae bc they dont have QTL terms, and add in loci
  mean.form <- deparse(sov[['meta']][['formulae']][['mean.null.formula']])
  var.form <- deparse(sov[['meta']][['formulae']][['var.null.formula']])
  for (marker.name in marker.names) {
    mean.form <- paste(mean.form, '+', marker.name)
    var.form <- paste(var.form, '+', marker.name)
  }

  dglm.fit <- dglm::dglm(formula = formula(mean.form),
                    dformula = formula(var.form),
                    data = modeling.df)

  mean.pred <- predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit

  var.pred <- predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  var.estim <- var.pred$fit/var.pred$residual.scale
  var.se <- var.pred$se.fit/var.pred$residual.scale

  prediction.tbl <- dplyr::bind_cols(modeling.df,
                                     data_frame(indiv.mean.estim = mean.estim,
                                                indiv.mean.lb = mean.estim - mean.se,
                                                indiv.mean.ub = mean.estim + mean.se,
                                                indiv.var.estim = exp(var.estim),
                                                indiv.var.lb = exp(var.estim - var.se),
                                                indiv.var.ub = exp(var.estim + var.se)))

  prediction.tbl %>%
    dplyr::group_by_(.dots = grouping.factor.names) %>%



  return(3)
}


validate.modeled_effects_mean_var_plot.input <- function(cross,
                                                         phenotype.name,
                                                         focal.covariate.names,
                                                         nuisance.covariates.names) {

  return(TRUE)
}



