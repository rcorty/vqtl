#' @title mean_var_plot_model_free
#' @rdname mean_var_plots
#'
#' @param cross the cross object from which phenotypes and genotypes are drawn
#' @param phenotype.name the name of the phenotype
#' @param grouping.factor.names the factors to group individuals by
#'
#' @return Nothing, just plot.
#' @export
#'
mean_var_plot_model_free <- function(cross,
                                     phenotype.name,
                                     grouping.factor.names) {

  validate.mean_var_plot_model_free.input(cross, phenotype.name, grouping.factor.names)

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



validate.mean_var_plot_model_free.input <- function(cross,
                                                    phenotype.name,
                                                    grouping.factor.names) {

  stopifnot(is.cross(cross))

  stopifnot(phenotype.name %in% names(qtl::pull.pheno(cross)))

  for (gf.name in grouping.factor.names) {

    stopifnot(xor((gf.name %in% colnames(qtl::pull.geno(cross))),
                  (gf.name %in% names(qtl::pull.pheno(cross)))))
  }
}



#' @title mean_var_plot_model_based
#' @rdname mean_var_plots
#'
#' @param cross the cross
#' @param phenotype.name the name of the phenotype of interest
#' @param focal.groups the focal covariates, whose effects will be plotted.  Markers or phenotypes.
#' @param nuisance.groups the nuisance covariates, whose effects will be modeled, then marginalized over.  Markers or phenotypes.
#'
#' @return nothing, just the plot.
#' @export
#'
mean_var_plot_model_based <- function(cross,
                                      phenotype.name,
                                      focal.groups = NULL,
                                      nuisance.groups = NULL,
                                      genotype.names = c('AA', 'AB', 'BB'),
                                      xlim = NULL,
                                      ylim = NULL,
                                      title = paste(phenotype.name, 'by', paste(focal.groups, collapse = ', '))) {

  validate.mean_var_plot_model_based.input(cross = cross,
                                               phenotype.name = phenotype.name,
                                               focal.groups = focal.groups,
                                               nuisance.groups = nuisance.groups)

  # use modeling.df from sov, and add in phen names and marker names
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))

  modeling.df[[phenotype.name]] <- cross[['pheno']][[phenotype.name]]

  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))],
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))],
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])

  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[,marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[,phen.name])
  }

  modeling.df[['placeholder']] <- NULL

  # modeling.df <- dplyr::bind_cols(response.df,
  #                                 marker.df,
  #                                 phen.df)



  # make formulae from covariate names and cross
  # problem with nuisance is NULL -- todo
  covar.form <- paste(focal.groups, collapse = '+')
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, '+', paste(nuisance.groups, collapse = '+'))
  }
  mean.form <- paste(phenotype.name, '~', covar.form)
  var.form <- paste('~', covar.form)

  # pull formulae from sov and adapt
  # use null formulae bc they dont have QTL terms, and add in loci
  # mean.form <- deparse(sov[['meta']][['formulae']][['mean.null.formula']])
  # var.form <- deparse(sov[['meta']][['formulae']][['var.null.formula']])
  # for (marker.name in marker.names) {
  #   mean.form <- paste(mean.form, '+', marker.name)
  #   var.form <- paste(var.form, '+', marker.name)
  # }

  dglm.fit <- dglm::dglm(formula = formula(mean.form),
                         dformula = formula(var.form),
                         data = modeling.df)

  mean.pred <- predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit

  sd.pred <- predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit

  indiv.prediction.tbl <- dplyr::bind_cols(modeling.df,
                                           dplyr::data_frame(indiv.mean.estim = mean.estim,
                                                             indiv.mean.lb = mean.estim - mean.se,
                                                             indiv.mean.ub = mean.estim + mean.se,
                                                             indiv.sd.estim = exp(sd.estim),
                                                             indiv.sd.lb = exp(sd.estim - sd.se),
                                                             indiv.sd.ub = exp(sd.estim + sd.se)))

  group.prediction.tbl <- indiv.prediction.tbl %>%
    dplyr::group_by_(.dots = c(focal.groups)) %>%
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim),
                     group.mean.lb = mean(indiv.mean.lb),
                     group.mean.ub = mean(indiv.mean.ub),
                     group.sd.estim = mean(indiv.sd.estim),
                     group.sd.lb = mean(indiv.sd.lb),
                     group.sd.ub = mean(indiv.sd.ub))


  p <- ggplot2::ggplot(data = group.prediction.tbl,
                       mapping = ggplot2::aes_string(color = focal.groups[1])) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = group.mean.estim, y = group.sd.estim)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = group.mean.lb, xend = group.mean.ub, y = group.sd.estim, yend = group.sd.estim)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = group.mean.estim, xend = group.mean.estim, y = group.sd.lb, yend = group.sd.ub)) +
    ggplot2::xlab('mean estimate +/- 1 SE') +
    ggplot2::ylab('SD estimate +/- 1 SE') +
    ggplot2::ggtitle(title)

  if (length(focal.groups) > 1) {
    p <- p + ggplot2::geom_point(mapping = ggplot2::aes_string(x = 'group.mean.estim',
                                                               y = 'group.sd.estim',
                                                               shape = focal.groups[2]),
                                 size = 3)
  }

  if (!is.null(xlim) & !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  if (!is.null(xlim) & is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  if (is.null(xlim) & !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }


  return(p)
}


validate.mean_var_plot_model_based.input <- function(cross,
                                                     phenotype.name,
                                                     focal.groups,
                                                     nuisance.groups) {

  return(TRUE)
}



