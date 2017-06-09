#' @title effects_over_genome_plot
#' @rdname effects_over_genome_plot
#'
#' @param sov the scanonevar
#' @param effect.names the effects to plot
#' @param mean.or.var mean or variance effects of covariates
#' @param se_ribbons Should a ribbon from estimate - se to estimate + se be plotted?
#'
#' @description Plots estimated effects and their standard errors at each locus in the genome.
#'
#' @return the plot
#' @export
#'
effects_over_genome_plot <- function(sov,
                                     effect.names = NULL,
                                     mean.or.var = c('both', 'mean', 'var'),
                                     se_ribbons = TRUE) {

  effect_name <- estimate <- chr.type <- pos <- se <- effect_type <- 'fake_global_for_CRAN'

  mean.or.var = match.arg(arg = mean.or.var)

  validate_effects_over_genome_plot_inputs(sov, effect.names, mean.or.var)

  wrangle_effects_over_genome_plot_inputs(sov, effect.names, mean.or.var)

  # get rid of association statistics
  result <- sov[['result']] %>%
    dplyr::select(-dplyr::matches('lod')) %>%
    dplyr::select(-dplyr::matches('asymp.p')) %>%
    dplyr::select(-dplyr::matches('empir.p'))

  ests <- result %>%
    dplyr::select(-dplyr::matches('mse|vse')) %>%
    tidyr::gather(key = effect_name, value = estimate, -(chr.type:pos)) %>%
    dplyr::mutate(effect_type = factor(x = ifelse(test = grepl(pattern = 'mef',
                                                               x = effect_name),
                                                  yes = 'mean',
                                                  no = ifelse(test = grepl(pattern = 'vef',
                                                                           x = effect_name),
                                                              yes = 'var',
                                                              no = 'unknown type')),
                                       levels = c('mean', 'var'))) %>%
    dplyr::mutate(effect_name = gsub(pattern = 'mef_', replacement = '', x = effect_name),
                  effect_name = gsub(pattern = 'vef_', replacement = '', x = effect_name))



  ses <- result %>%
    dplyr::select(-dplyr::matches('mef|vef')) %>%
    tidyr::gather(key = effect_name, value = se, -(chr.type:pos)) %>%
    dplyr::mutate(effect_type = factor(x = ifelse(test = grepl(pattern = 'mse',
                                                               x = effect_name),
                                                  yes = 'mean',
                                                  no = ifelse(test = grepl(pattern = 'vse',
                                                                           x = effect_name),
                                                              yes = 'var',
                                                              no = 'unknown type')),
                                       levels = c('mean', 'var'))) %>%
    dplyr::mutate(effect_name = gsub(pattern = 'mse_', replacement = '', x = effect_name),
                  effect_name = gsub(pattern = 'vse_', replacement = '', x = effect_name))

  to_plot <- dplyr::inner_join(x = ests, y = ses, by = c("chr.type", "chr", "loc.name", "pos", "effect_name", "effect_type"))


  # sort levels of effect_name so that mean effects come first
  # dplyr::mutate(effect_name = factor(x = effect_name),
  #               effect_name = factor(x = effect_name,
  #                                    levels = c(grep(pattern = 'mef', x = levels(effect_name), value = TRUE),
  #                                               grep(pattern = 'vef', x = levels(effect_name), value = TRUE)))) %>%

  if (!is.null(effect.names)) {
    to_plot <- to_plot %>%
      dplyr::filter(grepl(pattern = paste0(effect.names, collapse = '|'), x = effect_name))
  }
  to_plot <- to_plot %>%
    dplyr::mutate(effect_name = gsub(pattern = 'mean.QTL.', replacement = 'QTL.', x = effect_name),
                  effect_name = gsub(pattern = 'var.QTL.', replacement = 'QTL.', x = effect_name),
                  effect_name = factor(effect_name))

  if (!(mean.or.var == 'both')) {
    to_plot <- to_plot %>%
      dplyr::filter(effect_type == mean.or.var)
  }

  # make sure no crap got through
  stopifnot(!any(to_plot$effect_type == 'unknown type'))

  ggplot2::ggplot(data = to_plot,
                  mapping = ggplot2::aes(x = pos, y = estimate, group = interaction(effect_name, effect_type))) +
    ggplot2::geom_abline(slope = 0, intercept = 0) +
    ggplot2::facet_grid(~chr,
                        scales = 'free_x', space = 'free_x', switch = 'x') +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = estimate - se, ymax = estimate + se, fill = effect_name), alpha = 0.3) +
    ggplot2::geom_line(mapping = ggplot2::aes(color = effect_name, linetype = effect_type), size = 1) +
    ggplot2::scale_color_brewer(type = 'qual', palette = 2) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = 2) +
    ggplot2::scale_linetype_manual(breaks = c('mean', 'var'), values = c(1, 2), drop = FALSE)
}


validate_effects_over_genome_plot_inputs <- function(sov, effect.names, mean.or.var) {

  return(TRUE)
}



wrangle_effects_over_genome_plot_inputs <- function(sov, effect.names, mean.or.var) {

  return(TRUE)
}