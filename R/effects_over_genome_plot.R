#' @title effects_over_genome_plot
#' @rdname effects_over_genome_plot
#'
#' @param sov the scanonevar
#' @param se_ribbons Should a ribbon from estimate - se to estimate + se be plotted?
#' @param effect_type_regex regex that matches 'mean', 'var', or both
#' @param covar_name_regex regex that matches the covars we want to plot
#' @param transform_var_effects combine variance effects w intercept and exponentiate?
#'
#' @description Plots estimated effects and their standard errors at each locus in the genome.
#'
#' @return the plot
#' @export
#'
effects_over_genome_plot <- function(sov,
                                     covar_name_regex = '.',
                                     effect_type_regex = '(mean|var)',
                                     transform_var_effects = TRUE,
                                     se_ribbons = TRUE) {

  chr <- chr.type <- covar_name <- covar_name_other <- effect_name <- 'fake_global_for_CRAN'
  effect_type <- est <- est_int  <- est_or_se <- est_other <- 'fake_global_for_CRAN'
  estimate <- lb <- loc.name <- pos <- se <- se_int <- se_other <- ub <- 'fake_global_for_CRAN'

  # get rid of association statistics
  result <- sov[['result']] %>%
    dplyr::select(-dplyr::matches('lod')) %>%
    dplyr::select(-dplyr::matches('asymp.p')) %>%
    dplyr::select(-dplyr::matches('empir.p'))

  long_result <- result %>%
    tidyr::gather(key = effect_name, value = estimate, -(chr.type:pos)) %>%
    tidyr::separate(col = effect_name, into = c('effect_type', 'est_or_se', 'covar_name'), sep = '__') %>%
    dplyr::mutate(effect_type = factor(x = effect_type, levels = c('m', 'v'), labels = c('mean', 'var'))) %>%
    dplyr::mutate(covar_name = gsub(pattern = '(mean|var)\\.QTL', replacement = 'QTL', x = covar_name)) %>%
    tidyr::spread(key = est_or_se, value = estimate) %>%
    dplyr::filter(grepl(pattern = effect_type_regex, x = effect_type))


  mean_to_plot <- long_result %>%
    dplyr::filter(effect_type == 'mean') %>%
    dplyr::mutate(lb = est - se,
                  ub = est + se) %>%
    dplyr::select(-se)

  if (transform_var_effects) {

    intercepts <- long_result %>% dplyr::filter(effect_type == 'var', covar_name == '(Intercept)')
    the_rest <- long_result %>% dplyr::filter(effect_type == 'var', covar_name != '(Intercept)')

    var_to_plot <- dplyr::inner_join(x = intercepts, y = the_rest,
                                     by = c('chr.type', 'chr', 'loc.name', 'pos', 'effect_type'),
                                     suffix = c('_int', '_other')) %>%
      dplyr::mutate(est = exp(est_int + est_other),
                    lb = exp(est_int - se_int + est_other - se_other),
                    ub = exp(est_int + se_int + est_other + se_other)) %>%
      dplyr::select(chr.type, chr, loc.name, pos, effect_type, covar_name_other, est, lb, ub) %>%
      dplyr::rename(covar_name = covar_name_other)

  } else {

    var_to_plot <- long_result %>%
      dplyr::filter(effect_type == 'var') %>%
      dplyr::mutate(lb = est - se,
                    ub = est + se) %>%
      dplyr::select(-se)

  }

  to_plot <- dplyr::bind_rows(mean_to_plot, var_to_plot) %>%
    dplyr::filter(grepl(pattern = covar_name_regex, x = covar_name))

  p <- ggplot2::ggplot(data = to_plot,
                  mapping = ggplot2::aes(x = pos, group = interaction(covar_name, effect_type))) +
    ggplot2::geom_abline(slope = 0, intercept = 0) +
    ggplot2::facet_grid(~chr,
                        scales = 'free_x', space = 'free_x', switch = 'x')

  if (se_ribbons) {
    p <- p +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = lb, ymax = ub, fill = covar_name), alpha = 0.3)
  }

  p <- p +
    ggplot2::geom_line(mapping = ggplot2::aes(y = est, color = covar_name, linetype = effect_type), size = 1) +
    ggplot2::scale_color_brewer(type = 'qual', palette = 2) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = 2) +
    ggplot2::scale_linetype_manual(breaks = c('mean', 'var'), values = c(1, 2), drop = FALSE)

  return(p)
}
