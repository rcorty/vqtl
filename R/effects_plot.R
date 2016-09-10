#' @title effects_plot
#' @rdname effects_plot
#'
#' @param sov the scanonevar
#' @param effect.names the effects to plot
#' @param mean.or.var mean or variance effects of covariates
#'
#' @return the plot
#' @export
#'
effects_plot <- function(sov,
                         effect.names = NULL,
                         mean.or.var = c('both', 'mean', 'var')) {

  mean.or.var = match.arg(mean.or.var)

  validate_effects_plot_inputs(sov, effect.names, mean.or.var)

  wrangle_effects_plot_inputs(sov, effect.names, mean.or.var)

  # get rid of association statistics
  result <- sov[['result']] %>%
    dplyr::select(-dplyr::matches('lod')) %>%
    dplyr::select(-dplyr::matches('asymp.p')) %>%
    dplyr::select(-dplyr::matches('empir.p'))

  to.plot <- result %>%
    tidyr::gather(key = effect_name, value = effect_estimate, -(chr.type:pos))


    # sort levels of effect_name so that mean effects come first
    # dplyr::mutate(effect_name = factor(x = effect_name),
    #               effect_name = factor(x = effect_name,
    #                                    levels = c(grep(pattern = 'mef', x = levels(effect_name), value = TRUE),
    #                                               grep(pattern = 'vef', x = levels(effect_name), value = TRUE)))) %>%

  if (!is.null(effect.names)) {
    to.plot <- to.plot %>%
      dplyr::filter(grepl(pattern = effect.names, x = effect_name))
  }

    # generate effect_type based on whether effect_name has 'mef' or 'vef'
  to.plot <- to.plot %>%
    dplyr::mutate(effect_type = factor(x = ifelse(test = grepl(pattern = 'mef',
                                                               x = effect_name),
                                                  yes = 'mean',
                                                  no = ifelse(test = grepl(pattern = 'vef',
                                                                           x = effect_name),
                                                              yes = 'var',
                                                              no = 'unknown type')),
                                       levels = c('mean', 'var')))

  if (!(mean.or.var == 'both')) {
    to.plot <- to.plot %>%
      dplyr::filter(effect_type == mean.or.var)
  }

  # remove information that is now in effect_type from effect_name
  to.plot <- to.plot %>%
    dplyr::mutate(effect_name = gsub(pattern = '_mef', replacement = '', x = effect_name),
                  effect_name = gsub(pattern = '_vef', replacement = '', x = effect_name),
                  effect_name = gsub(pattern = 'mean.QTL.', replacement = 'QTL.', x = effect_name),
                  effect_name = gsub(pattern = 'var.QTL.', replacement = 'QTL.', x = effect_name),
                  effect_name = factor(effect_name))

  stopifnot(!any(to.plot$effect_type == 'unknown type'))

  ggplot2::ggplot(data = to.plot,
                  mapping = ggplot2::aes(x = pos, y = effect_estimate)) +
    ggplot2::facet_grid(~chr,
                        scales = 'free_x', space = 'free_x', switch = 'x') +
    ggplot2::geom_line(mapping = ggplot2::aes(color = effect_name, linetype = effect_type), size = 1) +
    ggplot2::scale_color_brewer(type = 'qual', palette = 2) +
    ggplot2::scale_linetype_manual(breaks = c('mean', 'var'), values = c(1, 2), drop = FALSE)
}


validate_effects_plot_inputs <- function(sov, effect.names, mean.or.var) {

  return(TRUE)
}



wrangle_effects_plot_inputs <- function(sov, effect.names, mean.or.var) {

  return(TRUE)
}