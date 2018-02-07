#' @title phenotype_at_marker_plot
#' @rdname plotting
#'
#' @description plotting functions for package vqtl
#'
#' @param phenotype_name The phenotype to plot
#' @param marker_name The marker to stratify observations by
#' @param color_by variable name to color the points by
#' @param point_alpha alpha value (see-throughness) of the plotted points
#' @param Ibars Should I bars be plotted showing the standard deviation of each group?
#' @param connectIbars Should the Ibars be connected horizontally?
#' @param genotype_labels plotting labels for genotype groups
#' @param shape_by a discrete phenotype to map to the shape aesthetic of the points
#'
#' @return nothing.  Just plots.
#' @export
#'
phenotype_at_marker_plot <- function(cross,
                                     phenotype_name,
                                     marker_name,
                                     color_by = NULL,
                                     shape_by = NULL,
                                     point_alpha = 1,
                                     point_size = 1,
                                     Ibars = TRUE,
                                     connectIbars = TRUE,
                                     genotype_labels = NULL) {

  sd <- ub <- val <- stat <- 'fake_global_for_CRAN'

  to_plot <- dplyr::data_frame(phenotype_name = qtl::pull.pheno(cross = cross)[[phenotype_name]],
                               marker_name = factor(x = qtl::pull.geno(cross = cross)[,marker_name]))

  if (!is.null(genotype_labels)) {
    to_plot[['marker_name']] <- factor(x = to_plot[['marker_name']],
                                       labels = genotype_labels)
  }

  the_plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = marker_name, y = phenotype_name))


  if (!is.null(color_by)) {

    to_add <- NULL
    if (color_by %in% names(qtl::pull.pheno(cross = cross))) {
      to_add <- qtl::pull.pheno(cross = cross)[[color_by]]
      if (length(unique(to_add)) < 10) { to_add <- factor(to_add) }
    } else if (color_by %in% colnames(qtl::pull.geno(cross = cross))) {
      to_add <- factor(qtl::pull.geno(cross = cross)[[color_by]])
    }

    if (is.null(to_add)) {
      stop('color_by argument, ', color_by, ', not found in phenotypes or genotypes')
    }

    to_add <- stats::setNames(object = dplyr::data_frame(placeholder_name = to_add), nm = color_by)
    to_plot <- dplyr::bind_cols(to_plot, to_add)
  }

  if (!is.null(shape_by)) {
    to_add <- NULL
    if (shape_by %in% names(qtl::pull.pheno(cross = cross))) {
      to_add <- qtl::pull.pheno(cross = cross)[[shape_by]]
      if (length(unique(to_add)) < 10) { to_add <- factor(to_add) }
    }

    if (is.null(to_add)) {
      stop('shape_by argument, ', shape_by, ', not found in phenotypes')
    }

    to_add <- stats::setNames(object = dplyr::data_frame(placeholder_name = to_add), nm = shape_by)
    to_plot <- dplyr::bind_cols(to_plot, to_add)
  }

  if (is.character(point_size)) {

    to_add <- NULL
    if (point_size %in% names(qtl::pull.pheno(cross = cross))) {
      to_add <- qtl::pull.pheno(cross = cross)[[point_size]]
      if (length(unique(to_add)) < 10) { to_add <- factor(to_add) }
    }

    if (is.null(to_add)) {
      stop('point_size argument, ', point_size, ', not found in phenotypes')
    }

    to_add <- stats::setNames(object = dplyr::data_frame(placeholder_name = to_add), nm = point_size)
    to_plot <- dplyr::bind_cols(to_plot, to_add)

    the_plot <- the_plot +
      ggplot2::geom_jitter(data = stats::na.omit(to_plot),
                           mapping = ggplot2::aes_string(color = color_by, shape = shape_by, size = point_size),
                           width = 0.3,
                           alpha = point_alpha)
  } else {

    the_plot <- the_plot +
      ggplot2::geom_jitter(data = stats::na.omit(to_plot),
                           mapping = ggplot2::aes_string(color = color_by, shape = shape_by),
                           width = 0.3,
                           size = point_size,
                           alpha = point_alpha)
  }





  if (Ibars) {

      genotype_summaries <- to_plot %>%
          dplyr::group_by(marker_name) %>%
          dplyr::summarise(mean = mean(phenotype_name, na.rm = TRUE),
                           sd = stats::sd(phenotype_name, na.rm = TRUE)) %>%
          stats::na.omit()

      the_plot <- the_plot +
          ggplot2::geom_point(data = genotype_summaries,
                              mapping = ggplot2::aes(x = marker_name, y = mean),
                              size = 3) +
          ggplot2::geom_segment(data = genotype_summaries,
                                mapping = ggplot2::aes(x = marker_name, xend = marker_name,
                                                       y = mean - sd, yend = mean + sd)) +
          ggplot2::geom_segment(data = genotype_summaries,
                                mapping = ggplot2::aes(x = as.numeric(marker_name) - 0.15, xend = as.numeric(marker_name) + 0.15,
                                                       y = mean - sd, yend = mean - sd)) +
          ggplot2::geom_segment(data = genotype_summaries,
                                mapping = ggplot2::aes(x = as.numeric(marker_name) - 0.15, xend = as.numeric(marker_name) + 0.15,
                                                       y = mean + sd, yend = mean + sd))


          # ggplot2::geom_errorbar(data = genotype_summaries,
          #                        mapping = ggplot2::aes(x = marker_name, ymin = mean - sd, ymax = mean + sd),
          #                        width = 0.3)
  }

  if (all(Ibars, connectIbars)) {

      genotype_summary_lines <- genotype_summaries %>%
          dplyr::mutate(lb = mean - sd, ub = mean + sd) %>%
          dplyr::select(-sd) %>%
          tidyr::gather(key = 'stat', value = 'val', mean:ub)

      # todo: figure out the right geom here
      the_plot <- the_plot +
          ggplot2::geom_line(data = stats::na.omit(genotype_summary_lines),
                             mapping = ggplot2::aes(x = marker_name, y = val, group = stat),
                             linetype = 2)
  }

  the_plot <- the_plot +
    ggplot2::xlab(label = marker_name) +
    ggplot2::ylab(label = phenotype_name) +
    ggplot2::theme_minimal()

  return(the_plot)
}
