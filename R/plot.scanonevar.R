#' @title plot.scanonevar
#'
#' @author Robert Corty \email{rcorty@@gmail.com}
#'
#' @description \code{plot.scanonevar} implements the plot generic for objects of class 'scanonevar'.
#' Because scanonevar objects can be viewed in terms of LODs or empirical p-values,
#' this plotting function checks the 'units' attribute to determine which to plot.
#'
#' @param x the \code{scanonevar} object to be plotted
#' @param y Optionally, a \code{scanone} object to be plotting for comparison to the \code{scanonevar} object.
#' @param chrs Optionally, the subset of the chromosomes to plot
#' @param tests_to_plot which one or ones of the three possible tests to plot ('mQTL', 'vQTL', and 'mvQTL')
#' @param plotting.units One of 'LOD', 'asymp.p', or 'empir.p', implying whether
#' LOD scores, asymptotic p-values, or empirical p-values should be plotted.
#' Defaults to 'LOD'
#' @param plot.title the title of the plot
#' @param marker.rug Should a marker rug be plotted? Defaults to TRUE.
#' @param ymax the top of the y axis
#' @param legend_pos the position of the legend
#' @param alpha_pos the position of the alpha values (false positive rate)
#' #@param show.equations Optionally, whether to write the modeling equations under the title.  Defaults to TRUE
#' @param ... additional plotting arguments
#'
#' @return Returns the plot.
#'
#' @details If such a strong signal was observed that the empirical p-value underflows R's
#'    float type, this function produces an error.  The author is open to suggestions on how
#'    to deal with this situation better.
#'
#'    These plots look better when both x (the scanonevar object) and y (optional scanone
#'    for comparison) are in units p values than when they are in LOD units.
#'
#' @importFrom dplyr %>%
#' @export
#'
#' @details none
#'
#' @examples
#' set.seed(27599)
#' test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), n.mar = 5), n.ind = 50)
#' test.sov <- scanonevar(cross = test.cross)
#' plot(x = test.sov)
#'
plot.scanonevar <- function(x,
                            y = NULL,
                            chrs = unique(x[['result']][['chr']]),
                            tests_to_plot = c('mQTL', 'vQTL', 'mvQTL'),
                            plotting.units = if(any(grepl(pattern = 'empir.p', x = names(x[['result']])))) {
                              'empir.p'
                            } else {
                              'LOD'
                            },
                            plot.title = x[['meta']][['formulae']][['mean.alt.formula']][[2]],
                            marker.rug = TRUE,
                            ymax = NULL,
                            legend_pos = NULL,
                            alpha_pos = c('left', 'right', 'none'),
                            ...) {

  chr <- lab <- pos <- val <- test <- loc.name <- 'fake_global_for_CRAN'

  alpha_pos <- match.arg(arg = alpha_pos)

  # can only plot if x is a scanonevar and y, if present, is a scanone
  stopifnot(is.scanonevar(x))
  if (!is.null(y)) {
    stopifnot('scanone' %in% class(y))
  }

  # filter down to requested chromosomes, and make y a tbl_df
  result <- dplyr::filter(.data = x[['result']], chr %in% chrs)
  if (!is.null(y)) {
    y <- dplyr::filter(.data = dplyr::tbl_df(y), chr %in% chrs)
  }

  # pull necessary columns into to.plot
  to.plot <- pull.plotting.columns_(sov = result,
                                    so = y,
                                    tests_to_plot = tests_to_plot,
                                    plotting.units = plotting.units)
  to.plot <- dplyr::mutate(.data = to.plot,
                           chr = factor(x = chr, levels = gtools::mixedsort(unique(chr)))) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(pos = pos - min(pos))

  p <- ggplot2::ggplot(data = to.plot)

  if (plotting.units == 'LOD') {
    p <- p +
      ggplot2::ylab(label = 'LOD') +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = 'lightgray'))
  } else {
    p <- p +
      ggplot2::ylab(label = '-log10(p)') +
      ggplot2::geom_hline(yintercept = -log10(c(0.05, 0.01)), color = 'gray')

    if (alpha_pos != 'none') {
      p <- p +
        ggplot2::geom_text(mapping = ggplot2::aes(x = x,
                                                  y = y,
                                                  label = lab),
                           data = data.frame(x = switch(EXPR = alpha_pos,
                                                        left = 0,
                                                        right = max(to.plot$pos[to.plot$chr == chrs[1]])),
                                             y = -log10(c(0.05, 0.01)),
                                             lab = c("alpha == 0.05",
                                                     "alpha == 0.01"),
                                             chr = chrs[1]),
                           vjust = 0,
                           hjust = switch(EXPR = alpha_pos,
                                          left = 0,
                                          right = 1),
                           parse = TRUE)
    }

    p <- p +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }

  p <- p +
    ggplot2::geom_line(ggplot2::aes(x = pos, y = val, col = test)) +
    ggplot2::scale_color_manual(name = 'test',
                                breaks = c('mQTL', 'vQTL', 'mvQTL', 'traditional'),
                                values = c('blue', 'red', 'black', 'darkgreen'),
                                drop = FALSE) +
    ggplot2::ggtitle(label = plot.title) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = 'white'),
                   panel.grid.major.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = 'lightgray'),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.title = ggplot2::element_blank())

  if (length(chrs) == 1) {
    p <- p +
      ggplot2::xlab(label = paste('Chromosome', chrs))
      # ggplot2::facet_grid(facets = ~chr,
      #                     scales = 'free_x', space = 'free_x', switch = 'x',
      #                     labeller = function(labels) ggplot2::label_both(labels = labels,
      #                                                                     sep = ''))
  } else {
    p <- p +
      ggplot2::facet_grid(facets = ~chr,
                          scales = 'free_x', space = 'free_x', switch = 'x') +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())

  }



  if (marker.rug) {
    # true.markers <- result %>% dplyr::filter(pos != round(pos))
    true.markers <- result %>%
      dplyr::filter(!grepl(pattern = '^chr[0-9]*_loc[0-9]*$', x = loc.name)) %>%
      dplyr::mutate(chr = factor(x = chr, levels = gtools::mixedsort(unique(chr)))) %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(pos = pos - min(pos))
      #dplyr::mutate(chr = stringr::str_pad(string = chr, width = 2, pad = '0'))
    p <- p + ggplot2::geom_rug(mapping = ggplot2::aes(x = pos),
                               data = true.markers)
  }

  if (!is.null(ymax)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, ymax))
  }

  if (!is.null(legend_pos)) {
    p <- p + ggplot2::theme(legend.position = legend_pos)
  }

  return(p)
}



pull.plotting.columns_ <- function(sov, so, tests_to_plot, plotting.units) {

  loc.name <- chr <- pos <- 'fake_global_for_CRAN'

  base.to.plot <- dplyr::select(.data = sov,
                                loc.name, chr, pos)

  to.plot <- dplyr::data_frame(test = NA, val = NA)

  # pull appropriate data into plotting columns
  if (plotting.units == 'LOD') {

    if (!is.null(sov[['mean.lod']]) & 'mQTL' %in% tests_to_plot) {
      to.plot <- dplyr::bind_rows(to.plot,
                                  dplyr::mutate(.data = base.to.plot,
                                                test = 'mQTL',
                                                val = sov[['mean.lod']]))
    }

    if (!is.null(sov[['var.lod']]) & 'vQTL' %in% tests_to_plot) {
      to.plot <- dplyr::bind_rows(to.plot,
                                  dplyr::mutate(.data = base.to.plot,
                                                test = 'vQTL',
                                                val = sov[['var.lod']]))
    }

    if (!is.null(sov[['joint.lod']]) & 'mvQTL' %in% tests_to_plot) {
      to.plot <- dplyr::bind_rows(to.plot,
                                  dplyr::mutate(.data = base.to.plot,
                                                test = 'mvQTL',
                                                val = sov[['joint.lod']]))
    }

    if (!is.null(so)) {
      to.plot <- dplyr::bind_rows(to.plot,
                                  dplyr::mutate(.data = base.to.plot,
                                                test = 'traditional',
                                                val = so[['lod']]))
    }
  }

  if (plotting.units == 'asymp.p') {

    to.plot <- dplyr::bind_rows(dplyr::mutate(.data = base.to.plot,
                                              test = 'mQTL',
                                              val = -log10(sov[['mean.asymp.p']])),
                                dplyr::mutate(.data = base.to.plot,
                                              test = 'vQTL',
                                              val = -log10(sov[['var.asymp.p']])),
                                dplyr::mutate(.data = base.to.plot,
                                              test = 'mvQTL',
                                              val = -log10(sov[['joint.asymp.p']])))
    if (!is.null(so)) {
      to.plot <- dplyr::bind_rows(to.plot,
                                  dplyr::mutate(.data = base.to.plot,
                                                test = 'traditional',
                                                val = -log10(lod2pval(so[['lod']], df = 2))))
    }
  }

  if (plotting.units == 'empir.p') {

    to.plot <- dplyr::bind_rows(dplyr::mutate(.data = base.to.plot,
                                              test = 'mQTL',
                                              val = -log10(sov[['mean.empir.p']])),
                                dplyr::mutate(.data = base.to.plot,
                                              test = 'vQTL',
                                              val = -log10(sov[['var.empir.p']])),
                                dplyr::mutate(.data = base.to.plot,
                                              test = 'mvQTL',
                                              val = -log10(sov[['joint.empir.p']])))

    if (!is.null(so)) {
      if ('empir.p' %in% names(so)) {
        to.plot <- dplyr::bind_rows(to.plot,
                                    dplyr::mutate(.data = base.to.plot,
                                                  test = 'traditional',
                                                  val = -log10(so[['empir.p']])))
      } else {
        stop('Plotting unit is empir.p, but provided scanone doesnt have empir.p')
      }
    }
  }

  # remove the top NA (placeholder) row
  to.plot <- to.plot[-1,]

  # straighten up factors
  to.plot$chr <- factor(to.plot$chr, levels = gtools::mixedsort(unique(to.plot$chr)))
  vqtl.test.names <- c('mQTL', 'vQTL', 'mvQTL')
  if (all(unique(to.plot$test) %in% vqtl.test.names)) {
    to.plot$test <- factor(to.plot$test, levels = vqtl.test.names)
  } else if (all(unique(to.plot$test) %in% c(vqtl.test.names, 'traditional'))) {
    to.plot$test <- factor(x = to.plot$test,
                           levels = c(vqtl.test.names, 'traditional'),
                           labels = c(vqtl.test.names, 'traditional'))
  } else {
    stop('unrecognized test type')
  }


  return(to.plot)
}



lod2pval <- function(lod, df) {
  return(stats::pchisq(q = lod, df = df, lower.tail = FALSE))
}

