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
#' @param plotting.units One of 'LOD', 'asymp.p', or 'empir.p', implying whether
#' LOD scores, asymptotic p-values, or empirical p-values should be plotted.
#' Defaults to 'LOD'
#' @param incl.markers Optionally, whether to draw a rug plot along the bottom indicating where the markers are.  Defaults to TRUE.
#' @param show.equations Optionally, whether to write the modeling equations under the title.  Defaults to TRUE.
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
#' @seealso  \code{\link{scanonevar}}, \code{\link{scanonevar.to.p.values}}
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
#' plot(test.sov)
#'
plot.scanonevar <- function(x,
                            y = NULL,
                            chrs = unique(x[['chr']]),
                            plotting.units = c('LOD', 'asymp.p', 'empir.p'),
                            plot.title = attr(x = x, which = 'mean.formula')[[2]],
                            marker.rug = TRUE)
{

  plotting.units <- match.arg(plotting.units)

  # can only plot if x is a scanonevar and y, if present, is a scanone
  stopifnot(is.scanonevar(x))
  if (!is.null(y)) {
    stopifnot('scanone' %in% class(y))
  }

  # filter down to requested chromosomes, and make y a tbl_df
  x <- dplyr::filter(.data = x, chr %in% chrs)
  if (!is.null(y)) {
    y <- dplyr::filter(.data = dplyr::tbl_df(y), chr %in% chrs)
  }

  # pull necessary columns into to.plot
  to.plot <- pull.plotting.columns_(sov = x, so = y, plotting.units = plotting.units)

  p <- ggplot2::ggplot(data = to.plot) +
    ggplot2::geom_line(ggplot2::aes(x = pos, y = val, col = test)) +
    ggplot2::facet_grid(facets = ~chr,
                        scales = 'free_x', space = 'free_x', switch = 'x',
                        labeller = function(labels) ggplot2::label_both(labels = labels,
                                                                        sep = '')) +
    ggplot2::scale_color_manual(name = c('mQTL', 'vQTL', 'mvQTL', 'traditional'),
                                values = c('blue', 'red', 'black', 'darkgreen')) +
    ggplot2::ggtitle(label = plot.title) +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = 'white'),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = 'lightgray'),
                   strip.background = ggplot2::element_rect(fill = 'lightgray'))

  if (marker.rug) {
    true.markers <- x %>% filter(pos != round(pos))
    p <- p + ggplot2::geom_rug(mapping = ggplot2::aes(x = pos),
                               data = true.markers)
  }

  return(p)
}



pull.plotting.columns_ <- function(sov, so, plotting.units) {

  base.to.plot <- dplyr::select(.data = sov,
                           loc.name, chr, pos)

  # pull appropriate data into plotting columns
  if (plotting.units == 'LOD') {
    to.plot <- dplyr::bind_rows(dplyr::mutate(.data = base.to.plot,
                                              test = 'mQTL',
                                              val = sov[['mean.lod']]),
                                dplyr::mutate(.data = base.to.plot,
                                              test = 'vQTL',
                                              val = sov[['var.lod']]),
                                dplyr::mutate(.data = base.to.plot,
                                              test = 'mvQTL',
                                              val = sov[['joint.lod']]))
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
    to.plot[['mean.to.plot']] <- -log10(sov[['mean.empir.p']])
    to.plot[['var.to.plot']] <- -log10(sov[['var.empir.p']])
    to.plot[['joint.to.plot']] <- -log10(sov[['joint.empir.p']])
    if (is.null(so)) {
      to.plot[['so.to.plot']] <- NA
    } else {
      to.plot[['so.to.plot']] <- -log10(so[['empir.p']])
    }
  }

  # straighten up factors
  to.plot$chr <- factor(to.plot$chr, levels = gtools::mixedsort(unique(to.plot$chr)))
  vqtl.test.names <- c('mQTL', 'vQTL', 'mvQTL')
  if (all(unique(to.plot$test) %in% vqtl.test.names)) {
    to.plot$test <- factor(to.plot$test, levels = vqtl.test.names)
  } else if (all(unique(to.plot$test) %in% c(vqtl.test.names, 'traditional'))) {
    to.plot$test <- factor(to.plot$test, levels = c(vqtl.test.names, 'traditional'))
  } else {
    stop('unrecognized test type')
  }


  return(to.plot)
}



lod2pval <- function(lod, df) {
  return(pchisq(q = lod, df = df, lower.tail = FALSE))
}

