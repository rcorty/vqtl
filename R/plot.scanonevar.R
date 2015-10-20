#'  @title plot.scanonevar
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{plot.scanonevar} implements the plot generic for objects of class 'scanonevar'.
#'    Because scanonevar objects can be viewed in terms of LODs or empirical p-values,
#'    this plotting function checks the 'units' attribute to determine which to plot.
#'
#'  @param x the \code{scanonevar} object to be plotted
#'  @param y Optionally, a \code{scanone} object to be plotting for comparison to the \code{scanonevar} object.
#'  @param chrs Optionally, the subset of the chromosomes to plot
#'  @param col Optionally, a vector specifying the colors of the scan lines.  Defaults to \code{c("black", "blue", "red", "darkgreen")}.
#'  @param bandcol Optionally, a background color for the even-index chromosomes in this scan.
#'  @param legends Optionally, the name to put for each scan in the legend.  Defaults to \code{c('mean or var', 'mean', 'var', 'scanone for comparison')}.
#'  @param legend.pos Optionally, the corner/edge where the legend should be drawn.  Defaults to \code{'top'}.
#'  @param gap Optionally, The space between chromosomes in Mb.  Defaults to 25.
#'  @param incl.markers Optionally, whether to draw a rug plot along the bottom indicating where the markers are.  Defaults to TRUE.
#'  @param main Optionally, a title for the plot.  Defaults to the name of the phenotype from the scanonevar.
#'  @param ylim Optionally, the y limits for the plot.  Defaults to \code{c(0, 1.05 * highest.point)}.
#'  @param show.equations Optionally, whether to write the modeling equations under the title.  Defaults to TRUE.
#'  @param ... Additionaly graphical parameters.
#'
#'  @return Returns nothing.  Only makes the plot.
#'
#'  @details If such a strong signal was observed that the empirical p-value underflows R's
#'    float type, this function produces an error.  The author is open to suggestions on how
#'    to deal with this situation better.
#'
#'  @seealso  \code{\link{scanonevar}}, \code{\link{convert.scanonevar.to.p.values}}
#'
#'  @details none
#'

plot.scanonevar <- function(x,
                            y = NULL,
                            chrs = unique(x$chr),
                            col = c("black", "blue", "red", "darkgreen"),
                            bandcol = 'lightgray',
                            legends = c('mean or var', 'mean', 'var', 'scanone for comparison'),
                            legend.pos = 'top',
                            gap = 25,
                            incl.markers = TRUE,
                            main = attr(x, 'pheno'),
                            ylim = c(0, 1.05*max(coords.y.locus, na.rm = TRUE)),
                            show.equations = (length(chrs) != 1),
                            ...)
{

  # hack to get R CMD CHECK to run without NOTEs that these globals are undefined
  chr <- pos <- len <- matches <- 'fake.global'

  # validate scanonevar object
  if (!is.scanonevar(x)) {
    stop(paste('scanonevar argument is not a valid scanonevar object:', attr(is.scanonevar(x), 'why.not')))
  }

  # store current graphical parameters and customize them for this plot
  start.pars <- par(no.readonly = TRUE)
  if (show.equations) {
    par(mar = c(3,4,5,1))
  } else {
    par(mar = c(3,4,2,1))
  }

  # convert scanone.for.comparison to tbl_df if needed
  if (!any(is.null(y), is.tbl(y))) {
    y <- tbl_df(y)
  }

  # subset scanonevar to necessary chrs only
  if (!identical(chrs, unique(x$chr))) {

    temp <- dplyr::filter(x, chr %in% chrs)
    y <- dplyr::filter(y, chr %in% chrs)

    class(temp) <- class(x)
    attr(temp, 'method') <- attr(x, 'method')
    attr(temp, 'type') <- attr(x, 'type')
    attr(temp, 'model') <- attr(x, 'model')
    attr(temp, 'mean.null.formula') <- attr(x, 'mean.null.formula')
    attr(temp, 'mean.alt.formula') <- attr(x, 'mean.alt.formula')
    attr(temp, 'var.null.formula') <- attr(x, 'var.null.formula')
    attr(temp, 'var.alt.formula') <- attr(x, 'var.alt.formula')
    attr(temp, 'pheno') <- attr(x, 'pheno')
    attr(temp, 'null.effects') <- attr(x, 'null.effects')
    attr(temp, 'units') <- attr(x, 'units')
    attr(temp, 'null.fit') <- attr(x, 'null.fit')
    temp$chr <- factor(temp$chr)

    x <- temp
  }

  # x coordinates for plotting
  x$chr <- factor(x$chr)
  levels(x$chr) <- mixedsort(levels(x$chr))
  coords.x.chr <- x %>%
    group_by(chr) %>%
    summarise(len = max(pos) - min(pos)) %>%
    mutate(start = cumsum(dplyr::lag(len + gap, default = 0))) %>%
    mutate(end = start + len) %>%
    mutate(middle = (start + end)/2)

  coords.x.locus <- left_join(coords.x.chr, x, by = 'chr') %>%
    mutate(coord.x = start + pos)

  # y coordinates for plotting
  if (attr(x, 'units') == 'lods') {
    coords.y.locus <- x %>% select(matches('lod'))
    ylab <- 'LOD'
  }
  if (attr(x, 'units') == 'emp.ps') {
    coords.y.locus <- -log10(select(x, matches('emp.p')))
    ylab = '-log10(empirical p)'
  }

  # make plotting area
  xlim <- c(-gap/2, max(coords.x.chr$end) + gap/2)
  plot(-42, -42, xlim = xlim, ylim = ylim,
       type = 'n', xaxs = 'i',
       xlab = NA, ylab = NA, axes = FALSE)

  # shade in bg for every other chr
  if (!is.null(bandcol) & length(chrs) > 1) {
    names.chrs.even <- chrs[seq(from = 2, to = length(chrs), by = 2)]
    xs.chrs.even <- dplyr::filter(coords.x.chr, chr %in% names.chrs.even)

    rect(xleft = xs.chrs.even$start - gap/2,
         xright = xs.chrs.even$end + gap/2,
         ybottom = ylim[1],
         ytop = ylim[2],
         border = bandcol, col = bandcol)
  }

  # draw x axis and label chrs
  segments(x0 = xlim[1], x1 = xlim[2],
           y0 = 0, y1 = 0)
  mtext(side = 1, text = chrs, at = coords.x.chr$middle, line = 0)
  mtext(side = 1, text = 'Chromosome', at = mean(xlim), line = 1)

  # draw y axis and label chrs
  axis(side = 2)
  mtext(side = 2, text = ylab, line = 2, at = mean(ylim))

  # plot lines
  for (test.idx in 1:ncol(coords.y.locus)) {
    test.name <- names(coords.y.locus)[test.idx]
    segments(x0 = coords.x.locus$coord.x,
             x1 = lead(coords.x.locus$coord.x),
             y0 = coords.y.locus[[test.name]],
             y1 = lead(coords.y.locus[[test.name]]),
             lwd = 2*with(coords.x.locus, pos != len),
             col = col[test.idx])
  }

  # plot scanone for comparison
  if (attr(x, 'units') == 'lods') {
    if (!is.null(y)) {
      segments(x0 = coords.x.locus$coord.x,
               x1 = lead(coords.x.locus$coord.x),
               y0 = y$lod,
               y1 = lead(y$lod),
               lwd = 2*with(coords.x.locus, pos != len),
               col = col[length(col)])
    }
  }
  if (attr(x, 'units') == 'emp.ps') {
    if (!is.null(y)) {
      segments(x0 = coords.x.locus$coord.x,
               x1 = lead(coords.x.locus$coord.x),
               y0 = -log10(y$emp.p.scanone),
               y1 = -log10(lead(y$emp.p.scanone)),
               lwd = 2*with(coords.x.locus, pos != len),
               col = col[length(col)])
    }
  }

  # note: I got rid of thresholds on the LOD scale, because it is visually confusing to
  # look at 3 or 4 different alpha05 thresholds and/or 3 or 4 different alpha 01's  -RWC

  # add thresholds on p-value scale
  if (attr(x, 'units') == 'emp.ps') {
    abline(h = -log10(c(0.05, 0.01)), lty = c(1, 2))
    text(x = coords.x.locus$coord.x[1], y = -log10(0.05), labels = 'alpha=0.05', adj = c(0.5, -0.2))
  }

  # plot lines at marker positions (rug plot)
  if (incl.markers) {
    marker.idxs <- !grepl(pattern = 'loc', x = x$marker.name)
    segments(x0 = coords.x.locus$coord.x[marker.idxs],
             x1 = coords.x.locus$coord.x[marker.idxs],
             y0 = 0,
             y1 = ylim[2]*0.03,
             col = 'gray50')
  }

  # draw the legend
  legend(x = legend.pos,
         legend = legends,
         fill = col, cex = 0.8, bty = 'n',
         x.intersp = 0.3, y.intersp = 1, xjust = 0.5, yjust = 0)

  # add the title
  if (show.equations) {
    title <- paste(attr(x = x, 'pheno'),
                   '\n', 'mean null:',
                   paste(as.character(attr(x, 'mean.null.formula'))[c(2,1,3)], collapse = ' '),
                   '\n', 'mean alt:',
                   paste(as.character(attr(x, 'mean.alt.formula'))[c(2,1,3)], collapse = ' '),
                   '\n', 'var null:',
                   paste(as.character(attr(x, 'var.null.formula')), collapse = ' '),
                   '\n', 'var alt:',
                   paste(as.character(attr(x, 'var.alt.formula')), collapse = ' '))
    mtext(text = title, side = 3, line = 0)
  } else {
    mtext(text = attr(x, 'pheno'), side = 3, line = 0)
  }

  # reset graphical parameteers to how they were on start
  par(start.pars)

  # return nothing
  invisible()
}
