#'  @title Plot Phenotype of interest Averaged (Marginalized) Across Specified Markers and Phenotypes
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{margin.plot} should be used to visually investigate the relationship
#'    between the phenotype of interest and other phenotypes.  \code{margin.plot} can also
#'    be used to visualize the relationship between the phenotype of interest and genetic
#'    loci of interest, but \code{predictive.plot} is usually preferrable.
#'
#'  @param cross The cross object to be plotted
#'  @param focal.phenotype.name the phenotype to put on the y-axis
#'  @param marginal.phen.names a list of phenotypes to average over (put on the x-axis).
#'  @param marginal.marker.names a list of marker names, whose values will be averaged over (put on the x-axis).
#'  @param genotype.plotting.names Labels for the genotype groups.  Defaults to \code{c('AA', 'AB', 'BB')}.
#'  @param subset the subset of individuals to use
#'  @param col optionally, color of dots, as in base R graphics.  Defaults to gray.
#'  @param pch optionally, plotting character, as in base R graphics.  Defaults to 19 (disc).
#'  @param my.xlab optionally, x axis label, as in base R graphics.  Defaults to the name of the marginal marker.
#'  @param my.ylab optionally, y axis label, as in base R graphics.  Defaults to focal phenotype name.
#'  @param title optionally, plot title, as in base R graphics.  Defaults to 'focal phenotype name by marginal phenotype name'.
#'  @param title.cex optionally, character expansion for title, as in base R graphics.  Defaults to 1.5.
#'  @param circle.alpha optionally, alpha (transparency) of discs.  Defaults to 0.2.
#'
#'  @return None.  Only makes plot.
#'
#'  @details none
#'
#'  @examples
#'  \dontrun{
#'    margin.plot(cross = my.cross,
#'                focal.phenotype.name = my.phenotype,
#'                marginal.phen.name = list('sex', 'age'),
#'                marginal.marker.name = list('chrA_markerB', 'chrC_markerD'))
#'
#'  }
#'
#'
margin.plot <- function(cross,
                        focal.phenotype.name,
                        marginal.phen.names = NULL,
                        marginal.marker.names = NULL,
                        genotype.plotting.names = c('A', 'H', 'B'),
                        subset = 1:nind(cross),
                        col = rep(rgb(0.5, 0.5, 0.5, 0.5), nind(cross)),
                        pch = 19,
                        my.xlab = marginal.marker.name,
                        my.ylab = focal.phenotype.name,
                        title = paste(focal.phenotype.name, 'by', marginal.phen.name),
                        title.cex = 1.5,
                        circle.alpha = 0.2) {

  if (any(missing(cross), missing(focal.phenotype.name), !(focal.phenotype.name %in% names(cross$pheno)))) {
    stop('Must provide a cross and a focal phenotype in that cross.')
  }

  num.plots <- sum(length(marginal.phen.names), length(marginal.marker.names))
  if (num.plots == 0) { stop('Must provide a marginal phenotype or marker.')}

  focal.phen <- cross$pheno[[focal.phenotype.name]][subset]
  if (!missing(col) & length(col) == nind(cross)) {
    col <- col[subset]
  }

  for (marginal.phen.name in marginal.phen.names) {

    marginal.phen <- cross$pheno[[marginal.phen.name]][subset]
    if (is.factor(marginal.phen)) {
      lev.names <- levels(marginal.phen)
      plotting.phen <- as.numeric(marginal.phen)
    } else {
      plotting.phen <- marginal.phen
    }

    plot(x = jitter(plotting.phen),
         y = focal.phen,
         xlab = marginal.phen.name,
         ylab = NA,
         xaxt = 'n',
         main = NA,
         col = alpha(col, 0.5),
         pch = pch,
         axes = FALSE)
    mtext(text = title, side = 3, line = 1, cex = title.cex)

    if (is.factor(marginal.phen)) {
      axis(side = 1, at = unique(marginal.phen), labels = lev.names, tick = TRUE)
    }

  }

  for (marginal.marker.name in marginal.marker.names) {

    chr.of.interest <- which(sapply(X = cross$geno, FUN = function(chr) { marginal.marker.name %in% colnames(chr$data)}))

    genotypes <- cross$geno[[chr.of.interest]]$data[subset,marginal.marker.name]
    plot(x = jitter(genotypes),
         y = focal.phen,
         xaxt = 'n',
         xlab = NA,
         ylab = NA,
         axes = FALSE,
         col = alpha(col, circle.alpha),
         pch = pch,
         cex.main = title.cex)
    axis(side = 2)
    mtext(text = genotype.plotting.names, side = 1, line = 0, at = 1:3)
    mtext(side = 1, text = my.xlab, line = 2)
    mtext(side = 2, text = my.ylab, line = 2)
    mtext(side = 3, text = title, line = 1, cex = title.cex)

    means <- aggregate(x = focal.phen, by = list(genotypes), FUN = mean)[,2]
    sds <- aggregate(x = focal.phen, by = list(genotypes), FUN = sd)[,2]

    x.start.wide <- 1:3 - 0.2
    x.end.wide <- 1:3 + 0.2
    x.start.narrow <- 1:3 - 0.1
    x.end.narrow <- 1:3 + 0.1

    # horizontal lines at means and means +/- 1SD
    segments(x0 = c(x.start.narrow, x.start.wide, x.start.wide),
             y0 = c(means, means - sds, means + sds),
             x1 = c(x.end.narrow, x.end.wide, x.end.wide),
             y1 = c(means, means - sds, means + sds),
             lwd = rep(c(4, 2, 2), each = 3))

    # vertical line down the middle of each genotype group
    segments(x0 = 1:3,
             y0 = means - sds,
             x1 = 1:3,
             y1 = means + sds,
             lwd = 2)

    # dotted horizontal lines connecting means and SD's
    segments(x0 = rep(1:2, 3),
             y0 = c(means[1:2], means[1:2] + sds[1:2], means[1:2] - sds[1:2]),
             x1 = rep(2:3, 3),
             y1 = c(means[2:3], means[2:3] + sds[2:3], means[2:3] - sds[2:3]),
             lty = 2)
  }

  # reset graphical parameteers to how they were on start
  # par(start.pars)

  # return nothing
  invisible()
}