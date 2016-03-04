#'  @title Summary of Peaks in Scanonevar
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{summary.scanonevar} prints out the loci in a scanonevar object
#'    that exceed \code{thresh}.  It is an S3 generic for summary().  It handles scanonevar
#'    objects in both LOD units and empirical p value units.
#'
#'  @param object the scanonevar object to be summarized
#'  @param thresh the threshold over which (for LODs) or under which (for emprirical p values)
#'    a locus will be printed.
#'  @param ... additional arguments controlling the summary
#'
#'  @return None.  Only prints results to screen.
#'
#'  @details none
summary.scanonevar <- function(object, thresh, ...) {

  # hack to get R CMD CHECK to run without NOTEs that these globals are undefined
  full.peak <- full.lod <- matches <- mean.peak <- mean.lod <- var.peak <- var.lod <- 'fake.global'
  emp.p.full.lod <- emp.p.mean.lod <- emp.p.var.lod <- 'fake.global'
  chr <- pos <- marker.name <- 'fake.global'

  if (!any(class(object) == "scanonevar")) {
    stop("Input should have class \"scanonevar\".")
  }

  if (missing(thresh)) {
    if (units(object) == 'lods') { thresh <- 3 }
    if (units(object) == 'emp.ps') { thresh <- 0.05 }
  }


  peaks <- get.peaks.from.scanonevar(object, thresh)

  if (units(object) == 'lods') {
    message('Full Model Peaks:')
    print(peaks %>%
            dplyr::filter(full.peak == TRUE, full.lod > thresh) %>%
            select(chr, pos, marker.name, full.lod))
    message('Mean Model Peaks:')
    print(peaks %>%
            dplyr::filter(mean.peak == TRUE, mean.lod > thresh) %>%
            select(chr, pos, marker.name, mean.lod))
    message('Var Model Peaks:')
    print(peaks %>%
            dplyr::filter(var.peak == TRUE, var.lod > thresh) %>%
            select(chr, pos, marker.name, var.lod))
  }

  if (units(object) == 'emp.ps') {
    message('Full Model Peaks:')
    print(peaks %>%
            dplyr::filter(full.peak == TRUE, emp.p.full.lod < thresh) %>%
            select(chr, pos, marker.name, emp.p.full.lod))
    message('Mean Model Peaks:')
    print(peaks %>%
            dplyr::filter(mean.peak == TRUE, emp.p.mean.lod < thresh) %>%
            select(chr, pos, marker.name, emp.p.mean.lod))
    message('Var Model Peaks:')
    print(peaks %>%
            dplyr::filter(var.peak == TRUE, emp.p.var.lod < thresh) %>%
            select(chr, pos, marker.name, emp.p.var.lod))
  }
}
