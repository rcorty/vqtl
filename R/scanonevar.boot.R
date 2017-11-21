#' @title scanonevar.boot
#' @name scanonevar.boot
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @description \code{scanonevar.boot} conducts a nonparametric bootstrap of one chromosome
#' to establish a confidence interval on any peaks
#'
#' @param sov the scanonevar whose significance should be assessed empirically in an FWER-controlling method
#' @param n.resamples the number of resamples
#' @param chr which chromosome to focus on
#' @param random.seed value to start the random number generator at, for reproducibility
#' @param n.cores number of cores to use for the permutations
#' @param silent Should all messaging be suppressed?
#' @param qtl_type which type of QTL did you detect and want a CI for?  mQTL, vQTL, or mvQTL.
#'
#' @return 27599
#' @export
#'
#' @importFrom foreach %dopar%
#'
#' @examples
#' set.seed(27599)
#' test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5), n.mar = 5), n.ind = 50)
#' sov <- scanonevar(cross = test.cross)
#'
scanonevar.boot <- function(sov,
                            n.resamples,
                            chr,
                            qtl_type = c('mQTL', 'vQTL', 'mvQTL'),
                            # bootstrap_type = c('simple', 'bayesian'),
                            random.seed = 27599,
                            n.cores = 1,
                            silent = FALSE) {

  stopifnot(is.scanonevar(sov))
  stopifnot(n.resamples > 0)
  stopifnot(chr %in% names(sov$meta$cross$geno))
  qtl_type <- match.arg(arg = qtl_type)
  # bootstrap_type <- match.arg(arg = bootstrap_type)

  # no way to split perms over cores...so never use more cores than there are perms
  n.cores <- min(n.cores, n.resamples)

  # execute the scan
  result <- scanonevar.boot_(sov = sov[['result']],
                             meta = sov[['meta']],
                             n.resamples = n.resamples,
                             chr = chr,
                             qtl_type = qtl_type,
                             # bootstrap_type = bootstrap_type,
                             seed = random.seed,
                             n.cores = n.cores,
                             silent = silent)

  sov <- list(meta = sov[['meta']],
              bootstrap_maxes = result)

  class(sov) <- c('scanonevar', class(sov))
  return(sov)
}



scanonevar.boot_ <- function(sov,
                             meta,
                             n.resamples,
                             chr,
                             qtl_type,
                             # bootstrap_type,
                             seed,
                             n.cores,
                             silent) {



  # if (n.cores != 1) {
  #   cl <- parallel::makeCluster(spec = n.cores)
  #   doParallel::registerDoParallel(cl = cl)
  # }

  max_positions <- rep(NA, n.resamples)
  c <- qtl:::subset.cross(x = meta$cross, chr = chr)

  for (resample_idx in 1:n.resamples) {

    bs_scan <- scanonevar(cross = qtl:::subset.cross(x = c, ind = sample(x = 1:qtl::nind(c), replace = TRUE)),
                          mean.formula = meta$scan.formulae$mean.alt.formula,
                          var.formula = meta$scan.formulae$var.alt.formula,
                          chrs = chr,
                          scan_types = qtl_type)

    max_positions[resample_idx] <- bs_scan$result$pos[which.max(x = bs_scan$resul[[paste0(qtl_type, '.lod')]])]

  }

  # if (n.cores != 1) {
  #   parallel::stopCluster(cl)
  # }

  # perms <- dplyr::bind_rows(perms)
  # sov <- calc.empir.ps(sov, perms)

  # return(list(sov = sov,
  #             perms = perms))

  return(max_positions)
}

