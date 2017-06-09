#' @title c.scanonevar
#' @name c.scanonevar
#'
#' @param ... the scanonevar objects with permutations to be combined
#'
#' @description combines scanonevar objects that have permutations to improve the precision of the p-value estimates.
#'
#' @return a scanonevar object that is the concatenation of the inputted
#' scanonevars
#'
#' @export
#'
c.scanonevar <- function(...) {
  sovs <- list(...)

  validate.c.scanonevar.input_(sovs)

  first.sov <- sovs[[1]]
  new.perms <- first.sov[['perms']]
  for (sov.idx in 2:(length(sovs))) {
    new.perms <- dplyr::bind_rows(new.perms,
                                  sovs[[sov.idx]][['perms']])
  }
  first.sov[['perms']] <- new.perms

  first.sov[['result']] <- calc.empir.ps(sov = first.sov[['result']],
                                         perms = new.perms)

  return(first.sov)

}


validate.c.scanonevar.input_ <- function(sovs) {

  stopifnot(all(sapply(X = sovs, FUN = is.scanonevar.w.perms)))

  first.sov <- sovs[[1]]
  for (sov.idx in 2:(length(sovs))) {
    # meta has to be the same
    stopifnot(identical(x = first.sov[['meta']], y = sovs[[sov.idx]][['meta']],
                        ignore.environment = TRUE))

    # result has to be the same
    # todo

    # perms have to have the same names...but shouldn't be the same
    stopifnot(all.equal(names(first.sov[['perms']]),
                        names(sovs[[sov.idx]][['perms']])))
  }


}
