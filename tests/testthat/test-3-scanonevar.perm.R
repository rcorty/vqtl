context("Testing scanonevar.perm")

set.seed(27599)
test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(40, 5), n.mar = 10, eq.spacing = TRUE, include.x = FALSE),
                             n.ind = 300,
                             type = 'f2')
test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                         size = qtl::nind(test.cross),
                                         replace = TRUE)
test.cross[['pheno']][['phenotype1']] <- rnorm(n = qtl::nind(test.cross))
test.cross[['pheno']][['phenotype2']] <- rnorm(n = qtl::nind(test.cross), mean = 0.3*test.cross$geno$`2`$data[,3])
test.cross[['pheno']][['phenotype4']] <- rnorm(n = qtl::nind(test.cross),
                                               mean = test.cross$pheno$sex + test.cross$geno$`2`$data[,3],
                                               sd = exp(test.cross$pheno$sex + test.cross$geno$`3`$data[,3]))


x <- scanonevar(cross = test.cross,
                mean.formula = phenotype2 ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
                var.formula = ~ sex + D2M3 + var.QTL.add + var.QTL.dom)

y1a <- scanonevar.perm(sov = x, n.perms = 10, random.seed = 27599)
y1b <- scanonevar.perm(sov = x, n.perms = 10, random.seed = 27599)
y1c <- scanonevar.perm(sov = x, n.perms = 10, random.seed = 27599, n.cores = 3)
y1d <- scanonevar.perm(sov = x, n.perms = 10, random.seed = 27599, n.cores = 3)


test_that(
  desc = 'testing scanonevar.perm',
  code = {

    # all the results should be scanones
    expect_true(is.scanonevar(y1a))
    expect_true(is.scanonevar(y1b))
    expect_true(is.scanonevar(y1c))
    expect_true(is.scanonevar(y1d))

    # perms done with the same number of cores should be identical
    expect_true(identical(y1a$perms, y1b$perms))
    expect_true(identical(y1c$perms, y1d$perms))

    # # joint lods should be pointwise be higher than mean and var lods over the genome, so its max should be higher too
    # # I no longer think this is true -- the models are not nested
    # all.joint.lod.maxes.greater.than.x <- function(sov, x) {
    #   all(sov$perms %>% dplyr::filter(test == 'joint') %>% .[['max.lod']] > sov$perms %>% dplyr::filter(test == x) %>% .[['max.lod']])
    # }
    #
    # expect_true(all.joint.lod.maxes.greater.than.x(y1a, 'mean'))
    # expect_true(all.joint.lod.maxes.greater.than.x(y1a, 'var'))
    #
    # expect_true(all.joint.lod.maxes.greater.than.x(y1b, 'mean'))
    # expect_true(all.joint.lod.maxes.greater.than.x(y1b, 'var'))
    #
    # expect_true(all.joint.lod.maxes.greater.than.x(y1c, 'mean'))
    # expect_true(all.joint.lod.maxes.greater.than.x(y1c, 'var'))
    #
    # expect_true(all.joint.lod.maxes.greater.than.x(y1d, 'mean'))
    # expect_true(all.joint.lod.maxes.greater.than.x(y1d, 'var'))


    # it should not be the case that all the perms of a give type are the same
    all.lod.maxes.same <- function(sov, type) {
      max.lods <- sov$perms %>% dplyr::filter(test == type) %>% .[['max.lod']]
      return(min(max.lods, na.rm = TRUE) == max(max.lods, na.rm = TRUE))
    }

    expect_false(all.lod.maxes.same(sov = y1a, 'mean'))
    expect_false(all.lod.maxes.same(sov = y1a, 'var'))
    expect_false(all.lod.maxes.same(sov = y1a, 'joint'))

    expect_false(all.lod.maxes.same(sov = y1c, 'mean'))
    expect_false(all.lod.maxes.same(sov = y1c, 'var'))
    expect_false(all.lod.maxes.same(sov = y1c, 'joint'))
  }
)



test_that(
  desc = 'testing c.scanonevar.perm',
  code = {

    expect_true(is.scanonevar(c(y1a, y1b)))
    expect_true(is.scanonevar(c(y1a, y1c)))
    expect_true(is.scanonevar(c(y1a, y1d)))

    expect_true(is.scanonevar(c(y1c, y1a)))
    expect_true(is.scanonevar(c(y1c, y1b)))
    expect_true(is.scanonevar(c(y1c, y1d)))

    expect_true(is.scanonevar(c(y1a, y1b, y1c, y1d)))
  }
)


test_that(
  desc = 'summary.scanonevar',
  code = {

    a <- summary(object = y1a, units = 'empir.p')
    b <- summary(object = y1b, units = 'empir.p')
    c <- summary(object = y1c, units = 'empir.p')
    d <- summary(object = y1d, units = 'empir.p')

    # summaries should be the same if the objects are the same
    expect_identical(object = a, expected = b)
    expect_identical(object = d, expected = d)

    # components of summary
    expect_true(all(c('intro', 'mean.peaks', 'var.peaks', 'joint.peaks') %in% names(a)))
    expect_true(all(c('intro', 'mean.peaks', 'var.peaks', 'joint.peaks') %in% names(c)))

    # size of summary
    expect_true(nrow(a$mean.peaks) <= nrow(y1a$result))
    expect_true(nrow(a$var.peaks) <= nrow(y1a$result))
    expect_true(nrow(a$joint.peaks) <= nrow(y1a$result))

    expect_true(nrow(c$mean.peaks) <= nrow(y1c$result))
    expect_true(nrow(c$var.peaks) <= nrow(y1c$result))
    expect_true(nrow(c$joint.peaks) <= nrow(y1c$result))

  }
)