context("Testing scanonevar")

test_that(
  desc = 'testing scanonevar',
  code = {
    set.seed(27599)
    test.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(20, 5)))
    test.cross[['pheno']][['sex']] <- sample(x = c(0, 1),
                                             size = qtl::nind(test.cross),
                                             replace = TRUE)


    sov <- scanonevar(cross = test.cross,
                      mean.formula = phenotype ~ sex + D1M2 + mean.QTL.add + mean.QTL.dom,
                      var.formula = ~ sex + D2M3 + var.QTL.add + var.QTL.dom)

    # should be a scanonevar object
    expect_true(is.scanonevar(sov))

    # joint lods should be pointwise higher than mean and variance lods
    expect_true(all(sov$result %>% .[['joint.lod']] >= sov$result %>% .[['mean.lod']]))
    expect_true(all(sov$result %>% .[['joint.lod']] >= sov$result %>% .[['var.lod']]))


    x <- summary(object = sov, units = 'LOD')

    # components of summary
    expect_true(all(c('intro', 'mean.peaks', 'var.peaks', 'joint.peaks') %in% names(x)))

    # absolute maximum size of each component of summary is size of result
    expect_true(nrow(x$mean.peaks) <= nrow(sov$result))
    expect_true(nrow(x$var.peaks) <= nrow(sov$result))
    expect_true(nrow(x$joint.peaks) <= nrow(sov$result))

    # should be no randomness in summary -- recalculate it and it's the same
    expect_identical(object = x, expected = summary(object = sov, units = 'LOD'))

  }
)


